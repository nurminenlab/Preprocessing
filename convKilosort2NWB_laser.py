from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path
import numpy as np
from pynwb import NWBHDF5IO, NWBFile
from pynwb.misc import Units
from pynwb.file import LabMetaData
from neuroconv.datainterfaces import KiloSortSortingInterface
import spikeinterface.extractors as se
from spikeinterface import read_binary
import pandas as pd
import scipy.io as sio
import glob

from datetime import date

def calculate_age_days(birthday: date, reference_date: date = None) -> int:
    """
    Calculate age in days from birthday to reference date.
    
    Parameters:
    -----------
    birthday : date
        Date of birth
    reference_date : date, optional
        Reference date for age calculation. If None, uses today's date.
        
    Returns:
    --------
    int
        Age in days
    """
    if reference_date is None:
        reference_date = date.today()
    
    return (reference_date - birthday).days


def parse_date_string(date_str: str) -> date:
    """
    Convert a date string in 'YYYY-MM-DD' format to a date object.
    
    Parameters:
    -----------
    date_str : str
        Date string in format 'YYYY-MM-DD' (e.g., '2025-12-05')
        
    Returns:
    --------
    date
        Python date object
    """
    year, month, day = map(int, date_str.split('-'))
    return date(year, month, day)


def dataframe_to_nwb_trials(nwbfile, df):
    """
    Add trials from a pandas DataFrame to an NWB file's trials table.
    
    Parameters:
    -----------
    nwbfile : NWBFile
        The NWB file object to which trials will be added.
    df : pd.DataFrame
        DataFrame containing trial information. Must have 'start_time' and 'stop_time' columns.
    """
    for col in df.columns:
        if col not in ['start_time', 'stop_time']:
            nwbfile.add_trial_column(name=col, description=f'Trial column: {col}')
    
    for _, row in df.iterrows():
        nwbfile.add_trial(
            start_time=row['start_time'],
            stop_time=row['stop_time'],
            **{col: row[col] for col in df.columns if col not in ['start_time', 'stop_time']}
        )

def linear_stimuli(expt_info,stim_type='orientation_ENV'):
    stims = np.array([])
    ENVS = np.array([])
    if stim_type=='orientation_ENV':
        for i in range(expt_info.trial_records.shape[0]):
            if ~np.isnan(expt_info.trial_records[i].tr_orientation).all():
                n_ori = expt_info.trial_records[i].tr_orientation.shape[0]
            else:
                n_ori = 1
            stims =  np.hstack((stims, expt_info.trial_records[i].tr_orientation))
            ENVS = np.hstack((ENVS, np.repeat(expt_info.trial_records[i].ENV, n_ori)))
    elif stim_type=='image':
        for i in range(expt_info.trial_records.shape[0]):
            stims = np.hstack((stims, expt_info.trial_records[i].trImage))
    return stims, ENVS

def parse_spikeglx_meta(meta_path):
    """
    Parse SpikeGLX .meta file and extract recording parameters.
    
    Parameters:
    -----------
    meta_path : str or Path
        Path to the .ap.meta file
        
    Returns:
    --------
    meta_dict : dict
        Dictionary containing metadata parameters
    """
    meta_dict = {}
    with open(meta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if '=' in line:
                key, value = line.split('=', 1)
                meta_dict[key.strip()] = value.strip()
    return meta_dict


def get_sampling_rate_from_meta(meta_path):
    """
    Extract sampling rate from SpikeGLX .meta file.
    
    Parameters:
    -----------
    meta_path : str or Path
        Path to the .ap.meta file
        
    Returns:
    --------
    sampling_rate : float
        Sampling rate in Hz
    """
    meta = parse_spikeglx_meta(meta_path)
    # SpikeGLX stores sampling rate as 'imSampRate' for imec probes
    if 'imSampRate' in meta:
        return float(meta['imSampRate'])
    # Fallback to 'niSampRate' for NI-DAQ
    elif 'niSampRate' in meta:
        return float(meta['niSampRate'])
    else:
        raise ValueError(f"Could not find sampling rate in meta file: {meta_path}")


def get_num_channels_from_meta(meta_path):
    """
    Extract number of saved channels from SpikeGLX .meta file.
    
    Parameters:
    -----------
    meta_path : str or Path
        Path to the .ap.meta or .nidq.meta file
        
    Returns:
    --------
    num_channels : int
        Number of saved channels
    """
    meta = parse_spikeglx_meta(meta_path)
    if 'nSavedChans' in meta:
        return int(meta['nSavedChans'])
    else:
        raise ValueError(f"Could not find nSavedChans in meta file: {meta_path}")


def get_dtype_from_meta(meta_path):
    """
    Extract data type from SpikeGLX .meta file.
    
    Parameters:
    -----------
    meta_path : str or Path
        Path to the .ap.meta or .nidq.meta file
        
    Returns:
    --------
    dtype : str
        Data type string (e.g., 'int16')
    """
    meta = parse_spikeglx_meta(meta_path)
    # SpikeGLX typically stores data as int16
    # Check if there's a typeThis field (some versions have it)
    if 'typeThis' in meta:
        type_code = meta['typeThis']
        # Type codes: 0=int16, 1=uint16, 2=int32, 3=uint32, 4=float32, 5=float64
        type_map = {0: 'int16', 1: 'uint16', 2: 'int32', 3: 'uint32', 4: 'float32', 5: 'float64'}
        return type_map.get(int(type_code), 'int16')
    # Default to int16 as SpikeGLX always saves raw data as int16
    return 'int16'


def get_uv_per_bit(meta_path):
    """
    Get the conversion factor from raw int16 values to microvolts for AP data.
    
    For Neuropixels probes, the formula is:
        voltage (uV) = raw_value * Vmax / (maxInt * AP_gain) * 1e6
    
    Parameters:
    -----------
    meta_path : str or Path
        Path to the .ap.meta file
        
    Returns:
    --------
    uv_per_bit : float
        Conversion factor: multiply raw int16 values by this to get microvolts
    """
    meta = parse_spikeglx_meta(meta_path)
    
    # Get voltage range
    Vmax = float(meta.get('imAiRangeMax', 0.6))
    
    # Determine maxInt from probe type
    # NP 1.0 (type 0): 10-bit ADC -> maxInt = 512
    # NP 2.0 (type 21, 24): 14-bit ADC -> maxInt = 8192
    prb_type = int(meta.get('imDatPrb_type', 0))
    if prb_type in [21, 24, 2013, 2014]:
        maxInt = 8192
    else:
        maxInt = 512  # NP 1.0 default
    
    # Get AP gain from imroTbl
    # Format: (probe_type,num_channels)(chan bank refid apgain lfgain)...
    imro_str = meta.get('imroTbl', '')
    entries = imro_str.split(')')
    ap_gain = 500  # default
    if len(entries) > 1:
        first_channel = entries[1].strip().lstrip('(')
        parts = first_channel.split()
        if len(parts) >= 4:
            ap_gain = float(parts[3])
    
    uv_per_bit = (Vmax / (maxInt * ap_gain)) * 1e6
    return uv_per_bit


def get_last_analog_channel_index(nidaq_meta_path):
    """
    Get the index of the last analog channel from the NIDAQ meta file.
    
    Parses the 'acqMnMaXaDw' field which gives counts of:
    MN (multiplexed neural), MA (multiplexed aux), XA (extra analog), DW (digital word).
    The last analog channel index = MN + MA + XA - 1.
    
    Parameters:
    -----------
    nidaq_meta_path : str or Path
        Path to the .nidq.meta file
        
    Returns:
    --------
    int
        Zero-based index of the last analog channel
    """
    meta = parse_spikeglx_meta(nidaq_meta_path)
    acq_str = meta.get('acqMnMaXaDw', '')
    if acq_str:
        parts = acq_str.split(',')
        mn = int(parts[0])
        ma = int(parts[1])
        xa = int(parts[2])
        total_analog = mn + ma + xa
        return total_analog - 1
    else:
        # Fallback: assume last saved channel is digital, second-to-last is last analog
        n_saved = int(meta.get('nSavedChans', 8))
        return n_saved - 2


def get_session_start_time_from_meta(meta_path):
    """
    Extract session start time from SpikeGLX .meta file.
    
    Parameters:
    -----------
    meta_path : str or Path
        Path to the .ap.meta or .nidq.meta file
        
    Returns:
    --------
    datetime
        Session start time as a timezone-aware datetime object
    """
    meta = parse_spikeglx_meta(meta_path)
    # SpikeGLX stores the file creation time in 'fileCreateTime' field
    # Format is typically: 2025-12-05T14:30:45
    if 'fileCreateTime' in meta:
        time_str = meta['fileCreateTime']
        # Parse the datetime string
        dt = datetime.fromisoformat(time_str)
        # Add timezone if not present (assume local time, use US/Central as default)
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=ZoneInfo("US/Central"))
        return dt
    else:
        raise ValueError(f"Could not find fileCreateTime in meta file: {meta_path}")


def load_spikeglx_binary(bin_path, num_channels, sampling_rate, dtype='int16'):
    """Load SpikeGLX binary file as a numpy memmap."""
    data = np.memmap(bin_path, dtype=dtype, mode='r')
    num_samples = len(data) // num_channels
    data = data.reshape((num_samples, num_channels))
    return data


def extract_waveforms(data, spike_times, channel_idx, sampling_rate, 
                      pre_samples=30, post_samples=60, max_waveforms=200):
    """
    Extract waveforms around spike times from the raw data.
    
    Parameters:
    -----------
    data : np.ndarray
        Raw data array (samples x channels)
    spike_times : np.ndarray
        Spike times in samples
    channel_idx : int
        Channel index to extract waveforms from
    sampling_rate : float
        Sampling rate in Hz
    pre_samples : int
        Number of samples before spike to include
    post_samples : int
        Number of samples after spike to include
    max_waveforms : int
        Maximum number of waveforms to extract for sample
        
    Returns:
    --------
    mean_waveform : np.ndarray
        Mean waveform across all spikes
    std_waveform : np.ndarray
        Standard deviation of waveform across all spikes
    sample_waveforms : np.ndarray
        Sample of waveforms (up to max_waveforms)
    snr : float
        Signal-to-noise ratio
    ap_type : str
        'somatic' or 'axonal' based on waveform characteristics
    trough_to_peak_ms : float
        Action potential duration (trough to peak time) in milliseconds
    """
    waveform_length = pre_samples + post_samples
    num_samples = data.shape[0]
    
    # Filter valid spike times (those that have enough samples before and after)
    valid_spikes = spike_times[(spike_times >= pre_samples) & 
                               (spike_times < num_samples - post_samples)]
    
    if len(valid_spikes) == 0:
        return np.zeros(waveform_length), np.zeros(waveform_length), np.zeros((1, waveform_length)), 0.0, 'unknown', 0.0
    
    # Extract all waveforms for mean calculation
    all_waveforms = []
    for spike_idx in valid_spikes:
        start_idx = int(spike_idx - pre_samples)
        end_idx = int(spike_idx + post_samples)
        waveform = data[start_idx:end_idx, channel_idx].astype(np.float32)
        all_waveforms.append(waveform)
    
    all_waveforms = np.array(all_waveforms)
    mean_waveform = np.mean(all_waveforms, axis=0)
    std_waveform = np.std(all_waveforms, axis=0)
    
    # Sample waveforms (random subset)
    if len(valid_spikes) > max_waveforms:
        sample_indices = np.random.choice(len(all_waveforms), max_waveforms, replace=False)
        sample_waveforms = all_waveforms[sample_indices]
    else:
        sample_waveforms = all_waveforms
    
    # Calculate SNR: peak-to-peak amplitude / (2 * std of noise)
    # Noise is estimated from the first few samples (baseline)
    baseline_samples = pre_samples // 2
    noise_std = np.std(all_waveforms[:, :baseline_samples])
    peak_to_peak = np.max(mean_waveform) - np.min(mean_waveform)
    snr = peak_to_peak / (2 * noise_std) if noise_std > 0 else 0.0
    
    # Calculate trough-to-peak duration and classify AP type
    # Find the trough (minimum) in the waveform
    trough_idx = np.argmin(mean_waveform)
    
    # Find the peak (maximum) after the trough
    # This is the repolarization peak that follows the negative deflection
    if trough_idx < len(mean_waveform) - 1:
        peak_idx = trough_idx + np.argmax(mean_waveform[trough_idx:])
    else:
        peak_idx = trough_idx
    
    # Calculate trough-to-peak duration in milliseconds
    trough_to_peak_samples = peak_idx - trough_idx
    trough_to_peak_ms = (trough_to_peak_samples / sampling_rate) * 1000.0
    
    # Classify AP type based on trough-to-peak duration
    # Somatic (regular spiking, pyramidal): typically > 0.4 ms (broader waveforms)
    # Axonal (fast spiking, interneurons): typically <= 0.4 ms (narrower waveforms)
    # This threshold is commonly used in the literature (e.g., Bartho et al., 2004)
    if trough_to_peak_ms > 0.4:
        ap_type = 'broad'
    else:
        ap_type = 'narrow'
    
    return mean_waveform, std_waveform, sample_waveforms, snr, ap_type, trough_to_peak_ms


def load_nidaq_binary(bin_path, num_channels=8, dtype='int16'):
    """
    Load NIDAQ binary file as a numpy memmap.
    
    Parameters:
    -----------
    bin_path : str or Path
        Path to the .nidq.bin file
    num_channels : int
        Number of channels in the recording (default 8)
    dtype : str
        Data type of the binary file
        
    Returns:
    --------
    data : np.ndarray
        Data array (samples x channels)
    """
    data = np.memmap(bin_path, dtype=dtype, mode='r')
    num_samples = len(data) // num_channels
    data = data.reshape((num_samples, num_channels))
    return data


def detect_ttl_onsets(ttl_signal, threshold=None):
    """
    Detect onset times (rising edges) of TTL signals.
    
    Parameters:
    -----------
    ttl_signal : np.ndarray
        1D array of TTL signal values
    threshold : float, optional
        Threshold for detecting high state. If None, uses midpoint of signal range.
        
    Returns:
    --------
    onset_indices : np.ndarray
        Indices where TTL signal transitions from low to high
    """
    if threshold is None:
        threshold = (np.max(ttl_signal) + np.min(ttl_signal)) / 2
    
    # Binarize the signal
    binary_signal = (ttl_signal > threshold).astype(int)
    
    # Find rising edges (0 to 1 transitions)
    diff_signal = np.diff(binary_signal)
    onset_indices = np.where(diff_signal == 1)[0] + 1
    
    return onset_indices


def detect_ttl_offsets(ttl_signal, threshold=None):
    """
    Detect offset times (falling edges) of TTL signals.
    
    Parameters:
    -----------
    ttl_signal : np.ndarray
        1D array of TTL signal values
    threshold : float, optional
        Threshold for detecting high state. If None, uses midpoint of signal range.
        
    Returns:
    --------
    offset_indices : np.ndarray
        Indices where TTL signal transitions from high to low
    """
    if threshold is None:
        threshold = (np.max(ttl_signal) + np.min(ttl_signal)) / 2
    
    # Binarize the signal
    binary_signal = (ttl_signal > threshold).astype(int)
    
    # Find falling edges (1 to 0 transitions)
    diff_signal = np.diff(binary_signal)
    offset_indices = np.where(diff_signal == -1)[0] + 1
    
    return offset_indices


def is_ttl_high(ttl_signal, sample_idx, threshold=None):
    """
    Check if TTL signal is high at a given sample index.
    
    Parameters:
    -----------
    ttl_signal : np.ndarray
        1D array of TTL signal values
    sample_idx : int
        Sample index to check
    threshold : float, optional
        Threshold for detecting high state
        
    Returns:
    --------
    bool
        True if TTL is high at the given sample
    """
    if threshold is None:
        threshold = (np.max(ttl_signal) + np.min(ttl_signal)) / 2
    
    if sample_idx < 0 or sample_idx >= len(ttl_signal):
        return False
    
    return ttl_signal[sample_idx] > threshold


def load_experiment_mat_file(nidaq_bin_path):
    """
    Find and load the .mat file in the same folder as the NIDAQ binary file.
    
    Parameters:
    -----------
    nidaq_bin_path : str or Path
        Path to the NIDAQ .bin file
        
    Returns:
    --------
    mat_data : dict
        Loaded .mat file contents
    """
    nidaq_folder = Path(nidaq_bin_path).parent
    mat_files = list(nidaq_folder.glob('*.mat'))
    
    if len(mat_files) == 0:
        raise FileNotFoundError(f"No .mat file found in {nidaq_folder}")
    elif len(mat_files) > 1:
        print(f"Warning: Multiple .mat files found, using: {mat_files[0]}")
    
    mat_path = mat_files[0]
    print(f"Loading experiment info from: {mat_path}")
    
    # Load mat file (try both old and new format)
    try:
        mat_data = sio.loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    except NotImplementedError:
        # For MATLAB v7.3 files, use h5py
        import h5py
        mat_data = {}
        with h5py.File(mat_path, 'r') as f:
            for key in f.keys():
                mat_data[key] = f[key][()]
    
    return mat_data


def expt_info_to_metadata(expt_info):
    """
    Convert expt_info structure to a dictionary suitable for NWB metadata.
    Excludes trial_records as those go in the trials table.
    
    Parameters:
    -----------
    expt_info : object
        Experiment info structure from mat file
        
    Returns:
    --------
    metadata_dict : dict
        Dictionary of experiment metadata
    """
    metadata_dict = {}
    
    # Get all attributes except trial_records
    if hasattr(expt_info, '_fieldnames'):
        # For scipy.io.loadmat with struct_as_record=False
        field_names = expt_info._fieldnames
    else:
        # Try to get attributes
        field_names = [attr for attr in dir(expt_info) 
                       if not attr.startswith('_') and not callable(getattr(expt_info, attr))]
    
    for field_name in field_names:
        if field_name == 'trial_records': # processed elsewhere
            continue
        
        try:
            value = getattr(expt_info, field_name)
            
            # Convert numpy arrays to lists for JSON serialization
            if isinstance(value, np.ndarray):
                value = value.tolist()
            elif hasattr(value, '_fieldnames'):
                # Nested struct - convert to dict
                value = str(value)  # Simple string representation
            
            metadata_dict[field_name] = value
        except Exception as e:
            print(f"Warning: Could not extract field {field_name}: {e}")
    
    return metadata_dict


def get_unit_depths(kilosort_folder):
    """
    Load unit depths from Kilosort output files.
    
    Parameters:
    -----------
    kilosort_folder : str or Path
        Path to the Kilosort output folder
        
    Returns:
    --------
    unit_depths : dict
        Dictionary mapping unit_id to depth
    """
    kilosort_folder = Path(kilosort_folder)
    
    # Load cluster info if available
    cluster_info_path = kilosort_folder / "cluster_info.tsv"
    if cluster_info_path.exists():
        cluster_info = pd.read_csv(cluster_info_path, sep='\t')
        if 'depth' in cluster_info.columns and 'cluster_id' in cluster_info.columns:
            return dict(zip(cluster_info['cluster_id'], cluster_info['depth']))
    
    # Alternative: compute from channel positions and spike templates
    channel_positions_path = kilosort_folder / "channel_positions.npy"
    spike_clusters_path = kilosort_folder / "spike_clusters.npy"
    spike_templates_path = kilosort_folder / "spike_templates.npy"
    templates_path = kilosort_folder / "templates.npy"
    
    if channel_positions_path.exists() and templates_path.exists():
        channel_positions = np.load(channel_positions_path)
        templates = np.load(templates_path)  # (n_templates, n_timepoints, n_channels)
        
        # Get the channel with maximum amplitude for each template
        template_amplitudes = np.max(templates, axis=1) - np.min(templates, axis=1)
        best_channels = np.argmax(template_amplitudes, axis=1)
        
        # Get depths for each template
        template_depths = channel_positions[best_channels, 1]  # y-coordinate is depth
        
        # Map clusters to templates
        if spike_clusters_path.exists():
            spike_clusters = np.load(spike_clusters_path)
            unique_clusters = np.unique(spike_clusters)
            
            # If we have spike_templates, use it to map clusters to templates
            if spike_templates_path.exists():
                spike_templates = np.load(spike_templates_path)
                cluster_depths = {}
                for cluster_id in unique_clusters:
                    cluster_mask = spike_clusters == cluster_id
                    cluster_templates = spike_templates[cluster_mask]
                    # Most common template for this cluster
                    if len(cluster_templates) > 0:
                        most_common_template = np.bincount(cluster_templates.flatten()).argmax()
                        if most_common_template < len(template_depths):
                            cluster_depths[cluster_id] = template_depths[most_common_template]
                return cluster_depths
            else:
                # Assume cluster_id == template_id
                return {i: template_depths[i] for i in range(len(template_depths))}
    
    return {}


def get_unit_channels(kilosort_folder):
    """
    Get the best channel for each unit from Kilosort output.
    
    Returns:
    --------
    unit_channels : dict
        Dictionary mapping unit_id to best channel index
    """
    kilosort_folder = Path(kilosort_folder)
    
    # Try cluster_info.tsv first
    cluster_info_path = kilosort_folder / "cluster_info.tsv"
    if cluster_info_path.exists():
        cluster_info = pd.read_csv(cluster_info_path, sep='\t')
        if 'ch' in cluster_info.columns and 'cluster_id' in cluster_info.columns:
            return dict(zip(cluster_info['cluster_id'], cluster_info['ch']))
    
    # Alternative: compute from templates
    templates_path = kilosort_folder / "templates.npy"
    spike_clusters_path = kilosort_folder / "spike_clusters.npy"
    spike_templates_path = kilosort_folder / "spike_templates.npy"
    
    if templates_path.exists():
        templates = np.load(templates_path)
        template_amplitudes = np.max(templates, axis=1) - np.min(templates, axis=1)
        best_channels = np.argmax(template_amplitudes, axis=1)
        
        if spike_clusters_path.exists() and spike_templates_path.exists():
            spike_clusters = np.load(spike_clusters_path)
            spike_templates = np.load(spike_templates_path)
            unique_clusters = np.unique(spike_clusters)
            
            cluster_channels = {}
            for cluster_id in unique_clusters:
                cluster_mask = spike_clusters == cluster_id
                cluster_templates = spike_templates[cluster_mask]
                if len(cluster_templates) > 0:
                    most_common_template = np.bincount(cluster_templates.flatten()).argmax()
                    if most_common_template < len(best_channels):
                        cluster_channels[cluster_id] = best_channels[most_common_template]
            return cluster_channels
        else:
            return {i: best_channels[i] for i in range(len(best_channels))}
    
    return {}


def _pad_true_runs(mask: np.ndarray, pad_before: int, pad_after: int) -> np.ndarray:
    """Expand each contiguous True-run in a bool mask by pad_before / pad_after samples."""
    if pad_before <= 0 and pad_after <= 0:
        return mask.copy()
    n = mask.size
    out = mask.copy()
    edges = np.diff(mask.astype(np.int8))
    starts = (np.where(edges == 1)[0] + 1).tolist()
    ends = (np.where(edges == -1)[0] + 1).tolist()  # exclusive
    if mask[0]:
        starts = [0] + starts
    if mask[-1]:
        ends = ends + [n]
    for s, e in zip(starts, ends):
        out[max(0, s - pad_before):min(n, e + pad_after)] = True
    return out


def _drop_short_valid_islands(bad_mask: np.ndarray, min_valid_samples: int) -> np.ndarray:
    """Mark contiguous valid (False) runs shorter than min_valid_samples as bad."""
    n = bad_mask.size
    if n == 0 or min_valid_samples <= 1:
        return bad_mask.copy()
    valid = ~bad_mask
    edges = np.diff(valid.astype(np.int8))
    starts = (np.where(edges == 1)[0] + 1).tolist()
    ends = (np.where(edges == -1)[0] + 1).tolist()
    if valid[0]:
        starts = [0] + starts
    if valid[-1]:
        ends = ends + [n]
    out = bad_mask.copy()
    for s, e in zip(starts, ends):
        if (e - s) < min_valid_samples:
            out[s:e] = True
    return out


def clean_pupil_blinks(
    pupil: np.ndarray,
    sampling_rate: float,
    speed_mad_n: float = 16.0,
    pad_before_ms: float = 50.0,
    pad_after_ms: float = 150.0,
    min_valid_ms: float = 20.0,
    invalid_low_frac: float = 0.10,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Remove blink artifacts from a pupil-diameter trace and linearly interpolate.

    Follows the standard pupillometry preprocessing recipe:

    * Kret & Sjak-Shie (2019), Behav. Res. Methods 51:1336-1352 -- "speed
      filter": per-sample dilation speed = max(|forward diff|, |backward diff|)
      / dt; samples with speed > median + n * MAD (n=16 by default) are marked
      as artifacts. MAD is used instead of SD because it is robust to the
      blinks themselves.
    * Mathot (2013) / Hershman et al. (2018) -- pad each blink run by
      pad_before_ms / pad_after_ms because the eyelid partially occludes the
      pupil around the blink core. Discard "recovered" islands shorter than
      min_valid_ms since they are usually still corrupted.
    * Signal-loss floor: full closures drive the analog eye-tracker channel
      close to zero; anything at or below invalid_low_frac * median(positive
      pupil) is flagged.
    * Gaps are then linearly interpolated. Leading / trailing edge gaps stay
      NaN (no extrapolation).

    Parameters
    ----------
    pupil : np.ndarray
        1-D pupil diameter samples.
    sampling_rate : float
        Sample rate of `pupil`, in Hz.
    speed_mad_n, pad_before_ms, pad_after_ms, min_valid_ms, invalid_low_frac
        Tuning knobs; defaults follow Kret & Sjak-Shie (2019) recommendations
        adjusted for a slightly more conservative post-blink pad.

    Returns
    -------
    cleaned : np.ndarray of float32
        Blink-cleaned pupil trace, same length as input.
    blink_mask : np.ndarray of bool
        True where the sample was flagged as blink and interpolated (or set to
        NaN at the edges).
    """
    pupil = np.asarray(pupil, dtype=np.float32).copy()
    n = pupil.size
    if n == 0:
        return pupil, np.zeros(0, dtype=bool)

    dt = 1.0 / float(sampling_rate)

    # 1) Signal-loss floor (full-closure blinks)
    positive = pupil[pupil > 0]
    if positive.size:
        low_thr = invalid_low_frac * np.median(positive)
    else:
        low_thr = 0.0
    low_mask = pupil <= low_thr

    # 2) Kret & Sjak-Shie dilation-speed filter
    if n >= 2:
        diffs = np.abs(np.diff(pupil)) / dt
        speed = np.zeros(n, dtype=np.float32)
        speed[:-1] = diffs
        speed[1:] = np.maximum(speed[1:], diffs)  # max(forward, backward)

        med_s = np.median(speed)
        mad_s = np.median(np.abs(speed - med_s))
        # If MAD collapses to 0 (flat regions dominate), fall back to std so
        # the threshold is still meaningful.
        scale = mad_s if mad_s > 0 else float(np.std(speed))
        speed_thr = med_s + speed_mad_n * scale
        speed_mask = speed > speed_thr
    else:
        speed_mask = np.zeros(n, dtype=bool)

    bad_mask = low_mask | speed_mask

    # 3) Discard tiny valid islands wedged inside blink clusters
    min_valid_samples = max(1, int(round(min_valid_ms * 1e-3 * sampling_rate)))
    bad_mask = _drop_short_valid_islands(bad_mask, min_valid_samples)

    # 4) Pad artifact runs
    pad_before = int(round(pad_before_ms * 1e-3 * sampling_rate))
    pad_after = int(round(pad_after_ms * 1e-3 * sampling_rate))
    bad_mask = _pad_true_runs(bad_mask, pad_before, pad_after)

    # 5) Linear interpolation over interior gaps
    good_idx = np.where(~bad_mask)[0]
    cleaned = np.full(n, np.nan, dtype=np.float32)
    if good_idx.size >= 2:
        first, last = int(good_idx[0]), int(good_idx[-1])
        interior = np.arange(first, last + 1)
        cleaned[interior] = np.interp(
            interior, good_idx, pupil[good_idx]
        ).astype(np.float32)
        # Preserve original values at good samples exactly
        cleaned[good_idx] = pupil[good_idx]
    elif good_idx.size == 1:
        cleaned[good_idx[0]] = pupil[good_idx[0]]

    return cleaned, bad_mask


def convert_session(
    monkey_name: str,
    date_str: str,
    birthday_str: str,
    sex: str,
    session_name: str = "BSD",
    gate: str = "g0",
    project: str = "Inhibition-neural-information-processing",
    data_root: str = "F:/localDATA/Electrophysiology",
    output_dir: str = None,
    species: str = "Callithrix jacchus",
    pre_samples: int = 30,
    post_samples: int = 60,
    max_sample_waveforms: int = 200,
    overwrite: bool = False,
) -> Path:
    """
    Run the full Kilosort -> NWB conversion for a single session.

    Parameters
    ----------
    monkey_name, date_str, birthday_str, sex
        Subject identifiers. `date_str` and `birthday_str` must be 'YYYY-MM-DD'.
    session_name, gate, project
        SpikeGLX session-name / gate / project identifiers used to construct
        both the source-data paths and the output file name.
    data_root
        Root of the raw Electrophysiology tree.
    output_dir
        If given, the .nwb file is written to this folder. Otherwise it is
        written next to the Kilosort output folder (original behavior).
    species
        Latin binomial for the NWB Subject metadata.
    pre_samples, post_samples, max_sample_waveforms
        Waveform extraction parameters.
    overwrite
        If True and the output file already exists, delete it before running.
        If False (default), an existing file will cause `run_conversion` to
        raise.

    Returns
    -------
    Path
        Full path to the written .nwb file.
    """
    import json

    file_name = date_str + "-" + monkey_name + "-" + project + ".nwb"

    recording_date = parse_date_string(date_str)
    birthday = parse_date_string(birthday_str)
    age_days = calculate_age_days(birthday, recording_date)

    subject_info = dict(
        subject_id=monkey_name,
        species=species,
        sex=sex,
        age=f"P{age_days}D",
    )

    # Kilosort output folder
    folder_path = (
        f"{data_root}/{monkey_name}/{date_str}/{session_name}_{gate}/"
        f"{session_name}_{gate}_imec0/kilosort4/"
    )

    # Paths to SpikeGLX .ap.bin / .ap.meta and NIDAQ .bin / .meta files
    spikeglx_bin_path = (
        f"{data_root}/{monkey_name}/{date_str}/{session_name}_{gate}/"
        f"{session_name}_{gate}_imec0/{session_name}_{gate}_t0.imec0.ap.bin"
    )
    spikeglx_meta_path = (
        f"{data_root}/{monkey_name}/{date_str}/{session_name}_{gate}/"
        f"{session_name}_{gate}_imec0/{session_name}_{gate}_t0.imec0.ap.meta"
    )
    NIDAQ_bin_path = (
        f"{data_root}/{monkey_name}/{date_str}/{session_name}_{gate}/"
        f"{session_name}_{gate}_t0.nidq.bin"
    )
    NIDAQ_meta_path = (
        f"{data_root}/{monkey_name}/{date_str}/{session_name}_{gate}/"
        f"{session_name}_{gate}_t0.nidq.meta"
    )

    # Read recording parameters from the .ap.meta file
    num_channels = get_num_channels_from_meta(spikeglx_meta_path)
    dtype = 'int16'  # SpikeGLX always saves raw data as int16
    sampling_rate = get_sampling_rate_from_meta(spikeglx_meta_path)

    # ============== CONVERSION ==============
    interface = KiloSortSortingInterface(folder_path=folder_path, verbose=False)

    metadata = interface.get_metadata()
    session_start_time = get_session_start_time_from_meta(spikeglx_meta_path)
    metadata["NWBFile"].update(
        session_start_time=session_start_time,
        institution="University of Houston, College of Optometry",
        lab="Neural information processing laboratory / nurminenlab",
    )

    # Add subject information (required for DANDI upload)
    metadata["Subject"] = subject_info

    # Decide where to write the .nwb file
    if output_dir is not None:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        nwbfile_path = str(Path(output_dir) / file_name)
    else:
        nwbfile_path = folder_path + file_name

    if overwrite and Path(nwbfile_path).exists():
        print(f"Overwrite requested; removing existing file: {nwbfile_path}")
        Path(nwbfile_path).unlink()

    # First, run the basic conversion to create the NWB file
    interface.run_conversion(nwbfile_path=nwbfile_path, metadata=metadata)

    # ============== ADD DEPTH, SNR, AND WAVEFORMS ==============
    print("Loading Kilosort data for additional processing...")

    unit_depths = get_unit_depths(folder_path)
    unit_channels = get_unit_channels(folder_path)

    spike_times_all = np.load(Path(folder_path) / "spike_times.npy").flatten()
    spike_clusters_all = np.load(Path(folder_path) / "spike_clusters.npy").flatten()

    raw_data = load_spikeglx_binary(spikeglx_bin_path, num_channels, sampling_rate, dtype)

    uv_per_bit = get_uv_per_bit(spikeglx_meta_path)

    with NWBHDF5IO(nwbfile_path, 'r+') as io:
        nwbfile = io.read()
        unit_ids = nwbfile.units.id[:]

        depths = []
        snrs = []
        mean_waveforms = []
        std_waveforms = []
        sample_waveforms_list = []
        ap_types = []
        trough_to_peak_durations = []

        for unit_id in unit_ids:
            depth = unit_depths.get(unit_id, np.nan)
            depths.append(depth)

            unit_spike_mask = spike_clusters_all == unit_id
            unit_spike_times = spike_times_all[unit_spike_mask]

            best_channel = unit_channels.get(unit_id, 0)

            if len(unit_spike_times) > 0 and best_channel < num_channels:
                mean_wf, std_wf, sample_wf, snr, ap_type, trough_to_peak_ms = extract_waveforms(
                    raw_data,
                    unit_spike_times,
                    best_channel,
                    sampling_rate,
                    pre_samples=pre_samples,
                    post_samples=post_samples,
                    max_waveforms=max_sample_waveforms,
                )
            else:
                waveform_length = pre_samples + post_samples
                mean_wf = np.zeros(waveform_length)
                std_wf = np.zeros(waveform_length)
                sample_wf = np.zeros((1, waveform_length))
                snr = 0.0
                ap_type = 'unknown'
                trough_to_peak_ms = 0.0

            mean_wf = mean_wf * uv_per_bit
            std_wf = std_wf * uv_per_bit
            sample_wf = sample_wf * uv_per_bit

            snrs.append(snr)
            mean_waveforms.append(mean_wf)
            std_waveforms.append(std_wf)
            sample_waveforms_list.append(sample_wf)
            ap_types.append(ap_type)
            trough_to_peak_durations.append(trough_to_peak_ms)

        nwbfile.units.add_column(
            name='depth',
            description='Depth of the unit from the cortical surface (um)',
            data=depths,
        )
        nwbfile.units.add_column(
            name='snr',
            description='Signal-to-noise ratio of the unit waveform',
            data=snrs,
        )
        nwbfile.units.add_column(
            name='waveform_mean',
            description='Mean waveform of the unit (microvolts)',
            data=mean_waveforms,
        )
        nwbfile.units.add_column(
            name='waveform_std',
            description='Standard deviation of the waveform across all spikes (microvolts)',
            data=std_waveforms,
        )
        nwbfile.units.add_column(
            name='ap_type',
            description='Action potential type classification: somatic (>0.4ms trough-to-peak) or axonal (<=0.4ms)',
            data=ap_types,
        )
        nwbfile.units.add_column(
            name='trough_to_peak_ms',
            description='Action potential duration from waveform trough to peak (milliseconds)',
            data=trough_to_peak_durations,
        )

        io.write(nwbfile)

    # ============== ADD TRIALS FROM NIDAQ DATA ==============
    print("\nProcessing NIDAQ data for trial extraction...")

    nidaq_sampling_rate = get_sampling_rate_from_meta(NIDAQ_meta_path)

    nidaq_num_channels = get_num_channels_from_meta(NIDAQ_meta_path)
    nidaq_data = load_nidaq_binary(NIDAQ_bin_path, num_channels=nidaq_num_channels)

    mat_data = load_experiment_mat_file(NIDAQ_bin_path)
    expt_info = mat_data['expt_info']

    # Get TTL signals
    fixation_ttl = nidaq_data[:, 4]     # Channel 4: fixation maintained
    stim_onset_ttl = nidaq_data[:, 5]   # Channel 5: stimulus onset

    laser_ch = 7
    laser_ttl = nidaq_data[:, laser_ch]

    stim_onsets = detect_ttl_onsets(stim_onset_ttl, threshold=15000)
    fixation_offsets = detect_ttl_offsets(fixation_ttl, threshold=15000)
    laser_onsets = detect_ttl_onsets(laser_ttl, threshold=15000)

    # Remove trials in which stimulus onset is less than 500 ms from the end
    # of the trial (fixstop)
    short_stim_inds = []
    stim_onsets_sec = stim_onsets / nidaq_sampling_rate
    fixation_offsets_sec = fixation_offsets / nidaq_sampling_rate
    laser_onsets_sec = laser_onsets / nidaq_sampling_rate
    for i, val in enumerate(stim_onsets_sec):
        fixstop_tmp = fixation_offsets >= val
        if np.min(np.abs(fixstop_tmp - val)) < expt_info.image_duration:
            short_stim_inds.append(i)

    # Check whether the laser was on for each trial
    laser_on = []
    for st_on in stim_onsets_sec:
        if np.any(np.abs(laser_onsets_sec - st_on) < 0.1):
            laser_on.append(True)
        else:
            laser_on.append(False)

    natural_images, ENVS = linear_stimuli(expt_info, stim_type='image')
    natural_images = natural_images[~np.isnan(natural_images)]

    trials = pd.DataFrame({
        'start_time': stim_onsets_sec,
        'stop_time': stim_onsets_sec + expt_info.image_duration,
        'natural_images': natural_images,
        'laser_on': laser_on,
    })

    # ---- Per-trial gaze traces (NIDAQ ch0 = x, ch1 = y) ----
    # Raw voltage samples spanning [stim_onset, stim_onset + image_duration],
    # sampled at nidaq_sampling_rate. Stored as ragged VectorData so trials
    # that hit the end of the recording (and are zero-padded shorter) still
    # round-trip cleanly.
    gaze_x_ch = 0
    gaze_y_ch = 1
    n_gaze_samples = int(round(float(expt_info.image_duration) * nidaq_sampling_rate))
    n_nidaq_samples = nidaq_data.shape[0]

    gaze_x_traces = []
    gaze_y_traces = []
    for onset_idx in stim_onsets:
        start = int(onset_idx)
        end = min(start + n_gaze_samples, n_nidaq_samples)
        gaze_x_traces.append(
            np.asarray(nidaq_data[start:end, gaze_x_ch], dtype=np.int16)
        )
        gaze_y_traces.append(
            np.asarray(nidaq_data[start:end, gaze_y_ch], dtype=np.int16)
        )

    gaze_description_suffix = (
        f' Raw NIDAQ int16 samples starting at stimulus onset, spanning '
        f'image_duration={float(expt_info.image_duration):.4f} s at '
        f'{nidaq_sampling_rate:.2f} Hz ({n_gaze_samples} samples/trial nominal).'
    )

    # ---- Per-trial pupil-diameter traces (NIDAQ ch2), blink-cleaned ----
    # Cleaning follows Kret & Sjak-Shie (2019, Behav. Res. Methods) speed-filter
    # + Mathot (2013) style pad-and-interpolate. See clean_pupil_blinks() docstring.
    pupil_ch = 2
    pupil_raw_traces = []
    pupil_clean_traces = []
    pupil_blink_masks = []
    for onset_idx in stim_onsets:
        start = int(onset_idx)
        end = min(start + n_gaze_samples, n_nidaq_samples)
        raw = np.asarray(nidaq_data[start:end, pupil_ch], dtype=np.float32)
        cleaned, blink_mask = clean_pupil_blinks(raw, nidaq_sampling_rate)
        pupil_raw_traces.append(np.asarray(nidaq_data[start:end, pupil_ch], dtype=np.int16))
        pupil_clean_traces.append(cleaned)
        pupil_blink_masks.append(blink_mask)

    pupil_description_suffix = (
        f' NIDAQ channel {pupil_ch}, starting at stimulus onset, spanning '
        f'image_duration={float(expt_info.image_duration):.4f} s at '
        f'{nidaq_sampling_rate:.2f} Hz ({n_gaze_samples} samples/trial nominal).'
    )

    with NWBHDF5IO(nwbfile_path, 'r+') as io:
        nwbfile = io.read()
        dataframe_to_nwb_trials(nwbfile, trials)
        nwbfile.add_trial_column(
            name='gaze_x',
            description='Horizontal gaze trace from NIDAQ channel 0.'
                        + gaze_description_suffix,
            data=gaze_x_traces,
            index=True,
        )
        nwbfile.add_trial_column(
            name='gaze_y',
            description='Vertical gaze trace from NIDAQ channel 1.'
                        + gaze_description_suffix,
            data=gaze_y_traces,
            index=True,
        )
        nwbfile.add_trial_column(
            name='pupil_diameter_raw',
            description='Raw pupil-diameter trace (int16 NIDAQ counts) BEFORE '
                        'blink cleaning, provided for QC.'
                        + pupil_description_suffix,
            data=pupil_raw_traces,
            index=True,
        )
        nwbfile.add_trial_column(
            name='pupil_diameter',
            description='Blink-cleaned pupil-diameter trace (float32). Blinks '
                        'detected via Kret & Sjak-Shie (2019) speed filter '
                        '(median + 16*MAD of |d(pupil)/dt|) plus signal-loss '
                        'floor; artifact runs padded by 50 ms before / 150 ms '
                        'after (Mathot 2013 / Hershman et al. 2018); short '
                        '(<20 ms) valid islands discarded; interior gaps '
                        'linearly interpolated; edge gaps left as NaN. See '
                        'pupil_blink_mask for which samples were '
                        'interpolated.'
                        + pupil_description_suffix,
            data=pupil_clean_traces,
            index=True,
        )
        nwbfile.add_trial_column(
            name='pupil_blink_mask',
            description='Boolean mask, True where the corresponding sample in '
                        'pupil_diameter was flagged as blink/artifact and '
                        'interpolated (or set to NaN at trial edges).'
                        + pupil_description_suffix,
            data=pupil_blink_masks,
            index=True,
        )
        io.write(nwbfile)

    # Store expt_info as a scratch entry (assigning to
    # nwbfile.experiment_description on an r+ file is silently dropped).
    expt_metadata = expt_info_to_metadata(expt_info)
    expt_metadata_str = json.dumps(expt_metadata, indent=2, default=str)

    with NWBHDF5IO(nwbfile_path, 'r+') as io:
        nwbfile = io.read()
        nwbfile.add_scratch(
            expt_metadata_str,
            name='expt_info_json',
            description='JSON-serialized experiment metadata from the session .mat '
                        'file (excludes trial_records, which are in the trials table).',
        )
        io.write(nwbfile)

    print("\nConversion complete!")
    print(f"NWB file saved to: {nwbfile_path}")

    # Quick verification
    from pynwb import read_nwb
    nwbfile = read_nwb(nwbfile_path)
    print(f"\nNWB file contains:")
    print(f"  - {len(nwbfile.units)} units")
    print(f"  - {len(nwbfile.trials)} trials")

    return Path(nwbfile_path)


if __name__ == "__main__":
    # Default single-session invocation (original hardcoded configuration).
    convert_session(
        monkey_name="MM005-Desta",
        date_str="2026-03-24",
        birthday_str="2021-04-23",
        sex="F",
        session_name="BSD",
        gate="g0",
        project="Inhibition-neural-information-processing",
    )