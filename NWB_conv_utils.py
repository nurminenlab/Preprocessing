from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path
import numpy as np
from pynwb import NWBHDF5IO, NWBFile
from pynwb.misc import Units
from neuroconv.datainterfaces import KiloSortSortingInterface
import spikeinterface.extractors as se
from spikeinterface import read_binary
import pandas as pd
import scipy.io as sio
import glob

"""Utility functions for converting Kilorted spikes and our experiment description files to NWB."""



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
        if field_name == 'trial_records':
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


