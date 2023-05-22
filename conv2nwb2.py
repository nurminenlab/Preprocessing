from pynwb import NWBFile, get_manager
from datetime import datetime
from pynwb.form.backends.hdf5 import HDF5IO
from pynwb.ecephys import ElectricalSeries, TimeSeries, ElectrodeTable, ElectrodeTableRegion
# import psutil
import scipy.io as scio
import numpy as np
import datapreprocesslib as dpl
import sys
sys.path.append('/home/lauri/code/brPy/')
from brpylib import NsxFile

# name of data file
datafile = '20170518-223225-001.ns5'
pulsfile = '20170518-223225-001.ns4'
stimfile = 'stimParamsOR.mat'

# data lives here
bs_fldr  = '/opt3/MM385/P1/20170518-223225/'

datafile = bs_fldr + datafile
pulsfile = bs_fldr + pulsfile
stimfile = bs_fldr + stimfile

###################################
# init NWB-file
f = NWBFile('N', 'Anesthetized marmoset recording', 'MM666', datetime.now(),
            experimenter='Dr. Lauri Nurminen',
            lab='Angelucci-lab',
            institution='University of Utah',
            experiment_description='Orientation tuning recording using Utah array in monkey V1',
            session_id='MM666')

# create device object
device = f.create_device(name='128-ch Cerebus', source="a source")

# create electrode group
electrode_name = 'Utah-array'
source = "Blackrock Microsystems"
description = "96-channel Utah array"
location = "V1"
electrode_group = f.create_electrode_group(electrode_name,
                                           source=source,
                                           description=description,
                                           location=location,
                                           device=device)

###################################
# open file and extract headers
nsx_file = NsxFile(bs_fldr + datafile)
Data     = nsx_file.getdata()

# get continuous pulse timeseries
pls_file = NsxFile(pulsfile)
pls_Data = pls_file.getdata()
pls_Data = pls_Data['data'][0]

# extract data parameters
NS5_numsamples = Data['data_headers'][0]['NumDataPoints']
sample_rate = Data['samp_per_s']
channel_num =  nsx_file.basic_header['ChannelCount']
nsx_file.close()

# set electrode table for the file
electrode_table = ElectrodeTable('V-Probe_table')
# add electrode entries
for idx in np.arange(channel_num):
    electrode_table.add_row(idx,
                    x=1.0, y=2.0, z=0.0,
                    imp=float(-idx),
                    location='V1', filtering='0.3-7500Hz',
                    description='channel %s' % idx, group=electrode_group)


# set electrode table
f.set_electrode_table(electrode_table)
# include all electrode at once
electrode_table_region = ElectrodeTableRegion(electrode_table, list(range(channel_num)), 'All electrodes')
# write data for all electrodes
ephys_ts = ElectricalSeries('orientation1',
                            'an hypothetical source',
                            Data['data'],
                            electrode_table_region,
                            timestamps=np.arange(NS5_numsamples)/sample_rate,
                            resolution=1.0,
                            comments="This data was randomly generated with numpy, using 1234 as the seed",
                            description="Random numbers generated with numpy.random.rand")

# read in stimulus info 
stim_mat = scio.loadmat(stimfile)
stim_mat = dpl.getStimParams(stim_mat)
# read pulses
pulse_samples= dpl.PulseCounter(pls_Data,1000)

stimulus_data = dpl.modStim4NWB2(stim_mat['stim_value'],pulse_samples,NS5_numsamples,stim_mat['on_time_ms'],stim_mat['off_time_ms'])
# stimulus_ts = TimeSeries('Stimulus timeseries',
#                          'VSG',
#                          stimulus_data,
#                          'degrees',
#                          timestamps=ephys_timestamps,
#                          resolution=0.001,
#                          comments="This data was randomly generated with numpy, using 1234 as the seed",
#                          description="Imagined orientation timeseries")
                         



filename = 'orientation-data.h5'
f.add_acquisition(ephys_ts)
io = HDF5IO(filename, manager=get_manager(), mode='w')
io.write(f)
io.close()
#f.add_stimulus(stimulus_ts)

