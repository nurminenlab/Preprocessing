import numpy as np
import tables as tb
import pandas as pd
import sys
sys.path.append('/home/lauri/code/Preprocessing/')
import datapreprocesslib as prepro

results_root = '/home/lauri/projects/CorrelatedVariability/results/'
file_name    = results_root+'collated_correlation_data_combinedsorting_400ms_MK-July2020.h5'

# information about penetration is stored in a csv file
notes_dir = '/home/lauri/projects/CorrelatedVariability/notes/'
notes_fle = 'penetrationinfo_macaque.csv'
df = pd.read_csv(notes_dir+notes_fle)

# we want to drop bad runs, as judged qualitatively, using voodoo you know
df = df[df.QUALITY != 'bad']
# and runs with visible drift
df = df[df.OTHER != 'drift']

animal = df['ANIMAL'].values
penetr = df['PENETRATION'].values
fldrs  = df['FOLDER'].values

# open file for read/write
data_file  = tb.open_file(file_name,'a')
# create group
data_group = data_file.create_group('/','data_group',"Data tables for each penetration") 

# loop through penetrations and animals
# we need dummy variable for skipping folders that are read in the loop
a = -1
while a < animal.shape[0]-1:
    a += 1
    # cluster_group = pd.read_csv('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'cluster_groups.csv',sep='\t')
    # cluster_id = cluster_group['cluster_id'].values

    # this part in a loop if condition met, that's it
    spkC_NoL = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spkCnts_NoL_400ms.npy')
    baseLine = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_baseLine_400ms.npy')
    spkR_NoL = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spkraster_NoL_400ms.npy')
    RF_NoL   = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_fieldParams_NoL_400ms.npy')
    SNR      = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_SNR_400ms.npy')
    depth    = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spikeDepths_400ms.npy')[:,0]
    spkDur   = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spikeWidths_400ms.npy')
    spkWav   = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_waveForms_400ms.npy')
    maxAmpC  = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_maxamp_contact_400ms.npy')
    
    # gather stimulus info
    diams              = np.squeeze(np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_diams_400ms.npy'))
    contrast           = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'contrast_400ms.npy')
    spatial_frequency  = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'TF_400ms.npy')
    temporal_frequency = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'SF_400ms.npy')
    
    # prepare table
    # needs to be done in the loop to accommodate different
    # array sizes in different penetrations
    dataTable = {'spkC_NoL':tb.Float64Col(shape=(spkC_NoL.shape[1],spkC_NoL.shape[2])),
                 'baseLine':tb.Float64Col(shape=(baseLine.shape[1],baseLine.shape[2])),
                 'spkR_NoL':tb.Float64Col(shape=(spkR_NoL.shape[1],spkR_NoL.shape[2],spkR_NoL.shape[3])),
                 'RF_NoL':tb.Float64Col(shape=(RF_NoL.shape[1])),
                 'SNR':tb.Float64Col(1),
                 'waveForms':tb.Float64Col(shape=(spkWav.shape[1],spkWav.shape[2])),
                 'maxAmpC':tb.Float64Col(1),
                 'depth':tb.Float64Col(1),
                 'layer':tb.StringCol(3),
                 'animal':tb.StringCol(5),
                 'penetration':tb.StringCol(2),
                 'spkDur':tb.Float64Col(1),
                 'diams':tb.Float64Col(shape=(diams.shape[0])),
                 'contrast':tb.Float64Col(shape=contrast.shape),
                 'spatial_frequency':tb.Float64Col(shape=(spatial_frequency.shape[1])),
                 'temporal_frequency':tb.Float64Col(shape=(temporal_frequency.shape[1]))}
    
    # create data table
    # replace fldrs with contrast
    table = data_file.create_table(data_group, df.ANIMAL.values[a]+df.PENETRATION.values[a], dataTable, 'Preprocessed spike-sorted data')

    layer_borders = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+'layer_borders.npy')
    
    # write unit data to table
    unit = table.row
    for u in range(spkC_NoL.shape[0]):
        unit['spkC_NoL'] = spkC_NoL[u,:,:]
        unit['baseLine'] = baseLine[u,:,:]
        unit['spkR_NoL'] = spkR_NoL[u,:,:,:]
        unit['RF_NoL']   = RF_NoL[u,:]
        unit['SNR']      = SNR[u]
        unit['depth']    = depth[u]
        unit['spkDur']   = spkDur[u]
        unit['diams']    = diams
        unit['contrast'] = contrast
        unit['spatial_frequency']  = spatial_frequency[0,:]
        unit['temporal_frequency'] = temporal_frequency[0,:]
        unit['waveForms'] = spkWav[u,:,:]
        unit['maxAmpC']   = maxAmpC[u,0]
        unit['layer']     = prepro.return_layer(layer_borders,depth[u])
        unit['animal']      = df.ANIMAL.values[a]
        unit['penetration'] = df.PENETRATION.values[a]
        
        # append row
        unit.append()

    # flush data
    table.flush()

# close file
data_file.close()
