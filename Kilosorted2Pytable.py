import numpy as np
import tables as tb
import pandas as pd

# information about penetration is stored in a csv file
notes_dir = '/home/lauri/projects/CorrelatedVariability/notes/'
notes_fle = 'penetrationinfo.csv'
df = pd.read_csv(notes_dir+notes_fle)

# we want to drop bad runs, as judged qualitatively, using voodoo you know
df = df[df.QUALITY != 'bad']
# and runs with visible drift
df = df[df.OTHER != 'drift']

animal = df['ANIMAL'].values
penetr = df['PENETRATION'].values
fldrs  = df['FOLDER'].values

results_root = '/home/lauri/projects/CorrelatedVariability/results/'
file_name    = results_root+'collated_correlation_data.h5'

# open file for read/write
data_file  = tb.open_file(file_name,'a')
# create group
data_group = data_file.create_group('/','data_group',"Data tables for each penetration") 

# loop through penetrations and animals
# we need dummy variable for skipping folders that are read in the loop
a = -1
while a < animal.shape[0]-1:
    a += 1
    cluster_group = pd.read_csv('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'cluster_groups.csv',sep='\t')
    cluster_id = cluster_group['cluster_id'].values
    # some runs were separated to different folders (separuns = yes)
    if df.SEPARUNS.values[a] == 'yes':
        
        # which unit numbers correspond across folders is stored in a dictionary file
        with open('/home/lauri/projects/CorrelatedVariability/notes/unit_corresp_'+df.ANIMAL.values[a]+df.PENETRATION.values[a]+'_c'+str(df.CONTRAST.values[a])+'.py') as corresp_dict_file:
            corresp_dict = eval(corresp_dict_file.read())

        # 
        unit_numbers = [corresp_dict[i] for i in corresp_dict]
        unit_numbers = np.array(unit_numbers)
        
        for i in range(unit_numbers.shape[1]):
            #
            if i == 0:
                # no need to load cluster id here as it was done before entering the if clause
                corresp_indx = [np.where(cluster_id == unit_numbers[j,i])[0][0] for j in range(unit_numbers[:,i].shape[0])]
                spkC_L   = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spkCnts_L.npy')[corresp_indx,:,:]
                spkC_NoL = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spkCnts_NoL.npy')[corresp_indx,:,:]
                baseLine = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_baseLine.npy')[corresp_indx,:]
                spkR_L   = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spkraster_L.npy')[corresp_indx,:,:,:]
                spkR_NoL = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spkraster_NoL.npy')[corresp_indx,:,:,:]
            else:
                a +=1
                #  cluster id needs to be updated for the new unit assignment
                cluster_group = pd.read_csv('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'cluster_groups.csv',sep='\t')
                cluster_id = cluster_group['cluster_id'].values
                corresp_indx = [np.where(cluster_id == unit_numbers[j,i])[0][0] for j in range(unit_numbers[:,i].shape[0])]
                    
                spkC_L   = np.append(spkC_L,   np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spkCnts_L.npy')[corresp_indx,:,:],axis=2)
                spkC_NoL = np.append(spkC_NoL, np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spkCnts_NoL.npy')[corresp_indx,:,:],axis=2)
                baseLine = np.append(baseLine, np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_baseLine.npy')[corresp_indx,:],axis=1)
                spkR_L   = np.append(spkR_L,   np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spkraster_L.npy')[corresp_indx,:,:,:], axis=2)
                spkR_NoL = np.append(spkR_NoL, np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spkraster_NoL.npy')[corresp_indx,:,:,:], axis=2)
                RF_L     = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_fieldParams_L.npy')
                RF_NoL   = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_fieldParams_NoL.npy')
                SNR      = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_SNR.npy')
                depth    = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spikeDepths.npy')
                spkDur   = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spikeWidths.npy')
                # gather stimulus info
                diams              = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'diams.npy')
                contrast           = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'contrast.npy')
                spatial_frequency  = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'TF.npy')
                temporal_frequency = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'SF.npy')
                
    else:
        # this part in a loop if condition met, that's it
        spkC_L   = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spkCnts_L.npy')
        spkC_NoL = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spkCnts_NoL.npy')
        baseLine = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_baseLine.npy')
        spkR_L   = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spkraster_L.npy')
        spkR_NoL = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spkraster_NoL.npy')
        RF_L     = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_fieldParams_L.npy')
        RF_NoL   = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_fieldParams_NoL.npy')
        SNR      = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_SNR.npy')
        depth    = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spikeDepths.npy')
        spkDur   = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'Kilosorted_spikeWidths.npy')
        # gather stimulus info
        diams              = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'diams.npy')
        contrast           = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'contrast.npy')
        spatial_frequency  = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'TF.npy')
        temporal_frequency = np.load('/opt3/'+df.ANIMAL.values[a]+'/'+df.PENETRATION.values[a]+'/'+df.FOLDER.values[a]+'/'+'SF.npy')
    
    
    # prepare table
    # needs to be done in the loop to accommodate different
    # array sizes in different penetrations
    dataTable = {'spkC_L':tb.Float64Col(shape=(spkC_L.shape[1],spkC_L.shape[2])),
                 'spkC_NoL':tb.Float64Col(shape=(spkC_NoL.shape[1],spkC_NoL.shape[2])),
                 'baseLine':tb.Float64Col(shape=(baseLine.shape[1])),
                 'spkR_L':tb.Float64Col(shape=(spkR_L.shape[1],spkR_L.shape[2],spkR_L.shape[3])),
                 'spkR_NoL':tb.Float64Col(shape=(spkR_NoL.shape[1],spkR_NoL.shape[2],spkR_NoL.shape[3])),
                 'RF_L':tb.Float64Col(shape=(RF_L.shape[1])),
                 'RF_NoL':tb.Float64Col(shape=(RF_NoL.shape[1])),
                 'SNR':tb.Float64Col(1),
                 'depth':tb.Float64Col(1),
                 'layer':tb.StringCol(3),
                 'spkDur':tb.Float64Col(1),
                 'diams':tb.Float64Col(shape=(diams.shape[1])),
                 'contrast':tb.Float64Col(shape=(contrast.shape[1])),
                 'spatial_frequency':tb.Float64Col(shape=(spatial_frequency.shape[1])),
                 'temporal_frequency':tb.Float64Col(shape=(temporal_frequency.shape[1]))}
    
    # create data table
    # replace fldrs with contrast
    table = data_file.create_table(data_group, df.ANIMAL.values[a]+df.PENETRATION.values[a]+str(df.CONTRAST.values[a]), dataTable, 'Preprocessed spike-sorted data')

    # write unit data to table
    unit = table.row
    for u in range(spkC_L.shape[0]):
        unit['spkC_L']   = spkC_L[u,:,:]
        unit['spkC_NoL'] = spkC_NoL[u,:,:]
        unit['baseLine'] = baseLine[u]
        unit['spkR_L']   = spkR_L[u,:,:,:]
        unit['spkR_NoL'] = spkR_NoL[u,:,:,:]
        unit['RF_L']     = RF_L[u,:]
        unit['RF_NoL']   = RF_NoL[u,:]
        unit['SNR']      = SNR[u]
        unit['depth']    = depth[u]
        unit['spkDur']   = spkDur[u]
        unit['diams']    = diams[0,:]
        unit['contrast'] = contrast[0,:]
        unit['spatial_frequency']  = spatial_frequency[0,:]
        unit['temporal_frequency'] = temporal_frequency[0,:]
        # append row
        unit.append()

    # flush data
    table.flush()

# close file
data_file.close()
