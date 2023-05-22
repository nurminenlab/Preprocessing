import numpy as np
import tables as tb
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
import sys
import pandas as pd
sys.path.append('/home/lauri/code/Libs/')
import spikelib as spk

# information about penetration is stored in a csv file
notes_dir = '/home/lauri/projects/CorrelatedVariability/notes/'
notes_fle = 'penetrationinfo.csv'
df = pd.read_csv(notes_dir+notes_fle)

animal = df['ANIMAL'].values
penetr = df['PENETRATION'].values
fldrs  = df['FOLDER'].values

# for debugging
# animal = animal[0:2]
# penetr = penetr[0:2]
# fldrs  = fldrs[0:2]

results_root = '/home/lauri/projects/CorrelatedVariability/results/pass_1/'

for p,val in enumerate(animal):
    file_name    = results_root+'qualityPlots/'+animal[p]+penetr[p]+fldrs[p]+'.pdf'
    pdf          = PdfPages(file_name)

    print('Now processing '+animal[p]+penetr[p]+fldrs[p])
    
    # load pre-processed Kilosort outputs
    spkC_L   = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_spkCnts_L.npy')
    spkC_NoL = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_spkCnts_NoL.npy')
    baseLine = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_baseLine.npy')
    spkR_L   = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_spkraster_L.npy')
    spkR_NoL = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_spkraster_NoL.npy')
    RF_L     = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_fieldParams_L.npy')
    RF_NoL   = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_fieldParams_NoL.npy')
    SNR      = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_SNR.npy')
    depth    = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_spikeDepths.npy')
    spkDur   = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_spikeWidths.npy')
    diams    = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'diams.npy')[0]

    # load manual unit classification
    df = pd.read_csv('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'cluster_groups.csv',sep='\t')
    cluster_group = df['group'].values

    # loop units and plot
    gs = GridSpec(5,2)
    for i in range(spkC_NoL.shape[0]):

        # don't plot noise clusters 
        if cluster_group[i] == 'noise' or np.max(np.mean(spkC_NoL[i,:,:],axis=1)) < 5 or SNR[i] < 2.0:
            print('Noise cluster encountered')
        else:
            plt.subplot(gs[0:2,0])
            plt.errorbar(diams, np.mean(spkC_NoL[i,:,:],axis=1), yerr=np.std(spkC_NoL[i,:,:],axis=1) / np.sqrt(spkC_NoL.shape[2]), fmt='ko-')
            plt.errorbar(diams, np.mean(spkC_L[i,:,:],axis=1), yerr=np.std(spkC_L[i,:,:],axis=1) / np.sqrt(spkC_NoL.shape[2]), fmt='go-')
            plt.title('unit '+str(i))
    
            a = -1
            for cond in range(spkR_NoL.shape[1]-1):
                a = a + 1
                max_psth = 0
                if a > 2:
                    a = 0

                psth   = spk.PSTH(spkR_NoL[i,cond,:,:],bin_width=40)
                psth_L = spk.PSTH(spkR_L[i,cond,:,:],bin_width=40)

                if cond <=2:
                    plt.subplot(gs[2+a,0])
                    psth   = np.mean(psth,axis=0) / 0.001
                    psth_L = np.mean(psth_L,axis=0) / 0.001
                    plt.plot(psth,'k-')
                    plt.plot(psth_L,'g-')
                    plt.plot([100, 100,], [0,np.max(psth)], 'k-')
                    plt.plot([100, 100,], [0,np.max(psth_L)], 'k-')
                    plt.axis('off')
                    if np.max([np.max(psth), np.max(psth_L)]) > max_psth:
                        max_psth = np.max([np.max(psth), np.max(psth_L)])
                
                else:
                    plt.subplot(gs[2+a,1])
                    psth   = np.mean(psth,axis=0) / 0.001
                    psth_L = np.mean(psth_L,axis=0) / 0.001
                    plt.plot(psth,'k-')
                    plt.plot(psth_L,'g-')
                    plt.plot([100, 100,], [0,np.max(psth)], 'k-')
                    plt.plot([100, 100,], [0,np.max(psth_L)], 'k-')
                    plt.axis('off')
                    if np.max([np.max(psth), np.max(psth_L)]) > max_psth:
                        max_psth = np.max([np.max(psth), np.max(psth_L)])

            # set equal y-range
            a = -1
            for cond in range(spkR_NoL.shape[1]-1):
                a = a + 1
                if a > 2:
                    a = 0

            if cond <=2:
                plt.subplot(gs[2+a,0])
                plt.ylim(0,max_psth)
            else:
                plt.subplot(gs[2+a,1])
                plt.subplot(gs[2+a,0])
                plt.ylim(0,max_psth)


            pdf.savefig()
            plt.clf()
    
    pdf.close()
