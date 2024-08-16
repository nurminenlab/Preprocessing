from lib import readSGLX
from utils import utils
import numpy as np
import subprocess 
import pandas as pd
from time import sleep

CatGT_dir  = ' C:/Users/lonurmin/Desktop/code/CatGT-win'
TPrime_dir = ' C:/Users/lonurmin/Desktop/code/TPrime-win'

# catGT parameters
DIR = ' -dir=C:\\localDATA\\Electrophysiology\\MM001-Sansa\\2024-07-03'
RUN = ' -run=BSD'
prs = ' -g=1 -t=0 -ap -ni -prb_fld -prb=0'
trialstart = ' -xa=0,0,3,2.5,1,0' # trial start
trialstop = ' -xia=0,0,3,2.5,3.5,0' # trial stop
fixstart = ' -xa=0,0,4,2.5,1,0' # fix start
fixstop = ' -xia=0,0,4,2.5,3.5,0' # fix stop
stimstart  = ' -xa=0,0,5,2.5,1,0' # stim start
stimstop  = ' -xia=0,0,5,2.5,3.5,0' # stim start
laserstart  = ' -xa=0,0,7,2.5,1,0' # laser start
laserstop  = ' -xia=0,0,7,2.5,3.5,0' # laser start

# run CatGT
print('CatGT is running, please wait for the process to finish')
#subprocess.run('cd '+CatGT_dir, shell=True)
subprocess.run('CatGT'+DIR+RUN+prs+trialstart+trialstop+fixstart+fixstop+stimstart+stimstop+laserstart+laserstop, shell=True)

# convert spike times to seconds
binFullPath = utils.getFilePath(windowTitle="Select binary ap file",filetypes=[("sGLX binary","*.bin")])
spikesFullPath = utils.getFilePath(windowTitle="Select spike times in samples",filetypes=[("KS output spikes_times","*.npy")])
meta = readSGLX.readMeta(binFullPath)
sRate = readSGLX.SampRate(meta)
spike_times_smp = np.load(spikesFullPath)
spike_times_sec = np.around(np.divide(spike_times_smp,sRate,dtype=float),decimals=6)
np.save(spikesFullPath.with_stem('spike_times_sec'),spike_times_sec)

# TPrime parameters
tostream = utils.getFilePath(windowTitle="SYNC tostream file (IMEC0 edgefile.txt)",filetypes=[("CatGT output","*xd_384_6_500.txt")])
fromstream = utils.getFilePath(windowTitle="SYNC fromstream file (usually NIDAQ edgefile.txt)",filetypes=[("CatGT output","*xa_6_500.txt")])
trialstart = utils.getFilePath(windowTitle="SYNC trial start (usually xa3)",filetypes=[("CatGT output","*xa_3_0.txt")])
trialstop  = utils.getFilePath(windowTitle="SYNC trial stop (usually xia3)",filetypes=[("CatGT output","*xia_3_0.txt")])
fixstart = utils.getFilePath(windowTitle="SYNC fixation start (usually xa4)",filetypes=[("CatGT output","*xa_4_0.txt")])
fixstop  = utils.getFilePath(windowTitle="SYNC fixation stop (usually xia4)",filetypes=[("CatGT output","*xia_4_0.txt")])
stimstart  = utils.getFilePath(windowTitle="SYNC stimulus (usually xa5)",filetypes=[("CatGT output","*xa_5_0.txt")])
stimstop   = utils.getFilePath(windowTitle="SYNC stimulus (usually xia5)",filetypes=[("CatGT output","*xia_5_0.txt")])
laserstart  = utils.getFilePath(windowTitle="SYNC laser (usually xa7)",filetypes=[("CatGT output","*xa_7_0.txt")])
laserstop  = utils.getFilePath(windowTitle="SYNC laser (usually xia7)",filetypes=[("CatGT output","*xia_7_0.txt")])

syncperiod = ' -syncperiod=1.0'
tostream = ' -tostream='+str(tostream)
fromstream = ' -fromstream=1,'+str(fromstream)
trialstart = ' -events=1,'+str(trialstart)+','+str(trialstart)[0:len(str(trialstart))-len(str(trialstart.stem)+'.txt')]+'trialstart.txt'
trialstop  = ' -events=1,'+str(trialstop)+','+str(trialstop)[0:len(str(trialstop))-len(str(trialstop.stem)+'.txt')]+'trialstop.txt'
fixstart   = ' -events=1,'+str(fixstart)+','+str(fixstart)[0:len(str(fixstart))-len(str(fixstart.stem)+'.txt')]+'fixstart.txt'
fixstop    = ' -events=1,'+str(fixstop)+','+str(fixstop)[0:len(str(fixstop))-len(str(fixstop.stem)+'.txt')]+'fixstop.txt'
stimstart  = ' -events=1,'+str(stimstart)+','+str(stimstart)[0:len(str(stimstart))-len(str(stimstart.stem)+'.txt')]+'stimstart.txt'
stimstop   = ' -events=1,'+str(stimstop)+','+str(stimstop)[0:len(str(stimstop))-len(str(stimstop.stem)+'.txt')]+'stimstop.txt'
laserstart = ' -events=1,'+str(laserstart)+','+str(laserstart)[0:len(str(laserstart))-len(str(laserstart.stem)+'.txt')]+'laserstart.txt'
laserstop  = ' -events=1,'+str(laserstop)+','+str(laserstop)[0:len(str(laserstop))-len(str(laserstop.stem)+'.txt')]+'laserstop.txt'

# run TPrime 
subprocess.run('cd '+TPrime_dir, shell=True)
subprocess.run('TPrime'+syncperiod+tostream+fromstream+trialstart+trialstop+fixstart+fixstop+stimstart+stimstop+laserstart+laserstop, shell=True) 

# get paths to the pulse files 
trialstartFullPath = utils.getFilePath(windowTitle="Select trialstart file",filetypes=[("TPrime output","trialstart.txt")])
trialstopFullPath  = utils.getFilePath(windowTitle="Select trialstop file",filetypes=[("TPrime output","trialstop.txt")])
stimstartFullPath  = utils.getFilePath(windowTitle="Select stimstart file",filetypes=[("TPrime output","stimstart.txt")])
stimstopFullPath   = utils.getFilePath(windowTitle="Select stimstop file",filetypes=[("TPrime output","stimstop.txt")])
fixstartFullPath   = utils.getFilePath(windowTitle="Select fixstart file",filetypes=[("TPrime output","fixstart.txt")])
fixstopFullPath    = utils.getFilePath(windowTitle="Select fixstop file",filetypes=[("TPrime output","fixstop.txt")])
laserstartFullPath  = utils.getFilePath(windowTitle="Select laserstart file",filetypes=[("TPrime output","laserstart.txt")])
laserstopFullPath   = utils.getFilePath(windowTitle="Select laserstop file",filetypes=[("TPrime output","laserstop.txt")])

# read the pulse files and convert to dataframe
## trials
trialstartDF = pd.read_csv(trialstartFullPath.absolute(),sep=" ",header=None)
trialstartDF.columns = ['trialstart']
trialstopDF = pd.read_csv(trialstopFullPath.absolute(),sep=" ",header=None)
trialstopDF.columns = ['trialstop']
trialsDF = trialstartDF.join(trialstopDF)

## stimuli
stimstartDF = pd.read_csv(stimstartFullPath.absolute(),sep=" ",header=None)
stimstartDF.columns = ['stimstart']
stimstopDF = pd.read_csv(stimstopFullPath.absolute(),sep=" ",header=None)
stimstopDF.columns = ['stimstop']
stimsDF = stimstartDF.join(stimstopDF)

# fixation
fixstartDF = pd.read_csv(fixstartFullPath.absolute(),sep=" ",header=None)
fixstartDF.columns = ['fixstart']
fixstopDF = pd.read_csv(fixstopFullPath.absolute(),sep=" ",header=None)
fixstopDF.columns = ['fixstop']
fixationDF = fixstartDF.join(fixstopDF)

## lasers
laserstartDF = pd.read_csv(laserstartFullPath.absolute(),sep=" ",header=None)
laserstartDF.columns = ['laserstart']
laserstopDF = pd.read_csv(laserstopFullPath.absolute(),sep=" ",header=None)
laserstopDF.columns = ['laserstop']
laserDF = laserstartDF.join(laserstopDF)

# save to dataframes
trialsDF.to_csv(trialstartFullPath.with_suffix('.csv'),index=False)
stimsDF.to_csv(stimstartFullPath.with_suffix('.csv'),index=False)
fixationDF.to_csv(fixstartFullPath.with_suffix('.csv'),index=False)
laserDF.to_csv(laserstartFullPath.with_suffix('.csv'),index=False)