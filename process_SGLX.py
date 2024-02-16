from lib import readSGLX
from utils import utils
import numpy as np
import subprocess 

CatGT_dir  = ' C:/Users/lonurmin/Desktop/code/CatGT-win'
TPrime_dir = ' C:/Users/lonurmin/Desktop/code/TPrime-win'

# catGT parameters
DIR = ' -dir=C:/localDATA/Electrophysiology/MM001-Sansa/2024-02-01'
RUN = ' -run=RF-mapping-wedge'
prs = ' -g=0 -t=0 -ap -ni -prb_fld -prb=0'
xa2 = ' -xa=0,0,2,2.5,1,0' # trial start
xia = ' -xia=0,0,2,2.5,3.5,0' # trial stop
xa  = ' -xa=0,0,4,2.5,1,25' # stim start

# TPrime parameters
tostream = utils.getFilePath(windowTitle="SYNC tostream file (IMEC0)")
fromstream = utils.getFilePath(windowTitle="SYNC fromstream file (usually NIDAQ)")
trialstart = utils.getFilePath(windowTitle="SYNC trial start (usually xa2)")
trialstop  = utils.getFilePath(windowTitle="SYNC trial stop (usually xia2)")
stimstart  = utils.getFilePath(windowTitle="SYNC stimulus (usually xa4)")

syncperiod = ' -syncperiod=1.0'
tostream = ' -tostream='+str(tostream)
fromstream = ' -fromstream=1,'+str(fromstream)
trialstart = ' -events=1,'+str(trialstart)+','+str(trialstart)[0:len(str(trialstart))-len(str(trialstart.stem)+'.txt')]+'trialstart.txt'
trialstop = ' -events=1,'+str(trialstop)+','+str(trialstop)[0:len(str(trialstop))-len(str(trialstop.stem)+'.txt')]+'trialstop.txt'
stimstart = ' -events=1,'+str(stimstart)+','+str(stimstart)[0:len(str(stimstart))-len(str(stimstart.stem)+'.txt')]+'stimstart.txt'

# run CatGT
subprocess.run('cd '+CatGT_dir, shell=True)
subprocess.run('CatGT'+DIR+RUN+prs+xia+xa+xa2, shell=True)

binFullPath = utils.getFilePath(windowTitle="Select binary file")
spikesFullPath = utils.getFilePath(windowTitle="Select spikes file")

# convert spike times to seconds
meta = readSGLX.readMeta(binFullPath)
sRate = readSGLX.SampRate(meta)
spike_times_smp = np.load(spikesFullPath)
spike_times_sec = np.around(np.divide(spike_times_smp,sRate,dtype=float),decimals=6)
np.save(spikesFullPath.with_stem('spike_times_sec'),spike_times_sec)

# run TPrime 
subprocess.run('cd '+TPrime_dir, shell=True)
subprocess.run('TPrime'+syncperiod+tostream+fromstream+trialstart+trialstop+stimstart, shell=True)