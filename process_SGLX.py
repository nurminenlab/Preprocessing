from lib import readSGLX
from utils import utils
import numpy as np
from pathlib import Path
from tkinter import Tk
from tkinter import filedialog
import subprocess 

#https://billkarsh.github.io/SpikeGLX/help/syncEdges/Sync_edges/

# run catGT to extract sync edges
subprocess.run("CatGT -dir=C:\localDATA\Electrophysiology\MM001-Sansa\2024-02-01\ -run=RF-mapping-wedge -g=1 -t=0 -ap -ni -prb_fld -prb=0", shell=True)

binFullPath = utils.getFilePath(windowTitle="Select binary file")
spikesFullPath = utils.getFilePath(windowTitle="Select spikes file")

meta = readSGLX.readMeta(binFullPath)
sRate = readSGLX.SampRate(meta)
spike_times_smp = np.load(spikesFullPath)
spike_times_sec = np.around(np.divide(spike_times_smp,sRate,dtype=float),decimals=6)
np.save(spikesFullPath.with_stem('spike_times_sec'),spike_times_sec)

# run TPrime align spikes to trials
subprocess.run("TPrime -syncperiod=1.0 -tostream=C:\localDATA\Electrophysiology\MM001-Sansa\2024-02-01\RF-mapping-wedge_g1\RF-mapping-wedge_g1_imec0\RF-mapping-wedge_g1_tcat.imec0.ap.xd_384_6_500.txt
 -fromstream=1,D:/CGT_OUT/catgt_demo_g0/demo_g0_imec1/demo_g0_tcat.imec1.ap.xd_384_6_500.txt ^
 -fromstream=2,D:/CGT_OUT/catgt_demo_g0/demo_g0_tcat.nidq.xd_1_3_500.txt ^
 -events=1,D:/CGT_OUT/catgt_demo_g0/demo_g0_imec1/spike_seconds.npy,D:/CGT_OUT/catgt_demo_g0/demo_g0_imec1/spike_seconds_adj.npy ^
 -events=2,D:/CGT_OUT/catgt_demo_g0/demo_g0_tcat.nidq.xa_0_25.txt,D:/CGT_OUT/catgt_demo_g0/go_cue.txt ^
 -events=2,D:/CGT_OUT/catgt_demo_g0/demo_g0_tcat.nidq.xd_1_2_0.txt,D:/CGT_OUT/catgt_demo_g0/nose_poke.txt
