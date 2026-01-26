from pynwb import read_nwb

folder_path = "F:/localDATA/Electrophysiology/MM004-Wolfjaw/2025-12-05/ori_contrast_Lipshutz_g0/ori_contrast_Lipshutz_g0_imec0/kilosort4/"
nwbfile_path = folder_path+"myNWB.nwb"

nwbfile = read_nwb(nwbfile_path)

