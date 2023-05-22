import sys
import numpy as np
sys.path.append('/opt/mountainlab/packages/pyms/mlpy/')
sys.path.append('/home/lauri/code/brPy/')
import brpylib as brlib 
import mdaio as mdaio

# name of data file
datafile = '20160329-231153-001'
stimfile = 'stimParamsOR.mat'

# data lives here
bs_fldr  = '/opt3/MM378/mountainsort_test/'

stimfile = bs_fldr + stimfile

nsx_file = brlib.NsxFile(bs_fldr + datafile + '.ns5')
Data     = nsx_file.getdata()
nsx_file.datafile.close()

# convert back to uint16
digital_factor = brlib.getdigfactor(nsx_file.extended_headers, 0)
Data['data'] /= brlib.getdigfactor(nsx_file.extended_headers, 0)
Data['data']  = Data['data'].astype(np.uint16)

mdaio.writemda16ui(Data['data'][0:24,:], bs_fldr + datafile + '.mda')
