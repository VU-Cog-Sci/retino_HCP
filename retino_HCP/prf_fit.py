"""
-----------------------------------------------------------------------------------------
prf_fit
-----------------------------------------------------------------------------------------
Goal of the script:
Create pRF estimates
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject name
sys.argv[2]: start voxel index
sys.argv[3]: end voxel index
sys.argv[4]: data file path
sys.argv[5]: main directory
-----------------------------------------------------------------------------------------
Output(s):
Nifti image files with fitting parameters per voxels
-----------------------------------------------------------------------------------------
"""

# General imports
from __future__ import division
import sys
import ctypes
import multiprocessing
import numpy as np
import tables
import platform
from math import *
import os
import glob
import json
import ipdb
import warnings
warnings.filterwarnings('ignore')

# MRI analysis imports
import nibabel as nb
import popeye.utilities as utils
from popeye.visual_stimulus import VisualStimulus
import popeye.css as css

# Get inputs
subject = sys.argv[1]
start_idx = sys.argv[2]
end_idx = sys.argv[3]
data_file = sys.argv[4]
base_dir = sys.argv[5]

# Define analysis parameters
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define cluster/server specific parameters
if 'lisa' in platform.uname()[1]:
    N_PROCS = 4
elif 'aeneas' in platform.uname()[1]:
    N_PROCS = 31
elif 'local' in platform.uname()[1]:
    N_PROCS = 8
else:
    N_PROCS = 25

# Define output file path and directories
base_file_name = os.path.split(data_file)[-1][:-7]
opfn = os.path.join(base_dir,'pp',subject,'prf',base_file_name + '_est_%s_to_%s.gii' %(start_idx,end_idx))
print('output file will be ' + opfn) 

try:
    os.makedirs(os.path.join(base_dir,'pp',subject,'prf'))
except:
    pass

# Load data
data = []
data_file_load = nb.load(data_file)
data.append(np.array([data_file_load.darrays[i].data for i in range(len(data_file_load.darrays))]))
data = np.vstack(data)
data_to_analyse = data[:,int(start_idx):int(end_idx)]

# Create stimulus design
visual_dm = []
file = tables.open_file(os.path.join(base_dir, 'retinotopysmall5.mat'))
visual_dm.append(file.get_node('/stim')[:])
file.close()
visual_dm = np.vstack(visual_dm).transpose((1,2,0))

stimulus = VisualStimulus(  stim_arr=visual_dm,
                            viewing_distance=analysis_info["screen_distance"],
                            screen_width=analysis_info["screen_width"],
                            scale_factor=1/5.0,
                            tr_length=analysis_info["TR"],
                            dtype=np.short)

# Initialize css model
css_model = css.CompressiveSpatialSummationModel(   stimulus = stimulus,
                                                    hrf_model = utils.spm_hrf)
css_model.hrf_delay = 0
print('models and stimulus loaded')

# Fit: define search grids
x_grid = (-12, 12)
y_grid = (-12, 12)
sigma_grid = (0.05, 10)
n_grid =  (0.01, 1)
css_grids =  (x_grid, y_grid, sigma_grid, n_grid)

# Fit: define search bounds
x_bound = (-30.0, 30.0)
y_bound = (-30.0, 30.0)
sigma_bound = (0.001, 70.0)
n_bound = (0.001, 1.5)
beta_bound = (-1e3, 1e3)
baseline_bound = (-1e3, 1e3)
css_bounds =   (x_bound, y_bound, sigma_bound, n_bound, beta_bound, baseline_bound)    # boundaries of the final gradient-descent

# Fit: define empty estimate and voxel indeces
estimates = np.zeros((7,data.shape[1]))
voxel_indices = [(xx, 0, 0) for xx in np.arange(int(start_idx),int(end_idx),1)]

# Define multiprocess bundle
bundle = utils.multiprocess_bundle( Fit = css.CompressiveSpatialSummationFit, 
                                    model = css_model, 
                                    data = data_to_analyse.T,
                                    grids = css_grids, 
                                    bounds = css_bounds, 
                                    indices = voxel_indices, 
                                    auto_fit = True, 
                                    verbose = 1, 
                                    Ns = 12)
# Run fitting
pool = multiprocessing.Pool(processes = N_PROCS)
output = pool.map(  func = utils.parallel_fit, 
                    iterable = bundle)


for fit in output:
    estimates[:6,fit.voxel_index[0]] = fit.estimate
    estimates[6,fit.voxel_index[0]] = fit.rsquared

# Free up memory
pool.close()
pool.join()

# Save data
darrays = [nb.gifti.gifti.GiftiDataArray(d) for d in estimates]

gii_out = nb.gifti.gifti.GiftiImage(header = data_file_load.header, 
                                    extra = data_file_load.extra,
                                    darrays = darrays)

nb.save(gii_out, opfn)
