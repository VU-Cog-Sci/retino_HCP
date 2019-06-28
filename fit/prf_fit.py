"""
-----------------------------------------------------------------------------------------
prf_fit.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create pRF estimates
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject name
sys.argv[2]: start voxel index -1
sys.argv[3]: end voxel index   -2
sys.argv[4]: data file path
sys.argv[5]: main directory   
-----------------------------------------------------------------------------------------
Output(s):
Gifti image files with fitting parameters per vertex
-----------------------------------------------------------------------------------------
"""

# General imports
from __future__ import division
import sys
import ctypes
import multiprocessing
import numpy as np
import scipy.io
import platform
from math import *
import os
import glob
import json
import ipdb
deb = ipdb.set_trace
opj = os.path.join
import warnings
warnings.filterwarnings('ignore')

# MRI analysis imports
import cifti
import popeye.utilities as utils
from popeye.visual_stimulus import VisualStimulus
import popeye.css as css
import popeye.og as og

# Get inputs
fit_model = sys.argv[1]
subject = sys.argv[2]
start_idx = sys.argv[3]
end_idx =  sys.argv[4]
data_file = sys.argv[5]
base_dir = sys.argv[6]

# Define analysis parameters
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define cluster/server specific parameters
if 'lisa' in platform.uname()[1]:
    N_PROCS = 16
elif 'aeneas' in platform.uname()[1]:
    N_PROCS = 2
else: # cartesius
    N_PROCS = 16
Ns = analysis_info["fit_step"]

# Define output file path and directories
base_file_name = os.path.split(data_file)[-1][:-13]
opfn_est = opj(base_dir,'pp_data',subject,fit_model,'fit', "{base_file_name}_est_{start}_to_{end}.dtseries.nii".format(  base_file_name = base_file_name,
                                                                                                                start = int(start_idx),
                                                                                                                end = int(end_idx)))

try: os.makedirs(opj(base_dir,'pp_data',subject,fit_model,'fit'))
except: pass

# Load data
data = []
data_file_load = cifti.read(data_file)
data = data_file_load[0]
data_to_analyse = data[:,int(start_idx):int(end_idx)]

# Create stimulus design
visual_dm_file = scipy.io.loadmat(opj(base_dir,'raw_data','stim','retinotopysmall5.mat'))
visual_dm = visual_dm_file['stim']

stimulus = VisualStimulus(  stim_arr = visual_dm,
                            viewing_distance = analysis_info["screen_distance"],
                            screen_width = analysis_info["screen_width"],
                            scale_factor = 1/10.0,
                            tr_length = analysis_info["TR"],
                            dtype = np.short)

# Initialize css model
if fit_model == 'gauss':
    fit_func = og.GaussianFit
    num_est = 6
    model_func = og.GaussianModel(  stimulus = stimulus,
                                    hrf_model = utils.spm_hrf)
elif fit_model == 'css':
    fit_func = css.CompressiveSpatialSummationFit
    num_est = 7
    model_func = css.CompressiveSpatialSummationModel(  stimulus = stimulus,
                                                        hrf_model = utils.spm_hrf)

model_func.hrf_delay = 0
print('models and stimulus loaded')

# Fit: define search grids
x_grid = (-12, 12)
y_grid = (-12, 12)
sigma_grid = (0.05, 15)
n_grid =  (0.01, 1.5)

# Fit: define search bounds
x_bound = (-30.0, 30.0)
y_bound = (-30.0, 30.0)
sigma_bound = (0.001, 70.0)
n_bound = (0.01, 3)
beta_bound = (-1e3, 1e3)
baseline_bound = (-1e3, 1e3)

if fit_model == 'gauss':
    fit_model_grids =  (x_grid, y_grid, sigma_grid)
    fit_model_bounds = (x_bound, y_bound, sigma_bound, beta_bound, baseline_bound)
elif fit_model == 'css':
    fit_model_grids =  (x_grid, y_grid, sigma_grid, n_grid)
    fit_model_bounds = (x_bound, y_bound, sigma_bound, n_bound, beta_bound, baseline_bound)

# Fit: define empty estimate and voxel indeces
estimates = np.zeros((num_est,data.shape[1]))
vertex_indices = [(xx, 0, 0) for xx in np.arange(int(start_idx),int(end_idx),1)]


# Define multiprocess bundle
bundle = utils.multiprocess_bundle( Fit = fit_func,
                                    model = model_func,
                                    data = data_to_analyse.T,
                                    grids = fit_model_grids, 
                                    bounds = fit_model_bounds, 
                                    indices = vertex_indices, 
                                    auto_fit = True, 
                                    verbose = 1, 
                                    Ns = Ns)
# Run fitting
pool = multiprocessing.Pool(processes = N_PROCS)
output = pool.map(  func = utils.parallel_fit, 
                    iterable = bundle)

for fit in output:
    estimates[:num_est-1,fit.voxel_index[0]] = fit.estimate
    estimates[num_est-1,fit.voxel_index[0]] = fit.rsquared

# Free up memory
pool.close()
pool.join()

# Save estimates data
bm_full = data_file_load[1][1]
series = cifti.Series(start=0, step=1, size=num_est)
cifti.write(opfn_est, estimates, (series, bm_full))