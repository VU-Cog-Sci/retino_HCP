"""
-----------------------------------------------------------------------------------------
prf_fit.py
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
import cifti

# MRI analysis imports
import nibabel as nb
import popeye.utilities as utils
from popeye.visual_stimulus import VisualStimulus
import popeye.css as css
import popeye.og as og

# Get inputs
fit_model = 'gauss'
subject = '999999'
start_idx = 1
end_idx = 5
data_file = '/Users/macbook/disks/ae_Home/hcp_code/tfMRI_RETALL_Atlas_1.6mm_MSMAll_hp2000_clean_sg_psc_all.dtseries.nii'
base_dir = '/Users/macbook/disks/ae_Shared/2018/visual/subcortical_hcp'


# Define analysis parameters
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define cluster/server specific parameters
if 'lisa' in platform.uname()[1]:
    N_PROCS = 16
elif 'aeneas' in platform.uname()[1]:
    N_PROCS = 8 #31
elif 'local' in platform.uname()[1]:
    N_PROCS = 8

# Define output file path and directories
base_file_name = os.path.split(data_file)[-1][:-7]
opfn_est = opj(base_dir, base_file_name + '_est_%s_to_%s.nii.' %(start_idx,end_idx))
#opfn_est = opj(base_dir,'pp_data',subject,fit_model,'fit',base_file_name + '_est_%s_to_%s.nii.' %(start_idx,end_idx))

try : os.makedirs(opj(base_dir))
#try: os.makedirs(opj(base_dir,'pp_data',subject,fit_model,'fit'))
except: pass

#loading with cifti
data = cifti.read(data_file)
data_file_shape = data.shape

# Load data
#data = []
#data_file_load = nb.load(data_file)
#data_file_dat = data_file_load.get_data()
#data_file_shape = data_file_load.shape

# load mask
#maskfn = opj(base_dir,'raw_data','RETBAR_ALL_tfMRI_data_sub_mask.nii.gz')
#data_mask = nb.load(maskfn).get_data()

#data_file_masked = data_file_dat[data_mask==1.0]
#data_to_analyse = data_file_masked[int(start_idx):int(end_idx),:]
#data_to_analyse_idx = np.where(data_mask==1.0)
#voxel_indices = [(xx, yy, zz) for xx,yy,zz in  zip(data_to_analyse_idx[0][:],data_to_analyse_idx[1][:],data_to_analyse_idx[2][:])]
#voxel_indices = voxel_indices[int(start_idx):int(end_idx)]

vertex_indices = [(xx, 0, 0) for xx in np.arange(int(start_idx),int(end_idx),1)]

# Create stimulus design
visual_dm_file = scipy.io.loadmat(opj(base_dir,'stim','retinotopysmall_all.mat'))
#visual_dm_file = scipy.io.loadmat(opj(base_dir,'raw_data','retinotopysmall_all.mat'))
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
#estimates = np.zeros((data_file_shape[0], data_file_shape[1], data_file_shape[2], num_est))

#data_to_analyse
index_start = 1
index_end = 5
data_to_analyse = data.get_data()[:,int(index_start):int(index_end)]
data_to_analyse.shape

# Define multiprocess bundle
bundle = utils.multiprocess_bundle( Fit = fit_func,
                                    model = model_func,
                                    data = data_to_analyse.T,
                                    grids = fit_model_grids, 
                                    bounds = fit_model_bounds, 
                                    indices = vertex_indices,
                                    auto_fit = True,
                                    verbose = 1,
                                    Ns = 6)
# Run fitting
pool = multiprocessing.Pool(processes = N_PROCS)
output = pool.map(  func = utils.parallel_fit, 
                    iterable = bundle)

# Save estimates data

for fit in output:
    estimates[:num_est-1,fit.voxel_index[0]] = fit.estimate
    estimates[num_est-1,fit.voxel_index[0]] = fit.rsquared

#for fit in output:
    # Fill estimates matrix
    #estimates[fit.voxel_index][:num_est-1] = fit.estimate
    #estimates[fit.voxel_index][num_est-1] = fit.rsquared

# Free up memory
pool.close()
pool.join()

#saving with cifti
#gii_out = nb.cifti2.cifti2.Cifti2Image(dataobj = estimates,
                                        #header = data.header, 
                                        #extra = data.extra)
#nb.cifti2.cifti2.save(gii_out, opfn_est)

series = cifti.Series(start=0, step=analysis_info["TR"]*1000, size=data_file_shape[0])
bm_full = data[1][1]
cifti.write(opfn_est, estimates, (series, bm_full))
    
