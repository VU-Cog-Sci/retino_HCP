"""
-----------------------------------------------------------------------------------------
prf_makegrid.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create pRF grid estimation
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: model
sys.argv[2]: grid prediction file
sys.argv[3]: main directory
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
import nibabel as nb
from prf_utils import *

# Get inputs
fit_model = sys.argv[1]
subject = sys.argv[2]
start_idx = sys.argv[3]
end_idx = sys.argv[4]
data_file = sys.argv[5]
base_dir = sys.argv[6]

# Define analysis parameters
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
    
# Create stimulus design
visual_dm_file = scipy.io.loadmat(opj(base_dir,'raw_data','vis_design.mat'))
visual_dm = visual_dm_file['stim']

# Load data
data = []
data_file_load = nb.load(data_file)
data.append(np.array([data_file_load.darrays[i].data for i in range(len(data_file_load.darrays))]))
data = np.vstack(data)
data_to_analyse = data[:,int(start_idx):int(end_idx)]

# Define output file path and directories
base_file_name = os.path.split(data_file)[-1][:-4]
opfn_est = opj(base_dir,'pp_data',subject,fit_model,'gridfit',base_file_name + '_est_%s_to_%s.gii' %(start_idx,end_idx))
try: os.makedirs(opj(base_dir,'pp_data',subject,fit_model,'gridfit'))
except: pass

# Fit: define search grids
x_grid_bound = (-13.4, 13.4)
y_grid_bound = (-7.5, 7.5)
sigma_grid_bound = (0.05, 9.5)
n_grid_bound = (0.01, 1.5)

# Fit: define search bounds
x_fit_bound = (-30.0, 30.0)
y_fit_bound = (-30.0, 30.0)
sigma_fit_bound = (0.001, 70.0)
n_fit_bound = (0.01, 3)
beta_fit_bound = (-1e3, 1e3)
baseline_fit_bound = (-1e3, 1e3)

if fit_model == 'gauss' or fit_model == 'gauss_sg':
    bound_grids  = (x_grid_bound, y_grid_bound, sigma_grid_bound)
    bound_fits = (x_fit_bound, y_fit_bound, sigma_fit_bound, beta_fit_bound, baseline_fit_bound)
elif fit_model == 'css' or fit_model == 'css_sg':
    bound_grids  = (x_grid_bound, y_grid_bound, sigma_grid_bound, n_grid_bound)
    bound_fits = (x_fit_bound, y_fit_bound, sigma_fit_bound, n_fit_bound, beta_fit_bound, baseline_fit_bound)
    
# Initialize the prf model
prf = prf_fit(  fit_model = fit_model, 
                visual_design = visual_dm, 
                screen_distance = analysis_info["screen_distance"],
                screen_width = analysis_info["screen_width"],
                tr =  analysis_info["TR"],
                grid_steps = analysis_info["grid_steps"],
                bound_grids = bound_grids,
                bound_fits = bound_fits,
                hrf_delay = analysis_info["hrf_delay"],
                sg_filter_window_length = analysis_info["sg_filt_window_length"],
                sg_filter_polyorder = analysis_info["sg_filt_polyorder"],
                sg_filter_deriv = analysis_info["sg_filt_deriv"], 
                )

# Load the grid
grid_prediction_file = opj(base_dir,'pp_data','makegrid',fit_model,'grid_predictions.hdf5')
prf.load_grid(save_file = grid_prediction_file)

# Fit the grid
prf.fit_grid(data = data_to_analyse)

# Save estimates data
darrays = [nb.gifti.gifti.GiftiDataArray(d) for d in prf.gridsearch_all]
gii_out = nb.gifti.gifti.GiftiImage(header = data_file_load.header, 
                                    extra = data_file_load.extra,
                                    darrays = darrays)
nb.save(gii_out, opfn_est)