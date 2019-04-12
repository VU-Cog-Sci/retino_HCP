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
grid_prediction_file = sys.argv[2]
base_dir = sys.argv[3]

# Define analysis parameters
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define cluster/server specific parameters
if 'lisa' in platform.uname()[1]:
    N_PROCS = 16
elif 'aeneas' in platform.uname()[1]:
    N_PROCS = 2
    
# Create stimulus design
visual_dm_file = scipy.io.loadmat(opj(base_dir,'raw_data','vis_design.mat'))
visual_dm = visual_dm_file['stim']

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
                sg_filter_window_length = analysis_info["sg_filt_window_length"],
                sg_filter_polyorder = analysis_info["sg_filt_polyorder"],
                sg_filter_deriv = analysis_info["sg_filt_deriv"], 
                )


# Make and save the grid
prf.make_grid(save_file = grid_prediction_file)