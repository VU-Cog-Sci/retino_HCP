from __future__ import division
import sys
import ctypes
import multiprocessing
import numpy as np
import popeye.utilities as utils
from popeye.visual_stimulus import VisualStimulus, simulate_bar_stimulus
import popeye.css_nohrf as css
from skimage.transform import rotate

from skimage.morphology import disk
import nibabel as nb
import matplotlib.pyplot as pl
from math import *
import json

import os
import glob
import gc
from IPython import embed as shell
from joblib import Parallel, delayed

from hrf_estimation.hrf import spmt  # , dspmt, ddspmt

from utils import *


with open('../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

base_dir = analysis_info['cluster_base_folder'] 
subject = sys.argv[1]

subject_folder = os.path.join(base_dir, subject)

############################################################################################################################################
#
#   create stimulus design matrices 
#
############################################################################################################################################

# visual_dm = []
# for i, d in enumerate(analysis_info["direction_order"]):
#     if i in (0,1):
#         visual_dm.append(design_matrix_wedge(direction=d))
#     if i in (2,3):
#         visual_dm.append(design_matrix_ring(direction=d))
#     if i in (4,5):
#         visual_dm.append(design_matrix_prf())

# using only the even runs which are the 'standards'
visual_dm = []
for i, d in enumerate(analysis_info["direction_order"]):
    if i == 0:
        visual_dm.append(design_matrix_wedge(direction=d))
    if i == 2:
        visual_dm.append(design_matrix_ring(direction=d))
    if i == 4:
        visual_dm.append(design_matrix_prf())

visual_dm = np.vstack(visual_dm).transpose((1,2,0))

stimulus = VisualStimulus(stim_arr=visual_dm, 
                            viewing_distance=analysis_info["screen_distance"], 
                            screen_width=analysis_info["screen_width"], 
                            scale_factor=0.05, 
                            tr_length=analysis_info["TR"], 
                            dtype=np.short)


############################################################################################################################################
#
#   load cii data, timepoints by grayordinates
#
############################################################################################################################################

averaged_runs = ['CCW','EXP','BOTHBARS']
cii_files = [glob.glob(os.path.join(subject_folder, '*%s*_sg_psc_av.nii'%run))[0] for run in averaged_runs]

data = []
for cii_file in cii_files:
    cii_in = nb.load(cii_file)
    data.append(cii_in.get_data())

data = np.vstack(data)

data = data.squeeze()
data = np.vstack(data)

estimates = np.zeros((7,data.shape[1]))

############################################################################################################################################
#
#   then, creating the model to be fit
#
############################################################################################################################################


def my_spmt(delay, tr):
    return spmt(np.arange(0, 33, tr))

# MODEL
# initialize css model
css_model = css.CompressiveSpatialSummationModel(stimulus, my_spmt)
css_model.hrf_delay = 0


# FIT
# define search grids
# these define min and max of the edge of the initial brute-force search.
x_grid = (-12.5, 12.5)
y_grid = (-12.5, 12.5)
s_grid = (0.25, 7.25)
n_grid = (0.45, 1.05)      # nonlinearity for css
b_grid = (-2.5, 2.5)
bas_grid = (-1.0, 1.0)       # baseline for css

# define search bounds
# these define the boundaries of the final gradient-descent search.
x_bound = (-100.0, 100.0)
y_bound = (-100.0, 100.0)
s_bound = (0.001, 70.0)
n_bound = (0.001, 1.5)     # nonlinearity for css
b_bound = (-1e3, 1e3)
bas_bound = (-1e3, 1e3)    # baseline for css

# order of css estimate parameters:
# [self.theta,self.rho,self.sigma_size,self.n,self.beta,self.baseline]
css_grids = (x_grid, y_grid, s_grid, n_grid, b_grid, bas_grid)
css_bounds = (x_bound, y_bound, s_bound, n_bound, b_bound, bas_bound)


############################################################################################################################################
#
#   actual fitting
#
############################################################################################################################################

voxel_indices = [(xx, 0, 0) for xx in np.arange(data.shape[1])]

bundle = utils.multiprocess_bundle(Fit=css.CompressiveSpatialSummationFit, model=css_model, data=data.T,
                                   grids=css_grids, bounds=css_bounds, indices=voxel_indices, auto_fit=True, verbose=1, Ns=5)

# run analysis
pool = multiprocessing.Pool(23)
output = pool.map(utils.parallel_fit, bundle)

for fit in output:
    estimates[:6,fit.voxel_index[0]] = fit.estimate
    estimates[6,fit.voxel_index[0]] = fit.rsquared

# try to free up memory by closing the pool and joining them with the main thread
pool.close()
pool.join()


############################################################################################################################################
#
#   outputs
#
############################################################################################################################################


cii_out = nb.Cifti2Image(dataobj=estimates, 
                        header=cii_in.header, 
                        nifti_header=cii_in.nifti_header, 
                        extra=cii_in.extra)

out_name = os.path.splitext(cii_file)[0] + '_est.nii'
out_file = os.path.abspath(out_name)
nb.save(cii_out, out_file)