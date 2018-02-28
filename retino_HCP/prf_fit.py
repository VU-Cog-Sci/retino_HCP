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

# from .utils import *

subject_folder = '/home/shared/2018/visual/HCP7TFIXED/671855/'
# subject_folder = sys.argv[1]


TR = 1.0
screen_distance = 101
screen_width = 29
n_timepoints_per_run = 300
run_order = ['RETCCW','RETCW','RETEXP','RETCON','RETBAR1','RETBAR2']
direction_order = ['CCW','CW','EXP','CON','PRF','PRF']

############################################################################################################################################
#
#   create stimulus design matrices 
#
############################################################################################################################################

visual_dm = []
for i, d in enumerate(direction_order):
    if i in (0,1):
        visual_dm.append(design_matrix_wedge(direction=d))
    if i in (2,3):
        visual_dm.append(design_matrix_ring(direction=d))
    if i in (4,5):
        visual_dm.append(design_matrix_prf())

visual_dm = np.vstack(visual_dm).transpose((1,2,0))

stimulus = VisualStimulus(stim_arr=visual_dm, 
                            viewing_distance=screen_distance, 
                            screen_width=screen_width, 
                            scale_factor=1.0/5.0, 
                            tr_length=TR, 
                            dtype=np.short)#, interp='nearest'




############################################################################################################################################
#
#   load cii data, timepoints by grayordinates
#
############################################################################################################################################

cii_files = [glob.glob(os.path.join(subject_folder, '*%s*_sg_psc.nii'%run))[0] for run in run_order]

data = []
for cii_file in cii_files:
    cii_in = nb.load(cii_file)
    data.append(cii_in.get_data())

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
x_grid = (-10, 10)
y_grid = (-10, 10)
s_grid = (0.25, 7.25)
n_grid = (0.45, 1.05)      # nonlinearity for css
b_grid = (-5.0, 5.0)
bas_grid = (-1.5, 1.5)       # baseline for css

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
                                   grids=css_grids, bounds=css_bounds, indices=voxel_indices, auto_fit=True, verbose=1, Ns=6)

output = Parallel(n_jobs=3)(delayed(css.CompressiveSpatialSummationFit)(model=css_model,
                                    data=d,
                                    grids=css_grids, 
                                    bounds=css_bounds, 
                                    voxel_index=vi, 
                                    auto_fit=True, 
                                    verbose=1, 
                                    Ns=6) for vi, d in zip(voxel_indices, data.T))


# run analysis
pool = multiprocessing.Pool(3)
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