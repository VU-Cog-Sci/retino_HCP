from __future__ import division
import sys
import ctypes
import multiprocessing
import numpy as np
import popeye.utilities as utils
from popeye.visual_stimulus import VisualStimulus, simulate_bar_stimulus
import popeye.css_nohrf as css
from skimage.transform import rotate
import tables
import platform

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

if 'lisa' in platform.uname()[1]:
    base_dir = analysis_info['lisa_cluster_base_folder'] 
    N_PROCS = 15
    print('on lisa')
elif 'localhost' in platform.uname()[1]:
    base_dir = analysis_info['lisa_cluster_base_folder'] 
    N_PROCS = 250
    print('on ascanius')
else:
    base_dir = analysis_info['cartesius_cluster_base_folder'] 
    N_PROCS = 23
    print('on cartesius')
subject = str(sys.argv[1])
hemi = str(sys.argv[2])

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
# visual_dm = []
# for i, d in enumerate(analysis_info["direction_order"]):
#     if i == 0:
#         visual_dm.append(design_matrix_wedge(direction=d))
#     if i == 2:
#         visual_dm.append(design_matrix_ring(direction=d))
#     if i == 4:
#         visual_dm.append(design_matrix_prf())

# visual_dm = np.vstack(visual_dm).transpose((1,2,0))

# the code above recreates Kendrick's design matrices, we'll now just load them.
visual_dm = []
# for i in [1,3,5]:
for i in [5]:
    file = tables.open_file(os.path.join(base_dir, 'retinotopysmall{i}.mat'.format(i=i)))
    visual_dm.append(file.get_node('/stim')[:])
    file.close()
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

# averaged_runs = ['CCW','EXP','RETBAR1']
averaged_runs = ['RETBAR1']
gii_files = [glob.glob(os.path.join(subject_folder, '*{run}*_{hemi}.func_bla_psc_av.gii'.format(run=run, hemi=hemi)))[0] for run in averaged_runs]

data = []
for gii_file in gii_files:
    gii_in = nb.load(gii_file)
    data.append(np.array([gii_in.darrays[i].data for i in range(len(gii_in.darrays))]))

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

print("starting fitting of {subject}, hemi {hemi}".format(subject=subject, hemi=hemi))

bundle = utils.multiprocess_bundle(Fit=css.CompressiveSpatialSummationFit, model=css_model, data=data.T,
                                   grids=css_grids, bounds=css_bounds, indices=voxel_indices, auto_fit=True, verbose=1, Ns=12)

# run analysis
pool = multiprocessing.Pool(N_PROCS)
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

darrays = [nb.gifti.gifti.GiftiDataArray(d) for d in estimates]

gii_out = nb.gifti.gifti.GiftiImage(header=gii_in.header, 
                        extra=gii_in.extra,
                        darrays=darrays)

out_name = os.path.splitext(gii_files[0])[0] + '_est.gii'
out_file = os.path.abspath(out_name)
nb.save(gii_out, out_file)
