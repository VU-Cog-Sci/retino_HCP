import numpy as np
import nibabel as nb
import matplotlib.pyplot as pl
from math import *
import json

import os
import glob
import gc

import tables
from joblib import Parallel, delayed
from sklearn.linear_model import Ridge

import popeye.utilities as utils
from popeye.spinach import generate_og_timeseries
from popeye.visual_stimulus import VisualStimulus, simulate_bar_stimulus
import popeye.css as css
from hrf_estimation.hrf import spmt  # , dspmt, ddspmt


#################################################################################
#####
#####   gii workflow
#####
#################################################################################

with open('../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

base_dir = '/home/shared/2018/visual/HCP7TFIXED/'

subject = str(sys.argv[1])
hemi = str(sys.argv[2])

subject_folder = os.path.join(base_dir, subject)


#################################################################################
#####
#####   create stimulus
#####
#################################################################################


visual_dm = []
# for i in [1,3,5]:
for i in [5]:
    file = tables.open_file(os.path.join(base_dir, 'retinotopysmall{i}.mat'.format(i=i)))
    visual_dm.append(file.get_node('/stim')[:])
    file.close()
visual_dm = np.vstack(visual_dm).transpose((1,2,0))

stimulus = VisualStimulus(stim_arr=visual_dm[::2,::2],
                            viewing_distance=analysis_info["screen_distance"],
                            screen_width=analysis_info["screen_width"],
                            scale_factor=0.05,
                            tr_length=analysis_info["TR"],
                            dtype=np.short)


#################################################################################
#####
#####   get data
#####
#################################################################################

averaged_runs = ['RETBAR1']
gii_files = [glob.glob(os.path.join(subject_folder, '*{run}*_{hemi}.func_bla_psc_av.gii'.format(run=run, hemi=hemi)))[0] for run in averaged_runs]

data = []
for gii_file in gii_files:
    gii_in = nb.load(gii_file)
    data.append(np.array([gii_in.darrays[i].data for i in range(len(gii_in.darrays))]))

data = np.vstack(data)

data = data.squeeze()
data = np.vstack(data)

data = np.nan_to_num(data)
estimates = np.zeros((16,data.shape[1]))

#################################################################################
#####
#####   parameters for modeling
#####
#################################################################################

def my_spmt(delay, tr):
    return spmt(np.arange(0, 33, tr))

# MODEL
# initialize css model
css_model = css.CompressiveSpatialSummationModel(stimulus, my_spmt)
css_model.hrf_delay = 0

xy_parspace = np.sort(np.r_[np.logspace(0,1.7,20)-1, -(np.logspace(0,1.4,20)-1)])
x_pars, y_pars, sigma_pars, n_pars = np.meshgrid(xy_parspace, xy_parspace, np.logspace(0,1.4,10)-0.8, np.linspace(0.1,1.2,5))

x_pars = x_pars.ravel()
y_pars = y_pars.ravel()
sigma_pars = sigma_pars.ravel()
n_pars = n_pars.ravel()

x_space, y_space = np.meshgrid(np.linspace(-30,30,100), np.linspace(-30,30,100))

#################################################################################
#####
#####   loading/creating DM
#####
#################################################################################

if os.path.isfile('dm.npz'):
    regs = np.load('dm.npz')['arr_0']
else:
    regs = np.ones((data.shape[0], n_pars.shape[0]+1))

    i=1
    for ixr, iyr, isigmar, inr in zip(x_pars, y_pars, sigma_pars, n_pars):
        regs[:,i] = css_model.generate_prediction(x=ixr, y=iyr, n=inr, sigma=isigmar, beta=1, baseline=0, hrf_delay=0)
        i += 1
        if i%1000 == 0:
            print('bla %i'%i)

#################################################################################
#####
#####   ridge step, with peak/trough detection
#####
#################################################################################

clf = Ridge(alpha=1e15)
clf.fit(regs, data)
peak_pars = np.argmax(clf.coef_, axis = 1)
trough_pars = np.argmin(clf.coef_, axis = 1)

#################################################################################
#####
#####   per-voxel part of fitting
#####
#################################################################################

for i, d in enumerate(data.T):
    peak_x = np.array([np.ones(data.shape[0]), regs[:,peak_pars[i]]]).T
    peak_lr = Ridge(alpha=0)
    peak_lr.fit(peak_x, d.reshape(-1, 1))
    rsq_peak = peak_lr.score(peak_x, d.reshape(-1, 1))

    trough_x = np.array([np.ones(data.shape[0]), regs[:,trough_pars[i]]]).T
    trough_lr = Ridge(alpha=0)
    trough_lr.fit(trough_x, d.reshape(-1, 1))
    rsq_trough = trough_lr.score(trough_x, d.reshape(-1, 1))

    estimates[:,i] = np.array([x_pars[peak_pars[i]], 
                    y_pars[peak_pars[i]],
                    np.angle(x_pars[peak_pars[i]] + 1j * y_pars[peak_pars[i]]), 
                    sqrt(x_pars[peak_pars[i]]**2 + y_pars[peak_pars[i]]**2),
                    sigma_pars[peak_pars[i]], 
                    n_pars[peak_pars[i]], 
                    peak_lr.coef_[0,0], 
                    peak_lr.coef_[0,1], 
                    rsq_peak,
                    x_pars[trough_pars[i]], 
                    y_pars[trough_pars[i]], 
                    np.angle(x_pars[trough_pars[i]] + 1j * y_pars[trough_pars[i]]), 
                    sqrt(x_pars[trough_pars[i]]**2 + y_pars[trough_pars[i]]**2),
                    sigma_pars[trough_pars[i]], 
                    n_pars[trough_pars[i]], 
                    trough_lr.coef_[0,0], 
                    trough_lr.coef_[0,1], 
                    rsq_trough])

    if i%1000 == 0:
        print('fitting %i'%i)

#################################################################################
#####
#####   save out
#####
#################################################################################

darrays = [nb.gifti.gifti.GiftiDataArray(d) for d in estimates]

gii_out = nb.gifti.gifti.GiftiImage(header=gii_in.header, 
                        extra=gii_in.extra,
                        darrays=darrays)

out_name = os.path.splitext(gii_file)[0] + '_est.gii'
out_file = os.path.abspath(out_name)
nb.save(gii_out, out_file)









