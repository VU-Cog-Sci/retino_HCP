"""
-----------------------------------------------------------------------------------------
extract_sum.py
-----------------------------------------------------------------------------------------
Goal of the script:
Extract statistics summary for each participants
# 0. ecc vs. size correlation r2
# 1. ecc vs. size slope
# 2. ecc vs. size intercept
# 3. contra-laterality index,
# 4. median of r2 of the fit
# 5. mean of r2 of the fit
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject number
sys.argv[2]: fit model ('gauss','css')
-----------------------------------------------------------------------------------------
Output(s):
hdf5_file summary.h5
-----------------------------------------------------------------------------------------
To run:
source activate i27
cd /home/szinte/projects/retino_HCP
python post_fit/extract_sum.py gauss
-----------------------------------------------------------------------------------------
"""

# Stop warnings
# -------------
import warnings
warnings.filterwarnings("ignore")

# General imports
# ---------------
import os
import sys
import json
import glob
import numpy as np
import ipdb
import time
import platform
import h5py
import scipy
from scipy.optimize import curve_fit
opj = os.path.join
deb = ipdb.set_trace

# Define analysis parameters
# --------------------------
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
    
# Define cluster/server specific parameters
# -----------------------------------------
if 'aeneas' in platform.uname()[1]:
    base_dir = analysis_info['aeneas_base_folder'] 
elif 'lisa' in platform.uname()[1]:
    base_dir = analysis_info['lisa_base_folder'] 
    
sign_idx, rsq_idx, ecc_idx, polar_real_idx, polar_imag_idx , size_idx, \
            non_lin_idx, amp_idx, baseline_idx, cov_idx, x_idx, y_idx = 0,1,2,3,4,5,6,7,8,9,10,11

# Get inputs
# ----------
subject = sys.argv[1]
fit_model = sys.argv[2]


# Create h5 for all subject
# -------------------------
summary_hdf5 = opj(base_dir,'pp_data',subject,fit_model,'h5')
try: os.makedirs(summary_hdf5)
except: pass
summary_hdf5_file = opj(base_dir,'pp_data',subject,fit_model,'h5','summary.h5')
try: 
    os.system('rm '+ summary_hdf5_file)
except: 
    pass
h5file = h5py.File(summary_hdf5_file, "a")


# Define r pearson weighted by fit rsquare
# ----------------------------------------
def m(x, w):
    return np.sum(x * w) / np.sum(w)

def cov(x, y, w):
    return np.sum(w * (x - m(x, w)) * (y - m(y, w))) / np.sum(w)

def weighted_corr(x, y, w):
    return cov(x, y, w) / np.sqrt(cov(x, x, w) * cov(y, y, w))


# Get summary values
# ------------------
# 0. ecc vs. size correlation r2
# 1. ecc vs. size slope
# 2. ecc vs. size intercept
# 3. contra-laterality index,
# 4. median of r2 of the fit
# 5. mean of r2 of the fit


print('summary processing of: %s'%subject)
for roi_num, roi in enumerate(analysis_info['rois']):
    for mask_dir in ['pos','neg']:
        deriv_mat_hemi = []
        deriv_mat = np.array([])
        summary_mat = np.array([])
        val_hemi = 0
        for hemi in ['L', 'R','LR']:
            h5_dir = opj(base_dir,'pp_data',subject,fit_model,'h5')

            # load data
            if hemi == 'LR':
                deriv_mat = np.row_stack((deriv_mat_hemi[0],deriv_mat_hemi[1]))
                draw = True
            else:
                # load derivatives
                if hemi == 'L': val_hemi = 1
                elif hemi == 'R': val_hemi = 2
                folder_alias = '{hemi}_{mask_dir}'.format(hemi = hemi,mask_dir = mask_dir)
                h5_file = h5py.File(opj(h5_dir,'{roi}.h5'.format(roi = roi)), "r")
                in_file = opj("prf_deriv_{hemi}_{mask_dir}".format(hemi = hemi, mask_dir = mask_dir))
                deriv_mat = h5_file['{folder_alias}/{in_file}'.format(folder_alias=folder_alias,in_file=in_file)]
                deriv_mat = np.vstack((deriv_mat,val_hemi*np.ones((1,deriv_mat.shape[1]))))
                deriv_mat = deriv_mat[:,:].T
                deriv_mat_hemi.append(deriv_mat)

        vertex_ini = deriv_mat.shape[0]
        if vertex_ini > 0:
            data4mask = deriv_mat
            deriv_mat = deriv_mat[np.logical_and(np.logical_and( data4mask[:,rsq_idx]>=analysis_info['rsq_threshold'],
                                                                 data4mask[:,cov_idx]>=analysis_info['cov_threshold']),
                                                                 data4mask[:,size_idx]>=analysis_info['size_threshold'])]
            vertex = deriv_mat.shape[0]

            if np.round(vertex) == 0:
                summary_mat = np.array([np.nan,np.nan,np.nan,np.nan,np.nan,np.nan])
            else:
                # compute correlation size/ecc
                if vertex >= 2:
                    linear_function = lambda x, a, b: a * x + b

                    ecc_data = deriv_mat[:,ecc_idx]
                    size_data = deriv_mat[:,size_idx]
                    weight_data = deriv_mat[:,rsq_idx]
                    coeffs, matcov = curve_fit( f = linear_function,
                                                xdata = ecc_data,
                                                ydata = size_data,
                                                sigma = weight_data)
                    size_fit = linear_function(ecc_data, coeffs[0], coeffs[1])

                    ecc_size_r = weighted_corr(ecc_data,size_data,weight_data)
                    ecc_size_slope = coeffs[0]
                    ecc_size_intercept = coeffs[0]
                else:
                    ecc_size_r = np.nan
                    ecc_size_slope = np.nan
                    ecc_size_intercept = np.nan

                # compute laterality index                    
                data_L = deriv_mat[deriv_mat[:,12]==1,:]
                contra_lat_L = np.nansum(data_L[data_L[:,x_idx] > 0,rsq_idx])/np.nansum(data_L[:,rsq_idx])
                data_R = deriv_mat[deriv_mat[:,12]==2,:]
                contra_lat_R = np.nansum(data_R[data_R[:,x_idx] < 0,rsq_idx])/np.nansum(data_R[:,rsq_idx])
                contra_lat = np.nanmean([contra_lat_L,contra_lat_R])

                # median r2 fit
                median_fit_r2 = np.nanmedian(deriv_mat[:,rsq_idx])
                mean_fit_r2   = np.nanmean(deriv_mat[:,rsq_idx])

                summary_mat = np.array([ecc_size_r,ecc_size_slope,ecc_size_intercept,contra_lat,median_fit_r2,mean_fit_r2])
        else:
            summary_mat = np.array([np.nan,np.nan,np.nan,np.nan,np.nan,np.nan])        


        # creat subject folder
        h5file.create_dataset('{subject}/{roi}/{mask_dir}_summary'.format(subject = subject, roi = roi, mask_dir = mask_dir),data = summary_mat,dtype='float32')
        time.sleep(0.1)