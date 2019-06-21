#To run:
#cd /home/ada/projects/retino_HCP
#python post_fit/pp_roi_cifti.py 999999 gauss 2500

import warnings
warnings.filterwarnings("ignore")

# General imports
# ---------------
import os
import sys
import json
import glob
import numpy as np
import matplotlib.pyplot as pl
import ipdb
import platform
import cifti
opj = os.path.join
deb = ipdb.set_trace
import cifti
# MRI imports
# -----------
import nibabel as nb
import cortex

# Functions import
# ----------------
import utils 

sys.exit('Popeye error with Python 2. Use Python 3 Aborting.') if sys.version_info[0] < 3 else None

subject = sys.argv[1]
fit_model = sys.argv[2]
job_vox = sys.argv[3]

if fit_model == 'gauss': fit_val = 6
elif fit_model == 'css': fit_val = 7
base_file_name = 'tfMRI_RETALL_Atlas_1.6mm_MSMAll_hp2000_clean_sg_psc'

# Define analysis parameters
# --------------------------
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define cluster/server specific parameters
# -----------------------------------------
if 'aeneas' in platform.uname()[1]:
    base_dir = analysis_info['aeneas_base_folder']
 
elif 'local' in platform.uname()[1]:
    base_dir = analysis_info['local_base_folder'] 
#Output directory:
deriv_dir = opj(base_dir,'pp_data',subject,fit_model,'deriv')

# determine number of vertex and time_series
data_file = opj(base_dir,'raw_data', subject, 'tfMRI_RETALL_Atlas_1.6mm_MSMAll_hp2000_clean_sg_psc.dtseries.nii')
data_file_load = cifti.read(data_file)
data = data_file_load[0]

# Check if all slices are present
start_idx =  np.arange(0,vox_num,job_vox)
end_idx = start_idx+job_vox
end_idx[-1] = vox_num
num_miss_part = 0

fit_est_files = []
for iter_job in np.arange(0,start_idx.shape[0],1):
    fit_est_file = opj(base_dir,'pp_data',subject,fit_model,'fit', '%s_est_%s_to_%s.dtseries.nii' %(base_file_name,str(int(start_idx[iter_job])),str(int(end_idx[iter_job]))))
    if os.path.isfile(fit_est_file):
        if os.path.getsize(fit_est_file) == 0:
            num_miss_part += 1 
        else:
            fit_est_files.append(fit_est_file)
    else:
        num_miss_part += 1

if num_miss_part != 0:
    sys.exit('%i missing files, analysis stopped'%num_miss_part)


#Combine fit files
print('combining fit files')
data_combined = np.zeros((fit_val, vox_num))
for fit_filename in fit_est_files:
    data_fit_file = cifti.read(fit_filename)
    data_fit = data_fit_file[0]
    data_fit.shape
    data_combined = data_combined + data_fit

prf_filename = opj(base_dir,'pp_data',subject,fit_model,'fit',"{bfn}.dtseries.nii".format(bfn= base_file_name))
bm_full = data_fit_file[1][1]
series = cifti.Series(start=0, step=1, size=fit_val)
cifti.write(prf_filename, data_combined, (series, bm_full)) 

# Compute derived measures from prfs
# ----------------------------------
print('extracting pRF derivatives')

utils.convert_fit_results(prf_filename = prf_filename,
                    output_dir = deriv_dir,
                    stim_radius = analysis_info['stim_radius'],
                    fit_model = fit_model)












