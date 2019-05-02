"""
-----------------------------------------------------------------------------------------
pp_roi.py
-----------------------------------------------------------------------------------------
Goal of the script:
Region of interests pre-processing
Compute pRF derivatives and plot on pycortex overlay.svg to determine visual ROI
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject number
sys.argv[2]: fit model ('gauss','css')
sys.argv[3]: voxels per fit (e.g 2500)
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
source activate i27
cd /home/szinte/projects/retino_HCP
python post_fit/pp_roi.py 999999 gauss 2500
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
import matplotlib.pyplot as pl
import ipdb
import platform
opj = os.path.join
deb = ipdb.set_trace

# MRI imports
# -----------
import nibabel as nb
import cortex

# Functions import
# ----------------
from utils import convert_fit_results

# Get inputs
# ----------
subject = sys.argv[1]
fit_model = sys.argv[2]
job_vox = float(sys.argv[3])
if fit_model == 'gauss': fit_val = 6
elif fit_model == 'css': fit_val = 7
base_file_name = 'RETBAR_ALL_tfMRI_data_sub'

# Define analysis parameters
# --------------------------
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define cluster/server specific parameters and folders
# -----------------------------------------------------
if 'aeneas' in platform.uname()[1]:
    base_dir = analysis_info['aeneas_base_folder'] 
elif 'local' in platform.uname()[1]:
    base_dir = analysis_info['local_base_folder'] 
deriv_dir = opj(base_dir,'pp_data',subject,fit_model,'deriv')
prf_filename = opj(base_dir,"pp_data",subject,fit_model,"fit","{bfn}_est.nii.gz".format(bfn = base_file_name))

# Check if all slices are present
# -------------------------------
maskfn = opj(base_dir,'raw_data','RETBAR_ALL_tfMRI_data_sub_mask.nii.gz')
data_mask = nb.load(maskfn).get_data()
start_idx =  np.arange(0,np.sum(data_mask),job_vox)
end_idx = start_idx+job_vox
end_idx[-1] = int(np.sum(data_mask))

num_miss_part = 0
fit_est_files = []
for iter_job in np.arange(0,start_idx.shape[0],1):
    fit_est_file = opj(base_dir,'pp_data',subject,fit_model,'fit', '%s_est_%s_to_%s.nii.gz' %(base_file_name,str(int(start_idx[iter_job])),str(int(end_idx[iter_job]))))
    if os.path.isfile(fit_est_file):
        if os.path.getsize(fit_est_file) == 0:
            num_miss_part += 1 
        else:
            fit_est_files.append(fit_est_file)
    else:
        num_miss_part += 1

if num_miss_part != 0:
    sys.exit('%i missing files, analysis stopped'%num_miss_part)

# Combine fit files
# -----------------
print('combining fit files')
data_combined = np.zeros((nb.load(fit_est_files[0]).shape[0],nb.load(fit_est_files[0]).shape[1],nb.load(fit_est_files[0]).shape[2],fit_val))

for fit_filename in fit_est_files:
    data_fit = []
    data_fit = nb.load(fit_filename).get_data()
    data_combined = data_combined + data_fit
img = nb.Nifti1Image(   dataobj = data_combined,
                        affine = nb.load(fit_est_files[0]).affine,
                        header =nb.load(fit_est_files[0]).header)
img.to_filename(prf_filename)

# Compute derived measures from prfs
# ----------------------------------
print('extracting pRF derivatives')
convert_fit_results(prf_filename = prf_filename,
                    output_dir = deriv_dir,
                    stim_radius = analysis_info['stim_radius'],
                    fit_model = fit_model,
                    mask_filename = maskfn)

# Compute masked data
# -------------------
for mask_dir in ['all','pos','neg']:
    exec('derivfn = opj(deriv_dir,"{mask_dir}","prf_deriv_{mask_dir}.nii.gz"))'.format(mask_dir = mask_dir))
    prf_deriv_load = nb.load(prf_filename)
    prf_deriv = prf_data_load.get_data()

    # threshold with R2

    # threshold with cov

    # threhsold with cov and R2




