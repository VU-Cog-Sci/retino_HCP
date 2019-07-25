"""
-----------------------------------------------------------------------------------------
launch_pp_roi.py
-----------------------------------------------------------------------------------------
Goal of the script:
run pp_roi codes for each subjects of hcp dataset
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: fit model ('gauss','css')
sys.argv[2]: max tmux at the time
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
source activate i27
cd /home/szinte/projects/retino_HCP
python post_fit/launch_pp_roi.py gauss 8
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
import time
import nibabel as nb
import cifti

opj = os.path.join
deb = ipdb.set_trace


# Get inputs
# ----------
fit_model = sys.argv[1]
# max_tmux = int(sys.argv[2])
job_vox = 5000

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
base_file_name = 'tfMRI_RETBAR1_7T_AP_Atlas_MSMAll_hp2000_clean.dtseries'


# Get subject list
# ----------------
subject_list = analysis_info['subject_list']

# Run pp_roi.py
# -------------
num_run = 0
for subject_fit in subject_list:
    # if num_run < max_tmux:
    #     if os.path.isdir(opj(base_dir,'pp_data',subject_fit,fit_model,'deriv'))==0:
    #         session_name = "{subject_fit}_{fit_model}_pp_roi".format(subject_fit = subject_fit, fit_model = fit_model)
    #         print("tmux new-session -d -s {session_name} 'python post_fit/pp_roi.py {subject_fit} {fit_model} {job_vox} 0'".format(\
    #             session_name = session_name, subject_fit = subject_fit, fit_model = fit_model, job_vox = job_vox))
    #         os.system("tmux new-session -d -s {session_name} 'python post_fit/pp_roi.py {subject_fit} {fit_model} {job_vox} 0'".format(\
    #             session_name = session_name, subject_fit = subject_fit, fit_model = fit_model, job_vox = job_vox))
    #         time.sleep(2)
    #         num_run = num_run + 1
    
    print("python post_fit/pp_roi.py {subject_fit} {fit_model} {job_vox} 0".format(subject_fit = subject_fit, fit_model = fit_model, job_vox = job_vox))
    os.system("python post_fit/pp_roi.py {subject_fit} {fit_model} {job_vox} 0".format(subject_fit = subject_fit, fit_model = fit_model, job_vox = job_vox))

