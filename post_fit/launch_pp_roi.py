"""
-----------------------------------------------------------------------------------------
launch_pp_roi.py
-----------------------------------------------------------------------------------------
Goal of the script:
run pp_roi codes for each subjects of hcp dataset
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: fit model ('gauss','css')
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
source activate i27
cd /home/szinte/projects/retino_HCP
python post_fit/launch_pp_roi.py gauss
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

# Check system
# ------------
sys.exit('Drawing Flatmaps only works with Python 2. Aborting.') if sys.version_info[0] > 2 else None

# Get inputs
# ----------
fit_model = sys.argv[1]
vox_num = 32492
job_vox = 2500

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

# Get subject list with fit ready
# -------------------------------
subject_fit_list = []

for subject in subject_list:
	
	if os.path.isdir(opj(base_dir,'pp_data',subject,fit_model,'fit')):
		start_idx =  np.arange(0,vox_num,job_vox)
		end_idx = start_idx+job_vox
		end_idx[-1] = vox_num
		num_miss_part = 0
		fit_est_files_L = []
		fit_est_files_R = []
		for hemi in ['L','R']:
		    for iter_job in np.arange(0,start_idx.shape[0],1):
		        fit_est_file = opj(base_dir,'pp_data',subject,fit_model,'fit', '%s_%s.func_bla_psc_est_%s_to_%s.gii' %(base_file_name,hemi,str(int(start_idx[iter_job])),str(int(end_idx[iter_job]))))
		        if os.path.isfile(fit_est_file):
		            if os.path.getsize(fit_est_file) == 0:
		                num_miss_part += 1 
		            else:
		                exec('fit_est_files_{hemi}.append(fit_est_file)'.format(hemi = hemi))
		        else:
		            num_miss_part += 1
		
		if num_miss_part == 0:
			subject_fit_list.append(subject)

# Run pp_roi.py
# -------------
for subject_fit in subject_fit_list:
	if os.path.isdir(opj(base_dir,'pp_data',subject_fit,fit_model,'figs','prf'))==0:
		session_name = "{subject_fit}_{fit_model}".format(subject_fit = subject_fit, fit_model = fit_model)
		print("tmux new-session -d -s {session_name} 'python post_fit/pp_roi.py {subject_fit} {fit_model} {job_vox}'".format(session_name = session_name, subject_fit = subject_fit, fit_model = fit_model, job_vox = job_vox))
		os.system("tmux new-session -d -s {session_name} 'python post_fit/pp_roi.py {subject_fit} {fit_model} {job_vox}'".format(session_name = session_name, subject_fit = subject_fit, fit_model = fit_model, job_vox = job_vox))
		