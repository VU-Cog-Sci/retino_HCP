"""
-----------------------------------------------------------------------------------------
launch_roi_plots.py
-----------------------------------------------------------------------------------------
Goal of the script:
run roi_plots codes for each subjects of hcp dataset
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: fit model ('gauss','css')
sys.argv[2]: max tmux at the time
sys.argv[3]: do single hemifield plot (1 = YES, 0 = NO)
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
source activate i27
cd /home/szinte/projects/retino_HCP
python post_fit/launch_roi_plots.py gauss 8 0
python post_fit/launch_roi_plots.py gauss 8 1
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
opj = os.path.join
deb = ipdb.set_trace

# Get inputs
# ----------
fit_model = sys.argv[1]
max_tmux = int(sys.argv[2])
draw_hemi = int(sys.argv[3])

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
	if os.path.isdir(opj(base_dir,'pp_data',subject,fit_model,'h5')):
		subject_fit_list.append(subject)

# Run roi_plots.py
# -------------
num_run = 0
for subject_fit in subject_fit_list:
	if num_run < max_tmux:
		if os.path.isdir(opj(base_dir,'pp_data',subject_fit,fit_model,'figs','prf'))==0:
			session_name = "{subject_fit}_{fit_model}_post_pp_roi".format(subject_fit = subject_fit, fit_model = fit_model)
			print("tmux new-session -d -s {session_name} 'python post_fit/roi_plots.py {subject_fit} {fit_model} {draw_hemi}'".format(session_name = session_name, subject_fit = subject_fit, fit_model = fit_model, draw_hemi = draw_hemi))
			os.system("tmux new-session -d -s {session_name} 'python post_fit/roi_plots.py {subject_fit} {fit_model} {draw_hemi}'".format(session_name = session_name, subject_fit = subject_fit, fit_model = fit_model, draw_hemi = draw_hemi))
			time.sleep(2)
			num_run = num_run + 1