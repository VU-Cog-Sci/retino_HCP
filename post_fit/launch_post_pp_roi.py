"""
-----------------------------------------------------------------------------------------
launch_post_pp_roi.py
-----------------------------------------------------------------------------------------
Goal of the script:
run post_pp_roi codes for each subjects of hcp dataset
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
python post_fit/launch_post_pp_roi.py gauss
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

# Check system
# ------------
sys.exit('Drawing Flatmaps only works with Python 2. Aborting.') if sys.version_info[0] > 2 else None

# Get inputs
# ----------
fit_model = sys.argv[1]

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
	if os.path.isdir(opj(base_dir,'pp_data',subject,fit_model,'figs','roi')):
		subject_fit_list.append(subject)

# Run post_pp_roi.py
# ------------------
for subject_fit in subject_fit_list:
	if os.path.isdir(opj(base_dir,'pp_data',subject_fit,fit_model,'h5'))==0:
		print("python post_fit/post_pp_roi.py {subject_fit} {fit_model}".format(subject_fit = subject_fit, fit_model = fit_model))
		os.system("python post_fit/post_pp_roi.py {subject_fit} {fit_model}".format(subject_fit = subject_fit, fit_model = fit_model))
		time.sleep(2)
		
