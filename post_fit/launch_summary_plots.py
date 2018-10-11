"""
-----------------------------------------------------------------------------------------
launch_summary_plots.py
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
python post_fit/launch_summary_plots.py gauss
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
analysed_subject = []

for subject in subject_list:
	if os.path.isdir(opj(base_dir,'pp_data',subject,fit_model,'figs','roi')):
		analysed_subject.append(subject)

analysed_subject.append('000000')
analysed_subject.append('999999')
analysed_subject.append('999999_ns13')
analysed_subject.append('999999_mm16')
analysed_subject.append('999999_mm16_ns13')

# Run post_pp_roi.py
# ------------------
for subject_fit in analysed_subject:	
	
	if subject == '000000':
		for eb in ['std','sem','95ci']:
			# print("python post_fit/summary_plots.py {subject_fit} {fit_model} {eb}".format(subject_fit = subject_fit, fit_model = fit_model, eb = eb))
			os.system("python post_fit/summary_plots.py {subject_fit} {fit_model} {eb} 0".format(subject_fit = subject_fit, fit_model = fit_model, eb = eb))
	else:
		# print("python post_fit/summary_plots.py {subject_fit} {fit_model} std".format(subject_fit = subject_fit, fit_model = fit_model))
		os.system("python post_fit/summary_plots.py {subject_fit} {fit_model} std 0".format(subject_fit = subject_fit, fit_model = fit_model))
	
	
