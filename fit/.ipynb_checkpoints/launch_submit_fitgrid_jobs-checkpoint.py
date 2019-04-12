"""
-----------------------------------------------------------------------------------------
launch_submit_fitgrid_jobs.py
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
cd /home/szinte/projects/retino_HCP
python fit/launch_submit_fitgrid_jobs.py gauss_sg
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

# Define analysis parameters
# --------------------------
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

fit_model = sys.argv[1]
num_vox = 5
dur_fit = 10

subs = analysis_info['subject_list']

# Subject to analyse
# ------------------
index_start = 0
index_end = 2
for subject in subs[index_start:index_end]:
    print("python fit/submit_fitgrid_jobs.py {subject} {fit_model} {num_vox} {dur_fit}".format(subject = subject, 
                                                                                               fit_model = fit_model, 
                                                                                               num_vox = num_vox, 
                                                                                               dur_fit = dur_fit))
    
    os.system("python fit/submit_fitgrid_jobs.py {subject} {fit_model} {num_vox} {dur_fit}".format(subject = subject, 
                                                                                               fit_model = fit_model, 
                                                                                               num_vox = num_vox, 
                                                                                               dur_fit = dur_fit))