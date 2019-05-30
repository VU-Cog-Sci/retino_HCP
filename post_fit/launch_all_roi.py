"""
-----------------------------------------------------------------------------------------
launch_pp_roi.py
-----------------------------------------------------------------------------------------
Goal of the script:
run pp_roi codes for each subjects of hcp dataset
-----------------------------------------------------------------------------------------
Input(s):
None
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
source activate i27
cd /home/szinte/projects/retino_HCP
python post_fit/launch_all_roi.py
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

opj = os.path.join
deb = ipdb.set_trace

# Define analysis parameters
# --------------------------
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Get subject list
# ----------------
subjects = analysis_info['subject_list']

# Run all codes
# -------------

for subject_num,subject in enumerate(subjects):
    print("\nSub: {} - {:1.0f}/{:1.0f}\n".format(subject,subject_num+1,len(subjects)))
    # pp_roi.py
    os.system("python post_fit/pp_roi.py {subject} gauss 2500 0".format(subject = subject))
    # post_pp_roi.py
    os.system("python post_fit/post_pp_roi.py {subject} gauss".format(subject = subject))
    # roi plots
    # os.system("python post_fit/roi_plots.py {subject} gauss 0 0".format(subject = subject))