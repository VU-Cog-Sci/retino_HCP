"""
-----------------------------------------------------------------------------------------
pre_fit.py
-----------------------------------------------------------------------------------------
Goal of the script:
Load individual subjects, filter data per run and concatenate the files, then create the
average subject across the HCP dataset time series
-----------------------------------------------------------------------------------------
Input(s):
None
-----------------------------------------------------------------------------------------
Output(s):
Cifti files
-----------------------------------------------------------------------------------------
To run:
python pre_fit/pre_fit.py
-----------------------------------------------------------------------------------------
"""

# General imports
import sys
import numpy as np
import platform
from math import *
import os
import json
import ipdb
deb = ipdb.set_trace
opj = os.path.join
import warnings
warnings.filterwarnings('ignore')

# MRI analysis imports
import cifti
from scipy.signal import savgol_filter

# Define analysis parameters
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define cluster/server specific parameters
raw_dir = analysis_info["{platform_name}_raw_folder".format(platform_name = platform.uname()[1])]
raw_shared_dir = analysis_info["{platform_name}_raw_shared_folder".format(platform_name = platform.uname()[1])]

# Filter and normalized data per run and combine it
# print('filtering and pscing data per run -using initial blanks- and combine runs together')
# for subject in analysis_info['subject_list']:
#     print(subject)
#     for task_num,task in enumerate(analysis_info['task_list']):

#         file_in_name = opj(raw_dir,subject,'MNINonLinear','Results',"{task}".format(task = task),"{task}_Atlas_1.6mm_MSMAll_hp2000_clean.dtseries.nii".format(task = task))

#         data_load = cifti.read(file_in_name)
#         data = data_load[0]
#         run_duration = data.shape[0]
        
#         data_filt = savgol_filter(  x = data.T, 
#                                     window_length = analysis_info['sg_filt_window_length'],
#                                     polyorder = analysis_info['sg_filt_polyorder'],
#                                     deriv = analysis_info['sg_filt_deriv'],
#                                     axis = 1, 
#                                     mode = 'nearest').T

#         data_filt = data - data_filt + data_filt.mean(axis=0)
        
#         data_mean_blank = data[0:analysis_info["blanks"][task_num]].mean(axis=0)
#         data_filt_psc = 100.0 * (data_filt - data_mean_blank)/data_mean_blank
        
#         if task_num == 0:
#             all_data_filt_psc = data_filt_psc
#         else:
#             all_data_filt_psc = np.concatenate((all_data_filt_psc,data_filt_psc),axis = 0)
            
    
#     # Save the 6 tasks data
#     try: os.makedirs(opj(raw_shared_dir,subject))
#     except: pass
#     file_out_name = opj(raw_shared_dir,subject,'tfMRI_RETALL_Atlas_1.6mm_MSMAll_hp2000_clean_sg_psc.dtseries.nii')

#     series = cifti.Series(start=0, step=analysis_info["TR"]*1000, size=all_data_filt_psc.shape[0])
#     bm_full = data_load[1][1]
#     cifti.write(file_out_name, all_data_filt_psc, (series, bm_full))
        
# Average across participant _sg
print('averaging subjects data')
prep_mean_data = np.zeros((1800,170494))

for subject_num, subject in enumerate(analysis_info['subject_list']):
    print(subject,subject_num)
    file_in_name = opj(raw_shared_dir,subject,'tfMRI_RETALL_Atlas_1.6mm_MSMAll_hp2000_clean_sg_psc.dtseries.nii')
    data_load = cifti.read(file_in_name)
    data = data_load[0]
    prep_mean_data += data/len(analysis_info['subject_list'])
    

try: os.makedirs(opj(raw_shared_dir,'999999'))
except: pass
file_out_name = opj(raw_shared_dir,'999999','tfMRI_RETALL_Atlas_1.6mm_MSMAll_hp2000_clean_sg_psc.dtseries.nii')

series = cifti.Series(start=0, step=analysis_info["TR"]*1000, size=prep_mean_data.shape[0])
bm_full = data_load[1][1]
cifti.write(file_out_name, prep_mean_data, (series, bm_full))