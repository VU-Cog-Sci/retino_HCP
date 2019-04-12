"""
-----------------------------------------------------------------------------------------
submit_makegrid_jobs
-----------------------------------------------------------------------------------------
Goal of the script:
create jobscript to run locally, in a cluster (LISA) or server (AENEAS)
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[2]: fit model ('gauss','css')
sys.argv[4]: job duration requested in hours (used 20h on lisa)
-----------------------------------------------------------------------------------------
Output(s):
.sh file to execute in server
-----------------------------------------------------------------------------------------
Exemple:
cd /home/szinte/projects/retino_HCP/
python fit/submit_makegrid_jobs.py gauss_sg 20
-----------------------------------------------------------------------------------------
"""

# General imports
import numpy as np
import os
import glob
import json
import sys
import nibabel as nb
import platform
import ipdb
deb = ipdb.set_trace
opj = os.path.join

# Get subject number and hemisphere to analyse
fit_model = sys.argv[1]
job_dur_req = float(sys.argv[2])

# Load the analysis parameters from json file
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define server or cluster settings
if 'lisa' in platform.uname()[1]:
    jobscript_template_file = opj(os.getcwd(),'fit','lisa_makegrid_jobscript_template.sh')
    base_dir = analysis_info['lisa_cluster_base_folder'] 
    sub_command = 'sbatch '
    print('analysis running on lisa')
elif 'aeneas' in platform.uname()[1]:
    jobscript_template_file     =   opj(os.getcwd(),'fit','aeneas_makegrid_jobscript_template.sh')
    base_dir = analysis_info['aeneas_base_folder'] 
    sub_command = 'sh '
    print('analysis running on aeneas')

fun_script = 'fit/prf_makegrid.py'

# Create job and log output folders
try:
    os.makedirs(opj(base_dir, 'pp_data', 'makegrid', fit_model, 'jobs'))
    os.makedirs(opj(base_dir, 'pp_data', 'makegrid', fit_model, 'log_outputs'))
except:
    pass
    
job_input = []

# Define output file
opfn = opj(base_dir,'pp_data','makegrid',fit_model,'grid_predictions.hdf5')
if os.path.isfile(opfn):
    if os.path.getsize(opfn) != 0:
        sys.exit('%s already exists and is non-empty. \naborting analysis \ndelete manually to run new grid'%opfn)

# create job shell
jobscript = open(jobscript_template_file)
working_string = jobscript.read()
jobscript.close()
job_dur = '%i:00:00'%job_dur_req
    
re_dict = { '---job_dur---':job_dur,
            '---fun_script---':fun_script,
            '---fit_model---':fit_model,
            '---save_file---':opfn,
            '---base_dir---':base_dir}

for e in re_dict.keys():
    working_string  =   working_string.replace(e, re_dict[e])

js_name =  opj(base_dir, 'pp_data', 'makegrid', fit_model, 'jobs', 'makegrid.sh')
of = open(js_name, 'w')
of.write(working_string)
of.close()

# Submit jobs
print('submitting ' + js_name + ' to queue')

if 'lisa' in platform.uname()[1]:
    os.chdir(opj(base_dir,'pp_data','makegrid',fit_model,'log_outputs'))
    os.system('sbatch ' + js_name)

elif 'aeneas' in platform.uname()[1]:
    os.system('sh ' + js_name)