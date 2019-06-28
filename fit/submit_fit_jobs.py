"""
-----------------------------------------------------------------------------------------
submit_fit_jobs
-----------------------------------------------------------------------------------------
Goal of the script:
create jobscript to run locally, in a cluster (LISA) or server (AENEAS)
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject name (e.g. 'sub-001')
sys.argv[2]: fit model ('gauss','css')
sys.argv[3]: voxel per jobs (used 400 on lisa)
sys.argv[4]: job duration requested in hours (used 10h on lisa)
-----------------------------------------------------------------------------------------
Output(s):
.sh file to execute in server
-----------------------------------------------------------------------------------------
Exemple:
#1: sub 192641 : ran 18/08/2018 at 20:15
python fit/submit_fit_jobs.py 999999 gauss 5000 20
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
import cifti
deb = ipdb.set_trace
opj = os.path.join

# Get subject number and hemisphere to analyse
subject = sys.argv[1]
fit_model = sys.argv[2]
job_vox = float(sys.argv[3])
job_dur_req = float(sys.argv[4])

# Load the analysis parameters from json file
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define server or cluster settings
if 'lisa' in platform.uname()[1]:
    jobscript_template_file = opj(os.getcwd(),'fit','lisa_jobscript_template.sh')
    base_dir = analysis_info['lisa_base_folder'] 
    sub_command = 'sbatch '
    print('analysis running on lisa')
elif 'aeneas' in platform.uname()[1]:
    jobscript_template_file     =   opj(os.getcwd(),'fit','aeneas_jobscript_template.sh')
    base_dir = analysis_info['aeneas_base_folder'] 
    sub_command = 'sh '
    print('analysis running on aeneas')
else: # cartesius
    jobscript_template_file = opj(os.getcwd(),'fit','cartesius_jobscript_template.sh')
    base_dir = analysis_info['cartesius_base_folder'] 
    sub_command = 'sbatch '
    print('analysis running on cartesius')


fit_script = 'fit/prf_fit.py'

# Create job and log output folders
try:
    os.makedirs(opj(base_dir, 'pp_data', subject, fit_model, 'jobs'))
    os.makedirs(opj(base_dir, 'pp_data', subject, fit_model, 'log_outputs'))
except:
    pass

data = []
    
# Determine data to analyse
data_file  =  "{basedir}/raw_data/{sub}/tfMRI_RETBARS_Atlas_1.6mm_MSMAll_hp2000_clean_sg_psc.dtseries.nii".format(basedir = base_dir,
                                                                                                                sub = subject)

# Cut it in small pieces of voxels
data_file_load = cifti.read(data_file)
data = data_file_load[0]
data_size = data.shape

start_idx =  np.arange(0,data_size[1],job_vox)
end_idx = start_idx+job_vox
end_idx[-1] = data_size[1]


print('%i jobs of %1.1fh each will be run/send to %s'%(start_idx.shape[0],job_dur_req,platform.uname()[1]))

job_input = []
for iter_job in np.arange(0,start_idx.shape[0],1):
    job_input = data[:,int(start_idx[iter_job]):int(end_idx[iter_job])]

    print('input data vox num: %i to %i'%(int(start_idx[iter_job]),int(end_idx[iter_job])))

    # Define output file
    base_file_name = os.path.split(data_file)[-1][:-13]
    opfn = opj(base_dir,'pp_data',subject,fit_model,'fit', "{base_file_name}_est_{start}_to_{end}.dtseries.nii".format(  base_file_name = base_file_name,
                                                                                                                         start = int(start_idx[iter_job]),
                                                                                                                         end = int(end_idx[iter_job])))
    
    if os.path.isfile(opfn):
        if os.path.getsize(opfn) != 0:
            print('output file %s already exists and is non-empty. aborting analysis of voxels %s to %s'%(opfn,str(int(start_idx[iter_job])),str(int(end_idx[iter_job]))))
            continue

    # create job shell
    jobscript = open(jobscript_template_file)
    working_string = jobscript.read()
    jobscript.close()
    job_dur = '%i:00:00'%job_dur_req
    
    re_dict = { '---job_dur---':job_dur,
                '---fit_file---':fit_script,
                '---fit_model---':fit_model,
                '---subject---':subject,
                '---start_idx---':str(int(start_idx[iter_job])),
                '---end_idx---':str(int(end_idx[iter_job])),
                '---data_file---':data_file,
                '---base_dir---':base_dir}

    for e in re_dict.keys():
        working_string  =   working_string.replace(e, re_dict[e])

    js_name =  opj(base_dir, 'pp_data', subject, fit_model, 'jobs', '%s_vox_%s_to_%s.sh'%(subject,str(int(start_idx[iter_job])),str(int(end_idx[iter_job]))))

    of = open(js_name, 'w')
    of.write(working_string)
    of.close()

    # Submit jobs
    print('submitting ' + js_name + ' to queue')

    if 'lisa' in platform.uname()[1]:
        os.chdir(opj(base_dir,'pp_data',subject,fit_model,'log_outputs'))
        os.system(sub_command + js_name)

    elif 'aeneas' in platform.uname()[1]:
        os.system(sub_command + js_name)
 
    else: # cartesius
        os.chdir(opj(base_dir,'pp_data',subject,fit_model,'log_outputs'))
        os.system(sub_command + js_name)
    