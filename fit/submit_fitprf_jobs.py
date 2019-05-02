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
sys.argv[3]: job voxel in a job (eg. 2500)
sys.argv[4]: job duration requested in hours (used 10h on lisa)
-----------------------------------------------------------------------------------------
Output(s):
.sh file to execute in server
-----------------------------------------------------------------------------------------
Exemple:
cd /home/szinte/projects/retino_HCP/
python fit/submit_fit_jobs.py 999999 gauss 2500 10
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
subject = sys.argv[1]
fit_model = sys.argv[2]
job_vox = int(sys.argv[3])
job_dur_req = float(sys.argv[4])

# Load the analysis parameters from json file
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define server or cluster settings
if 'lisa' in platform.uname()[1]:
    jobscript_template_file = opj(os.getcwd(),'fit','lisa_jobscript_template.sh')
    base_dir = analysis_info['lisa_cluster_base_folder'] 
    sub_command = 'qsub '
    print('analysis running on lisa')
elif 'aeneas' in platform.uname()[1]:
    jobscript_template_file     =   opj(os.getcwd(),'fit','aeneas_jobscript_template.sh')
    base_dir = analysis_info['aeneas_base_folder'] 
    sub_command = 'sh '
    print('analysis running on aeneas')
elif 'local' in platform.uname()[1]:
    jobscript_template_file     =   opj(os.getcwd(),'fit','local_jobscript_template.sh')
    base_dir = analysis_info['local_base_folder'] 
    sub_command = 'sh '
    print('analysis running on local')

fit_script = 'fit/prf_fit.py'

# Create job and log output folders
try:
    os.makedirs(opj(base_dir, 'pp_data', subject, fit_model, 'jobs'))
    os.makedirs(opj(base_dir, 'pp_data', subject, fit_model, 'log_outputs'))
except:
    pass

data = []
    
# Determine data to analysed
data_file  =  opj(base_dir,'raw_data',subject,'RETBAR_ALL_tfMRI_data_sub.nii.gz')
data_file_load = nb.load(data_file)
data_file_shape = data_file_load.shape

# load or create mask
maskfn = opj(base_dir,'raw_data','RETBAR_ALL_tfMRI_data_sub_mask.nii.gz')
if os.path.isfile(maskfn)==False:
    data_file_dat = data_file_load.get_data()
    data_mask = np.zeros(data_file_shape[0:3])
    data_file_tc_std = np.std(data_file_dat,3)
    data_mask[data_file_tc_std>0] = 1

    img = nb.Nifti1Image(   dataobj = data_mask,
                            affine = data_file_load.affine,
                            header = data_file_load.header)

    img.to_filename(maskfn)
else:
    data_mask = nb.load(maskfn).get_data()

# get idx of non empty voxels
start_idx =  np.arange(0,np.sum(data_mask),job_vox)
end_idx = start_idx+job_vox
end_idx[-1] = int(np.sum(data_mask))

print('%i jobs of %1.1fh each will be run/send to %s'%(start_idx.shape[0],job_dur_req,platform.uname()[1]))

for iter_job in np.arange(0,start_idx.shape[0],1):
    

    print('input data vox num: %i to %i'%(int(start_idx[iter_job]),int(end_idx[iter_job])))

    # Define output file
    base_file_name = os.path.split(data_file[0])[-1][:-7]
    opfn = opj(base_dir,'pp_data',subject,fit_model,'fit',base_file_name + '_est_%s_to_%s.nii.gz' %(str(int(start_idx[iter_job])),str(int(end_idx[iter_job]))))


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
                '---data_file---':data_file[0],
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
        os.system('qsub ' + js_name)

    elif 'aeneas' in platform.uname()[1]:
        os.system('sh ' + js_name)
        
    elif 'local' in platform.uname()[1]:
        os.system('sh ' + js_name)
    