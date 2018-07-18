import numpy as np
from IPython import embed as shell
import re
import os
import glob
import json
import sys
import nibabel as nb
import platform

with open('../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

if 'lisa' in platform.uname()[1]:
    jobscript_template_file = os.path.join(os.getcwd(), 'lisa_jobscript_template.sh')
    base_dir = analysis_info['lisa_cluster_base_folder'] 
    sub_command = 'qsub '
    print('on lisa')
else:
    jobscript_template_file = os.path.join(os.getcwd(), 'cartesius_jobscript_template.sh')
    base_dir = analysis_info['cartesius_cluster_base_folder'] 
    sub_command = 'sbatch '
    print('on cartesius')

# only run the last, mean subject.
subject_directories = [sorted([os.path.split(fn)[1] for fn in glob.glob(os.path.join(base_dir, '*')) if os.path.isdir(fn)])[-1]]


for sd in subject_directories:
    for hemi in ["L","R"]: # gifti files are separated into different hemis

        jobscript = open(jobscript_template_file)
        working_string = jobscript.read()
        jobscript.close()

        RE_dict = {
            '---SUBJECT---':                sd,
            '---HEMI---':                   hemi,
        }

        for e in RE_dict.keys():
            working_string = working_string.replace(e, RE_dict[e])

        js_name = os.path.expanduser(os.path.join('~', 'jobs', '{sd}_{hemi}.sh'.format(sd=sd, hemi=hemi)))
        of = open(js_name, 'w')
        of.write(working_string)
        of.close()

        print('submitting ' + js_name + ' to queue')
        print(working_string)
        os.system(sub_command + js_name)
