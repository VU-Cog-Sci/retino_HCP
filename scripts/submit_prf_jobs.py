import numpy as np
from IPython import embed as shell
import re
import os
import glob
import json
import sys
import nibabel as nb

jobscript_template_file = os.path.join(os.getcwd(), 'jobscript_template.sh')

with open('../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

base_dir = analysis_info['cluster_base_folder'] 

subject_directories = [os.path.split(fn)[1] for fn in glob.glob(os.path.join(base_dir, '*')) if os.path.isdir(fn)]


for sd in subject_directories[:10]:
    jobscript = open(jobscript_template_file)
    working_string = jobscript.read()
    jobscript.close()

    RE_dict = {
        '---SUBJECT---':                sd
    }

    for e in RE_dict.keys():
        working_string = working_string.replace(e, RE_dict[e])

    js_name = os.path.expanduser(os.path.join('~', 'jobs', sd + '.sh'))
    of = open(js_name, 'w')
    of.write(working_string)
    of.close()

    print('submitting ' + js_name + ' to queue')
    print(working_string)
    os.system('sbatch ' + js_name)
