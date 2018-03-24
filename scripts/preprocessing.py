import os
import glob
import json
from joblib import Parallel, delayed

from ..retino_HCP.utils import *

with open('../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

base_dir = analysis_info['cluster_base_folder'] 

subject_directories = [fn for fn in glob.glob(os.path.join(base_dir, '*')) if os.path.isdir(fn)]

for sd in subject_directories:
    # get ordered files
    rawest_cii_files = [glob.glob(os.path.join(sd, '*%s*dtseries.nii'%run))[0] for run in analysis_info["run_order"]]
    # temporal preprocessing

    pp_out_files = Parallel(n_jobs=6)(delayed(sg_psc_cii)(cii_file) for cii_file in rawest_cii_files)

    # each pair averaged:
    av_wedge_file = average_phase_encoded_ciis(pp_out_files[0],pp_out_files[1])
    av_ring_file = average_phase_encoded_ciis(pp_out_files[2],pp_out_files[3])
    av_bar_file = average_bar_ciis(pp_out_files[4],pp_out_files[5])





#################################################################################

with open('../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

base_dir = analysis_info['cluster_base_folder'] 

subject_directories = [fn for fn in glob.glob(os.path.join(base_dir, '*')) if os.path.isdir(fn)]

for sd in subject_directories:
    os.chdir(sd)
    cii_files = glob.glob('*.dseries.nii')

    wbc = """wb_command -cifti-separate {cii} \
        COLUMN -volume-all data_sub.nii \
        -metric CORTEX_LEFT {cii_n}_L.func.gii \
        -metric CORTEX_RIGHT {cii_n}_R.func.gii"""

    for cii in cii_files:
        wbc_c = wbc.format(cii=cii, cii_n=cii[-4])

        print(wbc_c)
        # os.system(wbc_c)