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
    for hemi in ["L","R"]: # gifti files are separated into different hemis
        # get ordered files
        rawest_gii_files = [glob.glob(os.path.join(sd, '*%s*_%s.func.gii'%(run, hemi)))[0] for run in analysis_info["run_order"]]
        # temporal preprocessing

        pp_out_files = Parallel(n_jobs=6)(delayed(sg_psc_gii)(gii_file) for gii_file in rawest_gii_files)

        # each pair averaged:
        av_wedge_file = average_phase_encoded_giis(pp_out_files[0],pp_out_files[1])
        av_ring_file = average_phase_encoded_giis(pp_out_files[2],pp_out_files[3])
        av_bar_file = average_bar_giis(pp_out_files[4],pp_out_files[5])





#################################################################################
#####
#####   conversion to GII
#####
#################################################################################

with open('../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

base_dir = analysis_info['cluster_base_folder'] 

subject_directories = sorted([fn for fn in glob.glob(os.path.join(base_dir, '*')) if os.path.isdir(fn)])

for sd in subject_directories[:10]:
    cii_files = glob.glob(os.path.join(sd,'*av_est.nii'))

    wbc = """wb_command -cifti-separate {cii} \
 COLUMN -volume-all {cii_n}_data_sub.nii \
 -metric CORTEX_LEFT {cii_n}_L.func.gii \
 -metric CORTEX_RIGHT {cii_n}_R.func.gii &"""

    for cii in cii_files:
        wbc_c = wbc.format(cii=cii, cii_n=cii[:-4])

        # if cii == cii_files[-1]:
        #     wbc_c = wbc_c[:-1]
        print(wbc_c)
        os.system(wbc_c)