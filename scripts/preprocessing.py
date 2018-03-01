import os
import glob
import json
from joblib import Parallel, delayed

from ..retino_HCP.utils import *

with open('../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

base_dir = '/home/shared/2018/visual/HCP7TFIXED'

subject_directories = [fn for fn in glob.glob(os.path.join(base_dir, '*')) if os.path.isdir(fn)]

for sd in subject_directories:
    # get ordered files
    rawest_cii_files = [glob.glob(os.path.join(sd, '*%s*dseries.nii'%run))[0] for run in analysis_info["run_order"]]
    # temporal preprocessing

    pp_out_files = Parallel(n_jobs=6)(delayed(sg_psc_cii)(cii_file) for cii_file in rawest_cii_files)

    pp_pairs = [[pp_out_files[0],pp_out_files[1]],
                [pp_out_files[2],pp_out_files[3]]
                [pp_out_files[4],pp_out_files[5]]]

    # first & second pair:
    av_wedge_file = average_phase_encoded_ciis(pp_out_files[0],pp_out_files[1])
    av_ring_file = average_phase_encoded_ciis(pp_out_files[2],pp_out_files[3])
    av_bar_file = average_bar_ciis(pp_out_files[4],pp_out_files[5])


