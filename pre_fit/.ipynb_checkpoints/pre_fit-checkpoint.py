"""
-----------------------------------------------------------------------------------------
pre_fit.py
-----------------------------------------------------------------------------------------
Goal of the script:
SG filter, PSC, AVG runs and combine data of both hemisphere
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject name
sys.argv[2]: start voxel index
sys.argv[3]: end voxel index
sys.argv[4]: data file path
sys.argv[5]: main directory
-----------------------------------------------------------------------------------------
Output(s):
# Preprocessed timeseries files
-----------------------------------------------------------------------------------------
To run:
cd /home/szinte/projects/retino_HCP
python pre_fit/pre_fit.py
-----------------------------------------------------------------------------------------
"""

# Stop warnings
# -------------
import warnings
warnings.filterwarnings("ignore")

# General imports
# ---------------
import json
import sys
import os
import glob
import ipdb
import platform
import numpy as np
opj = os.path.join
deb = ipdb.set_trace

# MRI analysis imports
# --------------------
import nibabel as nb
from scipy.signal import savgol_filter
from nipype.interfaces.freesurfer import SurfaceTransform

with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

with open('select_block.json') as f:
    json_s = f.read()
    select_block = json.loads(json_s)
    
trans_cmd = 'rsync -avuz --progress'

# Define cluster/server specific parameters
# -----------------------------------------
if 'aeneas' in platform.uname()[1]:
    base_dir = analysis_info['aeneas_base_folder'] 
elif 'lisa' in platform.uname()[1]:
    base_dir = analysis_info['lisa_base_folder'] 
    
# Copy files in raw_data folder
# ----------------------------
for sub_name in analysis_info['subject_list'] :
    sub_session = select_block["{sub}_sessions".format(sub = sub_name)]
    dest_folder = "{base_dir}/raw_data/{sub}".format(base_dir = base_dir, sub = sub_name)
    try: os.makedirs(dest_folder)
    except: pass
      
    for session in sub_session:
        sub_session_run = select_block["{sub}_{ses}_run".format(sub = sub_name, ses = session)]

        for run in sub_session_run:
            orig_folder = "{base_dir}/derivatives/fmriprep/{sub}/{ses}/func".format(base_dir = base_dir, sub = sub_name, ses=session)
            
            for hemi in ['L','R']:
                orig_file = "{orig_fold}/{sub}_{ses}_task-prf_{run}_space-fsaverage6_hemi-{hemi}.func.gii".format(orig_fold = orig_folder, 
                                                                                                                  sub = sub_name, ses=session, 
                                                                                                                  run = run, hemi = hemi)
                dest_file = "{dest_fold}/{sub}_{ses}_task-prf_{run}_space-fsaverage6_hemi-{hemi}.func.gii".format(dest_fold = dest_folder, 
                                                                                                                  sub = sub_name, ses=session, 
                                                                                                                  run = run, hemi = hemi)

                os.system("{cmd} {orig} {dest}".format(cmd = trans_cmd, orig = orig_file, dest = dest_file))
                
                
# SG + PSC + AVG + COMBINE HEMI
# -----------------------------
# sxfm = SurfaceTransform()
# sxfm.inputs.source_subject = "fsaverage6"
# sxfm.inputs.target_subject = "fsaverage"
# sxfm.terminal_output = 'none'
# sxfm.inputs.subjects_dir = opj(base_dir,'derivatives','freesurfer')

for sub_name in analysis_info['subject_list'] :
    
    # SG + PSC
    # --------
    print(sub_name+': sg + psc')
    file_list = sorted(glob.glob("{base_dir}/raw_data/{sub}/*func.gii".format(base_dir = base_dir, sub = sub_name)))
    
    for file in file_list:
        
        
        print(file)
        # load
        pp_hemi = []
        pp_hemi_file = nb.load(file)
        pp_hemi.append(np.array([pp_hemi_file.darrays[i].data for i in range(len(pp_hemi_file.darrays))]))
        pp_hemi = np.vstack(pp_hemi)

        # sg filter
        pp_hemi_filt = savgol_filter( x = pp_hemi.T,
                                      window_length = analysis_info['sg_filt_window_length'],
                                      polyorder = analysis_info['sg_filt_polyorder'],
                                      deriv = analysis_info['sg_filt_deriv'],
                                      axis = 1, 
                                      mode = 'nearest').T

        pp_hemi_sg = pp_hemi - pp_hemi_filt + pp_hemi_filt.mean(axis=0)

        # percent signal change
        pp_hemi_sg_median = np.median(pp_hemi_sg, axis=0)
        pp_hemi_sg_psc = 100.0 * (pp_hemi_sg - pp_hemi_sg_median)/pp_hemi_sg_median

        # save
        gii_out_name = file[:-4] + '_sg_psc.gii'
        darrays_pp_hemi_sg_psc = [nb.gifti.gifti.GiftiDataArray(d) for d in pp_hemi_sg_psc]
        gii_out = nb.gifti.gifti.GiftiImage(header = pp_hemi_file.header, 
                                            extra = pp_hemi_file.extra,
                                            darrays = darrays_pp_hemi_sg_psc)

        nb.save(img = gii_out,filename = gii_out_name)
        
#         # convert to fsaverage
#         if 'hemi-L' in gii_out_name: 
#             sxfm.inputs.hemi = "lh"
#             hemi = 'L'
#         elif 'hemi-R' in gii_out_name: 
#             sxfm.inputs.hemi = "rh"
#             hemi = 'R'

#         sxfm.inputs.source_file = gii_out_name
#         sxfm.inputs.out_file = gii_out_name[:-33]+"fsaverage_hemi-{hemi}.func_sg_psc.gii".format(hemi = hemi)
#         print(sxfm.inputs.out_file)
#         sxfm.run()
    
    
    # AVERAGE RUNS
    # ------------
    print(sub_name+': average runs')
    for hemi in ['L','R']:
        file_list = sorted(glob.glob("{base_dir}/raw_data/{sub}/*fsaverage6_hemi-{hemi}.func_sg_psc.gii".format(base_dir = base_dir, sub = sub_name, hemi = hemi)))

        pp_hemi_sg_psc_avg = np.zeros((120,40962))
        for file in file_list:
            print(file)
            # load
            pp_hemi_sg_psc = []
            pp_hemi_sg_psc_file = nb.load(file)
            pp_hemi_sg_psc.append(np.array([pp_hemi_sg_psc_file.darrays[i].data for i in range(len(pp_hemi_sg_psc_file.darrays))]))
            pp_hemi_sg_psc = np.vstack(pp_hemi_sg_psc)

            
            # avg
            pp_hemi_sg_psc_avg += pp_hemi_sg_psc/len(file_list)
            
        # save
        gii_out_name = "{base_dir}/raw_data/{sub}/{sub}_task-prf_space-fsaverage6_hemi-{hemi}.func_sg_psc.gii".format(base_dir = base_dir, 
                                                                                                                     sub = sub_name, 
                                                                                                                     hemi = hemi)

        darrays_pp_hemi_sg_psc_avg = [nb.gifti.gifti.GiftiDataArray(d) for d in pp_hemi_sg_psc_avg]
        gii_out = nb.gifti.gifti.GiftiImage(header = pp_hemi_file.header,
                                            extra = pp_hemi_file.extra,
                                            darrays = darrays_pp_hemi_sg_psc_avg)

        nb.save(img = gii_out,filename = gii_out_name)
        
#         # convert to fsaverage
#         if 'hemi-L' in gii_out_name: 
#             sxfm.inputs.hemi = "lh"
#             hemi = 'L'
#         elif 'hemi-R' in gii_out_name: 
#             sxfm.inputs.hemi = "rh"
#             hemi = 'R'

#         sxfm.inputs.source_file = gii_out_name
#         sxfm.inputs.out_file = gii_out_name[:-33]+"fsaverage_hemi-{hemi}.func_sg_psc.gii".format(hemi = hemi)
#         print(sxfm.inputs.out_file)
#         sxfm.run()
        
        
    # COMBINE HEMISPHERE
    # ------------------
    print(sub_name+': combine hemisphere')
    
#     sub_session = select_block["{sub}_sessions".format(sub = sub_name)]
#     dest_folder = "{base_dir}/raw_data/{sub}".format(base_dir = base_dir, sub = sub_name)
    
#     for session in sub_session:
#         sub_session_run = select_block["{sub}_{ses}_run".format(sub = sub_name, ses = session)]

#         for run in sub_session_run:
#             hemi_mat = []
#             for hemi in ['L','R']:
#                 hemi_filename = "{base_dir}/raw_data/{sub}/{sub}_{ses}_task-prf_{run}_space-fsaverage_hemi-{hemi}.func_sg_psc.gii".format(base_dir = base_dir, 
#                                                                                                         sub = sub_name, ses=session, hemi = hemi,
#                                                                                                         run = run)
#                 hemi_file = nb.load(hemi_filename)
#                 hemi_mat.append(np.array([hemi_file.darrays[i].data for i in range(len(hemi_file.darrays))]))
            
#             hemi_mat = np.hstack(hemi_mat)
#             # save
#             gii_out_name = "{base_dir}/raw_data/{sub}/{sub}_{ses}_task-prf_{run}_space-fsaverage.func_sg_psc.gii".format(base_dir = base_dir, 
#                                                                                                      sub = sub_name,ses=session, 
#                                                                                                         run = run)
#             darrays_hemi_mat = [nb.gifti.gifti.GiftiDataArray(d) for d in hemi_mat]
#             gii_out = nb.gifti.gifti.GiftiImage(header = hemi_file.header,
#                                         extra = hemi_file.extra,
#                                         darrays = darrays_hemi_mat)
    
#             nb.save(img = gii_out,filename = gii_out_name)
            
    
    hemi_mat=[]
    for hemi in ['L','R']:
        # load
        hemi_filename = "{base_dir}/raw_data/{sub}/{sub}_task-prf_space-fsaverage6_hemi-{hemi}.func_sg_psc.gii".format(base_dir = base_dir, 
                                                                                                                      sub = sub_name, 
                                                                                                                      hemi = hemi)
        hemi_file = nb.load(hemi_filename)

        # load
        hemi_mat.append(np.array([hemi_file.darrays[i].data for i in range(len(hemi_file.darrays))]))

    hemi_mat = np.hstack(hemi_mat)
    
    # save
    gii_out_name = "{base_dir}/raw_data/{sub}/{sub}_task-prf_space-fsaverage6.func_sg_psc.gii".format(base_dir = base_dir, 
                                                                                                     sub = sub_name)
    darrays_hemi_mat = [nb.gifti.gifti.GiftiDataArray(d) for d in hemi_mat]
    gii_out = nb.gifti.gifti.GiftiImage(header = hemi_file.header,
                                        extra = hemi_file.extra,
                                        darrays = darrays_hemi_mat)
    
    nb.save(img = gii_out,filename = gii_out_name)
    
    
    