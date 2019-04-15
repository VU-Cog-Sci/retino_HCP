"""
-----------------------------------------------------------------------------------------
pp_roi.py
-----------------------------------------------------------------------------------------
Goal of the script:
Region of interests pre-processing
Compute pRF derivatives and plot on pycortex overlay.svg to determine visual ROI
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject number
sys.argv[2]: fit model ('gauss','css')
sys.argv[3]: voxels per fit (e.g 2500)
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
cd /home/szinte/projects/retino_HCP
python post_fit/pp_roi.py sub-01 gauss_sg 2500
python post_fit/pp_roi.py sub-02 gauss_sg 2500
-----------------------------------------------------------------------------------------
"""

# Stop warnings
# -------------
import warnings
warnings.filterwarnings("ignore")

# General imports
# ---------------
import os
import sys
import json
import glob
import numpy as np
import matplotlib.pyplot as pl
import ipdb
import platform
opj = os.path.join
deb = ipdb.set_trace

# MRI imports
# -----------
import nibabel as nb
import cortex
from nipype.interfaces.freesurfer import SurfaceTransform

# Functions import
# ----------------
from utils import set_pycortex_config_file, convert_fit_results, draw_cortex_vertex

# Check system
# ------------
sys.exit('Drawing Flatmaps only works with Python 2. Aborting.') if sys.version_info[0] > 2 else None

# Get inputs
# ----------
subject = sys.argv[1]
fit_model = sys.argv[2]
job_vox = float(sys.argv[3])
if fit_model == 'gauss' or fit_model == 'gauss_sg': fit_val = 6
elif fit_model == 'css' or fit_model == 'css_sg': fit_val = 7
base_file_name = "{sub}_task-prf_space-fsaverage6.func_sg_psc".format(sub = subject)

# Define analysis parameters
# --------------------------
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define cluster/server specific parameters
# -----------------------------------------
if 'aeneas' in platform.uname()[1]:
    base_dir = analysis_info['aeneas_base_folder'] 
elif 'lisa' in platform.uname()[1]:
    base_dir = analysis_info['lisa_base_folder'] 
    
# Change cortex database folder
# -----------------------------
pycortex_folder     =   opj(base_dir,'pp_data','cortex')
set_pycortex_config_file(project_folder = pycortex_folder)

# Conversion settings
# -------------------
sxfm = SurfaceTransform()
sxfm.inputs.source_subject = "fsaverage6"
sxfm.inputs.target_subject = "fsaverage"
sxfm.terminal_output = 'none'
sxfm.inputs.subjects_dir = opj(base_dir,'derivatives','freesurfer')

# Determine number of vertex and time_serie
# -----------------------------------------
data = []
data_file  =  sorted(glob.glob(opj(base_dir,'raw_data',subject,"{bfn}.gii".format(bfn = base_file_name))))
data_file_load = nb.load(data_file[0])
data.append(np.array([data_file_load.darrays[i].data for i in range(len(data_file_load.darrays))]))
data = np.vstack(data) 
ts_num,vox_num = data.shape[0],data.shape[1]

for type_data in ["gridfit","fit"]:
    
    # Determine derivative folder
    # ---------------------------
    deriv_dir = opj(base_dir,'pp_data',subject,fit_model,'deriv',type_data)
    
    # Check if all slices are present
    # -------------------------------
    start_idx =  np.arange(0,vox_num,job_vox)
    end_idx = start_idx+job_vox
    end_idx[-1] = vox_num
    num_miss_part = 0
    fit_est_files_job = []
    for iter_job in np.arange(0,start_idx.shape[0],1):
        fit_est_file = opj(base_dir,'pp_data',subject,fit_model,type_data,
                           "{bfn}_est_{start_idx}_to_{end_idx}.gii".format(  bfn = base_file_name,
                                                                             start_idx = str(int(start_idx[iter_job])),
                                                                             end_idx = str(int(end_idx[iter_job]))))
        if os.path.isfile(fit_est_file):
            if os.path.getsize(fit_est_file) == 0:  num_miss_part += 1 
            else: fit_est_files_job.append(fit_est_file)
        else:
            num_miss_part += 1

    if num_miss_part != 0:
        sys.exit('%i missing files of %s, analysis stopped'%(num_miss_part,type_data))

    # Combine fit files
    # -----------------
    print('%s: combining fit files'%type_data)

    data = np.zeros((fit_val,vox_num))
    idx_start = 0
    idx_end = 0
    for fit_filename in fit_est_files_job:
        data_fit = []
        data_fit_file = nb.load(fit_filename)
        data_fit.append(np.array([data_fit_file.darrays[i].data for i in range(len(data_fit_file.darrays))]))
        data_fit = np.vstack(data_fit)
        idx_end += data_fit.shape[0] 
        data[:,idx_start:idx_end] = data_fit.T
        idx_start += data_fit.shape[0]

    # Seperate hemi files
    # -------------------
    data_L = data[:,0:vox_num/2]
    data_R = data[:,vox_num/2:vox_num]

    for hemi in ['L','R']:
        exec("darrays_est_{hemi} = [nb.gifti.gifti.GiftiDataArray(d) for d in data_{hemi}]".format(hemi = hemi))
        exec("gii_out_{hemi} = nb.gifti.gifti.GiftiImage(header = data_fit_file.header, extra = data_fit_file.extra, darrays = darrays_est_{hemi})".format(hemi = hemi))
        exec("prf_filename_{hemi} = opj(base_dir,'pp_data',subject,fit_model,'{type_data}','{bfn}_est_{hemi}.gii')".format(bfn =base_file_name, hemi = hemi, type_data = type_data))
        exec("nb.save(gii_out_{hemi}, prf_filename_{hemi})".format(hemi = hemi))

    # Compute derived measures from prfs
    # ----------------------------------
    print('%s: extracting pRF derivatives'%type_data)
    for hemi in ['L','R']:
        exec("prf_filename = prf_filename_{hemi}".format(hemi = hemi)) 
        convert_fit_results(prf_filename = prf_filename,
                            output_dir = deriv_dir,
                            hemi = hemi,
                            stim_width = analysis_info['stim_width'],
                            stim_height = analysis_info['stim_height'],
                            fit_model = fit_model)

    # Resample gii to fsaverage
    # -------------------------
    print('%s: converting derivative files to fsaverage'%type_data)

    for hemi in ['L','R']:        
        if hemi == 'L': sxfm.inputs.hemi = "lh"
        elif hemi == 'R': sxfm.inputs.hemi = "rh"
            
        for mask_dir in ['all','pos','neg']:
            sxfm.inputs.source_file = opj(deriv_dir,mask_dir,"prf_deriv_{hemi}_{mask_dir}.gii".format(hemi = hemi, mask_dir = mask_dir))
            sxfm.inputs.out_file = opj(deriv_dir,mask_dir,"prf_deriv_{hemi}_{mask_dir}_fsaverage.gii".format(hemi = hemi, mask_dir = mask_dir))
            print(sxfm.inputs.out_file)
            sxfm.run()

    # Create derivatives flatmaps
    # ---------------------------
    print('%s: draw deriv maps'%type_data)
    cmap_neg_pos = 'RdBu_r'
    cmap_polar = 'hsv'
    cmap_gain = 'viridis'
    col_offset = 1/14.0
    cmap_pos = 'Reds'
    polar_col_steps = [4.0, 8.0, 16.0, 255.0]
    cmap_ecc_size = 'Spectral'
    sign_idx, rsq_idx, ecc_idx, polar_real_idx, polar_imag_idx , size_idx, \
                non_lin_idx, amp_idx, baseline_idx, cov_idx, x_idx, y_idx = 0,1,2,3,4,5,6,7,8,9,10,11

    for mask_dir in ['all','pos','neg']:

        # Create figure folders
        vertex_names = []
        all_vertex   = []
        exec('fig_roi_dir_{mask_dir} = opj(base_dir,"pp_data",subject,fit_model,"figs","{type_data}","roi","{mask_dir}")'.format(mask_dir=mask_dir,type_data = type_data))
        try: exec('os.makedirs(fig_roi_dir_{mask_dir})'.format(mask_dir=mask_dir))
        except: pass

        # Combine hemispheres
        deriv_mat=[]
        for hemi in ['L','R']:
            deriv_file = nb.load(opj(deriv_dir,mask_dir,"prf_deriv_{hemi}_{mask_dir}_fsaverage.gii".format(hemi = hemi, mask_dir = mask_dir)))
            deriv_mat.append(np.array([deriv_file.darrays[i].data for i in range(len(deriv_file.darrays))]))
        deriv_mat = np.hstack(deriv_mat)

        # R-square
        rsq_data = deriv_mat[rsq_idx,:]
        alpha = rsq_data
        param_rsq = {'subject': 'fsaverage', 'data': rsq_data.T, 'cmap': cmap_pos, 'alpha': alpha.T, 'vmin': 0,'vmax': 1,'cbar': 'discrete'}
        vertex_names.append('rsq')

        # Polar angle
        pol_comp_num = deriv_mat[polar_real_idx,:] + 1j * deriv_mat[polar_imag_idx,:]
        polar_ang = np.angle(pol_comp_num)
        ang_norm = (polar_ang + np.pi) / (np.pi * 2.0)

        for cmap_steps in polar_col_steps:
            param_polar = {'data': ang_norm.T, 'cmap': cmap_polar, 'alpha': alpha.T, 'vmin': 0, 'vmax': 1, 'cmap_steps': cmap_steps,\
                           'curv_brightness': 0.05, 'curv_contrast': 0.1, 'cbar': 'polar', 'col_offset': col_offset}
            exec('param_polar_{csteps} = param_polar'.format(csteps = int(cmap_steps)))
            exec('vertex_names.append("polar_{csteps}")'.format(csteps = int(cmap_steps)))

        # Eccentricity
        ecc_data = deriv_mat[ecc_idx,:]
        param_ecc = {'data': ecc_data.T, 'cmap': cmap_ecc_size, 'alpha': alpha.T, 'vmin': 0, 'vmax': 8,'cbar': 'ecc'}
        vertex_names.append('ecc')

        # Sign
        sign_data = deriv_mat[sign_idx,:]
        param_sign = {'data': sign_data.T, 'cmap': cmap_neg_pos, 'alpha': alpha.T, 'vmin': -1, 'vmax': 1, 'cbar': 'discrete'}
        vertex_names.append('sign')

        # Size
        size_data = deriv_mat[size_idx,:]
        param_size = {'data': size_data.T, 'cmap': cmap_ecc_size, 'alpha': alpha.T, 'vmin': 0, 'vmax': 8, 'cbar': 'discrete'}
        vertex_names.append('size')

        # Amplitude
        amp_data = deriv_mat[amp_idx,:]
        param_amp = {'data': amp_data.T, 'cmap': cmap_neg_pos, 'alpha': alpha.T, 'vmin': -3, 'vmax': 3, 'cbar': 'discrete'}
        vertex_names.append('amp')

        # Baseline
        baseline_data = deriv_mat[baseline_idx,:]
        param_baseline = {'data': baseline_data.T, 'cmap': cmap_neg_pos, 'alpha': alpha.T, 'vmin': -1, 'vmax': 1,\
                          'curv_brightness': 0.05, 'curv_contrast': 0.1,'cbar': 'discrete'}
        vertex_names.append('baseline')

        # Non-linearity
        non_lin_data = deriv_mat[non_lin_idx,:]
        param_non_lin = {'data': non_lin_data.T, 'cmap': cmap_pos, 'alpha': alpha.T, 'vmin': 0, 'vmax': 1.5, 'cbar': 'discrete'}
        vertex_names.append('non_lin')

        # Coverage
        cov_data = deriv_mat[cov_idx,:]
        param_cov = {'data': cov_data.T, 'cmap': cmap_pos, 'alpha': alpha.T,'vmin': 0, 'vmax': 1, 'cbar': 'discrete'}
        vertex_names.append('cov')

        # Draw figures

        # if fit_model == 'gauss' and subject == '999999': 
        #     dataset_name = 'dataset_{mask_dir}.hdf'.format(mask_dir = mask_dir)
        #     dataset_webgl = cortex.Dataset()

        for vertex_name in vertex_names:
            roi_name = '{vertex_name}_{mask_dir}'.format(vertex_name = vertex_name, mask_dir = mask_dir)
            roi_param = {   'subject': 'fsaverage',
                            'add_roi': False,
                            'roi_name': roi_name}

            exec('param_{vertex_name}.update(roi_param)'.format(vertex_name = vertex_name))
            exec('vertex_rgb = draw_cortex_vertex(**param_{vertex_name})'.format(vertex_name=vertex_name))
            exec('pl.savefig(opj(fig_roi_dir_{mask_dir}, "{vertex_name}_{mask_dir}.pdf"),facecolor="w")'.format(mask_dir=mask_dir,vertex_name = vertex_name))
            # if fit_model == 'gauss' and subject == '999999': 
            #     dataset_webgl.append(**{vertex_name:vertex_rgb})

        # if fit_model == 'gauss' and subject == '999999': 
        #     print('saving dataset: {dataset_name}'.format(dataset_name = dataset_name))
        #     print('cortex.webgl.make_static(outpath = os.path.join(fig_roi_dir_{mask_dir}, data = dataset_name, recache = True))'.format(mask_dir = mask_dir))

        pl.close()