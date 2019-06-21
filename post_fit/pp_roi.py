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
#To run:
#cd /home/ada/projects/retino_HCP
#python post_fit/pp_roi_cifti.py 999999 gauss 2500
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
import cifti

# MRI imports
# -----------
import nibabel as nb
import cortex

# Functions import
# ----------------
import utils

# Check system
# ------------
sys.exit('Popeye error with Python 2. Use Python 3 Aborting.') if sys.version_info[0] < 3 else None

#Get inputs
subject = sys.argv[1]
fit_model = sys.argv[2]
job_vox = sys.argv[3]

if fit_model == 'gauss': fit_val = 6
elif fit_model == 'css': fit_val = 7
base_file_name = 'tfMRI_RETALL_Atlas_1.6mm_MSMAll_hp2000_clean_sg_psc'

# Define analysis parameters
# --------------------------
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define cluster/server specific parameters
# -----------------------------------------
if 'aeneas' in platform.uname()[1]:
    base_dir = analysis_info['aeneas_base_folder']
 
elif 'local' in platform.uname()[1]:
    base_dir = analysis_info['local_base_folder'] 
#Output directory:
deriv_dir = opj(base_dir,'pp_data',subject,fit_model,'deriv')

# determine number of vertex and time_series
data_file = opj(base_dir,'raw_data', subject, 'tfMRI_RETALL_Atlas_1.6mm_MSMAll_hp2000_clean_sg_psc.dtseries.nii')
data_file_load = cifti.read(data_file)
data = data_file_load[0]

# Check if all slices are present
start_idx =  np.arange(0,vox_num,job_vox)
end_idx = start_idx+job_vox
end_idx[-1] = vox_num
num_miss_part = 0

fit_est_files = []
for iter_job in np.arange(0,start_idx.shape[0],1):
    fit_est_file = opj(base_dir,'pp_data',subject,fit_model,'fit', '%s_est_%s_to_%s.dtseries.nii' %(base_file_name,str(int(start_idx[iter_job])),str(int(end_idx[iter_job]))))
    if os.path.isfile(fit_est_file):
        if os.path.getsize(fit_est_file) == 0:
            num_miss_part += 1 
        else:
            fit_est_files.append(fit_est_file)
    else:
        num_miss_part += 1

if num_miss_part != 0:
    sys.exit('%i missing files, analysis stopped'%num_miss_part)


#Combine fit files
print('combining fit files')
data_combined = np.zeros((fit_val, vox_num))
for fit_filename in fit_est_files:
    data_fit_file = cifti.read(fit_filename)
    data_fit = data_fit_file[0]
    data_combined = data_combined + data_fit

prf_filename = opj(base_dir,'pp_data',subject,fit_model,'fit',"{bfn}.dtseries.nii".format(bfn= base_file_name))
bm_full = data_fit_file[1][1]
series = cifti.Series(start=0, step=1, size=fit_val)
cifti.write(prf_filename, data_combined, (series, bm_full)) 

# Compute derived measures from prfs
# ----------------------------------
print('extracting pRF derivatives')

utils.convert_fit_results(prf_filename = prf_filename,
                    output_dir = deriv_dir,
                    stim_radius = analysis_info['stim_radius'],
                    fit_model = fit_model)

# Resample gii to fsaverage
# -------------------------
print('converting derivative files to fsaverage')
resample_cmd = """{main_cmd} -metric-resample {metric_in} {current_sphere} {new_sphere} ADAP_BARY_AREA {metric_out} -area-metrics {current_area} {new_area}"""
for hemi in ['L','R']:

    current_sphere = opj(base_dir,'raw_data/surfaces/resample_fsaverage','fs_LR-deformed_to-fsaverage.{hemi}.sphere.{num_vox_k}k_fs_LR.surf.gii'.format(hemi=hemi, num_vox_k = int(np.round(vox_num/1000))))
    new_sphere = opj(base_dir,'raw_data/surfaces/resample_fsaverage','fsaverage_std_sphere.{hemi}.164k_fsavg_{hemi}.surf.gii'.format(hemi=hemi))
    current_area = opj(base_dir,'raw_data/surfaces/resample_fsaverage','fs_LR.{hemi}.midthickness_va_avg.{num_vox_k}k_fs_LR.shape.gii'.format(hemi=hemi,num_vox_k = int(np.round(vox_num/1000))))
    new_area = opj(base_dir,'raw_data/surfaces/resample_fsaverage','fsaverage.{hemi}.midthickness_va_avg.164k_fsavg_{hemi}.shape.gii'.format(hemi=hemi))

    for mask_dir in ['all','pos','neg']:
        
        metric_in = opj(deriv_dir,mask_dir,"prf_deriv_{hemi}_{mask_dir}.gii".format(hemi = hemi, mask_dir = mask_dir))
        metric_out = opj(deriv_dir,mask_dir,"prf_deriv_{hemi}_{mask_dir}_fsaverage.func.gii".format(hemi = hemi, mask_dir = mask_dir))

        os.system(resample_cmd.format(  main_cmd = main_cmd,
                                        metric_in = metric_in, 
                                        current_sphere = current_sphere, 
                                        new_sphere = new_sphere, 
                                        metric_out = metric_out, 
                                        current_area = current_area, 
                                        new_area = new_area))




if draw_roi:
    
    # Change cortex database folder
    # -----------------------------
    pycortex_folder     =   opj(base_dir,'pp_data','cortex')
    set_pycortex_config_file(project_folder = pycortex_folder)
    
    # Create derivatives flatmaps
    # ---------------------------
    print('draw deriv maps')
    cmap_neg_pos = 'RdBu_r'
    cmap_polar = 'hsv'
    col_offset = 1/14.0
    polar_col_steps = [4.0, 8.0, 16.0, 255.0]
    cmap_ecc_size = 'Spectral'
    cmap_pos = 'Reds'
    sign_idx, rsq_idx, ecc_idx, polar_real_idx, polar_imag_idx , size_idx, \
                non_lin_idx, amp_idx, baseline_idx, cov_idx, x_idx, y_idx = 0,1,2,3,4,5,6,7,8,9,10,11

    for mask_dir in ['all','pos','neg']:

        # Create figure folders
        vertex_names = []
        all_vertex   = []
        exec('fig_roi_dir_{mask_dir} = opj(base_dir,"pp_data",subject,fit_model,"figs","roi","{mask_dir}")'.format(mask_dir=mask_dir))
        try: exec('os.makedirs(fig_roi_dir_{mask_dir})'.format(mask_dir=mask_dir))
        except: pass

        # Combine hemispheres
        deriv_mat=[]
        for hemi in ['L','R']:
            deriv_file = nb.load(opj(deriv_dir,mask_dir,"prf_deriv_{hemi}_{mask_dir}_fsaverage.func.gii".format(hemi = hemi, mask_dir = mask_dir)))
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
        param_amp = {'data': amp_data.T, 'cmap': cmap_neg_pos, 'alpha': alpha.T, 'vmin': -0.005, 'vmax': 0.005, 'cbar': 'discrete'}
        vertex_names.append('amp')

        # Baseline
        baseline_data = deriv_mat[baseline_idx,:]
        param_baseline = {'data': baseline_data.T, 'cmap': cmap_neg_pos, 'alpha': alpha.T, 'vmin': -10, 'vmax': 10,\
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