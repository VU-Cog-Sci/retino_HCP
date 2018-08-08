"""
-----------------------------------------------------------------------------------------
pp_roi.py
-----------------------------------------------------------------------------------------
Goal of the script:
Region of interests pre-processing
Compute pRF parameters and plot on pycortex overlay to determine ROI
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject number
sys.argv[2]: voxels per fit (e.g 400)
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
source activate i27
cd /home/szinte/projects/retino_HCP
python retino_HCP/pp_roi.py 192641 400
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import os
import sys
import json
import glob
import numpy as np
import matplotlib.pyplot as pl
import ipdb
import platform

# MRI imports
import nibabel as nb
import cortex

# Function import
from utils import set_pycortex_config_file, convert_fit_results

# from pRF_gazeMod.utils.utils import combine_cv_prf_fit_results_all_runs
# from pRF_gazeMod.utils.prf import convert_fit_results,draw_cortex_volume


# Get inputs
subject = sys.argv[1]
job_vox = float(sys.argv[2])
fit_val = 7
vox_num = 32492
base_file_name = 'tfMRI_RETBAR1_7T_AP_Atlas_MSMAll_hp2000_clean.dtseries'
# Get clock to rename old overlay.svg file
import datetime
now = datetime.datetime.now()

## Check for python version:
sys.exit('Drawing Flatmaps only works with Python 2.x. The current environment seems to have a higher version. Aborting.') if sys.version_info[0] > 2 else None

# Define analysis parameters
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define cluster/server specific parameters
if 'aeneas' in platform.uname()[1]:
    base_dir = analysis_info['aeneas_base_folder'] 
elif 'local' in platform.uname()[1]:
    base_dir = analysis_info['local_base_folder'] 

# change cortex database folder
pycortex_folder     =   os.path.join(base_dir,'pp','cortex')
set_pycortex_config_file(project_folder     =   pycortex_folder)

# Check if all slices are present. If not, the script will abort.
start_idx =  np.arange(0,vox_num,job_vox)
end_idx = start_idx+job_vox
end_idx[-1] = vox_num
num_miss_part = 0
fit_files_L = []
fit_files_R = []
for hemi in ['L','R']:
    for iter_job in np.arange(0,start_idx.shape[0],1):
        fit_file = os.path.join(base_dir,'pp',subject,'prf', '%s_%s.func_bla_psc_est_%s_to_%s.gii' %(base_file_name,hemi,str(int(start_idx[iter_job])),str(int(end_idx[iter_job]))))
        if os.path.isfile(fit_file):
            if os.path.getsize(fit_file) == 0:
                num_miss_part += 1 
            else:
                exec('fit_files_{hemi}.append(fit_file)'.format(hemi = hemi))
        else:
            num_miss_part += 1

if num_miss_part != 0:
    # sys.exit('%i missing files, analysis stopped'%num_miss_part)
    print('%i missing files, partial analysis'%num_miss_part)

# Post-fit pRF analysis
# ---------------------

# combine fit files and save
# left hemi
data_L = []
data_L = np.zeros((fit_val,vox_num))
for fit_filename_L in fit_files_L:
    data_fit_L = []
    data_fit_file_L = nb.load(fit_filename_L)
    data_fit_L.append(np.array([data_fit_file_L.darrays[i].data for i in range(len(data_fit_file_L.darrays))]))
    data_fit_L = np.vstack(data_fit_L)
    data_L = data_L + data_fit_L

darrays_L = [nb.gifti.gifti.GiftiDataArray(d) for d in data_L]
gii_out_L = nb.gifti.gifti.GiftiImage(  header = data_fit_file_L.header, 
                                        extra = data_fit_file_L.extra,
                                        darrays = darrays_L)
nb.save(gii_out_L, os.path.join(base_dir,'pp',subject,'prf','%s_L.func_bla_psc_est.gii'%(base_file_name)))

# right hemi
data_R = []
data_R = np.zeros((fit_val,vox_num))
for fit_filename_R in fit_files_R:
    data_fit_R = []
    data_fit_file_R = nb.load(fit_filename_R)
    data_fit_R.append(np.array([data_fit_file_R.darrays[i].data for i in range(len(data_fit_file_R.darrays))]))
    data_fit_R = np.vstack(data_fit_R)
    data_R = data_R + data_fit_R

darrays_R = [nb.gifti.gifti.GiftiDataArray(d) for d in data_L]
gii_out_R = nb.gifti.gifti.GiftiImage(  header = data_fit_file_R.header, 
                                        extra = data_fit_file_R.extra,
                                        darrays = darrays_R)

nb.save(gii_out_R, os.path.join(base_dir,'pp',subject,'prf','%s_R.func_bla_psc_est.gii'%(base_file_name)))



# compute pRF measures

deriv_dir = os.path.join(base_dir,'pp',subject,'deriv')

for hemi in ['L','R']:
    prf_filename = sorted(glob.glob(os.path.join(base_dir,'pp',subject,'prf','%s_%s.func_bla_psc_est.gii'%(base_file_name, hemi))))
    convert_fit_results(prf_filename = prf_filename,
                        output_dir = deriv_dir,
                        stim_radius = analysis_info['stim_radius'])


# Renaming .svg file before creating a new one below
# --------------------------------------------------
if analysis_info['keep_svg'] == 0:
    try: 
        overlays_file = sorted(glob.glob(os.path.join(analysis_info['FS_subject_dir'], 'cortex', 'db', sub_id,'overlays.svg')))[0]          # "/home/shared/2017/visual/pRF_gazeMod/derivatives/freesurfer/"
        os.rename(overlays_file, overlays_file[:-4] + '_' + now.strftime("%Y_%m_%d_%H_%M") + '.svg')
    except:pass


# Calculate and save x gain ratio and amplitude gain ratio
# --------------------------------------------------------
for mask_dir in ['all']:
    data                                =   {'gazeLeft' : [], 'gazeRight' : [], 'retinal_x_gain' : [], 'screen_x_gain' : [],
                                             'retinal_gain_index' : [], 'amplitude_change' : [], 'y_change' : []}

    exec('deriv_dir_{mask_dir} = os.path.join(analysis_info["aeneas_project_directory"],"pp", sub_id, "deriv_gazeAll","{mask_dir}")'.format(mask_dir = mask_dir))
    try: exec('os.makedirs(deriv_dir_{mask_dir})'.format(mask_dir = mask_dir))                              # make averaging folder
    except OSError: pass


    file_name_gazeLeft                  =   sorted(glob.glob(os.path.join(\
                                                        analysis_info['aeneas_project_directory'],\
                                                        'pp', sub_id,'deriv_gazeLeft',\
                                                        mask_dir,'loo_prf_mean_all.nii.gz')))                              # get all files derived for gaze left
    file_name_gazeRight                 =   sorted(glob.glob(os.path.join(\
                                                        analysis_info['aeneas_project_directory'],\
                                                        'pp', sub_id,'deriv_gazeRight',\
                                                        mask_dir,'loo_prf_mean_all.nii.gz')))                              # get all files derived for gaze right
    data['gazeLeft']                    =   np.array(nb.load(file_name_gazeLeft[0]).get_data())
    data['gazeRight']                   =   nb.load(file_name_gazeRight[0]).get_data()

    data['retinal_x_gain']              =   (data['gazeRight'][..., 0] - data['gazeLeft'][..., 0]) / analysis_info['gazeShift'] * 100
    
    data['screen_x_gain']               =   ((data['gazeRight'][..., 0] + analysis_info['x_grid'][1]) - \
                                            (data['gazeLeft'][..., 0] + analysis_info['x_grid'][0])) / analysis_info['gazeShift'] * 100
    
    data['retinal_gain_index']          =   (analysis_info['gazeShift'] - (data['gazeRight'][..., 0] - data['gazeLeft'][..., 0]))\
                                             / analysis_info['gazeShift']
    
    data['amplitude_change']            =   (data['gazeRight'][..., 3] - data['gazeLeft'][..., 3]) * 100
    data['y_change']                    =   (data['gazeRight'][..., 1] - data['gazeLeft'][..., 1]) * 100




    # Saving the calculated variables
    aff                                 =   nb.load(file_name_gazeLeft[0]).affine                                  # get image affine
    hdr                                 =   nb.load(file_name_gazeLeft[0]).header                                  # get image header
    
    for key in ['retinal_x_gain', 'screen_x_gain', 'retinal_gain_index', 'amplitude_change', 'y_change']:

        nii_file                            =   nb.Nifti1Image(                                                     # create nifti image file
                                                        dataobj         =   data[key],                       # array containing image data
                                                        affine          =   aff,                                # define affine
                                                        header          =   hdr)                                # define header
        file_name                           =   'loo_prf_{0}_all.nii.gz'.format(key)
        exec('nii_file.to_filename(os.path.join(deriv_dir_{mask_dir}, file_name))'.\
                                                            format(mask_dir = mask_dir))                        # write image to filename





for mask_dir in ['all','pos','neg']:

    exec('deriv_dir_{mask_dir} = os.path.join(analysis_info["aeneas_project_directory"],"pp", sub_id, "deriv_gazeAll","{mask_dir}")'.format(mask_dir = mask_dir))
    try: exec('os.makedirs(deriv_dir_{mask_dir})'.format(mask_dir = mask_dir))                              # make averaging folder
    except: pass

    all_files_names_gazeLeft        =   sorted(glob.glob(os.path.join(\
                                                        analysis_info['aeneas_project_directory'],\
                                                        'pp', sub_id,'deriv_%s'%gaze_conditions[0],\
                                                        mask_dir,'*.nii.gz')))                              # get all files derived for gaze left
    all_files_names_gazeRight       =   sorted(glob.glob(os.path.join(\
                                                        analysis_info['aeneas_project_directory'],\
                                                        'pp', sub_id,'deriv_%s'%gaze_conditions[1],\
                                                        mask_dir,'*.nii.gz')))                              # get all files derived for gaze right

    for file_name_gazeLeft, file_name_gazeRight in zip(all_files_names_gazeLeft,\
                                                       all_files_names_gazeRight):                          # define group of files to average
        avg_all_cond                =   np.nanmean(np.array([nb.load(file_name_gazeLeft).get_data(),nb.load\
                                                        (file_name_gazeRight).get_data()]),axis = 0)        # load files to average
        aff                         =   nb.load(file_name_gazeLeft).affine                                  # get image affine
        hdr                         =   nb.load(file_name_gazeLeft).header                                  # get image header
        avg_all_cond_img            =   nb.Nifti1Image(                                                     # create nifti image file
                                                    dataobj         =   avg_all_cond,                       # array containing image data
                                                    affine          =   aff,                                # define affine
                                                    header          =   hdr)                                # define header
        file_name                   =   file_name_gazeLeft.split('/')[-1]
        exec('avg_all_cond_img.to_filename(os.path.join(deriv_dir_{mask_dir}, file_name))'.\
                                                        format(mask_dir = mask_dir))                        # write image to filename



# Create ROI overlay
# ------------------
print('draw roi maps')

# define volume parameters
cv_rsq_th                       =   analysis_info["rsq_threshold_roi"]                                      # cv rsq threshold
cmap_neg_pos                    =   'RdBu_r'                                                                # colormap for negative to positive data
cmap_polar                      =   'hsv'                                                                   # colormap for polar angle (e.g. 'gist_rainbow','Retinotopy_RYBCR')                                                  
cmap_gain                       =   'viridis'
col_offset                      =   1/14.0                                                                  # angle offset necessary to avoid left hemisphere to go to right hemisphere (just for simple visualization)
polar_col_steps                 =   [4.0, 8.0, 16.0, 255.0]                                                 # colormap steps for polar maps    
cmap_ecc_size                   =   'Spectral'                                                              # colormap for eccentricity and size
cmap_pos                        =   'Reds'                                                                  # colormap for positive only data
data_names                      =   {'cv_rsq','ecc','sign','size','amp','baseline',\
                                     'stim_ratio','polar_real','polar_imag',\
                                     'retinal_x_gain', 'screen_x_gain', 'retinal_gain_index'}        # list of volumes to plot
# data_names                      =   {'cv_rsq','size','polar_real','polar_imag','retinal_x_gain', 'screen_x_gain',
#                                      'retinal_gain_index'}        # list of volumes to plot


# create figs_roi folders
figs_roi = os.path.join(analysis_info["aeneas_project_directory"],"pp", sub_id, "figs_roi")                 # create figure folder
try: os.makedirs(figs_roi)                                                                                  # create folder
except OSError: pass

for mask_dir in ['all']:#,'pos','neg']:
    
    vol_names = []
    all_vol   = []

    # create folder
    exec('fig_roi_dir_{mask_dir} = os.path.join(analysis_info["aeneas_project_directory"],"pp", sub_id, "figs_roi","{mask_dir}")'.format(mask_dir=mask_dir)) # create figure folder
    try: exec('os.makedirs(fig_roi_dir_{mask_dir})'.format(mask_dir=mask_dir))                              # create folder
    except: pass

    # load data 
    for data_name in data_names:
        exec('{data_name}_f = os.path.join(deriv_dir_{mask_dir},"loo_prf_{data_name}_{mask_dir}.nii.gz")'.format(mask_dir=mask_dir,data_name = data_name)) # define data to load
        exec('{data_name}_d = nb.load({data_name}_f).get_data()'.format(data_name = data_name))             # load data

    # Cross-validated r-square (loo) or r-square (all)
    cv_rsq_th_mask                  =   cv_rsq_d >= cv_rsq_th                                               # create threshold mask
    rsq_val                         =   np.copy(cv_rsq_d)                                                   # define cv_rsq val for alpha
    rsq_val                         =   np.nan_to_num(rsq_val)                                              # put nan to 0
    rsq_val[rsq_val<0.0]            =   0                                                                   # get rid of negative cv_rsq
    alpha                           =   np.zeros(cv_rsq_d.shape)                                            # pre-define zero (transparent) alpha map
    alpha[cv_rsq_th_mask]           =   1                                                                   # mask it wiht threshold mask
    alpha                           =   alpha * rsq_val                                                     # define non-transparency by cv_rsq
    vmin,vmax                       =   -0.7,0.7                                                            # specify colormap range for cross-validated r-square values
    cmap_rsq                        =   cmap_neg_pos

    param_cv_rsq                    =   {   'subject':          sub_id,                                     # subject in cortex database
                                            'xfmname':          xfmname,                                    # transform name in database
                                            'data':             cv_rsq_d.T,                                 # data to draw
                                            'cmap':             cmap_pos,                                   # colormap to use
                                            'alpha':            alpha.T,                                    # alpha map
                                            'vmin':             vmin,                                       # minimal value in colormap
                                            'vmax':             vmax,                                       # maximal value in colormap 
                                            'cbar':             'discrete'}                                 # color bar layout
    vol_names.append('cv_rsq')
    
    # Polar angle
    pol_comp_num                    =   polar_real_d + 1j * polar_imag_d                                    # define complex number
    polar_ang                       =   np.angle(pol_comp_num)                                              # convert in angle (radian) from -pi to pi
    ang_norm                        =   (polar_ang + np.pi) / (np.pi * 2.0)                                 # add pi and divided by 2pi to have normalized to 1
    
    for cmap_steps in polar_col_steps:
        param_polar                 =   {       'subject':          sub_id,                                 # subject in cortex database
                                                'xfmname':          xfmname,                                # transform name in database
                                                'data':             ang_norm.T,                             # data to draw
                                                'cmap':             cmap_polar,                             # colormap to use
                                                'alpha':            alpha.T,                                # alpha map
                                                'vmin':             0,                                      # minimal value iqn colormap
                                                'vmax':             1,                                      # maximal value in colormap
                                                'cmap_steps':       cmap_steps,                             # steps of colormap colors
                                                'curv_brightness':  0.05,                                   # curvature brightness
                                                'curv_contrast':    0.1,                                    # curvature contrast
                                                'cbar':             'polar',                                # color bar layout
                                                'col_offset':       col_offset}
        exec('param_polar_{csteps} = param_polar'.format(csteps = int(cmap_steps)))                         # get param for each cmap_steps
        exec('vol_names.append("polar_{csteps}")'.format(csteps = int(cmap_steps)))                         # save each volume name
        
    # eccentricity
    param_ecc                       =   {       'subject':          sub_id,                                 # subject in cortex database
                                                'xfmname':          xfmname,                                # transform name in database
                                                'data':             ecc_d.T,                                # data to draw
                                                'cmap':             cmap_ecc_size,                          # colormap to use
                                                'alpha':            alpha.T,                                # alpha map
                                                'vmin':             0,                                      # minimal value iqn colormap
                                                'vmax':             6,                                      # maximal value in colormap
                                                'curv_brightness':  0.05,                                   # curvature brightness
                                                'curv_contrast':    0.1,                                    # curvature contrast
                                                'cbar':             'ecc'}                                  # color bar layout
    vol_names.append('ecc')                                                                                 # save volume name

    # 4. Sign
    param_sign                      =   {       'subject':          sub_id,                                 # subject in cortex database
                                                'xfmname':          xfmname,                                # transform name in database
                                                'data':             sign_d.T,                               # data to draw
                                                'cmap':             cmap_neg_pos,                           # colormap to use
                                                'alpha':            alpha.T,                                # alpha map
                                                'vmin':             -1,                                     # minimal value iqn colormap
                                                'vmax':             1,                                      # maximal value in colormap
                                                'cbar':             'discrete'}                             # color bar layout  
    vol_names.append('sign')                                                                                # save volume name
    
    # Size
    param_size                      =   {       'subject':          sub_id,                                 # subject in cortex database
                                                'xfmname':          xfmname,                                # transform name in database
                                                'data':             size_d.T,                               # data to draw
                                                'cmap':             cmap_ecc_size,                          # colormap to use
                                                'alpha':            alpha.T,                                # alpha map
                                                'vmin':             0,                                      # minimal value iqn colormap
                                                'vmax':             15,                                     # maximal value in colormap
                                                'cbar':             'discrete'}                             # color bar layout
    vol_names.append('size')                                                                                # save volume name

    # Amplitude
    amp_d                           =   np.abs(amp_d)                                                       # take absolute value
    param_amp                       =   {       'subject':          sub_id,                                 # subject in cortex database
                                                'xfmname':          xfmname,                                # transform name in database
                                                'data':             amp_d.T,                                # data to draw
                                                'cmap':             cmap_pos,                               # colormap to use
                                                'alpha':            alpha.T,                                # alpha map
                                                'vmin':             0,                                      # minimal value iqn colormap
                                                'vmax':             1,                                      # maximal value in colormap
                                                'cbar':             'discrete'}                             # color bar layout
    vol_names.append('amp')                                                                                 # save volume name

    # Baseline
    param_baseline                  =   {       'subject':          sub_id,                                 # subject in cortex database
                                                'xfmname':          xfmname,                                # transform name in database
                                                'data':             baseline_d.T,                           # data to draw
                                                'cmap':             cmap_neg_pos,                           # colormap to use
                                                'alpha':            alpha.T,                                # alpha map
                                                'vmin':             -0.5,                                   # minimal value iqn colormap
                                                'vmax':             0.5,                                    # maximal value in colormap
                                                'curv_brightness':  0.05,                                   # curvature brightness
                                                'curv_contrast':    0.1,                                    # curvature contrast
                                                'cbar':             'discrete'}                             # color bar layout
    vol_names.append('baseline')                                                                            # save volume name
    
    # # Non-linearity
    # param_non_lin                   =   {       'subject':          sub_id,                                 # subject in cortex database
    #                                             'xfmname':          xfmname,                                # transform name in database
    #                                             'data':             non_lin_d.T,                            # data to draw
    #                                             'cmap':             cmap_pos,                               # colormap to use
    #                                             'alpha':            alpha.T,                                # alpha map
    #                                             'vmin':             0,                                      # minimal value iqn colormap
    #                                             'vmax':             1.5,                                    # maximal value in colormap
    #                                             'cbar':             'discrete'}                             # color bar layout
        
    # vol_names.append('non_lin')                                                                             # save volume name

    # Stim ratio
    param_stim_ratio                =   {       'subject':          sub_id,                                 # subject in cortex database
                                                'xfmname':          xfmname,                                # transform name in database
                                                'data':             stim_ratio_d.T,                         # data to draw
                                                'cmap':             cmap_pos,                               # colormap to use
                                                'alpha':            alpha.T,                                # alpha map
                                                'vmin':             0,                                      # minimal value iqn colormap
                                                'vmax':             1,                                      # maximal value in colormap
                                                'cbar':             'discrete'}                             # color bar layout
    vol_names.append('stim_ratio')      

    # retinal_x_gain
    param_retinal_x_gain            =   {       'subject':          sub_id,                                 # subject in cortex database
                                                'xfmname':          xfmname,                                # transform name in database
                                                'data':             data['retinal_x_gain'].T,               # data to draw
                                                'cmap':             cmap_gain,                               # colormap to use
                                                'alpha':            alpha.T,                                # alpha map
                                                'vmin':             0,                                      # minimal value iqn colormap
                                                'vmax':             100,                                      # maximal value in colormap
                                                'cbar':             'discrete'}                             # color bar layout
    vol_names.append('retinal_x_gain')                                                                      # save volume name

    # screen_x_gain
    param_screen_x_gain             =   {       'subject':          sub_id,                                 # subject in cortex database
                                                'xfmname':          xfmname,                                # transform name in database
                                                'data':             data['screen_x_gain'].T,                # data to draw
                                                'cmap':             cmap_gain,                               # colormap to use
                                                'alpha':            alpha.T,                                # alpha map
                                                'vmin':             0,                                      # minimal value iqn colormap
                                                'vmax':             100,                                      # maximal value in colormap
                                                'cbar':             'discrete',
                                                }
    vol_names.append('screen_x_gain')                                                                       # save volume name
    
    # retinal gain index
    param_retinal_gain_index        =   {       'subject':          sub_id,                                 # subject in cortex database
                                                'xfmname':          xfmname,                                # transform name in database
                                                'data':             data['retinal_gain_index'].T,           # data to draw
                                                'cmap':             cmap_gain,                              # colormap to use
                                                'alpha':            alpha.T,                                # alpha map
                                                'vmin':             0,                                      # minimal value iqn colormap
                                                'vmax':             1,                                      # maximal value in colormap
                                                'curv_brightness':  0.05,                                   # curvature brightness
                                                'curv_contrast':    0.1,                                    # curvature contrast
                                                'cbar':             'discrete'}                             # color bar layout
    vol_names.append('retinal_gain_index')                                                                  # save volume name



    # draw volumes, add to the overlay, create dataset
    dataset_name                   =   'loo_dataset_{mask_dir}.hdf'.\
                                                        format(mask_dir = mask_dir)                         # define dataset name

    dataset_webgl                   =   cortex.Dataset()                                                    # pre-define dataset

    for vol_name in vol_names:
        roi_name                        =   'loo_{vol_name}_{mask_dir}'.\
                                                    format(vol_name = vol_name, mask_dir = mask_dir)        # define roi name
        
        if mask_dir == 'pos':
            roi_param                   =   {   'add_roi':          True,                                   # add roi to overlay.svg
                                                'roi_name':         roi_name}                               # name of the roi
            exec('param_{vol_name}.update(roi_param)'.format(vol_name = vol_name))
        else:
            roi_param                   =   {   'add_roi':          False,                                  # add roi to overlay.svg
                                            }
            exec('param_{vol_name}.update(roi_param)'.format(vol_name = vol_name))

        exec('volrgb = draw_cortex_volume(**param_{vol_name})'.format(vol_name=vol_name))                   # draw and save a quickshow
        exec('pl.savefig(os.path.join(fig_roi_dir_{mask_dir}, "loo_{vol_name}_{mask_dir}.pdf"),facecolor="w")'.format(mask_dir=mask_dir,vol_name = vol_name))
        dataset_webgl.append(**{vol_name:volrgb})                                                           # add created figure to dataset
    

    print('saving dataset: {dataset_name}'.format(dataset_name = dataset_name))                             # saving dataset
    exec('dataset_webgl.save(os.path.join(fig_roi_dir_{mask_dir}, dataset_name))'.format(mask_dir = mask_dir)) # saving dataset
    pl.close()

