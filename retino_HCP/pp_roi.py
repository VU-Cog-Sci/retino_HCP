"""
-----------------------------------------------------------------------------------------
pp_roi.py
-----------------------------------------------------------------------------------------
Goal of the script:
Region of interests pre-processing
Compute pRF parameters and plot on pycortex overlay to determine ROI
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject name (e.g. 'sub-001')
sys.argv[2]: gaze condition (e.g. 'gazeRight')
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
cd /home/szinte/projects/pRF_analysis/
python pp_roi.py 'sub-004'
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
from utils import set_pycortex_config_file
# from pRF_gazeMod.utils.utils import combine_cv_prf_fit_results_all_runs
# from pRF_gazeMod.utils.prf import convert_fit_results,draw_cortex_volume

# Get clock to rename old overlay.svg file
import datetime
now = datetime.datetime.now()

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
# ---------------------------------------------------------------
# determine iter job for left hemi
ipdb.set_trace()
subject = '536647'
job_vox_L = 240.0
start_idx =  np.arange(0,data_size[1],job_vox_L)
end_idx = start_idx+job_vox_L
for iter_job in np.arange(0,start_idx.shape[0],1):
    opfn = os.path.join(base_dir,'pp',subject,'prf', '*_est_%s_to_%s.gii' %(str(int(start_idx[iter_job])),str(int(end_idx[iter_job]))))
    if os.path.isfile(opfn):
            if os.path.getsize(opfn) != 0:
                print('missing')


project_dir                             =   analysis_info['aeneas_project_directory']
mask_data                               =   nb.load(os.path.join(project_dir, 'pp', sub_id, 'masks',        # get the cortical mask
                                                    'cortex_cortical.nii.gz')).get_data()

slices                                  =   np.arange(mask_data.shape[2])[mask_data.mean(axis=(0,1))>0]     # get only slice with voxels

for gaze_condition in gaze_conditions:

    base_file_names                 =   sorted(glob.glob(os.path.join(project_dir,'pp',sub_id, 'av_' + gaze_condition, 'loo', '*.nii.gz')))
    base_file_names                 =   [file_name.split('/')[-1] for file_name in base_file_names]

    for base_file_name in base_file_names:
        for slice_nr in slices:
            opfn                    =   os.path.join(project_dir,'pp',sub_id, 'cv_' + gaze_condition, 'prf', base_file_name[:-7] + '_est_{0}.nii.gz'.format(slice_nr))

            if os.path.isfile(opfn) == False:
                print(opfn)
                sys.exit('Condition: {gaze_condition}, slice {slice_nr} does not exist, can\'t combine individual slices. aborting.'.format(slice_nr=slice_nr, gaze_condition=gaze_condition))                      # quit code if data exist

## Check for python version:
sys.exit('Drawing Flatmaps only works with Python 2.x. The current environment seems to have a higher version. Aborting.') if sys.version_info[0] > 2 else None

# Renaming .svg file before creating a new one below
# --------------------------------------------------
if analysis_info['keep_svg'] == 0:
    try: 
        overlays_file = sorted(glob.glob(os.path.join(analysis_info['FS_subject_dir'], 'cortex', 'db', sub_id,'overlays.svg')))[0]          # "/home/shared/2017/visual/pRF_gazeMod/derivatives/freesurfer/"
        os.rename(overlays_file, overlays_file[:-4] + '_' + now.strftime("%Y_%m_%d_%H_%M") + '.svg')
    except:pass


# Post-fit pRF analysis
# ---------------------
for gaze_condition in gaze_conditions:

    # define folders
    base_dir                        =   os.path.join(analysis_info['aeneas_project_directory'], 'pp',\
                                                        sub_id,'av_%s'%gaze_condition)                      # avg analysis folder
    fit_dir                         =   os.path.join(analysis_info['aeneas_project_directory'], 'pp',\
                                                        sub_id,'cv_%s'%gaze_condition, 'prf')               # cv analysis folder
    deriv_dir                       =   os.path.join(analysis_info['aeneas_project_directory'],\
                                                        'pp', sub_id, 'deriv_%s'%gaze_condition)            # new derivatives folder
    try: os.makedirs(deriv_dir) 
    except: pass                                                                                            # folder where pRF summary analysis will be stored

    # # combine files
    # output_files                    =   combine_cv_prf_fit_results_all_runs(
    #                                                 basedir         =   base_dir,                           # pp analysis folder
    #                                                 fit_dir         =   fit_dir,                            # pRF analysis folders
    #                                                 avg_folder      =   'loo')                              # leave-one-out or all folder

    # prf_filenames                   =   sorted(glob.glob(os.path.join(fit_dir,'*all.nii.gz')))              # define combine analysis file name

    # # compute pRF measures
    # convert_fit_results(                            prf_filenames   =   prf_filenames,                      # file path to pRF analysis image
    #                                                 output_dir      =   deriv_dir,                          # output derivative folder
    #                                                 stim_radius     =   analysis_info['stim_radius'],       # stimulus radius
    #                                                 typeData        =   'all')                              # type of data to analyse to specify output file name

    # average across conditions
    print('making average across conditions')



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

