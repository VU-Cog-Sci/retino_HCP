"""
-----------------------------------------------------------------------------------------
post_pp_roi.py
-----------------------------------------------------------------------------------------
Goal of the script:
Draw basic plot of pRF analysis
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject number
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
source activate i27
cd /home/szinte/projects/retino_HCP
python retino_HCP/post_pp_roi.py 192641
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
import h5py

# MRI imports
import nibabel as nb
import cortex

# Define analysis parameters
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define cluster/server specific parameters
if 'aeneas' in platform.uname()[1]:
    base_dir = analysis_info['aeneas_base_folder'] 
elif 'local' in platform.uname()[1]:
    base_dir = analysis_info['local_base_folder'] 

# Define input
subject = sys.argv[1]
deriv_dir = os.path.join(base_dir,'pp',subject,'deriv')

# # Get index of region of interest
# altas_file = h5py.File('%s/atlas.mat'%base_dir, 'r')
# val_vox_roi_num = altas_file['/glasser2016'].value
# roi_ref = { 'V1':2,
#             'V2':5,
#             'V3':6,
#             'V4':7,
#             'V3A':14,
#             'V3B':20,
#             'LO1':21,
#             'LO2':22,
#             'MT':24,
#             'MST':3,
#             'FEF':11,
#             'IPS1':18}



# load data


# Draw main analysis figure
rois = ['V1','V2']


for hemi in ['L','R']:
    if hemi == 'L':
        altas_idx = np.arange(0,32492,1)
    elif hemi == 'R':
        altas_idx = np.arange(0,32492,1)

    for mask_dir in ['all','pos','neg']:

        deriv_mat = []
        deriv_file = os.path.join(deriv_dir,'%s'%mask_dir,'prf_deriv_%s_%s.gii'%(hemi,mask_dir))
        deriv_file_load = nb.load(deriv_file)
        deriv_mat.append(np.array([deriv_file_load.darrays[i].data for i in range(len(deriv_file_load.darrays))]))
        deriv_mat = np.vstack(deriv_mat)

        roi_mask = altas_file

        ipdb.set_trace()

        




# numArea = 0
# for roi_text in rois:
#     # Get data
    
#     pRF_fit_L                       =   np.zeros((pRF_mean_L.shape[0],9))                                           #  create empty ndarray
#     pRF_fit_L[:,:-2]                =   pRF_mean_L[:,:7]                                                            # first 6 columns (exclude 2 empty column)
#     pRF_fit_L[:,-2]                 =   pRF_stim_ratio_L                                                            # add last columns
#     pRF_fit_L[:,-1]                 =   pRF_ecc_L                                                                   # add last columns

#     pRF_fit_L_mask                  =   pRF_fit_L[:,6] >= analysis_info['rsq_threshold']                            # define rsq mask for left hemi
#     pRF_fit_L                       =   pRF_fit_L[pRF_fit_L_mask]                                                   # mask left hemi data

    
#     pRF_fit_R                       =   np.zeros((pRF_mean_R.shape[0],9))                                           # create empty ndarray
#     pRF_fit_R[:,:-2]                =   pRF_mean_R[:,:7]                                                            # first 6 columns (exclude 2 empty column)
#     pRF_fit_R[:,-2]                 =   pRF_stim_ratio_R                                                            # add last column
#     pRF_fit_R[:,-1]                 =   pRF_ecc_R                                                                   # add last column

#     pRF_fit_R_mask                  =   pRF_fit_R[:,6] >= analysis_info['rsq_threshold']                            # define rsq mask for right hemi
#     pRF_fit_R                       =   pRF_fit_R[pRF_fit_R_mask]                                                   # mask left hemi data

#     # get both hemisphere
#     pRF_fit_LR                      =   np.zeros((pRF_fit_L.shape[0]+pRF_fit_R.shape[0],pRF_fit_R.shape[1]))        # create empty ndarray
#     pRF_fit_LR[0:pRF_fit_L.shape[0],:] = pRF_fit_L                                                                  # fill it with left hemi.
#     pRF_fit_LR[pRF_fit_L.shape[0]:,:] =  pRF_fit_R                                                                  # and fill it with right hemi.

#     # define main parameters
#     param_all                       =   dict(
#                                             roi_t                   =   roi_text,                                   # current region of interest
#                                             p_width                 =   400,                                        # individual plot width
#                                             p_height                =   400,                                        # individual plot height
#                                             min_border_large        =   10,                                         # large border between figure and axis
#                                             min_border_small        =   5,                                          # small border between figure and axis
#                                             bg_color                =   tuple([229,229,229]),                       # background color         
#                                             stim_color              =   tuple([250,250,250]),                       # stimuli color         
#                                             hist_fill_color         =   tuple([255,255,255]),                       # histogram color
#                                             hist_line_color         =   tuple([0,0,0]),                             # histogram color
#                                             stim_radius             =   analysis_info['stim_radius'],               # stimulus radius
#                                             cmap                    =   'Spectral',                                 # colormap to use
#                                             cmap_steps              =   8,                                          # colormap steps
#                                             col_offset              =   0,                                          # colormap offset
#                                             vmin                    =   0,                                          # val min to draw with colormap
#                                             vmax                    =   4,                                          # val max to draw with colormap
#                                             leg_xy_max_ratio        =   1.8                                         # xy axis max based on vmax
#                                             )

#     # pRF map
#     # -------
#     params_pRFmap                   =   param_all                                                                   # get main params
#     _                               =   params_pRFmap.update(                                                       # update parameters
#                                             dict(
#                                             x_range                 =   (-8, 8),                                    # x axis limits
#                                             y_range                 =   (-8, 8),                                    # y axis limits
#                                             x_label                 =   'Horizontal coordinate (dva)',              # x axis label
#                                             y_label                 =   'Vertical coordinate (dva)',                # y axis label
#                                             x_tick_steps            =   2,                                          # x axis major tick steps
#                                             y_tick_steps            =   2,                                          # y axis major tick steps
#                                             hist_range              =   (0,0.5),                                    # histogram x/y range
#                                             hist_steps              =   0.5,                                        # x axis major tick steps
#                                             h_hist_bins             =   16,                                         # hor. histogram bins
#                                             v_hist_bins             =   16)                                         # ver. histogram bins
#                                             )

#     # Left hemisphere
#     title                           =   '{}-left hemisphere: pRF position and size (n={})'.format(
#                                                                             roi_text,pRF_fit_L.shape[0])            # define titlte
#     _                               =   params_pRFmap.update(                                                       # update parameters
#                                             dict(              
#                                             dataMat                 =   pRF_fit_L,                                  # main data matrix
#                                             main_fig_title          =   title)                                      # title
#                                             )
#     f_pRFmap_lh,_                   =   draw_pRFmap(                                                                # draw pRF map figure
#                                             params                  =   params_pRFmap)                              # figure parameter

#     # Right hemisphere
#     title                           =   '{}-right hemisphere: pRF position and size (n={})'.format(
#                                                                             roi_text,pRF_fit_R.shape[0])            # define title
#     _                               =   params_pRFmap.update(                                                       # update parameters
#                                             dict(
#                                             dataMat                 =   pRF_fit_R,                                  # main data matrix
#                                             main_fig_title          =   title)                                      # title
#                                             )
#     f_pRFmap_rh,_                   =   draw_pRFmap(                                                                # draw pRF map figure
#                                             params                  =   params_pRFmap)                              # figure parameter

#     # Both hemisphere
#     title                           =   '{}: pRF position and size (n={})'.format(roi_text,pRF_fit_LR.shape[0])     # define title
#     _                               =   params_pRFmap.update(                                                       # update parameters
#                                             dict(
#                                             dataMat                 =   pRF_fit_LR,                                 # main data matrix
#                                             main_fig_title          =   title)                                      # title
#                                             )
#     f_pRFmap_lrh,_                  =   draw_pRFmap(                                                                # draw pRF map figure
#                                             params                  =   params_pRFmap)                              # figure parameter

#     # pRF cov
#     # -------
#     params_pRFcov                   =   param_all                                                                   # get main params
#     _                               =   params_pRFcov.update(                                                       # update parameters
#                                             dict(
#                                             x_range                 =   (-8, 8),                                    # x axis limits
#                                             y_range                 =   (-8, 8),                                    # y axis limits
#                                             x_label                 =   'Horizontal coordinate (dva)',              # x axis label
#                                             y_label                 =   'Vertical coordinate (dva)',                # y axis label
#                                             x_tick_steps            =   2,                                          # x axis major tick steps
#                                             y_tick_steps            =   2,                                          # y axis major tick steps
#                                             smooth_factor           =   15,                                         # pixel per degree in image
#                                             cmap                    =   'viridis',                                  # colormap to use
#                                             cmap_steps              =   10,                                         # colormap steps
#                                             col_offset              =   0,                                          # colormap offset
#                                             vmin                    =   0,                                          # val min to draw with colormap
#                                             vmax                    =   1,                                          # val max to draw with colormap
#                                             colorbar_tick_steps     =   0.2,                                        # colorbar axis major tick steps
#                                             colorbar_label          =   'pRF coverage (norm.)')                     # colorbar label
#                                             )

#     # Left hemisphere
#     title                           =   '{}-left hemisphere: pRF coverage (n={})'.format(
#                                                                             roi_text,pRF_fit_L.shape[0])            # define titlte
#     _                               =   params_pRFcov.update(                                                       # update parameters
#                                             dict(              
#                                             dataMat                 =   pRF_fit_L,                                  # main data matrix
#                                             main_fig_title          =   title)                                      # title
#                                             )
#     f_pRFcov_lh,_                   =   draw_pRFcov(                                                                # draw pRF map figure
#                                             params                  =   params_pRFcov)                              # figure parameter

#     # Right hemisphere
#     title                           =   '{}-right hemisphere: pRF coverage (n={})'.format(
#                                                                             roi_text,pRF_fit_R.shape[0])            # define title
#     _                               =   params_pRFcov.update(                                                       # update parameters
#                                             dict(
#                                             dataMat                 =   pRF_fit_R,                                  # main data matrix
#                                             main_fig_title          =   title)                                      # title
#                                             )
#     f_pRFcov_rh,_                   =   draw_pRFcov(                                                                # draw pRF map figure
#                                             params                  =   params_pRFcov)                              # figure parameter

#     # Both hemisphere
#     title                           =   '{}: pRF coverage (n={})'.format(roi_text,pRF_fit_LR.shape[0])     # define title
#     _                               =   params_pRFcov.update(                                                       # update parameters
#                                             dict(
#                                             dataMat                 =   pRF_fit_LR,                                 # main data matrix
#                                             main_fig_title          =   title)                                      # title
#                                             )
#     f_pRFcov_lrh,_                  =   draw_pRFcov(                                                                # draw pRF map figure
#                                             params                  =   params_pRFcov)                              # figure parameter

#     # pRF ecc
#     # -------
#     f_pRFecc_lh                         =   []
#     f_pRFecc_rh                         =   []
#     f_pRFecc_lrh                        =   []
#     first_pRFecc_fig                    =   [] 
#     numData                             =   0

#     for typeData in ['Size','CV R2','Amplitude','Stim-ratio']:

#         params_pRFecc                   =   param_all                                                               # get main params     
#         _                               =   params_pRFecc.update(                                                   # update parameters
#                                                 dict(
#                                                 x_range                 =   (0, 10),                                # x axis limits
#                                                 x_label                 =   'Eccentricity (dva)',                   # x axis label
#                                                 x_tick_steps            =   2,                                      # x axis major tick steps
#                                                 hist_range              =   (0,0.5),                                # histogram x/y range
#                                                 hist_steps              =   0.5,                                    # x axis major tick steps
#                                                 h_hist_bins             =   20,                                     # hor. histogram bins
#                                                 cmap                    =   'Spectral',                                 # colormap to use
#                                                 cmap_steps              =   8,                                          # colormap steps
#                                                 col_offset              =   0,                                          # colormap offset
#                                                 vmin                    =   0,                                          # val min to draw with colormap
#                                                 vmax                    =   4)                                          # val max to draw with colormap
#                                                 )

#         if typeData == 'Size':
#             _                           =   params_pRFecc.update(                                                   # update parameters
#                                                 dict(
#                                                 y_range                 =   (0, 4),                                 # y axis limits (adapted for size)
#                                                 y_label                 =   'Size (dva)',                           # y axis label
#                                                 y_source_label          =   'sigma',                                # y axis source val in data source
#                                                 y_tick_steps            =   1,                                      # y axis major tick steps
#                                                 v_hist_bins             =   16)                                     # ver. histogram bins
#                                                 )
#         elif typeData == 'CV R2':
#             _                           =   params_pRFecc.update(                                                   # update parameters
#                                                 dict(
#                                                 y_range                 =   (0, 1),                                 # y axis limits (adapted for size)
#                                                 y_label                 =   'Cross-validated R2',                   # y axis label
#                                                 y_source_label          =   'cv_rsq',                               # y axis source val in data source
#                                                 y_tick_steps            =   0.2,                                    # y axis major tick steps
#                                                 v_hist_bins             =   20)                                     # ver. histogram bins
#                                                 )
#         elif typeData == 'Amplitude':
#             _                           =   params_pRFecc.update(                                                   # update parameters
#                                                 dict(
#                                                 y_range                 =   (0, 1),                                 # y axis limits (adapted for size)
#                                                 y_label                 =   'Amplitude (z-score)',                  # y axis label
#                                                 y_source_label          =   'beta',                                 # y axis source val in data source
#                                                 y_tick_steps            =   0.25,                                   # y axis major tick steps
#                                                 v_hist_bins             =   16)                                     # ver. histogram bins   
#                                                 )

#         elif typeData == 'Stim-ratio':
#             _                           =   params_pRFecc.update(                                                   # update parameters
#                                                 dict(
#                                                 y_range                 =   (-10, 110),                             # y axis limits (adapted for size)
#                                                 y_label                 =   'Stimulus ratio (%)',                   # y axis label
#                                                 y_source_label          =   'stim_ratio',                           # y axis source val in data source
#                                                 y_tick_steps            =   10,                                     # y axis major tick steps
#                                                 v_hist_bins             =   12)                                     # ver. histogram bins
#                                                 )
#         # Left hemisphere
#         title                           =   '{}-left hemisphere: Eccentricity vs. {}'.format(roi_text,typeData)     # define title
#         _                               =   params_pRFecc.update(                                                   # update parameters
#                                                 dict(
#                                                 dataMat                 =   pRF_fit_L,                              # main data matrix
#                                                 main_fig_title          =   title)                                  # title
#                                                 )
#         out1, out2                      =   draw_pRFecc(                                                            # draw pRF ecc figure
#                                                 params                  =   params_pRFecc)                          # figure parameter
#         _                               =   f_pRFecc_lh.append(out1)
#         _                               =   first_pRFecc_fig.append(out2)

#         # Right hemisphere
#         title                           =   '{}-right hemisphere: Eccentricity vs. {}'.format(roi_text,typeData)    # define title
#         _                               =   params_pRFecc.update(                                                   # update parameters
#                                                 dict(              
#                                                 dataMat                 =   pRF_fit_R,                              # main data matrix
#                                                 main_fig_title          =   title))                                 # title

#         out1, _                         =   draw_pRFecc(                                                            # draw pRF ecc figure
#                                             params                  =   params_pRFecc,                              # figure parameter
#                                             old_main_fig            =   first_pRFecc_fig[numData])                  # handle to first fig
#         _                               =   f_pRFecc_rh.append(out1)

#         # Both hemispheres
#         title                           =   '{}: Eccentricity vs. {}'.format(roi_text,typeData)                     # define title 
#         _                               =   params_pRFecc.update(                                                   # update parameters
#                                                 dict(
#                                                 dataMat                 =   pRF_fit_LR,                             # main data matrix
#                                                 main_fig_title          =   title))                                 # title
#         out1, _                         =   draw_pRFecc(                                                            # draw pRF ecc figure
#                                                 params                  =   params_pRFecc,                          # figure parameter
#                                                 old_main_fig            =   first_pRFecc_fig[numData])              # handle to first fig
#         _                               =   f_pRFecc_lrh.append(out1)

#         numData                         =   numData + 1

#     # Combine figures
#     # ---------------
#     # left hemisphere
#     all_fL                          =   gridplot([
#                                             [f_pRFmap_lh,       f_pRFecc_lh[0],     f_pRFecc_lh[1]],                # specify figure 1st row
#                                             [f_pRFcov_lh,       f_pRFecc_lh[2],    f_pRFecc_lh[3]]                  # specify figure 2nd row
#                                             ])
#     output_L_file_html              =   opj(pRFbasic_folder_L,'%sL_pRF.html'%roi_text)                              # define html file
#     _                               =   output_file(output_L_file_html, title='%s L pRF analysis'%roi_text)         # get html saved
#     _                               =   save(all_fL)


#     # right hemisphere
#     all_fR                          =   gridplot([
#                                             [f_pRFmap_rh,       f_pRFecc_rh[0],     f_pRFecc_rh[1]],                # specify figure 1st row
#                                             [f_pRFcov_rh,       f_pRFecc_rh[2],     f_pRFecc_rh[3]]                 # specify figure 2nd row
#                                             ])
#     output_R_file_html              =   opj(pRFbasic_folder_R,'%sR_pRF.html'%roi_text)                              # define html file
#     _                               =   output_file(output_R_file_html, title='%s R pRF analysis'%roi_text)         # get html saved
#     _                               =   save(all_fR)

#     # both hemisphere
#     all_fLR                         =   gridplot([
#                                             [f_pRFmap_lrh,      f_pRFecc_lrh[0],    f_pRFecc_lrh[1]],               # specify figure 1st row
#                                             [f_pRFcov_lrh,      f_pRFecc_lrh[2],    f_pRFecc_lrh[3]]                # specify figure 2nd row
#                                             ])
#     output_LR_file_html             =   opj(pRFbasic_folder,'%s_pRF.html'%roi_text)                                 # define html file
#     _                               =   output_file(output_LR_file_html, title='%s pRF analysis'%roi_text)          # get html saved
#     _                               =   save(all_fLR)