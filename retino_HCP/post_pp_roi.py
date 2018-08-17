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
import h5py
opj = os.path.join
deb = ipdb.set_trace

# MRI imports
# -----------
import nibabel as nb
import cortex

# Functions import
# ----------------
from plot_class import PlotOperator

# Bokeh imports
# ---------------
from bokeh.io import output_notebook, show, save, output_file, export_png, export_svgs
from bokeh.layouts import row, column, gridplot

# Get inputs
# ----------
subject = sys.argv[1]

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
deriv_dir = opj(base_dir,'pp',subject,'deriv')
h5_dir = opj(base_dir,'pp',subject,'h5')

# Draw main analysis figure
# -------------------------
print('creating bokeh plots')

# Initialize data dictionary that will save all data arrays

for roi in analysis_info['rois']:
    
    for mask_dir in ['all','pos','neg']:

        for hemi in ['L', 'R']:

            # create folders
            exec('fig_bokeh_dir_{mask_dir} = opj(base_dir,"pp",subject,"figs_bokeh","{mask_dir}")'.format(mask_dir=mask_dir))
            try: exec('os.makedirs(fig_bokeh_dir_{mask_dir})'.format(mask_dir=mask_dir))
            except: pass

            # load data
            folder_alias = '{hemi}_{mask_dir}'.format(hemi = hemi,mask_dir = mask_dir)
            h5_file = h5py.File(opj(h5_dir,'{roi}.h5'.format(roi = analysis_info['rois'][roi])), "r")
            in_file = opj("prf_deriv_{hemi}_{mask_dir}".format(hemi = hemi, mask_dir = mask_dir))
            data = h5_file['{folder_alias}/{in_file}'.format(folder_alias=folder_alias,in_file=in_file)]
            data = data[:,data[1,:]>=analysis_info['rsq_threshold']]
            data = data[:,data[9,:]>=analysis_info['cov_threshold']]
            data = data.T

            roi_text = analysis_info['rois'][roi]
            vertex = data.shape[0]
            
            param_all = {   'roi_t':            analysis_info['rois'][roi], 
                            'p_width':          400, 
                            'p_height':         400, 
                            'min_border_large': 10, 
                            'min_border_small': 5,
                            'bg_color':         tuple([229,229,229]), 
                            'stim_color':       tuple([250,250,250]), 
                            'hist_fill_color':  tuple([255,255,255]),
                            'hist_line_color':  tuple([0,0,0]), 
                            'stim_radius':      analysis_info['stim_radius'], 
                            'rsq_threshold':    analysis_info['rsq_threshold'],
                            'cmap':             'Spectral',
                            'cmap_steps':       8,
                            'col_offset':       0,
                            'vmin':             0,
                            'vmax':             8,
                            'leg_xy_max_ratio': 1.8,
                            # 'saving_figs':      False,
                            # 'svg_folder':       'fig_bokeh_dir_{mask_dir}'.format(mask_dir=mask_dir),
                            'dataMat':          data
                            }

            plotter = PlotOperator(**param_all)

            # pRFmap
            title = '{roi}{hemi} hemisphere: pRFs (n={vertex})'.format(roi = roi_text, hemi = hemi, vertex = vertex)
            params_pRFmap = param_all
            params_pRFmap.update(\
                        {   'x_range':          (-15, 15),
                            'y_range':          (-15, 15),
                            'x_label':          'Horizontal coordinate (dva)',
                            'y_label':          'Vertical coordinate (dva)',
                            'x_source_label':   'x',
                            'y_source_label':   'y',
                            'x_tick_steps':     5,
                            'y_tick_steps':     5,
                            'hist_range':       (0,0.5),
                            'hist_steps':       5,
                            'h_hist_bins':      16,
                            'v_hist_bins':      16
                            'main_fig_title':   title})

            f_pRFmap = plotter.draw_figure( parameters = params_pRFmap, 
                                            plot = 'map')
            
            # pRFecc
            old_main_fig = []
            for numData, type_comp in enumerate(['Size','R2','Non-Linearity','Amplitude','Coverage']):

                params_pRFecc = param_all
                params_pRFecc.update(   
                           {    'x_range':          (0, 10),
                                'x_label':          'Eccentricity (dva)',
                                'x_tick_steps':     2,
                                'x_source_label':   'ecc',
                                'vmin':             0,
                                'vmax':             4,
                                'hist_range':       (0,0.5),
                                'hist_steps':       0.5,
                                'h_hist_bins':      20})

                if type_comp == 'Size':
                    params_pRFecc.update(
                                {   'y_range':          (0, 4),
                                    'y_label':          'Size (dva)',
                                    'y_source_label':   'sigma',
                                    'y_tick_steps':     1,
                                    'v_hist_bins':      16})

                elif type_comp == 'R2':
                    params_pRFecc.update(
                                {   'y_range':          (0, 1)
                                    'y_label':          'R2',
                                    'y_source_label':   'rsq',
                                    'y_tick_steps':     0.2,
                                    'v_hist_bins':      10})

                elif type_comp == 'Non-Linearity'
                    params_pRFecc.update(
                                {   'y_range':          (0, 2)
                                    'y_label':          'Non-linearity',
                                    'y_source_label':   'non_lin',
                                    'y_tick_steps':     0.4,
                                    'v_hist_bins':      20})

                elif type_comp == 'Amplitude':
                    params_pRFecc.update(
                                {   'y_range':          (0, 1),
                                    'y_label':          'Amplitude (z-score)',
                                    "y_source_label":   'beta',
                                    'y_tick_steps':     0.25,
                                    'v_hist_bins':      10})

                elif type_comp == 'Coverage':
                    params_pRFecc.update(
                                {   'y_range':          (-10, 110),
                                    'y_label':          'Stimulus ratio (%)',
                                    'y_source_label':   'stim_ratio',
                                    'y_tick_steps':     10,
                                    'v_hist_bins':      12})

                title = '{roi}{hemi} hemisphere: Eccentricity vs. {typeData}'.format(roi = roi_text, hemi = hemi, type_comp = type_comp)
                params_pRFecc.update({'main_fig_title':   title})

                out1, old_main_fig  = plotter.draw_figure(  parameters = params_pRFecc,
                                                            plot = 'ecc',
                                                            old_main_fig = old_main_fig)
                f_pRFecc[numData].append(out1)
                



            all_fL = gridplot([
                                [f_pRFmap[0],f_pRFecc[0],None]                # specify figure 1st row
                                ])
            
            
            exec('output_file_html = opj(fig_bokeh_dir_{mask_dir}_{hemi},"{roi_text}_pRF.html")'.format(mask_dir = mask_dir,roi_text = roi_text, hemi = hemi))
            output_file(output_file_html, title='%s pRF analysis'%roi_text)
            save(all_fL)
            
            
            
    
    # # pRF ecc
    # # -------
    # f_pRFecc_basic                      =   {'lh' : [], 'rh' : [], 'lrh' : []}
    # first_pRFecc_fig_basic              =   []

    # 

    #     params_pRFecc                   =   param_all                                                               # get main params     
    #     _                               =   params_pRFecc.update(                                                   # update parameters
    #                                             dict(
    #                                             x_range                 =   (0, 10),                                # x axis limits
    #                                             x_label                 =   'Eccentricity (dva)',                   # x axis label
    #                                             x_tick_steps            =   2,                                      # x axis major tick steps
    #                                             x_source_label          =   'ecc',
    #                                             vmin                    =   0,
    #                                             vmax                    =   4,
    #                                             hist_range              =   (0,0.5),                                # histogram x/y range
    #                                             hist_steps              =   0.5,                                    # x axis major tick steps
    #                                             h_hist_bins             =   20),                                     # hor. histogram bins
    #                                             condition               =   'ecc'
    #                                             )

    #     if typeData == 'Size':
    #         _                           =   params_pRFecc.update(                                                   # update parameters
    #                                             dict(
    #                                             y_range                 =   (0, 4),                                 # y axis limits (adapted for size)
    #                                             y_label                 =   'Size (dva)',                           # y axis label
    #                                             y_source_label          =   'sigma',                                # y axis source val in data source
    #                                             y_tick_steps            =   1,                                      # y axis major tick steps
    #                                             v_hist_bins             =   16)                                      # ver. histogram bins
    #                                             )
    #     elif typeData == 'CV R2':
    #         _                           =   params_pRFecc.update(                                                   # update parameters
    #                                             dict(
    #                                             y_range                 =   (0, 1),                                 # y axis limits (adapted for size)
    #                                             y_label                 =   'Cross-validated R2',                   # y axis label
    #                                             y_source_label          =   'cv_rsq',                               # y axis source val in data source
    #                                             y_tick_steps            =   0.2,                                    # y axis major tick steps
    #                                             v_hist_bins             =   10)                                     # ver. histogram bins
    #                                             )
    #     elif typeData == 'R2':
    #         _                           =   params_pRFecc.update(                                                   # update parameters
    #                                             dict(
    #                                             y_range                 =   (0, 1),                                 # y axis limits (adapted for size)
    #                                             y_label                 =   'R2',                                   # y axis label
    #                                             y_source_label          =   'rsq',                                  # y axis source val in data source
    #                                             y_tick_steps            =   0.2,                                    # y axis major tick steps
    #                                             v_hist_bins             =   10)                                     # ver. histogram bins
    #                                             )
    #     elif typeData == 'Amplitude':
    #         _                           =   params_pRFecc.update(                                                   # update parameters
    #                                             dict(
    #                                             y_range                 =   (0, 1),                                 # y axis limits (adapted for size)
    #                                             y_label                 =   'Amplitude (z-score)',                  # y axis label
    #                                             y_source_label          =   'beta',                                 # y axis source val in data source
    #                                             y_tick_steps            =   0.25,                                   # y axis major tick steps
    #                                             v_hist_bins             =   10)                                     # ver. histogram bins   
    #                                             )

    #     elif typeData == 'Stim-ratio':
    #         _                           =   params_pRFecc.update(                                                   # update parameters
    #                                             dict(
    #                                             y_range                 =   (-10, 110),                             # y axis limits (adapted for size)
    #                                             y_label                 =   'Stimulus ratio (%)',                   # y axis label
    #                                             y_source_label          =   'stim_ratio',                           # y axis source val in data source
    #                                             y_tick_steps            =   10,                                     # y axis major tick steps
    #                                             v_hist_bins             =   12)                                     # ver. histogram bins
    #                                             )
        

    #     for hemi, list_name, ext in [('-left', 'lh', '_L'), ('-right', 'rh', '_R'), ('', 'lrh', '_LR')]:
    #         title                           =   '{roi}{hemi} hemisphere: Eccentricity vs. {typeData}'.\
    #                                              format(roi = roi_text, hemi = hemi, typeData = typeData)           # define title
    #         _                               =   params_pRFecc.update(                                                   # update parameters
    #                                                 dict(
    #                                                 dataMat                 =   data['pRF_fit%s'%ext],                              # main data matrix
    #                                                 ecc_gainRatio_counter   =   0,
    #                                                 main_fig_title          =   title,
    #                                                 hemisphere              =   list_name,
    #                                                 typeData                =   typeData))                                  # title
            
    #         old_fig = [] if list_name == 'lh' else first_pRFecc_fig_basic[numData]

    #         out1, out2                      =   plotter.draw_figure(parameters = params_pRFecc,
    #                                                                 # old_main_fig = old_fig,                                                                    
    #                                                                 plot = 'ecc')
    #         _                               =   f_pRFecc_basic[list_name].append(out1)
    #         _                               =   first_pRFecc_fig_basic.append(out2)


    # # pRF cov
    # # -------
    # f_pRFcov                        =   {'lh' : [], 'rh' : [], 'lrh' : []}
    # params_pRFcov                   =   param_all                                                                   # get main params
    # _                               =   params_pRFcov.update(                                                       # update parameters
    #                                         dict(
    #                                         x_range                 =   (-8, 8),                                    # x axis limits
    #                                         y_range                 =   (-8, 8),                                    # y axis limits
    #                                         x_label                 =   'Horizontal coordinate (dva)',              # x axis label
    #                                         y_label                 =   'Vertical coordinate (dva)',                # y axis label
    #                                         x_tick_steps            =   2,                                          # x axis major tick steps
    #                                         y_tick_steps            =   2,                                          # y axis major tick steps
    #                                         dataMat                 =   data['pRF_fit_LR'],
    #                                         smooth_factor           =   15,                                         # pixel per degree in image
    #                                         cmap                    =   'viridis',                                  # colormap to use
    #                                         cmap_steps              =   10,                                         # colormap steps
    #                                         col_offset              =   0,                                          # colormap offset
    #                                         vmin                    =   0,                                          # val min to draw with colormap
    #                                         vmax                    =   1,                                          # val max to draw with colormap
    #                                         colorbar_tick_steps     =   0.2,                                        # colorbar axis major tick steps
    #                                         condition               =   'cov',
    #                                         colorbar_label          =   'pRF coverage (norm.)')                     # colorbar label
    #                                         )
    
    # for hemi, list_name, ext in [('-left', 'lh', '_L'), ('-right', 'rh', '_R'), ('', 'lrh', '_LR')]:
    #     title                           =   '{roi}{hemi} hemisphere: pRF coverage (n={voxel})'.\
    #                                          format(roi = roi_text, hemi = hemi, voxel = data['pRF_fit_L'].shape[0])            # define titlte
    #     _                               =   params_pRFcov.update(dict(              
    #                                             dataMat                 =   data['pRF_fit%s'%ext],                                  # main data matrix
    #                                             main_fig_title          =   title,
    #                                             hemisphere              =   list_name,
    #                                             typeData                =   None))
    #     f_pRFcov[list_name],_           =   plotter.draw_figure(parameters = params_pRFcov, plot = 'cov')



    # # Combining Figures
    # # -----------------
    # for list_name in ['lh', 'rh', 'lrh']:
    #     all_fL                          =     gridplot([
    #                                             [f_pRFmap[list_name],                f_pRFecc_basic[list_name][0],     f_pRFcov[list_name]],                # specify figure 1st row
    #                                             [f_pRFecc_basic[list_name][1],       f_pRFecc_basic[list_name][2],     None],                               # specify figure 2nd row
    #                                             [f_pRFecc_basic[list_name][3],       f_pRFecc_basic[list_name][4],     None]])
    #     folder_name                     =   opj(pRFbasic_folder, list_name.upper()) if list_name != 'lrh' else pRFbasic_folder
    #     output_file_html                =   opj(folder_name,'%s_pRF.html'%roi_text)                                                          # define html file
    #     _                               =   output_file(output_file_html, title='%s pRF analysis'%roi_text)                                   # get html saved
    #     _                               =   save(all_fL)



