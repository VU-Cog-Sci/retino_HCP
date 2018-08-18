"""
-----------------------------------------------------------------------------------------
post_pp_roi.py
-----------------------------------------------------------------------------------------
Goal of the script:
Draw basic plot of pRF analysis
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject number
sys.argv[2]: fit model ('gauss','css')
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
source activate i36
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
fit_model = sys.argv[2]

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
deriv_dir = opj(base_dir,'pp_data',subject,fit_model,'deriv')
h5_dir = opj(base_dir,'pp_data',subject,fit_model,'h5')

# Draw main analysis figure
# -------------------------
print('creating bokeh plots')

# Initialize data dictionary that will save all data arrays
for roi in analysis_info['rois']:
    
    for mask_dir in ['all','pos','neg']:

        data_hemi = []
        for hemi in ['L', 'R', 'LR']:

            # create folders
            roi_text = analysis_info['rois'][roi]
            exec('fig_bokeh_dir_{mask_dir}_{hemi} = opj(base_dir,"pp_data",subject,fit_model,"figs","prf","{mask_dir}","{hemi}")'.format(mask_dir=mask_dir, hemi = hemi))
            try: exec('os.makedirs(fig_bokeh_dir_{mask_dir}_{hemi})'.format(mask_dir=mask_dir,hemi = hemi))
            except: pass
            
            # load data
            if hemi == 'LR':
                data = np.row_stack((data_hemi[0],data_hemi[1]))
            else:
                folder_alias = '{hemi}_{mask_dir}'.format(hemi = hemi,mask_dir = mask_dir)
                h5_file = h5py.File(opj(h5_dir,'{roi}.h5'.format(roi = analysis_info['rois'][roi])), "r")
                in_file = opj("prf_deriv_{hemi}_{mask_dir}".format(hemi = hemi, mask_dir = mask_dir))
                data = h5_file['{folder_alias}/{in_file}'.format(folder_alias=folder_alias,in_file=in_file)]
                data = data[:,:].T
                data_hemi.append(data)

            vertex_ini = data.shape[0]
            if vertex_ini != 1:
                data = data[data[:,1]>=analysis_info['rsq_threshold'],:]
                data = data[data[:,9]>=analysis_info['cov_threshold'],:]
                vertex = data.shape[0]
                

                print("drawing {roi}_{hemi}_{mask_dir} figures, n={vertex}".format(roi = analysis_info['rois'][roi],hemi = hemi,vertex = vertex, mask_dir = mask_dir)) 
                
                data_source = { 'sign':             data[:,0],
                                'rsq':              data[:,1],
                                'ecc':              data[:,2],
                                'sigma':            data[:,5],
                                'non_lin':          data[:,6],
                                'beta':             data[:,7],
                                'baseline':         data[:,8],
                                'cov':              data[:,9],
                                'x':                data[:,10],
                                'y':                data[:,11],
                                'colors_ref':       data[:, 2]}

                
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
                                'saving_figs':      False,
                                'svg_folder':       'fig_bokeh_dir_{mask_dir}'.format(mask_dir=mask_dir),
                                'dataMat':          data,
                                'data_source':      data_source,
                                'hist_range':       (0,0.5),
                                'hist_steps':       0.5,
                                'h_hist_bins':      16,
                                }

                plotter = PlotOperator(**param_all)

                # pRFmap
                title = '{roi}_{hemi}_{mask_dir}: pRF map (n={vertex})'.format(roi = roi_text, hemi = hemi, vertex = vertex, mask_dir = mask_dir)
                params_pRFmap = param_all
                params_pRFmap.update(\
                            {   'x_range':          (-16, 16),
                                'y_range':          (-16, 16),
                                'x_label':          'Horizontal coordinate (dva)',
                                'y_label':          'Vertical coordinate (dva)',
                                'x_source_label':   'x',
                                'y_source_label':   'y',
                                'x_tick_steps':     4,
                                'y_tick_steps':     4,
                                'v_hist_bins':      16,
                                'main_fig_title':   title})

                f_pRFmap = plotter.draw_figure(parameters = params_pRFmap, 
                                                            plot = 'map')
                
                # pRFecc
                old_main_fig = []
                f_pRFecc = []
                if fit_model == 'gauss':
                    type_comp_list = ['Size','R2','Amplitude','Coverage','Baseline']
                elif fit_model == 'css':
                    type_comp_list = ['Size','R2','Non-Linearity','Amplitude','Coverage','Baseline']

                for numData, type_comp in enumerate(type_comp_list):

                    params_pRFecc = param_all
                    params_pRFecc.update(   
                               {    'x_range':          (0, 16),
                                    'x_label':          'Eccentricity (dva)',
                                    'x_tick_steps':     2,
                                    'x_source_label':   'ecc'})

                    if type_comp == 'Size':
                        params_pRFecc.update(
                                    {   'y_range':          (0, 8),
                                        'y_label':          'Size (dva)',
                                        'y_source_label':   'sigma',
                                        'y_tick_steps':     1,
                                        'v_hist_bins':      16})

                    elif type_comp == 'R2':
                        params_pRFecc.update(
                                    {   'y_range':          (0, 1),
                                        'y_label':          'R2',
                                        'y_source_label':   'rsq',
                                        'y_tick_steps':     0.2,
                                        'v_hist_bins':      15})

                    elif type_comp == 'Non-Linearity':
                        params_pRFecc.update(
                                    {   'y_range':          (0, 1.5),
                                        'y_label':          'Non-linearity',
                                        'y_source_label':   'non_lin',
                                        'y_tick_steps':     0.25,
                                        'v_hist_bins':      18})

                    elif type_comp == 'Amplitude':
                        params_pRFecc.update(
                                    {   'y_range':          (-1, 1),
                                        'y_label':          'Amplitude (z-score)',
                                        "y_source_label":   'beta',
                                        'y_tick_steps':     0.5,
                                        'v_hist_bins':      16})

                    elif type_comp == 'Coverage':
                        params_pRFecc.update(
                                    {   'y_range':          (0, 1),
                                        'y_label':          'pRF coverage (%)',
                                        'y_source_label':   'cov',
                                        'y_tick_steps':     0.2,
                                        'v_hist_bins':      15})
                    elif type_comp == 'Baseline':
                        params_pRFecc.update(
                                    {   'y_range':          (-1, 1),
                                        'y_label':          'Baseline (z-score)',
                                        'y_source_label':   'baseline',
                                        'y_tick_steps':     0.5,
                                        'v_hist_bins':      16})

                    title = '{roi}_{hemi}_{mask_dir}: Eccentricity vs. {type_comp}'.format(roi = roi_text, hemi = hemi, type_comp = type_comp, mask_dir = mask_dir)
                    params_pRFecc.update({'main_fig_title':   title})

                    out1, old_main_fig  = plotter.draw_figure(  parameters = params_pRFecc,
                                                                plot = 'ecc',
                                                                old_main_fig = old_main_fig)
                    f_pRFecc.append(out1)
                

                # # pRF cov
                # # -------
                params_pRFcov = param_all
                params_pRFcov.update(
                            {   'dataMat':          data,
                                'x_range':          (-16, 16), 
                                'y_range':          (-16, 16),
                                'x_label':          'Horizontal coordinate (dva)',
                                'y_label':          'Vertical coordinate (dva)',
                                'x_tick_steps':     4,
                                'y_tick_steps':     4,
                                'smooth_factor':    15,
                                'cmap':             'viridis',
                                'cmap_steps':       10,
                                'col_offset':       0,
                                'vmin':             0,
                                'vmax':             1,
                                'cb_tick_steps':    0.2,
                                'condition':        'cov',
                                'cb_label':         'pRF coverage (norm.)'})

                f_pRFcov = plotter.draw_figure(parameters = params_pRFcov, plot = 'cov')

                # save files
                if fit_model == 'gauss':
                    all_f1 = gridplot([ [f_pRFecc[0],f_pRFecc[1],f_pRFecc[2]],
                                        [f_pRFecc[3],f_pRFecc[4],None]])
                elif fit_model == 'css':
                    all_f1 = gridplot([ [f_pRFecc[0],f_pRFecc[1],f_pRFecc[2]],
                                        [f_pRFecc[3],f_pRFecc[4],f_pRFecc[5]]])
                
                exec('output_file_html = opj(fig_bokeh_dir_{mask_dir}_{hemi},"{roi_text}_{hemi}_{mask_dir}_pRFecc.html")'.format(mask_dir = mask_dir,roi_text = roi_text, hemi = hemi))
                output_file(output_file_html, title='%s pRF analysis'%roi_text)
                save(all_f1)

                all_f2 = gridplot([ [f_pRFmap[0],f_pRFcov[0]]])
                exec('output_file_html = opj(fig_bokeh_dir_{mask_dir}_{hemi},"{roi_text}_{hemi}__{mask_dir}_pRFmap.html")'.format(mask_dir = mask_dir,roi_text = roi_text, hemi = hemi))
                output_file(output_file_html, title='%s pRF analysis'%roi_text)
                save(all_f2)        
            else:
                print("drawing of {roi}_{hemi}_{mask_dir} figures is not possible as n={vertex}".format(roi = analysis_info['rois'][roi],hemi = hemi,vertex = vertex,mask_dir = mask_dir)) 

