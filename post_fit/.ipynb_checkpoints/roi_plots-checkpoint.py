"""
-----------------------------------------------------------------------------------------
roi_plots.py
-----------------------------------------------------------------------------------------
Goal of the script:
Draw roi plots (maps, ecc vs. params, laterality, time course)
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject number
sys.argv[2]: fit model ('gauss','css')
sys.argv[3]: do single hemifield plot (1 = YES, 0 = NO)
sys.argv[4]: save svg (1 = YES, 0 = NO)
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
source activate i27
cd /home/szinte/projects/retino_HCP
python post_fit/roi_plots.py sub-01 gauss_sg 1 0
python post_fit/roi_plots.py sub-02 gauss_sg 1 0
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
import scipy.io
opj = os.path.join
deb = ipdb.set_trace

# MRI imports
# -----------
import nibabel as nb
import cortex

# Functions import
# ----------------
from plot_class import PlotOperator
from utils import set_pycortex_config_file, mask_gii_2_hdf5 
from prf_utils import *

# Bokeh imports
# ---------------
from bokeh.io import output_notebook, show, save, output_file, export_png, export_svgs
from bokeh.layouts import row, column, gridplot


# import popeye.utilities as utils
# from popeye.visual_stimulus import VisualStimulus
# import popeye.css as css
# import popeye.og as og

# Get inputs
# ----------
subject = sys.argv[1]
fit_model = sys.argv[2]
draw_hemi = int(sys.argv[3])
save_svg = int(sys.argv[4])

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

deriv_dir = opj(base_dir,'pp_data',subject,fit_model,'deriv')
roi_masks_dir = opj(base_dir,'pp_data',subject,fit_model,'roi_masks')
h5_dir = opj(base_dir,'pp_data',subject,fit_model,'h5')
try: os.makedirs(roi_masks_dir)
except OSError: pass

# Create stimulus design and define model
# ---------------------------------------
visual_dm_file = scipy.io.loadmat(opj(base_dir,'raw_data','vis_design.mat'))
visual_dm = visual_dm_file['stim']

# Fit: define search grids
x_grid_bound = (-13.4, 13.4)
y_grid_bound = (-7.5, 7.5)
sigma_grid_bound = (0.05, 9.5)
n_grid_bound = (0.01, 1.5)

# Fit: define search bounds
x_fit_bound = (-30.0, 30.0)
y_fit_bound = (-30.0, 30.0)
sigma_fit_bound = (0.001, 70.0)
n_fit_bound = (0.01, 3)
beta_fit_bound = (-1e3, 1e3)
baseline_fit_bound = (-1e3, 1e3)

if fit_model == 'gauss' or fit_model == 'gauss_sg':
    bound_grids  = (x_grid_bound, y_grid_bound, sigma_grid_bound)
    bound_fits = (x_fit_bound, y_fit_bound, sigma_fit_bound, beta_fit_bound, baseline_fit_bound)
elif fit_model == 'css' or fit_model == 'css_sg':
    bound_grids  = (x_grid_bound, y_grid_bound, sigma_grid_bound, n_grid_bound)
    bound_fits = (x_fit_bound, y_fit_bound, sigma_fit_bound, n_fit_bound, beta_fit_bound, baseline_fit_bound)
    
# Initialize the prf model
prf = prf_fit(  fit_model = fit_model, 
                visual_design = visual_dm, 
                screen_distance = analysis_info["screen_distance"],
                screen_width = analysis_info["screen_width"],
                tr =  analysis_info["TR"],
                grid_steps = analysis_info["grid_steps"],
                hrf_delay = analysis_info["hrf_delay"],
                bound_grids = bound_grids,
                bound_fits = bound_fits,
                sg_filter_window_length = analysis_info["sg_filt_window_length"],
                sg_filter_polyorder = analysis_info["sg_filt_polyorder"],
                sg_filter_deriv = analysis_info["sg_filt_deriv"], 
                )

step_r2 = [0,100/3.0,250/3.0,100]
list_r2_level = ['High','Low']
step_params = [0,100/3.0,200/3.0,100]
list_params_level = ['High','Low']
list_params = ['ecc','amp','size','cov']
num_plots = len(list_params)*len(step_params)*(len(step_r2)-1)   

# Draw main analysis figure
# -------------------------
print('creating bokeh plots')
sign_idx, rsq_idx, ecc_idx, polar_real_idx, polar_imag_idx , size_idx, \
            non_lin_idx, amp_idx, baseline_idx, cov_idx, x_idx, y_idx = 0,1,2,3,4,5,6,7,8,9,10,11

# Initialize data dictionary that will save all data arrays
for roi_num, roi in enumerate(analysis_info['rois']):
    
    for mask_dir in ['pos','neg']:    
        data_hemi = []
        val_hemi = 0
        tc_mat_hemi = []
        
        for hemi in ['L', 'R', 'LR']:
            # create folder
            roi_text = roi
            exec('fig_bokeh_dir_{mask_dir}_{hemi} = opj(base_dir,"pp_data",subject,fit_model,"figs","prf","{mask_dir}","{hemi}")'.format(mask_dir=mask_dir, hemi = hemi))
            try: exec('os.makedirs(fig_bokeh_dir_{mask_dir}_{hemi})'.format(mask_dir=mask_dir,hemi = hemi))
            except: pass

            if save_svg == 1:
                # figure svg main folder
                exec('fig_bokeh_dir_svg_{mask_dir}_{hemi} = opj(base_dir,"pp_data",subject,fit_model,"figs","svg","{mask_dir}","{hemi}")'.format(mask_dir=mask_dir, hemi = hemi))
                exec('svg_folder = fig_bokeh_dir_svg_{mask_dir}_{hemi}'.format(mask_dir=mask_dir,hemi=hemi))
                try: os.makedirs(svg_folder)
                except: pass

            # load data
            if hemi == 'LR':
                data = np.row_stack((data_hemi[0],data_hemi[1]))
                tc_mat = np.row_stack((tc_mat_hemi[0],tc_mat_hemi[1]))
                draw = True
            
            else:
                # derivative data
                if hemi == 'L': val_hemi = 1
                elif hemi == 'R': val_hemi = 2
                folder_alias = '{hemi}_{mask_dir}'.format(hemi = hemi,mask_dir = mask_dir)
                h5_file = h5py.File(opj(h5_dir,'{roi}.h5'.format(roi = roi_text)), "r")
                in_file = opj("prf_deriv_{hemi}_{mask_dir}".format(hemi = hemi, mask_dir = mask_dir))
                data = h5_file['{folder_alias}/{in_file}'.format(folder_alias=folder_alias,in_file=in_file)]
                data = np.vstack((data,val_hemi*np.ones((1,data.shape[1]))))
                data = data[:,:].T
                data_hemi.append(data)
                if draw_hemi == 1:
                    draw = True
                elif draw_hemi == 0:
                    draw = False

                # time course data
                if hemi == 'L' or hemi == 'R':
                    mask_file = opj(roi_masks_dir,"masks_{hemi}_fsaverage6.gii".format(hemi = hemi))
                    gii_in_mask = nb.load(mask_file)
                    mask_mat = np.array([gii_in_mask.darrays[i].data for i in range(len(gii_in_mask.darrays))])

                    # load time course of original data
                    
                    orig_data_file  =  sorted(glob.glob(opj(base_dir,'raw_data',subject,"{sub}_task-prf_space-fsaverage6_hemi-{hemi}.func_sg_psc.gii".format(sub = subject, hemi = hemi))))
                    orig_data_file_load = nb.load(orig_data_file[0])
                    orig_data = []
                    orig_data.append(np.array([orig_data_file_load.darrays[i].data for i in range(len(orig_data_file_load.darrays))]))
                    orig_data = np.vstack(orig_data)

                # select corresponding time course
                mask_mat_roi = mask_mat[roi_num,:]
                mask_mat_roi = np.round(mask_mat_roi)
                tc_mat = orig_data[:, mask_mat_roi==1].T
                tc_mat_hemi.append(tc_mat)

            vertex_ini = data.shape[0]

            if draw == True:
                
                if vertex_ini > 0:
                    data4mask = data
                    data = data[np.logical_and(np.logical_and( data4mask[:,rsq_idx]>=analysis_info['rsq_threshold'],
                                                                data4mask[:,cov_idx]>=analysis_info['cov_threshold']),
                                                                data4mask[:,size_idx]>=analysis_info['size_threshold'])]

                    tc_mat = tc_mat[np.logical_and(np.logical_and(  data4mask[:,rsq_idx]>=analysis_info['rsq_threshold'],
                                                                    data4mask[:,cov_idx]>=analysis_info['cov_threshold']),
                                                                    data4mask[:,size_idx]>=analysis_info['size_threshold'])]

                    vertex = data.shape[0]

                    
                    if vertex > 0:

                        print("drawing {roi}_{hemi}_{mask_dir} figures, n={vertex}".format(roi = roi_text,hemi = hemi,vertex = vertex, mask_dir = mask_dir)) 

                        data_source = { 'sign':             data[:,sign_idx],
                                        'rsq':              data[:,rsq_idx],
                                        'ecc':              data[:,ecc_idx],
                                        'sigma':            data[:,size_idx],
                                        'non_lin':          data[:,non_lin_idx],
                                        'beta':             data[:,amp_idx],
                                        'baseline':         data[:,baseline_idx],
                                        'cov':              data[:,cov_idx],
                                        'x':                data[:,x_idx],
                                        'y':                data[:,y_idx],
                                        'colors_ref':       data[:,ecc_idx]}

                        param_all = {   'roi_t':            roi_text, 
                                        'p_width':          400, 
                                        'p_height':         400, 
                                        'min_border_large': 10, 
                                        'min_border_small': 5,
                                        'bg_color':         tuple([229,229,229]), 
                                        'stim_color':       tuple([250,250,250]), 
                                        'hist_fill_color':  tuple([255,255,255]),
                                        'hist_line_color':  tuple([0,0,0]), 
                                        'stim_width':       analysis_info['stim_width'], 
                                        'stim_height':      analysis_info['stim_height'], 
                                        'rsq_threshold':    analysis_info['rsq_threshold'],
                                        'cmap':             'Spectral',
                                        'cmap_steps':       9,
                                        'col_offset':       0,
                                        'vmin':             0,
                                        'vmax':             9,
                                        'leg_xy_max_ratio': 1.8,
                                        'dataMat':          data,
                                        'data_source':      data_source,
                                        'hist_range':       (0,0.5),
                                        'hist_steps':       0.5,
                                        'h_hist_bins':      16,
                                        'link_x':           False,
                                        'link_y':           False,
                                        'save_svg':         save_svg
                                        }

                        if save_svg == 1:
                            param_all.update({
                                        'svg_folder':       svg_folder})

                        plotter = PlotOperator(**param_all)

                        # pRFmap
                        title = '{roi}{hemi} {mask_dir}: pRF map (n={vertex})'.format(roi = roi_text, hemi = hemi, vertex = vertex, mask_dir = mask_dir)
                        params_pRFmap = param_all
                        params_pRFmap.update(\
                                    {   'x_range':          (-12, 12),
                                        'y_range':          (-12, 12),
                                        'x_label':          'Horizontal coordinate (dva)',
                                        'y_label':          'Vertical coordinate (dva)',
                                        'x_source_label':   'x',
                                        'y_source_label':   'y',
                                        'x_tick_steps':     4,
                                        'y_tick_steps':     4,
                                        'v_hist_bins':      16,
                                        'main_fig_title':   title})
                        if save_svg == 1:
                            params_pRFmap.update(
                                    {   'svg_filename':     '{roi}_{hemi}_{mask_dir}_pRFmap'.format(mask_dir=mask_dir,hemi=hemi,roi= roi_text)})

                        f_pRFmap,old_main_fig1 = plotter.draw_figure(   parameters = params_pRFmap, 
                                                                        plot = 'map')

                        
                        # pRFecc
                        old_main_fig = []
                        f_pRFecc = []
                        if fit_model == 'gauss' or fit_model == 'gauss_sg':
                            type_comp_list = ['Size','R2','Amplitude','Coverage','Baseline']
                        elif fit_model == 'css' or fit_model == 'css_sg':
                            type_comp_list = ['Size','R2','Non-Linearity','Amplitude','Coverage','Baseline']

                        for numData, type_comp in enumerate(type_comp_list):

                            params_pRFecc = param_all
                            params_pRFecc.update(   
                                       {    'x_range':          (0, 10),
                                            'x_label':          'Eccentricity (dva)',
                                            'x_tick_steps':     2,
                                            'x_source_label':   'ecc',
                                            'draw_reg':         False,
                                            'link_x':           True})
                            if save_svg == 1:
                                params_pRFecc.update(
                                            {   'svg_subfolder':    '{roi}_{hemi}_{mask_dir}_pRFecc'.format(mask_dir=mask_dir,hemi=hemi,roi= roi_text)})

                            if type_comp == 'Size':
                                params_pRFecc.update(
                                            {   'y_range':          (0, 10),
                                                'y_label':          'Size (dva)',
                                                'y_source_label':   'sigma',
                                                'y_tick_steps':     2,
                                                'v_hist_bins':      20,
                                                'draw_reg':         True})
                                if save_svg == 1:
                                    params_pRFecc.update(
                                            {   'svg_filename':     '{roi}_{hemi}_{mask_dir}_pRFecc_size'.format(mask_dir=mask_dir,hemi=hemi,roi= roi_text)})

                            elif type_comp == 'R2':
                                params_pRFecc.update(
                                            {   'y_range':          (0, 1),
                                                'y_label':          'R2',
                                                'y_source_label':   'rsq',
                                                'y_tick_steps':     0.2,
                                                'v_hist_bins':      15})
                                if save_svg == 1:
                                    params_pRFecc.update(
                                            {   'svg_filename':     '{roi}_{hemi}_{mask_dir}_pRFecc_rsq'.format(mask_dir=mask_dir,hemi=hemi,roi= roi_text)})

                            elif type_comp == 'Non-Linearity':
                                params_pRFecc.update(
                                            {   'y_range':          (0, 1.5),
                                                'y_label':          'Non-linearity',
                                                'y_source_label':   'non_lin',
                                                'y_tick_steps':     0.25,
                                                'v_hist_bins':      18})
                                if save_svg == 1:
                                    params_pRFecc.update(
                                            {   'svg_filename':     '{roi}_{hemi}_{mask_dir}_pRFecc_non_lin'.format(mask_dir=mask_dir,hemi=hemi,roi= roi_text)})

                            elif type_comp == 'Amplitude':
                                params_pRFecc.update(
                                            {   'y_range':          (-1, 1),
                                                'y_label':          'Amplitude (z-score)',
                                                "y_source_label":   'beta',
                                                'y_tick_steps':     0.5,
                                                'v_hist_bins':      16})
                                if save_svg == 1:
                                    params_pRFecc.update(
                                            {   'svg_filename':     '{roi}_{hemi}_{mask_dir}_pRFecc_amp'.format(mask_dir=mask_dir,hemi=hemi,roi= roi_text)})

                            elif type_comp == 'Coverage':
                                params_pRFecc.update(
                                            {   'y_range':          (0, 1),
                                                'y_label':          'pRF coverage (%)',
                                                'y_source_label':   'cov',
                                                'y_tick_steps':     0.2,
                                                'v_hist_bins':      15})
                                if save_svg == 1:
                                    params_pRFecc.update(
                                            {   'svg_filename':     '{roi}_{hemi}_{mask_dir}_pRFecc_cov'.format(mask_dir=mask_dir,hemi=hemi,roi= roi_text)})

                            elif type_comp == 'Baseline':
                                params_pRFecc.update(
                                            {   'y_range':          (-1, 1),
                                                'y_label':          'Baseline (z-score)',
                                                'y_source_label':   'baseline',
                                                'y_tick_steps':     0.5,
                                                'v_hist_bins':      16})
                                if save_svg == 1:
                                    params_pRFecc.update(
                                            {   'svg_filename':     '{roi}_{hemi}_{mask_dir}_pRFecc_baseline'.format(mask_dir=mask_dir,hemi=hemi,roi= roi_text)})

                            title = '{roi}{hemi} {mask_dir}: Eccentricity vs. {type_comp}'.format(roi = roi_text, hemi = hemi, type_comp = type_comp, mask_dir = mask_dir)
                            params_pRFecc.update({'main_fig_title':   title})

                            out1, old_main_fig  = plotter.draw_figure(  parameters = params_pRFecc,
                                                                        plot = 'ecc',
                                                                        old_main_fig = old_main_fig)
                            f_pRFecc.append(out1)

                        # pRF cov
                        # -------
                        params_pRFcov = param_all
                        params_pRFcov.update(
                                    {   'dataMat':          data,
                                        'x_range':          (-12, 12), 
                                        'y_range':          (-12, 12),
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
                                        'cb_label':         'pRF coverage (norm.)',
                                        'link_x':           True,
                                        'link_y':           True})
                        if save_svg == 1:
                            params_pRFcov.update(
                                    {   'svg_filename':     '{roi}_{hemi}_{mask_dir}_pRFcov'.format(mask_dir=mask_dir,hemi=hemi,roi= roi_text)})
                        title = '{roi}{hemi} {mask_dir}: pRF density map'.format(roi = roi_text, hemi = hemi, mask_dir = mask_dir)
                        params_pRFcov.update({'main_fig_title':   title})
                        f_pRFcov = plotter.draw_figure(parameters = params_pRFcov, plot = 'cov',old_main_fig = old_main_fig1)

                        
                        # pRF lat
                        # -------
                        params_pRFlat = param_all
                        params_pRFlat.update(
                                    {   'p_width':          500, 
                                        'p_height':         500,
                                        'dataMat':          data,
                                        'x_range':          (-2.6, 2.6), 
                                        'y_range':          (-2.6, 2.6),
                                        'vmin':             0,
                                        'vmax':             0.2,
                                        'weighted_data':    True,
                                        'main_fig_title':   title,
                                        'cmap':             'hsv',
                                        'cmap_steps':       16,
                                        'ang_bins':         36,
                                        'hemi':             hemi})
                        if save_svg == 1:
                            params_pRFlat.update(
                                    {   'svg_filename':     '{roi}_{hemi}_{mask_dir}_pRFlat'.format(mask_dir=mask_dir,hemi=hemi,roi= roi_text)})
                        title = '{roi}{hemi} {mask_dir}: pRF laterality histogram'.format(roi = roi_text, hemi = hemi, type_comp = type_comp, mask_dir = mask_dir)
                        params_pRFcov.update({'main_fig_title':   title})
                        f_pRFlat = plotter.draw_figure(parameters = params_pRFlat, plot = 'lat')

                        # pRF tc
                        # -------
                        params_pRFtc = param_all
                        params_pRFtc.update(
                                    {   'p_width':          500, 
                                        'p_height':         500,
                                        'x_range_map':      (-12,12),
                                        'y_range_map':      (-12,12),
                                        'x_label_map':      'Horizontal coord. (dva)',
                                        'y_label_map':      'Vertical coord. (dva)',
                                        'x_tick_map':       4,
                                        'y_tick_map':       4,
                                        'x_range_tc':       (0,analysis_info['n_timepoints_per_run']*analysis_info['TR']),
                                        'x_label_tc':       'Time (s)',
                                        'y_label_tc':       'BOLD signal change (%)',
                                        'x_tick_tc':        20,
                                        'tr_dur':           analysis_info['TR'],
                                        'model_line_color': tuple([254,51,10]),
                                        'model_fill_color': tuple([254,51,10])
                                    })

                        if save_svg == 1:
                            params_pRFtc.update(
                                        {   'svg_subfolder':    '{roi}_{hemi}_{mask_dir}_pRFtc'.format(mask_dir=mask_dir,hemi=hemi,roi= roi_text)})

                        f_pRFtc = []

                        # get index matrices
                        prct_r2 = np.nanpercentile(a = data[:,rsq_idx],q = step_r2)
                        idx_r2 = []
                        for r2_step in np.arange(0,len(step_r2)-1,2):

                            idx_r2.append(np.logical_and(data[:,rsq_idx]>=prct_r2[r2_step],data[:,rsq_idx]<=prct_r2[r2_step+1]))


                        # get vertex number to draw
                        for param in list_params:
                            exec('prct_{param} = np.nanpercentile(a = data[:,{param}_idx],q = step_params)'.format(param = param))
                            exec('idx_{param} = []'.format(param = param))

                            for param_step in np.arange(0,len(step_params)-1,2):
                                exec('idx_{param}.append(np.logical_and(data[:,{param}_idx]>=prct_{param}[param_step],data[:,{param}_idx]<=prct_{param}[param_step+1]))'.format(param = param))

                            exec('num_{param} = []'.format(param = param))
                            exec('num_{param} = []'.format(param = param))
                            for r2_step in np.arange(0,2,1):

                                for param_step in np.arange(0,2,1):

                                    exec('mat = np.where(np.logical_and(idx_r2[r2_step],idx_{param}[param_step]))'.format(param = param))
                                    if mat[0].size == 0:
                                        exec('num_{param}.append(-1)'.format(param = param))
                                    else:
                                        exec('num_{param}.append(mat[0][np.random.choice(len(mat[0]))])'.format(param = param))

                        for param in list_params:
                            exec('num_vertex = num_{param}'.format(param = param))

                            for r2_level in list_r2_level:
                                if r2_level == 'Low':
                                    num_vertex2draw = num_vertex[0:2]
                                elif r2_level == 'High':
                                    num_vertex2draw = num_vertex[2:4]

                                params_pRFtc.update(   
                                                {   'params':               param,
                                                    'r2_level':             r2_level,
                                                    'deriv_mat':            data,
                                                    'tc_mat':               tc_mat,
                                                    'num_vertex':           num_vertex2draw,
                                                    'fit_model':            fit_model,
                                                    'model_func':           prf,
                                                    'mask_dir':             mask_dir,
                                                    'title':                '{roi} {hemi} {sign}'.format(sign = mask_dir,hemi=hemi,roi = roi_text)
                                                })
                                if save_svg == 1:
                                    params_pRFtc.update(
                                                    {   'svg_filename':     '{roi}_{hemi}_{mask_dir}_pRFtc_rsq_{r2_level}_{param}'.format(mask_dir=mask_dir,hemi=hemi,roi= roi_text,r2_level = r2_level, param = param)})

                                out4,main_fig4  = plotter.draw_figure(    parameters =    params_pRFtc,
                                                                plot =          'tc')

                                f_pRFtc.append(out4)

                        # save files
                        if fit_model == 'gauss' or fit_model == 'gauss_sg':
                            all_f1 = gridplot([ [f_pRFecc[0],f_pRFecc[1],f_pRFecc[2]],
                                                [f_pRFecc[3],f_pRFecc[4],None]],toolbar_location='right')

                            all_f4 = gridplot([ [f_pRFtc[0],f_pRFtc[1]],
                                                [f_pRFtc[2],f_pRFtc[3]],
                                                [f_pRFtc[4],f_pRFtc[5]],
                                                [f_pRFtc[6],f_pRFtc[7]]],toolbar_location = None)
                        elif fit_model == 'css' or fit_model == 'css_sg':
                            all_f1 = gridplot([ [f_pRFecc[0],f_pRFecc[1],f_pRFecc[2]],
                                                [f_pRFecc[3],f_pRFecc[4],f_pRFecc[5]]],toolbar_location='right')
                            
                            all_f4 = gridplot([ [f_pRFtc[0],f_pRFtc[1]],
                                                [f_pRFtc[2],f_pRFtc[3]],
                                                [f_pRFtc[2],f_pRFtc[3]],
                                                [f_pRFtc[4],f_pRFtc[5]]])

                        exec('output_file_html = opj(fig_bokeh_dir_{mask_dir}_{hemi},"{roi_text}_{hemi}_{mask_dir}_pRFecc.html")'.format(mask_dir = mask_dir,roi_text = roi_text, hemi = hemi))
                        output_file(output_file_html, title='Subject: %s | ROI: %s | Hemisphere: %s | Sign: %s | Vertex: %i | Figures: pRF parameters and density'%(subject,roi_text,hemi,mask_dir,vertex))
                        save(all_f1)

                        all_f2 = gridplot([ [f_pRFmap,f_pRFcov[0]]],toolbar_location='right')
                        exec('output_file_html = opj(fig_bokeh_dir_{mask_dir}_{hemi},"{roi_text}_{hemi}_{mask_dir}_pRFmap.html")'.format(mask_dir = mask_dir,roi_text = roi_text, hemi = hemi))
                        output_file(output_file_html, title='Subject: %s | ROI: %s | Hemisphere: %s | Sign: %s | Vertex: %i | Figures:pRF maps parameters and density'%(subject,roi_text,hemi,mask_dir,vertex))
                        save(all_f2)
                        
                        

                        all_f3 = gridplot([[f_pRFlat[0]]],toolbar_location = 'right')
                        exec('output_file_html = opj(fig_bokeh_dir_{mask_dir}_{hemi},"{roi_text}_{hemi}_{mask_dir}_pRFlat.html")'.format(mask_dir = mask_dir,roi_text = roi_text, hemi = hemi))
                        output_file(output_file_html, title='Subject: %s | ROI: %s | Hemisphere: %s | Sign: %s | Vertex: %i | Figures:pRF laterality histogram'%(subject,roi_text,hemi,mask_dir,vertex))
                        save(all_f3)

                        exec('output_file_html = opj(fig_bokeh_dir_{mask_dir}_{hemi},"{roi_text}_{hemi}_{mask_dir}_pRFtc.html")'.format(mask_dir = mask_dir,roi_text = roi_text, hemi = hemi))
                        output_file(output_file_html, title='Subject: %s | ROI: %s | Hemisphere: %s | Sign: %s | Figures: pRF time course'%(subject,roi_text,hemi,mask_dir))
                        save(all_f4)
                        

                    else:
                        print("drawing {roi}_{hemi}_{mask_dir} figures not possible: n={vertex}".format(roi = roi_text,hemi = hemi,vertex = vertex,mask_dir = mask_dir)) 

