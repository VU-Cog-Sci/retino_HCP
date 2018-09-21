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
sys.argv[3]: do single hemifield plot (1 = YES, 0 = NO)
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
source activate i27
cd /home/szinte/projects/retino_HCP
python post_fit/post_pp_roi.py 999999 css 1
python post_fit/post_pp_roi.py 999999 gauss 1
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
from utils import set_pycortex_config_file, mask_gii_2_hdf5 

# Bokeh imports
# ---------------
from bokeh.io import output_notebook, show, save, output_file, export_png, export_svgs
from bokeh.layouts import row, column, gridplot

# Get inputs
# ----------
subject = sys.argv[1]
fit_model = sys.argv[2]
draw_hemi = int(sys.argv[3])

# Define analysis parameters
# --------------------------
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define cluster/server specific parameters
# -----------------------------------------
if 'aeneas' in platform.uname()[1]:
    base_dir = analysis_info['aeneas_base_folder'] 
    main_cmd = '/home/szinte/software/workbench/bin_rh_linux64/wb_command'
elif 'local' in platform.uname()[1]:
    base_dir = analysis_info['local_base_folder'] 
    main_cmd = '/Applications/workbench/bin_macosx64/wb_command'

deriv_dir = opj(base_dir,'pp_data',subject,fit_model,'deriv')
roi_masks_dir = opj(base_dir,'pp_data',subject,fit_model,'roi_masks')
h5_dir = opj(base_dir,'pp_data',subject,fit_model,'h5')
try: os.makedirs(roi_masks_dir)
except OSError: pass

# Check system
# ------------
sys.exit('Drawing Flatmaps only works with Python 2. Aborting.') if sys.version_info[0] > 2 else None

# Change cortex database folder
# -----------------------------
pycortex_folder     =   opj(base_dir,'pp_data','cortex')
set_pycortex_config_file(   project_folder = pycortex_folder)

# Create mask from overlay.svg
# ----------------------------
print('creating roi masks from overlay.svg')
masks = cortex.utils.get_roi_verts( subject = 'fsaverage', 
                                    roi = analysis_info['rois'], 
                                    mask = True)
mat_masks = []
for roi in analysis_info['rois']:
    mat_masks.append(masks[roi])
mat_masks = np.vstack(mat_masks)
mat_masks = mat_masks.astype('float32')

prf_deriv_L_all_fsaverage = nb.load(opj(deriv_dir,'all','prf_deriv_L_all_fsaverage.func.gii'))
mat_masks_L = mat_masks[:,0:163842]
darrays = [nb.gifti.gifti.GiftiDataArray(d) for d in mat_masks_L]
gii_out = nb.gifti.gifti.GiftiImage(header = prf_deriv_L_all_fsaverage.header,
                                    extra = prf_deriv_L_all_fsaverage.extra,
                                    darrays = darrays)
nb.save(gii_out,opj(roi_masks_dir,"masks_L_fsaverage.func.gii"))

prf_deriv_R_all_fsaverage = nb.load(opj(deriv_dir,'all','prf_deriv_R_all_fsaverage.func.gii'))
mat_masks_R = mat_masks[:,163842:327684]
darrays = [nb.gifti.gifti.GiftiDataArray(d) for d in mat_masks_R]
gii_out = nb.gifti.gifti.GiftiImage(header = prf_deriv_R_all_fsaverage.header, 
                                    extra = prf_deriv_R_all_fsaverage.extra, 
                                    darrays = darrays)
nb.save(gii_out,opj(roi_masks_dir,"masks_R_fsaverage.func.gii"))

resample_cmd = """{main_cmd} -metric-resample {metric_in} {current_sphere} {new_sphere} ADAP_BARY_AREA {metric_out} -area-metrics {current_area} {new_area}"""
for hemi in ['L','R']:

    current_sphere = opj(base_dir,'raw_data/surfaces/resample_fsaverage','fsaverage_std_sphere.{hemi}.164k_fsavg_{hemi}.surf.gii'.format(hemi=hemi))
    new_sphere = opj(base_dir,'raw_data/surfaces/resample_fsaverage','fs_LR-deformed_to-fsaverage.{hemi}.sphere.32k_fs_LR.surf.gii'.format(hemi=hemi))
    current_area = opj(base_dir,'raw_data/surfaces/resample_fsaverage','fsaverage.{hemi}.midthickness_va_avg.164k_fsavg_{hemi}.shape.gii'.format(hemi=hemi))
    new_area = opj(base_dir,'raw_data/surfaces/resample_fsaverage','fs_LR.{hemi}.midthickness_va_avg.32k_fs_LR.shape.gii'.format(hemi=hemi))

    metric_in = opj(roi_masks_dir,"masks_{hemi}_fsaverage.func.gii".format(hemi = hemi))
    metric_out = opj(roi_masks_dir,"masks_{hemi}.func.gii".format(hemi = hemi))

    os.system(resample_cmd.format(  main_cmd = main_cmd,
                                    metric_in = metric_in, 
                                    current_sphere = current_sphere, 
                                    new_sphere = new_sphere, 
                                    metric_out = metric_out, 
                                    current_area = current_area, 
                                    new_area = new_area))

# Save ROIS data in hdf5
# ----------------------
print('creating h5 files')
for roi_num, roi in enumerate(analysis_info['rois']):
    try: os.makedirs(h5_dir)
    except OSError: pass

    h5_file = opj(h5_dir,'{roi}.h5'.format(roi = roi))
    try: os.system('rm '+ h5_file)
    except: pass

    for hemi in ['L','R']:

        mask_file = opj(roi_masks_dir,"masks_{hemi}.func.gii".format(hemi = hemi))
        
        for mask_dir in ['all','pos','neg']:
            
            in_file = opj(deriv_dir,mask_dir,"prf_deriv_{hemi}_{mask_dir}.gii".format(hemi = hemi, mask_dir = mask_dir))
            folder_alias = '{hemi}_{mask_dir}'.format(hemi = hemi,mask_dir = mask_dir)
            
            mask_gii_2_hdf5(in_file = in_file,
                            mask_file = mask_file,
                            hdf5_file = h5_file,
                            folder_alias = folder_alias,
                            roi_num = roi_num)


deb()
# Draw main analysis figure
# -------------------------
print('creating bokeh plots')
sign_idx, rsq_idx, ecc_idx, polar_real_idx, polar_imag_idx , size_idx, \
            non_lin_idx, amp_idx, baseline_idx, cov_idx, x_idx, y_idx = 0,1,2,3,4,5,6,7,8,9,10,11

# Initialize data dictionary that will save all data arrays
for roi in analysis_info['rois']:
    for mask_dir in ['all','pos','neg']:
        data_hemi = []
        val_hemi = 0
        for hemi in ['L', 'R', 'LR']:

            # create folder
            roi_text = roi
            exec('fig_bokeh_dir_{mask_dir}_{hemi} = opj(base_dir,"pp_data",subject,fit_model,"figs","prf","{mask_dir}","{hemi}")'.format(mask_dir=mask_dir, hemi = hemi))
            try: exec('os.makedirs(fig_bokeh_dir_{mask_dir}_{hemi})'.format(mask_dir=mask_dir,hemi = hemi))
            except: pass
            
            # load data
            if hemi == 'LR':
                data = np.row_stack((data_hemi[0],data_hemi[1]))
                draw = True
            else:

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

            vertex_ini = data.shape[0]

            if draw == True:
                if vertex_ini > 0:
                    data = data[data[:,rsq_idx]>=analysis_info['rsq_threshold'],:]
                    data = data[data[:,cov_idx]>=analysis_info['cov_threshold'],:]
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
                                        'link_x':           False,
                                        'link_y':           False,
                                        }

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

                        f_pRFmap,old_main_fig1 = plotter.draw_figure(parameters = params_pRFmap, 
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
                                       {    'x_range':          (0, 10),
                                            'x_label':          'Eccentricity (dva)',
                                            'x_tick_steps':     2,
                                            'x_source_label':   'ecc',
                                            'draw_reg':         False,
                                            'link_x':           True})

                            if type_comp == 'Size':
                                params_pRFecc.update(
                                            {   'y_range':          (0, 10),
                                                'y_label':          'Size (dva)',
                                                'y_source_label':   'sigma',
                                                'y_tick_steps':     2,
                                                'v_hist_bins':      20,
                                                'draw_reg':         True})

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
                        title = '{roi}{hemi} {mask_dir}: pRF density map'.format(roi = roi_text, hemi = hemi, type_comp = type_comp, mask_dir = mask_dir)
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
                        title = '{roi}{hemi} {mask_dir}: pRF laterality histogram'.format(roi = roi_text, hemi = hemi, type_comp = type_comp, mask_dir = mask_dir)
                        params_pRFcov.update({'main_fig_title':   title})
                        f_pRFlat = plotter.draw_figure(parameters = params_pRFlat, plot = 'lat')


                        # save files
                        if fit_model == 'gauss':
                            all_f1 = gridplot([ [f_pRFecc[0],f_pRFecc[1],f_pRFecc[2]],
                                                [f_pRFecc[3],f_pRFecc[4],None]],toolbar_location='right')
                        elif fit_model == 'css':
                            all_f1 = gridplot([ [f_pRFecc[0],f_pRFecc[1],f_pRFecc[2]],
                                                [f_pRFecc[3],f_pRFecc[4],f_pRFecc[5]]],toolbar_location='right')
                        
                        exec('output_file_html = opj(fig_bokeh_dir_{mask_dir}_{hemi},"{roi_text}_{hemi}_{mask_dir}_pRFecc.html")'.format(mask_dir = mask_dir,roi_text = roi_text, hemi = hemi))
                        output_file(output_file_html, title='Subject: %s | ROI: %s | Hemisphere: %s | Type: %s | Vertex: %i | Figures: pRF parameters and density'%(subject,roi_text,hemi,mask_dir,vertex))
                        save(all_f1)

                        all_f2 = gridplot([ [f_pRFmap,f_pRFcov[0]]],toolbar_location='right')
                        exec('output_file_html = opj(fig_bokeh_dir_{mask_dir}_{hemi},"{roi_text}_{hemi}_{mask_dir}_pRFmap.html")'.format(mask_dir = mask_dir,roi_text = roi_text, hemi = hemi))
                        output_file(output_file_html, title='Subject: %s | ROI: %s | Hemisphere: %s | Type: %s | Vertex: %i | Figures:pRF maps parameters and density'%(subject,roi_text,hemi,mask_dir,vertex))
                        save(all_f2)

                        all_f3 = gridplot([[f_pRFlat[0]]],toolbar_location = 'right')
                        exec('output_file_html = opj(fig_bokeh_dir_{mask_dir}_{hemi},"{roi_text}_{hemi}_{mask_dir}_pRFlat.html")'.format(mask_dir = mask_dir,roi_text = roi_text, hemi = hemi))
                        output_file(output_file_html, title='Subject: %s | ROI: %s | Hemisphere: %s | Type: %s | Vertex: %i | Figures:pRF laterality histogram'%(subject,roi_text,hemi,mask_dir,vertex))
                        save(all_f3)

                    else:
                        print("drawing {roi}_{hemi}_{mask_dir} figures not possible: n={vertex}".format(roi = roi_text,hemi = hemi,vertex = vertex,mask_dir = mask_dir)) 
        
