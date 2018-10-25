"""
-----------------------------------------------------------------------------------------
polar_prog.py
-----------------------------------------------------------------------------------------
Goal of the script:
Draw polar progression for ANG
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject number
sys.argv[2]: fit model ('gauss','css')
sys.argv[3]: save svg (1 = YES, 0 = NO)
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
source activate i27
cd /home/szinte/projects/retino_HCP
python post_fit/polar_prog.py 999999 css 1
python post_fit/polar_prog.py 999999 gauss 1
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
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colors as colors
import ipdb
import platform
import wquantiles
opj = os.path.join
deb = ipdb.set_trace

# MRI imports
# -----------
import nibabel as nb
import cortex

# Bokeh import
# ------------
from bokeh.plotting import figure 
from bokeh.models import ColumnDataSource
from bokeh.models.tools import HoverTool
from bokeh.models.annotations import Span, Label
from bokeh.layouts import row, column, gridplot
from bokeh.models import BoxZoomTool, BoxSelectTool, Spacer, WheelZoomTool, PanTool, ResetTool
from bokeh.models.glyphs import Text
from bokeh.models.mappers import LinearColorMapper
from bokeh.io import output_notebook, show,save, output_file, export_png, export_svgs
from bokeh.models import Arrow, OpenHead, NormalHead, VeeHead
import matplotlib.image as mpimg

# Settings
# --------
subject = sys.argv[1]
fit_model = sys.argv[2]
save_svg = int(sys.argv[3])
rot_left = 10
rot_right = 20
nbins = 20
small_size_dot = 5
big_size_dot = small_size_dot*3


# Functions import
# ----------------
from utils import set_pycortex_config_file, draw_cortex_vertex, rotate_pts, get_colors, roi_coord_mask, rot_coord

# Check system
# ------------
sys.exit('Drawing Flatmaps only works with Python 2. Aborting.') if sys.version_info[0] > 2 else None

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
if save_svg == 1:
    svg_folder = opj(base_dir,"pp_data",subject,fit_model,"figs","svg","neg","LR","ANG_LR_neg_pRFpol")
    try: os.makedirs(svg_folder)
    except: pass


# Change cortex database folder
# -----------------------------
pycortex_folder     =   opj(base_dir,'pp_data','cortex')
set_pycortex_config_file(project_folder = pycortex_folder)

# Create derivatives flatmaps
# ---------------------------
roi = 'ANG'
cmap_polar = 'hsv'
col_offset = 1/14.0
sign_idx, rsq_idx, ecc_idx, polar_real_idx, polar_imag_idx , size_idx, \
            non_lin_idx, amp_idx, baseline_idx, cov_idx, x_idx, y_idx = 0,1,2,3,4,5,6,7,8,9,10,11
mask_dir  = 'neg'

# Get data and combine hemispheres
deriv_mat=[]
for hemi in ['L','R']:
    deriv_file = nb.load(opj(deriv_dir,mask_dir,"prf_deriv_{hemi}_{mask_dir}_fsaverage.func.gii".format(hemi = hemi, mask_dir = mask_dir)))
    deriv_mat.append(np.array([deriv_file.darrays[i].data for i in range(len(deriv_file.darrays))]))
deriv_mat = np.hstack(deriv_mat)
deriv_mat = deriv_mat.T
threshold_mask = np.logical_and(np.logical_and( deriv_mat[:,rsq_idx]>=analysis_info['rsq_threshold'],
                                                deriv_mat[:,cov_idx]>=analysis_info['cov_threshold']),
                                                deriv_mat[:,size_idx]>=analysis_info['size_threshold'])

# R-square
rsq_data = deriv_mat[:,rsq_idx]
rsq_data[~threshold_mask] = np.nan
alpha = rsq_data

# Polar angle
real_num = deriv_mat[:,polar_real_idx]
real_num[~threshold_mask] = np.nan
imag_num = deriv_mat[:,polar_imag_idx]
imag_num[~threshold_mask] = np.nan
pol_comp_num = real_num + 1j * imag_num
polar_ang = np.angle(pol_comp_num)
ang_norm = (polar_ang + np.pi) / (np.pi * 2.0)

cmap_steps = 255
param_polar = {'data': ang_norm.T, 'cmap': cmap_polar, 'alpha': alpha.T, 'vmin': 0, 'vmax': 1, 'cmap_steps': cmap_steps,\
                       'curv_brightness': 0.05, 'curv_contrast': 0.1, 'cbar': 'polar', 'col_offset': col_offset,\
                        'subject': 'fsaverage','add_roi': False, 'zoom_roi': 'ANG', 'zoom_hem': 'left', 'zoom_margin': 1.0}
vertex_polar = draw_cortex_vertex(**param_polar)
plt.savefig(opj(svg_folder,'ANG_left_polar.png'),bbox_inches = 'tight')
plt.savefig(opj(svg_folder,'ANG_left_polar.pdf'),bbox_inches = 'tight')

param_polar = {'data': ang_norm.T, 'cmap': cmap_polar, 'alpha': alpha.T, 'vmin': 0, 'vmax': 1, 'cmap_steps': cmap_steps,\
                       'curv_brightness': 0.05, 'curv_contrast': 0.1, 'cbar': 'polar', 'col_offset': col_offset,\
                        'subject': 'fsaverage','add_roi': False, 'zoom_roi': 'ANG', 'zoom_hem': 'right', 'zoom_margin': 1.0}
vertex_polar = draw_cortex_vertex(**param_polar)
plt.savefig(opj(svg_folder,'ANG_right_polar.png'),bbox_inches = 'tight')
plt.savefig(opj(svg_folder,'ANG_right_polar.pdf'),bbox_inches = 'tight')

# -----------------
# Zoom ROI in bokeh
# -----------------

# Settings
arrow_start = [0.5,0.2]
arrow_end = [0.5,0.8]
arrow_ctr = [0.5,0.5]
p_width = 400
p_height = 400
min_border_large = 10
stim_color = tuple([250,250,250])
bg_color = tuple([229,229,229])
valtick = (np.linspace(0,1,11))

x_range_zoom = (0,1)
y_range_zoom = (0,1)
x_label_zoom = 'ROI short axis (norm.)'
y_label_zoom = 'ROI long axis (norm.)'

# Load images
left_filename = opj(svg_folder,'ANG_left_polar.png')
img_left = mpimg.imread(left_filename)
img_left = img_left*255.0
img_left = np.flipud(img_left)
right_filename = opj(svg_folder,'ANG_right_polar.png')
img_right = mpimg.imread(right_filename)
img_right = img_right*255.0
img_right = np.flipud(img_right)

# Left ANG
# --------
# Define figure
zoom_left = figure( plot_width = p_width,
                    plot_height = p_height,
                    x_range = x_range_zoom,
                    y_range = y_range_zoom,
                    min_border_left = min_border_large,
                    min_border_right = min_border_large,
                    min_border_bottom = min_border_large,
                    min_border_top = min_border_large,
                    toolbar_location = None,
                    title = '')

zoom_left.xaxis.axis_label = x_label_zoom
zoom_left.yaxis.axis_label = y_label_zoom
zoom_left.grid.grid_line_color = None
zoom_left.axis.minor_tick_in = False
zoom_left.axis.minor_tick_out = False
zoom_left.axis.major_tick_in = False
zoom_left.outline_line_alpha = 0
xaxis_ticker = np.linspace(0.075,0.95,10)
zoom_left.xaxis.ticker = xaxis_ticker

dict_label = {}
for val_num, val in enumerate(xaxis_ticker):
    dict_label.update({val: '%1.0g'%valtick[val_num]})
zoom_left.xaxis.ticker = xaxis_ticker    
zoom_left.xaxis.major_label_overrides = dict_label
yaxis_ticker = np.linspace(0.025,0.975,11)
valtick = (np.linspace(0,1,11))
dict_label = {}
for val_num, val in enumerate(yaxis_ticker):
    dict_label.update({val: '%1.0g'%valtick[val_num]})
zoom_left.yaxis.ticker = yaxis_ticker    
zoom_left.yaxis.major_label_overrides = dict_label
zoom_left.background_fill_color = bg_color
zoom_left.axis.axis_label_standoff = 10
zoom_left.axis.axis_label_text_font_style = 'normal'

# Draw roi
zoom_left.image_rgba(image=[img_left.astype('uint8')], x=-0.05, y=-0.06, dw=1.1,dh=1.1)

# Draw arrow
arrow_start_left = rotate_pts(arrow_start,arrow_ctr,-rot_left)
arrow_end_left = rotate_pts(arrow_end,arrow_ctr,-rot_left)
zoom_left.add_layout(Arrow(end = NormalHead(size = 10,fill_color = 'black', line_color = 'black',line_width = 6), 
                            x_start = arrow_start_left[0], y_start = arrow_start_left[1],
                            x_end = arrow_end_left[0], y_end = arrow_end_left[1], line_color = 'black',line_width = 6))
zoom_left.add_layout(Arrow(end = NormalHead(size = 10,fill_color = 'white', line_color= 'white',line_width = 2), 
                            x_start = arrow_start_left[0], y_start = arrow_start_left[1],
                            x_end = arrow_end_left[0], y_end = arrow_end_left[1],line_color = 'white',line_width = 2))

# Draw title
zoom_left.text(x = 0.05, y = 0.95, text = ['ANG LH'],text_color = 'white',text_font_size = '10pt', text_align = 'left', text_baseline = 'middle')

# Right ANG
# ---------
# Define figure
zoom_right = figure(plot_width = p_width,
                    plot_height = p_height,
                    x_range = x_range_zoom,
                    y_range = y_range_zoom,
                    min_border_left = min_border_large,
                    min_border_right = min_border_large,
                    min_border_bottom = min_border_large,
                    min_border_top = min_border_large,
                    toolbar_location = None,
                    title = '')

zoom_right.xaxis.axis_label = x_label_zoom
zoom_right.yaxis.axis_label = y_label_zoom
zoom_right.grid.grid_line_color = None
zoom_right.axis.minor_tick_in = False
zoom_right.axis.minor_tick_out = False
zoom_right.axis.major_tick_in = False
zoom_right.outline_line_alpha = 0

xaxis_ticker = np.linspace(0.18,0.85,11)
zoom_right.xaxis.ticker = xaxis_ticker
dict_label = {}
for val_num, val in enumerate(xaxis_ticker):
    dict_label.update({val: '%1.0g'%valtick[val_num]})
zoom_right.xaxis.ticker = xaxis_ticker    
zoom_right.xaxis.major_label_overrides = dict_label
yaxis_ticker = np.linspace(0.025,0.975,11)
valtick = (np.linspace(0,1,11))
dict_label = {}
for val_num, val in enumerate(yaxis_ticker):
    dict_label.update({val: '%1.0g'%valtick[val_num]})
zoom_right.yaxis.ticker = yaxis_ticker    
zoom_right.yaxis.major_label_overrides = dict_label
zoom_right.background_fill_color = bg_color
zoom_right.axis.axis_label_standoff = 10
zoom_right.axis.axis_label_text_font_style = 'normal'

# Draw roi
zoom_right.image_rgba(image=[img_right.astype('uint8')], x=-0.05, y=-0.06, dw=1.1,dh=1.1)

# Draw arrow
arrow_start_right = rotate_pts(arrow_start,arrow_ctr,-rot_right)
arrow_end_right = rotate_pts(arrow_end,arrow_ctr,-rot_right)
zoom_right.add_layout(Arrow(end = NormalHead(size = 10,fill_color = 'black', line_color= 'black',line_width=6), 
                            x_start = arrow_start_right[0], y_start = arrow_start_right[1],
                            x_end = arrow_end_right[0], y_end = arrow_end_right[1],line_color = 'black',line_width=6))
zoom_right.add_layout(Arrow(end = NormalHead(size = 10,fill_color = 'white', line_color= 'white',line_width=2), 
                            x_start = arrow_start_right[0], y_start = arrow_start_right[1],
                            x_end = arrow_end_right[0], y_end = arrow_end_right[1],line_color = 'white',line_width=2))

# Draw title
zoom_right.text(x = 0.05, y = 0.95, text = ['ANG RH'],text_color = 'white',text_font_size = '10pt', text_align = 'left', text_baseline = 'middle')

# ROI rotation
# ------------
# define left coordinates
ang_coord_left, ang_mask_left = roi_coord_mask('ANG','left')
alpha_left = alpha[ang_mask_left]
ang_norm_left = ang_norm[ang_mask_left]
real_num_left = real_num[ang_mask_left]
imag_num_left = imag_num[ang_mask_left]
ang_norm_left = ang_norm_left[~np.isnan(alpha_left)]
ang_norm_col_left = get_colors(ang_norm_left, 'hsv', 255, 1/14.0, 0, 1)
real_num_left = real_num_left[~np.isnan(alpha_left)]
imag_num_left = imag_num_left[~np.isnan(alpha_left)]
ang_coord_left = ang_coord_left[~np.isnan(alpha_left)]
ang_coord_rot_left = rot_coord(ang_coord_left,rot_left)
alpha_left = alpha_left[~np.isnan(alpha_left)]

# define right coordinates
ang_coord_right, ang_mask_right = roi_coord_mask('ANG','right')
alpha_right = alpha[ang_mask_right]
ang_norm_right = ang_norm[ang_mask_right]
real_num_right = real_num[ang_mask_right]
imag_num_right = imag_num[ang_mask_right]
ang_norm_right = ang_norm_right[~np.isnan(alpha_right)]
ang_norm_col_right = get_colors(ang_norm_right, 'hsv', 255, 1/14.0, 0, 1)
real_num_right = real_num_right[~np.isnan(alpha_right)]
imag_num_right = imag_num_right[~np.isnan(alpha_right)]
ang_coord_right = ang_coord_right[~np.isnan(alpha_right)]
ang_coord_rot_right = rot_coord(ang_coord_right,rot_right)
alpha_right = alpha_right[~np.isnan(alpha_right)]

# Polar angle progression computation
# -----------------------------------

# ANG Left
# --------
# Rename values
roi_long_left = ang_coord_rot_left[:,1]
weights_left = alpha_left
real_left = real_num_left
imag_left = imag_num_left
comp_num_left = real_left+ 1j*imag_left

# Define value for scatter plot
scatter_roi_long_left_norm = (roi_long_left - np.min(roi_long_left))/(np.max(roi_long_left)-np.min(roi_long_left))
scatter_angle_left = np.angle(real_left + 1j * imag_left)
scatter_angle_norm_left = (scatter_angle_left + np.pi) / (np.pi * 2.0)
scatter_angle_norm_left_col = get_colors(scatter_angle_norm_left, 'hsv', 255, 1/14.0, 0, 1)
scatter_rsq_left_norm = (weights_left - np.min(weights_left))/(np.max(weights_left)-np.min(weights_left))
scatter_size_dots_left = np.round(small_size_dot+scatter_rsq_left_norm*big_size_dot)

# Define binned values
bin_left = np.zeros(nbins+1)
roi_long_axis = np.linspace(0,1,nbins+1)
for tbin,bin_val in enumerate(roi_long_axis):
    bin_left[tbin] = wquantiles.quantile(data = roi_long_left, weights = weights_left, quantile = bin_val)

bin_rsq_left, bin_roi_long_left, bin_angle_left, bin_angle_norm_left = np.zeros(nbins), np.zeros(nbins), np.zeros(nbins), np.zeros(nbins)
bin_angle_left_ebx, bin_angle_left_std = [],[]
for bin_left_num in np.arange(0,nbins,1):
    val2pick = np.logical_and(roi_long_left >= bin_left[bin_left_num],roi_long_left <= bin_left[bin_left_num+1])
    # Average
    bin_comp_num_left = np.average(a = comp_num_left[val2pick],weights = weights_left[val2pick])
    bin_angle_left[bin_left_num] = np.angle(bin_comp_num_left)
    bin_angle_norm_left[bin_left_num] = (bin_angle_left[bin_left_num] + np.pi) / (np.pi * 2.0)
    bin_rsq_left[bin_left_num] = np.average(a = weights_left[val2pick],weights = weights_left[val2pick])
    bin_roi_long_left[bin_left_num] = np.average(a = roi_long_left[val2pick], weights = weights_left[val2pick])
    
    # Error bar
    bin_angle_left_ebx.append([ bin_roi_long_left[bin_left_num], bin_roi_long_left[bin_left_num]])
    bin_angle_left_std.append([ bin_angle_left[bin_left_num]-np.nanmean(weights_left[val2pick])/np.abs(bin_comp_num_left),\
                                bin_angle_left[bin_left_num]+np.nanmean(weights_left[val2pick])/np.abs(bin_comp_num_left)])
    
    # bin_angle_left_std.append([ bin_angle_left[bin_left_num]-np.nanstd(np.angle(comp_num_left[val2pick])),\
    #                             bin_angle_left[bin_left_num]+np.nanstd(np.angle(comp_num_left[val2pick]))])

# Normalization 
bin_angle_norm_left_col = get_colors(bin_angle_norm_left, 'hsv', 255, 1/14.0, 0, 1)
bin_rsq_left_norm = (bin_rsq_left - np.min(bin_rsq_left))/(np.max(bin_rsq_left)-np.min(bin_rsq_left))
bin_size_dots_left = np.round(small_size_dot+bin_rsq_left_norm*big_size_dot)
bin_roi_long_left_norm = (bin_roi_long_left - np.min(roi_long_left))/(np.max(roi_long_left)-np.min(roi_long_left))
bin_angle_left_ebx_norm = np.ndarray.tolist((bin_angle_left_ebx - np.min(roi_long_left))/(np.max(roi_long_left)-np.min(roi_long_left)))

# ANG Right
# ---------
# Rename values
roi_long_right = ang_coord_rot_right[:,1]
weights_right = alpha_right
real_right = real_num_right
imag_right = imag_num_right

# Add rotation to have it centered on 0 for plot
comp_num_right = real_right + 1j*imag_right
comp_num_rot_right = comp_num_right * (np.cos(np.pi)+1j*np.sin(np.pi))
real_right = np.real(comp_num_rot_right)
imag_right = np.imag(comp_num_rot_right)
comp_num_right = real_right+ 1j*imag_right

# Define value for scatter plot
scatter_roi_long_right_norm = (roi_long_right - np.min(roi_long_right))/(np.max(roi_long_right)-np.min(roi_long_right))
scatter_angle_right = np.angle(real_right + 1j * imag_right)
scatter_angle_norm_right = (scatter_angle_right + np.pi) / (np.pi * 2.0)
scatter_angle_norm_right_col = get_colors(scatter_angle_norm_right, 'hsv', 255, 8/14.0, 0, 1)
scatter_rsq_right_norm = (weights_right - np.min(weights_right))/(np.max(weights_right)-np.min(weights_right))
scatter_size_dots_right = np.round(small_size_dot+scatter_rsq_right_norm*big_size_dot)

# Define binned values
bin_right = np.zeros(nbins+1)
roi_long_axis = np.linspace(0,1,nbins+1)
for tbin,bin_val in enumerate(roi_long_axis):
    bin_right[tbin] = wquantiles.quantile(data = roi_long_right, weights = weights_right, quantile = bin_val)

bin_rsq_right, bin_roi_long_right, bin_angle_right, bin_angle_norm_right = np.zeros(nbins), np.zeros(nbins), np.zeros(nbins), np.zeros(nbins)
bin_angle_right_ebx, bin_angle_right_std = [],[]
for bin_right_num in np.arange(0,nbins,1):
    val2pick = np.logical_and(roi_long_right >= bin_right[bin_right_num],roi_long_right <= bin_right[bin_right_num+1])
    
    # Average
    bin_comp_num_right = np.average(a = comp_num_right[val2pick],weights = weights_right[val2pick])
    bin_angle_right[bin_right_num] = np.angle(bin_comp_num_right)
    bin_angle_norm_right[bin_right_num] = (bin_angle_right[bin_right_num] + np.pi) / (np.pi * 2.0)
    bin_rsq_right[bin_right_num] = np.average(a = weights_right[val2pick],weights = weights_right[val2pick])
    bin_roi_long_right[bin_right_num] = np.average(a = roi_long_right[val2pick], weights = weights_right[val2pick])
    
    # Error bar
    bin_angle_right_ebx.append([ bin_roi_long_right[bin_right_num], bin_roi_long_right[bin_right_num]])

    bin_angle_right_std.append([bin_angle_right[bin_right_num]-np.nanmean(weights_right[val2pick])/np.abs(bin_comp_num_right),\
                                bin_angle_right[bin_right_num]+np.nanmean(weights_right[val2pick])/np.abs(bin_comp_num_right)])

    # bin_angle_right_std.append([ bin_angle_right[bin_right_num]-np.nanstd(np.angle(comp_num_right[val2pick])),\
    #                             bin_angle_right[bin_right_num]+np.nanstd(np.angle(comp_num_right[val2pick]))])

# Normalization
bin_angle_norm_right_col = get_colors(bin_angle_norm_right, 'hsv', 255, 8/14.0, 0, 1)
bin_rsq_right_norm = (bin_rsq_right - np.min(bin_rsq_right))/(np.max(bin_rsq_right)-np.min(bin_rsq_right))
bin_size_dots_right = np.round(small_size_dot+bin_rsq_right_norm*big_size_dot)
bin_roi_long_right_norm = (bin_roi_long_right - np.min(roi_long_right))/(np.max(roi_long_right)-np.min(roi_long_right))
bin_angle_right_ebx_norm = np.ndarray.tolist((bin_angle_right_ebx - np.min(roi_long_right))/(np.max(roi_long_right)-np.min(roi_long_right)))

# --------
# Drawings
# --------

# General settings
x_range_pp = (0,1)
x_axis_ticker = np.linspace(x_range_pp[0],x_range_pp[1],11)
x_label_pp = 'Main direction (norm.)'
y_range_pp = (-np.pi, np.pi)
y_axis_ticker = np.linspace(y_range_pp[0],y_range_pp[1],5)
y_label_pp = 'Polar angle relative to opposite hemifield (rad)'
bin_angle_left_eb = bin_angle_left_std
bin_angle_right_eb = bin_angle_right_std

# Define ANG polar progression figure
pp_left = figure(   plot_width = p_width, plot_height = p_height, x_range = x_range_pp, y_range = y_range_pp, min_border_left = min_border_large, min_border_right = min_border_large,
                    min_border_bottom = min_border_large, min_border_top = min_border_large, toolbar_location = None, title = '')

pp_left.xaxis.axis_label = x_label_pp
pp_left.yaxis.axis_label = y_label_pp
pp_left.grid.grid_line_color = None
pp_left.axis.minor_tick_in = False
pp_left.axis.minor_tick_out = False
pp_left.axis.major_tick_in = False
pp_left.outline_line_alpha = 0
pp_left.xaxis.ticker = x_axis_ticker
pp_left.yaxis.ticker = y_axis_ticker
pp_left.background_fill_color = bg_color
pp_left.axis.axis_label_standoff = 10
pp_left.axis.axis_label_text_font_style = 'normal'
dict_label = {}
txt_label = ['-pi','-pi/2','0','pi/2','pi']
for val_num, val in enumerate(y_axis_ticker):
    dict_label.update({val: txt_label[val_num]})
pp_left.yaxis.major_label_overrides = dict_label

# Plot contra-lateral hemifield
pp_left.quad(bottom = -np.pi/2, left = 0, right = 1, top = np.pi/2, color = stim_color)

# Plot values
pp_left.circle(x = scatter_roi_long_left_norm, y = scatter_angle_left, line_color = scatter_angle_norm_left_col,line_alpha = 0.5,fill_color = None,size = scatter_size_dots_left,line_width = 0.5)
pp_left.multi_line(bin_angle_left_ebx_norm,bin_angle_left_eb,line_color = 'black',line_width = 1)
pp_left.line(x = bin_roi_long_left_norm, y = bin_angle_left, color = 'black',line_width = 1)
pp_left.circle(x = bin_roi_long_left_norm, y = bin_angle_left, line_color = 'black',fill_color = bin_angle_norm_left_col,size = bin_size_dots_left,line_width = 1)

# Define ANG polar progression figure
pp_right = figure(  plot_width = p_width, plot_height = p_height, x_range = x_range_pp, y_range = y_range_pp, min_border_left = min_border_large, min_border_right = min_border_large,
                    min_border_bottom = min_border_large, min_border_top = min_border_large, toolbar_location = None, title = '')

pp_right.xaxis.axis_label = x_label_pp
pp_right.yaxis.axis_label = y_label_pp
pp_right.grid.grid_line_color = None
pp_right.axis.minor_tick_in = False
pp_right.axis.minor_tick_out = False
pp_right.axis.major_tick_in = False
pp_right.outline_line_alpha = 0
pp_right.xaxis.ticker = x_axis_ticker
pp_right.yaxis.ticker = y_axis_ticker
pp_right.background_fill_color = bg_color
pp_right.axis.axis_label_standoff = 10
pp_right.axis.axis_label_text_font_style = 'normal'
dict_label = {}
txt_label = ['-pi','-pi/2','0','pi/2','pi']
for val_num, val in enumerate(y_axis_ticker):
    dict_label.update({val: txt_label[val_num]})
pp_right.yaxis.major_label_overrides = dict_label

# Plot contra-lateral hemifield
pp_right.quad(bottom = -np.pi/2, left = 0, right = 1, top = np.pi/2, color = stim_color)

# Plot values
pp_right.circle(x = scatter_roi_long_right_norm, y = scatter_angle_right, line_color = scatter_angle_norm_right_col,line_alpha = 0.5,fill_color = None,size = scatter_size_dots_right,line_width = 0.5)
pp_right.multi_line(bin_angle_right_ebx_norm,bin_angle_right_eb,line_color = 'black',line_width = 1)
pp_right.line(x = bin_roi_long_right_norm, y = bin_angle_right, color = 'black',line_width = 1)
pp_right.circle(x = bin_roi_long_right_norm, y = bin_angle_right, line_color = 'black',fill_color = bin_angle_norm_right_col,size = bin_size_dots_right,line_width = 1)

f = column(row(zoom_left,pp_left),row(zoom_right,pp_right))

# create folder
fig_bokeh_dir_neg_LR = opj(base_dir,"pp_data",subject,fit_model,"figs","prf","neg","LR")
try: os.makedirs(fig_bokeh_dir_neg)
except: pass
output_file_html = opj(fig_bokeh_dir_neg_LR,"ANG_LR_neg_pRFpol.html")

output_file(output_file_html, title='Subject: %s | ROI: ANG | Hemisphere: LR | Sign: neg | Figures: pRF polar progression'%subject)
save(f)


# save in svg
# -----------
if save_svg == 1:

    fig_dict = {'pp_left': pp_left,'pp_right': pp_right}
    # save figure
    for fig in fig_dict:
        fig_dict[fig].output_backend = 'svg'
        output_file_svg = opj(svg_folder,"ANG_LR_neg_pRFpol_{fig}.svg".format(fig = fig))
        export_svgs(fig_dict[fig], filename = output_file_svg)
