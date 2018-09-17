# check on figure details to insert with the other plot
# insert it in the main code with weighted values
import warnings
warnings.filterwarnings("ignore")

import h5py
import os
import numpy as np
import ipdb
opj = os.path.join

from bokeh.io import output_notebook, show,save, output_file, export_png, export_svgs
output_notebook()
from bokeh.plotting import figure, show, output_file
from bokeh.models import ColumnDataSource
from bokeh.models.tools import HoverTool

import cortex
import matplotlib.colors as colors

mask_dir = 'pos'
roi_text = 'SUP_PAR'

data_hemi = []
val_hemi = 0
rsq_idx, polar_real_idx, polar_imag_idx, x_idx, hemi_idx = 1, 3, 4, 10, 12
for hemi in ['L', 'R', 'LR']:

    if hemi == 'LR':
        data = np.column_stack((data_hemi[0],data_hemi[1]))
    else:
        if hemi == 'L': val_hemi = 1
        elif hemi == 'R': val_hemi = 2
        
        h5_dir = '/Users/martin/disks/ae_S/2018/visual/nprf_hcp/pp_data/999999/gauss/h5'
        folder_alias = '{hemi}_{mask_dir}'.format(hemi = hemi,mask_dir = mask_dir)
        h5_file = h5py.File(opj(h5_dir,'{roi}.h5'.format(roi = roi_text)), "r")
        in_file = opj("prf_deriv_{hemi}_{mask_dir}".format(hemi = hemi, mask_dir = mask_dir))
        data = h5_file['{folder_alias}/{in_file}'.format(folder_alias=folder_alias,in_file=in_file)]
        data = np.vstack((data,val_hemi*np.ones((1,data.shape[1]))))
        data_hemi.append(data)
        
def convert_on_axis(val_in,min_val,max_val,min_axis,max_axis):
    range_val = max_val - min_val
    range_axis = max_axis - min_axis
    val_out = (val_in/range_axis)*range_val + min_val
    return val_out

# figure parameters
min_val, max_val = 1, 2                          # drawn maximum and minimum
min_axis, max_axis = 0, 0.20                     # axis minimum and maximum
axis_tick_num = 5                                # axis tick number
bin_num = 24                                     # annular histogram bin number
hemi_col_L,hemi_col_R = '#ff6a00','#009dff'      # colors of hemisphere data
bg_col = tuple([250,250,250])                     # colors of center of the plot
weighted_data = False

data = data[:,~np.isnan(data[rsq_idx,:])]
weighted_text = 'R2-weighted '
output_file_png = "{roi}_{mask_dir}_r2weighted.png".format(roi = roi_text,mask_dir = mask_dir)
output_file_html = "{roi}_{mask_dir}_r2weighted.html".format(roi = roi_text,mask_dir = mask_dir)
if weighted_data == False:
    weighted_text = ''
    data[rsq_idx,:] = np.ones((1,data.shape[1]))
    output_file_png = "{roi}_{mask_dir}.png".format(roi = roi_text,mask_dir = mask_dir)
    output_file_html = "{roi}_{mask_dir}.html".format(roi = roi_text,mask_dir = mask_dir)

main_fig = figure( plot_width = 1000,
            plot_height = 1000,
            # title = "{roi}-{hemi} {weighted_text}laterality".format(roi = roi_text,hemi = hemi,weighted_text=weighted_text),
            x_axis_type = None, 
            y_axis_type = None,
            x_range = (-max_val*1.1, max_val*1.1), 
            y_range = (-max_val*1.1, max_val*1.1),
            min_border = 0, 
            outline_line_color = "white",
            background_fill_color = "white",
            tools =   "")

# figure background
cmap = 'hsv'
cmap_steps = 16
col_offset = 1/14.0
base = cortex.utils.get_cmap(cmap)
val = np.fmod(np.linspace(0+col_offset, 1+col_offset,cmap_steps,endpoint=False),1.0)
colmap = colors.LinearSegmentedColormap.from_list('my_colmap',base(val),N = cmap_steps)
vmin, vmax = 0, 2*np.pi
data_col = np.linspace(vmin,vmax,cmap_steps)
vrange = float(vmax) - float(vmin)
norm_data = ((data_col-float(vmin))/vrange)*cmap_steps
col_mat_rgb = colmap(norm_data.astype(int)) * 255.0
colors_val_rgb = ["#%02x%02x%02x" % (int(r), int(g), int(b)) for r, g, b in zip(col_mat_rgb[:,0], col_mat_rgb[:,1], col_mat_rgb[:,2])]

wedge_ang = np.linspace(0,2*np.pi,cmap_steps+1) + 2*np.pi/cmap_steps/2 + np.pi
wedge_ang = wedge_ang[:-1]
wedge_ang_step = wedge_ang[1]-wedge_ang[0]

main_fig.annular_wedge(x = 0,
                y = 0,
                inner_radius = 1, 
                outer_radius = 2, 
                start_angle = wedge_ang - wedge_ang_step/2, 
                end_angle = wedge_ang + wedge_ang_step/2, 
                direction = 'anticlock',
                fill_color = colors_val_rgb,
                fill_alpha = 0.25,
                line_color = None)

# histogram ticks
ticks_axis = np.linspace(min_axis,max_axis,axis_tick_num)
ticks_val = convert_on_axis(ticks_axis,min_val,max_val,min_axis,max_axis)
main_fig.circle(x = 0,
         y = 0,
         radius = ticks_val,
         fill_color = None,
         line_color = 'black',
         line_width = 1,
         line_dash = 'dashed')

# minor axes
for angLine_minor in np.arange(0,2*np.pi,2*np.pi/cmap_steps):
    line_x0_min,line_y0_min = min_val * np.cos(angLine_minor), min_val * np.sin(angLine_minor)
    line_x1_min,line_y1_min = max_val * np.cos(angLine_minor), max_val * np.sin(angLine_minor)
    main_fig.segment(x0 = line_x0_min,
              y0 = line_y0_min,
              x1 = line_x1_min,
              y1 = line_y1_min,
              line_color = 'black',
              line_width = 1,
              line_dash = 'dashed')
            
    tick_val = 0.05
    line_x0_min,line_y0_min = max_val * np.cos(angLine_minor), max_val * np.sin(angLine_minor)
    line_x1_min,line_y1_min = (max_val+tick_val) * np.cos(angLine_minor), (max_val+tick_val) * np.sin(angLine_minor)
    main_fig.segment(x0 = line_x0_min,
              y0 = line_y0_min,
              x1 = line_x1_min,
              y1 = line_y1_min,
              line_color = 'black',
              line_width = 1)

# major axes
for angLine_major in np.arange(0,2*np.pi,np.pi/2):
    line_x0_maj,line_y0_maj = min_val * np.cos(angLine_major), min_val * np.sin(angLine_major)
    line_x1_maj,line_y1_maj = (max_val+tick_val) * np.cos(angLine_major), (max_val+tick_val) * np.sin(angLine_major)
    main_fig.segment(x0 = line_x0_maj,
              y0 = line_y0_maj,
              x1 = line_x1_maj,
              y1 = line_y1_maj,
              line_color = "black")

# angular histogram
bins = 36
bin_angle = 2*np.pi/bins


if hemi == 'L' or hemi == 'LR':
    data_L = data[:,data[hemi_idx,:] == 1]
    weights_val_L = data_L[rsq_idx,:]
    pol_comp_num_L = data_L[polar_real_idx,:] + 1j * data_L[polar_imag_idx,:]
    polar_ang_L = np.angle(pol_comp_num_L)
    
    hist_L, bin_edges_L = np.histogram( a = polar_ang_L,
                                        range = (-np.pi-bin_angle/2,np.pi-bin_angle/2),
                                        bins = bins,
                                        weights = weights_val_L)
    hist_perc_L = hist_L/np.nansum(hist_L)
    hist_percent_L = hist_perc_L*100
    
    hist_val_L = convert_on_axis(hist_perc_L,min_val,max_val,min_axis,max_axis)
    start_angle_hist_L = bin_edges_L[:-1]
    end_angle_hist_L = bin_edges_L[1:]

    start_angle_hist_deg_L = np.degrees(start_angle_hist_L)
    end_angle_hist_deg_L = np.degrees(end_angle_hist_L)
    
    hist_data_source_L = {  'hist_L': hist_L,
                            'hist_percent_L': hist_percent_L,
                            'hist_val_L': hist_val_L,
                            'start_angle_L': start_angle_hist_L,
                            'end_angle_L': end_angle_hist_L,
                            'start_angle_deg_L': start_angle_hist_deg_L,
                            'end_angle_deg_L': end_angle_hist_deg_L}
    hist_source_L = ColumnDataSource(data = hist_data_source_L)
    
    an_wedges_L = main_fig.annular_wedge(  x = 0,
                                    y = 0,
                                    inner_radius = min_val,
                                    outer_radius = 'hist_val_L',
                                    start_angle = 'start_angle_L',
                                    end_angle = 'end_angle_L',
                                    fill_color = hemi_col_L,
                                    source = hist_source_L,
                                    line_width = 1,
                                    direction = 'anticlock',
                                    line_color = 'black',
                                    fill_alpha = 0.6,
                                    hover_fill_color = 'black',
                                    hover_line_color = 'black',
                                    hover_fill_alpha = 0.5,
                                    hover_line_alpha = 0.5)
    
    hist_tooltips_L = [ ('LH vertex', 'n = @hist_L{0}'),
                        ('Prop.', '@hist_percent_L{0.0}%'),
                        ('Edges','(@start_angle_deg_L{0} deg,@end_angle_deg_L{0}) deg')]
    hist_hover_L = HoverTool( tooltips = hist_tooltips_L,
                        mode = 'mouse',
                        renderers = [an_wedges_L])
    main_fig.add_tools(hist_hover_L)


if hemi == 'R' or hemi == 'LR':
    data_R = data[:,data[hemi_idx,:] == 2]
    weights_val_R = data_R[rsq_idx,:]
    pol_comp_num_R = data_R[polar_real_idx,:] + 1j * data_R[polar_imag_idx,:]
    polar_ang_R = np.angle(pol_comp_num_R)
    polar_ang_R = polar_ang_R[~np.isnan(polar_ang_R)]
    
    hist_R, bin_edges_R = np.histogram( a = polar_ang_R,
                                        range = (-np.pi-bin_angle/2,np.pi-bin_angle/2),
                                        bins = bins,
                                        weights = weights_val_R)

    hist_perc_R = hist_R/np.nansum(hist_R)
    hist_percent_R = hist_perc_R*100

    hist_val_R = convert_on_axis(hist_perc_R,min_val,max_val,min_axis,max_axis)
    start_angle_hist_R = bin_edges_R[:-1]
    end_angle_hist_R = bin_edges_R[1:]

    start_angle_hist_deg_R = np.degrees(start_angle_hist_R)
    end_angle_hist_deg_R = np.degrees(end_angle_hist_R)

    hist_data_source_R = {  'hist_R': hist_R,
                            'hist_percent_R': hist_percent_R,
                            'hist_val_R': hist_val_R,
                            'start_angle_R': start_angle_hist_R,
                            'end_angle_R': end_angle_hist_R,
                            'start_angle_deg_R': start_angle_hist_deg_R,
                            'end_angle_deg_R': end_angle_hist_deg_R}
    hist_source_R = ColumnDataSource(data = hist_data_source_R)
    
    an_wedges_R = main_fig.annular_wedge(  x = 0,
                                    y = 0,
                                    inner_radius = min_val,
                                    outer_radius = 'hist_val_R',
                                    start_angle = 'start_angle_R',
                                    end_angle = 'end_angle_R',
                                    fill_color = hemi_col_R,
                                    source = hist_source_R,
                                    line_width = 1,
                                    direction = 'anticlock',
                                    line_color = 'black',
                                    fill_alpha = 0.6,
                                    hover_fill_color = 'black',
                                    hover_line_color = 'black',
                                    hover_fill_alpha = 0.5,
                                    hover_line_alpha = 0.5)
    hist_tooltips_R = [ ('RH vertex', 'n = @hist_R{0}'),
                        ('Prop.', '@hist_percent_R{0.0}%'),
                        ('Edges','(@start_angle_deg_R{0} deg,@end_angle_deg_R{0} deg)')]
    hist_hover_R = HoverTool( tooltips = hist_tooltips_R,
                        mode = 'mouse',
                        renderers = [an_wedges_R])

    main_fig.add_tools(hist_hover_R)
    


# major axis values
main_fig.text(x = 0.125,
       y = ticks_val , 
       text = np.round(ticks_axis*100), 
       text_font_size = "16pt", 
       text_align = "center", 
       text_baseline = "bottom")

# axis label
main_fig.text(x = -0.12, 
       y = (max_val-min_val)/2+min_val, 
       text = ['Prop. (%)'],
       angle = np.pi/2, 
       text_font_size = "20pt", 
       text_align = "center")

# central plot
main_fig.circle(   x = 0,
            y = 0,
            radius = min_val,
            line_width = 1,
            fill_color = bg_col,
            line_color = 'black')

main_fig.circle(   x = 0,
            y = 0,
            radius = max_val,
            fill_color = None,
            line_color = 'black',
            line_width = 1)

# central plot axis
tick_val = 0.05
bar_height = 1
bar_width = 0.3
bar_ctr = [0,0]

# y axis
main_fig.segment(bar_ctr[0]-bar_width,bar_ctr[1]-bar_height/2,bar_ctr[0]-bar_width,bar_ctr[1]+bar_height/2,line_color = 'black',line_width = 2)
main_fig.segment(bar_ctr[0]-bar_width,bar_ctr[1]-bar_height/2,bar_ctr[0]-bar_width-tick_val,bar_ctr[1]-bar_height/2,line_color = 'black',line_width = 2)
main_fig.segment(bar_ctr[0]-bar_width,bar_ctr[1],bar_ctr[0]-bar_width-tick_val,bar_ctr[1],line_color = 'black',line_width = 2)
main_fig.segment(bar_ctr[0]-bar_width,bar_ctr[1]+bar_height/2,bar_ctr[0]-bar_width-tick_val,bar_ctr[1]+bar_height/2,line_color = 'black',line_width = 2)
main_fig.text(bar_ctr[0]-bar_width-0.2,bar_ctr[1],['Contra-laterality'],angle = np.pi/2,text_font_size = "20pt", text_align = "center")
main_fig.text(bar_ctr[0]-bar_width-0.075,bar_ctr[1],['index (%)'],angle = np.pi/2,text_font_size = "20pt", text_align = "center")

# x axis
main_fig.segment(bar_ctr[0]-bar_width,bar_ctr[1]-bar_height/2,bar_ctr[0]+bar_width,bar_ctr[1]-bar_height/2,line_color = 'black',line_width = 2)
main_fig.segment(bar_ctr[0]-bar_width,bar_ctr[1]-bar_height/2,bar_ctr[0]-bar_width,bar_ctr[1]-bar_height/2-tick_val,line_color = 'black',line_width = 2)
main_fig.segment(bar_ctr[0],bar_ctr[1]-bar_height/2,bar_ctr[0],bar_ctr[1]-bar_height/2-tick_val,line_color = 'black',line_width = 2)
main_fig.segment(bar_ctr[0]+bar_width,bar_ctr[1]-bar_height/2,bar_ctr[0]+bar_width,bar_ctr[1]-bar_height/2-tick_val,line_color = 'black',line_width = 2)
main_fig.text(bar_ctr[0]-bar_width/2,bar_ctr[1]-bar_height/2-0.2,['RH'],text_font_size = "20pt",text_align = "center")
main_fig.text(bar_ctr[0]+bar_width/2,bar_ctr[1]-bar_height/2-0.2,['LH'],text_font_size = "20pt",text_align = "center")

# plots
if hemi == 'R' or hemi == 'LR':
    val_R = np.sum(data_R[rsq_idx,data_R[x_idx,:] < 0])/np.sum(data_R[rsq_idx,:])
    val_text_R = '%1.1f %%'%(val_R*100)
    main_fig.quad( left = bar_ctr[0]-bar_width, 
            right = bar_ctr[0], 
            top = bar_ctr[1]-bar_height/2+val_R*bar_height, 
            bottom = bar_ctr[1]-bar_height/2,
            fill_color = hemi_col_R, 
            line_width = 2,
            line_color = 'black',
            fill_alpha = 0.8)
    main_fig.text(x = bar_ctr[0]-bar_width/2,
           y = bar_ctr[1]-bar_height/2 +(val_R*bar_height*0.5),
           text = [val_text_R],
            angle = np.pi/2,
           text_font_size = "20pt",
           text_align = "center",
           text_baseline = "middle",
           text_color = 'black')

if hemi == 'L' or hemi == 'LR':
    val_L = np.sum(data_L[rsq_idx,data_L[x_idx,:] > 0])/np.sum(data_L[rsq_idx,:])
    val_text_L = '%1.1f %%'%(val_L*100)
    main_fig.quad( left = bar_ctr[0], 
            right = bar_ctr[0]+bar_width, 
            top = bar_ctr[1]-bar_height/2+val_L*bar_height, 
            bottom = bar_ctr[1]-bar_height/2,
            fill_color = hemi_col_L, 
            line_width = 2,
            line_color = 'black',
            fill_alpha = 0.8)
    main_fig.text(x = bar_ctr[0]+bar_width/2,
           y = bar_ctr[1]-bar_height/2 +(val_L*bar_height*0.5),
           text = [val_text_L],
           angle = np.pi/2,
           text_font_size = "20pt",
           text_align = "center",
           text_baseline = "middle",
           text_color = 'black')

output_file(output_file_html, title="{roi}-{hemi} {weighted_text}laterality".format(roi = roi_text,hemi = hemi,weighted_text=weighted_text))
export_png(main_fig, filename=output_file_png)

