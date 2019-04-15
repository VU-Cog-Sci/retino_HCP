"""
-----------------------------------------------------------------------------------------
summary_plots.py
-----------------------------------------------------------------------------------------
Goal of the script:
Draw summary plots (ecc vs. size, lat per roi)
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject number
sys.argv[2]: fit model ('gauss','css')
sys.argv[3]: error-bar = 'std', '95ci' or 'sem'
sys.argv[4]: save_svg (1 = YES, 0 = NO)
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
source activate i27
cd /home/szinte/projects/retino_HCP
python post_fit/summary_plots.py sub-01 gauss_sg 95ci 0
python post_fit/summary_plots.py sub-02 gauss_sg 95ci 0
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
import ipdb
import platform
import h5py
from scipy.optimize import curve_fit
opj = os.path.join
deb = ipdb.set_trace

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

# Other imports
# -------------
import numpy as np
import cortex
import matplotlib.colors as colors
import scipy as sy
import numpy as np
import ipdb
deb = ipdb.set_trace

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

# Define figure parameters
# ------------------------
es_r_idx, es_slop_idx, es_intercept_idx, lat_idx, median_rsq_idx, mean_rsq_idx = 0, 1, 2, 3, 4, 5
subject = sys.argv[1]
fit_model = sys.argv[2]
subject = sys.argv[1]
fit_model = sys.argv[2]
eb = sys.argv[3]
save_svg = int(sys.argv[4])
avg = '_mean' # or '_median'

num_sum_val = 6
avg_subject = False

# Get data
# --------

for type_data in ["gridfit","fit"]:
    if avg_subject:
        end_file_avg = avg
    else:
        end_file_avg = ''
    h5file = h5py.File(opj(base_dir,'pp_data',subject,fit_model,'h5',type_data,'summary.h5'), "r")

    for type_fig in ['early_vis','late_vis','dmn']:
        exec('summary_mat_{type_fig} = np.zeros((len(analysis_info["{type_fig}_rois"]),num_sum_val))*np.nan'.format(type_fig=type_fig))
        exec('summary_eb_top_{type_fig} = np.zeros((len(analysis_info["{type_fig}_rois"]),num_sum_val))*np.nan'.format(type_fig=type_fig))
        exec('summary_eb_bot_{type_fig} = np.zeros((len(analysis_info["{type_fig}_rois"]),num_sum_val))*np.nan'.format(type_fig=type_fig))

        if type_fig == 'early_vis':
            mask_dir = 'pos'
        elif type_fig == 'late_vis':
            mask_dir = 'pos'
        elif type_fig == 'dmn':
            mask_dir = 'neg'

        for roi_num,roi in enumerate(analysis_info['{type_fig}_rois'.format(type_fig = type_fig)]):

            folder_alias = '{subject}/{roi}/{mask_dir}_summary{end_file}'.format(subject = subject, roi = roi, mask_dir = mask_dir, end_file = end_file_avg)

            exec('summary_mat_{type_fig}[roi_num,:] = h5file["{folder_alias}"]'.format(type_fig = type_fig,folder_alias=folder_alias))

            if avg_subject:
                folder_alias = '{subject}/{roi}/{mask_dir}_summary{end_file}'.format(subject = subject, roi = roi, mask_dir = mask_dir, end_file = '_{eb}'.format(eb=eb))
                if eb == '95ci':
                    exec('summary_eb_{type_fig} = h5file["{folder_alias}"]'.format(type_fig = type_fig,folder_alias=folder_alias))
                    exec('summary_eb_bot_{type_fig}[roi_num,:] = summary_eb_{type_fig}[0,:]'.format(type_fig = type_fig))
                    exec('summary_eb_top_{type_fig}[roi_num,:] = summary_eb_{type_fig}[1,:]'.format(type_fig = type_fig))
                else:
                    exec('summary_eb_{type_fig} = h5file["{folder_alias}"]'.format(type_fig = type_fig,folder_alias=folder_alias))
                    exec('summary_eb_bot_{type_fig}[roi_num,:] = summary_mat_{type_fig}[roi_num,:] - summary_eb_{type_fig}'.format(type_fig = type_fig))
                    exec('summary_eb_top_{type_fig}[roi_num,:] = summary_mat_{type_fig}[roi_num,:] + summary_eb_{type_fig}'.format(type_fig = type_fig))


    # Ecc/size figure
    # ---------------
    p_width = 400
    p_height = 400
    x_range_es = (0,10)
    y_range_es = (0,10)
    x_tick_steps_es = 2
    y_tick_steps_es = 2
    min_border_large = 10
    stim_width = analysis_info['stim_width']
    stim_color = tuple([250,250,250])
    bg_color = tuple([229,229,229])
    x_label_es = 'Eccentricity (dva)'
    y_label_es = 'Size (dva)'
    cmap = 'J4'

    # define colors
    cmap_steps = len(analysis_info['rois'])+1
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
    early_vis_num = len(analysis_info['early_vis_rois'])
    late_vis_num = len(analysis_info['late_vis_rois'])
    dmn_num = len(analysis_info['dmn_rois'])
    col_early_vis = colors_val_rgb[0:early_vis_num]
    col_late_vis = colors_val_rgb[early_vis_num:early_vis_num+late_vis_num]
    col_dmn = colors_val_rgb[early_vis_num+late_vis_num:early_vis_num+late_vis_num+dmn_num]

    # linear function
    linear_function = lambda x, a, b: a * x + b

    # define figure
    es_fig_all = []
    for type_fig in ['early_vis','late_vis','dmn']:
        if type_fig == 'early_vis':
            title = '{sub}: Low level vision ROIs: pRF sign [+]'.format(sub = subject)
        elif type_fig == 'late_vis':
            title = '{sub}: High level vision ROIs: pRF sign [+]'.format(sub = subject)
        elif type_fig == 'dmn':
            title = '{sub}: Default Network ROIs: pRF sign [-]'.format(sub = subject)
        exec('col = col_{type_fig}'.format(type_fig = type_fig))
        list_roi = analysis_info['{type_fig}_rois'.format(type_fig = type_fig)]

        es_fig = figure(plot_width = p_width,
                    plot_height = p_height,
                    x_range = x_range_es,
                    y_range = y_range_es,
                    min_border_left = min_border_large,
                    min_border_right = min_border_large,
                    min_border_bottom = min_border_large,
                    min_border_top = min_border_large,
                    toolbar_location = None,
                    title = title,
                    tools = "")
        es_fig.xaxis.axis_label = x_label_es
        es_fig.yaxis.axis_label = y_label_es
        es_fig.grid.grid_line_color = None
        es_fig.axis.minor_tick_in = False
        es_fig.axis.minor_tick_out = False
        es_fig.axis.major_tick_in = False
        es_fig.outline_line_alpha = 0
        es_fig.xaxis.ticker = np.arange(x_range_es[0],x_range_es[1] + x_tick_steps_es, x_tick_steps_es)
        es_fig.yaxis.ticker = np.arange(y_range_es[0],y_range_es[1] + y_tick_steps_es, y_tick_steps_es)
        es_fig.background_fill_color = bg_color
        es_fig.axis.axis_label_standoff = 10
        es_fig.axis.axis_label_text_font_style = 'normal'

        # plot stimulus area
        es_fig.quad(bottom = y_range_es[0], left = 0, right = stim_width, top = y_range_es[1], color = stim_color)

        # plot ecc/size lines
        x_es = np.array((x_range_es[0],x_range_es[1]))

        for roi_num,roi in enumerate(list_roi):

            exec('summary_mat = summary_mat_{type_fig}'.format(type_fig = type_fig))

            if avg_subject:

                exec('summary_eb_bot = summary_eb_bot_{type_fig}'.format(type_fig = type_fig))
                exec('summary_eb_top = summary_eb_top_{type_fig}'.format(type_fig = type_fig))


            # get values
            rsq_es = summary_mat[roi_num,es_r_idx]
            slope_es = summary_mat[roi_num,es_slop_idx]
            intercept_es = summary_mat[roi_num,es_intercept_idx]
            y_es = linear_function(x_es,slope_es,intercept_es)

            # draw error bar areas
            if avg_subject:
                # bottom
                rsq_es_bot = summary_eb_bot[roi_num,es_r_idx]
                slope_es_bot = summary_eb_bot[roi_num,es_slop_idx]
                intercept_es_bot = summary_eb_bot[roi_num,es_intercept_idx]
                y_es_bot = linear_function(x_es,slope_es_bot,intercept_es_bot)

                # top
                rsq_es_top = summary_eb_top[roi_num,es_r_idx]
                slope_es_top = summary_eb_top[roi_num,es_slop_idx]
                intercept_es_top = summary_eb_top[roi_num,es_intercept_idx]
                y_es_top = linear_function(x_es,slope_es_top,intercept_es_top)

                # draw error bar
                x_eb = [x_es[0],x_es[1],x_es[1],x_es[0]]
                y_eb = [y_es_bot[0],y_es_bot[1],y_es_top[1],y_es_top[0]]

                if eb == '95ci':
                    legend_txt = '{roi} (r: {:1.2g} +/- [{:1.2g},{:1.2g}])'.format(rsq_es,rsq_es_bot,rsq_es_top,roi = roi)
                else:
                    rsq_eb = rsq_es_top - rsq_es
                    legend_txt = '{roi} (r: {:1.2g} +/- {:1.2g})'.format(rsq_es, rsq_eb, roi = roi)
                es_fig.patch(x_eb,y_eb,line_color = None, fill_color = col[roi_num], fill_alpha = 0.1,legend=legend_txt)
            else:
                legend_txt = '{roi} (r: {:1.2g})'.format(rsq_es, roi = roi)

            # draw line
            es_fig.line(x_es,y_es,line_color=col[roi_num],legend=legend_txt,line_width = 2)

            # legend
            es_fig.legend.location = "top_left"
            es_fig.legend.click_policy="hide"
            es_fig.legend.background_fill_alpha = 0
            es_fig.legend.label_text_font = '8pt'
            es_fig.legend.margin = 10
            es_fig.legend.spacing = -2
            es_fig.legend.glyph_width = 10
            es_fig.legend.label_text_baseline = 'middle'
            es_fig.legend.border_line_color = None

        es_fig_all.append(es_fig)

    f1 = column( row(es_fig_all[0],es_fig_all[1],es_fig_all[2]))

    # Ecc/size r vs. ROI figures
    # ---------------------------
    es_roi_fig_all = []
    for type_fig in ['early_vis','late_vis','dmn']:
        p_width = 400
        p_height = 400
        x_range_es_roi = (-1,len(analysis_info["{type_fig}_rois".format(type_fig = type_fig)]))
        y_range_es_roi = (-0.5,1)
        x_tick_steps_es_roi = 1
        y_tick_steps_es_roi = 0.25
        x_label_es_roi = 'ROI'
        y_label_es_roi = 'Ecc. vs. Size r'

        exec('col = col_{type_fig}'.format(type_fig = type_fig))


        list_roi = analysis_info['{type_fig}_rois'.format(type_fig = type_fig)]

        es_roi_fig = figure(plot_width = p_width,
                        plot_height = int(p_height*0.4),
                        x_range = x_range_es_roi,
                        y_range = y_range_es_roi,
                        min_border_left = min_border_large,
                        min_border_right = min_border_large,
                        min_border_bottom = min_border_large,
                        min_border_top = min_border_large,
                        toolbar_location = None,
                        tools = "")
        es_roi_fig.xaxis.axis_label = ''
        es_roi_fig.yaxis.axis_label = y_label_es_roi
        es_roi_fig.grid.grid_line_color = None
        es_roi_fig.axis.minor_tick_in = False
        es_roi_fig.axis.minor_tick_out = False
        es_roi_fig.axis.major_tick_in = False
        es_roi_fig.outline_line_alpha = 0
        es_roi_fig.yaxis.ticker = np.arange(y_range_es_roi[0],y_range_es_roi[1] + y_tick_steps_es_roi, y_tick_steps_es_roi)
        es_roi_fig.xaxis.major_label_orientation = np.pi/3
        dict_label = {-1:''}
        for roi_num, roi in enumerate(list_roi):
            dict_label.update({roi_num: roi})
        dict_label.update({roi_num+1:''})
        es_roi_fig.xaxis.major_label_overrides = dict_label

        es_roi_fig.background_fill_color = bg_color
        es_roi_fig.yaxis.axis_label_standoff = 10
        es_roi_fig.xaxis.axis_label_standoff = 10
        es_roi_fig.axis.axis_label_text_font_style = 'normal'
        es_roi_fig.axis.axis_label_text_align = 'left'
        es_roi_fig.xaxis.major_label_text_font_size = '0pt'

        # span
        es_roi_fig.add_layout(Span(location = 0, dimension = 'width', line_alpha = 0.5, line_color = 'black', line_width = 1, line_dash = 'dashed'))

        # draw error bar
        if avg_subject:
            exec('summary_eb_bot = summary_eb_bot_{type_fig}'.format(type_fig = type_fig))
            exec('summary_eb_top = summary_eb_top_{type_fig}'.format(type_fig = type_fig))

            x_mat = []
            y_mat = []
            for roi_num,roi in enumerate(list_roi):
                x_mat.append([roi_num,roi_num])
                y_mat.append([summary_eb_bot[roi_num,es_r_idx],summary_eb_top[roi_num,es_r_idx]])

            es_roi_fig.multi_line(x_mat, y_mat,line_color = col,line_width = 2)

        # draw dots
        x_circle = np.zeros(len(list_roi))
        y_circle = np.zeros(len(list_roi))
        for roi_num,roi in enumerate(list_roi):
            exec('summary_mat = summary_mat_{type_fig}'.format(type_fig = type_fig))
            exec('col = col_{type_fig}'.format(type_fig = type_fig))
            x_circle[roi_num] = roi_num
            y_circle[roi_num] = summary_mat[roi_num,es_r_idx]

        es_roi_fig.circle(x = x_circle,y = y_circle, size = 10, fill_color = 'white', line_color = col,line_width = 2)


        es_roi_fig_all.append(es_roi_fig)

    f2 = column( row(es_roi_fig_all[0],es_roi_fig_all[1],es_roi_fig_all[2]))

    # Contra-laterality vs. ROI figures
    # ---------------------------------
    lat_roi_fig_all = []
    for type_fig in ['early_vis','late_vis','dmn']:
        p_width = 400
        p_height = 400
        x_range_lat_roi = (-1,len(analysis_info["{type_fig}_rois".format(type_fig = type_fig)]))
        y_range_lat_roi = (0,1)
        x_tick_steps_lat_roi = 1
        y_tick_steps_lat_roi = 0.25
        x_label_lat_roi = 'ROI'
        y_label_lat_roi = 'Contra-lat. index (%)'

        exec('col = col_{type_fig}'.format(type_fig = type_fig))


        list_roi = analysis_info['{type_fig}_rois'.format(type_fig = type_fig)]

        lat_roi_fig = figure(   plot_width = p_width,
                                plot_height = int(p_height*0.65),
                                x_range = x_range_lat_roi,
                                y_range = y_range_lat_roi,
                                min_border_left = min_border_large,
                                min_border_right = min_border_large,
                                min_border_bottom = min_border_large,
                                min_border_top = min_border_large,
                                toolbar_location = None,
                                tools = "")

        lat_roi_fig.xaxis.axis_label = x_label_lat_roi
        lat_roi_fig.yaxis.axis_label = y_label_lat_roi
        lat_roi_fig.grid.grid_line_color = None
        lat_roi_fig.axis.minor_tick_in = False
        lat_roi_fig.axis.minor_tick_out = False
        lat_roi_fig.axis.major_tick_in = False
        lat_roi_fig.outline_line_alpha = 0
        lat_roi_fig.yaxis.ticker = np.arange(y_range_lat_roi[0],y_range_lat_roi[1] + y_tick_steps_lat_roi, y_tick_steps_lat_roi)
        lat_roi_fig.xaxis.major_label_orientation = np.pi/3
        dict_label = {-1:''}
        for roi_num, roi in enumerate(list_roi):
            dict_label.update({roi_num: roi})
        dict_label.update({roi_num+1:''})
        lat_roi_fig.xaxis.major_label_overrides = dict_label

        lat_roi_fig.background_fill_color = bg_color
        lat_roi_fig.yaxis.axis_label_standoff = 10
        lat_roi_fig.xaxis.axis_label_standoff = 10
        lat_roi_fig.axis.axis_label_text_font_style = 'normal'
        lat_roi_fig.axis.axis_label_text_align = 'left'

        # draw error bar
        if avg_subject:
            exec('summary_eb_bot = summary_eb_bot_{type_fig}'.format(type_fig = type_fig))
            exec('summary_eb_top = summary_eb_top_{type_fig}'.format(type_fig = type_fig))

            x_mat = []
            y_mat = []
            for roi_num,roi in enumerate(list_roi):
                x_mat.append([roi_num,roi_num])
                y_mat.append([summary_eb_bot[roi_num,lat_idx],summary_eb_top[roi_num,lat_idx]])

            lat_roi_fig.multi_line(x_mat, y_mat,line_color = col,line_width = 2)

        # draw dots
        x_circle = np.zeros(len(list_roi))
        y_circle = np.zeros(len(list_roi))
        for roi_num,roi in enumerate(list_roi):
            exec('summary_mat = summary_mat_{type_fig}'.format(type_fig = type_fig))
            exec('col = col_{type_fig}'.format(type_fig = type_fig))
            x_circle[roi_num] = roi_num
            y_circle[roi_num] = summary_mat[roi_num,lat_idx]

        lat_roi_fig.circle(x = x_circle,y = y_circle, size = 10, fill_color = 'white', line_color = col,line_width = 2)


        lat_roi_fig_all.append(lat_roi_fig)

    f3 = column( row(lat_roi_fig_all[0],lat_roi_fig_all[1],lat_roi_fig_all[2]))

    # Combine plots together and save
    # -------------------------------
    all_f = gridplot([[f1],[f2],[f3]],toolbar_location = None)

    exec('fig_bokeh_dir = opj(base_dir,"pp_data","{subject}",fit_model,"figs","{type_data}","prf")'.format(subject = subject, type_data = type_data))
    try: os.makedirs(fig_bokeh_dir)
    except: pass
    if avg:
        output_file_html = opj(fig_bokeh_dir,"{subject}_summary_{eb}.html".format(subject = subject,eb = eb))
    else:
        output_file_html = opj(fig_bokeh_dir,"{subject}_summary.html".format(subject = subject,eb = eb))

    output_file(output_file_html, title="Subject: {subject} | Figures: Summary".format(subject = subject))
    print('saving {subject} summary figure'.format(subject=subject))
    save(all_f)


    # save in svg
    # -----------
    if save_svg == 1:
        from bokeh.io import export_svgs
        import os
        import time
        opj = os.path.join

        svg_folder = opj(base_dir,"pp_data",subject,fit_model,"figs",type_data,"svg","summary")
        try: os.makedirs(opj(svg_folder))
        except: pass

        fig_dict = {'es_fig_all_early':         es_fig_all[0], 
                    'es_fig_all_late':          es_fig_all[1], 
                    'es_fig_all_dmn':           es_fig_all[2],
                    'es_roi_fig_all_early':     es_roi_fig_all[0], 
                    'es_roi_fig_all_late':      es_roi_fig_all[1],
                    'es_roi_fig_all_dmn':       es_roi_fig_all[2],
                    'lat_roi_fig_all_early':    lat_roi_fig_all[0],
                    'lat_roi_fig_all_late':     lat_roi_fig_all[1],
                    'lat_roi_fig_all_dmn':      lat_roi_fig_all[2],
                    }

        # save figure
        for fig in fig_dict:
            fig_dict[fig].output_backend = 'svg'
            output_file_svg = opj(svg_folder,"{fig}.svg".format(fig = fig))

            export_svgs(fig_dict[fig], filename = output_file_svg)