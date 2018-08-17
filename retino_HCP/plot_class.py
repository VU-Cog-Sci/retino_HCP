from __future__ import division, print_function
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

# Other imports
# -------------
import numpy as np
import cortex
import matplotlib.colors as colors
import scipy as sy
import numpy as np
import ipdb
deb = ipdb.set_trace

from popeye.spinach import generate_og_receptive_fields
from IPython import embed as shell


class PlotOperator(object):
    """
    class docstring
    """
    def __init__(self, **params):
        self.discrete_cbar = False
        self.params = params
        
        for k,v in params.items():
            setattr(self, k, v)


    def get_stim_circle(self, main_fig):
        return main_fig.circle(x = 0, y = 0, radius = self.stim_radius, color = self.stim_color)


    def get_span(self, dimension, location = 0,color = 'black'):
        fig_span     =   Span(location            =   location,                                   # create a infinte line span                                           
                              dimension           =   dimension,                                  # define dimension
                              line_alpha          =   0.5,                                        # define alpha of the line
                              line_color          =   color,                                    # define color of the line
                              line_width          =   1,                                          # define width of the line
                              line_dash           =   'dashed')                                   # define dash of the line
        return fig_span
        

    def save_figures_svg(self, main_fig, h_hist = None, v_hist = None, leg = None):
        from bokeh.io import export_svgs
        import os
        
        # Create folders to save the svg files in
        basedir                         =   os.path.join(self.svg_folder, self.roi_t, self.condition)
        hemi_dir                        =   os.path.join(basedir, self.hemisphere) if self.hemisphere != 'lrh' else basedir
        try:   os.makedirs(hemi_dir)
        except OSError: pass
        
        # save figures in dict to loop over them
        figures         = dict( main_fig =  main_fig,
                                h_hist   =  h_hist,
                                v_hist   =  v_hist,
                                leg      =  leg)

        # save figure
        for fig in ['main_fig', 'h_hist', 'v_hist', 'leg']:
            if figures[fig] != None:
                
                figures[fig].output_backend = 'svg'
                condition                   = self.condition if self.typeData == None else '{cond}_{typeData}'.\
                                                                                            format(cond = self.condition, typeData = self.typeData)
                fName                       = os.path.join(hemi_dir,'pRF_{condition}_{fig}.svg'.\
                                                                 format(condition = condition, fig = fig))
                export_svgs(figures[fig], filename = fName)

        return None


    def get_colors(self, data):
        import cortex
        import numpy as np
        import matplotlib.colors as colors

        base                =   cortex.utils.get_cmap(self.cmap)                       
        val                 =   np.fmod(np.linspace(0 + self.col_offset, 
                                                    1 + self.col_offset,  
                                                    self.cmap_steps + 1,
                                                    endpoint=False), 1.0)              
        self.colmap         =   colors.LinearSegmentedColormap.from_list(
                                                    'my_colmap', base(val), N = self.cmap_steps)           

        self.vrange         = float(self.vmax) - float(self.vmin)                           # define data range                                                 
        norm_data           = (( data - float(self.vmin) ) / self.vrange) * self.cmap_steps 
        col_mat_rgb         = self.colmap(norm_data.astype(int)) * 255.0
        colors_val_rgb      = ["#%02x%02x%02x" % (int(r), int(g), int(b)) for r, g, b in zip(col_mat_rgb[:,0], col_mat_rgb[:,1], col_mat_rgb[:,2])]
        return colors_val_rgb


    def get_weighted_regression_line(self, main_fig, data_source, non_nan = False, rsq_string = 'cv_rsq'):
        import numpy as np
        from scipy.optimize import curve_fit
        
        linear_function = lambda x, a, b: a * x + b
        
        x_reg                           =   data_source[self.x_source_label]                                                      # x data for regression line
        y_reg                           =   data_source[self.y_source_label]                                   # y data for regression line
        weight_reg                      =   data_source[rsq_string]                                                   # weight values for regression line

        if non_nan:
            x_reg                       =   x_reg[(~np.isnan(x_reg) & ~np.isnan(y_reg))]
            y_reg                       =   y_reg[(~np.isnan(x_reg) & ~np.isnan(y_reg))]

        coeffs, matcov                  =   curve_fit(                                                              # Use non-linear least squares to fit a function, f, to data.
                                                f                   =   linear_function,                                       # fit function to use
                                                xdata               =   x_reg,                                      # x data for regression fit
                                                ydata               =   y_reg,                                      # y data for regression fit
                                                sigma               =   weight_reg)                                 # weight
        # import ipdb ; ipdb.set_trace()
        x_fit                           =   np.arange(self.stim_fig_xlim[0], self.stim_fig_xlim[1] + self.x_tick_steps, self.x_tick_steps) # define fitted line x values
        y_fit                           =   linear_function(x_fit, coeffs[0], coeffs[1])                                       # define fitted line y values

        plot_reg                        =   main_fig.line(                                                          # draw fitted line on main figure
                                                x                   =  x_fit,                                       # define x of the line
                                                y                   =  y_fit,                                       # define y of the line
                                                line_color          =  'black',                                     # define line color
                                                line_width          =  4,
                                                line_alpha          =  0.5)                                         # define line alpha
        return plot_reg


    def create_horizontal_histogram(self, data_source, main_fig):

        v_a         = self.x_source_label
        h_hist_val, h_hist_edges = np.histogram(a = data_source[self.x_source_label], bins = self.h_hist_bins, range = self.x_range)

        h_hist_val_norm                 =   h_hist_val / np.shape(data_source[v_a])
        h_hist_val_norm_prct            =   h_hist_val_norm*100.0                                                   # put value in prctage


        data_h_hist                     =   np.abs((h_hist_edges[:-1] + h_hist_edges[1:])/2)
        colors_val_h_hist               =   self.get_colors(data_h_hist)
        
        # create hor histogram plot source
        h_hist_data_source              =   {   'hist_val':             h_hist_val,                                 # histogram value
                                                'hist_val_norm':        h_hist_val_norm,                            # histogram value normalized
                                                'hist_val_norm_prct':   h_hist_val_norm_prct,                       # histogram value normalized in percentage
                                                'bin_edges_left':       h_hist_edges[:-1],                          # histogram bin edges starting from left edge
                                                'bin_edges_right':      h_hist_edges[1:],                           # histogram bin edges strating from right edge
                                                'colors':               colors_val_h_hist                           # histogram colors
                                            }
        
        h_hist_source                   =   ColumnDataSource(data = h_hist_data_source)                             # define ColumnDataSource
        

        # define settings
        h_hist                          =   figure(                                                                 # create hor. histogram figure
                                                plot_width          =   self.p_width,                               # define figure widht
                                                plot_height         =   int(self.p_height/4),                       # define figure height
                                                x_range             =   main_fig.x_range,                           # define figure x range
                                                y_range             =   self.hist_range,                            # define figure y range
                                                min_border_bottom   =   self.min_border_small,                      # define bottom border space
                                                min_border_top      =   self.min_border_large,                      # define top border space
                                                y_axis_location     =   'left',                                     # define location of y axis
                                                x_axis_location     =   'below',                                    # define location of x axis
                                                title               =   self.main_fig_title,                        # specify title
                                                toolbar_location    =   None)                                       # specify toolbar location
        
        # determining axis ticker based on histogram orientation
        h_hist.yaxis.ticker = np.arange(self.hist_range[0],self.hist_range[1]+self.hist_steps,self.hist_steps)      # set y axis ticks
        h_hist.xaxis.ticker = np.arange(self.stim_fig_xlim[0], self.stim_fig_xlim[1], self.x_tick_steps)

        h_hist.grid.grid_line_color     =   None                                                                    # set axis grid color  
        h_hist.axis.minor_tick_in       =   False                                                                   # set axis minor tick in
        h_hist.axis.minor_tick_out      =   False                                                                   # set axis minor tick out
        h_hist.axis.major_tick_in       =   False                                                                   # set axis major tick in
        h_hist.outline_line_alpha       =   0                                                                       # set contour box alpha
        h_hist.background_fill_color    =   self.bg_color                                                      # set background color
        h_hist.xaxis.major_label_text_font_size = '0pt'                                                             # set x axis font size (make it disappear)

        # draw stim
        h_hist_stim                     =   h_hist.quad(                                                            # create a quad glyph of the stimulus
                                                bottom              =   self.stim_fig_ylim[0],                           # define bottom value
                                                left                =   -self.stim_radius,                                          # define left value
                                                right               =   +self.stim_radius,                     # define right value
                                                top                 =   self.stim_fig_ylim[1],                           # define top value
                                                color               =   self.stim_color)                       # define color

        # draw plot
        h_hist_plot                     =   h_hist.quad(                                                            # create a quad glyph of the histogram
                                                bottom              =   0,                                          # define bottom value
                                                left                =   'bin_edges_left',                           # define left value
                                                right               =   'bin_edges_right',                          # define right value
                                                top                 =   'hist_val_norm',                            # define top value
                                                source              =   h_hist_source,                              # define source
                                                fill_color          =   'colors',                                   # define color
                                                line_color          =   'black',                                    # define line color
                                                fill_alpha          =   0.5,
                                                line_alpha          =   0.5,
                                                hover_fill_color    =   'black',                                    # specify hover fill color 
                                                hover_line_color    =   'black',                                    # specify hover line color
                                                hover_fill_alpha    =   0.5,                                        # specify hover fill alpha
                                                hover_line_alpha    =   0.5)                                        # specify hover line alpha
                                                

        # add a span
        hor_center_hist                =   Span(                                                                    # define horizontal infinite span line
                                                location            =   0,                                          # define location value
                                                dimension           =   'height',                                   # define dimension value
                                                line_alpha          =   0.5,                                        # define alpha of the line
                                                line_color          =   'black',                                    # define color of the line
                                                line_width          =   1,                                          # define width of the line
                                                line_dash           =   'dashed')                                   # define dash of the line

        h_hist.add_layout(hor_center_hist)                                                                          # add infinite span line to figure
        
        
        # add a hover tool
        h_hist_tooltips                 = [     ('Voxels',              'n = @hist_val{0}'),                        # number of voxels
                                                ('Prop.',               '@hist_val_norm_prct{0.0} %'),              # proportion
                                                ('Edges',               '(@bin_edges_left{0.00},@bin_edges_right{0.00})')]# edge of distribution

        h_hist_hover                    =   HoverTool(                                                              # create hover tool
                                                tooltips            =   h_hist_tooltips,                            # specify content
                                                mode                =   'vline',                                    # specify mode
                                                renderers           =   [h_hist_plot])                              # specify renderer
        h_hist.add_tools(h_hist_hover)
        return h_hist


    def create_vertical_histogram(self, data_source, main_fig, colors = False, draw_stim = False):
        

        # compute histogram
        v_hist_val, v_hist_edges        =   np.histogram(a = data_source[self.y_source_label], bins = self.v_hist_bins, range = self.y_range)
        v_hist_val_norm                 =   v_hist_val / np.shape(data_source[self.y_source_label])
        v_hist_val_norm_prct            =   v_hist_val_norm*100.0                                                   # put value in prctage

        # create ver histogram plot source
        v_hist_data_source              =   {   'hist_val':             v_hist_val,                                 # histogram value
                                                'hist_val_norm':        v_hist_val_norm,                            # histogram value normalized
                                                'hist_val_norm_prct':   v_hist_val_norm_prct,                       # histogram value normalized in percentage
                                                'bin_edges_left':       v_hist_edges[:-1],                          # histogram bin edges starting from left edge
                                                'bin_edges_right':      v_hist_edges[1:],                           # histogram bin edges strating from right edge
                                            }

        if colors:
            data_hist                     =   np.abs((v_hist_edges[:-1] + v_hist_edges[1:])/2)
            colors_val_hist               =   self.get_colors(data = data_hist)
            v_hist_data_source.update(        dict(colors = colors_val_hist))
            self.hist_fill_color          =   'colors'

        v_hist_source                     =   ColumnDataSource(data = v_hist_data_source)                            # define ColumnDataSource
        

        # define settings
        v_hist                            =   figure(                                                                 # create vertical histogram figure
                                                plot_width          =   int(self.p_width/4),                   # define figure width 
                                                plot_height         =   self.p_height,                         # define figure height
                                                x_range             =   self.hist_range,                       # define figure x range
                                                y_range             =   main_fig.y_range,                           # define figure y range
                                                min_border_left     =   self.min_border_small,                 # define left border space
                                                min_border_right    =   self.min_border_large,                 # define right border space 
                                                y_axis_location     =   'left',                                     # define y axis location
                                                x_axis_location     =   'below',                                    # define x axis location
                                                toolbar_location    =   'right',
                                                tools               =   main_fig.tools)
        

        # # determining axis ticker based on histogram orientation
        v_hist.xaxis.ticker = np.arange(self.hist_range[0],self.hist_range[1]+self.hist_steps,self.hist_steps)      # set x axis ticks
        v_hist.yaxis.ticker = np.arange(self.stim_fig_ylim[0], self.stim_fig_ylim[1], self.y_tick_steps)            # set y axis ticks

        v_hist.grid.grid_line_color     =   None                                                                    # set both axis grid line color
        v_hist.axis.minor_tick_in       =   False                                                                   # set both axis minor tick in
        v_hist.axis.minor_tick_out      =   False                                                                   # set both axis minor tick out
        v_hist.axis.major_tick_in       =   False                                                                   # set both axis major tick in
        v_hist.yaxis.major_label_text_font_size = '0pt'                                                             # set y axis font size (make it disappear)
        v_hist.outline_line_alpha       =   0                                                                       # set box contour alpha
        v_hist.background_fill_color    =   self.bg_color                                                      # define background color

        if draw_stim:
            # draw stim
            v_hist_stim                     =   v_hist.quad(                                                            # create quad glyph of the stimulus
                                                left                =   self.stim_fig_ylim[0],                           # define left value
                                                bottom              =   +self.stim_radius,                     # define bottom value
                                                top                 =   -self.stim_radius,                     # define top value
                                                right               =   self.stim_fig_ylim[1],                           # define right value
                                                color               =   self.stim_color)                       # define color
        

        # draw plot
        v_hist_plot                     =   v_hist.quad(                                                            # create quad glyph of the vertical histogram
                                                left                =   0,                                          # define left value
                                                bottom              =   'bin_edges_left',                           # define bottom value
                                                top                 =   'bin_edges_right',                          # define top value
                                                right               =   'hist_val_norm',                            # define right value
                                                source              =   v_hist_source,                              # define source
                                                fill_color          =   self.hist_fill_color,                       # define fill color
                                                line_color          =   'black',                                    # define line color
                                                fill_alpha          =   0.5,                                        # define fill alpha
                                                line_alpha          =   0.5,                                        # define fill alpha
                                                hover_fill_color    =   'black',                                    # specify hover fill color 
                                                hover_line_color    =   'black',                                    # specify hover line color
                                                hover_fill_alpha    =   0.5,                                        # specify hover fill alpha
                                                hover_line_alpha    =   0.5)                                        # specify hover line alpha
        
        ver_center_hist                 =   Span(                                                                   # create vertical infinite span line
                                                location            =   0,                                          # define location
                                                dimension           =   'width',                                    # define dimension
                                                line_alpha          =   0.5,                                        # define line alpha
                                                line_color          =   'black',                                    # define line color
                                                line_width          =   1,                                          # define line width
                                                line_dash           =   'dashed')                                   # define line dash

        v_hist.add_layout(ver_center_hist)                                                                          # add infinite line span to figure
        
        # add a hover tool
        v_hist_tooltips                 = [     ('Voxels',              'n = @hist_val{0}'),                        # number of voxels
                                                ('Prop.',               '@hist_val_norm_prct{0.0} %'),              # proportion
                                                ('Edges',               '(@bin_edges_left{0.00},@bin_edges_right{0.00})')] # edge of distribution

        v_hist_hover                    =   HoverTool(                                                              # create hover tool
                                                tooltips            =   v_hist_tooltips,                            # specify content
                                                mode                =   'hline',                                    # specify mode
                                                renderers           =   [v_hist_plot])                              # specify renderer

        v_hist.add_tools(v_hist_hover)                                                                              # add hover tool to vhist plot
        
        return v_hist


    def create_legend_figure(self):
        from bokeh.models import ColumnDataSource
        from bokeh.plotting import figure
        
        # legend figure
        xy_max                          =   self.leg_xy_max_ratio*self.vmax
        leg                             =   figure(                                                                 # create vertical histogram figure
                                                plot_width          =   int(self.p_width/4),                   # define figure width 
                                                plot_height         =   int(self.p_height/4),                  # define figure height
                                                x_range             =   (-xy_max,xy_max),                           # define figure x range
                                                y_range             =   (-xy_max,xy_max),                           # define figure y range
                                                toolbar_location    =   None)                                       # define toolbar location
        leg.grid.grid_line_color        =   None                                                                    # set both axis grid line color
        leg.axis.major_label_text_font_size = '0pt'                                                                 # set y axis font size (make it disappear)
        leg.axis.axis_line_color        =   None                                                                    # set axis line color
        leg.axis.minor_tick_in          =   False                                                                   # set both axis minor tick in
        leg.axis.minor_tick_out         =   False                                                                   # set both axis minor tick out
        leg.axis.major_tick_in          =   False                                                                   # set both axis major tick in
        leg.axis.major_tick_out         =   False                                                                   # set both axis major tick in
        leg.outline_line_alpha          =   0                                                                       # set box contour alpha
        leg.background_fill_color       =   tuple([255,255,255])                                                    # define background color
        
        # define data leg
        data_leg                        =   np.linspace(self.vmax,self.vmin,self.cmap_steps+1)       # define colorbar values
        radius_val                      =   data_leg[:-1]
        data_leg                        =   (data_leg[:-1]+data_leg[1:])/2
        colors_val_leg                  =   self.get_colors(data = data_leg)

        # create data source
        leg_data_source                 =   {   'x':                    np.zeros(data_leg.shape),                   # coord x
                                                'y':                    np.zeros(data_leg.shape),                   # coord y
                                                'radius':               radius_val,                                 # radius
                                                'colors':               colors_val_leg                              # color
                                            }
        
        leg_source                      =   ColumnDataSource(data = leg_data_source)                                # define ColumnDataSource
        
        # draw circles
        plot_leg_circles                =   leg.circle(                                                             # create circle of the main plots
                                                x                   =   'x',                                        # define x coord
                                                y                   =   'y',                                        # define y coord
                                                radius              =   'radius',                                   # define radius
                                                fill_color          =   'colors',                                   # define fill color
                                                line_color          =   None,                                       # define line color
                                                fill_alpha          =   1,                                          # define fill alpha
                                                source              =   leg_source)                                 # specify data source
        plot_leg_line                   =   leg.line(
                                                x                   =  [self.vmax*1.2,self.vmax*1.2],
                                                y                   =  [0,self.vmax], 
                                                line_width          =  1.5,
                                                line_color          =  'black',
                                                line_cap            =  'round')
        
        text_leg                        =   Text(
                                                x                   =  self.vmax*1.2, 
                                                y                   =  self.vmax/2, 
                                                text                =  ['%0.f dva'%self.vmax], 
                                                text_color          = 'black',
                                                text_font_size      = '8pt',
                                                text_align          = 'center',
                                                text_baseline       = 'top',
                                                angle               = np.pi/2)
        leg.add_glyph(text_leg)
        return leg
        
    def initialize_main_fig_attributes(self, main_fig):

        main_fig.xaxis.axis_label       =   self.x_label                                                       # define x axis label
        main_fig.yaxis.axis_label       =   self.y_label                                                       # define y axis label
        main_fig.grid.grid_line_color   =   None                                                                    # define color of the grids for both axis
        main_fig.axis.minor_tick_in     =   False                                                                   # set minor tick in
        main_fig.axis.minor_tick_out    =   False                                                                   # set minor tick out
        main_fig.axis.major_tick_in     =   False                                                                   # set major tick in
        main_fig.outline_line_alpha     =   0                                                                       # change alpha of box contour
        main_fig.yaxis.ticker           =   np.arange(self.stim_fig_ylim[0], self.stim_fig_ylim[1], self.y_tick_steps)     # define y axis ticks
        main_fig.background_fill_color  =   self.bg_color                                                      # define backgroud color
        main_fig.axis.axis_label_standoff = 10                                                                      # set space between axis and label
        main_fig.axis.axis_label_text_font_style = 'normal' 
        if self.condition != 'roi': 
            main_fig.xaxis.ticker       = np.arange(self.stim_fig_xlim[0], self.stim_fig_xlim[1], self.x_tick_steps)     # define x axis ticks
        
        return main_fig


    def initialize_main_fig(self, old_main_fig =[], colors = True, gainRatio = None):
        # Bokeh import
        # ------------
        from bokeh.models import ColumnDataSource
        from bokeh.plotting import figure
        import numpy as np

        # Main figure
        # -----------
        if not old_main_fig:
            x_range                     =    self.x_range                                                   # define x range for the first time based on self.params
            y_range                     =    self.y_range                                                   # define y range for the first time based on self.params
        else:
            x_range                     =    old_main_fig.x_range                                           # define x range based on first figure to have shared axis
            y_range                     =    self.y_range
            # y_range                     =    old_main_fig.y_range                                           # define y range based on first figure to have shared axis
        
        data_source                     =   self.data_source
        if colors:
            color                       =   self.get_colors(data = data_source['colors_ref'])               # get the colors
        else:
            color                       =   []
        data_source.update({'color':color})

        # create the main plot source
        main_source                     =   ColumnDataSource(data = data_source)                            # define ColumnDataSource
        
        # define stimuli settings
        self.stim_fig_xlim              =   (self.x_range[0] - 5 * self.x_tick_steps, \
                                                        self.x_range[1] + 5 * self.x_tick_steps)            # define stimuli max axis
        self.stim_fig_ylim              =   (self.y_range[0] - 5 * self.y_tick_steps, \
                                                         self.y_range[1] + 5 * self.y_tick_steps)           # define stimuli max axis
        

        # figure settings
        main_fig                        =   figure(                                                         # create a figure in bokeh
                                                plot_width          =   self.p_width,                       # define figure width in pixel
                                                plot_height         =   self.p_height,                      # define figure height in pixel
                                                min_border_top      =   self.min_border_large,              # define top border size
                                                min_border_right    =   self.min_border_large,              # define right border size
                                                toolbar_location    =   None,                               # define toolbar location
                                                x_range             =   x_range,                            # define x limits
                                                y_range             =   y_range,                            # define y limits
                                                tools               =   "pan,wheel_zoom,box_zoom,reset")    # define tools   

        main_fig                        =   self.initialize_main_fig_attributes(main_fig = main_fig)
        
        main_fig.add_layout(self.get_span('width'))                                                         # add width span
        main_fig.add_layout(self.get_span('height'))                                                        # add height span
        
        return main_fig, main_source, data_source

    def draw_pRFmap(self, params, old_main_fig =[]):
        """
        -----------------------------------------------------------------------------------------
        draw_pRFmap(self, params,old_main_fig =[])
        -----------------------------------------------------------------------------------------
        Goal of the script:
        Create a graph with pRF position and size as well as side historgrams
        -----------------------------------------------------------------------------------------
        Input(s):
        params: dict containing a set of parameters for the figure
        old_main_fig: handle to the central figure to conserve same axis property across plots
        -----------------------------------------------------------------------------------------
        Output(s):
        none
        -----------------------------------------------------------------------------------------
        """
        
        self.condition                  =   'map'
        self.left_hist_lim              =   -self.stim_radius
        main_fig,main_source,data_source=   self.initialize_main_fig(old_main_fig)
        _                               =   main_fig.circle(x = 0, y = 0, radius = self.stim_radius, color = self.stim_color)

        # plot data
        plot_data                       =   main_fig.circle(                                                        # create circle of the main plots
                                                x                   =   'x',                                        # define x coord
                                                y                   =   'y',                                        # define y coord
                                                radius              =   'sigma',                                    # define radius
                                                fill_color          =   'color',                                    # define fill color
                                                line_color          =   'black',                                    # define fill color
                                                fill_alpha          =   0.5,                                        # define fill alpha
                                                line_alpha          =   0.5,                                        # define line alpha
                                                source              =   main_source,                                # specify data source
                                                hover_fill_color    =   'black',                                    # specify hover fill color 
                                                hover_line_color    =   'black',                                    # specify hover line color
                                                hover_fill_alpha    =   0.5,                                        # specify hover fill alpha
                                                hover_line_alpha    =   0.5)                                        # specify hover line alpha

        h_hist                          =   self.create_horizontal_histogram(data_source = data_source, main_fig = main_fig)
        v_hist                          =   self.create_vertical_histogram(data_source = data_source, main_fig = main_fig, colors = True, draw_stim = True)
        leg                             =   self.create_legend_figure()

        # add a hover tool
        main_fig_tooltips               =   [   ('R2',                  '@rsq{0.00}'),                              # rsq of hover tool
                                                ('[X,Y]',               '[@x{0.0}, @y{0.00}]'),                     # x,y coord of hover tool
                                                ('Eccentricity',        '@ecc{0.00} dva'),                          # eccentricity of hover tool
                                                ('Size',                '@sigma{0.00} dva'),                        # size of hover tool
                                                ('Baseline',            '@baseline{0.00}'),                         # baseline of hover tool
                                                ('Amplitude',           '@beta{0.00}'),                             # amplitude of hover tool
                                                ('Non-linarity',        '@non_lin{0.00}'),                          # non-linearity of hover tool
                                                ('Coverage',            '@cov{0.0}')]                               # coverage
        main_fig_hover                  =   HoverTool(                                                              # create hover tool
                                                tooltips            =   main_fig_tooltips,                          # specify content
                                                mode                =   'mouse',                                    # specify mode
                                                renderers           =   [plot_data])                                # specify renderer
        main_fig.add_tools(main_fig_hover)                                                                          # add hover tool to main plot

        # Put figure together
        # -------------------
        f                               =   column(                                                                 # define figures coluns
                                                row(h_hist,   leg),              # define figure first row
                                                row(main_fig, v_hist))                                              # define figure second row
        
        # save figures as .svg files
        if self.saving_figs: self.save_figures_svg(main_fig = main_fig, h_hist = h_hist, v_hist = v_hist, leg = leg)

        return (f,main_fig)

    def draw_pRFecc(self, params, old_main_fig =[]):
        """
        -----------------------------------------------------------------------------------------
        draw_pRFecc(self.params,old_main_fig =[])
        -----------------------------------------------------------------------------------------
        Goal of the script:
        Create a graph with pRF eccentricity as a function of ...
        -----------------------------------------------------------------------------------------
        Input(s):
        self.params: dict containing a set of parameters for the figure
        old_main_fig: handle to the central figure to conserve same axis property across plots
        -----------------------------------------------------------------------------------------
        Output(s):
        none
        -----------------------------------------------------------------------------------------
        """
        main_fig,main_source,data_source=   self.initialize_main_fig(old_main_fig)
        self.left_hist_lim              =   0
        
        # plot stimulus circle
        plot_stim                       =   main_fig.quad(                                                          # create a quad glyph of the stimulus
                                                bottom              =   self.stim_fig_ylim[0],                      # define bottom value
                                                left                =   0,                                          # define left value
                                                right               =   self.stim_radius,                           # define right value
                                                top                 =   self.stim_fig_ylim[1],                      # define top value
                                                color               =   self.stim_color)                            # define color
        # plot data
        
        plot_data                       =   main_fig.circle(                                                        # create circle of the main plots
                                                x                   =   'ecc',                                      # define x coord
                                                y                   =   self.y_source_label,                        # define y coord
                                                size                =   10,                                         # define radius
                                                fill_color          =   'color',                                    # define fill color
                                                line_color          =   'black',                                    # define fill color
                                                fill_alpha          =   0.5,                                        # define fill alpha
                                                line_alpha          =   0.5,                                        # define line alpha
                                                source              =   main_source,                                # specify data source
                                                hover_fill_color    =   'black',                                    # specify hover fill color 
                                                hover_line_color    =   'black',                                    # specify hover line color
                                                hover_fill_alpha    =   0.5,                                        # specify hover fill alpha
                                                hover_line_alpha    =   0.5)                                        # specify hover line alpha

        

        ### regression line weighted by cv rsq ###
        try:
            plot_reg                        =   self.get_weighted_regression_line(data_source = data_source, main_fig = main_fig)
        except: pass

        h_hist                          =   self.create_horizontal_histogram(data_source = data_source, main_fig = main_fig)
        v_hist                          =   self.create_vertical_histogram(data_source = data_source, main_fig = main_fig, colors = False, draw_stim = False)
        leg                             =   self.create_legend_figure()
        
        
        # add a hover tool
        main_fig_tooltips               =   [   ('R2',                  '@rsq{0.00}'),                              # rsq of hover tool
                                                ('[X,Y]',               '[@x{0.0}, @y{0.00}]'),                     # x,y coord of hover tool
                                                ('Eccentricity',        '@ecc{0.00} dva'),                          # eccentricity of hover tool
                                                ('Size',                '@sigma{0.00} dva'),                        # size of hover tool
                                                ('Baseline',            '@baseline{0.00}'),                         # baseline of hover tool
                                                ('Amplitude',           '@beta{0.00}'),                             # amplitude of hover tool
                                                ('Non-linarity',        '@non_lin{0.00}'),                          # non-linearity of hover tool
                                                ('Coverage',            '@cov{0.0}')]                               # coverage

        main_fig_hover                  =   HoverTool(                                                              # create hover tool
                                                tooltips            =   main_fig_tooltips,                          # specify content
                                                mode                =   'mouse',                                    # specify mode
                                                renderers           =   [plot_data])                                # specify renderer
        main_fig.add_tools(main_fig_hover)                                                                          # add hover tool to main plot

        # Put figure together
        # -------------------
        f                               =   column(                                                                 # define figures coluns
                                                row(h_hist,   leg),              # define figure first row
                                                row(main_fig, v_hist))                                              # define figure second row

        
        
        # save figures as .svg files
        # if self.saving_figs and self.condition == 'ecc': 
        #     self.save_figures_svg(main_fig = main_fig, h_hist = h_hist, v_hist = v_hist, leg = leg)
        # else:
        #     self.save_figures_svg(main_fig = main_fig, h_hist = h_hist, v_hist = v_hist, leg = None)

        
        return (f,main_fig)


    def draw_pRFroi(self, params, old_main_fig = [], condition = 'roi'):
        """
        -----------------------------------------------------------------------------------------
        draw_pRFroi(self.params)
        -----------------------------------------------------------------------------------------
        Goal of the script:
        Create a graph with parameters as a function of ROI
        -----------------------------------------------------------------------------------------
        Input(s):
        self.params: dict containing a set of parameters for the figure
        -----------------------------------------------------------------------------------------
        Output(s):
        none
        -----------------------------------------------------------------------------------------
        """
        # Additional imports
        # ------------------
        from bokeh.transform import dodge
        from bokeh.models import ColumnDataSource, Whisker
        self.condition                  = 'roi'

        # Main figure
        # -----------
        y_range                         =    self.y_range                                                      # define y range for the first time based on self.params
        x_range                         =    self.x_range                                                      # define x range for the first time based on self.params
        rois                            =    self.rois
        y_source_label                  =    self.y_source_label
        self.stim_fig_ylim              =   (self.y_range[0] - 5 * self.y_tick_steps, self.y_range[1] + 5 * self.y_tick_steps) # define stimuli max axis


        # define main data source
        data_source                     =   get_data_source_dict(data = self.dataMat, condition = self.condition)


        if y_source_label in ['retinal_x_gain', 'screen_x_gain', 'retinal_gain_index', 'amplitude_change', 'y_change']:
            upper       = [x+e for x,e in zip(data_source['%s'%self.y_source_label], data_source['%s_std'%self.y_source_label])]
            lower       = [x-e for x,e in zip(data_source['%s'%self.y_source_label], data_source['%s_std'%self.y_source_label])]
            data        = { 'rois' : rois,
                            y_source_label  : data_source[y_source_label],
                            '%s_std'%y_source_label : data_source['%s_std'%self.y_source_label],
                            'upper': upper,
                            'lower': lower}
        else:
            if self.y_source_label == 'voxel':
                upper_left  = []
                lower_left  = []
                upper_right = []
                lower_right = []
                data        = {'voxel_left_std'  : [],
                               'voxel_right_std' : []}
            else:
                upper_left  = [x+e for x,e in zip(data_source['%s_left'%self.y_source_label], data_source['%s_left_std'%self.y_source_label])]
                lower_left  = [x-e for x,e in zip(data_source['%s_left'%self.y_source_label], data_source['%s_left_std'%self.y_source_label])]
                upper_right = [x+e for x,e in zip(data_source['%s_right'%self.y_source_label], data_source['%s_right_std'%self.y_source_label])]
                lower_right = [x-e for x,e in zip(data_source['%s_right'%self.y_source_label], data_source['%s_right_std'%self.y_source_label])]
                data        = {'%s_left_std' : data_source['%s_left_std'%self.y_source_label],
                               '%s_right_std' : data_source['%s_right_std'%self.y_source_label]}

            data.update(      { 'rois' : rois,
                                '%s_left'%y_source_label  : data_source['%s_left'%y_source_label],
                                '%s_right'%y_source_label : data_source['%s_right'%y_source_label],
                                'upper_left' : upper_left,
                                'lower_left' : lower_left,
                                'upper_right': upper_right,
                                'lower_right': lower_right})
        

        # # define data source
        # data_source.update({'xs':          [[x1,x2] for x1, x2 in zip(data_source['x_left'], data_source['x_right'])], #dataMat_gazeLeft[:,0],dataMat_gazeRight[:,0])],   # list of list of xLeft/xRight
        #                     'ys':          [[y1,y2] for y1, y2 in zip(data_source['y_left'], data_source['y_right'])]})
        # # plot data
        # plot_data_line                  =   main_fig.multi_line(                                                    # create lines between gazeLeft/Right pRF
        #                                         xs                  =   self.xs,                               # list of list of x coord.
        #                                         ys                  =   self.ys,                               # list of list of y coord.
        #                                         line_color          =   'gray',                                    # define fill color
        #                                         line_alpha          =   0.2,                                        # define line alpha
        #                                         source              =   main_source,                                # define datasources
        #                                         hover_line_color    =   'black',                                    # specify hover line color
        #                                         hover_line_alpha    =   1                                           # define hover line alpha
        #                                         )

        # if self.y_source_label == 'gaze_gain':
        #     upper_gain = [x+e for x,e in zip(data_source[self.y_source_label], data_source['%s_std'%self.y_source_label])]
        #     lower_gain = [x-e for x,e in zip(data_source[self.y_source_label], data_source['%s_std'%self.y_source_label])]
        # elif self.y_source_label == 'voxel':
        #     upper_left = []
        #     lower_left = []
        #     upper_right = []
        #     lower_right = []
        # else:
        #     upper_left = [x+e for x,e in zip(data_source['%s_left'%self.y_source_label], data_source['%s_left_std'%self.y_source_label])]
        #     lower_left = [x-e for x,e in zip(data_source['%s_left'%self.y_source_label], data_source['%s_left_std'%self.y_source_label])]
        #     upper_right = [x+e for x,e in zip(data_source['%s_right'%self.y_source_label], data_source['%s_right_std'%self.y_source_label])]
        #     lower_right = [x-e for x,e in zip(data_source['%s_right'%self.y_source_label], data_source['%s_right_std'%self.y_source_label])]

            # data_source.update(dict(upper_left              =   upper_left,
            #                         lower_left              =   lower_left,
            #                         upper_right             =   upper_right,
            #                         lower_right             =   lower_right))


        # if self.y_source_label == 'retinal_x_gain':
        #     data_source.update({self.y_source_label : self.dataMat[-1][:,-5]})
        # elif self.y_source_label == 'screen_x_gain':
        #     data_source.update({self.y_source_label : self.dataMat[-1][:,-4]})
        # elif self.y_source_label == 'retinal_gain_index':
        #     data_source.update({self.y_source_label : self.dataMat[-1][:,-3]})
        # elif self.y_source_label == 'amplitude_change':
        #     data_source.update({self.y_source_label : self.dataMat[-1][:,-2]})
        # elif self.y_source_label == 'y_change':
        #     data_source.update({self.y_source_label : self.dataMat[-1][:,-1]})

       
        source      =  ColumnDataSource(data = data)
            
        # figure settings
        main_fig                        =   figure(                                                        # create a figure in bokeh
                                                plot_width          =   self.p_width,                          # define figure width in pixel
                                                plot_height         =   self.p_height,                         # define figure height in pixel
                                                min_border_left     =   self.min_border_large,                 # define left border size
                                                min_border_top      =   self.min_border_large,                 # define top border size
                                                min_border_right    =   self.min_border_large,                 # define right border size
                                                toolbar_location    =   None,                                  # define toolbar location
                                                x_range             =   rois,                          # define x limits
                                                y_range             =   y_range,                               # define y limits
                                                title               =   self.main_fig_title,                   # specify title
                                                tools               =   "")                                    # define tools
        main_fig                        =   self.initialize_main_fig_attributes(main_fig)
            

        if y_source_label in ['retinal_x_gain', 'screen_x_gain', 'retinal_gain_index', 'amplitude_change', 'y_change']:
            
            avg_data_all                    =   main_fig.vbar(
                                                    x                   =   'rois',                                 # define x value
                                                    top                 =   y_source_label,                       # define top value
                                                    bottom              =   0,                                      # define bottom value
                                                    width               =   self.bar_width_gain,                    # define bar widht
                                                    source              =   source,                                 # specify data source
                                                    fill_color          =   '#B774C1',                              # define fill color
                                                    line_color          =   'black',                                # define line color
                                                    fill_alpha          =   0.5,                                    # define fill alpha
                                                    line_alpha          =   1,                                      # define line alpha
                                                    hover_fill_color    =   'black',                                # specify hover fill color 
                                                    hover_line_color    =   'black',                                # specify hover line color
                                                    hover_fill_alpha    =   0.5,                                    # specify hover fill alpha
                                                    hover_line_alpha    =   1)                                      # specify hover line alpha
            # import ipdb ; ipdb.set_trace()
            renderers                       =   [avg_data_all]

            if self.across_subjects:
                main_fig_tooltips               =   [('%s'%self.tootips_txt, '@%s'%y_source_label+'{0.00}'+' %s'%self.tootips_end),
                                                     ('Std:', '@%s_std'%y_source_label+'{0.00}'+' %s'%self.tootips_end)]
            else:
                main_fig_tooltips               =   [('%s'%self.tootips_txt, '@%s'%y_source_label+'{0.00}'+' %s'%self.tootips_end)]
        
        else:

            avg_data_left                   =   main_fig.vbar(
                                                x                   =   dodge('rois',-0.125,range=main_fig.x_range), # define x value
                                                top                 =   '%s_left'%y_source_label,     # define top value
                                                width               =   self.bar_width,                    # define bar widht
                                                source              =   source,                            # specify data source
                                                fill_color          =   '#FF4136',                              # define fill color
                                                line_color          =   'black',                                # define line color
                                                fill_alpha          =   0.5,                                    # define fill alpha
                                                line_alpha          =   1,                                    # define line alpha
                                                legend              =   'gaze left',                           # define legend
                                                hover_fill_color    =   'black',                                # specify hover fill color 
                                                hover_line_color    =   'black',                                # specify hover line color
                                                hover_fill_alpha    =   0.5,                                    # specify hover fill alpha
                                                hover_line_alpha    =   1)                                    # specify hover line alpha
        
            avg_data_right                  =   main_fig.vbar(
                                                x                   =   dodge('rois', 0.125,range=main_fig.x_range), # define x value
                                                top                 =   '%s_right'%y_source_label,    # define top value
                                                width               =   self.bar_width,                    # define bar widht
                                                source              =   source,                            # specify data source
                                                fill_color          =   '#0074D9',                              # define fill color
                                                line_color          =   'black',                                # define line color
                                                fill_alpha          =   0.5,                                      # define fill alpha
                                                line_alpha          =   1,                                    # define line alpha
                                                legend              =   'gaze right',                          # define legend
                                                hover_fill_color    =   'black',                                # specify hover fill color 
                                                hover_line_color    =   'black',                                # specify hover line color
                                                hover_fill_alpha    =   0.5,                                    # specify hover fill alpha
                                                hover_line_alpha    =   1)                                    # specify hover line alpha
        
            renderers                       =   [avg_data_left,avg_data_right]

            if self.across_subjects:
                main_fig_tooltips           =   [ ('Left: %s'%self.tootips_txt, '@%s_left'%y_source_label+'{0.00}'+' %s'%self.tootips_end), 
                                                  ('Right: %s'%self.tootips_txt, '@%s_right'%y_source_label+'{0.00}'+' %s'%self.tootips_end),
                                                  ('Std Left:', '@%s_left_std'%y_source_label+'{0.00}'+' %s'%self.tootips_end),
                                                  ('Std Right:','@%s_right_std'%y_source_label+'{0.00}'+' %s'%self.tootips_end)]
            else:
                main_fig_tooltips           =   [ ('Left: %s'%self.tootips_txt, '@%s_left'%y_source_label+'{0.00}'+' %s'%self.tootips_end), 
                                                  ('Right: %s'%self.tootips_txt, '@%s_right'%y_source_label+'{0.00}'+' %s'%self.tootips_end)]
        # add some span
        hor_center_main_fig             =   Span(location            =   0,                                          # define location
                                                 dimension           =   'width',                                    # define dimension
                                                 line_alpha          =   0.5,                                        # define alpha of the line
                                                 line_color          =   'black',                                    # define color of the line
                                                 line_width          =   1)                                          # define width of the line
        main_fig.add_layout(hor_center_main_fig)                                                                    # add the horizontal line span to the figure
        main_fig.legend.border_line_alpha  = 0                                                                      # make legend box line transparent
        main_fig.legend.background_fill_alpha = 0                                                                   # make legend box background transparent

        if self.across_subjects:
            if y_source_label in ['retinal_x_gain', 'screen_x_gain', 'retinal_gain_index', 'amplitude_change', 'y_change']:
                main_fig.add_layout(Whisker(source=source, base="rois", upper="upper", lower="lower", level="overlay"))
            else:
                main_fig.add_layout(Whisker(source=source, base=dodge('rois',-0.125,range=main_fig.x_range), upper="upper_left", lower="lower_left", level="overlay"))
                main_fig.add_layout(Whisker(source=source, base=dodge('rois', 0.125,range=main_fig.x_range), upper="upper_right",lower="lower_right",level="overlay"))



        # add a hover tool
        main_fig_hover                  =   HoverTool(                                                              # create hover tool
                                                tooltips            =   main_fig_tooltips,                          # specify content
                                                mode                =   'mouse',                                    # specify mode
                                                renderers           =   renderers)                                  # specify renderer
        main_fig.add_tools(main_fig_hover)                                                                          # add hover tool to main plot


        # Put figure together
        # -------------------
        f                               =   column(row(main_fig))                                                   # define figure second row
        
        #import ipdb ; ipdb.set_trace()
        # save figures as .svg files
        if self.saving_figs: self.save_figures_svg(main_fig = main_fig)

        return (f,main_fig)


    def draw_pRFcov(self, params, old_main_fig =[]):
        """
        -----------------------------------------------------------------------------------------
        draw_pRFcov(params,old_main_fig =[])
        -----------------------------------------------------------------------------------------
        Goal of the script:
        Create a graph with pRF position and size as well as side historgrams
        -----------------------------------------------------------------------------------------
        Input(s):
        params: dict containing a set of parameters for the figure
        old_main_fig: handle to the central figure to conserve same axis property across plots
        -----------------------------------------------------------------------------------------
        Output(s):
        none
        -----------------------------------------------------------------------------------------
        """

        data_leg                        =   np.linspace(self.vmin,self.vmax,self.cmap_steps+1)       # define colorbar values
        data_leg_val                    =   (data_leg[:-1]+data_leg[1:])/2
        colors_val_leg                  =   self.get_colors(data_leg_val)
        dataMat                         =   self.dataMat
        smooth_factor                   =   self.smooth_factor
        
        deg_x, deg_y                    =   np.meshgrid(np.linspace(self.x_range[0]*1.5, self.x_range[1]*1.5, (-self.x_range[0]*1.5+self.x_range[1]*1.5)*smooth_factor), 
                                                        np.linspace(self.x_range[0]*1.5, self.x_range[1]*1.5, (-self.x_range[0]*1.5+self.x_range[1]*1.5)*smooth_factor))                     # define prfs in visual space


        pRFs                            =   generate_og_receptive_fields(                                           
                                                                          dataMat[:,10].astype(np.float64),         # coordinate of the center x of the Gaussian
                                                                          dataMat[:,11].astype(np.float64),         # coordinate of the center y of the Gaussian
                                                                          dataMat[:,5].astype(np.float64),          # dispersion of the Gaussian
                                                                          np.ones((dataMat.shape[0])),              # amplitude of the Gaussian
                                                                          deg_x,                                    # coordinate matrix along the horizontal dimension of the display (degrees)
                                                                          deg_y)                                    # coordinate matrix along the vertical dimension of the display (degrees)

        pRFs_rsq                        =   pRFs*dataMat[:,1]                                                       # weighted by cv rsq
        pRFs_rsq_sum                    =   np.nansum(pRFs_rsq,axis = 2)                                            # sum of pRF
        maxVal                          =   np.nanmax(pRFs_rsq_sum)                                                 # define max pixel
        minVal                          =   np.nanmin(pRFs_rsq_sum)
        valRange                        =   maxVal - minVal                                                         # define data range                                                 
        pRFs_rsq_sum_norm               =   ((pRFs_rsq_sum-minVal)/valRange)

        # define main data source
        data_source                     =   {'image':           pRFs_rsq_sum_norm,                               # pRF coverage image 
                                             'data_leg':        data_leg,                                           # scale for colorbar
                                             'colors':          colors_val_leg}
        # create the main plot source
        main_source                     =   ColumnDataSource(data = data_source)                                    # define ColumnDataSource
        
        # Main figure
        # -----------
        
        if not old_main_fig:
            x_range                         =    self.x_range                                                  # define x range for the first time based on params
            y_range                         =    self.y_range                                                  # define y range for the first time based on params
        else:
            x_range                         =    old_main_fig.x_range                                               # define x range based on first figure to have shared axis
            y_range                         =    old_main_fig.y_range                                               # define y range based on first figure to have shared axis
        


        # define stimuli settings
        stim_fig_xlim                   =   (self.x_range[0]-5*self.x_tick_steps,self.x_range[1]+5*self.x_tick_steps) # define stimuli max axis
        stim_fig_ylim                   =   (self.y_range[0]-5*self.y_tick_steps,self.y_range[1]+5*self.y_tick_steps) # define stimuli max axis
        
        # figure settings
        main_fig                        =   figure(                                                                 # create a figure in bokeh
                                                plot_width          =   self.p_width,                          # define figure width in pixel
                                                plot_height         =   self.p_height,                         # define figure height in pixel
                                                min_border_top      =   self.min_border_large,                 # define top border size
                                                min_border_right    =   self.min_border_large,                 # define right border size
                                                toolbar_location    =   None,                                       # define toolbar location
                                                x_range             =   x_range,                                    # define x limits
                                                y_range             =   y_range,                                    # define y limits
                                                tools               =   "pan,wheel_zoom,box_zoom,reset")            # define tools to show


        main_fig.xaxis.axis_label       =   self.x_label                                                       # define x axis label
        main_fig.yaxis.axis_label       =   self.y_label                                                       # define y axis label
        main_fig.grid.grid_line_color   =   None                                                                    # define color of the grids for both axis
        main_fig.axis.minor_tick_in     =   False                                                                   # set minor tick in
        main_fig.axis.minor_tick_out    =   False                                                                   # set minor tick out
        main_fig.axis.major_tick_in     =   False                                                                   # set major tick in
        main_fig.outline_line_alpha     =   0                                                                       # change alpha of box contour
        
        main_fig.yaxis.ticker           =   np.linspace(self.x_range[0], self.x_range[1], (-self.x_range[0]+self.x_range[1])/self.y_tick_steps+1)     # define y axis ticks
        main_fig.xaxis.ticker           =   np.linspace(self.y_range[0], self.y_range[1], (-self.y_range[0]+self.y_range[1])/self.y_tick_steps+1)     # define y axis ticks

        
        main_fig.background_fill_color  =   self.bg_color                                                           # define backgroud color
        main_fig.axis.axis_label_standoff = 10                                                                      # set space between axis and label
        main_fig.axis.axis_label_text_font_style = 'normal' 


        main_fig.add_layout(self.get_span('height',color = 'white'))                                                # add the vertical line span to the figure
        main_fig.add_layout(self.get_span('width',color = 'white'))                                                 # add the horizontal line span to the figure

        # colormap definition
        color_mapper                    =  LinearColorMapper(
                                                palette             =   colors_val_leg,                             # palette of color to use in Hex format
                                                low                 =   self.vmin,                                  # min of the colormap
                                                high                =   self.vmax)                                  # max of the colormap
                                                
        # plot data
        plot_data                       =  main_fig.image(
                                                image               =   [pRFs_rsq_sum_norm],                        # define image to plot
                                                x                   =   self.x_range[0]*1.5,                       # define x left bottom position
                                                y                   =   self.y_range[0]*1.5,                       # define y left bottom position
                                                dw                  =   [(-self.x_range[0]*1.5+self.x_range[1]*1.5)], # define width size
                                                dh                  =   [(-self.y_range[0]*1.5+self.y_range[1]*1.5)], # define height size
                                                color_mapper        =   color_mapper)                               # define colormap

        # plot stimulus circle
        plot_stim                       =   main_fig.circle(                                                        # create circle plots for the stim circle
                                                x                   =   0,                                          # define x coord
                                                y                   =   0,                                          # define y coord
                                                radius              =   self.stim_radius,                           # define radius
                                                line_alpha          =   0.5,                                        # define line alpha
                                                fill_color          =   None,                                       # define color
                                                line_color          =   'white',                                    # define line color
                                                line_dash           =   'dashed')

        # Colorbar
        # --------
        colorbar_fig                    =   figure(                                                                 # create vertical histogram figure
                                                plot_width          =   int(self.p_width/4),                   # define figure width 
                                                plot_height         =   self.p_height,                         # define figure height
                                                x_range             =  (0,1),                                       # define x range
                                                y_range             =  (self.vmin,self.vmax),             # define y range
                                                min_border_left     =   self.min_border_small,                 # define left border space
                                                min_border_right    =   self.min_border_large,                 # define right border space 
                                                toolbar_location    =   None,                                       # define toolbar location
                                                y_axis_location     =   'right',                                    # define y axis location
                                                )

        colorbar_fig.grid.grid_line_color     =   None                                                              # set both axis grid line color
        colorbar_fig.axis.minor_tick_in       =   False                                                             # set both axis minor tick in
        colorbar_fig.axis.minor_tick_out      =   False                                                             # set both axis minor tick out
        colorbar_fig.axis.major_tick_in       =   False                                                             # set both axis major tick in
        colorbar_fig.xaxis.major_tick_out     =   False                                                             # set both axis major tick in
        colorbar_fig.xaxis.major_label_text_font_size = '0pt'                                                       # set y axis font size (make it disappear)
        colorbar_fig.outline_line_alpha       =   0                                                                 # set box contour alpha
        colorbar_fig.axis.axis_label_standoff =   10                                                                # set space between axis and label
        colorbar_fig.yaxis.axis_label         =   self.cb_label
        colorbar_fig.axis.axis_label_text_font_style = 'normal'                                                     # set axis label font style
        colorbar_fig.yaxis.ticker             =   np.arange(self.vmin, self.vmax + self.cb_tick_steps,
                                                            self.cb_tick_steps) # define x axis ticks
        
        plot_colorbar                   =   colorbar_fig.quad(
                                                left                =   0,                                          # define left value
                                                bottom              =   data_leg[:-1],                              # define bottom value
                                                top                 =   data_leg[1:],                               # define top value
                                                right               =   1,                                          # define right value
                                                color               =   colors_val_leg)                             # define color
                                                
        
        # up space left
        s1 = Spacer(width=int(self.p_width), height=int(self.p_height/4))
        s2 = Spacer(width=int(self.p_width/4), height=int(self.p_height/4))

        # Put figure together
        # -------------------
        f                               =   column(                                                                 # define figures coluns
                                                row(s1, s2),                                                        # define figure second row
                                                row(main_fig, colorbar_fig))                                        # define figure second row

        # save figures as .svg files
        # if self.saving_figs: self.save_figures_svg(main_fig = main_fig, h_hist = None, v_hist = None, leg = colorbar_fig)

        return (f,main_fig)


    def draw_figure(self, parameters, plot, old_main_fig = []):
        
        # update params and atttributes to include plot specific params
        # current_params = self.params 
        # new_params     = current_params.update(parameters)
        # self.params    = new_params.update(parameters)
            
        for k,v in parameters.items():
            setattr(self, k, v)
        
        if plot == 'map':
            f, main_fig = self.draw_pRFmap(params = parameters, old_main_fig = old_main_fig)
        elif plot == 'ecc':
            f, main_fig = self.draw_pRFecc(params = parameters, old_main_fig = old_main_fig)
        elif plot == 'cov':
            f, main_fig = self.draw_pRFcov(params = parameters, old_main_fig = old_main_fig)
        
        return (f, main_fig)

