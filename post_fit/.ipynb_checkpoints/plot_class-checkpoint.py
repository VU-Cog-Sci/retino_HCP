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
from decimal import Decimal
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


    def get_span(self, dimension, location = 0,color = 'black'):
        fig_span     =   Span(location            =   location,                                   # create a infinte line span                                           
                              dimension           =   dimension,                                  # define dimension
                              line_alpha          =   0.5,                                        # define alpha of the line
                              line_color          =   color,                                    # define color of the line
                              line_width          =   1,                                          # define width of the line
                              line_dash           =   'dashed')                                   # define dash of the line
        return fig_span
        

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


    def get_weighted_regression_line(self, main_fig, data_source, rsq_string = 'rsq'):
        import numpy as np
        from scipy.optimize import curve_fit
        
        
        linear_function = lambda x, a, b: a * x + b
        
        x_reg                           =   data_source[self.x_source_label]                                                      # x data for regression line
        y_reg                           =   data_source[self.y_source_label]                                   # y data for regression line
        weight_reg                      =   data_source[rsq_string]                                                   # weight values for regression line
        
        x_reg                       =   x_reg[(~np.isnan(x_reg) & ~np.isnan(y_reg))]
        y_reg                       =   y_reg[(~np.isnan(x_reg) & ~np.isnan(y_reg))]

        coeffs, matcov                  =   curve_fit(                                                              # Use non-linear least squares to fit a function, f, to data.
                                                f                   =   linear_function,                                       # fit function to use
                                                xdata               =   x_reg,                                      # x data for regression fit
                                                ydata               =   y_reg,                                      # y data for regression fit
                                                sigma               =   weight_reg)                                 # weight
        
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
                                                toolbar_location    =   None)                                       # specify toolbar location
        
        # determining axis ticker based on histogram orientation
        h_hist.yaxis.ticker = np.arange(self.hist_range[0],self.hist_range[1]+self.hist_steps,self.hist_steps)      # set y axis ticks
        h_hist.xaxis.ticker = np.arange(self.stim_fig_xlim[0], self.stim_fig_xlim[1], self.x_tick_steps)

        h_hist.grid.grid_line_color     =   None                                                                    # set axis grid color  
        h_hist.axis.minor_tick_in       =   False                                                                   # set axis minor tick in
        h_hist.axis.minor_tick_out      =   False                                                                   # set axis minor tick out
        h_hist.axis.major_tick_in       =   False                                                                   # set axis major tick in
        h_hist.outline_line_alpha       =   0                                                                       # set contour box alpha
        h_hist.background_fill_color    =   self.bg_color                                                           # set background color
        h_hist.xaxis.major_label_text_font_size = '0pt'                                                             # set x axis font size (make it disappear)

        # draw stim
        h_hist_stim                     =   h_hist.quad(                                                            # create a quad glyph of the stimulus
                                                bottom              =   self.stim_fig_ylim[0],                      # define bottom value
                                                left                =   -self.stim_width/2.0,                           # define left value
                                                right               =   +self.stim_width/2.0,                           # define right value
                                                top                 =   self.stim_fig_ylim[1],                      # define top value
                                                color               =   self.stim_color)                            # define color

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
        v_hist                            =   figure(                                                               # create vertical histogram figure
                                                plot_width          =   int(self.p_width/4),                        # define figure width 
                                                plot_height         =   self.p_height,                              # define figure height
                                                x_range             =   self.hist_range,                            # define figure x range
                                                y_range             =   main_fig.y_range,                           # define figure y range
                                                min_border_left     =   self.min_border_small,                      # define left border space
                                                min_border_right    =   self.min_border_large,                      # define right border space 
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
        v_hist.background_fill_color    =   self.bg_color                                                           # define background color

        if draw_stim:
            # draw stim
            v_hist_stim                     =   v_hist.quad(                                                        # create quad glyph of the stimulus
                                                left                =   self.stim_fig_ylim[0],                      # define left value
                                                bottom              =   +self.stim_height/2.0,                          # define bottom value
                                                top                 =   -self.stim_height/2.0,                          # define top value
                                                right               =   self.stim_fig_ylim[1],                      # define right value
                                                color               =   self.stim_color)                            # define color
        

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
            if self.link_x == True and self.link_y == False:
                x_range                     =    old_main_fig.x_range                                       # define x range based on first figure to have shared axis
                y_range                     =    self.y_range                                               # define y range based on last give value
            elif self.link_x == False and self.link_y == True:
                x_range                     =    self.x_range                                               # define x range based on last give value
                y_range                     =    old_main_fig.y_range                                       # define y range based on first figure to have shared axis
            elif self.link_y == True and self.link_y == True:
                x_range                     =    old_main_fig.x_range                                       # define x range based on first figure to have shared axis
                y_range                     =    old_main_fig.y_range                                       # define y range based on first figure to have shared axis
            
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
        main_fig,main_source,data_source=   self.initialize_main_fig(old_main_fig)
        
        _                               =  main_fig.quad(                                                           # create stimulus frame
                                                left                =   -self.stim_width/2.0,                       # define left value
                                                bottom              =   -self.stim_height/2.0,                      # define bottom value
                                                top                 =   +self.stim_height/2.0,                      # define top value
                                                right               =   +self.stim_width/2.0,                       # define right value
                                                color               =   self.stim_color)                            # define color
        
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

                                                
        # add title
        title_txt                       =   main_fig.text(
                                                x                   =   self.x_range[1] - (self.x_range[1]-self.x_range[0])*0.05,
                                                y                   =   self.y_range[0] + (self.y_range[1]-self.y_range[0])*0.05,
                                                text                =   [self.main_fig_title],
                                                text_font_size      =   "8pt",
                                                text_align          = "right"
                                                )

        # Put figure together
        # -------------------
        f                               =   column(                                                                 # define figures coluns
                                                row(h_hist,   leg),              # define figure first row
                                                row(main_fig, v_hist))                                              # define figure second row
        
        # save in svg
        # -----------
        if self.save_svg == 1:
            from bokeh.io import export_svgs
            import os
            opj = os.path.join
                  
            try:   os.makedirs(opj(self.svg_folder,self.svg_filename))
            except: pass
        
            fig_dict = {'h_hist': h_hist, 'leg': leg, 'main_fig':main_fig, 'v_hist':v_hist}
            # save figure
            for fig in fig_dict:
                fig_dict[fig].output_backend = 'svg'
                output_file_svg = opj(self.svg_folder,self.svg_filename,"{svg_filename}_{fig}.svg".format(svg_filename = self.svg_filename, fig = fig))
                export_svgs(fig_dict[fig], filename = output_file_svg)


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
        
        # plot stimulus area
        plot_stim                       =   main_fig.quad(                                                          # create a quad glyph of the stimulus
                                                bottom              =   self.stim_fig_ylim[0],                      # define bottom value
                                                left                =   0,                                          # define left value
                                                right               =   self.stim_width/2.0,                        # define right value
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

        

        # regression line weighted by r2
        if self.draw_reg:
            if self.dataMat.shape[0] > 10:
                plot_reg                        =   self.get_weighted_regression_line(data_source = data_source, main_fig = main_fig)
        

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

        # add title
        title_txt                       =   main_fig.text(
                                                x                   =   self.x_range[1] - (self.x_range[1]-self.x_range[0])*0.05,
                                                y                   =   self.y_range[0] + (self.y_range[1]-self.y_range[0])*0.05,
                                                text                =   [self.main_fig_title],
                                                text_font_size      =   "8pt",
                                                text_align          = "right"
                                                )

        # Put figure together
        # -------------------
        f                               =   column(                                                                 # define figures coluns
                                                row(h_hist,   leg),              # define figure first row
                                                row(main_fig, v_hist))                                              # define figure second row

        
        # save in svg
        # -----------
        if self.save_svg == 1:
            from bokeh.io import export_svgs
            import os
            opj = os.path.join
                  
            try:   os.makedirs(opj(self.svg_folder,self.svg_subfolder,self.svg_filename))
            except: pass
        
            fig_dict = {'h_hist': h_hist, 'leg': leg, 'main_fig':main_fig, 'v_hist':v_hist}
            # save figure
            for fig in fig_dict:
                fig_dict[fig].output_backend = 'svg'
                output_file_svg = opj(self.svg_folder,self.svg_subfolder,self.svg_filename,"{svg_filename}_{fig}.svg".format(svg_filename = self.svg_filename, fig = fig))
                export_svgs(fig_dict[fig], filename = output_file_svg)

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
            x_range                         =    self.x_range                                                       # define x range for the first time based on params
            y_range                         =    self.y_range                                                       # define y range for the first time based on params
        else:
            x_range                         =    old_main_fig.x_range                                               # define x range based on first figure to have shared axis
            y_range                         =    old_main_fig.y_range                                               # define y range based on first figure to have shared axis
        
        # define stimuli settings
        stim_fig_xlim                   =   (self.x_range[0]-5*self.x_tick_steps,self.x_range[1]+5*self.x_tick_steps) # define stimuli max axis
        stim_fig_ylim                   =   (self.y_range[0]-5*self.y_tick_steps,self.y_range[1]+5*self.y_tick_steps) # define stimuli max axis
        

        main_fig,main_source,data_source=   self.initialize_main_fig(old_main_fig)                                                

        main_fig.xaxis.axis_label       =   self.x_label                                                            # define x axis label
        main_fig.yaxis.axis_label       =   self.y_label                                                            # define y axis label
        main_fig.grid.grid_line_color   =   None                                                                    # define color of the grids for both axis
        main_fig.axis.minor_tick_in     =   False                                                                   # set minor tick in
        main_fig.axis.minor_tick_out    =   False                                                                   # set minor tick out
        main_fig.axis.major_tick_in     =   False                                                                   # set major tick in
        main_fig.outline_line_alpha     =   0                                                                       # change alpha of box contour
        
        main_fig.yaxis.ticker           =   np.linspace(self.x_range[0], self.x_range[1], (-self.x_range[0]+self.x_range[1])/self.y_tick_steps+1)     # define y axis ticks
        main_fig.xaxis.ticker           =   np.linspace(self.y_range[0], self.y_range[1], (-self.y_range[0]+self.y_range[1])/self.y_tick_steps+1)     # define y axis ticks

        
        main_fig.background_fill_color  =   colors_val_leg[0]                                                       # define backgroud color
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
        _                               =  main_fig.quad(                                                           # create stimulus frame
                                                left                =   -self.stim_width/2.0,                       # define left value
                                                bottom              =   -self.stim_height/2.0,                      # define bottom value
                                                top                 =   +self.stim_height/2.0,                      # define top value
                                                right               =   +self.stim_width/2.0,                       # define right value
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
                                                
        # add title
        title_txt                       =   main_fig.text(
                                                x                   =   self.x_range[1] - (self.x_range[1]-self.x_range[0])*0.05,
                                                y                   =   self.y_range[0] + (self.y_range[1]-self.y_range[0])*0.05,
                                                text                =   [self.main_fig_title],
                                                text_font_size      =   "8pt",
                                                text_align          =   "right",
                                                text_color          =   'white'
                                                )


        # up space left
        s1 = Spacer(width=int(self.p_width), height=int(self.p_height/4))
        s2 = Spacer(width=int(self.p_width/4), height=int(self.p_height/4))

        # Put figure together
        # -------------------
        f                               =   column(                                                                 # define figures coluns
                                                row(s1, s2),                                                        # define figure second row
                                                row(main_fig, colorbar_fig))                                        # define figure second row


        # save in svg
        # -----------
        if self.save_svg == 1:
            from bokeh.io import export_svgs
            import os
            opj = os.path.join
                  
            try:   os.makedirs(opj(self.svg_folder,self.svg_filename))
            except: pass
        
            fig_dict = {'main_fig': main_fig, 'colorbar_fig': colorbar_fig}
            # save figure
            for fig in fig_dict:
                fig_dict[fig].output_backend = 'svg'
                output_file_svg = opj(self.svg_folder,self.svg_filename,"{svg_filename}_{fig}.svg".format(svg_filename = self.svg_filename, fig = fig))
                export_svgs(fig_dict[fig], filename = output_file_svg)

        return (f,main_fig)

    

    def draw_pRFlat(self, params, old_main_fig =[]):
        """
        -----------------------------------------------------------------------------------------
        draw_pRFlat(params,old_main_fig =[])
        -----------------------------------------------------------------------------------------
        Goal of the script:
        Create a graph with pRF laterality index
        -----------------------------------------------------------------------------------------
        Input(s):
        params: dict containing a set of parameters for the figure
        old_main_fig: handle to the central figure to conserve same axis property across plots
        -----------------------------------------------------------------------------------------
        Output(s):
        none
        -----------------------------------------------------------------------------------------
        """
        
        def convert_on_axis(val_in,min_val,max_val,min_axis,max_axis):
            range_val = max_val - min_val
            range_axis = max_axis - min_axis
            val_out = (val_in/range_axis)*range_val + min_val
            return val_out

        # figure parameters
        min_val, max_val = 1, 2                                                                                     # drawn maximum and minimum
        min_axis, max_axis = self.vmin, self.vmax                                                                   # axis minimum and maximum
        axis_tick_num = 5                                                                                           # axis tick number
        bin_num = 24                                                                                                # annular histogram bin number
        hemi_col_L,hemi_col_R = '#ff6a00','#009dff'                                                                 # colors of hemisphere data
        bg_col = tuple([250,250,250])                                                                               # colors of center of the plot
        weighted_data = self.weighted_data
        

        # get data
        rsq_idx, polar_real_idx, polar_imag_idx, x_idx, hemi_idx = 1, 3, 4, 10, 12
        dataMat = self.dataMat
        data = dataMat[~np.isnan(dataMat[:,rsq_idx]),:]


        if weighted_data == False:
            data[:,rsq_idx] = np.ones((data.shape[0]))

        # figure settings
        main_fig                        =   figure(                                                                 # create a figure in bokeh
                                                plot_width          =   self.p_width,                               # define figure width in pixel
                                                plot_height         =   self.p_height,                              # define figure height in pixel
                                                x_range             =   self.x_range,                               # define x limits
                                                y_range             =   self.y_range,                               # define y limits
                                                min_border_left     =   self.min_border_large,                      # define left border size
                                                min_border_right    =   self.min_border_large,                      # define right border size
                                                min_border_bottom   =   self.min_border_large,                      # define bottom border space
                                                min_border_top      =   self.min_border_large,                      # define top border space
                                                x_axis_type         =   None,                                       # no main x axis
                                                y_axis_type         =   None,                                       # no main y axis
                                                outline_line_color  =   "white",                                    # ?
                                                tools               =   "pan,wheel_zoom,box_zoom,reset")            # define tools to show

        # figure background
        cmap =  self.cmap
        cmap_steps = self.cmap_steps
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

        bg_wedge                        =   main_fig.annular_wedge(
                                                x                   =   0,                                          # wedge center x
                                                y                   =   0,                                          # wedge center y
                                                inner_radius        =   min_val,                                    # wedge inner radius val
                                                outer_radius        =   max_val,                                    # wedge outer radius val
                                                start_angle         =   wedge_ang - wedge_ang_step/2,               # wedge starting angles
                                                end_angle           =   wedge_ang + wedge_ang_step/2,               # wedge ending angles
                                                direction           =   'anticlock',                                # wedge direction
                                                fill_color          =   colors_val_rgb,                             # wedges colors
                                                fill_alpha          =   0.25,                                       # wedge colors fill alpha
                                                line_color          =   None)                                       # wedge line colors
        
        
        # histogram ticks
        ticks_axis = np.linspace(min_axis,max_axis,axis_tick_num)
        ticks_val = convert_on_axis(ticks_axis,min_val,max_val,min_axis,max_axis)
        main_fig.circle(                        x                   =   0,
                                                y                   =   0,
                                                radius              =   ticks_val,
                                                fill_color          =   None,
                                                line_color          =   'black',
                                                line_width          =   0.5,
                                                line_dash           =   'dashed')

        # minor axes
        for angLine_minor in np.arange(0,2*np.pi,2*np.pi/cmap_steps):
            line_x0_min,line_y0_min = min_val * np.cos(angLine_minor), min_val * np.sin(angLine_minor)
            line_x1_min,line_y1_min = max_val * np.cos(angLine_minor), max_val * np.sin(angLine_minor)
            main_fig.segment(                       x0                  =   line_x0_min,
                                                    y0                  =   line_y0_min,
                                                    x1                  =   line_x1_min,
                                                    y1                  =   line_y1_min,
                                                    line_color          =   'black',
                                                    line_width          =   0.5,
                                                    line_dash           =   'dashed')
                    
            tick_val = 0.05
            line_x0_min,line_y0_min = max_val * np.cos(angLine_minor), max_val * np.sin(angLine_minor)
            line_x1_min,line_y1_min = (max_val+tick_val) * np.cos(angLine_minor), (max_val+tick_val) * np.sin(angLine_minor)
            main_fig.segment(                       x0                  =   line_x0_min,
                                                    y0                  =   line_y0_min,
                                                    x1                  =   line_x1_min,
                                                    y1                  =   line_y1_min,
                                                    line_color          =   'black',
                                                    line_width          =   0.5)

        # major axes
        for angLine_major in np.arange(0,2*np.pi,np.pi/2):
            line_x0_maj,line_y0_maj = min_val * np.cos(angLine_major), min_val * np.sin(angLine_major)
            line_x1_maj,line_y1_maj = (max_val+tick_val) * np.cos(angLine_major), (max_val+tick_val) * np.sin(angLine_major)
            main_fig.segment(                       x0                  =   line_x0_maj,
                                                    y0                  =   line_y0_maj,
                                                    x1                  =   line_x1_maj,
                                                    y1                  =   line_y1_maj,
                                                    line_color          =   "black")

        # angular histogram
        bins = self.ang_bins
        bin_angle = 2*np.pi/bins

        if self.hemi == 'L' or self.hemi == 'LR':
            data_L = data[data[:,hemi_idx]== 1,:]
            if data_L.shape[0] > 0:
                weights_val_L = data_L[:,rsq_idx]
                pol_comp_num_L = data_L[:,polar_real_idx] + 1j * data_L[:,polar_imag_idx]
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
                
                an_wedges_L = main_fig.annular_wedge(   x                   =   0,
                                                        y                   =   0,
                                                        inner_radius        =   min_val,
                                                        outer_radius        =   'hist_val_L',
                                                        start_angle         =   'start_angle_L',
                                                        end_angle           =   'end_angle_L',
                                                        fill_color          =   hemi_col_L,
                                                        source              =   hist_source_L,
                                                        line_width          =   0.5,
                                                        direction           =   'anticlock',
                                                        line_color          =   'black',
                                                        fill_alpha          =   0.6,
                                                        hover_fill_color    =   'black',
                                                        hover_line_color    =   'black',
                                                        hover_fill_alpha    =   0.5,
                                                        hover_line_alpha    =   0.5)
                
                hist_tooltips_L = [ ('LH vertex', 'n = @hist_L{0}'),
                                    ('Prop.', '@hist_percent_L{0.0}%'),
                                    ('Edges','(@start_angle_deg_L{0} deg,@end_angle_deg_L{0} deg)')]
                hist_hover_L = HoverTool(               tooltips            =   hist_tooltips_L,
                                                        mode                = 'mouse',
                                                        renderers           =   [an_wedges_L])
                main_fig.add_tools(hist_hover_L)

        if self.hemi == 'R' or self.hemi == 'LR':

            data_R = data[data[:,hemi_idx] == 2,:]

            if data_R.shape[0] > 0:
                weights_val_R = data_R[:,rsq_idx]
                pol_comp_num_R = data_R[:,polar_real_idx] + 1j * data_R[:,polar_imag_idx]
                polar_ang_R = np.angle(pol_comp_num_R)
                
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
                
                an_wedges_R = main_fig.annular_wedge(          x                   =   0,
                                                        y                   =   0,
                                                        inner_radius        =   min_val,
                                                        outer_radius        =   'hist_val_R',
                                                        start_angle         =   'start_angle_R',
                                                        end_angle           =   'end_angle_R',
                                                        fill_color          =   hemi_col_R,
                                                        source              =   hist_source_R,
                                                        line_width          =   0.5,
                                                        direction           =   'anticlock',
                                                        line_color          =   'black',
                                                        fill_alpha          =   0.6,
                                                        hover_fill_color    =   'black',
                                                        hover_line_color    =   'black',
                                                        hover_fill_alpha    =   0.5,
                                                        hover_line_alpha    =   0.5)
                
                hist_tooltips_R = [ ('RH vertex', 'n = @hist_R{0}'),
                                    ('Prop.', '@hist_percent_R{0.0}%'),
                                    ('Edges','(@start_angle_deg_R{0} deg,@end_angle_deg_R{0} deg)')]
                hist_hover_R = HoverTool( tooltips = hist_tooltips_R,
                                    mode = 'mouse',
                                    renderers = [an_wedges_R])

                main_fig.add_tools(hist_hover_R)
    
        # major axis values
        main_fig.text(x = 0.125,y = ticks_val , text = np.round(ticks_axis*100), text_font_size = "8pt", text_align = "center", text_baseline = "bottom")

        # axis label
        main_fig.text(x = -0.12, y = (max_val-min_val)/2+min_val, text = ['Prop. (%)'],angle = np.pi/2, text_font_size = "8pt", text_align = "center")

        # central plot
        main_fig.circle(   x = 0,y = 0,radius = min_val,line_width = 0.5,fill_color = bg_col,line_color = 'black')
        main_fig.circle(   x = 0,y = 0,radius = max_val,fill_color = None,line_color = 'black',line_width = 0.5)

        # central plot axis
        tick_val = 0.05
        bar_height = 1
        bar_width = 0.3
        bar_ctr = [0,0]

        # y axis
        main_fig.segment(bar_ctr[0]-bar_width,bar_ctr[1]-bar_height/2,bar_ctr[0]-bar_width,bar_ctr[1]+bar_height/2,line_color = 'black',line_width = 1)
        main_fig.segment(bar_ctr[0]-bar_width,bar_ctr[1]-bar_height/2,bar_ctr[0]-bar_width-tick_val,bar_ctr[1]-bar_height/2,line_color = 'black',line_width = 1)
        main_fig.segment(bar_ctr[0]-bar_width,bar_ctr[1],bar_ctr[0]-bar_width-tick_val,bar_ctr[1],line_color = 'black',line_width = 1)
        main_fig.segment(bar_ctr[0]-bar_width,bar_ctr[1]+bar_height/2,bar_ctr[0]-bar_width-tick_val,bar_ctr[1]+bar_height/2,line_color = 'black',line_width = 1)
        main_fig.text(bar_ctr[0]-bar_width-0.2,bar_ctr[1],['Contra-laterality'],angle = np.pi/2,text_font_size = "8pt", text_align = "center")
        main_fig.text(bar_ctr[0]-bar_width-0.075,bar_ctr[1],['index (%)'],angle = np.pi/2,text_font_size = "8pt", text_align = "center")

        # x axis
        main_fig.segment(bar_ctr[0]-bar_width,bar_ctr[1]-bar_height/2,bar_ctr[0]+bar_width,bar_ctr[1]-bar_height/2,line_color = 'black',line_width = 1)
        main_fig.segment(bar_ctr[0]-bar_width,bar_ctr[1]-bar_height/2,bar_ctr[0]-bar_width,bar_ctr[1]-bar_height/2-tick_val,line_color = 'black',line_width = 1)
        main_fig.segment(bar_ctr[0],bar_ctr[1]-bar_height/2,bar_ctr[0],bar_ctr[1]-bar_height/2-tick_val,line_color = 'black',line_width = 1)
        main_fig.segment(bar_ctr[0]+bar_width,bar_ctr[1]-bar_height/2,bar_ctr[0]+bar_width,bar_ctr[1]-bar_height/2-tick_val,line_color = 'black',line_width = 1)
        main_fig.text(bar_ctr[0]-bar_width/2,bar_ctr[1]-bar_height/2-0.2,['RH'],text_font_size = "8pt",text_align = "center")
        main_fig.text(bar_ctr[0]+bar_width/2,bar_ctr[1]-bar_height/2-0.2,['LH'],text_font_size = "8pt",text_align = "center")

        # plots
        if self.hemi == 'R' or self.hemi == 'LR':
            if data_R.shape[0] > 0:
                val_R = np.sum(data_R[data_R[:,x_idx] < 0,rsq_idx])/np.sum(data_R[:,rsq_idx])
                val_text_R = '%1.1f %%'%(val_R*100)
                main_fig.quad( left = bar_ctr[0]-bar_width, 
                        right = bar_ctr[0], 
                        top = bar_ctr[1]-bar_height/2+val_R*bar_height, 
                        bottom = bar_ctr[1]-bar_height/2,
                        fill_color = hemi_col_R, 
                        line_width = 1,
                        line_color = 'black',
                        fill_alpha = 0.8)
                main_fig.text(x = bar_ctr[0]-bar_width/2,
                       y = bar_ctr[1]-bar_height/2 +(val_R*bar_height*0.5),
                       text = [val_text_R],
                        angle = np.pi/2,
                       text_font_size = "8pt",
                       text_align = "center",
                       text_baseline = "middle",
                       text_color = 'black')

        if self.hemi == 'L' or self.hemi == 'LR':
            if data_L.shape[0] > 0:
                val_L = np.sum(data_L[data_L[:,x_idx] > 0,rsq_idx])/np.sum(data_L[:,rsq_idx])
                val_text_L = '%1.1f %%'%(val_L*100)
                main_fig.quad( left = bar_ctr[0], 
                        right = bar_ctr[0]+bar_width, 
                        top = bar_ctr[1]-bar_height/2+val_L*bar_height, 
                        bottom = bar_ctr[1]-bar_height/2,
                        fill_color = hemi_col_L, 
                        line_width = 1,
                        line_color = 'black',
                        fill_alpha = 0.8)
                main_fig.text(x = bar_ctr[0]+bar_width/2,
                       y = bar_ctr[1]-bar_height/2 +(val_L*bar_height*0.5),
                       text = [val_text_L],
                       angle = np.pi/2,
                       text_font_size = "8pt",
                       text_align = "center",
                       text_baseline = "middle",
                       text_color = 'black')


        # title
        title_txt                       =   main_fig.text(
                                                x                   =   self.x_range[1] - (self.x_range[1]-self.x_range[0])*0.05,
                                                y                   =   self.y_range[0] + (self.y_range[1]-self.y_range[0])*0.05,
                                                text                =   [self.main_fig_title],
                                                text_font_size      =   "8pt",
                                                text_align          =   "right",
                                                )
        # up space
        s1 = Spacer(width=int(self.p_width), height=int(self.p_height/4))




        # Put figure together
        # -------------------
        f                               =   column(                                                                 # define figures coluns
                                                row(main_fig))                                                  # define figure second row

        # save in svg
        # -----------
        if self.save_svg == 1:
            from bokeh.io import export_svgs
            import os
            opj = os.path.join
                  
            try:   os.makedirs(opj(self.svg_folder,self.svg_filename))
            except: pass
        
            fig_dict = {'main_fig': main_fig}
            # save figure
            for fig in fig_dict:
                fig_dict[fig].output_backend = 'svg'
                output_file_svg = opj(self.svg_folder,self.svg_filename,"{svg_filename}_{fig}.svg".format(svg_filename = self.svg_filename, fig = fig))
                export_svgs(fig_dict[fig], filename = output_file_svg)

        return (f,main_fig)

    

    def draw_pRFtc(self, params):
        """
        -----------------------------------------------------------------------------------------
        draw_pRFtc(self, params)
        -----------------------------------------------------------------------------------------
        Goal of the script:
        Create a graph with pRF timecourse
        -----------------------------------------------------------------------------------------
        Input(s):
        params: dict containing a set of parameters for the figure
        old_main_fig: handle to the central figure to conserve same axis property across plots
        -----------------------------------------------------------------------------------------
        Output(s):
        none
        -----------------------------------------------------------------------------------------
        """
        sign_idx, rsq_idx, ecc_idx, polar_real_idx, polar_imag_idx , size_idx, \
            non_lin_idx, amp_idx, baseline_idx, cov_idx, x_idx, y_idx = 0,1,2,3,4,5,6,7,8,9,10,11

        
        def get_prediction(fit_model,num_vertex):
            
            # get data time course
            tc_data_mat = self.tc_mat[num_vertex,:]
            

            # # get model time course
            if fit_model == 'gauss' or fit_model == 'gauss_sg':
                tc_model_mat = self.model_func.generate_prediction( 
                                                    x = self.deriv_mat[num_vertex,x_idx], 
                                                    y = self.deriv_mat[num_vertex,y_idx], 
                                                    sigma = self.deriv_mat[num_vertex,size_idx],
                                                    beta = self.deriv_mat[num_vertex,amp_idx], 
                                                    baseline = self.deriv_mat[num_vertex,baseline_idx])

            elif fit_model == 'css' or fit_model == 'css_sg':
                tc_model_mat = self.model_func.generate_prediction( 
                                                    x = self.deriv_mat[num_vertex,x_idx], 
                                                    y = self.deriv_mat[num_vertex,y_idx], 
                                                    sigma = self.deriv_mat[num_vertex,size_idx],
                                                    beta = self.deriv_mat[num_vertex,amp_idx], 
                                                    n = self.deriv_mat[num_vertex,non_lin_idx], 
                                                    baseline = self.deriv_mat[num_vertex,baseline_idx])
        
            deriv_model_mat = self.deriv_mat[num_vertex,:]

            return (tc_data_mat, tc_model_mat,deriv_model_mat)

        # Time course - high parameter
        # ---------------------------
        
        
        # get model and data time course
        if self.num_vertex[1] != -1:
            deb()
            tc_data_mat,tc_model_mat,deriv_model_mat = get_prediction(fit_model = self.fit_model,num_vertex = self.num_vertex[1])
            low_val = np.nanpercentile(tc_data_mat,5)*1.5
            if np.round(low_val,0): low_val_dec_round = 1
            elif np.round(low_val,1): low_val_dec_round = 2
            else: low_val_dec_round = 0

            high_val = np.nanpercentile(tc_data_mat,95)*1.5
            if np.round(high_val,0): high_val_dec_round = 1
            elif np.round(high_val,1): high_val_dec_round = 2
            else: high_val_dec_round = 0
            y_range_tc = (round(Decimal(low_val),low_val_dec_round),round(Decimal(high_val),high_val_dec_round))
        
            # get data
            x_data = np.arange(1,tc_data_mat.shape[0]+1,1)*self.tr_dur
            y_data = tc_data_mat
            x_model = np.arange(1,tc_model_mat.shape[0]+1,1)*self.tr_dur
            y_model = tc_model_mat

            high_param_tc_data_source = { 'x_data':x_data,
                                          'y_data':y_data,
                                          'x_model':x_model,
                                          'y_model':y_model}
            high_param_tc_source = ColumnDataSource(data = high_param_tc_data_source)
        else:
            y_range_tc = (-1,5)

        high_param_tc_fig              =   figure(
                                                plot_width          =   self.p_width,
                                                plot_height         =   int(self.p_height*0.42),
                                                x_range             =   self.x_range_tc,
                                                y_range             =   y_range_tc,
                                                min_border_left     =   self.min_border_large,
                                                min_border_right    =   self.min_border_large,
                                                min_border_bottom   =   self.min_border_large,
                                                min_border_top      =   self.min_border_large,
                                                toolbar_location    =   None,
                                                tools               =   "")
        high_param_tc_fig.xaxis.axis_label       =   ''
        high_param_tc_fig.yaxis.axis_label       =   self.y_label_tc
        high_param_tc_fig.grid.grid_line_color   =   None
        high_param_tc_fig.axis.minor_tick_in     =   False
        high_param_tc_fig.axis.minor_tick_out    =   False
        high_param_tc_fig.axis.major_tick_in     =   False
        high_param_tc_fig.outline_line_alpha     =   0
        high_param_tc_fig.xaxis.ticker           =   np.arange(self.x_range_tc[0],self.x_range_tc[1] + self.x_tick_tc, self.x_tick_tc)
        high_param_tc_fig.background_fill_color  =   self.bg_color
        high_param_tc_fig.axis.axis_label_standoff = 10
        high_param_tc_fig.axis.axis_label_text_font_style = 'normal'
        high_param_tc_fig.xaxis.major_label_text_font_size = '0pt'

        deb()
        # span
        high_param_tc_fig.add_layout(Span(location = 0, dimension = 'width', line_alpha = 0.5, line_color = 'black', line_width = 1, line_dash = 'dashed'))

        if self.num_vertex[1] != -1:
            # plot data
            high_param_tc_fig.circle(x = 'x_data', y = 'y_data',fill_color = 'black',line_color = 'black',line_width = 1,size = 2,source=high_param_tc_source, legend = "data")

            # plot data model
            high_param_tc_plot_model = high_param_tc_fig.line(x = 'x_model', y = 'y_model', line_width = 2, line_color = self.model_line_color, source = high_param_tc_source, legend = "model")

            # data hover
            high_param_tc_fig_tooltips = [  ('data: ',  ' @x_data{0} s, @y_data{0.0} %'),
                                            ('model: ',  '@x_model{0} s, @y_model{0.0} %')
                                        ]                               # coverage
            high_param_tc_fig_hover = HoverTool(tooltips = high_param_tc_fig_tooltips,
                                                mode = 'vline',
                                                renderers = [high_param_tc_plot_model])
            high_param_tc_fig.add_tools(high_param_tc_fig_hover)

            # legend
            high_param_tc_fig.legend.location = "top_right"
            high_param_tc_fig.legend.click_policy="hide"
            high_param_tc_fig.legend.background_fill_alpha = 0
            high_param_tc_fig.legend.label_text_font = '8pt'
            high_param_tc_fig.legend.margin = 5
            high_param_tc_fig.legend.padding = 5
            high_param_tc_fig.legend.spacing = -2
            high_param_tc_fig.legend.glyph_width = 10
            high_param_tc_fig.legend.label_text_baseline = 'middle'
            high_param_tc_fig.legend.border_line_color = None

            


        # text
        x_text = self.x_range_tc[0] + (self.x_range_tc[1]-self.x_range_tc[0])*0.03
        y_text = y_range_tc[1] - (y_range_tc[1]-y_range_tc[0])*0.11
        
        text = '{val_r2} r2 + High {params}: {title}'.format(val_r2 = self.r2_level, params = self.params,title = self.title)
        high_param_tc_fig.text(x=x_text,y=y_text,text = [text],text_font_size = '8pt',text_font_style = 'bold')

        # pRF map - low parameter
        # -----------------------
        high_param_map_fig              =   figure(
                                                plot_width          =   int(self.p_height/2),
                                                plot_height         =   int(self.p_height*0.42),
                                                x_range             =   self.x_range_map,
                                                y_range             =   self.y_range_map,
                                                min_border_left     =   self.min_border_large,
                                                min_border_right    =   self.min_border_large,
                                                min_border_bottom   =   self.min_border_large,
                                                min_border_top      =   self.min_border_large,
                                                toolbar_location    =   None,
                                                tools               =   "")

        high_param_map_fig.xaxis.axis_label       =   ''
        high_param_map_fig.yaxis.axis_label       =   self.y_label_map
        high_param_map_fig.grid.grid_line_color   =   None
        high_param_map_fig.axis.minor_tick_in     =   False
        high_param_map_fig.axis.minor_tick_out    =   False
        high_param_map_fig.axis.major_tick_in     =   False
        high_param_map_fig.outline_line_alpha     =   0
        high_param_map_fig.yaxis.ticker           =   np.arange(self.y_range_map[0],self.y_range_map[1] + self.y_tick_map, self.y_tick_map)
        high_param_map_fig.xaxis.ticker           =   np.arange(self.x_range_map[0],self.x_range_map[1] + self.x_tick_map, self.x_tick_map)
        high_param_map_fig.background_fill_color  =   self.bg_color
        high_param_map_fig.axis.axis_label_standoff = 10
        high_param_map_fig.axis.axis_label_text_font_style = 'normal'
        high_param_map_fig.xaxis.major_label_text_font_size = '0pt'

        if self.num_vertex[1] != -1:
            # stimulus circle
            high_param_map_fig.quad(left = -self.stim_width/2.0, bottom = -self.stim_height/2.0, top = +self.stim_height/2.0, right = +self.stim_width/2.0, color = self.stim_color)

            # spans
            high_param_map_fig.add_layout(Span(location = 0, dimension = 'width', line_alpha = 0.5, line_color = 'black', line_width = 1, line_dash = 'dashed'))
            high_param_map_fig.add_layout(Span(location = 0, dimension = 'height', line_alpha = 0.5, line_color = 'black', line_width = 1, line_dash = 'dashed'))

            # plot rf
            high_param_map_fig.circle(  x                   =   deriv_model_mat[x_idx],
                                        y                   =   deriv_model_mat[y_idx],
                                        radius              =   deriv_model_mat[size_idx],
                                        fill_color          =   self.model_fill_color,
                                        line_color          =   'black',
                                        fill_alpha          =   1,
                                        line_alpha          =   1)

            # text
            x_text1 = self.x_range_map[0] + (self.x_range_map[1]-self.x_range_map[0])*0.05
            x_text2 = self.x_range_map[0] + (self.x_range_map[1]-self.x_range_map[0])*0.55
            y_text = self.y_range_map[1] - (self.y_range_map[1]-self.y_range_map[0])*0.025
            if self.fit_model == 'gauss' or self.fit_model == 'gauss_sg':
                text1 = 'r2:     \t{:1.2f}\nEcc.: \t{:1.1f} dva\nSize: \t{:1.1f} dva'.format(
                                                                                                deriv_model_mat[rsq_idx],
                                                                                                deriv_model_mat[ecc_idx],
                                                                                                deriv_model_mat[size_idx],
                                                                                                 )
                text2 = 'Cov.: \t{:1.0f} %\nAmp.: \t{:1.2f}'.format(   deriv_model_mat[cov_idx]*100,
                                                                        deriv_model_mat[amp_idx],
                                                                        )

            elif self.fit_model == 'css' or self.fit_model == 'css_sg':
                text1 = 'r2:     \t{:1.2f}\nEcc.: \t{:1.1f} dva\nSize: \t{:1.1f} dva'.format(   deriv_model_mat[rsq_idx],
                                                                                                deriv_model_mat[ecc_idx],
                                                                                                deriv_model_mat[size_idx],
                                                                                                 )
                text2 = 'n:    \t{:1.1f}\nCov.: \t{:1.0f} %\nAmp.: \t{:1.2f}'.format(  deriv_model_mat[non_lin_idx],
                                                                                        deriv_model_mat[cov_idx]*100,
                                                                                        deriv_model_mat[amp_idx],
                                                                        )

            high_param_map_fig.text(x=x_text1,y=y_text,text = [text1],text_font_size = '8pt',text_baseline = 'top')
            high_param_map_fig.text(x=x_text2,y=y_text,text = [text2],text_font_size = '8pt',text_baseline = 'top')


        # Time course - low parameter
        # ----------------------------

        # get model and data time course        
        if self.num_vertex[0] != -1:
            
            tc_data_mat,tc_model_mat,deriv_model_mat = get_prediction(fit_model = self.fit_model,num_vertex = self.num_vertex[0])
            low_val = np.nanpercentile(tc_data_mat,5)*1.5
            if np.round(low_val,0): low_val_dec_round = 1
            elif np.round(low_val,1): low_val_dec_round = 2
            else: low_val_dec_round = 0

            high_val = np.nanpercentile(tc_data_mat,95)*1.5
            if np.round(high_val,0): high_val_dec_round = 1
            elif np.round(high_val,1): high_val_dec_round = 2
            else: high_val_dec_round = 0
            y_range_tc = (round(Decimal(low_val),low_val_dec_round),round(Decimal(high_val),high_val_dec_round))
            
            
            # get data
            x_data = np.arange(1,tc_data_mat.shape[0]+1,1)*self.tr_dur
            y_data = tc_data_mat
            x_model = np.arange(1,tc_model_mat.shape[0]+1,1)*self.tr_dur
            y_model = tc_model_mat
            
            
            low_param_tc_data_source = {  'x_data':x_data,
                                          'y_data':y_data,
                                          'x_model':x_model,
                                          'y_model':y_model}
            low_param_tc_source = ColumnDataSource(data = low_param_tc_data_source)
        
        else:
            y_range_tc = (-1,5)

        low_param_tc_fig              =   figure(
                                                plot_width          =   self.p_width,
                                                plot_height         =   int(self.p_height/2),
                                                x_range             =   self.x_range_tc,
                                                y_range             =   y_range_tc,
                                                min_border_left     =   self.min_border_large,
                                                min_border_right    =   self.min_border_large,
                                                min_border_bottom   =   self.min_border_large,
                                                min_border_top      =   self.min_border_large,
                                                toolbar_location    =   None,
                                                tools               =   "")

        low_param_tc_fig.xaxis.axis_label       =   self.x_label_tc
        low_param_tc_fig.yaxis.axis_label       =   self.y_label_tc
        low_param_tc_fig.grid.grid_line_color   =   None
        low_param_tc_fig.axis.minor_tick_in     =   False
        low_param_tc_fig.axis.minor_tick_out    =   False
        low_param_tc_fig.axis.major_tick_in     =   False
        low_param_tc_fig.outline_line_alpha     =   0
        low_param_tc_fig.xaxis.ticker           =   np.arange(self.x_range_tc[0],self.x_range_tc[1] + self.x_tick_tc, self.x_tick_tc)
        low_param_tc_fig.background_fill_color  =   self.bg_color
        low_param_tc_fig.axis.axis_label_standoff = 10
        low_param_tc_fig.axis.axis_label_text_font_style = 'normal'
        low_param_tc_fig.add_layout(Span(location = 0, dimension = 'width', line_alpha = 0.5, line_color = 'black', line_width = 1, line_dash = 'dashed'))

        if self.num_vertex[0] != -1:
            # plot data
            low_param_tc_fig.circle(x = 'x_data', y = 'y_data', fill_color = 'black', line_color = 'black', line_width = 1, size = 2, source = low_param_tc_source, legend="data" )

            # plot data model
            low_param_tc_plot_model = low_param_tc_fig.line(x = 'x_model', y = 'y_model', line_width = 2, line_color = self.model_line_color, source = low_param_tc_source, legend="model")


            # data hover
            low_param_tc_fig_tooltips = [   ('data: ',  ' @x_data{0} s, @y_data{0.0} %'),
                                            ('model: ',  '@x_model{0} s, @y_model{0.0} %')
                                        ]                               # coverage
            low_param_tc_fig_hover = HoverTool(tooltips = low_param_tc_fig_tooltips,
                                                mode = 'vline',
                                                renderers = [low_param_tc_plot_model])
            low_param_tc_fig.add_tools(low_param_tc_fig_hover)

            # legend
            low_param_tc_fig.legend.location = "top_right"
            low_param_tc_fig.legend.click_policy="hide"
            low_param_tc_fig.legend.background_fill_alpha = 0
            low_param_tc_fig.legend.label_text_font = '8pt'
            low_param_tc_fig.legend.margin = 5
            low_param_tc_fig.legend.padding = 5
            low_param_tc_fig.legend.spacing = -2
            low_param_tc_fig.legend.glyph_width = 10
            low_param_tc_fig.legend.label_text_baseline = 'middle'
            low_param_tc_fig.legend.border_line_color = None

        # text
        x_text = self.x_range_tc[0] + (self.x_range_tc[1]-self.x_range_tc[0])*0.03
        y_text = y_range_tc[1] - (y_range_tc[1]-y_range_tc[0])*0.11
        text = '{val_r2} r2 + Low {params}: {title}'.format(val_r2 = self.r2_level, params = self.params,title = self.title)
        low_param_tc_fig.text(x=x_text,y=y_text,text = [text],text_font_size = '8pt',text_font_style = 'bold')

        # pRF map - low parameter
        # -----------------------
        low_param_map_fig              =   figure(                                                                  # create a figure in bokeh
                                                plot_width          =   int(self.p_height/2),                            # define figure width in pixel
                                                plot_height         =   int(self.p_height/2),                            # define figure height in pixel
                                                x_range             =   self.x_range_map,                               # define x limits
                                                y_range             =   self.y_range_map,                               # define y limits
                                                min_border_left     =   self.min_border_large,                      # define left border size
                                                min_border_right    =   self.min_border_large,                      # define right border size
                                                min_border_bottom   =   self.min_border_large,                      # define bottom border space
                                                min_border_top      =   self.min_border_large,                      # define top border space
                                                toolbar_location    =   None,
                                                tools               =   "")                                         # define tools to show

        low_param_map_fig.xaxis.axis_label       =   self.x_label_map
        low_param_map_fig.yaxis.axis_label       =   self.y_label_map
        low_param_map_fig.grid.grid_line_color   =   None
        low_param_map_fig.axis.minor_tick_in     =   False
        low_param_map_fig.axis.minor_tick_out    =   False
        low_param_map_fig.axis.major_tick_in     =   False
        low_param_map_fig.outline_line_alpha     =   0
        low_param_map_fig.yaxis.ticker           =   np.arange(self.y_range_map[0],self.y_range_map[1] + self.y_tick_map, self.y_tick_map)
        low_param_map_fig.xaxis.ticker           =   np.arange(self.x_range_map[0],self.x_range_map[1] + self.x_tick_map, self.x_tick_map)
        low_param_map_fig.background_fill_color  =   self.bg_color
        low_param_map_fig.axis.axis_label_standoff = 10
        low_param_map_fig.axis.axis_label_text_font_style = 'normal'

        if self.num_vertex[0] != -1:
            # stimulus circle
            low_param_map_fig.quad(left = -self.stim_width/2.0, bottom = -self.stim_height/2.0, top = +self.stim_height/2.0, right = +self.stim_width/2.0, color = self.stim_color)

            # spans
            low_param_map_fig.add_layout(Span(location = 0, dimension = 'width', line_alpha = 0.5, line_color = 'black', line_width = 1, line_dash = 'dashed'))
            low_param_map_fig.add_layout(Span(location = 0, dimension = 'height', line_alpha = 0.5, line_color = 'black', line_width = 1, line_dash = 'dashed'))



            # plot rf
            low_param_map_fig.circle(   x                   =   deriv_model_mat[x_idx],
                                        y                   =   deriv_model_mat[y_idx],
                                        radius              =   deriv_model_mat[size_idx],
                                        fill_color          =   self.model_fill_color,
                                        line_color          =   'black',
                                        fill_alpha          =   1,
                                        line_alpha          =   1)


            # text
            x_text1 = self.x_range_map[0] + (self.x_range_map[1]-self.x_range_map[0])*0.05
            x_text2 = self.x_range_map[0] + (self.x_range_map[1]-self.x_range_map[0])*0.55
            y_text = self.y_range_map[1] - (self.y_range_map[1]-self.y_range_map[0])*0.025
            if self.fit_model == 'gauss' or self.fit_model == 'gauss_sg':
                text1 = 'r2:     \t{:1.2f}\nEcc.: \t{:1.1f} dva\nSize: \t{:1.1f} dva'.format(
                                                                                                deriv_model_mat[rsq_idx],
                                                                                                deriv_model_mat[ecc_idx],
                                                                                                deriv_model_mat[size_idx],
                                                                                                 )
                text2 = 'Cov.: \t{:1.0f} %\nAmp.: \t{:1.1f}'.format(   deriv_model_mat[cov_idx]*100,
                                                                       deriv_model_mat[amp_idx],
                                                                        )

            elif self.fit_model == 'css' or self.fit_model == 'css_sg':
                text1 = 'r2:     \t{:1.2f}\nEcc.: \t{:1.1f} dva\nSize: \t{:1.2f} dva'.format(   deriv_model_mat[rsq_idx],
                                                                                                deriv_model_mat[ecc_idx],
                                                                                                deriv_model_mat[size_idx],
                                                                                                 )
                text2 = 'n:    \t{:1.1f}\nCov.: \t{:1.0f} %\nAmp.: \t{:1.2f}'.format(  deriv_model_mat[non_lin_idx],
                                                                                        deriv_model_mat[cov_idx]*100,
                                                                                        deriv_model_mat[amp_idx],
                                                                        )

            low_param_map_fig.text(x=x_text1,y=y_text,text = [text1],text_font_size = '8pt',text_baseline = 'top')
            low_param_map_fig.text(x=x_text2,y=y_text,text = [text2],text_font_size = '8pt',text_baseline = 'top')


        # Time course stimuli legend
        # --------------------------
        time_leg_fig              =   figure(   plot_width          =   self.p_width,
                                                plot_height         =   int(self.p_height/8),
                                                x_range             =   self.x_range_tc,
                                                y_range             =   (0,1),
                                                x_axis_type         =   None,
                                                y_axis_type         =   None,
                                                outline_line_color  =   "white",
                                                min_border_left     =   self.min_border_large,
                                                min_border_right    =   self.min_border_large,
                                                min_border_bottom   =   self.min_border_large,
                                                min_border_top      =   self.min_border_large,
                                                toolbar_location    =   None)

        stim_on = ([16,30],[31,52],[68,89],[90,104])
        stim_off = ([0,15],[53,67],[105,119])
        stim_dir = (['down'],['left'],['right'],['up'])

        for t_stim_on in np.arange(0,len(stim_on),1):
            time_leg_fig.quad(left=stim_on[t_stim_on][0]*self.tr_dur, right=stim_on[t_stim_on][1]*self.tr_dur, top=1, bottom=0, fill_color="black",line_color = 'white',line_width = 1)
            x_dir_txt = (stim_on[t_stim_on][1]*self.tr_dur+stim_on[t_stim_on][0]*self.tr_dur)/2.0
            time_leg_fig.text(x = x_dir_txt, y=0.5,text = ['{text}'.format(text = stim_dir[t_stim_on][0])],text_align = 'center',text_font_size = '8pt',text_color = 'white',text_baseline = 'middle')

        for t_stim_off in np.arange(0,len(stim_off),1):
            time_leg_fig.quad(left=stim_off[t_stim_off][0]*self.tr_dur, right=stim_off[t_stim_off][1]*self.tr_dur, top=1, bottom=0, fill_color=self.bg_color,line_color = 'white',line_width = 1)

        # up-right space
        s2 = Spacer(width=int(self.p_height/2), height=int(self.p_height/10))


        # Put figure together
        # -------------------
        f   =   column( row(time_leg_fig,s2),
                        row(high_param_tc_fig,high_param_map_fig),
                        row(low_param_tc_fig,low_param_map_fig))
        main_fig = high_param_tc_fig

        # save in svg
        # -----------
        if self.save_svg == 1:
            from bokeh.io import export_svgs
            import os
            import time
            opj = os.path.join
                  
            try: os.makedirs(opj(self.svg_folder,self.svg_subfolder,self.svg_filename))
            except: pass

            fig_dict = {'time_leg_fig': time_leg_fig, 'high_param_tc_fig': high_param_tc_fig, 'high_param_map_fig':high_param_map_fig,
                        'low_param_tc_fig': low_param_tc_fig, 'low_param_map_fig': low_param_map_fig}

            # save figure
            for fig in fig_dict:
                fig_dict[fig].output_backend = 'svg'
                output_file_svg = opj(self.svg_folder,self.svg_subfolder,self.svg_filename,"{svg_filename}_{fig}.svg".format(svg_filename = self.svg_filename, fig = fig))
                
                try: export_svgs(fig_dict[fig], filename = output_file_svg)
                except: pass

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
        elif plot == 'lat':
            f, main_fig = self.draw_pRFlat(params = parameters, old_main_fig = old_main_fig)
        elif plot == 'tc':
            f, main_fig = self.draw_pRFtc(params = parameters)
        
        return (f, main_fig)

