def draw_cortex_vertex(subject,data,cmap,vmin,vmax,cbar = 'discrete',cmap_steps = 255,\
                        alpha = None,depth = 1,thick = 1,height = 1024,sampler = 'nearest',\
                        with_curvature = True,with_labels = False,with_colorbar = False,\
                        with_borders = False,curv_brightness = 0.95,curv_contrast = 0.05,add_roi = False,\
                        roi_name = 'empty',col_offset = 0, zoom_roi = None, zoom_hem = None, zoom_margin = 0.0):
    """
    Plot brain data onto a previously saved flatmap.

    Parameters
    ----------
    subject             : subject id (e.g. 'sub-001')
    data                : the data you would like to plot on a flatmap
    cmap                : colormap that shoudl be used for plotting
    vmin                : minimal value iqn colormap
    vmax                : maximal value in colormap
    cbar                : color bar layout
    cmap_steps          : number of colormap bins
    alpha               : alpha map
    depth               : Value between 0 and 1 for how deep to sample the surface for the flatmap (0 = gray/white matter boundary, 1 = pial surface)
    thick               : Number of layers through the cortical sheet to sample. Only applies for pixelwise = True
    height              : Height of the image to render. Automatically scales the width for the aspect of the subject's flatmap
    sampler             : Name of sampling function used to sample underlying volume data. Options include 'trilinear', 'nearest', 'lanczos'
    with_curvature      : Display the rois, labels, colorbar, annotated flatmap borders, or cross-hatch dropout?
    with_labels         : Display labels?
    with_colorbar       : Display pycortex' colorbar?
    with_borders        : Display borders?
    curv_brightness     : Mean brightness of background. 0 = black, 1 = white, intermediate values are corresponding grayscale values.
    curv_contrast       : Contrast of curvature. 1 = maximal contrast (black/white), 0 = no contrast (solid color for curvature equal to curvature_brightness).
    add_roi             : add roi -image- to overlay.svg
    roi_name            : roi name
    col_offset          : colormap offset between 0 and 1
    zoom_roi            : name of the roi on which to zoom on
    zoom_hem            : hemifield fo the roi zoom
    zoom_margin         : margin in mm around the zoom

    Returns
    -------
    vertex_rgb - pycortex vertex file
    """
    
    import cortex
    import pdb
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    from matplotlib import cm
    import matplotlib as mpl
    import ipdb
    deb = ipdb.set_trace
    
    # define colormap
    base = cortex.utils.get_cmap(cmap)
    val = np.fmod(np.linspace(0+col_offset, 1+col_offset,cmap_steps+1,endpoint=False),1.0)
    colmap = colors.LinearSegmentedColormap.from_list(  'my_colmap',
                                                        base(val),
                                                        N = cmap_steps)
    
    # convert data to RGB
    vrange = float(vmax) - float(vmin)
    norm_data = ((data-float(vmin))/vrange)*cmap_steps
    mat = colmap(norm_data.astype(int))
    
    
    mat_hsv = colors.rgb_to_hsv(mat[:,0:3])
    alpha[np.isnan(alpha)]=0
    mat_hsv[:,2] = mat_hsv[:,2]*alpha
    mat = colors.hsv_to_rgb(mat_hsv)


    vertex_rgb = cortex.VertexRGB(  red = mat[...,0],
                                    green = mat[...,1],
                                    blue = mat[...,2],
                                    alpha = alpha,
                                    subject = subject,
                                    )
    
    vertex_rgb_fig = cortex.quickshow(  braindata = vertex_rgb,
                                        depth = depth,
                                        thick = thick,
                                        height = height,
                                        sampler = sampler,
                                        with_curvature = with_curvature,
                                        with_labels = with_labels,
                                        with_colorbar = with_colorbar,
                                        with_borders = with_borders,
                                        curvature_brightness = curv_brightness,
                                        curvature_contrast = curv_contrast)
    
    zoom_plot = []
    if zoom_roi != None:
        roi_verts = cortex.get_roi_verts(subject, zoom_roi)[zoom_roi]
        roi_map = cortex.Vertex.empty(subject)
        roi_map.data[roi_verts] = 1

        (lflatpts, lpolys), (rflatpts, rpolys) = cortex.db.get_surf(subject, "flat", nudge=True)
        sel_pts = dict(left=lflatpts, right=rflatpts)[zoom_hem]
        roi_pts = sel_pts[np.nonzero(getattr(roi_map, zoom_hem))[0],:2]
        
        xmin, ymin = roi_pts.min(0) - zoom_margin
        xmax, ymax = roi_pts.max(0) + zoom_margin
        zoom_width = roi_pts.max(0)[0]- roi_pts.min(0)[0]
        zoom_height = roi_pts.max(0)[1]- roi_pts.min(0)[1]
        sqr_zoom_size = np.round(np.max([zoom_width,zoom_height]))*1.1
        zoom_ctr_x = np.mean([roi_pts.max(0)[0],roi_pts.min(0)[0]])
        zoom_ctr_y = np.mean([roi_pts.max(0)[1],roi_pts.min(0)[1]])
        mat = [zoom_ctr_x-sqr_zoom_size/2.0,zoom_ctr_x+sqr_zoom_size/2.0,zoom_ctr_y-sqr_zoom_size/2.0,zoom_ctr_y+sqr_zoom_size/2.0]
        zoom_plot = pl.axis(mat)
    else:
        
        if cbar == 'polar':

            # Polar angle color bar
            colorbar_location = [0.5, 0.07, 0.8, 0.2]
            n = 200
            cbar_axis = vertex_rgb_fig.add_axes(colorbar_location, projection='polar')
            norm = mpl.colors.Normalize(0, 2*np.pi)

            # Plot a color mesh on the polar plot
            # with the color set by the angle
            t = np.linspace(2*np.pi,0,n)
            r = np.linspace(1,0,2)
            rg, tg = np.meshgrid(r,t)
            c = tg
            im = cbar_axis.pcolormesh(t, r, c.T,norm= norm, cmap = colmap)
            cbar_axis.set_theta_zero_location("W")
            cbar_axis.set_yticklabels([])
            cbar_axis.set_xticklabels([])
            cbar_axis.spines['polar'].set_visible(False)

        elif cbar == 'ecc':
        
            # Ecc color bar
            colorbar_location = [0.5, 0.07, 0.8, 0.2]
            n = 200
            cbar_axis = vertex_rgb_fig.add_axes(colorbar_location, projection='polar')

            t = np.linspace(0,2*np.pi, n)
            r = np.linspace(0,1, n)
            rg, tg = np.meshgrid(r,t)
            c = tg
            
            im = cbar_axis.pcolormesh(t, r, c, norm = mpl.colors.Normalize(0, 2*np.pi), cmap = colmap)
            cbar_axis.tick_params(pad = 1,labelsize = 15)
            cbar_axis.spines['polar'].set_visible(False)
            
            # superimpose new axis for dva labeling
            box = cbar_axis.get_position()
            cbar_axis.set_yticklabels([])
            cbar_axis.set_xticklabels([])
            axl = vertex_rgb_fig.add_axes(  [1.8*box.xmin,
                                        0.5*(box.ymin+box.ymax),
                                        box.width/600,
                                        box.height*0.5],
                                        axisbg = None)
            axl.spines['top'].set_visible(False)
            axl.spines['right'].set_visible(False)
            axl.spines['bottom'].set_visible(False)
            axl.yaxis.set_ticks_position('right')
            axl.xaxis.set_ticks_position('none')
            axl.set_xticklabels([])
            axl.set_yticklabels(np.linspace(vmin,vmax,3),size = 'x-large')
            axl.set_ylabel('$dva$\t\t', rotation = 0, size = 'x-large')
            axl.yaxis.set_label_coords(box.xmax+30,0.4)
            axl.patch.set_alpha(0.5)

        elif cbar == 'discrete':

            # Discrete color bars
            # -------------------
            colorbar_location= [0.9, 0.05, 0.03, 0.25]
            cmaplist = [colmap(i) for i in range(colmap.N)]

            # define the bins and normalize
            bounds = np.linspace(vmin, vmax, cmap_steps + 1)
            bounds_label = np.linspace(vmin, vmax, 3)
            norm = mpl.colors.BoundaryNorm(bounds, colmap.N)
            
            cbar_axis = vertex_rgb_fig.add_axes(colorbar_location)
            cb = mpl.colorbar.ColorbarBase(cbar_axis,cmap = colmap,norm = norm,ticks = bounds_label,boundaries = bounds)

        # add to overalt
        if add_roi == True:
            cortex.utils.add_roi(   data = vertex_rgb,
                                    name = roi_name,
                                    open_inkscape = False,
                                    add_path = False,
                                    depth = depth,
                                    thick = thick,
                                    sampler = sampler,
                                    with_curvature = with_curvature,
                                    with_colorbar = with_colorbar,
                                    with_borders = with_borders,
                                    curvature_brightness = curv_brightness,
                                    curvature_contrast = curv_contrast)

    return vertex_rgb


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
opj = os.path.join
deb = ipdb.set_trace

# MRI imports
# -----------
import nibabel as nb
import cortex

subject = '999999'

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

fit_model = 'gauss'
deriv_dir = opj(base_dir,'pp_data',subject,fit_model,'deriv')

# Create derivatives flatmaps
# ---------------------------
cmap_neg_pos = 'BuBkRd'#'RdBu_r'
sign_idx, rsq_idx, ecc_idx, polar_real_idx, polar_imag_idx , size_idx, \
            non_lin_idx, amp_idx, baseline_idx, cov_idx, x_idx, y_idx = 0,1,2,3,4,5,6,7,8,9,10,11
mask_dir  = 'all'

# Get data and combine hemispheres
deriv_mat=[]
for hemi in ['L','R']:
    deriv_file = nb.load(opj(deriv_dir,mask_dir,"prf_deriv_{hemi}_{mask_dir}_fsaverage.func.gii".format(hemi = hemi, mask_dir = mask_dir)))
    deriv_mat.append(np.array([deriv_file.darrays[i].data for i in range(len(deriv_file.darrays))]))
deriv_mat = np.hstack(deriv_mat)
deriv_mat = deriv_mat.T

threshold_mask = np.logical_and(np.logical_and( deriv_mat[:,rsq_idx]>=0.0,
                                                deriv_mat[:,cov_idx]>=0.0),
                                                deriv_mat[:,size_idx]>=0.1)

# R-square
rsq_data = deriv_mat[:,rsq_idx]
rsq_data[~threshold_mask] = np.nan
alpha = rsq_data

# Sign
sign_data = deriv_mat[:,sign_idx]

param_sign = {'data': sign_data.T, 'cmap': cmap_neg_pos, 'alpha': alpha.T, 'vmin': -1, 'vmax': 1, 'cbar': 'discrete','height': 2160,
              'curv_brightness': 1, 'curv_contrast': 1, 'subject': 'fsaverage','add_roi': False, 'with_labels': False}
vertex_sign = draw_cortex_vertex(**param_sign)

plt.savefig('/home/szinte/projects/retino_HCP/flat_sign.pdf',facecolor = "w")
