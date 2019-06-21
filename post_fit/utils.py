import os
import nibabel as nb
import numpy as np
import matplotlib.pyplot as pl
import scipy as sp
from scipy.signal import savgol_filter
from skimage.transform import rotate
from math import *
import cortex
import cifti

def set_pycortex_config_file(project_folder):

    
    # Import necessary modules
    import os
    import cortex
    from pathlib import Path

    import ipdb

    # Define the new database and colormaps folder
    pycortex_db_folder = project_folder + '/db/'
    pycortex_cm_folder = project_folder + '/colormaps/'

    # Get pycortex config file location
    pycortex_config_file  = cortex.options.usercfg

    # Create name of new config file that will be written
    new_pycortex_config_file = pycortex_config_file[:-4] + '_new.cfg'

    # Create the new config file
    Path(new_pycortex_config_file).touch()

    # Open the config file in read mode and the newly created one in write mode.
    # Loop over every line in the original file and copy it into the new one.
    # For the lines containing either 'filestore' or 'colormap', it will
    # change the saved folder path to the newly created one above (e.g. pycortex_db_folder)
    with open(pycortex_config_file, 'r') as fileIn:
        with open(new_pycortex_config_file, 'w') as fileOut:

            for line in fileIn:

                if 'filestore' in line:
                    newline = 'filestore=' + pycortex_db_folder
                    fileOut.write(newline)
                    newline = '\n'

                elif 'Colormaps' in line:
                    newline = 'Colormaps=' + pycortex_cm_folder
                    fileOut.write(newline)
                    newline = '\n'

                else:
                    newline = line

                fileOut.write(newline)

    
    # Renames the original config file als '_old' and the newly created one to the original name
    os.rename(pycortex_config_file, pycortex_config_file[:-4] + '_old.cfg')
    os.rename(new_pycortex_config_file, pycortex_config_file)
    return None

def convert_fit_results(prf_filename, 
                        output_dir, 
                        stim_radius,
                        fit_model):
    """
    Convert pRF fitting value in different parameters for following analysis
   
    Parameters
    ----------
    prf_filename: absolute paths to prf result files.
    output_dir: absolute path to directory into which to put the resulting files.
    stim_radius: stimulus radius in deg
    fit_model: fit model ('gauss','css')

    Returns
    -------
    prf_deriv_all: derivative of pRF analysis for all pRF voxels
    prf_deriv_neg: derivative of pRF analysis for all negative pRF voxels
    prf_deriv_pos: derivative of pRF analysis for all positive pRF voxels

    stucture output:
    columns: 1->32492
    row00 : sign
    row01 : R2
    row02 : eccentricity in deg
    row03 : polar angle real component in deg
    row04 : polar angle imaginary component in deg
    row05 : size in deg
    row06 : non-linerity or nans
    row07 : amplitude
    row08 : baseline
    row09 : coverage
    row10 : x
    row11 : y
    
    ['prf_sign','prf_rsq','prf_ecc','prf_polar_real','prf_polar_imag','prf_size','prf_non_lin','prf_amp','prf_baseline','prf_cov','prf_x','prf_y']

    """

    # Imports
    # -------
    # General imports
    import os
    import nibabel as nb
    import glob
    import numpy as np
    import ipdb
    deb = ipdb.set_trace
  

    # Popeye imports
    from popeye.spinach import generate_og_receptive_fields


    # Create folders
    # --------------
    try:
        os.makedirs(os.path.join(output_dir,'all'))
        os.makedirs(os.path.join(output_dir,'pos'))
        os.makedirs(os.path.join(output_dir,'neg'))
    except:
        pass

    # Get data details
    # ----------------
    prf_data_file = cifti.read(prf_filename)
    prf_data = prf_data_file[0]

    # Compute derived measures from prfs
    # ----------------------------------
    # get data index
    if fit_model == 'gauss':
        x_idx, y_idx, sigma_idx, beta_idx, baseline_idx, rsq_idx = 0, 1, 2, 3, 4, 5
    elif fit_model == 'css':
        x_idx, y_idx, sigma_idx, non_lin_idx, beta_idx, baseline_idx, rsq_idx = 0, 1, 2, 3, 4, 5, 6
    
    # pRF sign
    prf_sign_all = np.sign((prf_data[beta_idx,:]))
    pos_mask = prf_sign_all > 0.0
    neg_mask = prf_sign_all < 0.0
    all_mask = pos_mask | neg_mask
    
    # r-square
    prf_rsq_all = prf_data[rsq_idx,:]

    # pRF eccentricity
    prf_ecc_all = np.nan_to_num(np.sqrt(prf_data[x_idx,:]**2 + prf_data[y_idx,:]**2))

    # pRF polar angle
    complex_polar = prf_data[x_idx,:] + 1j * prf_data[y_idx,:]
    normed_polar = complex_polar / np.abs(complex_polar)
    prf_polar_real_all = np.real(normed_polar)
    prf_polar_imag_all = np.imag(normed_polar)
    
    # pRF size
    prf_size_all = prf_data[sigma_idx,:].astype(np.float64)
    prf_size_all[prf_size_all<1e-4] = 1e-4

    # pRF non-linearity
    if fit_model == 'gauss':
        prf_non_lin_all = np.zeros((prf_size_all.shape))*np.nan
    elif fit_model == 'css':
        prf_non_lin_all = prf_data[non_lin_idx,:]

    # pRF amplitude
    if fit_model == 'gauss':
        prf_amp_all = prf_data[beta_idx,:]
    elif fit_model == 'css':
        prf_amp_all = prf_data[beta_idx,:]

    # pRF baseline
    if fit_model == 'gauss':
        prf_baseline_all = prf_data[baseline_idx,:]
    elif fit_model == 'css':
        prf_baseline_all = prf_data[baseline_idx,:]

    # pRF coverage
    deg_x, deg_y = np.meshgrid(np.linspace(-30, 30, 121), np.linspace(-30, 30, 121))         # define prfs in visual space
    
    rfs = generate_og_receptive_fields( prf_data[x_idx,:],
                                        prf_data[y_idx,:],
                                        prf_size_all,
                                        np.ones(np.prod(prf_data[0,:].shape[0])),
                                        deg_x,
                                        deg_y)
    if fit_model == 'css':
        rfs = rfs ** prf_data[non_lin_idx,:]

    total_prf_content = rfs.reshape((-1, prf_data.shape[1])).sum(axis=0)
    stim_vignet = np.sqrt(deg_x ** 2 + deg_y**2) < stim_radius    
    prf_cov_all = rfs[stim_vignet, :].sum(axis=0) / total_prf_content

    # pRF x
    prf_x_all = prf_data[x_idx,:]

    # pRF y
    prf_y_all = prf_data[y_idx,:]

    # Saving
    # ------
    for mask_dir in ['all','pos','neg']:
        print('saving: %s'%('os.path.join(output_dir,"{mask_dir}","prf_deriv_{mask_dir}.nii")'.format(mask_dir = mask_dir)))
        for output_type in ['prf_sign','prf_rsq','prf_ecc','prf_polar_real','prf_polar_imag','prf_size','prf_non_lin','prf_amp','prf_baseline','prf_cov','prf_x','prf_y']:
            exec('{output_type}_{mask_dir} = np.copy({output_type}_all)'.format(mask_dir = mask_dir, output_type = output_type))
            exec('{output_type}_{mask_dir}[~{mask_dir}_mask] = np.nan'.format(mask_dir = mask_dir, output_type = output_type))
        
        exec('prf_deriv_{mask_dir} = np.row_stack((prf_sign_{mask_dir},prf_rsq_{mask_dir},prf_ecc_{mask_dir},prf_polar_real_{mask_dir},\
                prf_polar_imag_{mask_dir},prf_size_{mask_dir},prf_non_lin_{mask_dir},prf_amp_{mask_dir},prf_baseline_{mask_dir},prf_cov_{mask_dir},\
                prf_x_{mask_dir},prf_y_{mask_dir}))'.format(mask_dir = mask_dir))
        
        exec('prf_deriv_{mask_dir} = prf_deriv_{mask_dir}.astype(np.float32)'.format(mask_dir = mask_dir))
        bm_full = prf_data_file[1][1]
        series = cifti.Series(start=0, step=1, size=12)
        exec('cifti.write(os.path.join(output_dir,"{mask_dir}","prf_deriv_{mask_dir}.nii"), prf_deriv_{mask_dir}, (series,bm_full))'.format(mask_dir = mask_dir))


    return None

def mask_gii_2_hdf5(in_file, mask_file, hdf5_file, folder_alias, roi_num):
    """masks data in in_file with mask in mask_file,
    to be stored in an hdf5 file

    Takes a list of 3D or 4D fMRI nifti-files and masks the
    data with all masks in the list of nifti-files mask_files.
    These files are assumed to represent the same space, i.e.
    that of the functional acquisitions. 
    These are saved in hdf5_file, in the folder folder_alias.

    Parameters
    ----------
    in_files : list
        list of absolute path to functional nifti-files.
        all nifti files are assumed to have the same ndim
    mask_file : list
        list of absolute path to mask nifti-files.
        mask_files are assumed to be 3D
    hdf5_file : str
        absolute path to hdf5 file.
    folder_alias : str
                name of the to-be-created folder in the hdf5 file.
    roi_num: roi row number

    Returns
    -------
    hdf5_file : str
        absolute path to hdf5 file.
    """

    import nibabel as nb
    import os.path as op
    import numpy as np
    import h5py
    import ipdb
    deb = ipdb.set_trace


    gii_in_data = nb.load(in_file)
    data_mat = np.array([gii_in_data.darrays[i].data for i in range(len(gii_in_data.darrays))])
    data_name = op.split(in_file)[-1].split('.gii')[0]

    gii_in_mask = nb.load(mask_file)
    mask_mat = np.array([gii_in_mask.darrays[i].data for i in range(len(gii_in_mask.darrays))])

    mask_mat = mask_mat[roi_num,:]
    mask_mat = np.round(mask_mat)
    
    roi_data = data_mat[:, mask_mat==1]
    
    try:
        h5file = h5py.File(hdf5_file, "r+")
    except:
        h5file = h5py.File(hdf5_file, "a")
    
    try:
        g_hemi = h5file.create_group(folder_alias)    
    except:
        None

    h5file.create_dataset('{folder_alias}/{data_name}'.format(folder_alias = folder_alias, data_name = data_name),data = roi_data,dtype='float32')

    return None

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


def roi_coord_mask(roi,hem,subject = 'fsaverage'):
    """
    -----------------------------------------------------------------------------------------
    roi_coord_mask(roi,hem,subject = 'fsaverage')
    -----------------------------------------------------------------------------------------
    Goal of the script:
    Ger ROI pts coordinates and mask
    -----------------------------------------------------------------------------------------
    Input(s):
    subject: subject name
    roi: roi name
    hem: hemifield
    margin: mm out of ROI
    -----------------------------------------------------------------------------------------
    Output(s):
    roi_pts: pts coordinates on flat map
    roi_mask: roi mask
    -----------------------------------------------------------------------------------------
    """

    # get data
    roi_verts = cortex.get_roi_verts(subject, roi)[roi]
    roi_map = cortex.Vertex.empty(subject)
    roi_map.data[roi_verts] = 1
    (lflatpts, lpolys), (rflatpts, rpolys) = cortex.db.get_surf(subject, "flat",nudge=True)
    sel_pts = dict(left=lflatpts, right=rflatpts)[hem]
    roi_pts = sel_pts[np.nonzero(getattr(roi_map, hem))[0],:2]
    if hem == 'left':
        roi_mask = np.hstack([roi_map.left[:],roi_map.right[:]*0])
    elif hem == 'right':
        roi_mask = np.hstack([roi_map.left[:]*0,roi_map.right[:]])
    roi_mask = roi_mask==1
    
    return(roi_pts,roi_mask)

def get_colors(data,cmap,cmap_steps,col_offset,vmin,vmax):
    """
    -----------------------------------------------------------------------------------------
    get_colors(data,cmap,cmap_steps,col_offset,vmin,vmax)
    -----------------------------------------------------------------------------------------
    Goal of the script:
    Return for a given dataset the corresponding colors in Bokeh
    -----------------------------------------------------------------------------------------
    Input(s):
    data: data from which to get colors on colormaps
    cmap: colromap
    cmpas: steps for the colormap
    col_offset: color offset on the colormap scale
    vmin: minimum corresponding value on the colormap
    vmax: maximum corresponding value on the colormap
    -----------------------------------------------------------------------------------------
    Output(s):
    colors_val_rgb: matrix of colors for bokeh
    -----------------------------------------------------------------------------------------
    """
    import matplotlib.colors as colors
    base = cortex.utils.get_cmap(cmap)
    val = np.fmod(np.linspace(0 + col_offset, 1 + col_offset, cmap_steps + 1, endpoint = False), 1.0)
    colmap = colors.LinearSegmentedColormap.from_list('my_colmap', base(val), N = cmap_steps)
    vrange = float(vmax) - float(vmin)
    norm_data = (( data - float(vmin) ) / vrange) * cmap_steps 
    col_mat_rgb = colmap(norm_data.astype(int)) * 255.0
    colors_val_rgb = ["#%02x%02x%02x" % (int(r), int(g), int(b)) for r, g, b in zip(col_mat_rgb[:,0], col_mat_rgb[:,1], col_mat_rgb[:,2])]
    return colors_val_rgb

def rotate_pts(pts,orig,rot_deg):
    """
    -----------------------------------------------------------------------------------------
    rotate_pts(pts,orig,rot_deg)
    -----------------------------------------------------------------------------------------
    Goal of the script:
    rotate pts aroud an origin
    -----------------------------------------------------------------------------------------
    Input(s):
    pts: coordinates to rotate
    rot_deg: rotation amount
    -----------------------------------------------------------------------------------------
    Output(s):
    rot_pts: rotated points
    -----------------------------------------------------------------------------------------
    """
    rot = np.radians(rot_deg)
    qx = orig[0] + np.cos(rot) * (pts[0] - orig[0]) + np.sin(rot) * (pts[1] - orig[1])
    qy = orig[1] + -np.sin(rot) * (pts[0] - orig[0]) + np.cos(rot) * (pts[1] - orig[1])
    rot_pts = [qx,qy]
    return rot_pts

def rot_coord(coord,rot_deg):
    """
    -----------------------------------------------------------------------------------------
    rot_coord(coord,rot_deg)
    -----------------------------------------------------------------------------------------
    Goal of the script:
    Rotate complex numbers with rotation matrix
    -----------------------------------------------------------------------------------------
    Input(s):
    coord: complex numbers coordinate set
    rot_deg: rotation in degrees
    -----------------------------------------------------------------------------------------
    Output(s):
    coord_rot: rotated coordinates
    -----------------------------------------------------------------------------------------
    """
    theta = np.radians(rot_deg)
    c, s = np.cos(theta), np.sin(theta)
    R = np.matrix([[c, -s], [s, c]])
    coord_rot = np.array(np.dot(coord,R))
    
    return coord_rot