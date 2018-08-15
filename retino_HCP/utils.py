import os
import nibabel as nb
import numpy as np
import matplotlib.pyplot as pl
import scipy as sp
from scipy.signal import savgol_filter
from skimage.transform import rotate
from math import *

#################################################################################
#####
#####   cii workflow
#####
#################################################################################

def sg_filter_cii(cii_file, polyorder=3, deriv=0, window_length=210, tr=1):
    """sg_filter_cii temporally filters cii file contents"""

    # cii data
    cii_in = nb.load(cii_file)

    window = np.int(window_length / tr)
    if window % 2 == 0:
        window += 1

    data = cii_in.get_data()
    data_filt = savgol_filter(data.T, window_length=window, polyorder=polyorder,
                              deriv=deriv, axis=1, mode='nearest').T
    data_filt = data - data_filt + data_filt.mean(axis=0)

    cii_out = nb.Cifti2Image(dataobj=data_filt, 
                            header=cii_in.header, 
                            nifti_header=cii_in.nifti_header, 
                            extra=cii_in.extra)

    out_name = os.path.splitext(cii_file)[0] + '_sg.nii'
    out_file = os.path.abspath(out_name)
    nb.save(cii_out, out_file)

    return out_file

def psc_cii(cii_file, method='median'):
    """psc_cii performs percent signal change conversion on cii file contents"""

    # cii data
    cii_in = nb.load(cii_file)

    data = cii_in.get_data()
    if method == 'mean':
        data_m = data.mean(axis=0)
    elif method == 'median':
        data_m = np.median(data, axis=0)

    data_conv = 100.0 * (data - data_m)/data_m

    cii_out = nb.Cifti2Image(dataobj=data_conv, 
                            header=cii_in.header, 
                            nifti_header=cii_in.nifti_header, 
                            extra=cii_in.extra)

    out_name = os.path.splitext(cii_file)[0] + '_psc.nii'
    out_file = os.path.abspath(out_name)
    nb.save(cii_out, out_file)

    return out_file

def sg_psc_cii(cii_file):
    sg_file = sg_filter_cii(cii_file)
    psc_file = psc_cii(sg_file)

    return psc_file

def average_phase_encoded_ciis(file_1, file_2, shift=9):
    """second data file will be reversed and shifted by twice the haemodynamic delay
    """
    # cii data
    ciis = [nb.load(cii_file) for cii_file in [file_1, file_2]]
    data = [cii_in.get_data() for cii_in in ciis]

    data[1] = data[1][::-1]
    data[1] = np.roll(data[1], shift, axis=0)

    m_data = np.mean(np.array(data), axis = 0)

    cii_out = nb.Cifti2Image(dataobj=m_data, 
                            header=ciis[0].header, 
                            nifti_header=ciis[0].nifti_header, 
                            extra=ciis[0].extra)

    out_name = os.path.splitext(file_1)[0] + '_av.nii'
    out_file = os.path.abspath(out_name)
    nb.save(cii_out, out_file)

    return out_file

def average_bar_ciis(file_1, file_2):
    """no reversal or shifting necessary for the bar files
    """
    # cii data
    ciis = [nb.load(cii_file) for cii_file in [file_1, file_2]]
    data = [cii_in.get_data() for cii_in in ciis]

    m_data = np.mean(np.array(data), axis = 0)

    cii_out = nb.Cifti2Image(dataobj=m_data, 
                            header=ciis[0].header, 
                            nifti_header=ciis[0].nifti_header, 
                            extra=ciis[0].extra)

    out_name = (os.path.splitext(file_1)[0] + '_av.nii').replace('BAR1', 'BOTHBARS')
    out_file = os.path.abspath(out_name)
    nb.save(cii_out, out_file)

    return out_file


#################################################################################
#####
#####   gii workflow
#####
#################################################################################

def sg_filter_gii(gii_file, polyorder=3, deriv=0, window_length=210, tr=1):
    """sg_filter_gii temporally filters gii file contents"""

    # gii data
    gii_in = nb.load(gii_file)

    window = np.int(window_length / tr)
    if window % 2 == 0:
        window += 1

    data = np.array([gii_in.darrays[i].data for i in range(len(gii_in.darrays))])

    data_filt = savgol_filter(data.T, window_length=window, polyorder=polyorder,
                              deriv=deriv, axis=1, mode='nearest').T
    data_filt = data - data_filt + data_filt.mean(axis=0)

    darrays = [nb.gifti.gifti.GiftiDataArray(d) for d in data_filt]

    gii_out = nb.gifti.gifti.GiftiImage(header=gii_in.header, 
                            extra=gii_in.extra,
                            darrays=darrays)

    out_name = os.path.splitext(gii_file)[0] + '_bla.gii'
    out_file = os.path.abspath(out_name)
    nb.save(gii_out, out_file)

    return out_file

def psc_gii(gii_file, method='median'):
    """psc_gii performs percent signal change conversion on gii file contents"""

    # gii data
    gii_in = nb.load(gii_file)

    data = np.array([gii_in.darrays[i].data for i in range(len(gii_in.darrays))])
    if method == 'mean':
        data_m = data.mean(axis=0)
    elif method == 'median':
        data_m = np.median(data, axis=0)

    data_conv = 100.0 * (data - data_m)/data_m

    darrays = [nb.gifti.gifti.GiftiDataArray(d) for d in data_conv]

    gii_out = nb.gifti.gifti.GiftiImage(header=gii_in.header, 
                            extra=gii_in.extra,
                            darrays=darrays)

    out_name = os.path.splitext(gii_file)[0] + '_psc.gii'
    out_file = os.path.abspath(out_name)
    nb.save(gii_out, out_file)

    return out_file

def sg_psc_gii(gii_file):
    sg_file = sg_filter_gii(gii_file)
    psc_file = psc_gii(sg_file)

    return psc_file

def average_phase_encoded_giis(file_1, file_2, shift=9):
    """second data file will be reversed and shifted by twice the haemodynamic delay
    """
    # gii data
    giis = [nb.load(gii_file) for gii_file in [file_1, file_2]]

    data = []
    for gii in giis:
        data.append(np.array([gii.darrays[i].data for i in range(len(gii.darrays))]))

    data[1] = data[1][::-1]
    data[1] = np.roll(data[1], shift, axis=0)

    m_data = np.mean(np.array(data), axis = 0)

    darrays = [nb.gifti.gifti.GiftiDataArray(d) for d in m_data]

    gii_out = nb.gifti.gifti.GiftiImage(header=giis[0].header, 
                            extra=giis[0].extra,
                            darrays=darrays)

    out_name = os.path.splitext(file_1)[0] + '_av.gii'
    out_file = os.path.abspath(out_name)
    nb.save(gii_out, out_file)

    return out_file

def average_bar_giis(file_1, file_2):
    """no reversal or shifting necessary for the bar files
    """
    # gii data
    giis = [nb.load(gii_file) for gii_file in [file_1, file_2]]

    data = []
    for gii in giis:
        data.append(np.array([gii.darrays[i].data for i in range(len(gii.darrays))]))

    m_data = np.mean(np.array(data), axis = 0)
    darrays = [nb.gifti.gifti.GiftiDataArray(d) for d in m_data]

    gii_out = nb.gifti.gifti.GiftiImage(header=giis[0].header, 
                            extra=giis[0].extra,
                            darrays=darrays)

    out_name = os.path.splitext(file_1)[0] + '_av.gii'
    out_file = os.path.abspath(out_name)
    nb.save(gii_out, out_file)

    return out_file


#########################################################################################################################
## 
## Design matrices for different runs
##
#########################################################################################################################

def design_matrix_wedge(direction='CW',
                        pre_post_blank=22, 
                        n_steps=32, 
                        n_reps=8, 
                        size=0.25, 
                        n_pix=100):
    """design_matrix_wedge creates a wedge design matrix"""

    X,Y = np.meshgrid(np.linspace(-1, 1, n_pix, endpoint=True), np.linspace(-1, 1, n_pix, endpoint=True))
    ecc_mask = np.sqrt(X**2+Y**2) <= 1.01

    dm = np.zeros((n_steps*n_reps+2*pre_post_blank, n_pix, n_pix), dtype=bool)

    start_wedge = (X >= 0) & (Y >= 0) * ecc_mask

    if direction == 'CW':
        dir_multiply = 1
    elif direction == 'CCW':
        dir_multiply = -1
    rotation_angles = dir_multiply * 360.0 * np.arange(0,n_steps) / n_steps
    one_cycle = np.array([rotate(start_wedge, rs) for rs in rotation_angles])

    for rep in range(n_reps):
        these_tps = pre_post_blank + rep*n_steps
        dm[these_tps:these_tps+n_steps,:,:] = one_cycle
    
    return dm.astype(bool)

def design_matrix_ring(direction='EXP',
                        pre_post_blank=22, 
                        n_stim_steps=28, 
                        n_blank_steps=4,
                        n_reps=8, 
                        size_slope=0.5, 
                        n_pix=100):
    """design_matrix_ring creates a ring design matrix"""

    total_steps = (n_stim_steps+n_blank_steps)

    X,Y = np.meshgrid(np.linspace(-1, 1, n_pix, endpoint=True), np.linspace(-1, 1, n_pix, endpoint=True))
    ecc = np.sqrt(X**2+Y**2)

    ecc_pos_steps = np.linspace(0.75/n_stim_steps,0.75,n_stim_steps, endpoint=True)
    ecc_size_steps = np.linspace(size_slope/n_stim_steps,size_slope,n_stim_steps, endpoint=True)

    dm = np.zeros((total_steps*n_reps+2*pre_post_blank, n_pix, n_pix), dtype=bool)

    if direction == 'EXP':
        dir_indices = np.arange(n_stim_steps)
    elif direction == 'CON':
        dir_indices = np.arange(n_stim_steps)[::-1]

    one_cycle = np.zeros((total_steps, n_pix, n_pix))
    for i, p,s in zip(np.arange(n_stim_steps), ecc_pos_steps, ecc_size_steps):
        s2 = s/2.0
        one_cycle[i] = (ecc >= (p-s2)) & (ecc <= (p+s2))

    # invert direction for CON
    one_cycle[:n_stim_steps] = one_cycle[dir_indices]

    for rep in range(n_reps):
        these_tps = pre_post_blank + rep*total_steps
        dm[these_tps:these_tps+total_steps,:,:] = one_cycle
    
    return dm.astype(bool)

def design_matrix_prf(pre_post_blank=16, 
                        n_stim_steps=28, 
                        n_blank_steps=4,
                        n_reps=8, 
                        bar_width=0.25, # of range 2 instead of 1
                        n_pix=100,
                        intermezzo_blank_steps=12,
                        directions=[0,-90,-180,-270,-45,-135,-225,-315]):

    s2 = bar_width/2.0
    total_steps = (n_stim_steps+n_blank_steps)

    X,Y = np.meshgrid(np.linspace(-1, 1, n_pix, endpoint=True), np.linspace(-1, 1, n_pix, endpoint=True))
    ecc_mask = np.sqrt(X**2+Y**2) <= 1.01

    dm = np.zeros((intermezzo_blank_steps+total_steps*n_reps+2*pre_post_blank, n_pix, n_pix), dtype=bool)

    one_cycle = np.zeros((total_steps, n_pix, n_pix))
    x_steps = np.linspace(-1, 1, n_stim_steps, endpoint=True)
    for i, x in enumerate(x_steps):
        one_cycle[i] = (X >= (x-s2)) & (X <= (x+s2)) & ecc_mask

    start_point = pre_post_blank
    for i, d in enumerate(directions):
        for t in np.arange(total_steps):
            dm[start_point+t] = rotate(one_cycle[t], d)
        if i == (len(directions)/2-1):
            print('halfway %i %i'%(d,i))
            start_point += intermezzo_blank_steps
        start_point += total_steps

    return np.round(dm).astype(bool)

def set_pycortex_config_file(project_folder):

    
    # Import necessary modules
    import os
    import cortex
    from pathlib import Path

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
                        hemi):
    """
    Convert pRF fitting value in different parameters for following analysis
   
    Parameters
    ----------
    prf_filename: absolute paths to prf result files.
    output_dir: absolute path to directory into which to put the resulting files.
    stim_radius: stimulus radius in deg
    hemi: brain hemisphere

    Returns
    -------
    prf_deriv_L_all and  prf_deriv_R_all: derivative of pRF analysis for all pRF voxels
    prf_deriv_L_neg and  prf_deriv_R_neg : derivative of pRF analysis for all negative pRF voxels
    prf_deriv_L_pos and  prf_deriv_R_pos : derivative of pRF analysis for all positive pRF voxels

    stucture:
    columns: 1->32492
    row00 : sign
    row01 : R2
    row02 : eccentricity in deg
    row03 : polar angle real component in deg
    row04 : polar angle imaginary component in deg
    row05 : size in deg
    row06 : non-linerity
    row07 : amplitude
    row08 : baseline
    row09 : coverage
    
    ['prf_sign','prf_rsq','prf_ecc','prf_polar_real','prf_polar_imag','prf_size','prf_non_lin','prf_amp','prf_baseline','prf_cov']

    None
    """

    # Imports
    # -------
    # General imports
    import os
    import nibabel as nb
    import glob
    import numpy as np
    import ipdb

    # Popeye imports
    from popeye.spinach import generate_og_receptive_fields

    try:
        os.makedirs(os.path.join(output_dir,'all'))
        os.makedirs(os.path.join(output_dir,'pos'))
        os.makedirs(os.path.join(output_dir,'neg'))
    except:
        pass

    # Get data details
    # ----------------
    prf_data = []
    prf_data_load = nb.load(prf_filename[0])
    prf_data.append(np.array([prf_data_load.darrays[i].data for i in range(len(prf_data_load.darrays))]))
    prf_data = np.vstack(prf_data)    
    ext = prf_data_load.extra
    hdr = prf_data_load.header
    

    # Compute derived measures from prfs
    # ----------------------------------
    # pRF sign
    prf_sign_all = np.sign((prf_data[4,:]))
    pos_mask = prf_sign_all > 0.0
    neg_mask = prf_sign_all < 0.0
    all_mask = pos_mask | neg_mask
    
    # r-square
    prf_rsq_all = prf_data[6,:]

    # pRF eccentricity
    prf_ecc_all = np.nan_to_num(np.sqrt(prf_data[0,:]**2 + prf_data[1,:]**2))

    # pRF polar angle
    complex_polar = prf_data[0,:] + 1j * prf_data[1,:]
    normed_polar = complex_polar / np.abs(complex_polar)
    prf_polar_real_all = np.real(normed_polar)
    prf_polar_imag_all = np.imag(normed_polar)
    
    # pRF size
    prf_size_all = prf_data[2,:].astype(np.float64)
    prf_size_all[prf_size_all<1e-4] = 1e-4

    # pRF non-linearity
    prf_non_lin_all = prf_data[3,:]


    # pRF amplitude
    prf_amp_all = prf_data[4,:]

    # pRF baseline
    prf_baseline_all = prf_data[5,:]


    # pRF coverage
    deg_x, deg_y = np.meshgrid(np.linspace(-30, 30, 121), np.linspace(-30, 30, 121))         # define prfs in visual space
    
    rfs = generate_og_receptive_fields( prf_data[0,:],
                                        prf_data[1,:],
                                        prf_size_all,
                                        np.ones(np.prod(prf_data[0,:].shape[0])),
                                        deg_x,
                                        deg_y)
    
    css_rfs = rfs ** prf_data[3,:]
    total_prf_content = css_rfs.reshape((-1, prf_data.shape[1])).sum(axis=0)
    stim_vignet = np.sqrt(deg_x ** 2 + deg_y**2) < stim_radius    
    prf_cov_all = css_rfs[stim_vignet, :].sum(axis=0) / total_prf_content



    # Saving
    # ------
    for mask_dir in ['all','pos','neg']:
        print('saving: %s'%('os.path.join(output_dir,"{mask_dir}","prf_deriv_{hemi}_{mask_dir}.gii")'.format(hemi = hemi, mask_dir = mask_dir)))
        for output_type in ['prf_sign','prf_rsq','prf_ecc','prf_polar_real','prf_polar_imag','prf_size','prf_non_lin','prf_amp','prf_baseline','prf_cov']:
            exec('{output_type}_{mask_dir} = np.copy({output_type}_all)'.format(mask_dir = mask_dir, output_type = output_type))
            exec('{output_type}_{mask_dir}[~{mask_dir}_mask] = np.nan'.format(mask_dir = mask_dir, output_type = output_type))
        
        exec('prf_deriv_{mask_dir} = np.row_stack((prf_sign_{mask_dir},prf_rsq_{mask_dir},prf_ecc_{mask_dir},prf_polar_real_{mask_dir},\
                prf_polar_imag_all,prf_size_all,prf_non_lin_all,prf_amp_all,prf_baseline_all,prf_cov_all))'.format(mask_dir = mask_dir))
        
        exec('prf_deriv_{mask_dir} = prf_deriv_{mask_dir}.astype(np.float32)'.format(mask_dir = mask_dir))
        exec('darrays = [nb.gifti.gifti.GiftiDataArray(d) for d in prf_deriv_{mask_dir}]'.format(mask_dir = mask_dir))
        exec('gii_out = nb.gifti.gifti.GiftiImage(header = hdr, extra = ext, darrays = darrays)')
        exec('nb.save(gii_out,os.path.join(output_dir,"{mask_dir}","prf_deriv_{hemi}_{mask_dir}.gii"))'.format(hemi = hemi, mask_dir = mask_dir))
            
    return None

def mask_gii_2_hdf5(in_file, mask_file, hdf5_file, folder_alias):
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

    success = True

    
    
     
    gii_in_data = nb.load(in_file)
    data_mat = np.array([gii_in_data.darrays[i].data for i in range(len(gii_in_data.darrays))])
    data_name = op.split(in_file)[-1].split('.gii')[0]
    

    gii_in_mask = nb.load(mask_file)
    mask_mat = np.array([gii_in_mask.darrays[i].data for i in range(len(gii_in_mask.darrays))])
    mask_mat = mask_mat[0,:]
    mask_name = op.split(mask_file)[-1].split('.')[3]

    roi_data = data_mat[:, mask_mat==1]

    try:
        h5file = h5py.File(hdf5_file, "r+")
    except OSError:
        h5file = h5py.File(hdf5_file, "a")
    
    g_hemi = h5file.create_group(folder_alias)
    
    
    dset = g_hemi.create_dataset(data_name,data = roi_data,dtype='float32')


    return None

def draw_cortex_vertex(subject,data,cmap,vmin,vmax,cbar = 'discrete',cmap_steps = 255,\
                        alpha = None,depth = 1,thick = 1,height = 1024,sampler = 'nearest',\
                        with_curvature = True,with_labels = False,with_colorbar = False,\
                        with_borders = False,curv_brightness = 0.95,curv_contrast = 0.05,add_roi = False,\
                        roi_name = 'empty',col_offset = 0):
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
                                    subject = subject)
    
    vertex_rgb_fig = cortex.quickflat.make_figure(  braindata = vertex_rgb,
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
    
    # Color bars
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
        axl.set_yticklabels([0,3,6], size = 'x-large')
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
    
