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

    import ipdb
    ipdb.set_trace()
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
