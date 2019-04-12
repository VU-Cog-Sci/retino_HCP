# General imports
from __future__ import division
import sys
import ctypes
import numpy as np
import scipy.io
import platform
from math import *
import os
import glob
import json
opj = os.path.join
import ipdb
deb = ipdb.set_trace
import warnings
warnings.filterwarnings('ignore')

# MRI analysis imports
import nibabel as nb
import popeye.utilities as utils
from popeye.visual_stimulus import VisualStimulus
import popeye.css as css
import popeye.og as og
import cifti
from joblib import Parallel, delayed
from scipy import signal

import h5py
from scipy.signal import fftconvolve, savgol_filter
from popeye.base import PopulationModel
from popeye.spinach import generate_og_receptive_field, generate_rf_timeseries_nomask

def fit_gradient_descent(model, data, ballpark, bounds, verbose=0):
    return utils.gradient_descent_search(data,
                                         utils.error_function,
                                         model.generate_prediction,
                                         ballpark,
                                         bounds,
                                         verbose)[0]

class CompressiveSpatialSummationModel(PopulationModel):
    
    r"""
    A Compressive Spatial Summation population receptive field model class    
    """
    
    def __init__(self, stimulus, hrf_model, cached_model_path=None, nuisance=None):
                
        PopulationModel.__init__(self, stimulus, hrf_model, cached_model_path, nuisance)

    # main method for deriving model time-series
    def generate_prediction(self, x, y, sigma, n, beta, baseline):
        
        # generate the RF
        rf = generate_og_receptive_field(
            x, y, sigma, self.stimulus.deg_x, self.stimulus.deg_y)
        
        # normalize by the integral
        rf /= ((2 * np.pi * sigma**2) * 1 /
               np.diff(self.stimulus.deg_x[0, 0:2])**2)
        
        # extract the stimulus time-series
        response = generate_rf_timeseries_nomask(self.stimulus.stim_arr, rf)
        
        # compression
        response **= n
        
        # convolve with the HRF
        hrf = self.hrf_model(self.hrf_delay, self.stimulus.tr_length)
        
        # convolve it with the stimulus
        model = fftconvolve(response, hrf)[0:len(response)]
        
        # units
        model /= np.max(model)
        
        # offset
        model += baseline
        
        # scale it by beta
        model *= beta

        return model

class CompressiveSpatialSummationModelFiltered(PopulationModel):
    
    r"""
    A Compressive Spatial Summation population receptive field model class
    Adapted to include a savgol_filter
    
    """
    
    def __init__(self, stimulus, hrf_model, cached_model_path=None, nuisance=None, sg_filter_window_length=120, sg_filter_polyorder=3, sg_filter_deriv = 0, tr=1.5):
                
        PopulationModel.__init__(self, stimulus, hrf_model, cached_model_path, nuisance)

        # sg filter
        self.sg_filter_window = np.int(sg_filter_window_length / tr)
        if self.sg_filter_window % 2 == 0:
            self.sg_filter_window += 1
        self.sg_filter_polyorder = sg_filter_polyorder
        self.sg_filter_deriv = sg_filter_deriv

    # main method for deriving model time-series
    def generate_prediction(self, x, y, sigma, n, beta, baseline):
        
        # generate the RF
        rf = generate_og_receptive_field(
            x, y, sigma, self.stimulus.deg_x, self.stimulus.deg_y)
        
        # normalize by the integral
        rf /= ((2 * np.pi * sigma**2) * 1 /
               np.diff(self.stimulus.deg_x[0, 0:2])**2)
        
        # extract the stimulus time-series
        response = generate_rf_timeseries_nomask(self.stimulus.stim_arr, rf)
        
        # compression
        response **= n
        
        # convolve with the HRF
        hrf = self.hrf_model(self.hrf_delay, self.stimulus.tr_length)
        
        # convolve it with the stimulus
        model = fftconvolve(response, hrf)[0:len(response)]
        
        # units
        model /= np.max(model)
        
        # at this point, add filtering with a savitzky-golay filter
        model_drift = savgol_filter(model, 
                                    window_length = self.sg_filter_window, 
                                    polyorder = self.sg_filter_polyorder,
                                    deriv = self.sg_filter_deriv, 
                                    mode = 'nearest')
        
        # demain model_drift, so baseline parameter is still interpretable
        model_drift_demeaned = model_drift-np.mean(model_drift)
        
        # and apply to data
        model -= model_drift_demeaned    
        
        # offset
        model += baseline
        
        # scale it by beta
        model *= beta

        return model
    
class GaussianModel(PopulationModel):
    
    r"""
    A Gaussian Spatial Summation population receptive field model class
    
    """
    
    def __init__(self, stimulus, hrf_model, cached_model_path=None, nuisance=None):
                
        PopulationModel.__init__(self, stimulus, hrf_model, cached_model_path, nuisance)

    # main method for deriving model time-series
    def generate_prediction(self, x, y, sigma, beta, baseline):
        
        # generate the RF
        rf = generate_og_receptive_field(
            x, y, sigma, self.stimulus.deg_x, self.stimulus.deg_y)
        
        # normalize by the integral
        rf /= ((2 * np.pi * sigma**2) * 1 /
               np.diff(self.stimulus.deg_x[0, 0:2])**2)
        
        # extract the stimulus time-series
        response = generate_rf_timeseries_nomask(self.stimulus.stim_arr, rf)
                
        # convolve with the HRF
        hrf = self.hrf_model(self.hrf_delay, self.stimulus.tr_length)
        
        # convolve it with the stimulus
        model = fftconvolve(response, hrf)[0:len(response)]
        
        # units
        model /= np.max(model)
        
        # offset
        model += baseline
        
        # scale it by beta
        model *= beta

        return model
    
class GaussianModelFiltered(PopulationModel):
    
    r"""
    A Gaussian Spatial Summation population receptive field model class
    Adapted to include a savgol_filter
    
    """
    
    def __init__(self, stimulus, hrf_model, cached_model_path=None, nuisance=None, sg_filter_window_length=120, sg_filter_polyorder=3, sg_filter_deriv = 0, tr=1.5):
                
        PopulationModel.__init__(self, stimulus, hrf_model, cached_model_path, nuisance)

        # sg filter
        self.sg_filter_window = np.int(sg_filter_window_length / tr)
        if self.sg_filter_window % 2 == 0:
            self.sg_filter_window += 1
        self.sg_filter_polyorder = sg_filter_polyorder
        self.sg_filter_deriv = sg_filter_deriv

    # main method for deriving model time-series
    def generate_prediction(self, x, y, sigma, beta, baseline):
        
        # generate the RF
        rf = generate_og_receptive_field(
            x, y, sigma, self.stimulus.deg_x, self.stimulus.deg_y)
        
        # normalize by the integral
        rf /= ((2 * np.pi * sigma**2) * 1 /
               np.diff(self.stimulus.deg_x[0, 0:2])**2)
        
        # extract the stimulus time-series
        response = generate_rf_timeseries_nomask(self.stimulus.stim_arr, rf)
                
        # convolve with the HRF
        hrf = self.hrf_model(self.hrf_delay, self.stimulus.tr_length)
        
        # convolve it with the stimulus
        model = fftconvolve(response, hrf)[0:len(response)]
        
        # units
        model /= np.max(model)
        
        # at this point, add filtering with a savitzky-golay filter
        model_drift = savgol_filter(model, 
                                    window_length = self.sg_filter_window, 
                                    polyorder = self.sg_filter_polyorder,
                                    deriv = self.sg_filter_deriv, 
                                    mode = 'nearest')
        
        # demain model_drift, so baseline parameter is still interpretable
        model_drift_demeaned = model_drift-np.mean(model_drift)
        
        # and apply to data
        model -= model_drift_demeaned    
        
        # offset
        model += baseline
        
        # scale it by beta
        model *= beta

        return model
    
class prf_fit(object):
    
    def __init__(self, fit_model, visual_design, screen_distance, screen_width, 
                 tr, bound_grids, grid_steps, bound_fits,
                 sg_filter_window_length = 210,sg_filter_polyorder = 3,sg_filter_deriv =0):

        self.stimulus = VisualStimulus( stim_arr = visual_design,
                                        viewing_distance = screen_distance, 
                                        screen_width =screen_width,
                                        scale_factor = 1,
                                        tr_length = tr,
                                        dtype = np.short)
        if fit_model == 'gauss':
            self.model_func = GaussianModel(stimulus = self.stimulus, 
                                                    hrf_model = utils.spm_hrf)    
        elif fit_model == 'gauss_sg':
            self.model_func = GaussianModelFiltered(stimulus = self.stimulus,
                                            hrf_model = utils.spm_hrf,
                                            sg_filter_window_length = sg_filter_window_length, 
                                            sg_filter_polyorder = sg_filter_polyorder,
                                            sg_filter_deriv = sg_filter_deriv,
                                            tr = tr)

        elif fit_model == 'css':
            self.model_func = CompressiveSpatialSummationModel( stimulus = self.stimulus,
                                                                    hrf_model = utils.spm_hrf)
        elif fit_model == 'css_sg':
            self.model_func = CompressiveSpatialSummationModelFiltered( stimulus = self.stimulus,
                                                                    hrf_model = utils.spm_hrf,
                                                                    sg_filter_window_length = sg_filter_window_length, 
                                                                    sg_filter_polyorder = sg_filter_polyorder,
                                                                    sg_filter_deriv = sg_filter_deriv,
                                                                    tr = tr)
        self.model_func.hrf_delay = 0
        self.predictions = None      
        self.fit_model =  fit_model
        self.bound_grids = bound_grids
        self.grid_steps = grid_steps
        self.bound_fits = bound_fits
        
    def make_grid(self,save_file):

        prf_xs = np.linspace(self.bound_grids[0][0],self.bound_grids[0][1],self.grid_steps)
        prf_ys = np.linspace(self.bound_grids[1][0],self.bound_grids[1][1],self.grid_steps)
        prf_sigma = np.linspace(self.bound_grids[2][0],self.bound_grids[2][1],self.grid_steps)
        
        if self.fit_model == 'gauss' or self.fit_model == 'gauss_sg':
            self.prf_xs, self.prf_ys, self.prf_sigma = np.meshgrid(prf_xs, prf_ys, prf_sigma)
                
        elif self.fit_model == 'css' or self.fit_model == 'css_sg':
            prf_n = np.linspace(self.bound_grids[3][0],self.bound_grids[3][1],self.grid_steps)
            self.prf_xs, self.prf_ys, self.prf_sigma, self.prf_n = np.meshgrid(prf_xs, prf_ys, prf_sigma, prf_n)
                
                
        self.predictions = np.zeros(list(self.prf_xs.shape) + [self.stimulus.run_length])
        self.predictions = self.predictions.reshape(-1, self.predictions.shape[-1]).T
        
        if self.fit_model == 'gauss' or self.fit_model == 'gauss_sg':
            for i, (x, y, s) in enumerate(zip(self.prf_xs.ravel(), self.prf_ys.ravel(), self.prf_sigma.ravel())):
                self.predictions[:, i] = self.model_func.generate_prediction(x, y, s, 1, 0)
        elif self.fit_model == 'css' or self.fit_model == 'css_sg':
            for i, (x, y, s, n) in enumerate(zip(self.prf_xs.ravel(), self.prf_ys.ravel(), self.prf_sigma.ravel(), self.prf_n.ravel())):
                self.predictions[:, i] = self.model_func.generate_prediction(x, y, s, n, 1, 0)
                
                
        # save
        if self.fit_model == 'gauss' or self.fit_model == 'gauss_sg':
            with h5py.File(save_file, 'a') as f:
                f.create_dataset('bound_grids', data = self.bound_grids)
                f.create_dataset('grid_steps', data = self.grid_steps)
                f.create_dataset('bound_fits', data = self.bound_fits)
                f.create_dataset('predictions', data = self.predictions)
                f.create_dataset('prf_xs', data = self.prf_xs)
                f.create_dataset('prf_ys', data = self.prf_ys)
                f.create_dataset('prf_sigma', data = self.prf_sigma)
                
        elif self.fit_model == 'css' or self.fit_model == 'css_sg':
            with h5py.File(save_file, 'a') as f:
                f.create_dataset('bound_grids', data = self.bound_grids)
                f.create_dataset('grid_steps', data = self.grid_steps)
                f.create_dataset('bound_fits', data = self.bound_fits)
                f.create_dataset('predictions', data = self.predictions)
                f.create_dataset('prf_xs', data = self.prf_xs)
                f.create_dataset('prf_ys', data = self.prf_ys)
                f.create_dataset('prf_sigma', data = self.prf_sigma)
                f.create_dataset('prf_n', data= self.prf_n)            
    
    def load_grid(self,save_file):
        grid_predictions = h5py.File(save_file, 'r')
        
        self.prf_xs = grid_predictions['prf_xs'][:]
        self.prf_ys = grid_predictions['prf_ys'][:]
        self.prf_sigma = grid_predictions['prf_sigma'][:]
        if self.fit_model == 'css' or self.fit_model == 'css_sg':
            self.prf_n = grid_predictions['prf_n'][:]  
        self.predictions = grid_predictions['predictions'][:]
    
    def fit_grid(self,data):
        
        self.data = data
        
        if self.fit_model == 'gauss' or self.fit_model == 'gauss_sg':
            fit_grid_params = np.ones((self.data.shape[1], 6))*np.nan
        elif self.fit_model == 'css' or self.fit_model == 'css_sg':
            fit_grid_params = np.ones((self.data.shape[1], 7))*np.nan
        
        for vox_num in range (0,self.data.shape[1]):
            data_ = self.data[:,vox_num]

            if self.fit_model == 'gauss' or self.fit_model == 'gauss_sg':
                prediction_params = np.ones((self.predictions.shape[1], 6))*np.nan
            elif self.fit_model == 'css' or self.fit_model == 'css_sg':
                prediction_params = np.ones((self.predictions.shape[1], 7))*np.nan

            for prediction_num in range(0,self.predictions.shape[1]):
                prediction_ = self.predictions[:,prediction_num]
                betas, residual, _, _ = np.linalg.lstsq((np.vstack([np.ones_like(prediction_),np.nan_to_num(prediction_)])).T, np.nan_to_num(data_).T)
                baseline = betas[0]
                beta = betas[1]
                if beta == 0:
                    rsqs = 0.0
                else:
                    rsqs = ((1 - residual / (data_.shape[0] * data_.var(axis=-1))))[0]

                if self.fit_model == 'gauss' or self.fit_model == 'gauss_sg':
                    prediction_params[prediction_num,:] = [ self.prf_xs.ravel()[prediction_num],
                                                            self.prf_ys.ravel()[prediction_num],
                                                            self.prf_sigma.ravel()[prediction_num],
                                                            beta,
                                                            baseline,
                                                            rsqs
                                                           ]
                elif self.fit_model == 'css' or self.fit_model == 'css_sg':
                    prediction_params[prediction_num,:] = [ self.prf_xs.ravel()[prediction_num],
                                                            self.prf_ys.ravel()[prediction_num],
                                                            self.prf_sigma.ravel()[prediction_num],
                                                            self.prf_n.ravel()[prediction_num],
                                                            beta,
                                                            baseline,
                                                            rsqs
                                                           ]
            best_prediction_num = np.argmax(prediction_params[:,-1], 0)
            fit_grid_params[vox_num,:] = prediction_params[best_prediction_num,:]

        self.gridsearch_params = fit_grid_params[:,:-1]
        self.gridsearch_r2 = fit_grid_params[:,-1]
        self.gridsearch_all = fit_grid_params
                
    def fit_prf(self, data, n_jobs, gridsearch_params):
        
        self.data = data
        self.n_jobs = n_jobs
        self.gridsearch_params = gridsearch_params
        if self.gridsearch_params is None:
            raise Exception('First use self.fit_grid!')
        
        
        prf_params = Parallel(self.n_jobs,verbose = 10,prefer='processes')(delayed(fit_gradient_descent)(self.model_func, data, ballpark, self.bound_fits)
                                       for (data,ballpark) in zip(self.data.T, self.gridsearch_params))
        prf_params = np.vstack(prf_params)
        if self.fit_model == 'gauss' or self.fit_model == 'gauss_sg':
            output = np.ones((self.data.shape[1],6))*np.nan
        elif self.fit_model == 'css' or self.fit_model == 'css_sg':
            output = np.ones((self.data.shape[1],7))*np.nan
            
        for vox in range(0,self.data.shape[1]):
            data_tc = self.data[:,vox]
            if self.fit_model == 'gauss' or self.fit_model == 'gauss_sg':
                model_tc = self.model_func.generate_prediction(prf_params[vox,0],prf_params[vox,1],prf_params[vox,2],prf_params[vox,3],prf_params[vox,4])
            elif self.fit_model == 'css' or self.fit_model == 'css_sg':
                model_tc = self.model_func.generate_prediction(prf_params[vox,0],prf_params[vox,1],prf_params[vox,2],prf_params[vox,3],prf_params[vox,4],prf_params[vox,5])
                
            output[vox,:] = np.hstack([prf_params[vox,:], utils.coeff_of_determination(data_tc,model_tc)/100.0])
        
        self.fit_output = output