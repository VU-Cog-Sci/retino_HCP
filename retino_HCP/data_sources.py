

def get_data_source_dict(data, condition, colors_val = None, x_gain = 0, y_shift = 0):
    # define main data source
    if (condition == 'map') or (condition == 'ecc') or (condition == 'ecc_gainRatio'):



        data_source                     =   {   
                                            'sign':         data[:,0],                                  # sign
                                            'rsq':          data[:,1],                                  # residual square root
                                            'ecc':          data[:,2],                                  # eccentricity
                                            'sigma':        data[:,5],                                  # size (sigma)
                                            'non_lin':      data[:,6],                                  # non-linearity
                                            'beta':         data[:,7],                                  # amplitude (beta)
                                            'baseline':     data[:,8],                                  # baseline
                                            'cov':          data[:,9]*100,                              # ratio of pRF in stim
                                            'x':            data[:,10],                                 # x coordinate
                                            'y':            data[:,11],                                 # y coordinate
                                            'color':        colors_val}
    elif condition == 'cor':
        dataMat_gazeLeft  = data[0]
        dataMat_gazeRight = data[1]
        dataMat_gazeAll   = data[2]
        data_source                     =   {
                                             'x_all':            dataMat_gazeAll[:,0],                    # x coordinate of gazeAll 
                                             'y_all':            dataMat_gazeAll[:,1],                    # y coordinate of gazeAll 
                                             'sigma_all':        dataMat_gazeAll[:,2],                    # size (sigma) of gazeAll 
                                             'beta_all':         dataMat_gazeAll[:,3],                    # amplitude (beta) of gazeAll 
                                             'baseline_all':     dataMat_gazeAll[:,4],                    # baseline of gazeAll 
                                             'rsq_all':          dataMat_gazeAll[:,5],                    # cross-validated residual square root of gazeAll 
                                             'cv_rsq_all':       dataMat_gazeAll[:,6],                    # cross-validated residual square root of gazeAll 
                                             'stim_ratio_all':   dataMat_gazeAll[:,7]*100,                # ratio of pRF in stim of gazeAll 
                                             'x_left':           dataMat_gazeLeft[:,0],                   # x coordinate of gazeLeft 
                                             'y_left':           dataMat_gazeLeft[:,1],                   # y coordinate of gazeLeft
                                             'sigma_left':       dataMat_gazeLeft[:,2],                   # size (sigma) of gazeLeft
                                             'beta_left':        dataMat_gazeLeft[:,3],                   # amplitude (beta) of gazeLeft
                                             'baseline_left':    dataMat_gazeLeft[:,4],                   # baseline of gazeLeft#                                         
                                             'rsq_left':         dataMat_gazeLeft[:,5],                   # cross-validated residual square root of gazeLeft
                                             'cv_rsq_left':      dataMat_gazeLeft[:,6],                   # cross-validated residual square root of gazeLeft
                                             'stim_ratio_left':  dataMat_gazeLeft[:,7]*100,               # ratio of pRF in stim of gazeLeft
                                             'x_right':          dataMat_gazeRight[:,0],                  # x coordinate of gazeRight 
                                             'y_right':          dataMat_gazeRight[:,1],                  # y coordinate of gazeRight
                                             'sigma_right':      dataMat_gazeRight[:,2],                  # size (sigma) of gazeRight
                                             'beta_right':       dataMat_gazeRight[:,3],                  # amplitude (beta) of gazeRight 
                                             'baseline_right':   dataMat_gazeRight[:,4],                  # baseline of gazeRight
                                             'rsq_right':        dataMat_gazeRight[:,5],                  # cross-validated residual square root of gazeRight
                                             'cv_rsq_right':     dataMat_gazeRight[:,6],                  # cross-validated residual square root of gazeRight
                                             'stim_ratio_right': dataMat_gazeRight[:,7]*100,              # ratio of pRF in stim of gazeRight
                                             'retinal_x_gain':   dataMat_gazeAll[:,-5],
                                             'screen_x_gain':    dataMat_gazeAll[:,-4],
                                             'retinal_gain_index':dataMat_gazeAll[:,-3],
                                             'amplitude_change': dataMat_gazeAll[:,-2],
                                             'y_change':         dataMat_gazeAll[:,-1]
                                            }
    elif condition == 'shift':
        # extract data from the data_params list
        dataMat_gazeLeft  = data[0]
        dataMat_gazeRight = data[1]
        dataMat_gazeAll   = data[2]
        # define data source
        data_source                     =   {
                                             'x_left':              dataMat_gazeLeft[:,0],                # x coordinate of gazeLeft 
                                             'y_left':              dataMat_gazeLeft[:,1],                # y coordinate of gazeLeft
                                             'sigma_left':          dataMat_gazeLeft[:,2],                # size (sigma) of gazeLeft
                                             'beta_left':           dataMat_gazeLeft[:,3],                # amplitude (beta) of gazeLeft
                                             'baseline_left':       dataMat_gazeLeft[:,4],                # baseline of gazeLeft#                                         
                                             'cv_rsq_left':         dataMat_gazeLeft[:,6],                # cross-validated residual square root of gazeLeft
                                             'stim_ratio_left':     dataMat_gazeLeft[:,7]*100,            # ratio of pRF in stim of gazeLeft
                                             'x_right':             dataMat_gazeRight[:,0],               # x coordinate of gazeRight 
                                             'y_right':             dataMat_gazeRight[:,1],               # y coordinate of gazeRight
                                             'sigma_right':         dataMat_gazeRight[:,2],               # size (sigma) of gazeRight
                                             'beta_right':          dataMat_gazeRight[:,3],               # amplitude (beta) of gazeRight 
                                             'baseline_right':      dataMat_gazeRight[:,4],               # baseline of gazeRight
                                             'cv_rsq_right':        dataMat_gazeRight[:,6],               # cross-validated residual square root of gazeRight
                                             'stim_ratio_right':    dataMat_gazeRight[:,7]*100,           # ratio of pRF in stim of gazeRight
                                             'retinal_x_gain':      dataMat_gazeAll[:,-5],
                                             'screen_x_gain':       dataMat_gazeAll[:,-4],
                                             'retinal_gain_index':  dataMat_gazeAll[:,-3],
                                             'amplitude_change':    dataMat_gazeAll[:,-2],
                                             'y_change':            dataMat_gazeAll[:,-1]
                                            }

    elif condition == 'shift_spatiotopic':
        # extract data from the data_params list
        dataMat_gazeLeft  = data[0]
        dataMat_gazeRight = data[1]
        # define data source
        data_source                     =   {
                                             'x_left':              dataMat_gazeLeft[:,0],                # x coordinate of gazeLeft 
                                             'y_left':              dataMat_gazeLeft[:,1],                # y coordinate of gazeLeft
                                             'sigma_left':          dataMat_gazeLeft[:,2],                # size (sigma) of gazeLeft
                                             'beta_left':           dataMat_gazeLeft[:,3],                # amplitude (beta) of gazeLeft
                                             'baseline_left':       dataMat_gazeLeft[:,4],                # baseline of gazeLeft#                                         
                                             'cv_rsq_left':         dataMat_gazeLeft[:,6],                # cross-validated residual square root of gazeLeft
                                             'stim_ratio_left':     dataMat_gazeLeft[:,7]*100,            # ratio of pRF in stim of gazeLeft
                                             'x_right':             dataMat_gazeRight[:,0],               # x coordinate of gazeRight 
                                             'y_right':             dataMat_gazeRight[:,1],               # y coordinate of gazeRight
                                             'sigma_right':         dataMat_gazeRight[:,2],               # size (sigma) of gazeRight
                                             'beta_right':          dataMat_gazeRight[:,3],               # amplitude (beta) of gazeRight 
                                             'baseline_right':      dataMat_gazeRight[:,4],               # baseline of gazeRight
                                             'cv_rsq_right':        dataMat_gazeRight[:,6],               # cross-validated residual square root of gazeRight
                                             'stim_ratio_right':    dataMat_gazeRight[:,7]*100,           # ratio of pRF in stim of gazeRight
                                             'retinal_x_gain':      dataMat_gazeAll[:,-5],
                                             'screen_x_gain':       dataMat_gazeAll[:,-4],
                                             'retinal_gain_index':  dataMat_gazeAll[:,-3],
                                             'amplitude_change':    dataMat_gazeAll[:,-2],
                                             'y_change':            dataMat_gazeAll[:,-1]
                                            }
    elif condition == 'gain':
        dataMat_gazeLeft  = data[0]
        dataMat_gazeRight = data[1]
        dataMat_gazeAll   = data[2]
        data_source                     =   {
                                             'x_all':            dataMat_gazeAll[:,0],                    # x coordinate of gazeAll 
                                             'y_all':            dataMat_gazeAll[:,1],                    # y coordinate of gazeAll 
                                             'sigma_all':        dataMat_gazeAll[:,2],                    # size (sigma) of gazeAll 
                                             'beta_all':         dataMat_gazeAll[:,3],                    # amplitude (beta) of gazeAll 
                                             'baseline_all':     dataMat_gazeAll[:,4],                    # baseline of gazeAll 
                                             'cv_rsq_all':       dataMat_gazeAll[:,6],                    # cross-validated residual square root of gazeAll 
                                             'stim_ratio_all':   dataMat_gazeAll[:,7]*100,                # ratio of pRF in stim of gazeAll 
                                             'ecc_all':          dataMat_gazeAll[:,8],                    # eccentricity of gazeAll
                                             'x_left':           dataMat_gazeLeft[:,0],                   # x coordinate of gazeLeft 
                                             'y_left':           dataMat_gazeLeft[:,1],                   # y coordinate of gazeLeft
                                             'sigma_left':       dataMat_gazeLeft[:,2],                   # size (sigma) of gazeLeft
                                             'beta_left':        dataMat_gazeLeft[:,3],                   # amplitude (beta) of gazeLeft
                                             'baseline_left':    dataMat_gazeLeft[:,4],                   # baseline of gazeLeft#                                         
                                             'cv_rsq_left':      dataMat_gazeLeft[:,6],                   # cross-validated residual square root of gazeLeft
                                             'stim_ratio_left':  dataMat_gazeLeft[:,7]*100,               # ratio of pRF in stim of gazeLeft
                                             'ecc_left':         dataMat_gazeLeft[:,8],                   # eccentricity of gazeLeft
                                             'x_right':          dataMat_gazeRight[:,0],                  # x coordinate of gazeRight 
                                             'y_right':          dataMat_gazeRight[:,1],                  # y coordinate of gazeRight
                                             'sigma_right':      dataMat_gazeRight[:,2],                  # size (sigma) of gazeRight
                                             'beta_right':       dataMat_gazeRight[:,3],                  # amplitude (beta) of gazeRight 
                                             'baseline_right':   dataMat_gazeRight[:,4],                  # baseline of gazeRight
                                             'cv_rsq_right':     dataMat_gazeRight[:,6],                  # cross-validated residual square root of gazeRight
                                             'stim_ratio_right': dataMat_gazeRight[:,7]*100,              # ratio of pRF in stim of gazeRight
                                             'ecc_right':        dataMat_gazeLeft[:,8],                   # eccentricity of gazeRight
                                             'x_gain':           x_gain,                                            # gaze gain
                                             'y_shift':          y_shift,                                           # y shift
                                             'gain_color':       colors_val,                                        # colors based on gain
                                             'retinal_x_gain':   dataMat_gazeAll[:,-5],
                                             'screen_x_gain':    dataMat_gazeAll[:,-4],
                                             'retinal_gain_index':dataMat_gazeAll[:,-3],
                                             'amplitude_change': dataMat_gazeAll[:,-2],
                                             'y_change':         dataMat_gazeAll[:,-1]
                                        }

    elif condition == 'roi':
        dataMat_gazeAll_avg     =   data[0]                          # data for gaze all averaged
        dataMat_gazeAll_std     =   data[1]                          # data for gaze all std
        dataMat_gazeAll_num     =   data[2]                          # data for gaze all num of voxel
        dataMat_gazeLeft_avg    =   data[3]                         # data for gaze left averaged
        dataMat_gazeLeft_std    =   data[4]                         # data for gaze left std
        dataMat_gazeLeft_num    =   data[5]                         # data for gaze left num of voxel
        dataMat_gazeRight_avg   =   data[6]                        # data for gaze right averaged
        dataMat_gazeRight_std   =   data[7]                        # data for gaze right std
        dataMat_gazeRight_num   =   data[8]                        # data for gaze right num of voxel
        dataMat_gazeGain_avg    =   data[9]                         # data for gaze gain avg
        dataMat_gazeGain_std    =   data[10]
        rois                    =   data[11]      
        dataMat_gazeAll         =   data[-1]
        # import ipdb ; ipdb.set_trace()

        data_source                     =   {
                                             'rois':                rois,                                 # rois
                                             # 'ecc_all':             self.dataMat_gazeAll_avg[:,8],             # avg eccentricity of gazeAll
                                             # 'sigma_all':           self.dataMat_gazeAll_avg[:,2],             # avg size (sigma) of gazeAll 
                                             # 'cv_rsq_all':          self.dataMat_gazeAll_avg[:,6],             # avg cross-validated residual square root of gazeAll 
                                             # 'beta_all':            self.dataMat_gazeAll_avg[:,3],             # avg amplitude (beta) of gazeAll 
                                             # 'stim_ratio_all':      self.dataMat_gazeAll_avg[:,7]*100,         # avg ratio of pRF in stim of gazeAll 
                                             # 'ecc_all_std':         self.dataMat_gazeAll_std[:,8],             # std eccentricity of gazeAll
                                             # 'sigma_all_std':       self.dataMat_gazeAll_std[:,2],             # std size (sigma) of gazeAll 
                                             # 'cv_rsq_all_std':      self.dataMat_gazeAll_std[:,6],             # std cross-validated residual square root of gazeAll 
                                             # 'beta_all_std':        self.dataMat_gazeAll_std[:,3],             # std amplitude (beta) of gazeAll 
                                             # 'stim_ratio_all_std':  self.dataMat_gazeAll_std[:,7]*100,         # std ratio of pRF in stim of gazeAll 
                                             # 'voxel_all':           self.dataMat_gazeAll_num,                  # number of voxels of gazeAll
                                             'ecc_left':            dataMat_gazeLeft_avg[:,8],            # avg eccentricity of gazeLeft
                                             'sigma_left':          dataMat_gazeLeft_avg[:,2],            # avg size (sigma) of gazeLeft 
                                             'cv_rsq_left':         dataMat_gazeLeft_avg[:,6],            # avg cross-validated residual square root of gazeLeft 
                                             'beta_left':           dataMat_gazeLeft_avg[:,3],            # avg amplitude (beta) of gazeLeft 
                                             'stim_ratio_left':     dataMat_gazeLeft_avg[:,7]*100,        # avg ratio of pRF in stim of gazeLeft 
                                             'ecc_left_std':        dataMat_gazeLeft_std[:,8],            # std eccentricity of gazeLeft
                                             'sigma_left_std':      dataMat_gazeLeft_std[:,2],            # std size (sigma) of gazeLeft 
                                             'cv_rsq_left_std':     dataMat_gazeLeft_std[:,6],            # std cross-validated residual square root of gazeLeft 
                                             'beta_left_std':       dataMat_gazeLeft_std[:,3],            # std amplitude (beta) of gazeLeft 
                                             'stim_ratio_left_std': dataMat_gazeLeft_std[:,7]*100,        # std ratio of pRF in stim of gazeLeft 
                                             'voxel_left':          dataMat_gazeLeft_num,                 # number of voxels of gazeLeft 
                                             'ecc_right':           dataMat_gazeRight_avg[:,8],           # avg eccentricity of gazeRight
                                             'sigma_right':         dataMat_gazeRight_avg[:,2],           # avg size (sigma) of gazeRight 
                                             'cv_rsq_right':        dataMat_gazeRight_avg[:,6],           # avg cross-validated residual square root of gazeRight 
                                             'beta_right':          dataMat_gazeRight_avg[:,3],           # avg amplitude (beta) of gazeRight 
                                             'stim_ratio_right':    dataMat_gazeRight_avg[:,7]*100,       # avg ratio of pRF in stim of gazeRight 
                                             'ecc_right_std':       dataMat_gazeRight_std[:,8],           # std eccentricity of gazeRight
                                             'sigma_right_std':     dataMat_gazeRight_std[:,2],           # std size (sigma) of gazeRight 
                                             'cv_rsq_right_std':    dataMat_gazeRight_std[:,6],           # std cross-validated residual square root of gazeRight 
                                             'beta_right_std':      dataMat_gazeRight_std[:,3],           # std amplitude (beta) of gazeRight 
                                             'stim_ratio_right_std':dataMat_gazeRight_std[:,7]*100,       # std ratio of pRF in stim of gazeRight 
                                             'voxel_right':         dataMat_gazeRight_num,                # number of voxels of gazeLeft 
                                             'gaze_gain':           dataMat_gazeGain_avg,                 # avg ratio of gaze gain  
                                             'gaze_gain_std':       dataMat_gazeGain_std,                 # std ratio of gaze gain  
                                             'retinal_x_gain':      dataMat_gazeAll_avg[:,-5],
                                             'screen_x_gain':       dataMat_gazeAll_avg[:,-4],
                                             'retinal_gain_index':  dataMat_gazeAll_avg[:,-3],
                                             'amplitude_change':    dataMat_gazeAll_avg[:,-2],
                                             'y_change':            dataMat_gazeAll_avg[:,-1],
                                             'retinal_x_gain_std':  dataMat_gazeAll_std[:,-5],
                                             'screen_x_gain_std':   dataMat_gazeAll_std[:,-4],
                                             'retinal_gain_index_std':dataMat_gazeAll_std[:,-3],
                                             'amplitude_change_std':dataMat_gazeAll_std[:,-2],
                                             'y_change_std':        dataMat_gazeAll_std[:,-1]
                                              }
    return data_source
