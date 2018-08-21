## retino_HCP

martin_dev branch
-----------------

## Analysis specifics
- HCP subjects were first pre-processed and averaged by Tomas
- Best subjects (based on R2 of V1 voxel selected by glasser2016 attlas) were chosen using data
  issued from Benson et al (http://dx.doi.org/10.1101/308247)
- pRF parameters are extracted using fit/submit_fit_jobs.py on Lisa for css and gaussian model
- pRF derivatives are analysed using pp_roi.py
- pRF derivatives data are resampled from 32k .gii to fsaverage using pp_roi.py
- pRF derivatives are plotted on flatmap with pycortex using pp_roi.py
- ROI are drawn in inkscape following the same procedure as in nPRF
- ROI masks are created from inkscape drawings and converted from fsaverage to 32k .gii
- pRF derivatives are saved as a function of 32k gii ROI mask as hdf5 files using post_pp_roi.py
- PRF derivatives summary are plotted with bokeh using post_pp_roi.py


# TO DO
# -----
# 0. understand why output of gauss different for baseline and amplitude??
# 0. make plot of comparison time course and model
# 1. make plot of comparison time course and model
# 2. Erase from database the time series of the plot and from the pp_roi and post_pp_roi

# 4. draw ROI in agreement with Tomas
# 5. make launcher of submit_fit_jobs pp_roi and post_pp_roi

