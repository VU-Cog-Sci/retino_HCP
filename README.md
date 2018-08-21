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
# 0. Erase from database the time series of the plot and from the pp_roi and post_pp_roi


# 1. Load Yeo, draw mask on flatmap of fsaverage
# 2. Draw with Inkscape DMN areas: ANG / MED_PAR / SUP_MED_FR / LAT_TEMP
# 3. Run averaged subject with gauss and css model
# 4. Draw Visual areas based on polarity/eccentricity flatmap on fsaverage subject: V1/V2/V3/LO/TO/DO/SUP_PAR/TPJ/sPCS/iPCS/mPCS/INS/DLPFC
# 5. Ran post_pp_roi based on all these regions in 10 subjects and look

# note: TPJ is the bottom section of SUP_PAR
# 4. draw ROI in agreement with Tomas
# 5. make launcher of submit_fit_jobs pp_roi and post_pp_roi

