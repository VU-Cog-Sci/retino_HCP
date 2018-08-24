## retino_HCP

martin_dev branch
-----------------

## Analysis specifics
---------------------
- HCP subjects were first pre-processed and averaged by Tomas
- Best subjects (based on R2 of V1 voxel selected by glasser2016 attlas) were chosen using data
  issued from Benson et al (http://dx.doi.org/10.1101/308247)
- pRF parameters are extracted using fit/submit_fit_jobs.py on Lisa for css and gaussian model
- get yeo atlas see others/get_yeo_atlas_gii.txt
- DMN regions (ANG/MED_PAR/SUP_MED_FR/LAT_TEMP) are put on fsaverage flatmap overlay.svg
- DMN regions of interest are drawn in inkscape of fsaverage
- pRF derivatives of subject 999999 are analysed and drawn on fsaverage using pp_roi.py
- Vision regions of interest (V1/V2/V3/VO/DO/LO/SUP_PAR/TPJ/sPCS/iPCS/mPCS/INS/DLPFC) are drawn for fsaverage in inkscape
- PRF derivatives summary for DMN and Vision ROI are plotted with bokeh using post_pp_roi.py for subject 999999
- pRF derivatives of all others subject are analysed using pp_roi.py and drawn using fsaverage and ROI of 999999
- PRF derivatives summary of all others subject for DMN and Vision ROI are plotted with bokeh using post_pp_roi.py

# TO DO
# -----
# 0. understand why output of gauss different for baseline and amplitude??
# 0. make plot of comparison time course and model

# 2. Ran post_pp_roi based on all these regions for suject 999999 => debug it (check how masks work before saving in h5file, how is it possible that with 0 vertex the code ran and bugged)



# Everyday
# - transfer data from lisa to aeneas
# - run launcher of pp_roi
# - run launcher of post_pp_roi



# 6. Check variability of individual subjects with 999999 ROIs
# 7. Decide what to draw as statistics across participants
