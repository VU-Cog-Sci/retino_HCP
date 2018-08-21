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
- Vision regions of interest (V1/V2/V3/LO/TO/DO/SUP_PAR/TPJ/sPCS/iPCS/mPCS/INS/DLPFC) are drawn for fsaverage in inkscape
- PRF derivatives summary for DMN and Vision ROI are plotted with bokeh using post_pp_roi.py for subject 999999
- pRF derivatives of all others subject are analysed using pp_roi.py and drawn using fsaverage and ROI of 999999
- PRF derivatives summary of all others subject for DMN and Vision ROI are plotted with bokeh using post_pp_roi.py

# TO DO
# -----
# 0. understand why output of gauss different for baseline and amplitude??
# 0. make plot of comparison time course and model
# 0. Erase from database the time series of the plot and from the pp_roi and post_pp_roi


# 3. Run averaged subject with gauss and css model
# 4. Draw Visual areas based on polarity/eccentricity flatmap on fsaverage subject: V1/V2/V3/LO/TO/DO/SUP_PAR/TPJ/sPCS/iPCS/mPCS/INS/DLPFC
# 5. Ran post_pp_roi based on all these regions in 10 subjects and look

# note: TPJ is the bottom section of SUP_PAR
# 4. draw ROI in agreement with Tomas
# 5. make launcher of submit_fit_jobs pp_roi and post_pp_roi

