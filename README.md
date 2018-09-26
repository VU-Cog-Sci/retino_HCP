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
- pRF derivatives of subject 999999 are analysed and drawn on fsaverage pycortex map using pp_roi.py
- Vision regions of interest (V1/V2/V3/VO/DO/LO/SUP_PAR/TPJ/sPCS/iPCS/mPCS/INS/DLPFC) are drawn manually for fsaverage in inkscape
- PRF derivatives summary for each ROI are put in h5 files with post_pp_roi.py
- roi plots for DMN and Vision ROI are plotted with bokeh using roi_plots.py for subject 999999
- pRF derivatives of all others subject are analysed using pp_roi.py and drawn using fsaverage and ROI of 999999
- roi plots for of all others subject for DMN and Vision ROI are plotted with bokeh using roi_plots.py
- summary statistics across participants are extracted in extract_sum.py
- summary statistics per participants are extracted in extract_sum.py

# TO DO
# -----
# 02. make statistics plot across participants
# 04. save webgl and make it html for 999999 subject
# 05. run launch_pp_roi on aeneas
# 06. run launch_post_pp_roi on aeneas
# 07. see if new version of pycortex work on i36 and reply email