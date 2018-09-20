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
# 01. make plot of comparison time course and model (understand why output of gauss different for baseline and amplitude)
# 02. statistics across participants - ouput correlations values of ecc vs. all per ROI
# 03. statistics across participants - output laterality index per ROI
# 04. change code to fit 999999 with more steps of bruteforce and run it
# 05. change code to fit 999999_16mm and run it
# 06. change code to draw roi/bokeh plots 999999_16mm
# 08. save webgl and make it html for 999999 subject
# 09. run launch_pp_roi on aeneas
# 10. run launch_post_pp_roi on aeneas
# 11. see if new version of pycortex work on i36 and reply email

# ORDER
# thursday: 09/07/
# friday 04/05