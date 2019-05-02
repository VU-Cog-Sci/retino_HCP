## retino_HCP

martin_dev_subcortical branch
-----------------------------

## Analysis specifics
---------------------
- HCP subcortical subjects were first pre-processed and averaged by Tomas
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


- pRF derivatives of each subjects are analysed and drawn on fsaverage pycortex map of nprf_hcp using pp_roi.py
- each ROI masks are and saved in hdf5 files using post_pp_roi.py
- roi plots for Vision and DMM rois are plotted with bokeh using roi_plots.py
- Compute summary statistics per participants using extract_sum.py
- Draw summary statistics per participants using summary_plots.py

- polar angle progression over ANG using polar_prog.py


# TO DO
# -----
# Run it on LISA to see if it works with 10 steps
# Compare with results of gauss obtain with martin_dev_indiv
