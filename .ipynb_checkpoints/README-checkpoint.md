## retino_HCP

martin_dev_indiv branch
-----------------------

## Analysis specificsd
---------------------
- individual subjects are first preprocessed using FMRIPREP on LISA, see pre_fit/fmriprep_lisa.py
- Well registrered runs were kept for analysis, see select_block.json
- data are sg filter, psc and averaged within subjects, see pre_fit/pre_fit.py
- pRF parameters are extracted using fit/submit_fit_jobs.py on Lisa for gaussian model

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
- polar angle progression over ANG using polar_prog.py

# TO DO
# -----
# 01. check summary plot use average weighted and correct error bars and see what is going on without threhsold
# 02. make error bar for subject 999999
# 04. change ecc plots to not have baseline and amplitude incorporated for all subjects
# 05. change map plot to have crosses instead of circles of different size
# 06. save individual viewers with colors changed directly to include the r2 in colors not transparency
# 07. make average of subject fit and plot every data of '000000'