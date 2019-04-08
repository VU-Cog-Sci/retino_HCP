## retino_HCP

martin_dev_indiv branch
-----------------------

## Analysis specificsd
---------------------
- individual subjects are first preprocessed using FMRIPREP on LISA, see pre_fit/fmriprep_lisa.py
- Well registrered runs were kept for analysis, see select_block.json
- data are sg filter, psc and averaged within subjects, see pre_fit/pre_fit.py
- pRF parameters are extracted using fit/submit_fit_jobs.py on Lisa for gaussian model

- pRF derivatives of subject 999999 are analysed and drawn on fsaverage pycortex map of nprf_hcp using pp_roi.py
- PRF derivatives summary for each ROI of nprf_hcp are put in h5 files with post_pp_roi.py
- roi plots for DMN and Vision ROI are plotted with bokeh using roi_plots.py for subject 999999
- summary statistics per participants are extracted in extract_sum.py
- polar angle progression over ANG using polar_prog.py

# TO DO
# -----
