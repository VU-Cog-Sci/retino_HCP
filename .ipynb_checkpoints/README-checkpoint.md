## retino_HCP

martin_dev_indiv branch
-----------------------
Branch analysing data of 2 individual subjects

## Analysis specifics
---------------------
- individual subjects are first preprocessed using FMRIPREP on LISA, see pre_fit/fmriprep_lisa.py
- Well registrered runs were kept for analysis, see select_block.json
- data are sg filter, psc and averaged within subjects, see pre_fit/pre_fit.py
- use pycortex_webgl to view the average time course on a pycortex server
- pRF parameters are extracted using fit/submit_fit_jobs.py on Lisa for gaussian model
- pRF derivatives of each subjects are analysed and drawn on fsaverage pycortex map of nprf_hcp using pp_roi.py
- each ROI masks are and saved in hdf5 files using post_pp_roi.py
- roi plots for Vision and DMM rois are plotted with bokeh using roi_plots.py
- Compute summary statistics per participants using extract_sum.py
- Draw summary statistics per participants using summary_plots.py
- polar angle progression over ANG using polar_prog.py
- make Figure S1 with post_fit/notebooks/MakeFigureS1.ipynb
