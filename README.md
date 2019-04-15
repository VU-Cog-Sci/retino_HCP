## retino_HCP

ada_dev_indiv branch
--------------------
version run 15/04/2019
with hrf_delay of 0.5*TR
with gaussian sg filtered
with 3 step fit (grid + gridfit + fit)
with sub-01 / sub-02

# Duration on Lisa
------------------
grid duration = ~2h
grid_fit  = ~5h for 2500 voxels (no multiprocessing)
fit = ~13h for 2500 voxels (16 mutliprocessing)

## Analysis specifics
---------------------
- individual subjects are first preprocessed using FMRIPREP on LISA, see pre_fit/fmriprep_lisa.py
- Well registrered runs were kept for analysis, see select_block.json
- data are sg filter, psc and averaged within subjects, see pre_fit/pre_fit.py
- use pycortex_webgl to view the average time course on a pycortex server
- pRF grid is make using fit/submit_makegrid_jobs.py on LISA
- prf grid is fit using fit/launch_submit_fitgrid_jobs on LISA
- prf are finally obtained using fit/launch_submit_fitprf_jobs on LISA
- pRF derivatives of each subjects are analysed and drawn on fsaverage pycortex map of nprf_hcp using pp_roi.py
- each ROI masks are and saved in hdf5 files using post_pp_roi.py
- roi plots for Vision and DMM rois are plotted with bokeh using roi_plots.py
- Compute summary statistics per participants using extract_sum.py
- Draw summary statistics per participants using summary_plots.py
- polar angle progression over ANG using polar_prog.py


# Conclusion
# ----------
# 1. No qualitative effect of hrf delay
# 2. No qualitative effect of sg filtered model
# 3. No qualitative effect of small grid and seperate fit (except that it took more time)
# => use normal popeye procedure (faster)
# => don't use hrf delay of 0.5 TR
# => don't use gauss_sg model

# Questions
# ---------
# Why the gridfit gives offset data (wrong baseline, amplitude not perfect)
# Why the fit actually take so long? maybe wrong baseline