## retino_HCP

ada_dev_subcortical_bars
------------------------

## Analysis specifics
---------------------
- data is pre-processed (sg filter + psc) using pre_fit/pre_fit.py psc is computed using blank period of each run
- data is averaged across subjects using pre_fit/pre_fit.py
- data are modeled in lisa using fit/prf_fit.py run with fit/launch_submit_fit_jobs.py and fit/submit_fit_jobs.py
- combine fits, convert to fsaverage, determine ROI and plot flatmaps with pp_roi.py
