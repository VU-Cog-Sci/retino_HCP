## retino_HCP
martin_dev
----------

## Analysis specifics
- HCP subjects were first pre-processed and averaged by Tomas
- Best subjects (based on R2 of V1 voxel selected by glasser2016 attlas) were chosen using data
  issued from Benson et al (http://dx.doi.org/10.1101/308247)
- pRF parameters are extracted using scripts/submit_prf_jobs.py on Lisa
- pRF derivatives are analysed using pp_roi.py
- pRF derivatives data are resampled to fsaverage using pp_roi.py
- pRF derivatives are plotted with pycortex using pp_roi.py
- ROI are determined from Conte69 atlas using pp_roi.py
- pRF derivatives are saved as a function of ROI in hdf5 files using pp_roi.py
- PRF derivatives summary are plotted with bokeh using post_pp_roi.py


