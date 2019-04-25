## retino_HCP

martin_dev branch
-----------------
code for analysing HCP dataset

## Analysis specifics
---------------------
- HCP subjects were first pre-processed and averaged
- pRF parameters are extracted using fit/submit_fit_jobs.py on Lisa with gaussian model
- get yeo atlas see others/get_yeo_atlas_gii.txt
- DMN regions (ANG/MED_PAR/SUP_MED_FR/LAT_TEMP) are put on fsaverage flatmap overlay.svg
- DMN regions of interest are drawn in inkscape of fsaverage post_fit/add_dmn_roi.py
- pRF derivatives of subject '999999' are analysed and drawn on fsaverage pycortex map using pp_roi.py
- Vision regions of interest (V1/V2/V3/VO/DO/LO/SUP_PAR/TPJ/sPCS/iPCS/mPCS/INS/DLPFC) are drawn manually in inkscape
- PRF derivatives summary for each ROI are put in h5 files with post_pp_roi.py
- roi plots for DMN and Vision ROI are plotted with bokeh using roi_plots.py for subject 999999
- pRF derivatives of all others subject are analysed using pp_roi.py and drawn using fsaverage and ROI of 999999
- roi plots for of all others subject for DMN and Vision ROI are plotted with bokeh using roi_plots.py
- polar angle progression over ANG using polar_prog.py
- make Figure 1C using post_fit/notebooks/MakeFigure1C.ipynb
- make Figure 2A using post_fit/notebooks/MakeFigure2A.ipynb
- make Figure 2B using post_fit/notebooks/MakeFigure2B.ipynb
- make Figure 3 and S3 using post_fit/notebooks/MakeFigure3.ipynb