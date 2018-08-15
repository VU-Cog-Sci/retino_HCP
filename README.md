## retino_HCP
[![Build Status](https://travis-ci.org/uwescience/retino_HCP.svg?branch=master)](https://travis-ci.org/uwescience/retino_HCP)

PRF fitting project for HCP data




### Cartesius cluster
```rsync -a --no-g --no-p -vzhe "ssh -T -c aes128-ctr -o Compression=no -x" --progress --include="*_av.gii" --include='*/' --exclude='*' /home/shared/2018/visual/HCP7TFIXED $USER@cartesius.surfsara.nl:/projects/0/pqsh283/```

### Lisa cluster
```rsync -a --no-g --no-p -vzhe "ssh -T -c aes128-ctr -o Compression=no -x" --progress --include="*_av.gii" --include='*/' --exclude='*' /home/shared/2018/visual/HCP7TFIXED $USER@lisa.surfsara.nl:/home/knapen/data/```


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


