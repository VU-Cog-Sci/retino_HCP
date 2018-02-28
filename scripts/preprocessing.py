import os
import glob
from joblib import Parallel, delayed

from retino_HCP.utils import *

base_dir = '/home/shared/2018/visual/HCP7TFIXED'

cii_files = glob.glob(os.path.join(base_dir, '*', '*dseries.nii'))

Parallel(n_jobs=12)(delayed(sg_psc_cii)(cii_file) for cii_file in cii_files)
