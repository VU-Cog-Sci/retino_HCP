from __future__ import division
import sys
import ctypes
import multiprocessing
import numpy as np
import popeye.utilities as utils
from popeye.visual_stimulus import VisualStimulus, simulate_bar_stimulus
import popeye.css_nohrf as css

from skimage.morphology import disk
import nibabel as nb
import matplotlib.pyplot as pl
from math import *
import json

import os
import glob
import gc
from IPython import embed as shell

from hrf_estimation.hrf import spmt  # , dspmt, ddspmt

TR = 1.0
screen_distance = 225
screen_height_degree = 8

