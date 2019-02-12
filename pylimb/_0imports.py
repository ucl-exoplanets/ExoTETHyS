from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import warnings
warnings.filterwarnings("ignore",
                        message='Matplotlib is building the font cache using fc-list. This may take a moment.')
warnings.filterwarnings("ignore",
                        message='The installed version of numexpr 2.4.4 is not supported in pandas and will be not be used')
warnings.filterwarnings("ignore",
                        message='\'second\' was found  to be \'60.0\', which is not in range [0,60). '
                                'Treating as 0 sec, +1 min [astropy.coordinates.angle_utilities]')

import os
import sys

if sys.version_info[0] > 2:
    from urllib.request import urlretrieve
else:
    from urllib import urlretrieve
    input = raw_input

import matplotlib
if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using non-interactive Agg backend')
    matplotlib.use('Agg')
else:
    matplotlib.use('TkAgg')

import glob
import time
import shutil
import astropy.io.fits as pyfits
import numpy as np
import os
from scipy.optimize import minimize
from scipy.interpolate import LinearNDInterpolator as interp

import copy
import pickle
