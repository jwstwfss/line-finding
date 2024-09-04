import copy as cp
import scipy.signal as si
import numpy as np
import matplotlib.pyplot as plt
import pickle

import os
import astropy.io.ascii as asciitable
import astropy.io.fits as fits
#import time
import mpfit
#import multiprocessing as mp
import math
from scipy import interpolate
from scipy import integrate

from trim_spec import trim_spec
from utilities import gaussian
from utilities import is_number
from utilities import read_config
#from specmodel import emissionline_model_spline
#from specmodel import model_resid_spline
from find_cwt import find_cwt
from find_cwt  import loop_field_cwt
#from find_cwt import test_obj_cwt # AA commented this out as it seems to be unused
from fitting import emissionline_model
from fitting import model_resid
from fitting import fit_obj
from fitting import get_ratio_indices
from fitting import get_fitpar_indices
#from fitting import fitandplot # MDR 2022/05/26 - Defined in fitting.py but not used so commented out.
from guis import *
from measure_z_interactive import *
import pickle
#from gather_secure_sample import *
#from measure_stack import *
#from dustfromstack import *
#from mle_stack import *



# try:
#     from stacking import *
# except ImportError:
#     pass
# #    print 'No stacking module. It is not needed for line finding. Skipping'

import warnings
warnings.filterwarnings('ignore')
