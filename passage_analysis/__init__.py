import copy as cp
import scipy.signal as si
import numpy as np
import matplotlib.pyplot as plt
import pickle

import os
import astropy.io.ascii as asciitable
import astropy.io.fits as fits

# import time
from passage_analysis import mpfit

# import multiprocessing as mp
import math
try:
    from scipy.signal import cwt as _cwt
except:
    raise ImportError(
        "SciPy version must be <=1.14.1, as required functions"
        " were removed in later versions. Check the README for "
        "detailed installation instructions."
    )
from scipy import interpolate
from scipy import integrate

from passage_analysis.trim_spec import trim_spec
from passage_analysis.utilities import gaussian, is_number, read_config

# from specmodel import emissionline_model_spline
# from specmodel import model_resid_spline
from passage_analysis.find_cwt import find_cwt, loop_field_cwt, test_obj_cwt
from passage_analysis.fitting import (
    emissionline_model,
    model_resid,
    fit_obj,
    get_ratio_indices,
    get_fitpar_indices,
    get_broad_indices,
)

# from fitting import fitandplot # MDR 2022/05/26 - Defined in fitting.py but not used so commented out.
from passage_analysis.guis import *
from passage_analysis.measure_z_interactive import *
import pickle

# from gather_secure_sample import *
# from measure_stack import *
# from dustfromstack import *
# from mle_stack import *


# try:
#     from stacking import *
# except ImportError:
#     pass
# #    print 'No stacking module. It is not needed for line finding. Skipping'
