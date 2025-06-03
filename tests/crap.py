#!/usr/bin/env python

"""
Compare transpiration sensitivity to PFT difference in g1 vs. inc. temp/VPD

This makes the plot in the sci reports paper, but with the low temp ramp down
turned off.

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (23.07.2015)"
__email__ = "mdekauwe@gmail.com"

import sys
import numpy as np
import os
import math
import matplotlib.pyplot as plt


from farq import FarquharC3
from solve_coupled_An_gs_leaf_temp_transpiration import CoupledModel
from utils import vpd_to_rh, get_dewpoint, calc_esat
import constants as c

tair = 30.0
rh = 95.0

esat = calc_esat(tair)
ea = rh / 100. * esat
vpd = (esat - ea) * c.PA_2_KPA


print(vpd)

tair = 40.0
rh = 90.0

esat = calc_esat(tair)
ea = rh / 100. * esat
vpd = (esat - ea) * c.PA_2_KPA


print(vpd)
