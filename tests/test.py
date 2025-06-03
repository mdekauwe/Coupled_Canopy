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

def get_values2(vpd, Ca, tair, par, pressure, C):
    kpa_2_pa = 1000.
    pa_2_kpa = 1.0 / kpa_2_pa

    #print rh, vpd
    gs_store = []
    et_store = []
    An_store = []
    tair_store = []
    Cs_store = []
    Ci_store = []
    et_conv = c.MOL_WATER_2_G_WATER * c.G_TO_KG * c.SEC_TO_DAY
    an_conv = c.UMOL_TO_MOL * c.MOL_C_TO_GRAMS_C * c.SEC_TO_DAY
    for i,cax in enumerate(Ca):

        (An, gsw, et, LE, Cs, Ci) = C.main(tair, par, vpd, wind, pressure, cax)
        gs_store.append(gsw) # mol H20 m-2 s-1
        et_store.append(et * et_conv) # mm d-1
        An_store.append(An) # umol m-2 s-1
        Cs_store.append(Cs)
        Ci_store.append(Ci)

    return gs_store, et_store, An_store, Cs_store, Ci_store

if __name__ == '__main__':




    # Parameters

    # A stuff
    JV_ratio = 2.0
    Vcmax25 = 50.0          # ENF CABLE value
    Jmax25 = Vcmax25 * JV_ratio
    Rd25 = 0.92
    Eaj = 30000.0
    Eav = 60000.0
    deltaSj = 650.0
    deltaSv = 650.0
    Hdv = 200000.0
    Hdj = 200000.0
    Q10 = 1.92
    D0 = None
    gamma = None
    g0 = 0.0
    g1 = 2
    # Misc stuff
    leaf_width = 0.01

    # Cambell & Norman, 11.5, pg 178
    # The solar absorptivities of leaves (-0.5) from Table 11.4 (Gates, 1980)
    # with canopies (~0.8) from Table 11.2 reveals a surprising difference.
    # The higher absorptivityof canopies arises because of multiple reflections
    # among leaves in a canopy and depends on the architecture of the canopy.
    SW_abs = 0.8 # use canopy absorptance of solar radiation

    # variables though obviously fixed here.
    par = 1500.0
    wind = 2.5
    pressure = 101325.0
    vpd = 6.
    tair = 28.
    Ca = np.linspace(350, 2000)

    CM = CoupledModel(g0, g1, D0, gamma, Vcmax25, Jmax25, Rd25,
                     Eaj, Eav,deltaSj, deltaSv, Hdv, Hdj, Q10, leaf_width,
                     SW_abs, gs_model="medlyn")


    (gs, et,
     an, Cs, Ci) = get_values2(vpd, Ca, tair, par, pressure, CM)

    print("vpd =", round(vpd,1), "tair =", round(tair,1), "An =", round(np.mean(an),2), "E =", round(np.mean(et),2))
