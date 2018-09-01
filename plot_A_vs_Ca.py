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
    SW_abs = 0.5 # absorptance to short_wave rad [0,1], typically 0.4-0.6

    # variables though obviously fixed here.
    par = 1500.0
    wind = 2.5
    pressure = 101325.0
    vpd = 1.5
    tair = 25
    Ca = np.linspace(0, 2000)

    CM = CoupledModel(g0, g1, D0, gamma, Vcmax25, Jmax25, Rd25,
                     Eaj, Eav,deltaSj, deltaSv, Hdv, Hdj, Q10, leaf_width,
                     SW_abs, gs_model="medlyn")


    (gs, et,
     an, Cs, Ci) = get_values2(vpd, Ca, tair, par, pressure, CM)

    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.3)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['font.size'] = 14
    plt.rcParams['legend.fontsize'] = 14
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    almost_black = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color'] = almost_black
    plt.rcParams['xtick.color'] = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color'] = almost_black

    # Change the default axis colors from black to a slightly lighter black,
    # and a little thinner (0.5 instead of 1)
    plt.rcParams['axes.edgecolor'] = almost_black
    plt.rcParams['axes.labelcolor'] = almost_black

    #colour_list = brewer2mpl.get_map('Accent', 'qualitative', 8).mpl_colors
    # CB palette  with grey:
    # from http://jfly.iam.u-tokyo.ac.jp/color/image/pallete.jpg
    colour_list = ["#CC79A7", "#E69F00", "#0072B2", "#009E73", "#F0E442",
                "#56B4E9", "#D55E00", "#000000"]

    ax1 = fig.add_subplot(111)

    ax1.plot(Ca, an, "r-", label="Ca")
    ax1.plot(Cs, an, "g-", label="Cs")
    ax1.plot(Ci, an, "b-", label="Ci")
    ax1.legend(numpoints=1, loc="best")
    ax1.set_xlim(0, 1250)

    ax1.set_ylabel("$A_{\mathrm{n}}$ ($\mathrm{\mu}$mol m$^{-2}$ s$^{-1}$)")
    ax1.set_xlabel("CO$_2$ ($\mathrm{\mu}$mol mol$^{-1}$)")
    fig.savefig("/Users/%s/Desktop/A_vs_CO2.pdf" % (os.getlogin()),
                bbox_inches='tight', pad_inches=0.1)
