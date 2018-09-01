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

def get_values(vpd, Ca, tair, par, pressure, C):
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
    for i,ta in enumerate(tair):

        #Td = get_dewpoint(ta, rh)
        #if Td > 0.0:
        (An, gsw, et, LE, Cs, Ci) = C.main(ta, par, vpd, wind, pressure, Ca)
        gs_store.append(gsw) # mol H20 m-2 s-1
        et_store.append(et * et_conv) # mm d-1
        An_store.append(An * an_conv) # g C m-2 d-1
        Cs_store.append(Cs)
        Ci_store.append(Ci)
        tair_store.append(ta)

    return gs_store, et_store, An_store, tair_store, Cs_store, Ci_store


if __name__ == '__main__':

    fig = plt.figure(figsize=(15,4))
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.3)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 10
    plt.rcParams['font.size'] = 10
    plt.rcParams['legend.fontsize'] = 8
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10

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

    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)


    # Parameters

    # A stuff
    JV_ratio = 1.67
    Vcmax25 = 40.0          # ENF CABLE value
    Jmax25 = Vcmax25 * JV_ratio
    Rd25 = 1.0
    Eaj = 30000.0
    Eav = 60000.0
    deltaSj = 650.0
    deltaSv = 650.0
    Hdv = 200000.0
    Hdj = 200000.0
    Q10 = 2.0

    # Misc stuff
    leaf_width = 0.01
    SW_abs = 0.5 # absorptance to short_wave rad [0,1], typically 0.4-0.6

    # variables though obviously fixed here.
    par = 1500.0
    wind = 2.5
    pressure = 101325.0
    vpd = 2.0
    tair = np.linspace(0.1, 40, 20)
    Ca1 = 500.
    Ca2 = 900.
    D0 = None
    gamma = None
    g0 = 0.0
    g1 = 2.35


    CM = CoupledModel(g0, g1, D0, gamma, Vcmax25, Jmax25, Rd25,
                     Eaj, Eav,deltaSj, deltaSv, Hdv, Hdj, Q10, leaf_width,
                     SW_abs, gs_model="medlyn")


    (gs_amb, et_amb,
     an_amb, tair_2plot,
     Cs_amb, Ci_amb) = get_values(vpd, Ca1, tair, par, pressure, CM)
    (gs_ele, et_ele,
     an_ele, tair_2plot,
     Cs_ele, Ci_ele) = get_values(vpd, Ca2, tair, par, pressure, CM)

    ax1.plot(tair_2plot, an_amb, "b-")
    ax1.plot(tair_2plot, an_ele, "r-")
    ax2.plot(tair_2plot, gs_amb, "b-")
    ax2.plot(tair_2plot, gs_ele, "r-")
    ax3.plot(tair_2plot, et_amb, "b-")
    ax3.plot(tair_2plot, et_ele, "r-")


    ax1.set_ylabel("$A_{\mathrm{n}}$ (g C m$^{-2}$ d$^{-1}$)")
    ax2.set_ylabel("$g_{\mathrm{s}}$ (mol m$^{-2}$ s$^{-1}$)")
    ax3.set_ylabel("$E$ (mm d$^{-1}$)")
    ax2.set_xlabel("Air temperature ($^{\circ}$C)")

    #ax1.set_ylim(0,4)
    #ax2.set_ylim(0,4)
    #ax3.set_ylim(0,4)
    #ax4.set_ylim(0,15)
    #ax5.set_ylim(0,15)
    #ax6.set_ylim(0,15)
    #ax7.set_ylim(0,0.18)
    #ax8.set_ylim(0,0.18)
    #ax9.set_ylim(0,0.18)


    ax1.set_xlim(0,40)
    ax2.set_xlim(0,40)
    ax3.set_xlim(0,40)


    ax1.locator_params(nbins=6)
    ax2.locator_params(nbins=6)
    ax3.locator_params(nbins=6)

    fig.savefig("/Users/%s/Desktop/A_gs_E_vs_Tair.pdf" % (os.getlogin()),
                bbox_inches='tight', pad_inches=0.1)
    plt.show()
