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
import pandas as pd

from farq import FarquharC3
from solve_coupled_An_gs_leaf_temp_transpiration import CoupledModel
from utils import vpd_to_rh, get_dewpoint, calc_esat
import constants as c



# A stuff
JV_ratio = 1.67
Vcmax25 = 58.2          # ENF CABLE value
Jmax25 = Vcmax25 * JV_ratio
Rd25 = 1.0
Eaj = 39676.89#30000.0
Eav = 82620.87#60000.0
deltaSj = 650.0
deltaSv = 650.0
Hdv = 200000.0
Hdj = 200000.0
Q10 = 2.0

# Misc stuff
leaf_width = 0.08
# Cambell & Norman, 11.5, pg 178
# The solar absorptivities of leaves (-0.5) from Table 11.4 (Gates, 1980)
# with canopies (~0.8) from Table 11.2 reveals a surprising difference.
# The higher absorptivityof canopies arises because of multiple reflections
# among leaves in a canopy and depends on the architecture of the canopy.
SW_abs = 0.8 # use canopy absorptance of solar radiation

# variables though obviously fixed here.

pressure = 101325.0
D0 = None
gamma = None
g0 = 0.0
g1 = 1.72

df = pd.read_csv("/Users/mdekauwe/Downloads/Rocca1_met_and_plant_data_drought_2003.csv")

CM = CoupledModel(g0, g1, D0, gamma, Vcmax25, Jmax25, Rd25,
                 Eaj, Eav,deltaSj, deltaSv, Hdv, Hdj, Q10, leaf_width,
                 SW_abs, gs_model="medlyn")

"""
gs_store = []
et_store = []
An_store = []
et_conv = c.MOL_WATER_2_G_WATER * c.G_TO_KG * c.SEC_TO_DAY
an_conv = c.UMOL_TO_MOL * c.MOL_C_TO_GRAMS_C * c.SEC_TO_DAY

Anx = 0.0
etx = 0.0
count = 0
for i in range(len(df)):


    (An, gsw, et, LE,
     Cs, Ci) = CM.main(df.TAir[i], df.PAR[i], df.VPD[i], df["u"][i], pressure, 400.)


    Anx = An #* c.UMOL_TO_MOL * c.MOL_C_TO_GRAMS_C * 1800.
    etx = et #* c.MOL_WATER_2_G_WATER * c.G_TO_KG * 1800.


    An_store.append(Anx) # g C m-2 d-1
    et_store.append(etx ) # mm d-1


    count += 1

fig = plt.figure(figsize=(14,5))
fig.subplots_adjust(hspace=0.1)
fig.subplots_adjust(wspace=0.1)
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.sans-serif'] = "Helvetica"
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['font.size'] = 10
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.plot(An_store[:47*4])
ax2.plot(et_store[:47*4])
plt.show()
sys.exit()
"""

gs_store = []
et_store = []
An_store = []
et_conv = c.MOL_WATER_2_G_WATER * c.G_TO_KG * c.SEC_TO_DAY
an_conv = c.UMOL_TO_MOL * c.MOL_C_TO_GRAMS_C * c.SEC_TO_DAY
Anx = 0.0
etx = 0.0
count = 0
for i in range(len(df)):


    (An, gsw, et, LE,
     Cs, Ci) = CM.main(df.TAir[i], df.PAR[i], df.VPD[i], df["u"][i], pressure, 400.)


    Anx += An * c.UMOL_TO_MOL * c.MOL_C_TO_GRAMS_C * 1800.
    etx += et * c.MOL_WATER_2_G_WATER * c.G_TO_KG * 1800.

    if count == 47:
        An_store.append(Anx * df.LAI[i]) # g C m-2 d-1
        et_store.append(etx * df.LAI[i]) # mm d-1
        Anx = 0.0
        etx = 0.0
        count = 0


    count += 1

fig = plt.figure(figsize=(14,5))
fig.subplots_adjust(hspace=0.1)
fig.subplots_adjust(wspace=0.1)
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.sans-serif'] = "Helvetica"
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['font.size'] = 10
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.plot(An_store)
ax2.plot(et_store)
plt.show()
