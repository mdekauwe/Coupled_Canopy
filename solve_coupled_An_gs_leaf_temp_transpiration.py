#!/usr/bin/env python

"""
Iteratively solve leaf temp, ci, gs, An and transpiration following Maetra
looping logic


That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (23.07.2015)"
__email__ = "mdekauwe@gmail.com"

import sys
import numpy as np
import os
import math

import constants as c
from farq import FarquharC3
from penman_monteith_leaf import PenmanMonteith

class CoupledModel(object):
    """Iteratively solve leaf temp, Ci, gs and An."""

    def __init__(self, g0, g1, D0, gamma, Vcmax25, Jmax25, Rd25, Eaj, Eav,
                 deltaSj, deltaSv, Hdv, Hdj, Q10, leaf_width, SW_abs,
                 gs_model, alpha=None, leaf_absorptance=0.5, iter_max=100):

        # set params
        self.g0 = g0
        self.g1 = g1
        self.D0 = D0
        self.gamma = gamma
        self.Vcmax25 = Vcmax25
        self.Jmax25 = Jmax25
        self.Rd25 = Rd25
        self.Eaj = Eaj
        self.Eav = Eav
        self.deltaSj = deltaSj
        self.deltaSv = deltaSv
        self.Hdv = Hdv
        self.Hdj = Hdj
        self.Q10 = Q10
        self.leaf_width = leaf_width
        self.alpha = alpha
        self.SW_abs = SW_abs
        self.gs_model = gs_model
        self.iter_max = iter_max

        self.emissivity_leaf = 0.99   # emissivity of leaf (-)
        self.leaf_absorptance = leaf_absorptance # leaf abs of solar rad [0,1]


    def main(self, tair, par, vpd, wind, pressure, Ca):
        """
        Parameters:
        ----------
        tair : float
            air temperature (deg C)
        par : float
            Photosynthetically active radiation (umol m-2 s-1)
        vpd : float
            Vapour pressure deficit (kPa, needs to be in Pa, see conversion
            below)
        wind : float
            wind speed (m s-1)
        pressure : float
            air pressure (using constant) (Pa)
        Ca : float
            ambient CO2 concentration

        Returns:
        --------
        An : float
            net leaf assimilation (umol m-2 s-1)
        gs : float
            stomatal conductance (mol m-2 s-1)
        et : float
            transpiration (mol H2O m-2 s-1)
        """

        F = FarquharC3(peaked_Jmax=True, peaked_Vcmax=True, model_Q10=True,
                       gs_model=self.gs_model, gamma=self.gamma, g0=self.g0,
                       g1=self.g1, D0=self.D0, alpha=self.alpha)
        P = PenmanMonteith(self.leaf_width, self.leaf_absorptance)

        # set initialise values
        dleaf = vpd
        dair = vpd
        Cs = Ca
        Tleaf = tair
        Tleaf_K = Tleaf + c.DEG_2_KELVIN

        #print "Start: %.3f %.3f %.3f" % (Cs, Tleaf, dleaf)
        #print


        iter = 0
        while True:
            (An, gsc) = F.calc_photosynthesis(Cs=Cs, Tleaf=Tleaf_K, Par=par,
                                              Jmax25=self.Jmax25,
                                              Vcmax25=self.Vcmax25,
                                              Q10=self.Q10, Eaj=self.Eaj,
                                              Eav=self.Eav,
                                              deltaSj=self.deltaSj,
                                              deltaSv=self.deltaSv,
                                              Rd25=self.Rd25, Hdv=self.Hdv,
                                              Hdj=self.Hdj, vpd=dleaf)



            # Calculate new Tleaf, dleaf, Cs
            (new_tleaf, et, le_et, gbH, gw) = self.calc_leaf_temp(P, Tleaf, tair, gsc,
                                                           par, vpd, pressure,
                                                           wind)


            gbc = gbH * c.GBH_2_GBC
            Cs = Ca - An / gbc # boundary layer of leaf
            if et == 0.0 or gw == 0.0:
                dleaf = dair
            else:
                dleaf = (et * pressure / gw) * c.PA_2_KPA # kPa


            #print "%f %f %f %f %f %f" %  (Cs, Tleaf, dleaf, An*12.*0.000001*86400., gs, et*18*0.001*86400.)

            # Check for convergence...?
            if math.fabs(Tleaf - new_tleaf) < 0.02:
                break

            if iter > self.iter_max:
                raise Exception('No convergence: %d' % (iter))

            # Update temperature & do another iteration
            Tleaf = new_tleaf
            Tleaf_K = Tleaf + c.DEG_2_KELVIN

            iter += 1

        gsw = gsc * c.GSC_2_GSW

        return (An, gsw, et, le_et)

    def calc_leaf_temp(self, P=None, tleaf=None, tair=None, gsc=None, par=None,
                       vpd=None, pressure=None, wind=None):
        """
        Resolve leaf temp

        Parameters:
        ----------
        P : object
            Penman-Montheith class instance
        tleaf : float
            leaf temperature (deg C)
        tair : float
            air temperature (deg C)
        gs : float
            stomatal conductance (mol m-2 s-1)
        par : float
            Photosynthetically active radiation (umol m-2 s-1)
        vpd : float
            Vapour pressure deficit (kPa, needs to be in Pa, see conversion
            below)
        pressure : float
            air pressure (using constant) (Pa)
        wind : float
            wind speed (m s-1)

        Returns:
        --------
        new_Tleaf : float
            new leaf temperature (deg C)
        et : float
            transpiration (mol H2O m-2 s-1)
        gbH : float
            total boundary layer conductance to heat for one side of the leaf
        gw : float
            total leaf conductance to water vapour (mol m-2 s-1)
        """
        tleaf_k = tleaf + c.DEG_2_KELVIN
        tair_k = tair + c.DEG_2_KELVIN

        air_density = pressure / (c.RSPECIFC_DRY_AIR * tair_k)

        # convert from mm s-1 to mol m-2 s-1
        cmolar = pressure / (c.RGAS * tair_k)

        # W m-2 = J m-2 s-1
        rnet = P.calc_rnet(par, tair, tair_k, tleaf_k, vpd, pressure)

        (grn, gh, gbH, gw) = P.calc_conductances(tair_k, tleaf, tair,
                                                 wind, gsc, cmolar)
        if gsc == 0.0:
            et = 0.0
            le_et = 0.0
        else:
            (et, le_et) = P.calc_et(tleaf, tair, vpd, pressure, wind, par,
                                    gh, gw, rnet)

        # D6 in Leuning
        Y = 1.0 / (1.0 + grn / gbH)

        # sensible heat exchanged between leaf and surroundings
        sensible_heat = Y * (rnet - le_et)

        # leaf-air temperature difference recalculated from energy balance.
        delta_T = (rnet - le_et) / (c.CP * c.AIR_MASS * gh)
        new_Tleaf = tair + delta_T

        return (new_Tleaf, et, le_et, gbH, gw)

if __name__ == '__main__':

    # Parameters

    # gs stuff
    g0 = 0.001
    g1 = 9.0
    D0 = 1.5 # kpa

    # A stuff
    Vcmax25 = 30.0
    Jmax25 = Vcmax25 * 2.0
    Rd25 = 2.0
    Eaj = 30000.0
    Eav = 60000.0
    deltaSj = 650.0
    deltaSv = 650.0
    Hdv = 200000.0
    Hdj = 200000.0
    Q10 = 2.0
    gamma = 0.0

    # Misc stuff
    leaf_width = 0.02

    # hack to get around doing seperate sunlit/shaded leaves
    leaf_absorptance = 0.8 # leaf absorptance of solar radiation [0,1]

    # variables though obviously fixed here.
    par = 1500.0
    tair = 25.0
    vpd = 1.0
    wind = 2.5
    pressure = 101325.0
    Ca = 400.0

    C = CoupledModel(g0, g1, D0, gamma, Vcmax25, Jmax25, Rd25, Eaj, Eav, deltaSj,
                     deltaSv, Hdv, Hdj, Q10, leaf_width, leaf_absorptance,
                     gs_model="leuning")
    (An, gsw, et, le_et) = C.main(tair, par, vpd, wind, pressure, Ca)
