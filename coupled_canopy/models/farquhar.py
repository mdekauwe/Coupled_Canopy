#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Model of C3 photosynthesis, this is passed to fitting function and we are
optimising Jmax25, Vcmax25, Rd, Eaj, Eav, deltaS

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (13.08.2012)"
__email__ = "mdekauwe@gmail.com"

import sys
import numpy as np
import os
import math
from coupled_canopy.utils import constants as c

class FarquharC3(object):
    """
    C3 photosynthesis model based on Farquhar et al. (1980)

    All calculations in Kelvins...

    References:
    -----------
    * De Pury and Farquhar, G. D. (1997) Simple scaling of photosynthesis from
      leaves to canopies without the errors of big-leaf models. Plant Cell and
      Environment, 20, 537-557.
    * Farquhar, G.D., Caemmerer, S. V. and Berry, J. A. (1980) A biochemical
      model of photosynthetic CO2 assimilation in leaves of C3 species. Planta,
      149, 78-90.
    * Medlyn, B. E., Dreyer, E., Ellsworth, D., Forstreuter, M., Harley, P.C.,
      Kirschbaum, M.U.F., Leroux, X., Montpied, P., Strassemeyer, J.,
      Walcroft, A., Wang, K. and Loustau, D. (2002) Temperature response of
      parameters of a biochemically based model of photosynthesis. II.
      A review of experimental data. Plant, Cell and Enviroment 25, 1167-1179.
    """

    def __init__(self, peaked_Jmax=False, peaked_Vcmax=False, Oi=210.0,
                 gamstar25=42.75, Kc25=404.9, Ko25=278.4, Ec=79430.0,
                 Eo=36380.0, Eag=37830.0, theta_hyperbol=0.9995,
                 theta_J=0., force_vcmax_fit_pts=None,
                 alpha=None, quantum_yield=0.3, absorptance=0.8,
                 change_over_pt=None, model_Q10=False,
                 gs_model=None, gamma=0.0, g0=None, g1=None, D0=1.5,
                 adjust_Vcmax_Jmax_for_low_temp=False):
        """
        Parameters
        ----------
        Oi : float
            intercellular concentration of O2 [mmol mol-1]
        gamstar25 : float
            CO2 compensation point - base rate at 25 deg C / 298 K [umol mol-1]
        Kc25 : float
            Michaelis-Menten coefficents for carboxylation by Rubisco at
            25degC [umol mol-1] or 298 K
        Ko25: float
            Michaelis-Menten coefficents for oxygenation by Rubisco at
            25degC [mmol mol-1]. Note value in Bernacchie 2001 is in mmol!!
            or 298 K
        Ec : float
            Activation energy for carboxylation [J mol-1]
        Eo : float
            Activation energy for oxygenation [J mol-1]
        Eag : float
            Activation energy at CO2 compensation point [J mol-1]
        RGAS : float
            Universal gas constant [J mol-1 K-1]
        theta_hyperbol : float
            Curvature of the light response.
            See Peltoniemi et al. 2012 Tree Phys, 32, 510-519
        theta_J : float
            Curvature of the light response
        alpha : float
            Leaf quantum yield (initial slope of the A-light response curve)
            [mol mol-1]
        peaked_Jmax : logical
            Use the peaked Arrhenius function (if true)
        peaked_Vcmax : logical
            Use the peaked Arrhenius function (if true)

        force_vcmax_fit_pts : None or npts
            Force Ac fit for first X points
        change_over_pt : None or value of Ci
            Explicitly set the transition point between Aj and Ac.

        gs_model : sting
            stomatal conductance model - Leuning/Medlyn
        gamma : float
            is the CO2 compensation point of photosynthesis (umol m-2 s-1)
        g0 : float
            residual stomatal conductance as net assimilation rate reaches
            zero (mol m-2 s-1)
        g1 : float
            and the slope of the sensitivity of stomatal conductance to
            assimilation (mol m-2 s-1)
        D0 : float
            the sensitivity of stomatal conductance to D (kPa)
        """
        self.peaked_Jmax = peaked_Jmax
        self.peaked_Vcmax = peaked_Vcmax
        self.Oi = Oi
        self.gamstar25 = gamstar25
        self.Kc25 = Kc25
        self.Ko25 = Ko25
        self.Ec = Ec
        self.Eo = Eo
        self.Eag = Eag
        self.theta_hyperbol = theta_hyperbol
        self.theta_J = theta_J
        if alpha is not None:
            self.alpha = alpha
        else:
            self.alpha = quantum_yield * absorptance # (Medlyn et al 2002)

        self.force_vcmax_fit_pts = force_vcmax_fit_pts
        self.change_over_pt = change_over_pt
        self.model_Q10 = model_Q10
        self.gs_model = gs_model
        self.gamma = gamma
        self.g0 = g0
        self.g1 = g1
        self.D0 = D0
        self.adjust_Vcmax_Jmax_for_low_temp = adjust_Vcmax_Jmax_for_low_temp

    def calc_photosynthesis(self, Cs=None, Tleaf=None, Par=None,
                            Jmax=None, Vcmax=None, Jmax25=None, Vcmax25=None,
                            Rd=None, Rd25=None, Q10=None, Eaj=None, Eav=None,
                            deltaSj=None, deltaSv=None, Hdv=200000.0,
                            Hdj=200000.0, Ear=None, vpd=None, mult=None):
        """
        Parameters
        ----------
        Cs : float
            leaf surface CO2 concentration [umol mol-1]
        Tleaf : float
            leaf temp [deg K]
        * Optional args:
        Jmax : float
            potential rate of electron transport at measurement temperature
            [deg K]
        Vcmax : float
            max rate of rubisco activity at measurement temperature [deg K]
        Jmax25 : float
            potential rate of electron transport at 25 deg or 298 K
        Vcmax25 : float
            max rate of rubisco activity at 25 deg or 298 K
        Rd : float
            Day "light" respiration [umol m-2 time unit-1]
        Q10 : float
            ratio of respiration at a given temperature divided by respiration
            at a temperature 10 degrees lower
        Eaj : float
            activation energy for the parameter [J mol-1]
        Eav : float
            activation energy for the parameter [J mol-1]
        deltaSj : float
            entropy factor [J mol-1 K-1)
        deltaSv : float
            entropy factor [J mol-1 K-1)
        HdV : float
            Deactivation energy for Vcmax [J mol-1]
        Hdj : float
            Deactivation energy for Jmax [J mol-1]
        Rd25 : float
            Estimate of respiration rate at the reference temperature 25 deg C
             or 298 K [deg K]
        Par : float
            photosynthetically active radiation [umol m-2 time unit-1]. Default
            is not to supply PAR, with measurements taken under light
            saturation.

        Returns:
        --------
        An : float
            Net leaf assimilation rate [umol m-2 s-1]
        gsc : float
            stomatal conductance to CO2 [mol m-2 s-1]
        gsw : float
            stomatal conductance to water vapour [mol H2O m-2 s-1]
        """
        self.check_supplied_args(Jmax, Vcmax, Rd, Jmax25, Vcmax25, Rd25)

        # Calculate temperature-dependent parameters
        (Km, gamma_star,
        Vcmax, Jmax,
        Rd) = self.calc_temperature_dependent_params(Tleaf,Q10, Rd25, Ear,
                                                     Vcmax25, Eav, deltaSv, Hdv,
                                                     Jmax25, Eaj, deltaSj, Hdj)

        (gs_over_a, g0) = self.calc_stomatal_ratio(Cs, vpd, mult)

        # Rate of electron transport, which is a function of absorbed PAR
        if Par is not None:
            J = self.calc_electron_transport_rate(Par, Jmax)
        # all measurements are calculated under saturated light!!
        else:
            J = Jmax
        Vj = J / 4.0


        if ( math.isclose(Par, 0.0) | math.isclose(Vj, 0.0) ):
            Cic = Cs
            Cij = Cs
        else:
            # Solution when Rubisco activity is limiting
            (Cic) = self.solve_ci(g0, gs_over_a, Rd, Cs, gamma_star, Vcmax, Km)

            # Solution when electron transport rate is limiting
            (Cij) = self.solve_ci(g0, gs_over_a, Rd, Cs, gamma_star, Vj,
                                  2.0*gamma_star)

        if Cic <= 0.0 or Cic > Cs:
            Ac = 0.0
        else:
            Ac = self.assim(Cic, gamma_star, a1=Vcmax, a2=Km)
        Aj = self.assim(Cij, gamma_star, a1=Vj, a2=2.0*gamma_star)

        # When below light-compensation points, assume Ci=Ca.
        if Aj <= Rd + 1E-09:
            Cij = Cs
            Aj = self.assim(Cij, gamma_star, a1=Vj, a2=2.0*gamma_star)

        # Hyperbolic minimum.
        A = -self.quadratic(a=1.0 - 1E-04, b=Ac + Aj, c=Ac * Aj, large=True)

        # Net photosynthesis
        An = A - Rd

        # Calculate conductance to CO2
        gsc = max(g0, g0 + gs_over_a * An)

        # Extra step here; GS can be negative
        if gsc < g0:
            gsc = g0

        # Calculate conductance to water
        gsw = gsc * c.GSC_2_GSW

        # calculate the real Ci
        if gsc > 0.0 and An > 0.0:
            Ci = Cs - An / gsc
        else:
            Ci = Cs

        if math.isclose(Cs, 0.0):
            An = 0.0 - Rd
            gsc = 0.0
            Ci = Cs

        return (An, gsc, Ci)

    def calc_stomatal_ratio(self, Cs, vpd, mult):
        """
        Calculate stomatal ratio: gs/A based on specified model.

        Parameters:
        --------
        Cs: : float
            CO2 concentration at leaf surface (µmol mol-1)
        VPD: float:
            vapor pressure deficit (kPa)

        Returns:
        --------
        gs_over_a: float
            Ratio of gs to A (mol m-2 s-2 per µmol m-2 s-1)
        g0: float
            Residual conductance to CO₂ (mol m-2 s-1)
        """
        # Apply numerical fix for g0
        g0 = self.g0 * c.GSW_2_GSC if not math.isclose(self.g0, 0.0) else 1E-09

        if self.gs_model == 'medlyn':

            # Prevent division by very low VPD
            vpd = max(vpd, 0.05)

            # 1.6 (from corrigendum to Medlyn et al 2011) is missing here,
            # because we are calculating conductance to CO2!
            if math.isclose(Cs, 0.0):
                gs_over_a = 0.0
            else:
                gs_over_a = (1.0 + self.g1 / math.sqrt(vpd)) / Cs
            ci_over_ca = self.g1 / (self.g1 + math.sqrt(vpd))

        elif self.gs_model == "leuning":

            gs_over_a = self.g1 / (Cs - self.gamma) / (1.0 + vpd / self.D0)

            # conductance to CO2
            gs_over_a *= c.GSW_2_GSC
            ci_over_ca = 1.0 - 1.6 * (1.0 + vpd / self.D0) / self.g1

        #elif self.gs_model == "ball_berry":
        #
        #    if RH is None:
        #        raise ValueError("Ball-Berry model requires relative_humidity (0–1).")
        #    gs_over_a = (self.g1 * RH) / Cs
        #    gs_over_a *= c.GSW_2_GSC
        #
        #elif self.gs_model == "tuzet":
        #    if RH is None:
        #        raise ValueError("Tuzet model requires relative_humidity (0–1).")
        #
        #    f = RH / (RH + math.exp(self.a * (self.RH0 - RH)))
        #    gs_over_a = (self.g1 * f) / Cs
        #    gs_over_a *= c.GSW_2_GSC


        elif self.gs_model == "user_defined":

            # Multiplier is user-defined.
            g0 = self.g0 / c.GSC_2_GSW
            gs_over_a = mult / c.GSC_2_GSW

        else:
            # Custom or unknown model
            raise NotImplementedError(
             f"Stomatal conductance model '{self.gs_model}' not implemented.")

        return gs_over_a, g0

    def calc_temperature_dependent_params(self, Tleaf, Q10, Rd25, Ear, Vcmax25,
                                          Eav, deltaSv, Hdv, Jmax25, Eaj,
                                          deltaSj, Hdj):
        """
        Calculate temperature-dependent photosynthetic parameters:
        """

        # Calculate Michaelis-Menten constant for O2/CO2
        Kc = self._arrhenius(self.Kc25, self.Ec, Tleaf)
        Ko = self._arrhenius(self.Ko25, self.Eo, Tleaf)
        Km = Kc * (1.0 + self.Oi / Ko)

        # Effect of temp on CO2 compensation point
        gamma_star = self._arrhenius(self.gamstar25, self.Eag, Tleaf)

        # Calculations at 25 degrees C or the measurement temperature
        if Rd25 is not None:
            Rd = self.calc_resp(Tleaf, Q10, Rd25, Ear)

        # Calculate Vcmax and Jmax temperature responses
        if Vcmax25 is not None:
            # Effect of temperature on Vcmax and Jamx
            if self.peaked_Vcmax:
                Vcmax = self._peaked_arrhenius(Vcmax25, Eav, Tleaf, deltaSv, Hdv)
            else:
                Vcmax = self._arrhenius(Vcmax25, Eav, Tleaf)

        if Jmax25 is not None:
            if self.peaked_Jmax:
                Jmax = self._peaked_arrhenius(Jmax25, Eaj, Tleaf, deltaSj, Hdj)
            else:
                Jmax = self._arrhenius(Jmax25, Eaj, Tleaf)

        if self.adjust_Vcmax_Jmax_for_low_temp:
            Jmax = self.adj_for_low_temp(Jmax, Tleaf)
            Vcmax = self.adj_for_low_temp(Vcmax, Tleaf)

        return Km, gamma_star, Vcmax, Jmax, Rd

    def adj_for_low_temp(self, param, Tk, lower_bound=0.0, upper_bound=10.0):
        """
        Function allowing Jmax/Vcmax to be forced linearly to zero at low T

        Parameters:
        ----------
        Tk : float
            air temperature (Kelvin)
        """
        Tc = Tk - c.DEG_2_KELVIN

        if Tc < lower_bound:
            param = 0.0
        elif Tc < upper_bound:
            param *= (Tc - lower_bound) / (upper_bound - lower_bound)

        return param

    def _arrhenius(self, k25, Ea, Tk):
        """ Temperature dependence of kinetic parameters is described by an
        Arrhenius function.

        Parameters:
        ----------
        k25 : float
            rate parameter value at 25 degC or 298 K
        Ea : float
            activation energy for the parameter [J mol-1]
        Tk : float
            leaf temperature [deg K]

        Returns:
        -------
        kt : float
            temperature dependence on parameter

        References:
        -----------
        * Medlyn et al. 2002, PCE, 25, 1167-1179.
        """
        return k25 * np.exp((Ea * (Tk - 298.15)) / (298.15 * c.RGAS * Tk))

    def _peaked_arrhenius(self, k25, Ea, Tk, deltaS, Hd):
        """ Temperature dependancy approximated by peaked Arrhenius eqn,
        accounting for the rate of inhibition at higher temperatures.

        Parameters:
        ----------
        k25 : float
            rate parameter value at 25 degC or 298 K
        Ea : float
            activation energy for the parameter [J mol-1]
        Tk : float
            leaf temperature [deg K]
        deltaS : float
            entropy factor [J mol-1 K-1)
        Hd : float
            describes rate of decrease about the optimum temp [J mol-1]

        Returns:
        -------
        kt : float
            temperature dependence on parameter

        References:
        -----------
        * Medlyn et al. 2002, PCE, 25, 1167-1179.

        """
        arg1 = self._arrhenius(k25, Ea, Tk)
        arg2 = 1.0 + np.exp((298.15 * deltaS - Hd) / (298.15 * c.RGAS))
        arg3 = 1.0 + np.exp((Tk * deltaS - Hd) / (Tk * c.RGAS))

        return arg1 * arg2 / arg3

    def calc_electron_transport_rate(self, Par, Jmax):
        """
        Electron transport rate for a given absorbed irradiance

        Reference:
        ----------
        * Farquhar G.D. & Wong S.C. (1984) An empirical model of stomatal
          conductance. Australian Journal of Plant Physiology 11, 191-210,
          eqn A but probably clearer in:
        * Leuning, R. et a., Leaf nitrogen, photosynthesis, conductance and
          transpiration: scaling from leaves to canopies, Plant Cell Environ.,
          18, 1183– 1200, 1995. Leuning 1995, eqn C3.
        """
        A = self.theta_J
        B = -(self.alpha * Par + Jmax);
        C = self.alpha * Par * Jmax;

        J = self.quadratic(a=A, b=B, c=C, large=False)

        return J

    def solve_ci(self, g0, gs_over_a, rd, Cs, gamma_star, gamma, beta):
        """
        Solve intercellular CO2 concentration using quadric equation, following
        Leuning 1990, see eqn 15a-c, solving simultaneous solution for Eqs 2, 12
        and 13

        Reference:
        ----------
        * Leuning (1990) Modelling Stomatal Behaviour and Photosynthesis of
          Eucalyptus grandis. Aust. J. Plant Physiol., 17, 159-75.
        """

        A = g0 + gs_over_a * (gamma - rd)

        arg1 = (1. - Cs * gs_over_a) * (gamma - rd)
        arg2 = g0 * (beta - Cs)
        arg3 = gs_over_a * (gamma * gamma_star + beta * rd)
        B = arg1 + arg2 - arg3

        arg1 = -(1.0 - Cs * gs_over_a)
        arg2 = (gamma * gamma_star + beta * rd)
        arg3 =  g0 * beta * Cs
        C = arg1 * arg2 - arg3

        Ci = self.quadratic(a=A, b=B, c=C, large=True)

        return Ci

    def assim(self, Ci, gamma_star, a1, a2):
        """calculation of photosynthesis with the limitation defined by the
        variables passed as a1 and a2, i.e. if we are calculating vcmax or
        jmax limited assimilation rates.

        Parameters:
        ----------
        Ci : float
            intercellular CO2 concentration.
        gamma_star : float
            CO2 compensation point in the abscence of mitochondrial respiration
        a1 : float
            variable depends on whether the calculation is light or rubisco
            limited.
        a2 : float
            variable depends on whether the calculation is light or rubisco
            limited.

        Returns:
        -------
        assimilation_rate : float
            assimilation rate assuming either light or rubisco limitation.
        """
        return a1 * (Ci - gamma_star) / (a2 + Ci)

    def calc_resp(self, Tleaf=None, Q10=None, Rd25=None, Ear=None, Tref=25.0):
        """ Calculate leaf respiration accounting for temperature dependence.

        Parameters:
        ----------
        Rd25 : float
            Estimate of respiration rate at the reference temperature 25 deg C
            or or 298 K
        Tref : float
            reference temperature
        Q10 : float
            ratio of respiration at a given temperature divided by respiration
            at a temperature 10 degrees lower
        Ear : float
            activation energy for the parameter [J mol-1]
        Returns:
        -------
        Rt : float
            leaf respiration

        References:
        -----------
        Tjoelker et al (2001) GCB, 7, 223-230.
        """
        if self.model_Q10:
            Rd = Rd25 * Q10**(((Tleaf - c.DEG_2_KELVIN) - Tref) / 10.0)
        else:
            Rd = self._arrhenius(Rd25, Ear, Tleaf)

        return Rd

    def check_supplied_args(self, Jmax, Vcmax, Rd, Jmax25, Vcmax25, Rd25):
        """ Check the user supplied arguments, either they supply the values
        at 25 deg C, or the supply Jmax and Vcmax at the measurement temp. It
        is of course possible they accidentally supply both or a random
        combination, so raise an exception if so

        Parameters
        ----------
        Jmax : float
            potential rate of electron transport at measurement temperature
            [deg K]
        Vcmax : float
            max rate of rubisco activity at measurement temperature [deg K]
        Rd : float
            Day "light" respiration [umol m-2 time unit-1]
        Jmax25 : float
            potential rate of electron transport at 25 deg or 298 K
        Vcmax25 : float
            max rate of rubisco activity at 25 deg or 298 K
        Rd25 : float
            Estimate of respiration rate at the reference temperature 25 deg C
             or 298 K [deg K]

        Returns
        -------
        Nothing
        """
        at_25C = all(x is not None for x in (Jmax25, Vcmax25, Rd25)) and \
                    all(x is None for x in (Jmax, Vcmax, Rd))

        at_measured = all(x is not None for x in (Jmax, Vcmax, Rd)) and \
                        all(x is None for x in (Jmax25, Vcmax25, Rd25))

        if not (at_25C or at_measured):
            raise ValueError("Supplied arguments are a mess!")

    def quadratic(self, a=None, b=None, c=None, large=False):
        """ minimilist quadratic solution as root for J solution should always
        be positive, so I have excluded other quadratic solution steps. I am
        only returning the smallest of the two roots

        Parameters:
        ----------
        a : float
            co-efficient
        b : float
            co-efficient
        c : float
            co-efficient

        Returns:
        -------
        val : float
            positive root
        """
        d = b**2.0 - 4.0 * a * c # discriminant
        if d < 0.0:
            raise ValueError('imaginary root found')
        #root1 = np.where(d>0.0, (-b - np.sqrt(d)) / (2.0 * a), d)
        #root2 = np.where(d>0.0, (-b + np.sqrt(d)) / (2.0 * a), d)

        if large:
            if math.isclose(a, 0.0) and b > 0.0:
                root = -c / b
            elif math.isclose(a, 0.0) and math.isclose(b, 0.0):
                root = 0.0
                if c != 0.0:
                    raise ValueError('Cant solve quadratic')
            else:
                root = (-b + np.sqrt(d)) / (2.0 * a)
        else:
            if math.isclose(a, 0.0) and b > 0.0:
                root = -c / b
            elif math.isclose(a, 0.0) and math.isclose(b, 0.0):
                root == 0.0
                if c != 0.0:
                    raise ValueError('Cant solve quadratic')
            else:
                root = (-b - np.sqrt(d)) / (2.0 * a)

        return root
