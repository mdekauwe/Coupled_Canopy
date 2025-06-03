
def apsim_co2_response(temp, co2, pathway="c3"):
    """
    Calculation of the CO2 modification on rue, based on Gifford and Morison
    (1993) who suggested that the increase in RUE of a wheat crop canopy under
    nhanced CO2 could be explained by the simple theoretical increase in
    light-limited photosynthesis when a small correction for the change in
    respiratory efficiency was made.

    Reference
    ---------
    * Reyenga, Howden, Meinke, Mckeon (1999), Modelling global change impact on
      wheat cropping in south-east Queensland, Australia. Enivironmental
      Modelling & Software 14:297-306
    * Gifford, R.M., Morison, J.I.L., 1993. Crop response to the global increase
      in atmospheric carbon dioxide concentration. In: Inter- national Crop
      Science I. CSSA, USA, pp. 325–331.
    * Bykov, O.D., Koshkin, V.A., Catsky, J., 1981. Carbon dioxide compensation
      concentration of C3 and C4 plants: dependence on temperature.
      Photosynthetica 15 (1), 114–121.

    """
    base_co2 = 350.

    if pathway == "c3":
        # co2 compensation point (ppm) based on Bykov et al. 1981
        gamma_star  = max( (163.0 - temp) / ( 5.0 - 0.1 * temp), 0.0)

        arg1 = (co2 - gamma_star) * (base_co2 + 2.0 *  gamma_star)
        arg2 = (co2 + 2.0 * gamma_star) * (base_co2 - gamma_star)
        co2_mod = arg1 / arg2
    else:
        # Mark Howden, personal communication
        co2_mod = (0.000143 * co2 + 0.95)
    return co2_mod
