
def apsim_co2_response(temp, co2, pathway="c3"):
    """
    Calculation of the CO2 modification on rue

    Reference
    ---------
    * Reyenga, Howden, Meinke, Mckeon (1999), Modelling global change impact on
      wheat cropping in south-east Queensland, Australia. Enivironmental
      Modelling & Software 14:297-306

    """
    base_co2 = 350.

    if pathway == "c3":
        # co2 compensation point (ppm)
        gamma_star  = max( (163.0 - temp) / ( 5.0 - 0.1 * temp), 0.0)

        arg1 = (co2 - gamma_star) * (base_co2 + 2.0 *  gamma_star)
        arg2 = (co2 + 2.0 * gamma_star) * (base_co2 - gamma_star)
        co2_mod = arg1 / arg2
    else:
        # Mark Howden, personal communication
        co2_mod = (0.000143 * co2 + 0.95)
    return co2_mod
