import math

def pvap_nh3(t: float) -> float:
    """
    Calculate saturation vapor pressure of NH3 (ammonia).

    Parameters
    ----------
    t : float
        Temperature in Kelvin.

    Returns
    -------
    pvap : float
        Saturation vapor pressure in dyne/cm^2.
    """
    pvap = math.exp(-86596.0 / t**2 - 2161.0 / t + 10.53)
    pvap *= 1e6  # convert from bars to dyne/cm^2
    return pvap
