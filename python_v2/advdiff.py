import math

def advdiff(qt: float, ad_qbelow: float, ad_qvs: float,
            ad_mixl: float, ad_dz: float, ad_rainf: float) -> tuple[float, float]:
    """
    Calculate divergence from advective-diffusive balance for condensate in a model layer.

    Parameters
    ----------
    qt : float
        Total mixing ratio of condensate + vapor (g/g).
    ad_qbelow : float
        Total mixing ratio of vapor in underlying layer (g/g).
    ad_qvs : float
        Saturation mixing ratio (g/g).
    ad_mixl : float
        Convective mixing length (cm).
    ad_dz : float
        Layer thickness (cm).
    ad_rainf : float
        Rain efficiency factor.

    Returns
    -------
    advdiff_val : float
        Difference from advective-diffusive balance.
    ad_qc : float
        Mixing ratio of condensed condensate (g/g).
    """
    # All vapor in excess of saturation condenses
    ad_qc = max(0.0, qt - ad_qvs)

    # Difference from advective-diffusive balance
    advdiff_val = ad_qbelow * math.exp(-ad_rainf * ad_qc * ad_dz / (qt * ad_mixl)) - qt

    return advdiff_val, ad_qc
