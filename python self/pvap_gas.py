# pvap_gas.py
from pvap_nh3 import pvap_nh3

def pvap_gas(gas_name: str, t_layer: float, p_layer: float = None) -> float:
    """
    Calculate vapor pressure (dyne/cm^2) for a condensate species.
    Only NH3 is currently implemented.

    Parameters
    ----------
    gas_name : str
        Chemical formula, e.g., "NH3".
    t_layer : float
        Temperature (K).
    p_layer : float, optional
        Layer pressure (dyne/cm^2), only needed for some species (not NH3).

    Returns
    -------
    float
        Vapor pressure in dyne/cm^2.

    Raises
    ------
    ValueError
        If gas_name is not recognized.
    """
    if gas_name == "NH3":
        return pvap_nh3(t_layer)
    else:
        # Match Fortran behavior: print and stop
        print(f"stop in pvap_gas(), bad gas_name = {gas_name}")
        raise ValueError(f"pvap_gas(): species '{gas_name}' not implemented")
