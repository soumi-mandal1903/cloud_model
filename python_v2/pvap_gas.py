# pvap_gas.py
from pvap_nh3 import pvap_nh3


def pvap_gas(gas_name: str, t_layer: float, p_layer: float = None) -> float:
    """
    Calculate vapor pressure (dyne/cm^2) for a given condensate species.

    Right now only NH3 is implemented (others raise NotImplementedError).

    Parameters
    ----------
    gas_name : str
        Chemical formula, e.g. "NH3".
    t_layer : float
        Temperature (K).
    p_layer : float, optional
        Pressure (dyne/cm^2), needed only for Mg2SiO4 (not implemented yet).

    Returns
    -------
    pvap : float
        Vapor pressure in dyne/cm^2
    """
    if gas_name == "NH3":
        return pvap_nh3(t_layer)
    else:
        raise NotImplementedError(f"pvap_gas: species '{gas_name}' not implemented yet.")
