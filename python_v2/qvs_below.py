# qvs_below.py
import numpy as np
from pvap_nh3 import pvap_nh3  # or general pvap_gas function if you later expand

# Global "common block" variables (to be set by initialization before use)
qv_dtdlnp = None
qv_p = None
qv_t = None
qv_factor = None
qv_gas_name = None


def init_qvs_below(dtdlnp, p_ref, t_ref, factor, gas_name):
    """
    Initialize the 'common block' values for qvs_below.

    Parameters
    ----------
    dtdlnp : float
        Temperature lapse rate wrt log(pressure) (K per ln(p)).
    p_ref : float
        Reference pressure (dyne/cm^2).
    t_ref : float
        Reference temperature (K).
    factor : float
        Scaling factor (typically epsilon = Mw/Ma).
    gas_name : str
        Name of gas (e.g. "NH3").
    """
    global qv_dtdlnp, qv_p, qv_t, qv_factor, qv_gas_name
    qv_dtdlnp = dtdlnp
    qv_p = p_ref
    qv_t = t_ref
    qv_factor = factor
    qv_gas_name = gas_name


def qvs_below(p_test):
    """
    Compute saturation mixing ratio for a gas below the model domain,
    extrapolated using the local temperature lapse rate.

    Parameters
    ----------
    p_test : float
        Test pressure (dyne/cm^2).

    Returns
    -------
    qvs : float
        Saturation mixing ratio (g/g).
    """
    global qv_dtdlnp, qv_p, qv_t, qv_factor, qv_gas_name

    if None in (qv_dtdlnp, qv_p, qv_t, qv_factor, qv_gas_name):
        raise RuntimeError("qvs_below called before initialization with init_qvs_below")

    # Extrapolate temperature to p_test
    t_test = qv_t + np.log(qv_p / p_test) * qv_dtdlnp

    # Call vapor pressure function
    if qv_gas_name.upper() == "NH3":
        pvap_test = pvap_nh3(t_test)
    else:
        raise NotImplementedError(f"pvap_gas not available for {qv_gas_name}")

    # Saturation mixing ratio
    qvs = qv_factor * pvap_test / p_test
    return qvs
