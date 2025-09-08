# qvs_below.py
import numpy as np
from pvap_gas import pvap_gas

# Globals (act like the COMMON block in Fortran)
qv_dtdlnp = None   # temperature lapse rate wrt log pressure
qv_p = None        # reference pressure (dyne/cm^2)
qv_t = None        # reference temperature (K)
qv_factor = None   # molecular weight ratio factor
qv_gas_name = None # string gas name ("NH3", "H2O", ...)

def qvs_below(p_test: float) -> float:
    """
    Compute saturation mixing ratio for a gas below the model domain.

    Parameters
    ----------
    p_test : float
        Test pressure [dyne/cm^2]

    Returns
    -------
    qvs : float
        Saturation mixing ratio [g/g]
    """
    global qv_dtdlnp, qv_p, qv_t, qv_factor, qv_gas_name

    if None in (qv_dtdlnp, qv_p, qv_t, qv_factor, qv_gas_name):
        raise RuntimeError("qvs_below globals not initialized")

    # Extrapolate T at this pressure
    t_test = qv_t + np.log(qv_p / p_test) * qv_dtdlnp

    # Compute vapor pressure [dyne/cm^2]
    pvap_test = pvap_gas(qv_gas_name, t_test, p_test)

    # Saturation mixing ratio
    return qv_factor * pvap_test / p_test


def init_qvs_below(dtdlnp, p_ref, t_ref, factor, gas_name):
    """
    Initialize qvs_below globals.

    Parameters
    ----------
    dtdlnp : float
        dT/dlnP (K)
    p_ref : float
        Reference pressure [dyne/cm^2]
    t_ref : float
        Reference temperature [K]
    factor : float
        Molecular weight ratio (M_atmos / M_gas)
    gas_name : str
        Gas species name
    """
    global qv_dtdlnp, qv_p, qv_t, qv_factor, qv_gas_name
    qv_dtdlnp = dtdlnp
    qv_p = p_ref
    qv_t = t_ref
    qv_factor = factor
    qv_gas_name = gas_name
