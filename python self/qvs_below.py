import numpy as np
from pvap_gas import pvap_gas

# Globals simulating the Fortran COMMON block
qv_dtdlnp = None   # dT/dlnP
qv_p = None        # reference pressure (dyne/cm^2)
qv_t = None        # reference temperature (K)
qv_factor = None   # molecular weight ratio (M_atmos / M_gas)
qv_gas_name = None # gas name string ("NH3", ...)

def qvs_below(p_test: float) -> float:
    """
    Compute saturation mixing ratio for a gas below the model domain.

    Parameters
    ----------
    p_test : float
        Test pressure (dyne/cm^2)

    Returns
    -------
    qvs : float
        Saturation mixing ratio (g/g)
    """
    global qv_dtdlnp, qv_p, qv_t, qv_factor, qv_gas_name

    if None in (qv_dtdlnp, qv_p, qv_t, qv_factor, qv_gas_name):
        raise RuntimeError("qvs_below globals not initialized")

    # Extrapolate temperature using lapse rate wrt log-pressure
    t_test = qv_t + np.log(qv_p / p_test) * qv_dtdlnp

    # Compute saturation vapor pressure
    pvap_test = pvap_gas(qv_gas_name, t_test, p_test)

    # Saturation mixing ratio (g/g)
    qvs = qv_factor * pvap_test / p_test
    return qvs


def init_qvs_below(dtdlnp: float, p_ref: float, t_ref: float, factor: float, gas_name: str):
    """
    Initialize global variables for qvs_below.

    Parameters
    ----------
    dtdlnp : float
        Temperature gradient w.r.t ln(pressure)
    p_ref : float
        Reference pressure (dyne/cm^2)
    t_ref : float
        Reference temperature (K)
    factor : float
        Molecular weight ratio (M_atmos / M_gas)
    gas_name : str
        Gas species name ("NH3", ...)
    """
    global qv_dtdlnp, qv_p, qv_t, qv_factor, qv_gas_name
    qv_dtdlnp = dtdlnp
    qv_p = p_ref
    qv_t = t_ref
    qv_factor = factor
    qv_gas_name = gas_name
