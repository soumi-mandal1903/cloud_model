# layer.py
import numpy as np
from advdiff import advdiff
from find_root import find_root

def layer(qt_guess, qvs, qbelow, mixl, dz, rainf, delta=1e-8):
    """
    Python port of the Fortran `layer` subroutine.

    Parameters
    ----------
    qt_guess : float
        Initial guess for total mixing ratio (g/g)
    qvs : float
        Saturation mixing ratio (g/g)
    qbelow : float
        Vapor mixing ratio in underlying layer (g/g)
    mixl : float
        Convective mixing length (cm)
    dz : float
        Layer thickness (cm)
    rainf : float
        Rain efficiency factor
    delta : float
        Convergence tolerance

    Returns
    -------
    qc : float
        Condensed condensate mixing ratio (g/g)
    qt : float
        Total mixing ratio (g/g)
    status : int
        Root-finding status flag
    """

    # Store state in globals for advdiff (matches Fortran COMMON block)
    from advdiff import set_advdiff_params
    set_advdiff_params(qbelow, qvs, mixl, dz, rainf)

    # Define wrapper for root finding
    f = lambda qt: advdiff(qt)

    # Bracket search range: must enclose solution
    qmin = max(qvs, 1e-12)
    qmax = max(qt_guess * 2.0, qmin * 2.0)

    # Call root solver
    qt_new, status = find_root(f, 0.0, qmin, qmax, delta)

    # Compute condensate mixing ratio
    qc = max(0.0, qt_new - qvs)

    return qc, qt_new, status
