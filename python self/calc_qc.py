# calc_qc.py

import numpy as np
from advdiff import advdiff
from find_root import find_root
from find_rootl import find_rootl
from vfall import vfall   # your Python translation of vfall.f

def calc_qc(q_below, q_vs, mixl, dz, rainf, rmin, rmax, delta=1e-8):
    """
    Calculate condensate mixing ratio in a model layer (Ackerman & Marley 2001).

    Parameters
    ----------
    q_below : float
        Total mixing ratio of vapor in the underlying layer (g/g).
    q_vs : float
        Saturation mixing ratio (g/g).
    mixl : float
        Convective mixing length (cm).
    dz : float
        Layer thickness (cm).
    rainf : float
        Rain efficiency factor.
    rmin, rmax : float
        Minimum and maximum trial radii (cm).
    delta : float
        Convergence tolerance.

    Returns
    -------
    qc : float
        Condensate mixing ratio (g/g).
    qt : float
        Total mixing ratio (g/g).
    rg : float
        Geometric mean particle radius (cm).
    reff : float
        Effective particle radius (cm).
    ndz : float
        Number of condensate particles per gram of atmosphere (1/g).
    kz : float
        Eddy diffusion coefficient (cm^2/s).
    status : int
        0 = success
        1 = max iterations reached
       -1 = no convergence
    """

    # --- Step 1: Solve advdiff balance for qt ---
    def f_advdiff(qt):
        return advdiff(qt, q_below, q_vs, mixl, dz, rainf)

    qt_guess_low = q_vs
    qt_guess_high = max(q_below, q_vs * 10)

    qt, status_qt = find_root(f_advdiff, 0.0, qt_guess_low, qt_guess_high, delta)

    if status_qt != 0:
        # fall back or warn
        return 0.0, qt, 0.0, 0.0, 0.0, 0.0, status_qt

    # condensate mixing ratio
    qc = max(0.0, qt - q_vs)

    # --- Step 2: Solve for particle radius using vfall balance ---
    def f_vfall(r):
        return vfall(r, qc, dz, rainf)

    r, status_r = find_rootl(f_vfall, 1.0, rmin, rmax, delta)

    if status_r != 0:
        return qc, qt, 0.0, 0.0, 0.0, 0.0, status_r

    # --- Step 3: Compute particle properties ---
    rg = r   # placeholder: geometric mean radius
    reff = rg * 1.5  # placeholder: effective radius scaling
    ndz = qc / (4/3 * np.pi * rg**3 * 1.0)  # assuming density=1 g/cm^3
    kz = mixl**2 / dz  # crude eddy diffusion estimate

    return qc, qt, rg, reff, ndz, kz, 0 # success                   
def compute_saturation_mixing_ratio(T, P, mw_cloud, mw_atmos):

    for iz in range(nz):
        # Initial guess for total mixing ratio
        if iz == 0:
            q_below = q_init
        else:
            q_below = qt[iz - 1]

        # Saturation mixing ratio at current layer
        q_vs = max(0.0, compute_saturation_mixing_ratio(t[iz], p[iz], mw_cloud, mw_atmos))

        # Convective mixing length
        mixl = max(mixl_min, compute_mixing_length(t, dz, iz))

        # Calculate condensate properties in the layer
        qc_layer, qt_layer, rg_layer, reff_layer, ndz_layer, kz_layer, status = calc_qc(
            q_below, q_vs, mixl, dz[iz], rainf, rmin=1e-4, rmax=1e2, delta=delta
        )

        if status != 0:
            print(f"Warning: calc_qc did not converge at layer {iz}, status {status}")

        qc[iz] = qc_layer
        qt[iz] = qt_layer
        rg[iz] = rg_layer
        reff[iz] = reff_layer
        ndz[iz] = ndz_layer
        cloudf[iz] = max(cloudf_min, min(1.0, qc_layer / (qc_layer + 1e-10)))  # crude cloud fraction       