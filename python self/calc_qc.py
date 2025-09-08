# calc_qc.py
import math
import numpy as np
from advdiff import AdvDiff
from find_root import find_root
from vfall import vfall
from globals import R_GAS, PI

def calc_qc(q_below, q_vs, mixl, dz, grav, mw_atmos, mw_cloud, rho_p, rainf,
            w_convect, t_layer, p_layer, sig_layer, delta=1e-8):
    """
    Calculate condensate properties in a single layer (Python version of Fortran calc_qc).

    Parameters
    ----------
    q_below : float
        Total mixing ratio in the layer below [g/g]
    q_vs : float
        Saturation mixing ratio [g/g]
    mixl : float
        Convective mixing length [cm]
    dz : float
        Layer thickness [cm]
    grav : float
        Gravity [cm/s^2]
    mw_atmos : float
        Atmospheric mean molecular weight [g/mol]
    mw_cloud : float
        Condensate molecular weight [g/mol]
    rho_p : float
        Condensate particle density [g/cm^3]
    rainf : float
        Rain efficiency factor
    w_convect : float
        Convective velocity [cm/s]
    t_layer : float
        Layer temperature [K]
    p_layer : float
        Layer pressure [dyne/cm^2]
    sig_layer : float
        Geometric std. dev. of particle distribution
    delta : float
        Convergence tolerance

    Returns
    -------
    qc_layer : float
        Condensate mixing ratio [g/g]
    qt_layer : float
        Total mixing ratio [g/g]
    rg_layer : float
        Geometric mean particle radius [cm]
    reff_layer : float
        Effective particle radius [cm]
    ndz_layer : float
        Column droplet number concentration [cm^-2]
    qt_top : float
        Total mixing ratio at top of layer [g/g]
    status_r : int
        Root-finder status for particle radius
    status_q : int
        Root-finder status for total mixing ratio
    """

    # Initialize outputs
    qc_layer = qt_layer = rg_layer = reff_layer = ndz_layer = qt_top = 0.0
    status_r = status_q = 0

    # Cloud-free layer
    if q_below < q_vs:
        qt_layer = qt_top = q_below
        return qc_layer, qt_layer, rg_layer, reff_layer, ndz_layer, qt_top, status_r, status_q

    # --- Step 1: Advective-diffusive solution for total mixing ratio ---
    adv = AdvDiff(ad_qbelow=q_below, ad_qvs=q_vs, ad_mixl=mixl, ad_dz=dz, ad_rainf=rainf)
    f_advdiff = lambda qt: adv.compute(qt)[0]

    qlo = q_vs
    qhi = max(q_below, qlo * 2)
    delta_q = q_below / 1000.0

    qt_top, status_q = find_root(f_advdiff, 0.0, qlo, qhi, delta_q)

    qt_layer = 0.5 * (q_below + qt_top)
    qc_layer = max(0.0, qt_layer - q_vs)

    # --- Step 2: Solve for particle radius using vfall balance ---
    rlo = 1e-10
    rhi = 10.0
    delta_v = w_convect / 1000.0

    f_vfall = lambda r: vfall(r, qc_layer, dz, rainf)
    rw_layer, status_r = find_root(f_vfall, w_convect, rlo, rhi, delta_v)

    # --- Step 3: Compute particle properties ---
    lnsig2 = 0.5 * math.log(sig_layer)**2
    sig_alpha = max(1.1, sig_layer)

    if rainf > 1.0:
        alpha = math.log(vfall(rw_layer*sig_alpha, qc_layer, dz, rainf)/w_convect) / math.log(sig_alpha)
    else:
        alpha = math.log(w_convect / vfall(rw_layer/sig_alpha, qc_layer, dz, rainf)) / math.log(sig_alpha)

    rg_layer = rw_layer * (rainf**(1.0/alpha)) * math.exp(-(alpha+6.0)*lnsig2)
    reff_layer = rg_layer * math.exp(5.0*lnsig2)

    # Atmospheric density
    rho_atmos = p_layer / (R_GAS / mw_atmos * t_layer)

    # Column droplet number concentration
    ndz_layer = 3.0 * rho_atmos * qc_layer * dz / (4.0 * PI * rho_p * rg_layer**3) * math.exp(-9.0*lnsig2)

    return qc_layer, qt_layer, rg_layer, reff_layer, ndz_layer, qt_top, status_r, status_q
