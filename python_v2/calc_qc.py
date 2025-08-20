# calc_qc.py
import math
import numpy as np

from pvap_nh3 import pvap_nh3
from advdiff import advdiff
from vfall import vfall, VFallEnv
from find_root import find_root

# Universal gas constant in cgs (erg mol^-1 K^-1 = dyne·cm mol^-1 K^-1)
R_GAS = 8.314462618e7
PI = math.pi


def calc_qc(
    gas_name: str,
    rainf_layer: float,
    rho_p: float,
    mw_cloud: float,
    q_below: float,
    supsat: float,
    w_convect: float,
    mixl: float,
    dz_layer: float,
    grav: float,
    mw_atmos: float,
    mfp: float,
    visc: float,
    t_layer: float,
    p_layer: float,
    sig_layer: float,
):
    """
    Python translation of the Fortran subroutine calc_qc (Ackerman & Marley style).
    All units are cgs.

    Parameters
    ----------
    gas_name : str
        Name of condensing vapor (only 'NH3' supported here).
    rainf_layer : float
        Rain factor for layer.
    rho_p : float
        Condensate density (g/cm^3).
    mw_cloud : float
        Molecular weight of condensing vapor (g/mol).
    q_below : float
        Total (vapor+condensate) mixing ratio below layer (g/g).
    supsat : float
        Fractional supersaturation persisting after condensation.
    w_convect : float
        Convective velocity scale (cm/s).
    mixl : float
        Convective mixing length scale (cm).
    dz_layer : float
        Layer thickness (cm).
    grav : float
        Gravitational acceleration (cm/s^2).
    mw_atmos : float
        Molecular weight of atmosphere (g/mol).
    mfp : float
        Mean free path (cm).
    visc : float
        Dynamic viscosity (dyne·s/cm^2).
    t_layer : float
        Temperature (K).
    p_layer : float
        Pressure (dyne/cm^2).
    sig_layer : float
        Geometric std. dev. of lognormal size distribution.

    Returns
    -------
    qc_layer : float
    qt_layer : float
    rg_layer : float
    reff_layer : float
    ndz_layer : float
    qt_top : float
    status_r : int   # status from rw solve (vfall)
    status_q : int   # status from qt_top solve (advdiff)
    """

    # --- Saturation vapor pressure (dyne/cm^2) ---
    if gas_name != "NH3":
        raise NotImplementedError("This calc_qc.py currently supports only NH3.")
    pvap = pvap_nh3(t_layer)

    # --- Saturation factor and atmospheric density ---
    fs = supsat + 1.0
    rho_atmos = p_layer / ((R_GAS / mw_atmos) * t_layer)

    # --- Saturated mass mixing ratio of vapor (g/g) ---
    # qvs = fs * pvap / ((R_GAS/mw_cloud)*t_layer) / rho_atmos
    qvs = fs * pvap / ((R_GAS / mw_cloud) * t_layer) / rho_atmos

    # Defaults (for cloud-free path) and status flags
    status_r = 0
    status_q = 0

    # ----------------------------
    # Cloud-free layer
    # ----------------------------
    if q_below < qvs:
        qt_layer = q_below
        qt_top = q_below
        qc_layer = 0.0
        rg_layer = 0.0
        reff_layer = 0.0
        ndz_layer = 0.0
        return (
            qc_layer,
            qt_layer,
            rg_layer,
            reff_layer,
            ndz_layer,
            qt_top,
            status_r,
            status_q,
        )

    # ----------------------------
    # Cloudy layer
    # 1) Solve for qt_top using advdiff balance
    # ----------------------------
    qhi = q_below
    qlo = qhi / 1e3
    delta_q = q_below / 1000.0

    # Residual for root finder: advdiff(qt, ...) - 0
    def adv_residual(qt):
        val, _ = advdiff(qt, q_below, qvs, mixl, dz_layer, rainf_layer)
        return val

    qt_top, status_q = find_root(
        adv_residual, target=0.0, xlo=qlo, xhi=qhi, tol=delta_q
    )

    # Trapezoid-rule layer mean (Fortran comment: could integrate exponential)
    qt_layer = 0.5 * (q_below + qt_top)

    # Diagnose condensate mixing ratio
    qc_layer = max(0.0, qt_layer - qvs)

    # ----------------------------
    # 2) Solve for rw_layer s.t. vfall(rw) = w_convect
    # ----------------------------
    rlo = 1.0e-10
    rhi = 1.0e1
    delta_v = w_convect / 1000.0

    def vfall_residual(r):
        # vfall must accept (r, grav, mw_atmos, mfp, visc, t_layer, p_layer, rho_p)
        vf_env = VFallEnv(grav, mw_atmos, mfp, visc, t_layer, p_layer, rho_p)
        return vfall(r, vf_env)

    rw_layer, status_r = find_root(
        vfall_residual, target=w_convect, xlo=rlo, xhi=rhi, tol=delta_v
    )

    # ----------------------------
    # 3) Size distribution and optical properties
    # ----------------------------
    lnsig2 = 0.5 * (math.log(sig_layer) ** 2)

    # sigma floor for alpha calculation
    sig_alpha = max(1.1, sig_layer)

    # Compute exponent alpha in vfall ∝ r^alpha near rw
    if rainf_layer > 1.0:
        # bulk precip at r > rw: exponent between rw and rw*sig
        alpha = math.log(
            vfall_residual(rw_layer * sig_alpha) / w_convect
        ) / math.log(sig_alpha)
    else:
        # bulk precip at r < rw: exponent between rw/sig and rw
        alpha = math.log(
            w_convect / vfall_residual(rw_layer / sig_alpha)
        ) / math.log(sig_alpha)

    # Geometric mean radius (cm)
    rg_layer = (rainf_layer ** (1.0 / alpha)) * rw_layer * math.exp(
        -(alpha + 6.0) * lnsig2
    )

    # Effective radius (cm)
    reff_layer = rg_layer * math.exp(5.0 * lnsig2)

    # Column number concentration (cm^-2)
    if rg_layer > 0.0:
        ndz_layer = (
            3.0
            * rho_atmos
            * qc_layer
            * dz_layer
            / (4.0 * PI * rho_p * (rg_layer ** 3))
            * math.exp(-9.0 * lnsig2)
        )
    else:
        ndz_layer = 0.0  # avoid division by zero if degenerate

    return (
        qc_layer,
        qt_layer,
        rg_layer,
        reff_layer,
        ndz_layer,
        qt_top,
        status_r,
        status_q,
    )
