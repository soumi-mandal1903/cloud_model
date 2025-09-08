# eddysed.py
import numpy as np
from layer import layer
from globals import STEFBOLTZ

def eddysed(
    p, t, dz, grav, mw_atmos, teff,
    mixl_min, cloudf_min, mw_cloud,
    rainf, rho_p, supsat, sig_layer,
    q_init, nsub_max=8, delta=1e-8
):
    """
    Python port of the Fortran `eddysed` main routine.
    Computes cloud condensate properties through atmosphere.

    Parameters
    ----------
    p : ndarray
        Pressure profile [dyne/cm^2]
    t : ndarray
        Temperature profile [K]
    dz : ndarray
        Layer thicknesses [cm]
    grav : float
        Gravitational acceleration [cm/s^2]
    mw_atmos : float
        Atmospheric mean molecular weight [g/mol]
    teff : float
        Effective temperature [K]
    mixl_min : float
        Minimum convective mixing length [cm]
    cloudf_min : float
        Minimum cloud fraction
    mw_cloud : float
        Cloud condensate molecular weight [g/mol]
    rainf : float
        Rainfall efficiency
    rho_p : float
        Particle density [g/cm^3]
    supsat : float
        Supersaturation (unused, typically 0)
    sig_layer : float
        Geometric standard deviation of lognormal size distribution
    q_init : float
        Initial mixing ratio at base [g/g]
    nsub_max : int
        Max number of sublayers for optical depth iteration
    delta : float
        Root-finding tolerance

    Returns
    -------
    qc : ndarray
        Condensate mixing ratio [g/g]
    qt : ndarray
        Total mixing ratio [g/g]
    ndz : ndarray
        Droplet number density profile
    rg : ndarray
        Geometric mean radius [cm]
    reff : ndarray
        Effective radius [cm]
    cloudf : ndarray
        Cloud fraction profile
    """

    nz = len(p)

    qc = np.zeros(nz)
    qt = np.zeros(nz)
    ndz = np.zeros(nz)
    rg = np.zeros(nz)
    reff = np.zeros(nz)
    cloudf = np.zeros(nz)

    # convective heat flux profile
    chf = np.full(nz, STEFBOLTZ * teff**4)

    # initial vapor-only mixing ratio
    q_below = q_init

    for k in range(nz):
        # Pressure midpoint
        p_layer = p[k]
        t_layer = t[k]

        # Call layer solver (already uses corrected find_root inside)
        qc_layer, qt_layer, rg_layer, reff_layer, ndz_layer, q_below, _, _ = layer(
            nsub_max=nsub_max,
            gas_name="NH3",          # or pass dynamically
            grav=grav,
            mw_atmos=mw_atmos,
            kz_min=mixl_min,
            cloudf_min=cloudf_min,
            mw_cloud=mw_cloud,
            rainf=rainf,
            rho_p=rho_p,
            supsat=supsat,
            sig_layer=sig_layer,
            cloudf=1.0,              # overwritten below
            q_below=q_below,
            t_layer=t_layer,
            p_layer=p_layer,
            kz=mixl_min,
            chf=chf[k],
            t_top=t_layer,
            t_bot=t_layer,
            p_top=p_layer,
            p_bot=p_layer + dz[k]*grav/p_layer,  # hydrostatic guess
        )

        qc[k] = qc_layer
        qt[k] = qt_layer
        rg[k] = rg_layer
        reff[k] = reff_layer
        ndz[k] = ndz_layer
        cloudf[k] = max(cloudf_min, min(1.0, 1.0 - qc_layer))

    return qc, qt, ndz, rg, reff, cloudf
