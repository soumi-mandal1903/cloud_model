# eddysed.py
"""
Python port of eddysed.f
Balance between eddy diffusion and sedimentation for condensates.

Requires:
    - pvap_gas(name, T, P) : saturation vapor pressure [dyne/cm^2]
    - layer(...)           : layer calculation routine (ported already)
    - find_rootl(func, target, p_lo, p_hi, delta_q) : root finder

Author: Ported from Ackerman & Marley (2001) Fortran implementation
"""

import numpy as np
from layer import layer
from pvap_nh3 import pvap_nh3  # currently only NH3 supported
from scipy.optimize import root_scalar
from qvs_below import qvs_below


# def qvs_below(p, q_target, dtdlnp, p_ref, t_ref, factor, gas_name):
#     """
#     Function for root finding: qvs(p) - q_target
#     """
#     T = t_ref + np.log(p_ref / p) * dtdlnp
#     pvap = pvap_gas(gas_name, T, p)
#     qvs = factor * pvap / p
#     return qvs - q_target


def eddysed(
    grav,
    teff,
    kz_min,
    cloudf_min,
    nsub_max,
    supsat,
    mw_atmos,
    do_virtual,
    nz,
    z,
    z_top,
    p,
    p_top,
    t,
    t_top,
    chf,
    ngas,
    gas_name,
    gas_mmr,
    gas_mw,
    rho_p,
    sig,
    rainf,
):
    """
    Main driver routine: calculate condensate size and concentration
    """
    # Allocate outputs
    kz = np.zeros(nz)
    qt = np.zeros((nz, ngas))
    qc = np.zeros((nz, ngas))
    ndz = np.zeros((nz, ngas))
    rg = np.zeros((nz, ngas))
    reff = np.zeros((nz, ngas))
    cloudf = np.zeros(nz)

    qc_path = np.zeros(ngas)
    q_below = np.array(gas_mmr, copy=True)

    # Determine bottom/top index ordering
    if z[1] > z[0]:
        ibot, itop, incr = 0, nz - 1, 1
    else:
        ibot, itop, incr = nz - 1, 0, -1

    # Loop over condensable gases
    for igas in range(ngas):
        t_bot = t_top[ibot - incr]
        p_bot = p_top[ibot - incr]

        qc_path[igas] = 0.0

        # --- Adjust mixing ratio below cloud base ---
        if do_virtual:
            qvs_factor = (supsat + 1.0) * gas_mw[igas] / mw_atmos
            pvap = pvap_nh3(t_bot)
            qvs = qvs_factor * pvap / p_bot

            if qvs <= q_below[igas]:
                # find cloud base pressure
                dtdlnp = (t_top[ibot] - t_bot) / np.log(p_bot / p_top[ibot])

                def f_root(p_test):
                    return qvs_below(
                        p_test,
                        q_below[igas],
                        dtdlnp,
                        p_bot,
                        t_bot,
                        qvs_factor,
                        gas_name[igas],
                    )

                sol = root_scalar(f_root, bracket=[p_bot, p_bot * 1e2], method="bisect")
                if not sol.converged:
                    raise RuntimeError(
                        f"Unable to find cloud base for {gas_name[igas]}"
                    )
                p_base = sol.root
                t_base = t_bot + np.log(p_bot / p_base) * dtdlnp

                # virtual layer midpoint
                p_layer = 0.5 * (p_bot + p_base)
                t_layer = t_bot + np.log(p_bot / p_layer) * dtdlnp

                # call layer() for virtual layer
                qc_layer, qt_layer, rg_layer, reff_layer, ndz_layer, kz_layer, report_status_r, report_status_q = layer(
                    nsub_max,
                    gas_name[igas],
                    grav,
                    mw_atmos,
                    kz_min,
                    cloudf_min,
                    gas_mw[igas],
                    rainf[ibot, igas],
                    rho_p[igas],
                    supsat,
                    sig[ibot, igas],
                    cloudf[ibot],
                    q_below[igas],
                    t_layer,
                    p_layer,
                    kz[ibot],
                    chf[ibot],
                    t_bot,
                    t_base,
                    p_bot,
                    p_base,
                )

                print("")
                print(f"eddysed(): condensing gas = {gas_name[igas]}")
                print(f" cloud base at p, T = {p_base/1e6:.3f} bar, {t_base:.2f} K")
                tau_below = (
                    1.5
                    * qc_layer
                    * (p_base - p_bot)
                    / grav
                    / (rho_p[igas] * reff_layer)
                )
                print(f" optical depth below domain = {tau_below:.3e}")
                print("")

        # --- Loop over atmospheric layers ---
        iz = ibot
        while True:
            qc_layer, qt_layer, rg_layer, reff_layer, ndz_layer, kz_layer, report_status_r, report_status_q = layer(
                nsub_max,
                gas_name[igas],
                grav,
                mw_atmos,
                kz_min,
                cloudf_min,
                gas_mw[igas],
                rainf[iz, igas],
                rho_p[igas],
                supsat,
                sig[iz, igas],
                cloudf[iz],
                q_below[igas],
                t[iz],
                p[iz],
                kz[iz],
                chf[iz],
                t_top[iz],
                t_top[iz - incr],
                p_top[iz],
                p_top[iz - incr],
            )

            qc[iz, igas] = qc_layer
            qt[iz, igas] = qt_layer
            rg[iz, igas] = rg_layer
            reff[iz, igas] = reff_layer
            ndz[iz, igas] = ndz_layer
            kz[iz] = kz_layer

            qc_path[igas] += qc[iz, igas] * (p_top[iz - incr] - p_top[iz]) / grav

            if iz == itop:
                break
            iz += incr

        print("")
        print(f"eddysed(): condensing gas = {gas_name[igas]}")
        print(f" condensate path = {qc_path[igas]*1e4:.3f} g/m^2")
        print("")

    return kz, qt, qc, ndz, rg, reff, cloudf
