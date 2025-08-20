import numpy as np
from calc_qc import calc_qc

# Globals
R_GAS = 8.314462618e7  # erg/mol/K = dyne*cm/mol/K
AVOGADRO = 6.02214076e23  # mol^-1
K_BOLTZ = R_GAS/AVOGADRO  # erg/K
PI = np.pi

def layer(nsub_max, gas_name, grav, mw_atmos, kz_min, cloudf_min, mw_cloud, rainf, rho_p, supsat, sig_layer, cloudf, q_below, t_layer, p_layer, kz, chf, t_top, t_bot, p_top, p_bot, report_status_r=True, report_status_q=True):
    """
    Calculate layer condensate properties by iterating on optical depth
    in one model layer (converging on optical depth over sublayers).

    Parameters
    ----------
    (See original Fortran docstring for input/output descriptions)

    Returns
    -------
    qc_layer, qt_layer, rg_layer, reff_layer, ndz_layer, q_below, report_status_r, report_status_q
    """

    # ---- Constants ----
    d_molecule = 2.827e-8  # H2 molecule diameter (cm)
    eps_k = 59.7           # Lennard-Jones parameter for viscosity (K)

    # ---- Atmospheric properties ----
    r_atmos = R_GAS / mw_atmos  # specific gas constant (erg/K/g)
    c_p = 7.0 / 2.0 * r_atmos  # specific heat (erg/K/g)
    dp_layer = p_bot - p_top
    dlnp = np.log(p_bot / p_top)

    dtdlnp = (t_top - t_bot) / dlnp
    lapse_ratio = (t_bot - t_top) / dlnp / (2.0 / 7.0 * t_layer)

    rho_atmos = p_layer / (r_atmos * t_layer)  # g/cm^3
    scale_h = r_atmos * t_layer / grav         # cm
    mixl = max(0.1, lapse_ratio) * scale_h     # convective mixing length
    scalef_kz = 1.0 / 3.0
    kz = scalef_kz * scale_h * (mixl / scale_h) ** (4.0 / 3.0) * (
        (r_atmos * chf) / (rho_atmos * c_p)
    ) ** (1.0 / 3.0)
    kz = max(kz, kz_min)

    w_convect = kz / mixl
    cloudf = cloudf_min + max(0.0, min(1.0, 1.0 - lapse_ratio)) * (1.0 - cloudf_min)

    n_atmos = p_layer / (K_BOLTZ * t_layer)
    mfp = 1.0 / (np.sqrt(2.0) * n_atmos * PI * d_molecule**2)
    visc = (
        5.0
        / 16.0
        * np.sqrt(PI * K_BOLTZ * t_layer * (mw_atmos / AVOGADRO))
        / (PI * d_molecule**2)
        / (1.22 * (t_layer / eps_k) ** (-0.16))
    )

    # ---- Iterative convergence loop ----
    converge = False
    nsub = 1
    opd_test = 0.0

    while not converge:
        qc_layer = 0.0
        qt_layer = 0.0
        ndz_layer = 0.0
        opd_layer = 0.0

        qt_bot_sub = q_below
        p_bot_sub = p_bot
        dp_sub = dp_layer / nsub

        for isub in range(nsub):
            qt_below = qt_bot_sub
            p_top_sub = p_bot_sub - dp_sub
            dz_sub = scale_h * np.log(p_bot_sub / p_top_sub)
            p_sub = 0.5 * (p_bot_sub + p_top_sub)
            t_sub = t_bot + np.log(p_bot / p_sub) * dtdlnp

            qc_sub, qt_sub, rg_sub, reff_sub, ndz_sub, qt_top, status_r, status_q = calc_qc(
                gas_name, rainf, rho_p, mw_cloud,
                qt_below, supsat, w_convect, mixl,
                dz_sub, grav, mw_atmos, mfp, visc,
                t_sub, p_sub, sig_layer
            )

            qc_layer += qc_sub * dp_sub / grav
            qt_layer += qt_sub * dp_sub / grav
            ndz_layer += ndz_sub

            if reff_sub > 0.0:
                opd_layer += 1.5 * qc_sub * dp_sub / grav / (rho_p * reff_sub)

            qt_bot_sub = qt_top
            p_bot_sub = p_top_sub

        # ---- Convergence check ----
        if nsub_max == 1:
            converge = True
        elif nsub == 1:
            opd_test = opd_layer
        elif opd_layer == 0.0 or nsub >= nsub_max:
            converge = True
        elif abs(1.0 - opd_test / opd_layer) <= 1e-2:
            converge = True
        else:
            opd_test = opd_layer

        nsub *= 2

    # ---- Error reporting ----
    if status_r != 0 and report_status_r:
        print(f"layer(): find_root(vfall) status = {status_r} for {gas_name} at p = {p_layer/1e6:.3e} bar")
        report_status_r = False

    if status_q != 0 and report_status_q:
        print(f"layer(): find_root(advdiff) status = {status_q} for {gas_name} at p = {p_layer/1e6:.3e} bar")
        report_status_q = False

    # ---- Update for next layer ----
    q_below = qt_top

    # ---- Layer averages ----
    if opd_layer > 0.0:
        reff_layer = 1.5 * qc_layer / (rho_p * opd_layer)
        lnsig2 = 0.5 * np.log(sig_layer) ** 2
        rg_layer = reff_layer * np.exp(-5.0 * lnsig2)
    else:
        reff_layer = 0.0
        rg_layer = 0.0

    qc_layer = qc_layer * grav / dp_layer
    qt_layer = qt_layer * grav / dp_layer

    return (
        qc_layer,
        qt_layer,
        rg_layer,
        reff_layer,
        ndz_layer,
        q_below,
        report_status_r,
        report_status_q,
    )
