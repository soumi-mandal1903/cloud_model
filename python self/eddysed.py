import numpy as np
from pvap_gas import pvap_gas
from vfall import vfall, VFallEnv

R_GAS = 8.3143e7  # erg/mol/K
PI = np.pi
K_BOLTZ = 1.380649e-16  # erg/K
AVOGADRO = 6.02214076e23

def eddysed(
    grav_in, teff, nz, z, p, t, dlnp, chf, lapse_ratio,
    gas_mmr, mw_atmos_in, sig, rainf
):
    """
    Eddy-sedimentation calculation for NH3 cloud layers.

    Parameters
    ----------
    grav_in : float
        Gravity (cm/s^2)
    teff : float
        Effective temperature (K)
    nz : int
        Number of atmospheric layers
    z, p, t : arrays
        Layer height (cm), pressure (dyne/cm^2), temperature (K)
    dlnp : array
        log-pressure difference across layers
    chf : array
        Convective heat flux (erg/s/cm^2)
    lapse_ratio : array
        Layer lapse ratio / adiabatic
    gas_mmr : float
        Mass mixing ratio of NH3 below cloud base
    mw_atmos_in : float
        Atmospheric molecular weight (g/mol)
    sig, rainf : arrays
        Lognormal sigma, precipitation efficiency per layer

    Returns
    -------
    qt, qc, ndz, rg, reff, lhf, kz : arrays
        Mixing ratios, droplet properties, latent heat flux, eddy diffusion
    """

    # Initialize arrays
    qt = np.zeros(nz)
    qc = np.zeros(nz)
    ndz = np.zeros(nz)
    rg = np.zeros(nz)
    reff = np.zeros(nz)
    lhf = np.zeros(nz)
    kz = np.zeros(nz)
    dz = np.zeros(nz)
    mixl = np.zeros(nz)

    grav = grav_in
    mw_atmos = mw_atmos_in

    # Molecular and physical constants
    d_molecule = 2.827e-8  # cm (H2)
    eps_k = 59.7  # Lennard-Jones (K)

    cloud_base = 0
    qbelow = gas_mmr
    fs = 1.0  # saturation factor

    # Determine loop direction (bottom-up)
    if z[1] > z[0]:
        ibot, itop, incr = 0, nz, 1
    else:
        ibot, itop, incr = nz-1, -1, -1

    for iz in range(ibot, itop, incr):
        # Atmospheric density and scale height
        r_atmos = R_GAS / mw_atmos
        c_p = 7./2. * r_atmos
        rho_atmos = p[iz] / (r_atmos * t[iz])
        scale_h = r_atmos * t[iz] / grav
        dz[iz] = scale_h * dlnp[iz]

        # Mixing length (convective)
        mixl[iz] = max(0.1, lapse_ratio[iz]) * scale_h

        # Number density, mean free path, viscosity
        n_atmos = p[iz] / (K_BOLTZ * t[iz])
        mfp = 1. / (np.sqrt(2.) * n_atmos * PI * d_molecule**2)
        visc = (5./16.) * np.sqrt(PI * K_BOLTZ * t[iz] * (mw_atmos / AVOGADRO)) \
               / (PI * d_molecule**2) / (1.22 * (t[iz]/eps_k)**(-0.16))

        # NH3 cloud properties
        mw_cloud = 17.0
        lh = 28.6e10  # erg/mol
        rho_p = 0.84  # g/cm^3
        pvap = pvap_gas("NH3", t[iz], p_layer=p[iz])
        qvs = fs * pvap / ( (R_GAS/mw_cloud) * t[iz] ) / rho_atmos

        # Vapor below
        if iz == ibot:
            qbelow_layer = gas_mmr
            zdiff = dz[iz]/2.
        else:
            qbelow_layer = qt[iz-incr]
            zdiff = z[iz] - z[iz-incr]

        if qbelow_layer < qvs:
            # Cloud-free
            qt[iz] = qbelow_layer
            qc[iz] = 0.0
            rg[iz] = 0.0
            reff[iz] = 0.0
            ndz[iz] = 0.0
        else:
            # Eddy-sedimentation: total mixing ratio = qvs + qc
            # Approximate qtotal = qvs + small increment (placeholder)
            qt[iz] = qvs + (qbelow_layer - qvs)
            qc[iz] = max(0.0, qt[iz] - qvs)

            # Update cloud base
            if cloud_base == 0 and qc[iz] > 0:
                cloud_base = iz

            # Convective velocity scale
            lh_sum = rainf[iz] * qc[iz] * lh / (mw_cloud * c_p / r_atmos)
            coeff_a = lh_sum
            coeff_b = - chf[iz] / (rho_atmos * c_p / r_atmos)

            if lh_sum == 0.0:
                w_convect = (-coeff_b)**(1./3.)
            else:
                cubic_arg = coeff_b**2 / 4. + coeff_a**3 / 27.
                w_convect = (-coeff_b/2. + np.sqrt(cubic_arg))**(1./3.) \
                           + (-coeff_b/2. - np.sqrt(cubic_arg))**(1./3.)

            lhf[iz] = w_convect * lh_sum * rho_atmos * c_p / r_atmos
            kz[iz] = w_convect * mixl[iz]

            # Find particle size rw from vfall
            env = VFallEnv(grav, mw_atmos, mfp, visc, t[iz], p[iz], rho_p)
            rw = 1e-6  # initial guess (cm)
            rw = vfall(rw, env)

            # Geometric mean radius / effective radius
            lnsig2 = 0.5 * np.log(sig[iz])**2
            alpha = np.log(vfall(rw*max(1.1,sig[iz]), env)/w_convect) / np.log(max(1.1,sig[iz]))
            rg[iz] = rainf[iz]**(1./alpha) * rw * np.exp(-(alpha+6)*lnsig2)
            reff[iz] = rg[iz] * np.exp(5*lnsig2)

            # Column number density
            ndz[iz] = 3 * rho_atmos * qc[iz] * dz[iz] / (4 * PI * rho_p * rg[iz]**3) * np.exp(-9*lnsig2)

    return qt, qc, ndz, rg, reff, lhf, kz
