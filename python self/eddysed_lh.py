# eddysed_nh3.py
import numpy as np
from vfall import vfall, VFallEnv

R_GAS = 8.3143e7  # erg/mol/K
K_BOLTZ = 1.3807e-16  # erg/K
AVOGADRO = 6.022e23
PI = np.pi
ZERO = 0.0

def pvap_nh3(T):
    """
    Vapor pressure of NH3 (dyne/cm^2)
    Placeholder function; replace with accurate Clausius-Clapeyron relation.
    """
    # Example: T in K, returns dyne/cm^2
    return 1e5 * np.exp(-4000/T)  

def eddysed_nh3(
    grav_in, teff, sig, rainf, nz, z, p, t, dlnp, chf, lapse_ratio,
    gas_mmr, mw_atmos_in
):
    """
    Eddy-sedimentation and cloud microphysics solver for NH3 only.
    
    Returns: kz, lhf, qt, qc, ndz, rg, reff
    """
    # Allocate arrays
    qc = np.zeros(nz)
    qt = np.zeros(nz)
    ndz = np.zeros(nz)
    rg = np.zeros(nz)
    reff = np.zeros(nz)
    kz = np.zeros(nz)
    lhf = np.zeros(nz)
    mixl = np.zeros(nz)
    dz = np.zeros(nz)
    
    # Common block for vfall
    grav = grav_in
    mw_atmos = mw_atmos_in
    
    d_molecule = 2.827e-8
    eps_k = 59.7
    r_atmos = R_GAS / mw_atmos
    c_p = 7./2. * r_atmos
    
    cloud_base = 0
    fs = 1.0  # saturation factor
    
    for iz in range(nz):
        # Atmospheric density
        rho_atmos = p[iz] / (r_atmos * t[iz])
        scale_h = r_atmos * t[iz] / grav
        dz[iz] = scale_h * dlnp[iz]
        mixl[iz] = max(0.1, lapse_ratio[iz]) * scale_h
        n_atmos = p[iz] / (K_BOLTZ * t[iz])
        mfp = 1. / (np.sqrt(2.)*n_atmos*PI*d_molecule**2)
        visc = (5./16.*np.sqrt(PI*K_BOLTZ*t[iz]*(mw_atmos/AVOGADRO)) /
                (PI*d_molecule**2) / (1.22 * (t[iz]/eps_k)**(-0.16)))
        
        t_layer = t[iz]
        p_layer = p[iz]
        
        # NH3 properties
        pvap = pvap_nh3(t[iz])
        mw_cloud = 17.
        lh = 28.6e10
        rho_p = 0.84
        
        # Saturation mass mixing ratio
        qvs = fs * pvap / ( (R_GAS/mw_cloud) * t[iz] ) / rho_atmos
        
        # Mass mixing ratio below layer
        if iz == 0:
            qbelow = gas_mmr
            zdiff = dz[iz]/2.
        else:
            qbelow = qt[iz-1]
            zdiff = z[iz] - z[iz-1]
        
        if qbelow < qvs:
            qt[iz] = qbelow
            qc[iz] = 0.
            rg[iz] = 0.
            reff[iz] = 0.
            ndz[iz] = 0.
        else:
            # Solve for total mixing ratio (adv-diff)
            qlo = qbelow / 1e3
            qhi = qbelow
            delta_q = qbelow / 1000.
            # Placeholder for root-finding; in practice, replace with advdiff
            qt[iz] = qhi
            qc[iz] = max(0., qt[iz] - qvs)
            
            lh_sum = rainf[iz] * qc[iz] * lh / (mw_cloud * c_p / r_atmos)
            
            # Convective velocity scale (cubic)
            coeff_a = lh_sum
            coeff_b = -chf[iz] / (rho_atmos * c_p / r_atmos)
            if lh_sum == 0.:
                w_convect = (-coeff_b)**(1./3.)
            else:
                cubic_arg = coeff_b**2 / 4. + coeff_a**3 / 27.
                w_convect = ((-coeff_b/2. + np.sqrt(cubic_arg))**(1./3.) +
                             (-coeff_b/2. - np.sqrt(cubic_arg))**(1./3.))
            
            # Eddy diffusion coefficient
            kz[iz] = w_convect * mixl[iz]
            
            # Find particle radius corresponding to convective speed
            rlo = 1e-10
            rhi = 10.
            delta_v = w_convect / 1000.
            env = VFallEnv(grav, mw_atmos, mfp, visc, t[iz], p[iz], rho_p)
            # Placeholder: vfall root-finding; here we just pick mid-radius
            rw = (rlo + rhi)/2. 
            
            # Lognormal microphysics
            sig_alpha = max(1.1, sig[iz])
            alpha = np.log(w_convect / vfall(rw/sig_alpha, env)) / np.log(sig_alpha)
            lnsig2 = 0.5*np.log(sig[iz])**2
            rg[iz] = rw * np.exp(-(alpha+6)*lnsig2)
            reff[iz] = rg[iz]*np.exp(5*lnsig2)
            ndz[iz] = 3*rho_atmos*qc[iz]*dz[iz]/(4*PI*rho_p*rg[iz]**3) * np.exp(-9*lnsig2)
            
            if cloud_base == 0:
                cloud_base = iz
            
            # Latent heat flux
            lhf[iz] = w_convect * lh_sum * rho_atmos * c_p / r_atmos
    
    # Condensate path
    cond_path = sum(qc[iz]*dlnp[iz]*p[iz]/grav for iz in range(nz))
    if cond_path > 0.:
        print(f"eddysed(): NH3 condensate path (g/m^2) = {cond_path*1e4:.3f}")
    
    return kz, lhf, qt, qc, ndz, rg, reff
