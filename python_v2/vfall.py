# vfall.py
import numpy as np

class VFallEnv:
    """
    Container for atmospheric + particle properties.
    Units: cgs (cm, g, dyne, s, K)
    """
    def __init__(self, grav, mw_atmos, mfp, visc, T, P, rhop):
        self.vf_grav = grav        # gravity (cm/s^2)
        self.vf_mw_atmos = mw_atmos  # molecular weight (g/mol)
        self.vf_mfp = mfp          # mean free path (cm)
        self.vf_visc = visc        # dynamic viscosity (dyne s/cm^2)
        self.vf_t = T              # temperature (K)
        self.vf_p = P              # pressure (dyne/cm^2)
        self.vf_rhop = rhop        # particle density (g/cm^3)

def vfall(r, env):
    """
    Calculate terminal fall speed of a particle.
    
    Parameters
    ----------
    r : float
        Particle radius (cm)
    env : VFallEnv
        Atmospheric/particle environment in cgs units.
    
    Returns
    -------
    vfall : float
        Terminal velocity (cm/s)
    """
    # constants
    R_GAS = 8.3143e7  # erg/mol/K = dyne*cm/mol/K
    b1 = 0.8          # Ackerman fit
    b2 = -0.01
    cdrag = 0.45
    
    # unpack env
    g = env.vf_grav
    mw = env.vf_mw_atmos
    mfp = env.vf_mfp
    visc = env.vf_visc
    T = env.vf_t
    P = env.vf_p
    rhop = env.vf_rhop
    
    # Knudsen number
    knudsen = mfp / r
    
    # atmospheric density (ideal gas law, cgs units)
    rho_atmos = P / ((R_GAS / mw) * T)
    drho = rhop - rho_atmos
    
    # slip correction
    slip = 1.0 + 1.26 * knudsen
    
    # Stokes terminal velocity
    v = slip * (2.0/9.0) * drho * g * r**2 / visc
    reynolds = 2.0 * r * rho_atmos * v / visc
    
    if reynolds > 1.0:
        # turbulent correction
        x = np.log(reynolds)
        y = b1 * x + b2 * x**2
        reynolds = np.exp(y)
        v = visc * reynolds / (2.0 * r * rho_atmos)
        
        if reynolds > 1e3:
            # high Re regime
            v = slip * np.sqrt(8.0 * drho * r * g / (3.0 * cdrag * rho_atmos))
    
    return v
