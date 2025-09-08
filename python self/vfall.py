# vfall_main.py
import numpy as np

class VFallEnv:
    """
    Container for atmospheric + particle properties (cgs units).
    """
    def __init__(self, grav, mw_atmos, mfp, visc, T, P, rhop, mw_cloud=None):
        self.vf_grav = grav        # gravity (cm/s^2)
        self.vf_mw_atmos = mw_atmos  # molecular weight (g/mol)
        self.vf_mfp = mfp          # mean free path (cm)
        self.vf_visc = visc        # dynamic viscosity (dyne s/cm^2)
        self.vf_t = T              # temperature (K)
        self.vf_p = P              # pressure (dyne/cm^2)
        self.vf_rhop = rhop        # particle density (g/cm^3)
        self.vf_mw_cloud = mw_cloud

def vfall(r, env):
    """
    Compute terminal fall speed of a particle in cgs units.
    """
    # Constants
    R_GAS = 8.3143e7  # erg/mol/K
    b1 = 0.8
    b2 = -0.01
    cdrag = 0.45

    g = env.vf_grav
    mw = env.vf_mw_atmos
    mfp = env.vf_mfp
    visc = env.vf_visc
    T = env.vf_t
    P = env.vf_p
    rhop = env.vf_rhop

    # Knudsen number
    knudsen = mfp / r

    # Atmospheric density
    rho_atmos = P / ((R_GAS / mw) * T)
    drho = rhop - rho_atmos

    # Slip correction
    slip = 1.0 + 1.26 * knudsen

    # Stokes velocity
    v = slip * (2.0/9.0) * drho * g * r**2 / visc
    reynolds = 2.0 * r * rho_atmos * v / visc

    if reynolds > 1.0:
        x = np.log(reynolds)
        y = b1 * x + b2 * x**2
        reynolds = np.exp(y)
        v = visc * reynolds / (2.0 * r * rho_atmos)
        if reynolds > 1e3:
            v = slip * np.sqrt(8.0 * drho * r * g / (3.0 * cdrag * rho_atmos))

    return v

# Old Fortran main loop translated
def main():
    grav = 980.0
    mw_atmos = 29.0
    mfp = 6.2e-6
    visc = 1.6e-4
    t_layer = 288.0
    p_layer = 1e6
    rho_p = 1.0
    mw_cloud = 18.0  # example, needed for completeness

    rad = np.array([10., 100., 300., 1000.])  # microns

    env = VFallEnv(grav, mw_atmos, mfp, visc, t_layer, p_layer, rho_p, mw_cloud)

    print(" r (micron) | vfall_old | rossow | vfall_new")
    for r_um in rad:
        r_cm = r_um * 1e-4  # convert microns to cm
        v_old = vfall(r_cm, env)           # placeholder for old vfall
        v_rossow = vfall(r_cm, env)        # placeholder for rossow
        v_new = vfall(r_cm, env)           # current vfall implementation
        print(f"{r_um:10.1f} | {v_old:10.4e} | {v_rossow:10.4e} | {v_new:10.4e}")

if __name__ == "__main__":
    main()
