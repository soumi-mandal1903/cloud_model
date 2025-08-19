import numpy as np

class VfallCalculator:
    def __init__(self, vf_grav, vf_mw_atmos, vf_mfp, vf_visc, vf_t, vf_p, vf_rhop):
        """
        Initialize with atmospheric and particle properties.
        All units are cgs.
        """
        self.vf_grav = vf_grav      # gravity (cm/s^2)
        self.vf_mw_atmos = vf_mw_atmos  # atmospheric molecular weight (g/mol)
        self.vf_mfp = vf_mfp        # mean free path (cm)
        self.vf_visc = vf_visc      # dynamic viscosity (dyne s/cm^2)
        self.vf_t = vf_t            # temperature (K)
        self.vf_p = vf_p            # pressure (dyne/cm^2)
        self.vf_rhop = vf_rhop      # particle density (g/cm^3)

    def vfall(self, r):
        """
        Calculate fallspeed for a spherical particle of radius r (cm).
        Returns fallspeed in cm/s.
        """
        # Parameters
        b1 = 0.8
        b2 = -0.01
        cdrag = 0.45
        R_GAS = 8.3143e7  # erg/mol/K

        knudsen = self.vf_mfp / r
        rho_atmos = self.vf_p / ((R_GAS / self.vf_mw_atmos) * self.vf_t)
        drho = self.vf_rhop - rho_atmos

        # Cunningham correction (slip factor)
        slip = 1. + 1.26 * knudsen

        # Stokes terminal velocity (low Reynolds number)
        vfall = slip * (2. / 9.) * drho * self.vf_grav * r ** 2 / self.vf_visc
        reynolds = 2. * r * rho_atmos * vfall / self.vf_visc

        if reynolds > 1.0:
            # Correct drag coefficient for turbulence
            x = np.log(reynolds)
            y = b1 * x + b2 * x ** 2
            reynolds = np.exp(y)
            vfall = self.vf_visc * reynolds / (2. * r * rho_atmos)

            if reynolds > 1e3:
                # Drag coefficient independent of Reynolds number
                vfall = slip * np.sqrt(8. * drho * r * self.vf_grav / (3. * cdrag * rho_atmos))

        return vfall

# Example usage:
# vfall_calc = VfallCalculator(vf_grav=980, vf_mw_atmos=29, vf_mfp=1e-5, vf_visc=1.8e-4,
#                              vf_t=300, vf_p=1e6, vf_rhop=2.5)
# speed = vfall_calc.vfall(1e-4)
# print(speed)