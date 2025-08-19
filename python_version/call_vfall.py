import numpy as np

class VFallCalculator:
    def __init__(self, grav, mw_atmos, mfp, visc, t_layer, p_layer, rho_p, mw_cloud):
        self.grav = grav
        self.mw_atmos = mw_atmos
        self.mfp = mfp
        self.visc = visc
        self.t_layer = t_layer
        self.p_layer = p_layer
        self.rho_p = rho_p
        self.mw_cloud = mw_cloud

    def vfall_save(self, r):
        # Stub: Replace with actual implementation
        return r * 1e3

    def vfall_rossow(self, r):
        # Stub: Replace with actual implementation
        return r * 2e3

    def vfall(self, r):
        # Stub: Replace with actual implementation
        return r * 3e3

def main():
    rad = np.array([10., 100., 300., 1000.])  # micron
    grav = 980.
    mw_atmos = 29.
    mfp = 6.2e-6
    visc = 1.6e-4
    t_layer = 288.
    p_layer = 1e6
    rho_p = 1.
    mw_cloud = 18.  # Example value

    vfall_calc = VFallCalculator(grav, mw_atmos, mfp, visc, t_layer, p_layer, rho_p, mw_cloud)

    print("r (micron), vfall_old, rossow, new =")
    for r_micron in rad:
        r = r_micron * 1e-4  # convert micron to cm
        print(f"{r_micron:8.1f} {vfall_calc.vfall_save(r):12.4e} {vfall_calc.vfall_rossow(r):12.4e} {vfall_calc.vfall(r):12.4e}")

if __name__ == "__main__":
    main()