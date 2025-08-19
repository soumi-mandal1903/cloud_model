import numpy as np

class VaporPressure:
    @staticmethod
    def pvap_h2o(t):
        """
        Calculate saturation vapor pressure (dyne/cm^2) over water.
        Uses Buck's expressions for liquid and ice.
        Input: t (temperature in K)
        Output: vapor pressure in dyne/cm^2
        """
        # Buck constants
        BAL = 6.1121e3
        BBL = 18.729
        BCL = 257.87
        BDL = 227.3

        BAI = 6.1115e3
        BBI = 23.036
        BCI = 279.82
        BDI = 333.7

        # Wexler constants (not used, but included for completeness)
        GG0 = -0.29912729e+4
        GG1 = -0.60170128e+4
        GG2 = 0.1887643854e+2
        GG3 = -0.28354721e-1
        GG4 = 0.17838301e-4
        GG5 = -0.84150417e-9
        GG6 = 0.44412543e-12
        GG7 = 0.28584870e+1

        HH0 = -0.58653696e+4
        HH1 = 0.2224103300e+2
        HH2 = 0.13749042e-1
        HH3 = -0.34031775e-4
        HH4 = 0.26967687e-7
        HH5 = 0.6918651

        DO_BUCK = True

        if t < 273.16:
            # Saturation vapor pressure over ice
            if DO_BUCK:
                tc = t - 273.16
                pvap_h2o = BAI * np.exp((BBI - tc / BDI) * tc / (tc + BCI))
            else:
                pvap_h2o = 10 * np.exp(1.0 / t *
                    (HH0 + (HH1 + HH5 * np.log(t) +
                    (HH2 + (HH3 + HH4 * t) * t) * t) * t))
        else:
            # Saturation vapor pressure over liquid water
            if DO_BUCK:
                if t < 1048.:
                    tc = t - 273.16
                    pvap_h2o = BAL * np.exp((BBL - tc / BDL) * tc / (tc + BCL))
                else:
                    pvap_h2o = 600e6
            else:
                pvap_h2o = 10 * np.exp((1.0 / (t * t)) *
                    (GG0 + (GG1 + (GG2 + GG7 * np.log(t) +
                    (GG3 + (GG4 + (GG5 + GG6 * t) * t) * t) * t) * t) * t))

        return pvap_h2o

# Example usage:
# vp = VaporPressure.pvap_h2o(300)
# print(vp)