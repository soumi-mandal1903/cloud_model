import numpy as np

class VaporPressure:
    @staticmethod
    def pvap_al2o3(t):
        """
        Calculate saturation vapor pressure (dyne/cm^2) over Al2O3.
        Kozasa et al. Ap J. 344 325
        Input: t (temperature in K)
        Output: vapor pressure in dyne/cm^2
        """
        return np.exp(-73503.0 / t + 22.01) * 1e6

# Example usage:
# vp = VaporPressure.pvap_al2o3(1800)
# print(vp)