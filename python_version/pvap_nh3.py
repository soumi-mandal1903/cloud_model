import numpy as np

class VaporPressure:
    @staticmethod
    def pvap_nh3(t):
        """
        Calculate saturation vapor pressure (dyne/cm^2) over NH3.
        Input: t (temperature in K)
        Output: vapor pressure in dyne/cm^2
        """
        pvap_nh3 = np.exp(-86596.0 / t**2 - 2161.0 / t + 10.53)
        pvap_nh3 = pvap_nh3 * 1e6  # convert from bars to dyne/cm^2
        return pvap_nh3

# Example usage:
# vp = VaporPressure.pvap_nh3(200)
# print(vp)