import numpy as np

class VaporPressure:
    @staticmethod
    def pvap_mgsio3(t, metallicityMH=0.0):
        """
        Calculate saturation vapor pressure (dyne/cm^2) over MgSiO3.
        Uses the latest expression with metallicity dependence.
        Input: t (temperature in K), metallicityMH (default 0.0)
        Output: vapor pressure in dyne/cm^2
        """
        pvap_mgsio3 = 10.0 ** (11.83 - 27250.0 / t - metallicityMH) * 1e6
        return pvap_mgsio3

# Example usage:
# vp = VaporPressure.pvap_mgsio3(1800)
# print(vp)