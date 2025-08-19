import numpy as np

class VaporPressure:
    @staticmethod
    def pvap_fe(t):
        """
        Calculate saturation vapor pressure (dyne/cm^2) over Fe.
        Uses the latest expression from Channon Visscher (2011).
        Input: t (temperature in K)
        Output: vapor pressure in dyne/cm^2
        """
        pvap_fe = 10.0 ** (7.09 - 20833.0 / t)
        pvap_fe = pvap_fe * 1e6  # convert from bars to dyne/cm^2
        return pvap_fe

# Example usage:
# vp = VaporPressure.pvap_fe(1800)
# print(vp)