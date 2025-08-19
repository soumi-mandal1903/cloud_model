import numpy as np

class VaporPressure:
    @staticmethod
    def pvap_cr(t):
        """
        Calculate saturation vapor pressure (dyne/cm^2) over Cr.
        Input: t (temperature in K)
        Output: vapor pressure in dyne/cm^2
        """
        pvap_cr_bars = 10.0 ** (7.2688 - 20353.0 / t)
        pvap_cr = pvap_cr_bars * 1e6  # convert from bars to dyne/cm^2
        return pvap_cr

# Example usage:
# vp = VaporPressure.pvap_cr(1800)
# print(vp)