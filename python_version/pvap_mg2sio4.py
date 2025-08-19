import numpy as np

class VaporPressure:
    @staticmethod
    def pvap_mg2sio4(t, p, metallicityMH=0.0):
        """
        Calculate saturation vapor pressure (dyne/cm^2) over Mg2SiO4.
        Uses the latest expression with pressure and metallicity dependence.
        Input: t (temperature in K), p (pressure in dyne/cm^2), metallicityMH (default 0.0)
        Output: vapor pressure in dyne/cm^2
        """
        # p is expected in dyne/cm^2, convert to bars for log10
        pvap_mg2sio4 = 10.0 ** (
            -32488.0 / t + 14.88 - 0.2 * np.log10(p * 1e6)
            - 1.4 * metallicityMH
        ) * 1e6  # convert from bars to dyne/cm^2
        return pvap_mg2sio4

# Example usage:
# vp = VaporPressure.pvap_mg2sio4(1800, 1e6)
# print(vp)