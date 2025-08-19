import numpy as np

class VaporPressure:
    @staticmethod
    def pvap_ch4(t):
        """
        Calculate saturation vapor pressure (dyne/cm^2) over CH4.
        Input: t (temperature in K)
        Output: vapor pressure in dyne/cm^2
        """
        AMR = 16.043 / 8.3143
        TCRIT = 90.68
        PCRIT = 0.11719
        AS = 2.213 - 2.650
        AL = 2.213 - 3.370
        ALS = 611.10
        ALV = 552.36

        # ic=1: temperature below triple point, ic=2: above
        ic = 1 if t <= TCRIT else 2

        C = np.zeros(2)
        B = np.zeros(2)
        A = np.zeros(2)

        C[0] = -AMR * AS
        C[1] = -AMR * AL
        B[0] = -AMR * (ALS + AS * TCRIT)
        B[1] = -AMR * (ALV + AL * TCRIT)
        A[0] = PCRIT * TCRIT ** (-C[0]) * np.exp(-B[0] / TCRIT)
        A[1] = PCRIT * TCRIT ** (-C[1]) * np.exp(-B[1] / TCRIT)

        idx = ic - 1  # Python uses 0-based indexing
        pvap_ch4 = A[idx] * t ** C[idx] * np.exp(B[idx] / t)
        pvap_ch4 = pvap_ch4 * 1e6  # convert from bars to dyne/cm^2

        return pvap_ch4

# Example usage:
# vp = VaporPressure.pvap_ch4(100)
# print(vp)