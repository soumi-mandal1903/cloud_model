class VaporPressure:
    @staticmethod
    def pvap_na2s(t, metallicityMH=0.0):
        """
        Calculate saturation vapor pressure (dyne/cm^2) over Na2S.
        Input: t (temperature in K), metallicityMH (default 0.0)
        Output: vapor pressure in dyne/cm^2
        """
        pvap_na2s_bars = 10.0 ** (8.5497 - 13889.0 / t - 0.5 * metallicityMH)
        pvap_na2s = pvap_na2s_bars * 1e6  # convert from bars to dyne/cm^2
        return pvap_na2s

# Example usage:
# vp = VaporPressure.pvap_na2s(1200)
# print(vp)