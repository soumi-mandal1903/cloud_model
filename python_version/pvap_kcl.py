class VaporPressure:
    @staticmethod
    def pvap_kcl(t):
        """
        Calculate saturation vapor pressure (dyne/cm^2) over KCl.
        Input: t (temperature in K)
        Output: vapor pressure in dyne/cm^2
        """
        pvap_kcl_bars = 10.0 ** (7.6106 - 11382.0 / t)
        pvap_kcl = pvap_kcl_bars * 1e6  # convert from bars to dyne/cm^2
        return pvap_kcl

# Example usage:
# vp = VaporPressure.pvap_kcl(1200)
# print(vp)