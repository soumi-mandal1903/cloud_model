class VaporPressure:
    @staticmethod
    def pvap_zns(t, metallicityMH=0.0):
        """
        Calculate saturation vapor pressure (dyne/cm^2) over ZnS.
        Input: t (temperature in K), metallicityMH (default 0.0)
        Output: vapor pressure in dyne/cm^2
        """
        pvap_zns_bars = 10.0 ** (12.8117 - 15873.0 / t - metallicityMH)
        pvap_zns = pvap_zns_bars * 1e6  # convert from bars to dyne/cm^2
        return pvap_zns

# Example usage:
# vp = VaporPressure.pvap_zns(1200)
# print(vp)