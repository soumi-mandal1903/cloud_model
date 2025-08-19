class VaporPressure:
    @staticmethod
    def pvap_mns(t, metallicityMH=0.0):
        """
        Calculate saturation vapor pressure (dyne/cm^2) over MnS.
        Input: t (temperature in K), metallicityMH (default 0.0)
        Output: vapor pressure in dyne/cm^2
        """
        pvap_mns_bars = 10.0 ** (11.5315 - 23810.0 / t - metallicityMH)
        pvap_mns = pvap_mns_bars * 1e6  # convert from bars to dyne/cm^2
        return pvap_mns

# Example usage:
# vp = VaporPressure.pvap_mns(1800)
# print(vp)