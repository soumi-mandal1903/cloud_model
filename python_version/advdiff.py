import numpy as np

class AdvDiffCalculator:
    def __init__(self):
        # These mimic the Fortran common block variables
        self.ad_qbelow = None
        self.ad_qvs = None
        self.ad_mixl = None
        self.ad_dz = None
        self.ad_rainf = None

    def set_params(self, qbelow, qvs, mixl, dz, rainf):
        self.ad_qbelow = qbelow
        self.ad_qvs = qvs
        self.ad_mixl = mixl
        self.ad_dz = dz
        self.ad_rainf = rainf

    def advdiff(self, qt):
        # All vapor in excess of saturation condenses (supsat=0)
        ad_qc = max(0.0, qt - self.ad_qvs)
        # Difference from advective-diffusive balance
        if qt == 0 or self.ad_mixl == 0:
            # Avoid division by zero
            return -qt
        exponent = -self.ad_rainf * ad_qc * self.ad_dz / (qt * self.ad_mixl)
        return self.ad_qbelow * np.exp(exponent) - qt

# Example usage:
# advcalc = AdvDiffCalculator()
# advcalc.set_params(qbelow, qvs, mixl, dz, rainf)
# result = advcalc.advdiff(qt)