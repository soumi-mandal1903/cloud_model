import math

import numpy as np

class AdvDiff:
    """
    Advective-diffusive balance calculator for condensate in a model layer.
    """

    def __init__(self, ad_qbelow, ad_qvs, ad_mixl, ad_dz, ad_rainf):
        self.ad_qbelow = ad_qbelow
        self.ad_qvs = ad_qvs
        self.ad_mixl = ad_mixl
        self.ad_dz = ad_dz
        self.ad_rainf = ad_rainf
        self.ad_qc = 0.0

    def compute(self, qt):
        """Return advdiff value and update ad_qc."""
        self.ad_qc = max(0.0, qt - self.ad_qvs)

        advdiff_val = (
            self.ad_qbelow
            * np.exp(-self.ad_rainf * self.ad_qc * self.ad_dz / (qt * self.ad_mixl))
            - qt
        )

        return advdiff_val, self.ad_qc
