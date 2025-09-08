# advdiff.py
import numpy as np

class AdvDiff:
    """
    Advective-diffusive balance calculator for condensate in a model layer,
    following the original Fortran algorithm from Ackerman & Marley (2000).
    """

    def __init__(self, ad_qbelow, ad_qvs, ad_mixl, ad_dz, ad_rainf):
        """
        Parameters
        ----------
        ad_qbelow : float
            Total vapor+condensate mixing ratio in underlying layer (g/g)
        ad_qvs : float
            Saturation mixing ratio (g/g)
        ad_mixl : float
            Convective mixing length (cm)
        ad_dz : float
            Layer thickness (cm)
        ad_rainf : float
            Rain efficiency factor
        """
        self.ad_qbelow = ad_qbelow
        self.ad_qvs = ad_qvs
        self.ad_mixl = ad_mixl
        self.ad_dz = ad_dz
        self.ad_rainf = ad_rainf
        self.ad_qc = 0.0

    def compute(self, qt):
        """
        Compute advective-diffusive residual for a given total mixing ratio.

        Parameters
        ----------
        qt : float
            Total mixing ratio (vapor + condensate) [g/g]

        Returns
        -------
        advdiff_val : float
            Residual from advective-diffusive balance (zero when balanced)
        ad_qc : float
            Condensed mixing ratio (g/g)
        """
        # all vapor in excess of saturation condenses (supsat = 0)
        self.ad_qc = max(0.0, qt - self.ad_qvs)

        # difference from advective-diffusive balance
        advdiff_val = (
            self.ad_qbelow
            * np.exp(-self.ad_rainf * self.ad_qc * self.ad_dz / (qt * self.ad_mixl))
            - qt
        )

        return advdiff_val, self.ad_qc
