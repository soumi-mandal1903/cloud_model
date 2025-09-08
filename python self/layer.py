# layer.py
import numpy as np
from advdiff import AdvDiff
from find_root import find_root

def layer(qt_guess, qvs, qbelow, mixl, dz, rainf, delta=1e-8):
    """
    Compute condensate in a single layer using AdvDiff class.
    """

    # Create an AdvDiff object for this layer
    adv = AdvDiff(ad_qbelow=qbelow, ad_qvs=qvs, ad_mixl=mixl, ad_dz=dz, ad_rainf=rainf)

    # Define root-finding wrapper
    f = lambda qt: adv.compute(qt)[0]  # return advdiff_val only

    # Bracket search range
    qmin = max(qvs, 1e-12)
    qmax = max(qt_guess * 2.0, qmin * 2.0)

    # Call root solver
    qt_new, status = find_root(f, 0.0, qmin, qmax, delta)

    # Compute condensate mixing ratio
    qc = max(0.0, qt_new - qvs)

    return qc, qt_new, status
