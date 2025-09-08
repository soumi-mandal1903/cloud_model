import numpy as np

import numpy as np

def find_root(func, y, xlow, xhigh, delta, max_iter=50, epsilon=1e-7):
    """
    Secant root finder for solving func(x) - y = 0 between [xlow, xhigh].

    Parameters
    ----------
    func : callable
        Function f(x).
    y : float
        Target value, solve for f(x) = y.
    xlow, xhigh : float
        Bracketing interval for root.
    delta : float
        Convergence tolerance for |f(x)-y|.
    max_iter : int
        Maximum iterations (default 50).
    epsilon : float
        Convergence tolerance for relative x changes.

    Returns
    -------
    xnew : float
        Estimated root.
    status : int
        0 = success
        1 = max iterations reached without convergence
       -1 = root not bracketed initially
    """
    x1, x2 = xlow, xhigh
    f1 = func(x1) - y
    f2 = func(x2) - y

    # Default return values
    xnew, fnew = x1, f1

    # Abort if not bracketed
    if f1 * f2 > 0:
        return xnew, -1

    iter_count = 0
    while abs(fnew) > delta and iter_count < max_iter and abs(x2 / x1 - 1.0) > epsilon:
        # Secant step
        slope = (f2 - f1) / (x2 - x1)
        xnew = x1 - f1 / slope
        fnew = func(xnew) - y

        # Update bracket
        if fnew * f1 <= 0 and abs(x2 / xnew - 1.0) > epsilon:
            x2, f2 = xnew, fnew
        else:
            x1, f1 = xnew, fnew

        iter_count += 1

    # Status flag
    if iter_count == max_iter and abs(x2 / x1 - 1.0) > epsilon:
        status = 1
    else:
        status = 0

    return xnew, status
