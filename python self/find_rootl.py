# find_rootl.py
import numpy as np

def find_rootl(f, y, xlow, xhigh, delta, max_iter=50, eps=1e-7):
    """
    Solve log(f(x)) - log(y) = 0 for x in [xlow, xhigh] using secant method.

    Parameters
    ----------
    f : callable
        Function f(x).
    y : float
        Target value.
    xlow, xhigh : float
        Bracketing interval.
    delta : float
        Convergence tolerance on log(f(x)).
    max_iter : int
        Maximum number of iterations (default: 50).
    eps : float
        Convergence tolerance on x interval.

    Returns
    -------
    xnew : float
        Estimated root.
    status : int
        0 = success
        1 = max iterations reached
       -1 = no convergence (root not bracketed)
    """

    x1, x2 = xlow, xhigh

    try:
        f1 = np.log(f(x1)) - np.log(y)
        f2 = np.log(f(x2)) - np.log(y)
    except (ValueError, ZeroDivisionError, FloatingPointError):
        return x1, -1

    # default
    xnew, fnew = x1, f1

    # check bracketing
    if f1 * f2 > 0:
        return xnew, -1

    # take logarithm of precision
    delta = np.log(delta)

    iter_count = 0
    while (abs(fnew) > delta and
           iter_count < max_iter and
           abs(x2 / x1 - 1.0) > eps):

        slope = (f2 - f1) / (x2 - x1)
        xnew = x1 - f1 / slope

        try:
            fnew = np.log(f(xnew)) - np.log(y)
        except (ValueError, ZeroDivisionError, FloatingPointError):
            return xnew, -1

        # keep bracketing
        if fnew * f1 <= 0 and abs(x2 / xnew - 1.0) > eps:
            x2, f2 = xnew, fnew
        else:
            x1, f1 = xnew, fnew

        iter_count += 1

    # status flag
    if iter_count == max_iter and abs(x2 / x1 - 1.0) > eps:
        return xnew, 1
    else:
        return xnew, 0
