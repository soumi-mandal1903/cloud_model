# find_root.py
import numpy as np

def find_root(f, target, xlo, xhi, tol):
    """
    Secant root finder, translated from Fortran find_root.f
    
    Solves f(x) - target = 0 for x between xlo and xhi to within tol.
    
    Returns
    -------
    xnew : float
        Root estimate
    status : int
        0 = converged
        1 = max iterations reached
       -1 = root not bracketed
    """
    
    MAX_ITER = 50
    EPSILON = 1e-7
    
    # initial endpoints
    x1, x2 = xlo, xhi
    f1 = f(x1) - target
    f2 = f(x2) - target
    
    # default outputs
    xnew = x1
    fnew = f1
    
    # check if root is bracketed
    if f1 * f2 > 0:
        return xnew, -1   # not bracketed
    
    iter_count = 0
    while (abs(fnew) > tol) and (iter_count < MAX_ITER) and (abs(x2/x1 - 1.0) > EPSILON):
        # secant step
        slope = (f2 - f1) / (x2 - x1)
        xnew = x1 - f1 / slope
        fnew = f(xnew) - target
        
        # update bracket
        if (fnew * f1 <= 0) and (abs(x2/xnew - 1.0) > EPSILON):
            x2, f2 = xnew, fnew
        else:
            x1, f1 = xnew, fnew
        
        iter_count += 1
    
    # set status flag
    if iter_count == MAX_ITER and (abs(x2/x1 - 1.0) > EPSILON):
        status = 1   # max iterations reached
    else:
        status = 0   # converged
    
    return xnew, status
