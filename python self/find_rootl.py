def find_rootl(f, y, xlow, xhigh, delta):
    """
    Secant root finder solving log(f(x)) - log(y) = 0
    between [xlow, xhigh], mimicking the Fortran algorithm.

    Returns:
        xnew: estimated root
        status: 0 = success
                1 = max iterations reached without convergence
               -1 = root not bracketed initially
    """
    MAX_ITER = 50
    EPSILON = 1e-7

    # copy input range
    x1 = xlow
    x2 = xhigh

    # evaluate function at endpoints in log space
    f1 = np.log(f(x1)) - np.log(y)
    f2 = np.log(f(x2)) - np.log(y)

    # default return values
    xnew = x1
    fnew = f1

    # abort if root not bracketed
    if f1 * f2 > 0:
        status = -1
        return xnew, status

    # take logarithm of precision
    delta_log = np.log(delta)

    iter_count = 0
    while abs(fnew) > delta_log and iter_count < MAX_ITER and (x2 / x1 - 1.0) > EPSILON:
        # secant estimate
        slope = (f2 - f1) / (x2 - x1)
        xnew = x1 - f1 / slope
        fnew = np.log(f(xnew)) - np.log(y)

        # update bracket
        if fnew * f1 <= 0 and (x2 / xnew - 1.0) > EPSILON:
            x2 = xnew
            f2 = fnew
        else:
            x1 = xnew
            f1 = fnew

        iter_count += 1

    # determine status
    if iter_count == MAX_ITER and (x2 / x1 - 1.0) > EPSILON:
        status = 1
    else:
        status = 0

    return xnew, status
