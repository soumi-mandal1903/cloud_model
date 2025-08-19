import numpy as np

class RootFinderLog:
    def __init__(self, max_iter=50, epsilon=1e-7):
        self.max_iter = max_iter
        self.epsilon = epsilon

    def find_rootl(self, func, y, xlow, xhigh, delta, *args):
        """
        Finds x such that log(func(x)) - log(y) = 0 using the secant method.
        Returns (xnew, status):
          status = 0: success
          status = 1: max iterations reached
          status = -1: root not bracketed
        """
        x1 = xlow
        x2 = xhigh

        try:
            f1 = np.log(func(x1, *args)) - np.log(y)
            f2 = np.log(func(x2, *args)) - np.log(y)
        except (ValueError, ZeroDivisionError, FloatingPointError):
            return None, -1

        xnew = x1
        fnew = f1

        if f1 * f2 > 0:
            return None, -1  # root not bracketed

        delta_log = np.log(delta)
        iter_count = 0

        while (abs(fnew) > delta_log and iter_count < self.max_iter and abs(x2 / x1 - 1.0) > self.epsilon):
            slope = (f2 - f1) / (x2 - x1)
            if slope == 0:
                return xnew, 1
            xnew = x1 - f1 / slope
            try:
                fnew = np.log(func(xnew, *args)) - np.log(y)
            except (ValueError, ZeroDivisionError, FloatingPointError):
                return xnew, -1

            if fnew * f1 <= 0 and abs(x2 / xnew - 1.0) > self.epsilon:
                x2 = xnew
                f2 = fnew
            else:
                x1 = xnew
                f1 = fnew

            iter_count += 1

        if iter_count == self.max_iter and abs(x2 / x1 - 1.0) > self.epsilon:
            status = 1
        else:
            status = 0

        return xnew, status

# Example usage:
# def my_func(x): return x**2 + 1
# finder = RootFinderLog()
# root, status = finder.find_rootl(my_func, y=2, xlow=0.5, xhigh=2.0, delta=1e-6)