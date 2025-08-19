import numpy as np

class RootFinder:
    def __init__(self, max_iter=50, epsilon=1e-7):
        self.max_iter = max_iter
        self.epsilon = epsilon

    def find_root(self, func, y, xlow, xhigh, delta, *args):
        """
        Finds x such that func(x) - y = 0 using the secant method.
        Returns (xnew, status):
          status = 0: success
          status = 1: max iterations reached
          status = -1: root not bracketed
        """
        x1 = xlow
        x2 = xhigh

        try:
            f1 = func(x1, *args) - y
            f2 = func(x2, *args) - y
        except Exception:
            return None, -1

        xnew = x1
        fnew = f1

        if f1 * f2 > 0:
            return None, -1  # root not bracketed

        iter_count = 0

        while (abs(fnew) > delta and iter_count < self.max_iter and abs(x2 / x1 - 1.0) > self.epsilon):
            slope = (f2 - f1) / (x2 - x1)
            if slope == 0:
                return xnew, 1
            xnew = x1 - f1 / slope
            try:
                fnew = func(xnew, *args) - y
            except Exception:
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
# def my_func(x): return x**2 - 2
# finder = RootFinder()
# root, status = finder.find_root(my_func, y=0, xlow=1, xhigh=2, delta=1e-6)