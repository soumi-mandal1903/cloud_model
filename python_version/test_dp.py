import numpy as np

class TestDewpoint:
    def __init__(self):
        self.mw_atmos = 29.66

    def run(self):
        # Water
        mw_cloud = 18.0
        gas_mmr = 10e-3
        press = 1e6

        # Uncomment for ammonia
        # mw_cloud = 17.0
        # gas_mmr = 3e-5 * mw_cloud / self.mw_atmos
        # press = 42250.0

        eps = mw_cloud / self.mw_atmos

        # Range of temperatures to search [K]
        tlo = 100.0
        thi = 300.0

        # rhs in pvap(t_dew) = rhs [dyne/cm^2]
        rhs = press * gas_mmr * self.mw_atmos / mw_cloud

        # Precision of search [dyne/cm^2]
        delta_p = rhs / 100.0

        # Find dewpoint using root finder and vapor pressure
        root_finder = RootFinder()
        vapor_pressure = VaporPressure()

        def pvap_h2o_func(t):
            return vapor_pressure.pvap_h2o(t)

        dp, status_t = root_finder.find_root(pvap_h2o_func, rhs, tlo, thi, delta_p)

        print(f'delta_p, dewp = {delta_p:.3e}, {dp:.2f}')

        if status_t != 0:
            print(f'find_root(pvap) status = {status_t}')

# Example usage:
# test_dp = TestDewpoint()
# test_dp.run()