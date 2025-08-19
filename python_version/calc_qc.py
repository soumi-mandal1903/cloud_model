import numpy as np
from scipy.optimize import root_scalar

# Constants from globals.h
R_GAS = 8.3143e7  # erg/mol/K
PI = 3.14159274   # exactly three in Indiana

class CloudLayerCalculator:
    def __init__(self, R_GAS, PI):
        self.R_GAS = R_GAS
        self.PI = PI

    def pvap_gas(self, gas_name, t_layer, p_layer):
        # Stub: Replace with actual vapor pressure calculation
        return 1.0

    def vfall(self, r, grav, mw_atmos, mfp, visc, t_layer, p_layer, rho_p):
        # Stub: Replace with actual terminal velocity calculation
        return r * 1e3

    def advdiff(self, qt, q_below, qvs, mixl, dz_layer, rainf_layer):
        # Stub: Replace with actual advective-diffusive calculation
        return qt - (q_below + qvs) / 2

    def find_root(self, func, target, lo, hi, tol, args=()):
        try:
            sol = root_scalar(lambda x: func(x, *args) - target, bracket=[lo, hi], xtol=tol, method='secant')
            if sol.converged:
                return sol.root, 0
            else:
                return None, 1
        except ValueError:
            # Root not bracketed
            return None, -1
        except Exception:
            return None, 2

    def calc_qc(self, gas_name, rainf_layer, rho_p, mw_cloud, 
                q_below, supsat, w_convect, mixl,
                dz_layer, grav, mw_atmos, mfp, visc, t_layer, p_layer,
                sig_layer):
        # Outputs
        qc_layer = qt_layer = rg_layer = reff_layer = ndz_layer = qt_top = 0.0
        status_r = status_q = 0

        pvap = self.pvap_gas(gas_name, t_layer, p_layer)
        fs = supsat + 1
        rho_atmos = p_layer / (self.R_GAS / mw_atmos * t_layer)
        qvs = fs * pvap / ((self.R_GAS / mw_cloud) * t_layer) / rho_atmos

        if q_below < qvs:
            qt_layer = q_below
            qt_top = q_below
            qc_layer = rg_layer = reff_layer = ndz_layer = 0.0
        else:
            qhi = q_below
            qlo = qhi / 1e3
            delta_q = q_below / 1000.

            # Find qt_top using advdiff
            def advdiff_func(qt):
                return self.advdiff(qt, q_below, qvs, mixl, dz_layer, rainf_layer)
            qt_top, status_q = self.find_root(advdiff_func, 0.0, qlo, qhi, delta_q)

            if status_q != 0 or qt_top is None:
                return {
                    "qc_layer": None, "qt_layer": None, "rg_layer": None,
                    "reff_layer": None, "ndz_layer": None, "qt_top": None,
                    "status_r": status_r, "status_q": status_q
                }

            qt_layer = 0.5 * (q_below + qt_top)
            qc_layer = max(0.0, qt_layer - qvs)

            # Find rw_layer using vfall
            rlo = 1e-10
            rhi = 10.0
            delta_v = w_convect / 1000.

            def vfall_func(r):
                return self.vfall(r, grav, mw_atmos, mfp, visc, t_layer, p_layer, rho_p)
            rw_layer, status_r = self.find_root(vfall_func, w_convect, rlo, rhi, delta_v)

            if status_r != 0 or rw_layer is None:
                return {
                    "qc_layer": qc_layer, "qt_layer": qt_layer, "rg_layer": None,
                    "reff_layer": None, "ndz_layer": None, "qt_top": qt_top,
                    "status_r": status_r, "status_q": status_q
                }

            lnsig2 = 0.5 * np.log(sig_layer) ** 2
            sig_alpha = max(1.1, sig_layer)


            # Compute alpha 
            if rainf_layer > 1:
                alpha = np.log(
                    self.vfall(rw_layer * sig_alpha, grav, mw_atmos, mfp, visc, t_layer, p_layer, rho_p) / w_convect
                ) / np.log(sig_alpha)
            else:
                alpha = np.log(
                    w_convect / self.vfall(rw_layer / sig_alpha, grav, mw_atmos, mfp, visc, t_layer, p_layer, rho_p)
                ) / np.log(sig_alpha)

            rg_layer = rainf_layer ** (1.0 / alpha) * rw_layer * np.exp(-(alpha + 6) * lnsig2)
            reff_layer = rg_layer * np.exp(5 * lnsig2)
            ndz_layer = 3 * rho_atmos * qc_layer * dz_layer / (4 * self.PI * rho_p * rg_layer ** 3) * np.exp(-9 * lnsig2)

        return {
            "qc_layer": qc_layer, "qt_layer": qt_layer, "rg_layer": rg_layer,
            "reff_layer": reff_layer, "ndz_layer": ndz_layer, "qt_top": qt_top,
            "status_r": status_r, "status_q": status_q
        }

# Example usage:
# calculator = CloudLayerCalculator(R_GAS, PI)
# result = calculator.calc_qc(...your arguments...)