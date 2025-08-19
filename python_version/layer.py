import numpy as np

class LayerCalculator:
    def __init__(self, R_GAS=8.3143e7, PI=3.14159274, AVOGADRO=6.02e23):
        self.R_GAS = R_GAS
        self.PI = PI
        self.AVOGADRO = AVOGADRO
        self.K_BOLTZ = R_GAS / AVOGADRO

    def calc_qc(self, gas_name, rainf, rho_p, mw_cloud,
                qt_below, supsat, w_convect, mixl,
                dz_sub, grav, mw_atmos, mfp, visc, t_sub, p_sub,
                sig_layer):
        # Stub: Replace with actual calc_qc logic
        # Returns qc_sub, qt_sub, rg_sub, reff_sub, ndz_sub, qt_top, status_r, status_q
        qc_sub = qt_sub = rg_sub = reff_sub = ndz_sub = qt_top = 0.0
        status_r = status_q = 0
        return qc_sub, qt_sub, rg_sub, reff_sub, ndz_sub, qt_top, status_r, status_q

    def layer(self, nsub_max, gas_name, grav, mw_atmos, kz_min, cloudf_min,
              mw_cloud, rainf, rho_p, supsat, sig_layer, cloudf, q_below,
              t_layer, p_layer, kz, chf, t_top, t_bot, p_top, p_bot,
              report_status_r=True, report_status_q=True):
        # Physical constants
        d_molecule = 2.827e-8
        eps_k = 59.7
        r_atmos = self.R_GAS / mw_atmos
        c_p = 7. / 2. * r_atmos

        dp_layer = p_bot - p_top
        dlnp = np.log(p_bot / p_top)
        dtdlnp = (t_top - t_bot) / dlnp
        lapse_ratio = (t_bot - t_top) / dlnp / (2. / 7. * t_layer)
        rho_atmos = p_layer / (r_atmos * t_layer)
        scale_h = r_atmos * t_layer / grav
        mixl = max(0.1, lapse_ratio) * scale_h
        scalef_kz = 1. / 3.
        kz_val = scalef_kz * scale_h * (mixl / scale_h) ** (4. / 3.) * ((r_atmos * chf) / (rho_atmos * c_p)) ** (1. / 3.)
        kz_val = max(kz_val, kz_min)
        w_convect = kz_val / mixl
        cloudf_val = cloudf_min + max(0., min(1., 1. - lapse_ratio)) * (1. - cloudf_min)
        n_atmos = p_layer / (self.K_BOLTZ * t_layer)
        mfp = 1. / (np.sqrt(2.) * n_atmos * self.PI * d_molecule ** 2)
        visc = (5. / 16.) * np.sqrt(self.PI * self.K_BOLTZ * t_layer * (mw_atmos / self.AVOGADRO)) / \
               (self.PI * d_molecule ** 2) / (1.22 * (t_layer / eps_k) ** (-0.16))

        status_r = 0
        status_q = 0
        nsub = 1
        converge = False

        while not converge:
            qc_layer = 0.
            qt_layer = 0.
            ndz_layer = 0.
            opd_layer = 0.
            qt_bot_sub = q_below
            p_bot_sub = p_bot
            dp_sub = dp_layer / nsub

            for isub in range(nsub):
                qt_below = qt_bot_sub
                p_top_sub = p_bot_sub - dp_sub
                dz_sub = scale_h * np.log(p_bot_sub / p_top_sub)
                p_sub = 0.5 * (p_bot_sub + p_top_sub)
                t_sub = t_bot + np.log(p_bot / p_sub) * dtdlnp

                qc_sub, qt_sub, rg_sub, reff_sub, ndz_sub, qt_top, status_r, status_q = self.calc_qc(
                    gas_name, rainf, rho_p, mw_cloud,
                    qt_below, supsat, w_convect, mixl,
                    dz_sub, grav, mw_atmos, mfp, visc, t_sub, p_sub,
                    sig_layer
                )

                qc_layer += qc_sub * dp_sub / grav
                qt_layer += qt_sub * dp_sub / grav
                ndz_layer += ndz_sub

                if reff_sub > 0.:
                    opd_layer += 1.5 * qc_sub * dp_sub / grav / (rho_p * reff_sub)

                qt_bot_sub = qt_top
                p_bot_sub = p_top_sub

            # Convergence check
            if nsub_max == 1:
                converge = True
            elif nsub == 1:
                opd_test = opd_layer
            elif opd_layer == 0. or nsub >= nsub_max:
                converge = True
            elif abs(1. - opd_test / opd_layer) <= 1e-2:
                converge = True
            else:
                opd_test = opd_layer

            nsub *= 2

        # Report problems finding root the first time it happens
        if status_r != 0 and report_status_r:
            print(f"layer(): find_root(vfall) status = {status_r} for {gas_name} at p = {p_layer / 1e6:.2f}")
            print(" there may be more instances not reported")
            print(f"status_r = {status_r}")
            report_status_r = False

        if status_q != 0 and report_status_q:
            print(f"layer(): find_root(advdiff) status = {status_q} for {gas_name} at p = {p_layer / 1e6:.2f}")
            print(" there may be more instances not reported")
            print(f"status_q = {status_q}")
            report_status_q = False

        q_below = qt_top

        if opd_layer > 0.:
            reff_layer = 1.5 * qc_layer / (rho_p * opd_layer)
            lnsig2 = 0.5 * np.log(sig_layer) ** 2
            rg_layer = reff_layer * np.exp(-5 * lnsig2)
        else:
            reff_layer = 0.
            rg_layer = 0.

        qc_layer = qc_layer * grav / dp_layer
        qt_layer = qt_layer * grav / dp_layer

        return {
            "qc_layer": qc_layer,
            "qt_layer": qt_layer,
            "rg_layer": rg_layer,
            "reff_layer": reff_layer,
            "ndz_layer": ndz_layer,
            "opd_layer": opd_layer,
            "status_r": status_r,
            "status_q": status_q,
            "q_below": q_below
        }

# Example usage:
# layer_calc = LayerCalculator()
# result = layer_calc.layer(nsub_max, gas_name, grav, mw_atmos, kz_min, cloudf_min,
#                           mw_cloud, rainf, rho_p, supsat, sig_layer, cloudf, q_below,
#                           t_layer, p_layer, kz, chf, t_top, t_bot, p_top, p_bot)