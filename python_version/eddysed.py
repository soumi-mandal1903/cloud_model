import numpy as np

# Constants from globals.h
R_GAS = 8.3143e7  # erg/mol/K
PI = 3.14159274

class EddySedimentationModel:
    def __init__(self, maxnz=50, maxngas=10):
        self.maxnz = maxnz
        self.maxngas = maxngas

    def pvap_gas(self, gas_name, t, p):
        # Stub: Replace with actual vapor pressure calculation
        return 1.0

    def qvs_below(self, p, t, factor, gas_name, dtdlnp):
        # Stub: Replace with actual qvs_below calculation
        return 1.0

    def layer(self, nsub_max, gas_name, grav, mw_atmos, kz_min, cloudf_min,
              gas_mw, rainf, rho_p, supsat, sig, cloudf, q_below,
              t_layer, p_layer, kz_layer, chf, t_top, t_base, p_top, p_base):
        # Stub: Replace with actual layer calculation
        # Returns qc_layer, qt_layer, rg_layer, reff_layer, ndz_layer, status_r, status_q
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0

    def eddysed(self, grav, teff, kz_min, cloudf_min, nsub_max, supsat, mw_atmos, do_virtual,
                nz, z, z_top, p, p_top, t, t_top, chf,
                ngas, gas_name, gas_mmr, gas_mw, rho_p, sig, rainf):
        # Allocate output arrays
        kz = np.zeros(nz)
        qt = np.zeros((nz, ngas))
        qc = np.zeros((nz, ngas))
        ndz = np.zeros((nz, ngas))
        rg = np.zeros((nz, ngas))
        reff = np.zeros((nz, ngas))
        cloudf = np.zeros(nz)

        # Internal arrays
        qc_path = np.zeros(ngas)
        q_below = np.array(gas_mmr).copy()
        report_status_r = [True] * ngas
        report_status_q = [True] * ngas
        cloud_base = [0] * ngas

        # Determine bottom/top indices and increment
        if z[1] > z[0]:
            ibot, itop, incr = 0, nz-1, 1
        else:
            ibot, itop, incr = nz-1, 0, -1

        # Loop over condensing gases
        for igas in range(ngas):
            t_bot = t_top[ibot-incr]
            p_bot = p_top[ibot-incr]
            qc_path[igas] = 0.0

            # Virtual layer adjustment at bottom
            if do_virtual:
                qvs_factor = (supsat+1)*gas_mw[igas]/mw_atmos
                pvap = self.pvap_gas(gas_name[igas], t_bot, p_bot)
                qvs = qvs_factor * pvap / p_bot

                if qvs <= q_below[igas]:
                    p_lo = p_bot
                    p_hi = p_bot * 1e2
                    delta_q = q_below[igas] / 1e2
                    dtdlnp = (t_top[ibot] - t_bot) / np.log(p_bot/p_top[ibot])

                    # Find cloud base pressure (stub)
                    p_base = p_bot * 0.9  # Replace with root finding
                    t_base = t_bot + np.log(p_bot/p_base) * dtdlnp

                    # Virtual layer properties
                    p_layer = 0.5 * (p_bot + p_base)
                    t_layer = t_bot + np.log(p_bot/p_layer) * dtdlnp

                    qc_layer, qt_layer, rg_layer, reff_layer, ndz_layer, status_r, status_q = self.layer(
                        nsub_max, gas_name[igas], grav, mw_atmos, kz_min, cloudf_min,
                        gas_mw[igas], rainf[ibot, igas], rho_p[igas], supsat, sig[ibot, igas],
                        cloudf[ibot], q_below[igas], t_layer, p_layer, 0, chf[ibot],
                        t_bot, t_base, p_bot, p_base
                    )

                    print(f"\neddysed(): condensing gas = {gas_name[igas]}")
                    print(f" cloud base at p, T = {p_base/1e6:.3f}, {t_base:.3f}")
                    print(f" optical depth below domain = {1.5*qc_layer*(p_base-p_bot)/grav/(rho_p[igas]*reff_layer):.3e}\n")

            # Loop over atmospheric layers
            for iz in range(ibot, itop+incr, incr):
                qc_layer, qt_layer, rg_layer, reff_layer, ndz_layer, status_r, status_q = self.layer(
                    nsub_max, gas_name[igas], grav, mw_atmos, kz_min, cloudf_min,
                    gas_mw[igas], rainf[iz, igas], rho_p[igas], supsat, sig[iz, igas],
                    cloudf[iz], q_below[igas], t[iz], p[iz], kz[iz], chf[iz],
                    t_top[iz], t_top[iz-incr], p_top[iz], p_top[iz-incr]
                )
                qc[iz, igas] = qc_layer
                qt[iz, igas] = qt_layer
                rg[iz, igas] = rg_layer
                reff[iz, igas] = reff_layer
                ndz[iz, igas] = ndz_layer

                # Accumulate vertical path
                qc_path[igas] += qc_layer * (p_top[iz-incr] - p_top[iz]) / grav

            print(f"\neddysed(): condensing gas = {gas_name[igas]}")
            print(f" condensate path = {qc_path[igas]*1e4:.3f} g/m^2\n")

        # Return all output arrays
        return {
            "kz": kz, "qt": qt, "qc": qc, "ndz": ndz,
            "rg": rg, "reff": reff, "cloudf": cloudf
        }