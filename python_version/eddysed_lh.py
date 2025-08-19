import numpy as np

# Constants from globals.h
R_GAS = 8.3143e7  # erg/mol/K
PI = 3.14159274
AVOGADRO = 6.02e23
K_BOLTZ = R_GAS / AVOGADRO

class EddySedLHModel:
    def __init__(self):
        pass

    # Vapor pressure stubs
    def pvap_ch4(self, T): return 1.0
    def pvap_nh3(self, T): return 1.0
    def pvap_h2o(self, T): return 1.0
    def pvap_fe(self, T): return 1.0
    def pvap_kcl(self, T): return 1.0
    def pvap_mgsio3(self, T): return 1.0

    # Microphysics stubs
    def advdiff(self, qt, qbelow, qvs, mixl_layer, zdiff, rainf_layer):
        # Replace with actual advdiff logic
        return qt - (qbelow + qvs) / 2

    def vfall(self, r, grav, mw_atmos, mfp, visc, t_layer, p_layer, rho_p, mw_cloud):
        # Replace with actual vfall logic
        return r * 1e3

    def find_root(self, func, target, lo, hi, tol, args=()):
        from scipy.optimize import root_scalar
        try:
            sol = root_scalar(lambda x: func(x, *args) - target, bracket=[lo, hi], xtol=tol, method='secant')
            if sol.converged:
                return sol.root, 0
            else:
                return None, 1
        except ValueError:
            return None, -1
        except Exception:
            return None, 2

    def eddysed(self, grav_in, teff, sig, rainf, 
                nz, z, p, t, dlnp, chf, lapse_ratio,
                ngas, gas_name, gas_mmr, mw_atmos_in):
        # Allocate output arrays
        kz = np.zeros(nz)
        lhf = np.zeros(nz)
        qt = np.zeros((nz, ngas))
        qc = np.zeros((nz, ngas))
        ndz = np.zeros((nz, ngas))
        rg = np.zeros((nz, ngas))
        reff = np.zeros((nz, ngas))

        # Local variables
        rw = np.zeros((nz, ngas))
        dz = np.zeros(nz)
        cloud_base = np.zeros(ngas, dtype=int)
        cond_path = np.zeros(ngas)

        # Physical constants
        d_molecule = 2.827e-8
        eps_k = 59.7
        grav = grav_in
        mw_atmos = mw_atmos_in
        r_atmos = R_GAS / mw_atmos
        c_p = 7. / 2. * r_atmos
        supsat = 0
        fs = supsat + 1

        # Determine bottom/top indices and increment
        if z[1] > z[0]:
            ibot, itop, incr = 0, nz-1, 1
        else:
            ibot, itop, incr = nz-1, 0, -1

        # Loop over layers
        for iz in range(ibot, itop+incr, incr):
            rho_atmos = p[iz] / (r_atmos * t[iz])
            scale_h = r_atmos * t[iz] / grav
            dz[iz] = scale_h * dlnp[iz]
            mixl = max(0.1, lapse_ratio[iz]) * scale_h
            n_atmos = p[iz] / (K_BOLTZ * t[iz])
            mfp = 1. / (np.sqrt(2.) * n_atmos * PI * d_molecule**2)
            visc = (5./16.) * np.sqrt(PI * K_BOLTZ * t[iz] * (mw_atmos / AVOGADRO)) / \
                   (PI * d_molecule**2) / (1.22 * (t[iz] / eps_k)**(-0.16))
            t_layer = t[iz]
            p_layer = p[iz]
            rlo = 1e-10
            rhi = 10.0
            lh_sum = 0.0

            # Loop over gases
            for igas in range(ngas):
                # Set vapor pressure, molecular weight, latent heat, density
                name = gas_name[igas]
                if name == 'CH4':
                    pvap = self.pvap_ch4(t[iz])
                    mw_cloud = 16.
                    lh = 0.0
                    rho_p = 0.49
                elif name == 'NH3':
                    pvap = self.pvap_nh3(t[iz])
                    mw_cloud = 17.
                    lh = 28.6e10
                    rho_p = 0.84
                elif name == 'H2O':
                    pvap = self.pvap_h2o(t[iz])
                    mw_cloud = 18.
                    lh = 46.7e10
                    rho_p = 0.93
                elif name == 'Fe':
                    pvap = self.pvap_fe(t[iz])
                    mw_cloud = 55.845
                    lh = 0.0
                    rho_p = 7.875
                elif name == 'KCl':
                    pvap = self.pvap_kcl(t[iz])
                    mw_cloud = 74.5
                    lh = 0.0
                    rho_p = 1.99
                elif name == 'MgSiO3':
                    pvap = self.pvap_mgsio3(t[iz])
                    mw_cloud = 100.4
                    lh = 0.0
                    rho_p = 3.192
                else:
                    raise ValueError(f"eddysed(): bad igas = {igas}")

                qvs = fs * pvap / ((R_GAS / mw_cloud) * t[iz]) / rho_atmos

                # Vapor mass mixing ratio in underlying layer
                if iz == ibot:
                    qbelow = gas_mmr[igas]
                    zdiff = dz[iz] / 2.
                else:
                    qbelow = qt[iz-incr, igas]
                    zdiff = z[iz] - z[iz-incr]

                # Cloud free layer
                if qbelow < qvs:
                    qt[iz, igas] = qbelow
                    qc[iz, igas] = 0.0
                    rg[iz, igas] = 0.0
                else:
                    # Cloudy layer: microphysics/dynamics solution
                    rainf_layer = rainf[iz, igas]
                    mixl_layer = mixl
                    qhi = qbelow
                    qlo = qhi / 1e3
                    delta_q = qbelow / 1000.

                    def adv_func(qt_val):
                        return self.advdiff(qt_val, qbelow, qvs, mixl_layer, zdiff, rainf_layer)

                    qt_val, status = self.find_root(adv_func, 0.0, qlo, qhi, delta_q)
                    qt[iz, igas] = qt_val if qt_val is not None else qbelow
                    qc[iz, igas] = max(0.0, qt[iz, igas] - qvs)

                    # Latent heat fluxes due to precipitation
                    if mw_cloud > 0 and c_p > 0 and r_atmos > 0:
                        lh_sum += rainf_layer * qc[iz, igas] * lh / (mw_cloud * c_p / r_atmos)

                # Optical properties for cloudy layers
                if qc[iz, igas] > 0.0:
                    if cloud_base[igas] == 0:
                        cloud_base[igas] = iz

                    delta_v = 1e-3 if lh_sum == 0 else lh_sum / 1000.
                    def vfall_func(r):
                        return self.vfall(r, grav, mw_atmos, mfp, visc, t_layer, p_layer, rho_p, mw_cloud)
                    rw_val, status = self.find_root(vfall_func, lh_sum, rlo, rhi, delta_v)
                    rw[iz, igas] = rw_val if rw_val is not None else 0.0

                    lnsig2 = 0.5 * np.log(sig[iz, igas]) ** 2
                    sig_alpha = max(1.1, sig[iz, igas])
                    if rainf[iz, igas] > 1:
                        alpha = np.log(self.vfall(rw[iz, igas] * sig_alpha, grav, mw_atmos, mfp, visc, t_layer, p_layer, rho_p, mw_cloud) / lh_sum) / np.log(sig_alpha)
                    else:
                        alpha = np.log(lh_sum / self.vfall(rw[iz, igas] / sig_alpha, grav, mw_atmos, mfp, visc, t_layer, p_layer, rho_p, mw_cloud)) / np.log(sig_alpha)

                    rg[iz, igas] = rainf[iz, igas] ** (1.0 / alpha) * rw[iz, igas] * np.exp(-(alpha + 6) * lnsig2)
                    reff[iz, igas] = rg[iz, igas] * np.exp(5 * lnsig2)
                    ndz[iz, igas] = 3 * rho_atmos * qc[iz, igas] * dz[iz] / (4 * PI * rho_p * rg[iz, igas] ** 3) * np.exp(-9 * lnsig2)
                else:
                    rg[iz, igas] = 0.0
                    reff[iz, igas] = 0.0
                    ndz[iz, igas] = 0.0

            # Convective velocity scale (cm/s)
            coeff_a = lh_sum
            coeff_b = -chf[iz] / (rho_atmos * c_p / r_atmos)
            if lh_sum == 0.0:
                w_convect = (-coeff_b) ** (1. / 3.)
            else:
                cubic_arg = coeff_b ** 2 / 4. + coeff_a ** 3 / 27.
                w_convect = ((-coeff_b / 2. + np.sqrt(cubic_arg)) ** (1. / 3.) +
                             (-coeff_b / 2. - np.sqrt(cubic_arg)) ** (1. / 3.))

            lhf[iz] = w_convect * lh_sum * rho_atmos * c_p / r_atmos
            kz[iz] = w_convect * mixl

        # Condensate path
        for igas in range(ngas):
            cond_path[igas] = np.sum(qc[:, igas] * dlnp * p / grav)
            # Print or log as needed

        # Return all output arrays
        return {
            "kz": kz,
            "lhf": lhf,
            "qt": qt,
            "qc": qc,
            "ndz": ndz,
            "rg": rg,
            "reff": reff,
            "cond_path": cond_path
        }