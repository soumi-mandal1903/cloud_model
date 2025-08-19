import numpy as np

# Import your previously defined classes
# from optics_initializer import OpticsInitializer
# from eddysed_model import EddySedimentationModel
# from optics_calculator import OpticsCalculator
# from root_finder import RootFinderLog

class MainModel:
    def __init__(self):
        # Constants
        self.R_GAS = 8.3143e7
        self.MAXNGAS = 10
        self.MAXNZ = 100
        self.MAXNWAVE = 196
        self.MAXNRAD = 100
        self.NCASE = 20

    def pvap_ch4(self, T): return 1.0
    def pvap_nh3(self, T): return 1.0
    def pvap_h2o(self, T): return 1.0
    def pvap_fe(self, T): return 1.0
    def pvap_kcl(self, T): return 1.0
    def pvap_mgsio3(self, T): return 1.0
    def pvap_al2o3(self, T): return 1.0
    def pvap_mg2sio4(self, T): return 1.0
    def pvap_mns(self, T): return 1.0
    def pvap_zns(self, T): return 1.0
    def pvap_na2s(self, T): return 1.0
    def pvap_cr(self, T): return 1.0

    def read_goyal(self):
        # Stub: Replace with actual file reading logic
        # Returns grav, teff, nz, z, z_top, p, p_top, t, t_top, chf
        nz = 50
        grav = 980.
        teff = 130.
        z = np.linspace(0, 1000, nz)
        z_top = np.linspace(0, 1000, nz+1)
        p = np.linspace(1e6, 1e3, nz)
        p_top = np.linspace(1e6, 1e3, nz+1)
        t = np.linspace(130, 100, nz)
        t_top = np.linspace(130, 100, nz+1)
        chf = np.ones(nz)
        return grav, teff, nz, z, z_top, p, p_top, t, t_top, chf

    def run(self):
        # Model parameters
        mw_atmos = 2.2
        ngas = 7
        gas_name = ['Al2O3', 'Fe', 'Na2S', 'KCl', 'Cr', 'Mg2SiO4', 'ZnS']
        gas_mw = np.zeros(self.MAXNGAS)
        gas_mmr = np.zeros(self.MAXNGAS)
        rho_p = np.zeros(self.MAXNGAS)

        # Assign molecular weights and densities
        for igas, name in enumerate(gas_name):
            if name == 'CH4':
                gas_mw[igas] = 16.
                rho_p[igas] = 0.49
            elif name == 'NH3':
                gas_mw[igas] = 17.
                rho_p[igas] = 0.84
            elif name == 'H2O':
                gas_mw[igas] = 18.
                rho_p[igas] = 0.93
            elif name == 'Fe':
                gas_mw[igas] = 55.845
                rho_p[igas] = 7.875
            elif name == 'KCl':
                gas_mw[igas] = 74.5
                rho_p[igas] = 1.99
            elif name == 'MgSiO3':
                gas_mw[igas] = 100.4
                rho_p[igas] = 3.192
            elif name == 'Mg2SiO4':
                gas_mw[igas] = 140.7
                rho_p[igas] = 3.214
            elif name == 'Al2O3':
                gas_mw[igas] = 101.961
                rho_p[igas] = 3.987
            elif name == 'Na2S':
                gas_mw[igas] = 78.05
                rho_p[igas] = 1.856
            elif name == 'MnS':
                gas_mw[igas] = 87.00
                rho_p[igas] = 4.0
            elif name == 'ZnS':
                gas_mw[igas] = 97.46
                rho_p[igas] = 4.04
            elif name == 'Cr':
                gas_mw[igas] = 51.996
                rho_p[igas] = 7.15
            else:
                raise ValueError(f"main(): bad gas_name = {name}")

        # Assign mixing ratios (stub values, fill as needed)
        for igas, name in enumerate(gas_name):
            gas_mmr[igas] = 1e-6 * (gas_mw[igas] / mw_atmos)

        # Model parameters for dynamics and microphysics
        kz_min = 1e5
        sig_all = 2.0
        rainf_all = 0.1
        cloudf_min = 0.75
        nsub_max = 64
        do_virtual = True
        do_subcloud = True
        supsat = 0
        do_cases = False
        rainf_case = np.array([
            0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
            1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0
        ])
        ncases = self.NCASE if do_cases else 1

        # Optics setup
        do_optics = True
        read_mie = True
        optics_init = OpticsInitializer()
        optics_result = optics_init.init_optics(do_optics, read_mie, ngas, gas_name)
        nwave = optics_result['nwave']
        nrad = optics_result['nrad']
        wave = optics_result['wave']
        radius = optics_result['radius']
        dr = optics_result['dr']
        qscat = optics_result['qscat']
        qext = optics_result['qext']
        cos_qscat = optics_result['cos_qscat']

        # Read atmospheric profiles
        grav, teff, nz, z, z_top, p, p_top, t, t_top, chf = self.read_goyal()
        lapse_ratio = np.ones(nz)

        # Allocate arrays
        sig = np.zeros((nz, ngas))
        rainf = np.zeros((nz, ngas))
        cloudf = np.zeros(nz)
        cloudf_gas = np.zeros(ngas)
        kz = np.zeros(nz)
        qt = np.zeros((nz, ngas))
        qc = np.zeros((nz, ngas))
        ndz = np.zeros((nz, ngas))
        rg = np.zeros((nz, ngas))
        reff = np.zeros((nz, ngas))
        opd = np.zeros((nz, nwave))
        w0 = np.zeros((nz, nwave))
        g0 = np.zeros((nz, nwave))
        opd_gas = np.zeros((nz, ngas))
        dewp = np.zeros((nz, ngas))
        rainf_case_arr = np.zeros(self.NCASE)
        reff_qc = np.zeros((self.NCASE, ngas))
        reff_lwp = np.zeros((self.NCASE, ngas))
        totpath_gas = np.zeros((self.NCASE, ngas))
        opd_case = np.zeros((self.NCASE, ngas))

        # Main case loop
        for icase in range(ncases):
            # Set sig and rainf for each layer/gas
            for igas in range(ngas):
                for iz in range(nz):
                    sig[iz, igas] = sig_all
                    rainf[iz, igas] = rainf_case[icase] if do_cases else rainf_all

            # Eddy sedimentation
            eddysed_model = EddySedimentationModel()
            eddy_result = eddysed_model.eddysed(
                grav, teff, kz_min, cloudf_min, nsub_max, supsat, mw_atmos, do_virtual,
                nz, z, z_top, p, p_top, t, t_top, chf,
                ngas, gas_name, gas_mmr, gas_mw, rho_p, sig, rainf,
            )
            kz = eddy_result['kz']
            qt = eddy_result['qt']
            qc = eddy_result['qc']
            ndz = eddy_result['ndz']
            rg = eddy_result['rg']
            reff = eddy_result['reff']
            cloudf = eddy_result['cloudf']

            # Optical properties
            optics_calc = OpticsCalculator()
            optics_out = optics_calc.calc_optics(
                do_subcloud, nz, ngas, nwave, nrad,
                z, gas_name, radius, dr, qscat, qext, cos_qscat,
                ndz, sig, rg
            )
            opd = optics_out['opd']
            w0 = optics_out['w0']
            g0 = optics_out['g0']
            opd_gas = optics_out['opd_gas']
            opd_tot = optics_out['opd_tot']

            # Diagnostics
            if z[1] > z[0]:
                itop = nz - 1
                ibot = 0
                incr = 1
            else:
                itop = 0
                ibot = nz - 1
                incr = -1

            for igas in range(ngas):
                rhoc_max = 0.
                lwp = 0.
                totpath = 0.
                opd_eff = 0.
                iz_max = ibot

                for iz in range(itop, ibot + incr, incr):
                    lwp_layer = qc[iz, igas] * (p_top[iz + incr] - p_top[iz]) / grav
                    lwp += lwp_layer
                    totpath += (qc[iz, igas] + qt[iz, igas]) * (p_top[iz + incr] - p_top[iz]) / grav
                    if reff[iz, igas] > 0.:
                        opd_eff += lwp_layer / reff[iz, igas]
                    rhoc = qc[iz, igas] * p[iz] / (self.R_GAS / mw_atmos * t[iz])
                    if rhoc > rhoc_max:
                        rhoc_max = rhoc
                        iz_max = iz

                reff_qc[icase, igas] = reff[iz_max, igas] * 1e4
                opd_case[icase, igas] = opd_gas[ibot, igas]
                if opd_eff > 0.:
                    reff_lwp[icase, igas] = lwp / opd_eff * 1e4
                    cloudf_gas[igas] = cloudf[iz_max]
                else:
                    reff_lwp[icase, igas] = 0.
                    cloudf_gas[igas] = 0.
                totpath_gas[icase, igas] = totpath

                # Dewpoint calculation (stubbed, fill with root finding as needed)
                for iz in range(nz):
                    eps = gas_mw[igas] / mw_atmos
                    rhs = p[iz] * gas_mmr[igas] / (eps + gas_mmr[igas])
                    delta_p = rhs / 100.
                    tlo, thi = 100., 400.
                    dewp[iz, igas] = 0.0  # Replace with root finding for pvap_xxx

            # Output (terminal and file)
            do_ttyo = not do_cases
            do_fileo = True

            if do_ttyo:
                print('\nmain():')
                print(f'Teff = {int(teff + 0.5)} K')
                print(f'grav = {int(grav / 1e2 + 0.5)} m/s^2')
                print(f'nz = {int(nz)}')
                for igas in range(ngas):
                    print(f'\ncondensing gas = {gas_name[igas]}')
                    print('iz  P(bar)   z(m)     T(K) kz(cm^2/s) reff(um) qc(g/g) col_opd cloudf dewp(K)')
                    for iz in range(itop, ibot + incr, incr):
                        print(f"{iz:3d} {p[iz]/1e6:9.1f} {z[iz]/1e2:9.2f} {t[iz]:7.1f} {kz[iz]:12.1f} "
                              f"{reff[iz, igas]*1e4:9.2f} {qc[iz, igas]:9.2e} {opd_gas[iz, igas]:9.2e} "
                              f"{cloudf[iz]:6.2f} {dewp[iz, igas]:7.2f}")
                    print(f'\n geometric scattering optical depth = {opd_gas[ibot, igas]:.3e}')
                    print(f' effective radius (um) at condensate conc. max. = {reff_qc[icase, igas]:.2f}')
                    print(f' effective radius (um) from lwp, opd = {reff_lwp[icase, igas]:.2f}')
                    print(f' fractional coverage at cloud base = {cloudf_gas[igas]:.2f}')
                print(f'\n total optical depth for all geometric scatterers = {opd_tot:.3e}')

            if do_fileo:
                with open('eddysed.out', 'w') as f:
                    f.write(f"{teff} {grav} {nz} {ngas} {sig_all} {rainf_all}\n")
                    for igas in range(ngas):
                        f.write(f"{gas_name[igas]}\n")
                        f.write(f"{gas_mmr[igas]} {opd_gas[ibot, igas]}\n")
                    for iz in range(nz):
                        f.write(f"{z[iz]} {p[iz]} {t[iz]} " +
                                " ".join(str(qt[iz, igas]) for igas in range(ngas)) + " " +
                                " ".join(str(qc[iz, igas]) for igas in range(ngas)) + " " +
                                " ".join(str(reff[iz, igas]) for igas in range(ngas)) + " " +
                                " ".join(str(rg[iz, igas]) for igas in range(ngas)) + " " +
                                " ".join(str(opd_gas[iz, igas]) for igas in range(ngas)) + " " +
                                " ".join(str(dewp[iz, igas]) for igas in range(ngas)) + "\n")

        # End of case loop

        if do_cases:
            igas = 0
            print("rainf    =", ", ".join(f"{rainf_case[icase]:8.2f}" for icase in range(self.NCASE)))
            print("reff lwp =", ", ".join(f"{reff_lwp[icase, igas]:8.2f}" for icase in range(self.NCASE)))
            print("opd      =", ", ".join(f"{opd_case[icase, igas]:8.2f}" for icase in range(self.NCASE)))
            print("totpath  =", ", ".join(f"{totpath_gas[icase, igas]:8.2f}" for icase in range(self.NCASE)))

# Example usage:
# main_model = MainModel()
# main_model.run()
#there are LARGE SECTIONS COMMENTED OUT IN THE ORIGINAL CODE,
# so the above code is a simplified version that focuses on the main structure and logic.