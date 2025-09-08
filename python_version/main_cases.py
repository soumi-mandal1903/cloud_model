import numpy as np

# Import your previously defined classes
# from optics_initializer import OpticsInitializer
# from eddysed_model import EddySedimentationModel
# from optics_calculator import OpticsCalculator

class MainCases:
    def __init__(self):
        # Constants
        self.R_GAS = 8.3143e7
        self.MAXNGAS = 10
        self.MAXNZ = 100
        self.MAXNWAVE = 196
        self.MAXNRAD = 100
        self.NCASE = 20

    def read_voyager(self, input_fname):
        # Stub: Replace with actual file reading logic
        # Returns grav, teff, nz, z, z_top, p, p_top, t, t_top, dlnp, chf, lapse_ratio
        nz = 50
        grav = 980.
        teff = 130.
        z = np.linspace(0, 1000, nz)
        z_top = np.linspace(0, 1000, nz+1)
        p = np.linspace(1e6, 1e3, nz)
        p_top = np.linspace(1e6, 1e3, nz+1)
        t = np.linspace(130, 100, nz)
        t_top = np.linspace(130, 100, nz+1)
        dlnp = np.ones(nz)
        chf = np.ones(nz)
        lapse_ratio = np.ones(nz)
        return grav, teff, nz, z, z_top, p, p_top, t, t_top, dlnp, chf, lapse_ratio

    def run(self):
        # Model parameters
        mw_atmos = 2.2
        ngas = 1
        gas_name = ['NH3']
        gas_mmr = np.zeros(self.MAXNGAS)
        gas_mmr[0] = 3.0e-5 * (17. / mw_atmos)

        kz_min = 1e5
        sig_all = 2.0
        rainf_all = 3.0
        cloudf_min = 0.75
        nsub_max = 64
        do_subcloud = False
        do_cases = True

        rainf_case = np.array([
            0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
            1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0
        ])
        ncases = self.NCASE if do_cases else 1

        # Optics setup
        do_optics =   True
        read_mie = False
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
        input_fname = 'input/profiles/voyager.input'
        grav, teff, nz, z, z_top, p, p_top, t, t_top, dlnp, chf, lapse_ratio = self.read_voyager(input_fname)

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
        rainf_case_arr = np.zeros(self.NCASE)
        reff_qc = np.zeros((self.NCASE, ngas))
        reff_lwp = np.zeros((self.NCASE, ngas))
        opd_case = np.zeros((self.NCASE, ngas))
        totpath_gas = np.zeros((self.NCASE, ngas))

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
                grav, teff, kz_min, cloudf_min, nsub_max, sig, rainf,
                nz, z, z_top, p, p_top, t, t_top, dlnp, chf, lapse_ratio,
                ngas, gas_name, gas_mmr, mw_atmos
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

            # Output (terminal and file)
            do_ttyo = not do_cases
            do_fileo = True

            if do_ttyo:
                print('\nmain():')
                print(f'Teff = {int(teff + 0.5)} K')
                print(f'grav = {int(grav / 1e2 + 0.5)} m/s^2')
                for igas in range(ngas):
                    print(f'\ncondensing gas = {gas_name[igas]}')
                    print(' P(bar)   z(m)     T(K) kz(cm^2/s) reff(um) qc(g/g) col_opd cloudf')
                    for iz in range(itop, ibot + incr, incr):
                        print(f"{p[iz]/1e6:9.1f} {z[iz]/1e2:9.2f} {t[iz]:7.1f} {kz[iz]:12.1f} "
                              f"{reff[iz, igas]*1e4:9.2f} {qc[iz, igas]:9.2e} {opd_gas[iz, igas]:9.2e} "
                              f"{cloudf[iz]:6.2f} {lapse_ratio[iz]:6.2f}")
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
                                " ".join(str(opd_gas[iz, igas]) for igas in range(ngas)) + "\n")

        # End of case loop

        if do_cases:
            igas = 0
            print("rainf    =", ", ".join(f"{rainf_case[icase]:8.2f}" for icase in range(self.NCASE)))
            print("reff lwp =", ", ".join(f"{reff_lwp[icase, igas]:8.2f}" for icase in range(self.NCASE)))
            print("opd      =", ", ".join(f"{opd_case[icase, igas]:8.2f}" for icase in range(self.NCASE)))
            print("totpath  =", ", ".join(f"{totpath_gas[icase, igas]:8.2f}" for icase in range(self.NCASE)))

# Example usage:
# main_cases = MainCases()
# main_cases.run()