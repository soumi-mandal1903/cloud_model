import numpy as np

class OpticsCalculator:
    def __init__(self, PI=3.14159274):
        self.PI = PI

    def calc_optics(self, do_subcloud, nz, ngas, nwave, nrad,
                    z, gas_name, radius, dr, qscat, qext, cos_qscat,
                    ndz, sig, rg):
        # Output arrays
        opd = np.zeros((nz, nwave))
        w0 = np.zeros((nz, nwave))
        g0 = np.zeros((nz, nwave))
        opd_gas = np.zeros((nz, ngas))
        opd_layer = np.zeros((nz, ngas))
        scat_gas = np.zeros((nz, nwave, ngas))
        ext_gas = np.zeros((nz, nwave, ngas))
        cqs_gas = np.zeros((nz, nwave, ngas))

        # Determine indices of top and bottom layers
        if z[1] > z[0]:
            itop, ibot, incr = nz-1, 0, -1
        else:
            itop, ibot, incr = 0, nz-1, 1

        # Initialize indices of bottoms of cloud layers for subcloud kludge
        ibot_cloud = np.full(ngas, ibot)
        if do_subcloud:
            ibot_cloud[:] = ibot

        # Loop over layers and condensing vapors
        for iz in range(itop, ibot+incr, incr):
            for igas in range(ngas):
                opd_layer[iz, igas] = 0.0
                scat_gas[iz, :, igas] = 0.0
                ext_gas[iz, :, igas] = 0.0
                cqs_gas[iz, :, igas] = 0.0

                if ndz[iz, igas] > 0.0:
                    r2 = rg[iz, igas]**2 * np.exp(2 * np.log(sig[iz, igas])**2)
                    opd_layer[iz, igas] = 2.0 * self.PI * r2 * ndz[iz, igas]

                    rsig = sig[iz, igas]
                    norm = 0.0
                    for irad in range(nrad):
                        rr = radius[irad, igas]
                        arg1 = dr[irad, igas] / (np.sqrt(2.0 * self.PI) * rr * np.log(rsig))
                        arg2 = -np.log(rr / rg[iz, igas])**2 / (2 * np.log(rsig)**2)
                        norm += arg1 * np.exp(arg2)
                    norm = ndz[iz, igas] / norm if norm != 0 else 0.0

                    for irad in range(nrad):
                        rr = radius[irad, igas]
                        arg1 = dr[irad, igas] / (np.sqrt(2.0 * self.PI) * np.log(rsig))
                        arg2 = -np.log(rr / rg[iz, igas])**2 / (2 * np.log(rsig)**2)
                        pir2ndz = norm * self.PI * rr * arg1 * np.exp(arg2)
                        for iwave in range(nwave):
                            scat_gas[iz, iwave, igas] += qscat[iwave, irad, igas] * pir2ndz
                            ext_gas[iz, iwave, igas] += qext[iwave, irad, igas] * pir2ndz
                            cqs_gas[iz, iwave, igas] += cos_qscat[iwave, irad, igas] * pir2ndz

                    if do_subcloud:
                        ibot_cloud[igas] = iz

        # Subcloud kludge to soften discontinuity at cloud base
        if do_subcloud:
            for igas in range(ngas):
                ibot_gas = ibot_cloud[igas]
                if ibot_gas != ibot:
                    if ibot_gas + incr == ibot:
                        ibot_subcloud = ibot_gas + incr
                    else:
                        ibot_subcloud = ibot_gas + 2 * incr
                    for iz in range(ibot_gas + incr, ibot_subcloud + incr, incr):
                        norm = 0.10 if iz == ibot_gas + incr else 0.05
                        opd_layer[iz, igas] += opd_layer[ibot_gas, igas] * norm
                        for iwave in range(nwave):
                            scat_gas[iz, iwave, igas] += scat_gas[ibot_gas, iwave, igas] * norm
                            ext_gas[iz, iwave, igas] += ext_gas[ibot_gas, iwave, igas] * norm
                            cqs_gas[iz, iwave, igas] += cqs_gas[ibot_gas, iwave, igas] * norm

        # Sum over gases and compute spectral optical depth profile etc
        for iz in range(itop, ibot+incr, incr):
            for iwave in range(nwave):
                opd_scat = np.sum(scat_gas[iz, iwave, :])
                opd_ext = np.sum(ext_gas[iz, iwave, :])
                cos_qs = np.sum(cqs_gas[iz, iwave, :])
                if opd_scat > 0.0:
                    opd[iz, iwave] = opd_ext
                    w0[iz, iwave] = opd_scat / opd_ext if opd_ext != 0 else 0.0
                    g0[iz, iwave] = cos_qs / opd_scat
                else:
                    opd[iz, iwave] = 0.0
                    w0[iz, iwave] = 0.0
                    g0[iz, iwave] = 0.0

        # Cumulative optical depths for conservative geometric scatterers
        opd_tot = 0.0
        for igas in range(ngas):
            opd_gas[itop, igas] = opd_layer[itop, igas]
            for iz in range(itop + incr, ibot + incr, incr):
                opd_gas[iz, igas] = opd_gas[iz - incr, igas] + opd_layer[iz, igas]
            opd_tot += opd_gas[ibot, igas]

        return {
            "opd": opd,
            "w0": w0,
            "g0": g0,
            "opd_gas": opd_gas,
            "opd_tot": opd_tot
        }