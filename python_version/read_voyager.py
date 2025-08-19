import numpy as np

class AtmosphereReader:
    STEFBOLTZ = 5.67051e-5  # erg/cm^2/s/K^4

    def __init__(self, mw_atmos=2.2):
        self.r_atmos = 8.3143e7 / mw_atmos

    def read_voyager(self, teff, grav, nzp1, t_top, p_top):
        """
        Calculate atmospheric profiles for Galileo probe sounding.
        Arguments:
            teff: effective temperature (K)
            grav: gravitational acceleration (cm/s^2)
            nzp1: number of layer edges
            t_top: array of temperatures at layer edges (K), length nzp1
            p_top: array of pressures at layer edges (bars), length nzp1
        Returns:
            nz: number of layers
            z: altitudes at layer midpoints (cm), length nz
            z_top: altitudes at layer edges (cm), length nzp1
            p: pressures at layer midpoints (dyne/cm^2), length nz
            t: temperatures at layer midpoints (K), length nz
            chf: convective heat flux (erg/cm^2/s), length nz
        """
        # Convert p_top from bars to dyne/cm^2
        p_top = np.array(p_top) * 1e6
        t_top = np.array(t_top)
        nz = nzp1 - 1

        z = np.zeros(nz)
        z_top = np.zeros(nzp1)
        p = np.zeros(nz)
        t = np.zeros(nz)
        chf = np.zeros(nz)

        z_top[nz] = 0.0

        for iz in range(nz - 1, -1, -1):
            itop = iz
            ibot = iz + 1
            dlnp = np.log(p_top[ibot] / p_top[itop])
            p[iz] = 0.5 * (p_top[itop] + p_top[ibot])
            scale_h = self.r_atmos * t_top[ibot] / grav
            dz_layer = scale_h * dlnp
            z_top[iz] = z_top[ibot] + dz_layer
            dtdz = (t_top[itop] - t_top[ibot]) / dz_layer
            dz_pmid = scale_h * np.log(p_top[ibot] / p[iz])
            z[iz] = z_top[ibot] + dz_pmid
            t[iz] = t_top[ibot] + dtdz * dz_pmid

        for iz in range(nz):
            chf[iz] = self.STEFBOLTZ * teff ** 4

        return nz, z, z_top, p, t, chf

# Example usage:
# reader = AtmosphereReader(mw_atmos=2.2)
# nz, z, z_top, p, t, chf = reader.read_voyager(teff, grav, nzp1, t_top, p_top)