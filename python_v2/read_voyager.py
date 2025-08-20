import numpy as np

# Stefan–Boltzmann constant (cgs units: erg cm^-2 s^-1 K^-4)
STEFBOLTZ = 5.670374419e-5  

def read_voyager(input_file, mw_atmos=2.2):
    """
    Read Voyager/Galileo-style atmospheric profile file.

    Parameters
    ----------
    input_file : str
        Path to input file with teff, grav, nzp1, and profile data.
    mw_atmos : float, optional
        Mean molecular weight of the atmosphere (amu).
        Default = 2.2 (Jupiter-like).

    Returns
    -------
    teff : float
        Effective temperature [K]
    grav : float
        Gravitational acceleration [cm/s^2]
    nz : int
        Number of model layers
    z : ndarray
        Altitude at pressure midpoints [cm]
    z_top : ndarray
        Altitude at layer tops [cm]
    p : ndarray
        Pressure at layer midpoints [dyne/cm^2]
    p_top : ndarray
        Pressure at layer tops [dyne/cm^2]
    t : ndarray
        Temperature at midpoints [K]
    t_top : ndarray
        Temperature at tops [K]
    chf : ndarray
        Convective heat flux [erg cm^-2 s^-1]
    """

    with open(input_file, "r") as f:
        teff = float(f.readline().strip())         # K
        grav = float(f.readline().strip())         # cm/s^2
        nzp1 = int(f.readline().strip())           # number of layer edges

        # skip one header line
        _ = f.readline()

        # read first top-of-atmosphere entry
        t_in, p_in = map(float, f.readline().split())
        t_top = [t_in]
        p_top = [p_in * 1e6]  # convert bar -> dyne/cm^2

        # remaining tops
        for _ in range(nzp1 - 1):
            t_in, p_in = map(float, f.readline().split())
            t_top.append(t_in)
            p_top.append(p_in * 1e6)

    t_top = np.array(t_top)
    p_top = np.array(p_top)

    # number of layers
    nz = nzp1 - 1

    # Gas constant per unit mass for atmosphere (erg g^-1 K^-1)
    R_atmos = 8.3143e7 / mw_atmos  

    # Allocate arrays
    z_top = np.zeros(nz + 1)
    z = np.zeros(nz)
    p = np.zeros(nz)
    t = np.zeros(nz)

    # bottom boundary condition
    z_top[nz] = 0.0  

    # integrate hydrostatically upward
    for iz in range(nz - 1, -1, -1):  # from bottom (nz-1) to top (0)
        itop = iz
        ibot = iz + 1
        dlnp = np.log(p_top[ibot] / p_top[itop])

        # midpoint pressure
        p[iz] = 0.5 * (p_top[itop] + p_top[ibot])

        scale_h = R_atmos * t_top[ibot] / grav
        dz_layer = scale_h * dlnp
        z_top[itop] = z_top[ibot] + dz_layer

        # lapse rate
        dtdz = (t_top[itop] - t_top[ibot]) / dz_layer
        dz_pmid = scale_h * np.log(p_top[ibot] / p[iz])
        z[iz] = z_top[ibot] + dz_pmid
        t[iz] = t_top[ibot] + dtdz * dz_pmid

    # convective heat flux (same for all layers)
    chf = np.full(nz, STEFBOLTZ * teff**4)

    return teff, grav, nz, z, z_top, p, p_top, t, t_top, chf
