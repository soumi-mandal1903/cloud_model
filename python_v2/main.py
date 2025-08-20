# main.py
import numpy as np
from math import sqrt

# Import your eddysed implementation (assumed to be in eddysed.py)
from eddysed import eddysed

# Constants / defaults (cgs units)
R_GAS = 8.314462618e7  # erg mol^-1 K^-1 (dyne*cm / mol K)

# User-changeable model parameters
INPUT_FILENAME = "/home/sansar1/Codes/cloud_model/input/profiles/voyager.input"
OUTPUT_FILENAME = "voyager.output"

# Physical defaults (change if you want)
GRAV = 980.0            # cm/s^2 (earth-like). Set to planet-appropriate value if needed.
TEFF = 124.4            # not used in this minimal script but kept for compatibility
MW_ATMOS = 2.3          # mean molecular weight of atmosphere (g/mol), H2-rich ~2.3
KZ_MIN = 1e5
CLOUDF_MIN = 0.75
N_SUB_MAX = 64
SUPSAT = 0.0
DO_VIRTUAL = True

SIG_ALL = 2.0
RAINF_ALL = 1

# NH3 properties (used here only)
GAS_NAME = ["NH3"]      # only one gas
NGAS = 1
MW_CLOUD = 17.0         # NH3 molecular weight (g/mol)
# Choose a plausible total mass mixing ratio below cloud (g/g).
# I used a typical value mentioned in comments earlier (you can change):
GAS_MMR = np.array([1.34e-4 * (MW_CLOUD / MW_ATMOS)])  # (g/g)
RHO_P = np.array([0.84])  # NH3 solid density (g/cm^3)

def read_voyager_input(fname):
    """
    Read a simple voyager.input format:
     - skip first 4 header lines
     - remaining lines: at least two columns: T  P
    Returns:
      T_mid (K) : 1D numpy array of temperatures
      P_mid (bar): 1D numpy array of pressures (bar)
    """
    with open(fname, "r") as f:
        lines = f.readlines()

    data_lines = lines[4:]
    T = []
    P = []
    for line in data_lines:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        T.append(float(parts[0]))
        # allow both 'E' and 'e'
        P.append(float(parts[1].replace("E", "e")))
    return np.array(T), np.array(P)


def build_layer_edges(p_mid):
    """
    Given mid-point pressures p_mid (bar), build p_top edges array length nz+1.
    We choose geometric means for interior edges, and extend edges at ends.
    Returns p_top (bar) array length nz+1.
    """
    nz = len(p_mid)
    p_top = np.zeros(nz + 1)
    if nz == 1:
        p_top[0] = p_mid[0] * 1.1
        p_top[1] = p_mid[0] / 1.1
        return p_top

    # interior edges: geometric mean between adjacent midpoints
    for i in range(1, nz):
        p_top[i] = 0.5*(p_mid[i - 1] + p_mid[i])
    # extrapolate top and bottom edges by geometric spacing
    p_top[0] = p_mid[0] + (p_mid[0] - p_mid[1]) * 0.5
    p_top[-1] = p_mid[-1] + (p_mid[-1] - p_mid[-2]) * (-0.5)
    return p_top


def main():
    # Read the voyager.input file: mid-point T (K) and P (bar)
    T_mid, P_mid_bar = read_voyager_input(INPUT_FILENAME)
    nz = len(T_mid)
    if nz == 0:
        raise RuntimeError("No profile levels read from input file.")

    # Convert pressures to dyne/cm^2 for the model (1 bar = 1e6 dyne/cm^2)
    p_mid = P_mid_bar * 1e6  # dyne/cm^2

    # Build p_top edges (bar -> convert to dyne too)
    p_top_bar = build_layer_edges(P_mid_bar)
    p_top = p_top_bar * 1e6  # dyne/cm^2

    # Build t_top array (length nz+1). Use simple linear/exponential interpolation:
    t_top = np.zeros(nz + 1)
    # interior top temps: average neighbor mid temps
    t_top[1:-1] = 0.5 * (T_mid[:-1] + T_mid[1:])
    # extrapolate ends using nearest-midpoint
    t_top[0] = T_mid[0]
    t_top[-1] = T_mid[-1]

    # z and z_top: not used in our simplified eddysed, but construct dummy altitudes
    # Use scale-height estimate: z = H * ln(p_ref/p)
    # Choose p_ref = p_mid[0]
    H = (R_GAS / MW_ATMOS) * T_mid.mean() / GRAV
    z_mid = H * np.log(p_mid[0] / p_mid)        # cm
    # z_top boundaries from pressures p_top
    p_top_array = p_top
    z_top = H * np.log(p_mid[0] / p_top_array)  # cm

    # chf (convective heat flux) - use small non-zero values to avoid zeros
    chf = np.full(nz, 1e5)  # erg/s/cm^2 (placeholder)

    # initialize arrays required by eddysed
    z = z_mid.copy()
    p = p_mid.copy()
    t = T_mid.copy()

    # allocate 2D arrays shaped (nz, ngas)
    sig = np.full((nz, NGAS), SIG_ALL)
    rainf = np.full((nz, NGAS), RAINF_ALL)
    kz = np.zeros(nz)
    qt = np.zeros((nz, NGAS))
    qc = np.zeros((nz, NGAS))
    ndz = np.zeros((nz, NGAS))
    rg = np.zeros((nz, NGAS))
    reff = np.zeros((nz, NGAS))
    cloudf = np.zeros(nz)

    # Call eddysed for NH3 only (it expects all arrays, we pass NGAS=1 arrays)
    kz_out, qt_out, qc_out, ndz_out, rg_out, reff_out, cloudf_out = eddysed(
        GRAV,
        TEFF,
        KZ_MIN,
        CLOUDF_MIN,
        N_SUB_MAX,
        SUPSAT,
        MW_ATMOS,
        DO_VIRTUAL,
        nz,
        z,
        z_top,
        p,
        p_top,
        t,
        t_top,
        chf,
        NGAS,
        GAS_NAME,
        GAS_MMR,
        np.array([MW_CLOUD]),
        RHO_P,
        sig,
        rainf,
    )

    # qc_out is shape (nz, 1) -> grab column 0
    qc_nh3 = qc_out[:, 0]

    # Output: write pressure (bar) and qc values to voyage.output
    with open(OUTPUT_FILENAME, "w") as fo:
        fo.write("# Pressure_bar    qc_g_per_g  (NH3)\n")
        for i in range(nz):
            fo.write(f"{P_mid_bar[i]:14.6e} {qc_nh3[i]:14.6e}\n")

    print(f"Wrote {OUTPUT_FILENAME} with {nz} levels (pressure in bar, qc in g/g).")


if __name__ == "__main__":
    main()
