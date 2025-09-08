# main.py
import numpy as np
from eddysed_lh import eddysed_nh3, pvap_nh3

# Constants / defaults (cgs units)
R_GAS = 8.314462618e7  # erg/mol/K

# User-changeable model parameters
INPUT_FILENAME = "/home/sansar1/Codes/cloud_model/input/profiles/voyager.input"
OUTPUT_FILENAME = "voyager.output"

# Physical defaults
GRAV = 980.0            # cm/s^2
TEFF = 124.4             # K (not used in this minimal NH3-only example)
MW_ATMOS = 2.3           # g/mol, H2-rich
SIG_ALL = 2.0
RAINF_ALL = 1
GAS_NAME = "NH3"         # only NH3
MW_CLOUD = 17.0
GAS_MMR = 1.34e-4 * (MW_CLOUD / MW_ATMOS)  # g/g
RHO_P = 0.84             # g/cm^3

def read_voyager_input(fname):
    """Read mid-layer T (K) and P (bar) from input file."""
    with open(fname, "r") as f:
        lines = f.readlines()
    data_lines = lines[4:]  # skip 4 header lines
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
        P.append(float(parts[1].replace("E", "e")))
    return np.array(T), np.array(P)

def build_layer_edges(p_mid):
    """Geometric mean for top/bottom layer edges."""
    nz = len(p_mid)
    p_top = np.zeros(nz + 1)
    if nz == 1:
        p_top[0] = p_mid[0] * 1.1
        p_top[1] = p_mid[0] / 1.1
        return p_top
    for i in range(1, nz):
        p_top[i] = 0.5*(p_mid[i-1] + p_mid[i])
    p_top[0] = p_mid[0] + 0.5*(p_mid[0]-p_mid[1])
    p_top[-1] = p_mid[-1] - 0.5*(p_mid[-1]-p_mid[-2])
    return p_top

def main():
    # Read profile
    T_mid, P_mid_bar = read_voyager_input(INPUT_FILENAME)
    nz = len(T_mid)
    if nz == 0:
        raise RuntimeError("No profile levels read from input file.")

    # Convert pressures to dyne/cm^2
    p_mid = P_mid_bar * 1e6

    # Build dummy z array for eddysed_nh3
    z_dummy = np.zeros(nz)

    # Compute dlnp
    dlnp = np.zeros(nz)
    dlnp[:-1] = np.log(p_mid[1:] / p_mid[:-1])
    dlnp[-1] = dlnp[-2]

    # Placeholder arrays
    sig = np.full(nz, SIG_ALL)
    rainf = np.full(nz, RAINF_ALL)
    chf = np.full(nz, 1e5)  # erg/cm^2/s
    lapse_ratio = np.ones(nz)

    # Call NH3 eddysed solver
    kz_out, lhf_out, qt_out, qc_out, ndz_out, rg_out, reff_out = eddysed_nh3(
        GRAV,
        TEFF,
        sig,
        rainf,
        nz,
        z_dummy,
        p_mid,
        T_mid,
        dlnp,
        chf,
        lapse_ratio,
        GAS_MMR,
        MW_ATMOS
    )

    # Output results
    with open(OUTPUT_FILENAME, "w") as fo:
        fo.write("# Pressure_bar    qc_g_per_g  (NH3)\n")
        for i in range(nz):
            fo.write(f"{P_mid_bar[i]:14.6e} {qc_out[i]:14.6e}\n")

    print(f"Wrote {OUTPUT_FILENAME} with {nz} levels (pressure in bar, qc in g/g).")

if __name__ == "__main__":
    main()
