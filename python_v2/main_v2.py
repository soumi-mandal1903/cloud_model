#!/usr/bin/env python3
import numpy as np
from eddysed import eddysed
from read_voyager import read_voyager

def main():
    # --------------------------------------------------
    # Read input atmosphere from Voyager/Galileo-style file
    # --------------------------------------------------
    input_file = "/home/sansar1/Codes/cloud_model/input/profiles/voyager.input"   # change to your actual file
    teff, grav, nz, z, z_top, p, p_top, t, t_top, chf = read_voyager(input_file)

    # --------------------------------------------------
    # Define other model inputs
    # --------------------------------------------------
    kz = np.full(nz, 1e6)   # Eddy diffusion coefficient (cm^2/s)

    ngas = 1
    gas_name = ["nh3"]
    gas_mmr = [0.001]       # Example value
    mw_atmos = 2.2

    lapse_ratio = np.ones(nz)
    lhf = np.zeros(nz)

    qt = np.zeros(nz)
    qc = np.zeros(nz)
    ndz = np.zeros(nz)
    rg = np.zeros(nz)
    reff = np.zeros(nz)
    cloudf = np.zeros(nz)

    # --------------------------------------------------
    # Run eddy sedimentation model
    # --------------------------------------------------
    kz_out, qt_out, qc_out, ndz_out, rg_out, reff_out, cloudf_out = eddysed(
        grav, teff,
        sig=5.67e-5,
        rainf=1.0,
        nz=nz,
        z=z,
        p=p,
        t=t,
        dlnp=np.gradient(np.log(p)),
        chf=chf,
        lapse_ratio=lapse_ratio,
        ngas=ngas,
        gas_name=gas_name,
        gas_mmr=gas_mmr,
        mw_atmos_in=mw_atmos,
        kz=kz,
        lhf=lhf,
        qt=qt,
        qc=qc,
        ndz=ndz,
        rg=rg,
        reff=reff,
    )

    # --------------------------------------------------
    # Save results to file
    # --------------------------------------------------
    header = (
        "z[km]    P[bar]      T[K]      qt[g/g]      qc[g/g]      "
        "ndz[#/cm3]      rg[cm]      reff[cm]      cloudf      kz[cm2/s]"
    )

    with open("results.out", "w") as f:
        f.write(header + "\n")
        for i in range(nz):
            f.write(
                f"{z[i]/1e5:8.3f} "     # km
                f"{p[i]*1e-6:10.5f} "   # bar
                f"{t[i]:10.3f} "
                f"{qt_out[i]:12.5e} "
                f"{qc_out[i]:12.5e} "
                f"{ndz_out[i]:12.5e} "
                f"{rg_out[i]:12.5e} "
                f"{reff_out[i]:12.5e} "
                f"{cloudf_out[i]:10.3f} "
                f"{kz_out[i]:12.5e}\n"
            )

    print("✅ Results written to results.out")

if __name__ == "__main__":
    main()
