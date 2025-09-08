# main.py
import numpy as np
import matplotlib.pyplot as plt
from read_voyager import read_voyager
from eddysed import eddysed

def main():
    # ---- Read Voyager input profile ----
    fname = "/home/sansar1/Codes/cloud_model/voyager.input"
    p, t, dz, grav, teff, nz, z, chf = read_voyager(fname)

    # ---- Model parameters ----
    mw_atmos = 2.2        # mean molecular weight of H2-He atmosphere [g/mol]
    mw_cloud = 17.0       # NH3 [g/mol]
    rho_p = 0.9           # particle density [g/cm^3]
    rainf = 3.0           # sedimentation efficiency factor
    supsat = 0.0          # assume no supersaturation
    sig_layer = 2.0       # width of size distribution
    mixl_min = 1e4        # minimum mixing length [cm]
    cloudf_min = 0.05     # minimum cloud fraction
    q_init = 1e-6         # initial vapor mixing ratio at base (g/g)

    # ---- Run eddysed ----
    qc, qt, ndz, rg, reff, cloudf = eddysed(
        p=p, t=t, dz=dz,
        grav=grav, mw_atmos=mw_atmos, teff=teff,
        mixl_min=mixl_min, cloudf_min=cloudf_min,
        mw_cloud=mw_cloud, rainf=rainf,
        rho_p=rho_p, supsat=supsat, sig_layer=sig_layer,
        q_init=q_init, nsub_max=8, delta=1e-8
    )

    # ---- Convert units ----
    p_bar = p * 1e-6            # dyne/cm^2 → bar
    qc = np.clip(qc, 1e-30, None)  # avoid log(0)
    log_qc = np.log10(qc)

    # ---- Plot ----
    plt.figure(figsize=(6, 8))
    plt.plot(log_qc, p_bar, marker="o", linestyle="-", markersize=4)
    plt.gca().invert_yaxis()
    plt.yscale("log")
    plt.xlabel("log10(qc) [g/g]")
    plt.ylabel("Pressure [bar]")
    plt.title("Voyager Profile — NH3 Cloud (EddySed)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("cloud_profile.png", dpi=200)
    print("Saved figure to cloud_profile.png")
    plt.show()

if __name__ == "__main__":
    main()
