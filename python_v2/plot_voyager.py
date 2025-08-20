# plot_voyager.py
import numpy as np
import matplotlib.pyplot as plt

INPUT = "/home/sansar1/Codes/cloud_model/voyager.output"   # change if filename differs

def read_voyager(fname):
    pressures = []
    qc = []
    with open(fname, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            p = float(parts[0])   # pressure in bar
            q = float(parts[1])
            pressures.append(p)
            qc.append(q)
    return np.array(pressures), np.array(qc)

def main():
    p_bar, qc = read_voyager(INPUT)

    # Sort by pressure (just in case). Many profiles are top->bottom or bottom->top.
    order = np.argsort(p_bar)   # ascending pressure
    p_bar = p_bar[order]
    qc = qc[order]

    # Prepare masked log10(qc) (mask zeros and negative values)
    qc_positive = qc > 0.0
    logqc = np.full_like(qc, np.nan)
    logqc[qc_positive] = np.log10(qc[qc_positive])

    fig, axes = plt.subplots(1, 2, figsize=(12, 7), sharey=True)

    # Left: linear qc
    ax = axes[0]
    ax.plot(qc, p_bar, marker='o', linestyle='-', markersize=4)
    ax.invert_yaxis()
    ax.set_xlabel('qc (g/g)')
    ax.set_ylabel('Pressure (bar)')
    ax.set_title('qc vs Pressure (linear)')
    ax.grid(True)

    # Right: log10(qc) where qc>0
    ax = axes[1]
    if np.any(qc_positive):
        ax.plot(logqc[qc_positive], p_bar[qc_positive], marker='o', linestyle='-',
                markersize=4)
        ax.set_xlabel('log10(qc)')
        ax.set_title('log10(qc) vs Pressure (qc>0)')
    else:
        ax.text(0.5, 0.5, 'No positive qc values to plot (all zeros)', ha='center')
        ax.set_xlabel('log10(qc)')
        ax.set_title('log10(qc) vs Pressure')
    ax.invert_yaxis()
    ax.grid(True)

    plt.tight_layout()
    plt.savefig("voyager_qc_profile.png", dpi=200)
    print("Saved figure to voyager_qc_profile.png")
    plt.show()

if __name__ == "__main__":
    main()
