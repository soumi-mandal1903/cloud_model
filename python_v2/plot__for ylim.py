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

    # Sort by pressure (just in case).
    order = np.argsort(p_bar)
    p_bar = p_bar[order]
    qc = qc[order]

    # Apply mask: only keep 0.1 <= p <= 1 bar
    mask = (p_bar >= 0.1) & (p_bar <= 1.0)
    p_bar = p_bar[mask]
    qc = qc[mask]

    # Prepare masked log10(qc) (mask zeros and negatives)
    qc_positive = qc > 0.0
    logqc = np.full_like(qc, np.nan)
    logqc[qc_positive] = np.log10(qc[qc_positive])

    fig, axes = plt.subplots(1, 2, figsize=(12, 7), sharey=True)

    # Left: linear qc
    ax = axes[0]
    ax.plot(qc, p_bar, marker='o', linestyle='-', markersize=4)
    ax.invert_yaxis()
    ax.set_yscale('log')
    ax.set_ylim(1.0, 0.1)  # restrict y-axis between 0.1 and 1 bar
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
        ax.set_ylim(1.0, 0.1)
        ax.set_yscale('log')
    else:
        ax.text(0.5, 0.5, 'No positive qc values in range (0.1–1 bar)', ha='center')
        ax.set_xlabel('log10(qc)')
        ax.set_title('log10(qc) vs Pressure')
    
    ax.grid(True)

    plt.tight_layout()
    plt.savefig("voyager_qc_profile_masked.png", dpi=200)
    print("Saved figure to voyager_qc_profile_masked.png")
    plt.show()

if __name__ == "__main__":
    main()
