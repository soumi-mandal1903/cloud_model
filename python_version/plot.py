import numpy as np
import matplotlib.pyplot as plt

filename = '/home/sansar1/Codes/cloud_code_nkl/cloud_code_nkl/input/profiles/voyager.input'

with open(filename, 'r') as f:
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
    P.append(float(parts[1].replace('E', 'e')))

T = np.array(T)
P = np.array(P)

# Apply filter: keep only values within the desired ranges
mask = (T >= 110) & (T <= 160) & (P >= 0.1) & (P <= 1.0)
T_filtered = T[mask]
P_filtered = P[mask]

plt.figure(figsize=(6, 8))
ax = plt.gca()
ax.plot(T_filtered, P_filtered, marker='o')

# Invert y-axis
ax.invert_yaxis()

# Axis labels and formatting
ax.set_xlim(110, 160)
ax.set_ylim(1.0, 0.1)
ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (bar)')
plt.yscale('log')
plt.title('Temperature vs Pressure (Voyager Profile)')
plt.grid(True)
plt.tight_layout()
plt.show()
