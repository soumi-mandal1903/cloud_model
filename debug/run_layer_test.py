import sys
sys.path.insert(0, "/home/sansar1/Codes/cloud_model/python_v2")
from read_voyager import read_voyager
from layer import layer
import numpy as np

input_file = "input/profiles/voyager.input"
teff, grav, nz, z, z_top, p, p_top, t, t_top, chf = read_voyager(input_file)

iz = nz // 2
print(f"Testing layer index {iz} / {nz}")

# Parameters
nsub_max = 8
gas_name = "NH3"
mw_atmos = 2.2
kz_min = 1e5
cloudf_min = 0.75
mw_cloud = 17.0
rainf = 1.0
rho_p = 0.84
supsat = 0.0
sig_layer = 2.0
cloudf = cloudf_min
q_below = 1.34e-4
kz = kz_min

# choose top/bot edges for this layer
t_top_edge = t_top[iz]
t_bot_edge = t_top[iz+1]

p_top_edge = p_top[iz]
p_bot_edge = p_top[iz+1]

outputs = layer(nsub_max, gas_name, grav, mw_atmos, kz_min, cloudf_min, mw_cloud, rainf, rho_p, supsat, sig_layer, cloudf, q_below, t[iz], p[iz], kz, chf[iz], t_top_edge, t_bot_edge, p_top_edge, p_bot_edge)

qc_layer, qt_layer, rg_layer, reff_layer, ndz_layer, q_below_next, report_status_r, report_status_q = outputs

print("qc_layer:", qc_layer)
print("qt_layer:", qt_layer)
print("rg_layer:", rg_layer)
print("reff_layer:", reff_layer)
print("ndz_layer:", ndz_layer)
print("status_r, status_q:", report_status_r, report_status_q)

# show nan counts in outputs arrays if available (none here)


