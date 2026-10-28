#!/usr/bin/env python3
"""R2' vs R0 differences of one WP0 gate, as documented for grl/ant (gate-log 'Stack S2 and baseline R2').
    wp0_r2diff.py <gate>   (module load python3)
Prints: the hub preparation log of R2', H_ice / H_ice_ref at t0 (R0 and R2' against their own reference and
against each other: cells differing, > 10 m, max), the V_ice difference over the run, the z_bed difference."""
import re, sys, numpy as np, netCDF4 as nc
G = "/work/ba1442/robinson/kryos-gate/ref"; g = sys.argv[1]
r0 = f"{G}/{g}/yx-918d996/v24/cont"; r2 = f"{G}/{g}/yx-5ac89c8/S4/cont"
for l in sorted({l.strip() for l in open(f"{r2}/out.out") if re.search(r"htopo_init|prepared|cleaning", l)}):
    print("R2 log:", l)
def f(d, v, k=0):
    return nc.Dataset(f"{d}/yelmo.nc")[v][k].filled(np.nan)
for lab, d in (("R0", r0), ("R2'", r2)):
    H, Hr = f(d, "H_ice"), f(d, "H_ice_ref"); x = np.abs(H - Hr)
    print(f"{lab} t0 H_ice vs H_ice_ref: cells {int((x>0).sum())}, max {np.nanmax(x):.4g}, cells >10 m {int((x>10).sum())}")
for v in ("H_ice", "H_ice_ref", "z_bed"):
    try:
        x = np.abs(f(r0, v) - f(r2, v)); print(f"R2'-R0 t0 {v}: cells {int((x>0).sum())}, max {np.nanmax(x):.4g}, >10 m {int((x>10).sum())}")
    except Exception as e: print(v, e)
a, b = nc.Dataset(f"{r0}/yelmo_ts.nc"), nc.Dataset(f"{r2}/yelmo_ts.nc")
n = min(a["V_ice"].shape[0], b["V_ice"].shape[0])
d = np.abs(a["V_ice"][:n] - b["V_ice"][:n]); print(f"V_ice diff over {n} records: max {d.max():.4g} (units {getattr(a['V_ice'],'units','?')}), final {d[-1]:.4g}")
