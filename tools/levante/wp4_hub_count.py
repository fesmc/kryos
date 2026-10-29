#!/usr/bin/env python3
"""Why does the hub log report 8385 cells removed by H_ice_min on NH-32KM, when the raw file has 1313
cells with 0 < H < 1 m (gate-log, phase 2 WP4 design Q8; offline, no kryos code)?

htopo_init (fesm-utils htopo.f90) prepares the reference in this order, in sp:
    remove_lakes   topoprep_remove_lakes(H_ice, z_bed, z_srf, 0, rho_ice, rho_sw)
    H_ice_min      n = count(H_ice < H_ice_min .and. H_ice /= 0)      ! negatives count; then H_ice < H_ice_min -> 0
Hypotheses: (1) remove_lakes runs first and rewrites H_ice = z_srf - z_bed in the cells where the
ice base is above or below the bed (thin, or negative where z_srf < z_bed), so the count is of the
rewritten field; (2) the count includes negative values, which `0 < H < 1` does not.

    wp4_hub_count.py <topo file> [H_ice_min=1.0] [rho_ice=910] [rho_sw=1028]

Prints, for the raw field and after remove_lakes: the cells with 0 < H < min, H < 0, and the
number htopo logs (H < min and H /= 0), and which cells remove_lakes changed. The arithmetic is
float32 in the order of the Fortran. Needs numpy and netCDF4 (Levante: `module load python3`).
"""
import sys
import numpy as np
import netCDF4 as nc


def load(ds, name):
    v = np.squeeze(ds[name][:])
    return np.ma.filled(v, np.nan).astype(np.float32) if np.ma.isMaskedArray(v) else v.astype(np.float32)


def remove_lakes(H, z_bed, z_srf, z_sl, rho_ice, rho_sw):
    """topoprep_remove_lakes (fesm-utils topoprep.f90), float32."""
    f32 = np.float32
    rho_sw_ice = f32(rho_sw) / f32(rho_ice)
    H = H.copy()
    H_grnd = H - rho_sw_ice * np.maximum(f32(z_sl) - z_bed, f32(0.0))
    gap = (z_srf - H) > z_bed
    lake = (H_grnd > 0) & gap
    thick = (~lake) & ((z_srf - H) < z_bed)
    H_new = np.where(lake | thick, z_srf - z_bed, H)
    return H_new.astype(np.float32), lake, thick


def counts(label, H, hmin):
    fin = np.isfinite(H)
    thin = fin & (H > 0) & (H < hmin)
    neg = fin & (H < 0)
    logged = fin & (H < hmin) & (H != 0)
    print(f"{label:34s} 0<H<{hmin:g}: {int(thin.sum()):7d}   H<0: {int(neg.sum()):7d}   htopo count (H<{hmin:g}, H/=0): {int(logged.sum()):7d}"
          f"   non-finite: {int((~fin).sum())}")


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    path = argv[1]
    hmin = np.float32(float(argv[2])) if len(argv) > 2 else np.float32(1.0)
    rho_ice = float(argv[3]) if len(argv) > 3 else 910.0
    rho_sw = float(argv[4]) if len(argv) > 4 else 1028.0
    ds = nc.Dataset(path)
    H, z_bed, z_srf = load(ds, "H_ice"), load(ds, "z_bed"), load(ds, "z_srf")
    print(f"{path}: {H.shape}, rho_ice = {rho_ice}, rho_sw = {rho_sw}")
    counts("raw", H, hmin)
    H2, lake, thick = remove_lakes(H, z_bed, z_srf, 0.0, rho_ice, rho_sw)
    counts("after remove_lakes (hub order)", H2, hmin)
    print(f"remove_lakes rewrote {int((lake | thick).sum())} cells: {int(lake.sum())} englacial lakes filled, "
          f"{int(thick.sum())} erroneously thick reduced")
    changed = (H2 != H)
    thin_new = np.isfinite(H2) & (H2 < hmin) & (H2 != 0)
    print(f"htopo-count cells that remove_lakes created (not in the raw count): "
          f"{int((thin_new & ~(np.isfinite(H) & (H < hmin) & (H != 0))).sum())}")
    print(f"  of them with H < 0: {int((thin_new & (H2 < 0) & changed).sum())}, with 0 < H < {float(hmin):g}: "
          f"{int((thin_new & (H2 > 0) & changed).sum())}")
    # the other order, for the record (H_ice_min first)
    H3 = np.where(np.isfinite(H) & (H < hmin), np.float32(0.0), H)
    H3, _, _ = remove_lakes(H3, z_bed, z_srf, 0.0, rho_ice, rho_sw)
    print("(H_ice_min before remove_lakes would leave "
          f"{int((np.isfinite(H3) & (H3 < hmin) & (H3 != 0)).sum())} cells below {float(hmin):g} m after both steps)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
