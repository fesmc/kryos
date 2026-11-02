# ant-friction: ANT friction structure without inversion (first step)

Branch `ant-friction` (from main b51cb81). Decision (Alex, setup-diagnosis open question 1): no cb_ref inversion; `[ymat_ant] enh_shear` 1 -> 3, `scale_zb` stays on (exp), explore `[ytill_ant] z0/z1` (-400/0 since calib-cf round 1, never tuned). A till-regions friction scheme follows separately.

## Setup

- Stack S9, ANT PD case of calib-cf (30 kyr cold start from the hub reference, PD climate, 64 threads), run root `$G/runs/calib/kx-624f686/S9/Antarctica-pd-af<cf>-z<-z0>_<-z1>` (`af05` = cf_ref 0.05, `af1` = 0.1; the folders carry the ".0" of the arguments, e.g. `Antarctica-pd-af05-z1600.0_0.0`). 8 runs, 55-58 min each on an exclusive node (7.6 node-hours). Maps `$G/calib/maps/ant-pd-af*-maps.png`.
- `tools/levante/ant_friction.sh <cf_ref> <z0> <z1> ...` -> `calib.sh` with overrides on `$G/calib/par/Antarctica-pd.nml`: `ymat_ant.enh_shear = 3`, `ytill_ant.cf_ref/z0/z1`; calving as the round-5 setup (calib-cf round 4 `-cf05-vm`: vm-m16, `tau_ice_flt` 85 kPa, `calv_grnd_method = "zero"`, `H_min_*` 0). Everything else as calib-cf (MARINE closure, cf_min 0.01, Shapiro GHF).
- Exp scaling (yelmo `calc_lambda_bed_exp`): cb_ref = max(cf_min, cf_ref exp((z_bed - z1) / (z1 - z0))), capped at cf_ref above z1. z1 is where the weakening starts, z1 - z0 the e-folding depth. -400/0: cb_ref = cf_ref / e at -400 m, cf_min (0.01) reached at -640 m (cf_ref 0.05) or -920 m (0.1).
- Grid: cf_ref 0.05, 0.1 x z0/z1 -400/0 (current), -800/0, -1600/0, -1000/-200 (8 runs).
- Metrics: `tools/levante/ant_friction_bands.sh` (regional bands, masks from BedMachine v4: grounded reference ice, land / marine bed, distance to its edge; WAIS / EAIS from region_3; velocity against Rignot 2011); `calib_metrics.sh` (V_af, area, grounding line).

## Results (PD, t = 30 kyr)

Bands (`$G/calib/out/af-bands.md`): H bias mean [m] on the grounded BedMachine ice; u ratio = median log10(u_s model / Rignot 2011); sliding = area with u_b/u_s > 0.5. V_af, GL from `calib_metrics` (`$G/calib/out/ant-pd-af*.md`; BedMachine V_af 56.83 m). All runs converged (|rate| <= 0.02 m/kyr), front mismatch not re-evaluated (calving as round 4 vm).

| run | cf_ref | z0/z1 | H bias all (rms) | WAIS | EAIS land | EAIS marine | int > 300 km | 100-300 km | margin < 100 km | u ratio WAIS / EAIS / margin | sliding WAIS / EAIS | dV_af [m] | GL median / p90 [km] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| cf05 (enh 1, equil) | 0.05 | -400/0 | +192 (430) | -160 | +362 | +118 | +225 | +144 | +167 | +0.16 / +0.05 / +0.12 | 0.62 / 0.09 | +6.05 | 0 / 32 |
| cf05-vm (enh 1) | 0.05 | -400/0 | +133 (426) | -213 | +314 | +34 | +174 | +84 | +83 | +0.16 / +0.07 / +0.13 | 0.61 / 0.09 | +4.57 | 0 / 45 |
| af05-z400_0 | 0.05 | -400/0 | -43 (353) | -264 | +67 | -109 | -74 | -38 | +39 | +0.20 / +0.15 / +0.18 | 0.56 / 0.07 | -0.35 | 0 / 32 |
| af05-z800_0 | 0.05 | -800/0 | -11 (335) | -172 | +70 | -59 | -33 | -4 | +46 | +0.15 / +0.15 / +0.18 | 0.41 / 0.05 | +0.46 | 0 / 32 |
| **af05-z1600_0** | 0.05 | -1600/0 | +22 (334) | -70 | +83 | -35 | -1 | +40 | +63 | +0.11 / +0.14 / +0.17 | 0.28 / 0.05 | +1.31 | 0 / 32 |
| af05-z1000_200 | 0.05 | -1000/-200 | +26 (334) | -103 | +97 | -25 | +3 | +37 | +77 | +0.12 / +0.14 / +0.17 | 0.33 / 0.05 | +1.45 | 0 / 32 |
| af1-z400_0 | 0.1 | -400/0 | +14 (345) | -176 | +92 | -20 | -22 | +23 | +107 | +0.14 / +0.14 / +0.14 | 0.41 / 0.04 | +1.20 | 0 / 32 |
| af1-z800_0 | 0.1 | -800/0 | +68 (344) | -40 | +121 | +29 | +35 | +91 | +128 | +0.08 / +0.12 / +0.14 | 0.20 / 0.03 | +2.62 | 0 / 45 |
| af1-z1600_0 | 0.1 | -1600/0 | +93 (348) | -1 | +137 | +62 | +61 | +126 | +140 | +0.09 / +0.11 / +0.14 | 0.13 / 0.03 | +3.31 | 0 / 45 |
| af1-z1000_200 | 0.1 | -1000/-200 | +71 (342) | -23 | +119 | +33 | +35 | +99 | +137 | +0.08 / +0.12 / +0.13 | 0.16 / 0.03 | +2.71 | 0 / 45 |

- **enh_shear 3** removes most of the EAIS land excess (+314-362 -> +67-137 m) and the rms drops from 426-430 to 334-353 m: finding 6 of the diagnosis holds. Cost: the EAIS is now faster than Rignot 2011 (u ratio +0.05 -> +0.11-0.15, a factor 1.3-1.4), the margins too (+0.13-0.18).
- **z0/z1 sets the WAIS**: at cf_ref 0.05 the WAIS bias goes -264 / -172 / -70 m for an e-folding depth of 400 / 800 / 1600 m, the sliding-dominated WAIS area 0.56 / 0.41 / 0.28 and the WAIS speed ratio +0.20 / +0.15 / +0.11. The EAIS land changes little (+67-83). Shifting z1 to -200 m (-1000/-200) is close to -1600/0, slightly worse on the WAIS.
- **cf_ref** shifts every band together (0.05 -> 0.1: +60-110 m on EAIS land / marine and the margin), as before. 0.1 balances the WAIS (-1 to -40 m) only with a thicker EAIS and margin (+120-140 m) and V_af +2.6-3.3 m.
- Recommended: **cf_ref 0.05, z0/z1 = -1600/0** (`af05-z1600_0`): lowest rms (334, as -1000/-200), all bands within -70..+83 m, interior -1 m, V_af +1.3 m, WAIS sliding area halved (0.56 -> 0.28) and WAIS speed closest to Rignot of the 0.05 runs.

## Per-basin cf_ref (till_cf)

Branch `till-cf` with Yelmo `till-cf-field`. Yelmo's `cf_ref` is a field (`dyn%now%cf_ref`, written as a 2D variable): `ytill.cf_ref` everywhere after `yelmo_init` and `yelmo_restart_init`, not read from the restart. The key `till_cf` of the ice instance names a group of the par file that gives ids their own value (Yelmo's `set_cf_ref_ids`, exact id match), applied after both calls, so a changed group takes effect on a restart:

```
&ice_ant
    till_cf = "till_cf_ant"   ! "" (default): ytill.cf_ref everywhere
/
&till_cf_ant
    ids    = "basins"         ! "basins" (Yelmo's bnd%basins) | "regions" (bnd%regions)
    codes  = 1 2              ! ids with their own cf_ref
    cf_ref = 0.05 0.1         ! [-- or deg]
/
```

`ids = "basins"` is Yelmo's basin set, `[yelmo_masks_ant] basins` (`[Antarctica] basins` does not change it; smoke runs 2026-11-02, `/work/ba1442/robinson/till-cf/runs`). init stops on an unknown key or `ids`, a different number of `codes` and `cf_ref`, or a code not in the id field. `cf_min` and the `scale_zb` / `scale_sed` scaling apply on top (the law runs from `cf_min` to the cell's `cf_ref`). No science par file has a `till_cf` group yet.

## Open

1. The EAIS is 30-40 % faster than Rignot 2011 with enh_shear 3 while the land beds are still +80 m thick: enh_shear 2 may balance speed and thickness better (not run).
2. WAIS still -70 m and the trend with the e-folding depth has not flattened: a weaker scaling (-3200/0, or scale_zb 0 with cf_ref 0.05) may close it; a per-basin cf_ref (`till_cf`, above) is the structural answer.
3. The calving front (vm-m16 at 85 kPa, untuned) is unchanged from round 4; the tau_ice_flt sweep (calib-cf round 5) comes on this friction setup.
4. Not in the generator yet: `tools/bin/kryos-calib-par.jl` still has enh_shear 1 (base) and z0/z1 -400/0 for ANT.
