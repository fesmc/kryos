# setup-diag: basal and flow setup of the science case (NH-32KM, ANT-32KM)

Analysis only, no new runs. Branch `setup-diag`. Question: what controls the thickness biases of the calib-cf states, and is `cf_ref` an effective knob there?

- Fields: Yelmo restarts (final state) of `$G/runs/calib/kx-624f686/S9/` North-pd-cf10 (and North-pd, cf 0.2), North-lgmrecon-cf10 (and North-lgmrecon), Antarctica-pd-cf05 (and -cf10); the standalone Greenland gate run with the yelmox_Greenland settings (`$G/runs/grl/kx-af8f1b9/S3/base`, GRL-16KM, 15 kyr).
- Obs (FesmData): Joughin 2018 velocity, MAR v3.14 1981-2010 SMB (ice-sheet mask; GRL-32KM mapped to NH-32KM by nearest neighbour), Rignot 2011, RACMO2.3p2 1981-2010, Dawson 2022 basal state, GHF datasets.
- Script `tools/levante/setup_diag.jl` (`sbatch tools/levante/setup_diag.sh`); full tables in `setup-diagnosis-tables.md`; maps `$G/setup-diag/out/setup_diag_*.png`.
- Regularized Coulomb, q = 0.2, u0 = 100 m/yr: at fixed stress u_b ~ u0 (taud/c_bed)^5, so sliding needs c_bed of order taud. c_bed = cb_ref N_eff.

## Diagnosis

### Greenland (North PD, cf_ref 0.1)

| | Greenland | interior > 300 km | 100-300 km | margin < 100 km |
|---|---|---|---|---|
| H bias, on reference ice [m] | +67 | **-257** | +56 | **+266** |
| SMB model / MAR, MAR ice-sheet cells [m ice/yr] | 0.36 / 0.20 | 0.23 / 0.20 | 0.41 / 0.35 | **0.39 / -0.17** |
| ablation area model / MAR | 0.04 / 0.14 | 0 / 0 | 0.01 / 0.04 | **0.17 / 0.59** |
| T_srf model / MAR [degC] | -17.1 / -22.5 | -24.4 / -29.2 | -19.2 / -22.5 | -10.2 / -14.8 |
| temperate bed / wet till (area) | 0.65 / 0.15 | 0.59 / 0.00 | 0.62 / 0.06 | 0.70 / 0.40 |
| basal melt, temperate, median [mm/yr] (till_rate 8) | 2.6 | 0.8 | 2.2 | 14.8 |
| N_eff/P0 > 0.9 (area) | **0.85** | 1.00 | 0.94 | 0.60 |
| c_bed/taud, median | **21** | 37 | 24 | 5.6 |
| cf for c_bed = taud (taud/N_eff), median | **0.004** | 0.002 | 0.004 | 0.016 |
| u_b/u_s > 0.5 (area) | 0.09 | 0 | 0.01 | 0.23 |
| log10 u_s model/obs, median | +0.01 | +0.05 | -0.06 | +0.15 |
| cf 0.2 -> 0.1: dH [m] | -55 | -78 | -47 | -63 |

1. **cf_ref is not an effective knob for Greenland.** 65 % of the bed is temperate but only 15 % has till water: basal melt (median 2.6 mm/yr) is below the North `bkt_till_rate` of 8 mm/yr on 75 % of the temperate bed, so the bucket stays dry and the TILL closure gives N = P0 on 85 % of the ice. c_bed is 20 times the driving stress, the flow is internal deformation (sliding-dominated on 9 %, at the margins only). For friction to engage on these dry beds cf_ref would have to be about 0.004: the sweep to 0.015 would be tuning the dry-till overburden, not till strength. This is the "tiny cf_ref".
2. **The V_af excess is at the margins and comes from the SMB, not the flow.** Interior -257 m, margins +266 m (east coast up to +800 m). Surface speeds match Joughin 2018 (median ratio 1.0, interior 1.1). The modelled SMB on MAR's ice sheet is 587 against 326 Gt/yr; in the 100 km margin band +0.39 against -0.17 m ice/yr, ablation area 17 % against 59 %. The ablation zone is mostly missing (chion PDD, ERA5 reference, 32 km; part of it is the thick-margin elevation feedback, a free-geometry run cannot separate them).
3. **The interior is too thin** with the right speed and SMB (interior SMB 0.23 against 0.20): too soft or too warm ice. T_srf is 3-5 K above MAR everywhere, 59 % of the interior bed is temperate. Lowering cf_ref thins the interior further (-78 m per halving): the calib-cf direction makes the interior worse while it trims the margins.
4. **The standalone yelmox Greenland has the same pattern** (GRL-16KM, cf_ref 0.8, till_rate 1 mm/yr, GHF 50, MAR+ITM): bias +194 m, interior -198, margins +444, SMB 601 against 340 Gt/yr, margin ablation area 3 % against 38 %, c_bed/taud 128, no sliding-dominated area. So cf_ref 0.8 there is not a calibrated value either; it sits in the same non-sliding regime. Neither setup is a friction-calibrated Greenland.

### North LGM (lgmrecon, cf_ref 0.1, frozen 21 ka)

| | North America | W Eurasia | Greenland |
|---|---|---|---|
| temperate / wet till (area) | 0.87 / 0.23 | 0.62 / 0.27 | 0.41 / 0.11 |
| N_eff/P0 > 0.9 (area) | 0.78 | 0.74 | 0.89 |
| cb_ref median | **0.007** | **0.007** | 0.067 |
| c_bed/taud median | 3.0 | 1.9 | 15.6 |
| u_b/u_s > 0.5 (area) | 0.32 | 0.44 | 0.16 |
| cf 0.2 -> 0.1: dH [m] | -98 | -129 | -70 |

5. The Laurentide bed is temperate (87 %) but the till is dry (N = P0 on 78 %); the sliding that exists comes from cb_ref collapsing to near cf_min (0.001) because `scale_zb` uses the current, isostatically depressed bed: under 2-3 km of ice the bed sits hundreds of m below z0 = -100 m. At the LGM the friction is set by the load and z0/z1, not by cf_ref (halving it: -100 m, a fifth of the +600-700 m overgrowth). The overgrowth is mostly SMB (round 3), and the hydrology has no wet-bed mode to soften it.

### Antarctica PD (cf_ref 0.05, scale_zb exp -400/0, MARINE closure)

| | WAIS | WAIS int. | EAIS | EAIS int. | EAIS bed >= 0 | EAIS bed < 0 |
|---|---|---|---|---|---|---|
| H bias [m] | -154 | -279 | +244 | +276 | **+348** | +142 |
| cb_ref median | 0.017 | 0.015 | 0.049 | 0.050 | 0.050 | 0.023 |
| N_eff/P0 median | 0.55 | 0.66 | 0.99 | 0.99 | **0.99** | 0.83 |
| c_bed/taud median | 2.8 | 3.2 | 12.8 | 22.7 | **21.7** | 5.5 |
| u_b/u_s > 0.5 (area) | **0.62** | 0.46 | 0.09 | 0.01 | 0.01 | 0.18 |
| log10 u_s model/obs, median | **+0.13** | +0.17 | +0.01 | -0.05 | -0.04 | +0.06 |
| temperate (model) / thawed (Dawson 2022) | 0.73 / 0.16 | 0.71 / 0.35 | 0.34 / 0.10 | 0.26 / 0.11 | 0.16 / 0.04 | 0.55 / 0.17 |
| beta median model / Dawson [Pa yr/m] | 2.4e3 / 1e4 | 1.6e4 / 1.6e4 | 7.7e7 / 1.3e4 | | | |
| Q_geo (Shapiro 2004) / range of 6 others [mW/m2] | 92 / 65-96 | | 53 / 38-64 | | | |
| enh_bar | 1.0 | | 1.0 | | | |

6. **EAIS too thick is a stiff-ice problem.** On land beds the MARINE closure gives N = P0 whatever the thermal or water state, c_bed is 22 times taud: no sliding at any plausible cf_ref. Speeds match Rignot 2011 (median ratio 0.9-1.0) and SMB matches RACMO (0.097 against 0.099), yet the ice is 348 m (16 %) too thick: the rate factor is too low. `[ymat_ant]` has `enh_shear = 1` (Greenland and North: 3). With the SIA scaling H ~ E^(-1/8), E = 3 thins by about 13 %, the size of the bias.
7. **WAIS too thin is a too-slippery bed.** The exp elevation scaling (round 1) drops cb_ref to 0.015 on the deep WAIS bed, the WAIS gets Shapiro's 92 mW/m2 (upper end of the seven datasets: Hazzard-Richards 96, Martos 87, Stal 81, Loesing 70), 73 % temperate against 16 % thawed in Dawson, sliding-dominated on 62 %, twice the observed speed, beta four times below Dawson's (linear-law inversion, indicative only).
8. The EAIS-thick / WAIS-thin dipole is therefore structural: scale_zb + MARINE closure + enh 1, with one cf_ref shifting both. Halving cf_ref thins both (WAIS -81 m, EAIS -91 m). The yelmox reference was never run like this: `yelmox_Antarctica.nml` has `sim.opt = True` (L21 cb_ref inversion) and `relax = True`, switched off in the conversion (`-p sim.relax=False sim.opt=False`); its ytill block (cf_ref 0.1, scale_zb 0, z0/z1 -400/0 inactive) and enh_shear 1 are what the inversion compensated.

## Setup audit

Kryos science par (`par/kryos_NH-ANT.nml`, `_calib`) against yelmox (`yelmox_North`, `_Greenland`, `_Antarctica`). All North Yelmo settings are converted 1:1 from `yelmox_North.nml` (yelmox 5ac89c8, conversion report `par/kryos_North_conversion.md`); its ytill/yhyd values come from the legacy yelmox North namelists (palma-ice yelmox v1.15.1, the Laurentide setup; ported in yelmox 6c9ebcc) and were never tuned for Greenland. Physics defaults (`input/yelmo_defaults.nml`: frz_scale on, frz_min 1e-3, capacity basal BC, gl_temperate) are the same in yelmox and kryos.

| setting | kryos North | yelmox Greenland | kryos ANT (calib) | yelmox ANT | effect / flag |
|---|---|---|---|---|---|
| `yhyd.bkt_N_closure` | 3 TILL | 3 TILL | 2 MARINE, p 1 | 2 MARINE | **ANT: N = P0 on every bed above sea level (finding 6)** |
| `yhyd.bkt_till_rate` [m/yr] | **0.008** | 0.001 | 0.001 | 0.001 | **North: temperate beds dry (findings 1, 5)** |
| `yhyd.W_til_max` [m] | 0.4 | 2.0 | 2.0 | 2.0 | North: smaller bucket, faster to fill and drain |
| `ytill.cf_ref` | 0.1 | 0.8 | 0.05 (0.1 base) | 0.1 + L21 opt | not a calibrated value in any of them |
| `ytill.scale_zb`, z0/z1 [m] | 2, -100/100 | 0 | **2 (round 1)**, -400/0 | 0 (z0/z1 unused) | **ANT: z0/z1 never tuned, creates the dipole (7, 8)**; both: uses the depressed bed (5) |
| `ytill.cf_min` | 0.001 | 0.01 | 0.01 | 0.01 | North: LGM floor 10x lower |
| `sim.opt` / `sim.relax` | off / off | off / off | **off / off** | **on / on** | **ANT reference relied on the inversion (8)** |
| `ymat.enh_shear` / `enh_shlf` | 3 / 0.7 | 3 / 0.7 | **1** / 0.5 | 1 / 0.5 | **EAIS too stiff without inversion (6)** |
| `ydyn` beta law, q, u0, beta_min | RC, 0.2, 100, 100 | same | same | same | - |
| `ydyn.beta_gl_scale` | 0 | 0 | 2 | 2 | - |
| `ydyn.ssa_iter_conv` / LIS tol | **0.1 / 0.1** | 0.01 / 0.01 | 0.01 / 0.01 | 0.01 / 0.01 | North SSA loose (minor) |
| GHF | Shapiro 2004 | const 50 | Shapiro 2004 | Shapiro 2004 | Greenland: Shapiro 58 (E/SE 70-90), Colgan 44; WAIS: Shapiro 92, upper end of the range 65-96. Units checked: W m-2 -> mW m-2 in the adapter, values plausible |
| `ytherm` | enth, capacity BC | same (defaults) | same | same | - |
| `nz_aa` | 10 | 6 | 10 | 10 | - |
| SMB | chion PDD, ERA5 1981-2010 | chion ITM (alb_ice 0.31 MAR-calibrated), MAR 1961-1990 | chion PDD, ERA5 | chion PDD | **Greenland margins: no ablation zone (2)** |
| calving | vm-m16 | equil | threshold (round 5) | equil | ANT front: calib-cf round 4 |
| `sed` / `scale_sed` | CRUST1 / 0 | none / 0 | none / 0 | none / 0 | sediment field unused |

calib-cf changes (`tools/bin/kryos-calib-par.jl`, `docs/dev/calib-cf.md`):

- **ANT `scale_zb` 0 -> 2 (round 1)**: questionable. It switched on z0/z1 = -400/0 that were inactive in yelmox (no tuning behind them); it makes marine beds 10-100 times weaker than land beds on top of the MARINE closure, which already weakens marine beds near flotation. The dipole followed.
- cf_ref per domain (ANT 0.05, North 0.1): judged on V_af, which mixes opposite-signed regional errors (Greenland interior/margin, EAIS/WAIS). The PD V_af gains came with worse interiors.
- Threshold calving for ANT: not related to the basal question.

## Correction plan (ranked by likely impact)

1. **Stop the cf_ref sweeps for Greenland and ANT** (the 0.015 sweep included): in both, the thickness bias is not in the friction-limited regime, and lowering cf_ref worsens the interiors.
2. **North SMB at the Greenland margin** (finding 2; Greenland V_af). Check chion on the North grid at fixed PD geometry against MAR (the method of the GRL ITM alb_ice calibration), to separate the scheme from the elevation feedback; then decide between the PDD parameters (ddf, sigma), the ITM with the calibrated `itm_chion` as in yelmox Greenland, or the reference climate. Target: MAR integral and ablation area per band. No friction tuning before this.
3. **ANT friction structure** (findings 6-8). Options, to choose before any run:
   a. restore the L21 cb_ref inversion for ANT PD (what yelmox ANT was built on; the calib-grl branch has `opt_cb_ref` work), with scale_zb 0 and enh_shear 1 as in yelmox;
   b. no inversion: `enh_shear` 1 -> 3, `scale_zb` 2 -> 0 (back to the base), then one cf_ref; optionally TILL closure for ANT so that land beds can slide when wet;
   c. GHF: replace Shapiro 2004 (92 mW/m2 over WAIS, upper end) by a newer dataset (Stal 2020 or Martos 2017), or an ensemble mean.
4. **North basal hydrology** (findings 1, 5): `bkt_till_rate` 8 -> 1 mm/yr and `W_til_max` 0.4 -> 2 m (Greenland/ANT values), so temperate beds with a few mm/yr of melt get wet and N drops to delta P0. This turns cf_ref into a real knob in NA and EUR (LGM overgrowth) and the Greenland outlet zones; the Greenland interior (melt 0.8 mm/yr) stays dry. Expect a large change of the LGM North; re-tune cf_ref afterwards (likely upward). Needs Alex: the 8 mm/yr is the Laurentide tuning.
5. **Greenland interior** (finding 3): check what T_srf chion hands to Yelmo (3-5 K above MAR; `firn_fac` refreezing warming?) and the interior GHF (Shapiro 49 vs Colgan 43 mW/m2) before touching enhancement.
6. **scale_zb on the depressed bed** (finding 5): `calc_lambda_bed_exp` takes the current z_bed; under LGM loads cb_ref goes to cf_min everywhere. Use a fixed reference bed (rebounded or PD) for the scaling (Yelmo change), or drop scale_zb in the North and use regions instead (the till-regions work).
7. Minor: North `ssa_iter_conv` / LIS tol 0.1 -> 0.01 as Greenland and ANT.

After 2-4, re-tune in this order: SMB (North, fixed geometry) -> basal hydrology (North) / friction structure (ANT) -> one cf_ref per domain on the PD state, judged by regional bands (interior / 100-300 km / margin, EAIS land / marine, WAIS) and velocity, not only V_af -> LGM.

## Open for Alex

1. ANT: inversion (3a) or structural fix (3b)? 3a reproduces the yelmox reference; 3b keeps the ANT free for the glacial cycle.
2. North till hydrology: accept 1 mm/yr / 2 m (one setting for Laurentide and Greenland), or region-dependent?
3. Greenland SMB: PDD tuning, or ITM with the MAR-calibrated parameters as in the standalone?
4. scale_zb on a reference bed: Yelmo change, or drop scale_zb in favour of till regions?
