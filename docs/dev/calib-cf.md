# calib-cf: basal friction (cf_ref) for the science par, NH-32KM and ANT-32KM

`$G = /work/ba1442/robinson/kryos-gate`. Branch `calib-cf`. Stack S9 (kryos 624f686, `kryos-v0`: vilma off, FastIsostasy), 64 threads, exclusive node.
Runs `$G/runs/calib/kx-624f686/S9/<Domain>-<case><TAG>`, metrics `$G/calib/out/<ant|nor>-<case><TAG>.md`, maps `$G/calib/maps/`.

## Method

No optimisation: Yelmo's elevation scaling of cf_ref, one cf_ref per domain, judged on the present-day (PD) and LGM states.

- `[ytill_ant]`: `scale_zb = 2` (exp), z0/z1 = -400/0 m, cf_min 0.01. `[ytill_nor]`: exp, -100/100 m, cf_min 0.001.
- Science variants `par/kryos_NH-ANT_calib.nml`, `par/kryos_NH-ANT_vilma2_calib.nml` (`tools/bin/kryos-calib-par.jl sci`): rounds 1-2 cf_ref 0.2 both domains; round 3 ANT 0.05, North 0.2; round 4 ANT 0.05, North 0.1 (`CF_ANT`, `CF_NOR` in the generator; North 0.1 is the value of the base); from round 5 also `[ycalv_ant]` threshold calving (see the decision below). `par/kryos_NH-ANT_vilma2.nml` unchanged.
- `par/calib/*.nml` are regenerated from the science variant, so they carry ANT 0.05 from round 3 and North 0.1 from round 4; earlier runs are reproduced with `ytill_ant.cf_ref=` / `ytill_nor.cf_ref=` overrides.
- Single domain, FastIsostasy, frozen forcing (`tools/bin/kryos-calib-par.jl run`, pars `par/calib/`):
  - pd: t 0 -> 30 kyr, alpha 0, sea level 0, orbit 0 ka, cold start from the hub reference.
  - lgm: t 30 -> 60 kyr, restarted from pd, alpha 1, sea level -120 m, orbit -21 ka.
  - lgmrecon (North): lgm forcing, cold start from PaleoMIST (`recon_nor`, whole domain, no bed filter, smoothing or float cut, as the science par at -80 ka; `recon_time = -21e3` selects the 20 ka record), t 0 -> 20 kyr.
  - lgm-ext (North): lgm restarted from the 60 kyr bundle of lgm, t 60 -> 90 kyr.
- Submit `tools/levante/calib.sh <Domain> <case> [set_nml ...]` (TAG, FROM, RSTT); cf_ref variants by `ytill_ant.cf_ref=...` / `ytill_nor.cf_ref=...`. Metrics `tools/levante/calib_metrics.sh` -> `calib_metrics.jl`; maps `tools/levante/calib_maps.sh` -> `tools/calib_maps.jl` (CairoMakie, plot environment `$G/julia-plot`). Everything runs on Levante (the generator too: `tools/levante/calib_gen.sh`); rounds 1-4 maps re-rendered there (the Mac copies in `$G/calib/maps/mac/`).
- V_af: grounded ice above flotation as m SLE (rho_i 917, rho_sw 1028, ocean 3.625e14 m2). References: PD BedMachine v4 (ANT), GEBCO 2026 (North); LGM PaleoMIST 20/22.5 ka, GLAC-1D 21 ka, ICE-6G_C 21 ka. North regions are lon/lat boxes (GRL, NA, EUR; see the metrics file).
- Rate: Yelmo V_sle over the last 5 kyr. Wall about 52 min per 30 kyr (ANT), 60-80 min (North).

## Antarctica (rounds 1 and 2)

| cf_ref | PD V_af [m] (BedMachine 56.8) | PD dV_af | PD A [1e6 km2] (12.30) | PD H bias / rms, overlap [m] | PD rate [m/kyr] | LGM V_af [m] (refs 69.3-71.8) | LGM A (refs 17.0-20.5) | LGM - PD V_af [m] (refs 12.5-15.0) | LGM rate |
|---|---|---|---|---|---|---|---|---|---|
| 0.2 | 67.99 | +11.2 | 12.18 | +406 / 515 | 0.06 | 74.69 | 13.78 | 6.7 | 0.01 |
| 0.1 | 65.73 | +8.9 | 12.26 | +316 / 466 | 0.05 | 71.23 | 13.79 | 5.5 | -0.01 |
| 0.05 | 62.88 | +6.0 | 12.32 | +207 / 420 | 0.03 | 67.92 | 13.84 | 5.0 | 0.00 |

PD grounding line: median distance 0 km, p90 45 / 32 / 32 km (0.2 / 0.1 / 0.05).
LGM against PaleoMIST 20 ka: dV_af +4.5 / +1.1 / -2.2 m, overlap bias +321 / +198 / +82 m, missing grounded area 3.3 / 3.3 / 3.2 M km2 (model-only area 0).

- PD: too thick in the EAIS interior and along the Transantarctic Mountains, too thin in the WAIS interior (`ant-pd*-maps.png`). Lower cf_ref thins the EAIS, about -2.6 m V_af per halving.
- LGM: no grounding-line advance onto the Ross and Weddell shelves or the outer continental shelf, for every cf_ref (area 13.8 against 17-20 M km2; `ant-lgm*-maps.png`). The glacial excess (LGM - PD) is about half the references'.

## North (rounds 1 and 2, cf_ref 0.2)

| run | t [kyr] | V_af [m] | rate [m/kyr] | A [1e6 km2] | NA | EUR | GRL | H bias overlap vs PaleoMIST 20 ka [m] |
|---|---|---|---|---|---|---|---|---|
| pd | 0-30 | 8.37 (GEBCO 7.36) | 0.00 | 2.21 (2.07) | 0.17 (0.07) | 0.07 (0.03) | 8.13 (7.26) | - |
| lgm | 30-60 | 66.14 | 1.38 | 17.22 | 37.0 | 14.8 | 14.3 | -43 |
| lgm-ext | 60-90 | 95.55 | 0.57 | 20.92 | 63.7 | 15.8 | 16.1 | +404 |
| lgmrecon | 0-20 | 105.0 (t 0) -> 126.44 | 0.71 | 24.67 | 85.2 | 24.8 | 16.5 | +724 |
| LGM refs | | 88.5-101.1 | | 25.0-35.9 | 60.1-68.5 | 15.6-18.4 | 12.3-14.5 | |

Series `nor-lgm-regions.png`, `vsle-series.png`; maps `nor-*-maps.png`.

- The LGM climate sustains an LGM-sized NH: from PaleoMIST 20 ka the ice grows, 21 m V_af in 20 kyr (NA +15, EUR +5.5, GRL +1.3), still at +0.5-0.7 m/kyr. Margins stay near the PaleoMIST ones (GL distance median 45 km), the interiors thicken by 700 m on average.
- From PD the NH reaches 95.5 m after 60 kyr, still growing (0.6 m/kyr, NA). EUR is level at 15.8 (lower end of the references), GRL at 16.1 (above them), NA at 63.7 (in the range, still growing).
- The iteration-1 deficit (66 m) is slow build-up, not a climate that cannot hold the ice. The equilibrium under frozen 21 ka forcing is larger than the references (above 126 m), most of it in NA and EUR, and thicker. The two runs have not met yet; the end state from PD may still differ from the one from PaleoMIST (hysteresis).

## Round 3

### North cf_ref 0.1 (pd, lgmrecon)

| run | cf_ref | V_af [m] | dV_af | rate [m/kyr] | A [1e6 km2] | NA | EUR | GRL | H bias / rms overlap [m] |
|---|---|---|---|---|---|---|---|---|---|
| pd | 0.2 | 8.37 | +1.01 (GEBCO) | 0.00 | 2.21 | 0.17 | 0.07 | 8.13 | +104 / 378 |
| pd | 0.1 | 8.05 | +0.69 | 0.00 | 2.20 | 0.17 | 0.07 | 7.82 | +52 / 379 |
| lgmrecon (20 kyr) | 0.2 | 126.44 | +25.4 (PaleoMIST 20 ka) | 0.71 | 24.67 | 85.2 | 24.8 | 16.5 | +724 / 1012 |
| lgmrecon (20 kyr) | 0.1 | 121.77 | +20.7 | 0.68 | 24.95 | 82.6 | 23.4 | 15.8 | +613 / 929 |
| refs | | | | | | PD 0.07 / LGM 60.1-68.5 | 0.03 / 15.6-18.4 | 7.26 / 12.3-14.5 | |

- PD: 0.1 is better (V_af +0.69 against +1.01 m, GRL +0.56 against +0.87, overlap bias +52 against +104 m), margin unchanged (GL median 0 km), converged.
- lgmrecon: halving cf_ref takes 4.7 m off the 21 m growth (NA -2.6, EUR -1.4, GRL -0.7); growth from t 0 is 17.0 m against 21.4 m. The overgrowth under frozen 21 ka forcing is mostly SMB (climate), about a fifth friction.
- Not run: lgm from pd at 0.1 (build-up from PD with less friction).

### ANT LGM: reduced shelf melt (cf_ref 0.05)

Run `Antarctica-lgm-cf05-melt01`: lgm from the cf05 pd bundle, `marine_shelf_ant.gamma_quad_nl` 14500 -> 1450 (x 0.1). Knob: in `quad-nl` the shelf melt is linear in gamma, so it scales the shelf melt by one factor and leaves the ocean temperature (`term_ocean`, dTo) and the thermal forcing untouched. The deep-ocean rate (`c_deep`) and the grounding-zone terms are separate parameters and stay, which is why the mean melt over floating ice drops by a third only (median 0 in both). The LGM ocean term is not a clean knob: dTo = 0.25 x the mean air anomaly enters the thermal forcing non-linearly and through the freezing point.
Shelf budget of the final states (`tools/levante/calib_shelf_diag.jl`, `$G/calib/out/ant-shelf-diag.log`; means over floating ice):

| run | ice area [M km2] | grounded | floating | bmb_shlf [m/yr] (Gt/yr) | cmb [m/yr] (Gt/yr) | tf_shlf mean [K] | dTo [K] |
|---|---|---|---|---|---|---|---|
| pd cf05 | 14.50 | 12.32 | 2.17 | -0.87 (-1734) | -0.49 (-982) | 1.62 | 0 |
| lgm cf05 | 14.77 | 13.84 | 0.93 | -0.30 (-253) | -1.81 (-1547) | 0.52 | -1.09 |
| lgm cf05 melt x 0.1 | 14.78 | 13.84 | 0.93 | -0.20 (-167) | -1.85 (-1584) | 0.52 | -1.09 |
| pd 0.2 / lgm 0.2 | 14.40 / 14.77 | 12.18 / 13.78 | 2.23 / 0.99 | -0.90 / -0.32 | -0.39 / -1.50 | 1.64 / 0.54 | 0 / -1.09 |

- **No advance with less melt**: grounded area 13.843 against 13.836 M km2, V_af 67.94 against 67.92 m, the maps are indistinguishable (`ant-lgm-cf05-melt01-maps.png`). At the LGM the melt is already small (the LGM ocean cools the shelf, tf 1.6 -> 0.5 K) and calving is six times larger.
- **The ice margin is pinned at its PD position.** The total ice area is 14.4-14.8 M km2 in every PD and LGM run; at the LGM the grounding line advances inside the PD calving front and the shelves shrink (2.2 -> 0.9 M km2). `[ycalv_ant]` has `calv_flt_method = calv_grnd_method = "equil"` with the level set: the calving rate is minus the ice velocity, the front does not move. The reference LGM sheets (17-20 M km2 grounded) lie beyond that front, so no friction, melt or sea level setting can reach them. (`[ycalv_nor]`: `vm-m16`, a moving front.)

## Round 4: ANT moving calving front (cf_ref 0.05)

pd (cold start, 30 kyr) and lgm (restart, 30 kyr) with `[ycalv_ant]` set as `[ycalv_nor]`: `calv_grnd_method = "zero"`, `H_min_grnd = H_min_flt = 0`, `tau_ice_flt = tau_ice_grnd = 85e3`, and `calv_flt_method = "vm-m16"` (TAG `-cf05-vm`) or, as the one alternative, `"threshold"` (`Hc_ref_flt` 200 m as in the file; TAG `-cf05-thr`). Overrides on `calib.sh`, every other key as round 2 (`equil`, TAG `-cf05`).
Ice front (`calib_shelf_diag.jl`, `$G/calib/out/ant-shelf-diag-calving.log`): ice area H > 1 m, grounded and floating; mismatch = model-only / BedMachine-only ice area. Maps: dashed lines on the model panel are the ice fronts (`ant-*-cf05{,-vm,-thr}-maps.png`).

| calving | PD V_af [m] (56.83) | PD grounded A (12.30) | PD ice area (14.41) | PD floating (2.10) | PD front mismatch | PD GL median / p90 [km] | PD rate |
|---|---|---|---|---|---|---|---|
| equil | 62.88 (+6.05) | 12.32 | 14.50 | 2.17 (H 261 m) | 0.32 / 0.23 | 0 / 32 | 0.03 |
| vm-m16 | 61.40 (+4.57) | 12.20 | 14.52 | 2.32 (H 158 m) | 0.92 / 0.80 | 0 / 45 | 0.02 |
| threshold | 63.07 (+6.24) | 12.35 | 13.96 | 1.61 (H 304 m) | 0.47 / 0.92 | 0 / 45 | 0.04 |

| calving | LGM V_af [m] (69.3-71.8) | LGM grounded A (17.0-20.5) | LGM ice area | LGM - PD V_af (12.5-15.0) | GL vs PaleoMIST 20 ka median / p90 [km] | overlap bias vs PaleoMIST [m] | LGM rate [m/kyr] |
|---|---|---|---|---|---|---|---|
| equil | 67.92 | 13.84 | 14.77 | 5.0 | 143 / 421 | +82 | 0.00 |
| vm-m16 | 63.15 | 13.34 | 15.51 | 1.8 | 172 / 699 | -44 | -0.02 |
| threshold | 73.69 | 14.87 | 16.26 | 10.6 | 91 / 298 | +254 | +0.15 |

- **vm-m16**: PD grounded state close to equil (V_af 1.5 m lower), but thin shelves (mean 158 m) and a rearranged front (Ross, Ronne fronts retreat, other coasts advance; mismatch 0.9 / 0.8 M km2). At the LGM the front moves out (ice area 15.5) but the grounded ice retreats in the Weddell sector (Filchner-Ronne inland, -1500 m) after a peak at 35 kyr (V_sle 65.5 -> 64.5): less grounded area and V_af than equil, 7 m below PaleoMIST.
- **threshold**: PD grounded state as equil (V_af +6.2, area 12.35, GL unchanged), shelves smaller (1.6 against 2.1 M km2; Ross and Ronne fronts behind BedMachine). At the LGM the grounding line advances, the Weddell embayment grounds, Ross does not (`ant-lgm-cf05-thr-maps.png`): grounded area +1.0 M km2 on equil, GL median 91 km against 143, LGM - PD V_af doubles (10.6 m against 5.0). V_af 73.7 is 1.9-4.4 m above the references and still growing (0.15 m/kyr); the interior is thicker (+254 m overlap).
- Not tuned: `Hc_ref_flt` (threshold) or `tau_ice` (vm-m16).

## Round 5: ANT vm-m16 tuning (on hold)

Decision (Alex, 2026-10-31): ANT calving is the level set with `vm-m16`, tuned; replaces the round-4 `threshold` choice. **On hold** (2026-10-31): the ANT friction setup of rounds 1-4 is off (`docs/dev/setup-diagnosis.md` on main); the sweep is redone on the corrected setup. The four PD runs were cancelled before they started (empty folders `Antarctica-pd-cf05-vm{60,120,170,240}` in the run root).

Knobs (yelmo 30361f96 of S9, `calc_ytopo_calving_lsf` in `yelmo_topography.f90`, `physics/calving/calving_ac.f90`):
- `vm-m16`: calving rate c = -u max(0, tau_vm / `tau_ice_flt`) on the floating front faces (Morlighem et al. 2016, Eq. 4), tau_vm = sqrt(3) B eps_t^(1/n) the tensile von Mises stress from the positive principal strain rates (Eq. 6-7), B = ATT_bar^(-1/n). The front moves at u (1 - tau_vm / tau_ice_flt): it stands where tau_vm = tau_ice_flt and advances where tau_vm is lower (thin, slowly spreading ice). One knob: `tau_ice_flt`.
- ATT_bar carries the enhancement factor, so `ymat_ant.enh_shlf` (0.5) scales tau_vm as enh^(-1/3): on floating ice it is the same knob as tau_ice_flt (not tuned).
- Marine-grounded faces (f_grnd_acx > 0, not land-land): `calv_grnd_method`; `"vm-m16"` exists (same law with `tau_ice_grnd`). Kept `"zero"` (tau_ice_grnd unused): the PD states have 41-48 marine-grounded front cells against 530-820 floating ones, and at the LGM a grounded law can only hold back the advance.
- Unused with the level set: `calv_tau`, `calv_thin`, `Hc_ref_*`, `k2`, `kt_*`, `sd_*`, `w2` (tau_eff only). Used, kept: `H_min_flt/grnd` 0 (thin-margin removal), `H_min_tau`, `[ytopo_ant] front_subgrid = "marine"`, `front_H_eff_min` 50 m, `lsf_method = "snap"`, `dt_lsf = -1`.

Planned range, `tau_ice_flt` 60 / 120 / 170 / 240 kPa (log steps of sqrt 2 around round 4's 85; overrides as round 4, `tools/levante/calib.sh` with `ycalv_ant.calv_flt_method="vm-m16"`, `calv_grnd_method="zero"`, `H_min_grnd=H_min_flt=0.0`, `tau_ice_flt=<t>e3`, `tau_ice_grnd=85.0e3`, TAG `-cf05-vm<t>`). tau_vm at the ice front of the round-4 PD states (`calib_shelf_diag.jl`, `front`, last restart; `$G/calib/out/ant-front-r4.log`), p10 / p25 / p50 / p75 / p90 [kPa]:

| PD run | all front cells | Ross | Filchner-Ronne | Ross / F-R model-only, BedMachine-only [M km2] |
|---|---|---|---|---|
| equil (front near BedMachine) | 46 / 67 / 97 / 137 / 174 | 29 / 40 / 59 / 86 / 96 | 50 / 56 / 99 / 134 / 140 | 0.02, 0.03 / 0.02, 0.00 |
| vm-m16 85 kPa | 24 / 43 / 87 / 185 / 232 | 16 / 20 / 42 / 77 / 152 | 27 / 45 / 107 / 157 / 188 | 0.06, 0.11 / 0.06, 0.27 |
| threshold | 89 / 105 / 126 / 150 / 181 | 88 / 97 / 106 / 129 / 153 | 84 / 104 / 114 / 127 / 139 | 0.00, 0.32 / 0.02, 0.23 |

- At a front near BedMachine (equil) tau_vm spans 46-174 kPa (p10-p90), median 97; Filchner-Ronne needs about 100-135 kPa to hold, Ross about 60-90. With 85 kPa the Filchner-Ronne front retreats (0.27 M km2 BedMachine-only) and thin unconfined shelves on other coasts advance (low tau_vm). The range brackets the BedMachine front stress; a single tau_ice_flt trades the Filchner-Ronne retreat against the advance elsewhere, which is the question of the sweep.

## Recommendation for the science par

- **Antarctica: cf_ref 0.05** (round 2; in the science variants from round 3). Halves the PD excess (+6.0 m against +11.2, overlap bias +207 m against +406), grounding line unchanged, converged. LGM V_af 67.9 m is 1.4-3.9 m below the references on a 3.2 M km2 smaller area: with the shelves grounded it would be in range. The LGM area deficit does not depend on cf_ref (13.78-13.84 M km2) nor on shelf melt (round 3): the calving front is fixed at its PD position (`equil` calving), so it is not a friction problem and should not be tuned with cf_ref (round 4: with a moving front the cf_ref 0.05 LGM V_af is 73.7 m, above the references). Lower than 0.05 would thin the PD further but the remaining PD pattern (EAIS thick, WAIS thin) points to the elevation scaling or the climate, not cf_ref alone.
- **North: cf_ref 0.1** (round 3; replaces "keep 0.2" of round 2). Better PD (+0.69 m, GRL +0.56 m, overlap bias halved), slightly less LGM overgrowth, no downside seen. The LGM overgrowth is mostly SMB and is not a reason to go lower. In the science variants from round 4.

- **ANT calving: threshold** (round 4; decided 2026-10-31, in the science variants from round 5; to be replaced by tuned vm-m16, round 5 on hold). The science case is the glacial cycle, so the LGM extent decides: threshold is the only setting with LGM advance (Weddell grounded, LGM - PD 10.6 m against 5.0 with equil), and the PD grounded ice is as with equil. vm-m16 is worse than equil at the LGM. `[ycalv_ant]` in the generator: `calv_flt_method = "threshold"` (`Hc_ref_flt` 200 m of the file), `calv_grnd_method = "zero"`, `H_min_grnd = H_min_flt = 0`, `tau_ice_flt = tau_ice_grnd = 85e3` (unused by threshold; as `[ycalv_nor]`), the setting of run `Antarctica-*-cf05-thr`.
  **PD cost** (accepted, documented): PD shelves 1.61 against 2.10 M km2 in BedMachine (-24 %), ice area 13.96 against 14.41, the Ross and Ronne fronts behind BedMachine (model-only / BedMachine-only ice area 0.47 / 0.92 M km2, equil 0.32 / 0.23). At the LGM V_af is 1.9-4.4 m above the references and not converged (0.15 m/kyr); Ross does not ground.

## Open questions

1. ANT: Ross does not ground at the LGM with any calving setting; the PD shelves are small with threshold (`Hc_ref_flt` not tuned).
2. Frozen 21 ka forcing overstates the LGM SMB over the glacial cycle; is the lgmrecon drift an acceptable criterion, or should the North be judged in a transient run (-80 ka start of the science par)?
3. GRL is above the LGM references in all North LGM runs (15.8-16.5 against 12.3-14.5 m) and above GEBCO at PD (+0.56 m at 0.1): `[ytill_nor]` scaling for Greenland?
