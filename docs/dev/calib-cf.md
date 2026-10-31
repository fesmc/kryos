# calib-cf: basal friction (cf_ref) for the science par, NH-32KM and ANT-32KM

`$G = /work/ba1442/robinson/kryos-gate`. Branch `calib-cf`. Stack S9 (kryos 624f686, `kryos-v0`: vilma off, FastIsostasy), 64 threads, exclusive node.
Runs `$G/runs/calib/kx-624f686/S9/<Domain>-<case><TAG>`, metrics `$G/calib/out/<ant|nor>-<case><TAG>.md`, maps `$G/calib/maps/`.

## Method

No optimisation: Yelmo's elevation scaling of cf_ref, one cf_ref per domain, judged on the present-day (PD) and LGM states.

- `[ytill_ant]`: `scale_zb = 2` (exp), z0/z1 = -400/0 m, cf_min 0.01. `[ytill_nor]`: exp, -100/100 m, cf_min 0.001.
- Science variants `par/kryos_NH-ANT_calib.nml`, `par/kryos_NH-ANT_vilma2_calib.nml` (`tools/bin/kryos-calib-par.jl sci`): rounds 1-2 cf_ref 0.2 both domains; from round 3 ANT 0.05, North 0.2 (`CF_ANT`, `CF_NOR` in the generator). `par/kryos_NH-ANT_vilma2.nml` unchanged.
- `par/calib/*.nml` are regenerated from the science variant, so they carry ANT 0.05 from round 3; the round 1-2 ANT runs had 0.2 (now `ytill_ant.cf_ref=0.2`).
- Single domain, FastIsostasy, frozen forcing (`tools/bin/kryos-calib-par.jl run`, pars `par/calib/`):
  - pd: t 0 -> 30 kyr, alpha 0, sea level 0, orbit 0 ka, cold start from the hub reference.
  - lgm: t 30 -> 60 kyr, restarted from pd, alpha 1, sea level -120 m, orbit -21 ka.
  - lgmrecon (North): lgm forcing, cold start from PaleoMIST (`recon_nor`, whole domain, no bed filter, smoothing or float cut, as the science par at -80 ka; `recon_time = -21e3` selects the 20 ka record), t 0 -> 20 kyr.
  - lgm-ext (North): lgm restarted from the 60 kyr bundle of lgm, t 60 -> 90 kyr.
- Submit `tools/levante/calib.sh <Domain> <case> [set_nml ...]` (TAG, FROM, RSTT); cf_ref variants by `ytill_ant.cf_ref=...` / `ytill_nor.cf_ref=...`. Metrics `tools/levante/calib_metrics.sh` -> `calib_metrics.jl`; maps `tools/calib_maps.jl` (CairoMakie).
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

## Recommendation for the science par

- **Antarctica: cf_ref 0.05** (round 2; in the science variants from round 3). Halves the PD excess (+6.0 m against +11.2, overlap bias +207 m against +406), grounding line unchanged, converged. LGM V_af 67.9 m is 1.4-3.9 m below the references on a 3.2 M km2 smaller area: with the shelves grounded it would be in range. The LGM area deficit does not depend on cf_ref (13.78-13.84 M km2) nor on shelf melt (round 3): the calving front is fixed at its PD position (`equil` calving), so it is not a friction problem and should not be tuned with cf_ref. Lower than 0.05 would thin the PD further but the remaining PD pattern (EAIS thick, WAIS thin) points to the elevation scaling or the climate, not cf_ref alone.
- **North: cf_ref 0.1** (round 3; replaces "keep 0.2" of round 2). Better PD (+0.69 m, GRL +0.56 m, overlap bias halved), slightly less LGM overgrowth, no downside seen. The LGM overgrowth is mostly SMB and is not a reason to go lower. Not yet in the science variant (`CF_NOR` still 0.2), pending agreement.

## Open questions

1. ANT calving: the LGM extent needs a moving front. `calv_flt_method` (and `calv_grnd_method`) other than `equil` for `[ycalv_ant]`, e.g. `vm-m16` as the North or `threshold`; this changes the PD front too and needs its own PD check. Not run (a decision outside this brief).
2. Frozen 21 ka forcing overstates the LGM SMB over the glacial cycle; is the lgmrecon drift an acceptable criterion, or should the North be judged in a transient run (-80 ka start of the science par)?
3. GRL is above the LGM references in all North LGM runs (15.8-16.5 against 12.3-14.5 m) and above GEBCO at PD (+0.56 m at 0.1): `[ytill_nor]` scaling for Greenland?
4. Apply North 0.1 to the science variants (`CF_NOR`)?
