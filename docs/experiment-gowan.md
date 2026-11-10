---
title: "Reproducing Gowan et al. (2021)"
---

The PaleoMIST 1.0 ice history (Gowan et al. 2021, minimal MIS 3 scenario) imposed on
VILMA2 on Gowan's 1-D Earth, 80 ka to present, and the resulting relative sea level (RSL)
scored against PaleoMIST itself and against the GAPSLIP 2.1 sea-level indicators. There is
no ice model: the ice is read from the history. This is run **B1** of
[docs/dev/vilma2-paleomist.md](https://github.com/fesmc/kryos/blob/main/docs/dev/vilma2-paleomist.md),
the reference against which every plastic run of [Simulating our own ice](experiment-plastic.md)
is scored.

- **Program:** standalone VILMA2, `$VILMA2/bin/vilma.x` (`VILMA2` of `tools/levante/env.sh`,
  default the vilma2 checkout of kryos; `make vilma` in it), not `kryos.x`. The published runs
  used `vilma.x` of the gate stack S11 (vilma2 c17d513, tree identical to 523b660). No kryos par
  file imposes the history everywhere: the namelist `pm.nml` is written by the job script from its
  environment.
- **Scripts** (`tools/levante/`):
  [pm_prep.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_prep.jl) (inputs),
  [pm_run.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_run.sh) (one run, submits itself),
  [pm_compare.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_compare.jl) (metrics and GAPSLIP scores),
  [pm_maps.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_maps.jl) (figures),
  [pm_iterref.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_iterref.sh) (iterated present-day reference).
- **Folders** (`env.sh`): your runs in `PM_DIR` (default `$KRYOS_WORK/paleomist`; `ROOT` of
  `pm_run.sh` overrides the parent folder of a run), the inputs in `PM_DATA` (default
  `$PM_DIR/data`, else the reference inputs `$PM_REF/data`). `PM_REF` (default
  `$KRYOS_DATA/paleomist`) holds the published runs `B1_l170`, `ref/I910_it0` and `B3_l170` (their
  `vilma_out.nc`, `pm.nml`, `run.log`) and the inputs `data/`, read-only: the references of every score.

Within kryos the same history enters every `geo = vilma2` par file as the ice outside the
domains' footprints (`[globe_ice] ext = "file"`, see [Global hub](global-hub.qmd)); the
coupling of a replayed PaleoMIST load through the kryos globe was tested by gate G2b
(`tests/g2b_replay.f90`, `tools/levante/archive/g2b.sh`), not as an experiment.

## What is imposed

| | Gowan et al. (2021), SELEN | B1, VILMA2 |
|---|---|---|
| Ice | ICESHEET 5 km grids, 2.5 kyr steps, 80-0 ka | FesmData `Global/GLOBAL-0.25DEG/GLOBAL-0.25DEG_IceHistory-PaleoMIST.nc` (`H_ice`), remapped conservatively to the Gauss grid by `vilma.x` (`remap_input = .true.`) |
| Load history | Heaviside steps | the 33 records as steps (`pm_force_B.nc`, `h_ice(lon, lat, time)`): record t_k enters over [t_k, t_k + 1 yr] and is held to t_{k+1} |
| Before 80 ka | equilibrium with the first record | equilibrium with the 80 ka load (spin-up) |
| Reference | present topography as final condition | PaleoMIST's 0 ka `z_bed` and `H_ice` (`pm_ref.nc`), relaxed (`i_eq = 1`) |
| Densities | 910 / 1000 (SELEN `data.inc`) | 931 / 1000 (B1); 910 / 1000 (`I910_it0`) |
| Sea-level equation | shoreline migration, grounding line, rotation | `rotation = .true.`, `deg1_frame = "cm"`, flotation, sloping coast (`sle_subgrid = .true.`) |

The step is a 1-yr ramp because `vilma.x` ramps the ice linearly between two forcing slices;
record t_k is compared with the state at t_k + 1 yr (`pm_prep.jl`, `pm_compare.jl` `trec`).

## Earth model

`pm_run.sh EARTH=1d`: Gowan's 1-D Earth on PREM (`earth = "PREM"`, `l_visc_3d = .false.`,
`pre_spinup_1d = .false.`), 12 layers:

| Layers | r_top [km] | Rheology | eta [Pa s] |
|---|---|---|---|
| lithosphere | 6371, 6291 | elastic (0) | 1e40 |
| upper mantle, 120-670 km | 6251, 6151, 5971, 5771 | Maxwell (1) | 4e20 |
| lower mantle, 670 km to the CMB | 5701, 5600, 4943, 4287, 3630 | Maxwell (1) | 4e22 |
| core | 3480 | fluid (2) | 0 |

So 120 km elastic lithosphere, upper mantle 4e20 to 670 km, lower mantle 4e22 Pa s; rho and
mu are VILMA2's PREM shell means. Spectral resolution `LMAX` 256 or 170 (170 differs by
0.4 m RMS and at most 2.4 m, at half the cost).

Other Earths of the same script, for the 3-D effect and its split: `3d` (the science setting
of `par/kryos_NH-ANT_vilma2.nml`: 11 layers, 80 km lithosphere, Bagge et al. 2021 structure 6,
`pre_spinup_1d`), `l80` (Gowan's viscosities on the 11 science layers), `bm` (the science
Earth collapsed to its 1-D pre-spin-up profile; needs `pm_visc1d.sh` first), `t3` (SELEN's
three-layer TABOO Earth of PaleoMIST).

## Timeline

- Spin-up: VILMA2 relaxes to equilibrium with the 80 ka load (`solid_earth_spinup`: passes of
  200 yr growing by 1.6, until the mean bed speed over a pass is below `RATE_TOL` = 1e-5 m/yr),
  capped at `TMAX` (B1: 2e6 yr; script default 5e5). B1 converges after 12 passes, 93 kyr.
- Transient: the whole forcing file (`time_init = -1.0e30`, `time_end = 1.0e30` in `pm.nml`),
  80 ka to 0 ka + 1 yr; explicit memory scheme `fe`, `cfl` 1, sub-steps of 100 yr on the 1-D
  Earth.
- Output: `vilma_out.nc` has one record per forcing record (65 under `pm_force_B.nc`): -80 ka,
  then each 2.5 kyr slice t and its 1-yr step t + 1, ending at 0 and +1 yr.
- Cost: B1_l170 3 min, B1_l256 6 min on 64 threads (shared partition); B3 (3-D) 1 h 44.

## Running

On Levante, with the environment of `tools/levante/env.sh` and the Julia project of
`tools/levante` (`tools/levante/setup.sh`, see [Getting started](getting-started.md#the-levante-tools)).
Inputs, once (writes `pm_ref.nc`, `pm_force_A.nc`, `pm_force_B.nc`, `pm_force_S.nc`, ...); without
them the runs read the reference inputs `$PM_REF/data`:

```bash
julia --project=tools/levante tools/levante/pm_prep.jl $FESMDATA $PM_DIR/data
```

B1 at lmax 170 and 256 (`pm_run.sh` submits itself: compute partition, one node, 64 threads, 8 h;
`PM_PARTITION`, `PM_TIME`, `NTHR` override; `RHOI` defaults to 910 now, B1 used 931):

```bash
RUN=B1_l170 LMAX=170 EARTH=1d FORCE=B TMAX=2.0e6 RHOI=931.0 tools/levante/pm_run.sh
RUN=B1_l256 LMAX=256 EARTH=1d FORCE=B TMAX=2.0e6 RHOI=931.0 tools/levante/pm_run.sh
```

B1 at 910 kg m-3 with the iterated present-day reference (iteration 0 = `I910_it0`, B1 at 910;
stops when the 0 ka bed misfit is below `TOL` m; submits itself, 3 h per run; runs in
`$PM_DIR/ref`):

```bash
BASE=I910 TOL=1.0 NMAX=8 tools/levante/pm_iterref.sh
```

`EARTH=bm` needs the 1-D profile of a 3-D run first:
`tools/levante/pm_visc1d.sh build` (once, in a login shell with the Intel modules), then
`tools/levante/pm_visc1d.sh $PM_DIR/B3_l170` (writes `$PM_DIR/data/pm_visc_bm.nc`). PaleoMIST's
own SELEN sea level at a few points, the check of `pm_prep.jl selen`:
`CASE=72_73_74_75/ehgr tools/levante/pm_selen.sh` (output in `$PM_DIR/selen/`).

## The switches

Environment of `pm_run.sh`:

| Variable | Values | |
|---|---|---|
| `RUN` | folder name | `$ROOT/$RUN`, must not exist |
| `LMAX` | 170, 256 | spectral degree |
| `EARTH` | `1d`, `3d`, `l80`, `bm`, `t3` | above |
| `FORCE` | `B` (80 -> 0 ka steps), `A` (22.5 ka ice held: equilibrium), `S` (smoke, 25 -> 20 ka); diagnosis variants `Bg`, `Bl`, `Bm`, ... | the forcing file `pm_force_<FORCE>.nc` |
| `TMAX`, `RATE_TOL` | 5e5, 1e-5 | spin-up cap [yr], tolerance [m/yr] |
| `RHOI`, `RHOW` | 910, 1000 | densities of the load |
| `ROT`, `DEG1` | `.true.`, `cm` | rotation, degree-1 frame |
| `REF`, `RSLREF` | `pm_ref`, empty | reference file; an rsl field gives `i_eq = 3` (iterated reference) |
| `NTHR` | 64 | OpenMP threads |
| `ROOT`, `PM_DATA` | `$PM_DIR`, `$PM_DIR/data` (else `$PM_REF/data`) | the parent folder of the run; the inputs |

The runs of the study: B1_l256, B1_l170 (1-D), B3_l170 (3-D), A1, A3 (equilibrium under the
22.5 ka ice), B1_L80, B1_bm (the split of the 3-D effect), `ref/I910_it*` (iterated reference). B1_l170,
B3_l170 and ref/I910_it0 are published in `PM_REF`; run the others with `pm_run.sh` into `PM_DIR`.

## Outputs

`$PM_DIR/<RUN>/` (the published runs: `$PM_REF/<RUN>/`): `pm.nml` (the namelist), `vilma_out.nc`
(`rsl`, `dz_bed`, `bsl`, `bsl_ive`, ... on the Gauss grid), `run.log` (spin-up passes, warnings,
memory), `input/` (links to `$VILMA2/input`), `fesmdata` (link), `maps/` (remapping weights, built
per run); the job log `$PM_DIR/<RUN>.out`.

VILMA2's `rsl` is relative to the relaxed reference and not zero at 0 ka (the Earth is still
rebounding); PaleoMIST's `dz_sl` is S(t) - S(0) by construction. The comparison therefore
takes rsl(t) - rsl(0 ka).

## Scoring and plots

Metrics against PaleoMIST (RMS, max, bias of RSL relative to present over the glaciated areas,
the ocean and the globe, on the 0.5 degree PaleoMIST grid) and the GAPSLIP 2.1 scores:

```bash
julia --project=tools/levante tools/levante/pm_compare.jl $FESMDATA $PM_DIR <out folder> B1_170=$PM_REF/B1_l170 B3=$PM_REF/B3_l170 [<key>=<folder> ...]
```

Keys: `B1` (B1_l256), `B1_170` (B1_l170), `B1_L80`, `B1_bm`, `B3` (B3_l170), `A1`, `A1_170`, `A3`;
each the folder of that name in the runs folder unless given as `key=folder`; a missing run is
skipped. `read_run` reads a standalone run (`vilma_out.nc`) or a kryos run (`geo_glob_vilma.nc`)
in place, so a kryos run can stand in for a key (this is how `plastic_score.sh` works).

Indicators: FesmData `Records/RSL/RSL-GAPSLIP_global.nc` (GAPSLIP 2.1, Gowan 2023), the primary
data (`alternative = 0`) of 0-80 ka: 8825 samples, 4752 index points. The model RSL relative to
present at each sample, bilinear in space, linear in time between the 2.5 kyr records.

- Index points: z = (model - obs) / sigma, sigma^2 = sigma_z^2 + (dRSL/dt sigma_t)^2, sigma_z
  half the 2-sigma elevation range on the side of the misfit, sigma_t half the 2-sigma age
  range. Reported: RMS of z, median |z|, fraction |z| <= 2, mean misfit.
- Limiting data: violated when no model RSL over the 2-sigma age range reaches the 2-sigma
  elevation range (marine: max < `rsl_min`; terrestrial: min > `rsl_max`).
- Groups of Gowan's regions: Laurentide (North America Arctic and Atlantic), Greenland,
  Fennoscandia/Barents (Europe, Eurasian Arctic), Antarctica, far field; and the marine-bed
  areas Hudson Bay, Hudson Strait, Ross Sea, James Ross Island.

Writes `pm_metrics.md`, `pm_metrics_B1_PM.csv`, `pm_sites.csv`, `pm_gmsl.csv`, `pm_fields.nc`,
`pm_indicators.md`, `pm_ind_samples.csv`, `pm_ind_sites.csv`, `pm_ind_curves.csv`. Figures
(CairoMakie) from that folder:

```bash
julia --project=tools/levante tools/levante/pm_maps.jl <compare folder>
```

`pm_B1_vs_PM.png`, `pm_effects.png`, `pm_decomp.png`, `pm_ind_sites.png`, `pm_ind_map.png`,
`pm_A_minus_B.png`, `pm_rsl0.png`, `pm_sites.png`, `pm_gmsl.png` ([Analysis and figures](experiment-analysis.md)).

### Results

From [docs/dev/vilma2-paleomist.md](https://github.com/fesmc/kryos/blob/main/docs/dev/vilma2-paleomist.md)
(figures in `docs/dev/figures/vilma2-paleomist/`):

- B1 reproduces PaleoMIST to about 10 m RMS globally and 24-29 m RMS over the formerly
  glaciated areas (22.5-10 ka); the far field is 5-8 m less deep at the LGM.
- Against GAPSLIP, B1 scores as PaleoMIST itself: median |z| 4.4 for both, 30-31 % of the index
  points within 2 sigma, RMS z 30.9 (B1_l256; B1_l170 30.4) against 29.5.
- Most of the remaining difference is the sea-level equation (PaleoMIST's SELEN keeps the water
  that grounded marine ice displaces), then the reference state and the ice density; B1 is not
  rebased on them.

## The bridge: PaleoMIST from the plastic ice

The plastic backend (see [PlasticIceSheet](backend-plastic.md)) rebuilds PaleoMIST's ice online
instead of reading it: at every 100 yr step a perfectly-plastic ice sheet on Gowan's margins
and tau_b domains, on the bed and sea level of VILMA2.
[par/kryos_North_plastic_vilma2_b1.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_North_plastic_vilma2_b1.nml)
is that configuration on B1's Earth, scored beside B1:

- ice: `ice_nor` (backend `plastic`, NH-32KM), margins `f_ice` > 0.5 of
  `{grid_name}_IceHistory-PaleoMIST.nc` interpolated by the signed distance (`margin_interp =
  "sdist"`); tau_b of Gowan's domains (`taub_method = "domains"`,
  `{grid_name}_IceHistory-PaleoMIST_taub.nc`, time dependent, `taub_default` 5000 Pa), unscaled
  (`taub_domain_factor_file = ""`, `taub_scale_file = ""`);
- Earth: the `[vilma]` group of `pm_run.sh EARTH=1d, LMAX=170`, the 80 ka relaxation capped at
  2e6 yr;
- reference: relaxed, no correction r: PaleoMIST's 0 ka bed and ice as the background X
  (`[X_z_bed]`, `[X_H_ice]`), today's topography (GEBCO2026) on NH-32KM;
- the rest of the world: the PaleoMIST history (`[globe_ice] ext = "file"`);
- densities: the run's, 910 and 1028 (B1: 931, 1000).

With Gowan's tau_b the plastic NH ice is 7-10 % larger in volume than PaleoMIST (offline check,
mostly the resolution); online, ESL at 20 ka is -131.7 m against PaleoMIST's -115.6 m, and over
all index points RMS z 25.9, median |z| 3.95 against 30.4 and 4.41 for B1_l170
([docs/dev/plastic-backend.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-backend.md),
"Baseline"). Scoring a kryos run beside B1:

```bash
tools/levante/plastic_score.sh <kryos run folder> <label> [<score folder>]
```

(a small job it submits itself: `pm_compare.jl` with `B1=$PM_REF/B1_l170 B1_170=$PM_REF/ref/I910_it0
B3=<run>`; output `<score folder>/out`, default `$KRYOS_WORK/score/<label>/out`; the run appears in
the tables under its folder name.) The variants with the Antarctic and Patagonian instances, the 3-D Earth
and the calibrated tau_b: [Simulating our own ice](experiment-plastic.md).
