---
title: "Tools"
---

The scripts in [tools/levante/](https://github.com/fesmc/kryos/tree/main/tools/levante) and the offline programs in
[tests/](https://github.com/fesmc/kryos/tree/main/tests) that the plastic-ice and VILMA2 work uses. Each one-line
description comes from the file's header, which also gives the full command line. The Julia environments
(`<julia-plot>`, `<tools>`) are described in [Analysis and figures](experiment-analysis.md#environment). Everything
runs on Levante.

## Setup

| File | What |
|---|---|
| [pm_prep.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_prep.jl) | Builds the inputs of the standalone `vilma.x` from the FesmData PaleoMIST product (reference, forcing variants A, B, S and the pm-diag variants). Subcommands: `budget` (ice volume against Gowan's), `margins` (registration against Gowan's polygons), `rslref` (the next reference correction), `taboo`, `selen` (PaleoMIST's SELEN sea level from its Stokes coefficients), `visc1d` |
| [pm_visc1d.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_visc1d.sh), [pm_visc1d.f90](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_visc1d.f90) | Computes the 1-D radial viscosity profile that VILMA2's `pre_spinup_1d` uses for a 3-D run, for `EARTH=bm` of `pm_run.sh` (build against `libvilma.a`, run, then `pm_prep.jl visc1d`) |
| [gia_refcorr.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/gia_refcorr.sh), [gia_refcorr.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/gia_refcorr.jl) | Computes the present-day reference correction r of a glacial-cycle run (standalone VILMA2 with the Earth of a kryos par file, iterated with `i_eq = 3` until its 0 ka bed is today's), then writes the kryos reference files (`_REFCORR.nc`, the domain topography with `z_bed + r`) |
| [pm_iterref.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_iterref.sh) | pm-ref: iterates the present-day reference of the standalone B1 run until its modelled 0 ka bed equals today's |
| [pm_selen.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_selen.sh) | Computes PaleoMIST's own SELEN sea level at points with Gowan's `calc_sl`, as a check of `pm_prep.jl selen` |
| [recon_t0.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/recon_t0.jl) | t0 check of the PaleoMIST start (`par/kryos_NH-ANT_vilma2.nml`): the ice of each domain at its first record against its PaleoMIST history at 80 ka |
| [vilma_link_test.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/vilma_link_test.sh) | Link test of `libvilma.a`: VILMA2's host-API tests built against the archive with the host link line |

## Runs

| File | What |
|---|---|
| [run_vilma2.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/run_vilma2.sh) | Stages and submits a run of the science case with `geo = vilma2` (stack S6, kryos built with `vilma=1`): `cont`, `cont2`, `chunk`, `mapgen`, with `set_nml` overrides |
| [pm_run.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_run.sh) | Runs the standalone VILMA2 (`vilma.x` of the S11 stack) under the PaleoMIST ice history: `RUN`, `LMAX`, `EARTH` (1d, l80, bm, 3d, t3), `FORCE` (A, B, S and the pm-diag variants) |
| [plastic_run.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_run.sh) | One kryos run of a plastic par file: stages the run folder (driver, data links, par file with the factor and scale files set, `group.key=value` overrides) and submits its Slurm job |
| [set_nml.py](https://github.com/fesmc/kryos/blob/main/tools/levante/set_nml.py) | Sets existing keys of a fully expanded par file in place (`group.key=value`) |
| [taub_fit/](https://github.com/fesmc/kryos/tree/main/tools/levante/taub_fit) | The Stage 1 fit namelists of the released factor files (`fit_NH.nml`, `fit_ANT.nml`, `fit_PAT.nml`) and their job `fit.sh` |
| [plastic_esmda.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_esmda.sh) `member` | One calibration or ensemble member: the frozen `kryos.x` on the member's `par.nml`, with `status` (exit code, seconds) |
| [jl.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/jl.sh) | Runs a Julia script of the tools project as a small job on the shared partition |
| [tests/plastic_slices.f90](https://github.com/fesmc/kryos/blob/main/tests/plastic_slices.f90) | Offline check against PaleoMIST: the plastic thickness at every slice on PaleoMIST's own bed, margin and tau_b (optionally with the scale and domain-factor files), against its `H_ice`; `plastic=1` |

## Calibration

| File | What |
|---|---|
| [tests/plastic_taub_fit.f90](https://github.com/fesmc/kryos/blob/main/tests/plastic_taub_fit.f90) | Stage 1: one tau_b factor per Gowan domain, fitted to the present-day thickness by Levenberg-Marquardt with a smoothness penalty, for a list of lambda; `plastic=1` |
| [plastic_taub_fit_report.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_taub_fit_report.jl) | Stage 1 report: L-curve and choice of lambda per grid, the `taub_domain_factor.nc` of all grids, present-day thickness per region before and after, maps |
| [plastic_taub_scale.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_taub_scale.jl) | Builds the scale groups of Gowan's domains from the `[scale]` table of a config, and the multiplier files (`taub_scale_file`) |
| [plastic_esmda.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_esmda.jl) | ES-MDA calibration of tau_b against the GAPSLIP 2.1 index points, VILMA2 in the loop: `init`, `update`, `selftest` |
| [plastic_esmda.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_esmda.sh) | The jobs of an ES-MDA calibration: `start` (init, then iteration 1), `submit` (an iteration's member array and its analysis) |
| [plastic_esmda.toml](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_esmda.toml) | Config of calibration v1: NH-32KM, 3 groups, independent model error |
| [plastic_esmda_v2.toml](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_esmda_v2.toml) | Config of calibration v2: NH-32KM and ANT-16KM, Stage 1 factors, present day pinned, correlated observation error |
| [plastic_ensemble.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_ensemble.jl) | Weighted ensemble: a Latin hypercube sample of the prior (`init`), and per-group correlated likelihood weights and ESS (`weights`) |
| [plastic_ensemble.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_ensemble.sh) | The jobs of a weighted ensemble: `start` (init, then wave 1), `submit` (a wave's members and the weights after them) |
| [plastic_ensemble_v1.toml](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_ensemble_v1.toml) | Config of ensemble v1, wave 1: NH-32KM, ANT-16KM and PAT-PAL-8KM on the 3-D Earth, 8 groups, per-group model error |
| [tests/plastic_closure_fit.f90](https://github.com/fesmc/kryos/blob/main/tests/plastic_closure_fit.f90) | Stage 1b of the implied mass balance: the flow closure (`yield` or `coulomb_n` parameter lists) against the observed present-day speed; with `prior_file`, the joint inversion on one set; `plastic=1` |

## Scoring

| File | What |
|---|---|
| [pm_compare.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_compare.jl) | The scorer: VILMA2 runs (standalone, or a kryos run's `geo_glob_vilma.nc`) against PaleoMIST and the GAPSLIP 2.1 indicators, by groups and marine-bed areas; included by the calibration and plot tools |
| [plastic_score.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_score.sh) | Scores one kryos run's VILMA2 stream against GAPSLIP 2.1 with `pm_compare.jl`, beside B1 and I910_it0 |
| [plastic_pd_regions.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_pd_regions.jl) | Present-day ice thickness of the plastic backend against the FesmData topography hubs, per region |
| [plastic_residual_corr.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_residual_corr.jl) | Correlation of the index-point misfit in distance and age, and the model error sigma_m, per group, region and sector |
| [plastic_earth_compare.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_earth_compare.jl) | Plastic runs on different Earths: GAPSLIP scores by group, region, marine-bed area and site, site maps, RSL curves, volume and ESL |

## Plots and reports

| File | What |
|---|---|
| [plastic_lgm_map.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_lgm_map.jl) | Global maps of a plastic run after Gowan et al. (2021, Fig. 2): `mode=age`, `margins`, `smb` |
| [plastic_mapkit.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_mapkit.jl) | The drawing module of the global maps (Equal Earth projection, VILMA2 and instance fields, coastlines, margins, scored sites) |
| [plastic_ensemble_sd_map.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_ensemble_sd_map.jl) | sd of RSL and ice thickness across ensemble members at one age, with that age's index points |
| [plastic_ensemble_spread.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_ensemble_spread.jl) | Surface-elevation spread of a calibration iteration: all members against the best by chi2/n |
| [plastic_esmda_report.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_esmda_report.jl) | A calibration's report: scores per iteration, factor curves, ice volume, present-day thickness per region, offline Greenland check |
| [plastic_esmda_diag.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_esmda_diag.jl) | A calibration's diagnostics: volume by group and ESL, prior-ensemble sensitivity, residuals by age, Antarctic sites, ANT maps |
| [plastic_implied_pd.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_implied_pd.jl) | Implied mass balance of a run at present day against MAR, RACMO and the observed surface speeds |
| [plastic_closure_maps.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_closure_maps.jl) | Maps of the implied SMB and the speed misfit of `plastic_closure_fit` sets (`:inv`: the inversion), with totals and ablation area |
| [plastic_maps.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_maps.jl) | Maps of the offline check (`plastic_slices`): H - PaleoMIST H_ice at given ages |
| [pm_maps.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_maps.jl) | Figures of the standalone VILMA2 PaleoMIST test from `pm_compare.jl` output |
| [gia_refcorr_plot.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/gia_refcorr_plot.jl) | The figure of a `gia_refcorr.sh` case: the converged r and the 0 ka bed misfit per iteration |

## Tests of the backend

| File | What |
|---|---|
| [tests/test_plastic_margin.f90](https://github.com/fesmc/kryos/blob/main/tests/test_plastic_margin.f90) | Margin chronology and tau_b (distance transform, interpolation, domains, scale groups, domain factors); runs with `plastic=0` |
| [tests/test_ice_plastic.f90](https://github.com/fesmc/kryos/blob/main/tests/test_ice_plastic.f90) | The backend through the ice wrappers on a toy grid (profile, masks, exports, restart, flotation); `plastic=1` only |

## Other scripts

- **Historical gate checks** (the kryos-against-yelmox and phase gates, see
  [Developer notes](developer-notes.md)): `compare_*.sh`, `gate_S*.sh`, `cmp_*.sh`, `cmpjob_*.sh`, `g2a*`, `g2b*`,
  `g2c*`, `g34_vilma2.sh`, `g3_budget.jl`, `g5*`, `g6*`, `g7*`, `hub_offline.jl`, `round_check.jl`, `register_k3.sh`,
  `run_kryos.sh`, `run_kryos_kk.sh`, `submit_ref.sh`, `submit_wp0.sh`, `check_kryos.sh`, `build_stack.sh`,
  `wp0_*.py`, `wp4_hub_count.py`.
- **Yelmo-backend calibration and diagnosis** (not used by the plastic work): `calib*`, `ant_friction*`, `setup_diag*`.
