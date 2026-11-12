---
title: "Tools"
---

The scripts in [tools/levante/](https://github.com/fesmc/kryos/tree/main/tools/levante) and the offline programs in
[tests/](https://github.com/fesmc/kryos/tree/main/tests) that the plastic-ice and VILMA2 work uses. Each one-line
description comes from the file's header, which also gives the full command line. Every script sources `env.sh`;
the Julia tools run with `julia --project=tools/levante` ([Getting started](getting-started.md#the-levante-tools),
[Analysis and figures](experiment-analysis.md#environment)). The model runs (`kryos.x`, `vilma.x`, the offline fits) are
staged and submitted with [runme](https://github.com/fesmc/runme) ([Getting started](getting-started.md)). Scripts that
run as Slurm jobs submit themselves: run them directly, not with `sbatch`. Everything runs on Levante and on AWI's
albedo (`KRYOS_MACHINE`, env.sh); [tools/albedo/stack.sh](https://github.com/fesmc/kryos/blob/main/tools/albedo/stack.sh)
builds a kryos stack on albedo.

## Environment

| File | What |
|---|---|
| [env.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/env.sh) | The environment of every tool: `KRYOS`, `KRYOS_WORK`, `FESMDATA`, `FESMDATA_SOURCES`, `KRYOS_DATA`, `VILMA2`, `PM_DIR`, `PM_REF`, `ACCOUNT`, `KRYOS_PART`, `KRYOS_PART_LONG`, `KRYOS_QUEUE`, `KRYOS_RUNME`, `KRYOS_JULIA`, with the defaults of `KRYOS_MACHINE` (levante, albedo) (overrides in the environment or `~/.kryos.env`); function `kryos_jl` |
| [setup.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/setup.sh) | The one-time setup: checks the paths and the builds; links `fesmdata` and `kryos_data` into the checkout and copies the shared grid files and remapping weights of `$KRYOS_DATA/maps` into its `maps/`; checks runme (>= 0.5.27) and creates `.runme/config.toml`; Julia (juliaup, 1.13) and the Julia project; creates `KRYOS_WORK`; safe to rerun |
| [.runme/](https://github.com/fesmc/kryos/tree/main/.runme) | The checkout as a runme project: `info.json` (executables `kryos`, `taub_fit`, `closure_fit`; defaults `input/kryos_defaults.nml`; the links `input`, `fesmdata`, `kryos_data`, `maps`, `ice_data`, `isostasy_data`), the job templates `submit_slurm`, `submit_slurm_omp` (the job writes `status`, exit code and seconds), `config.default.toml` (`dkrz_levante`, account `ba1442`, 16 threads) |
| [Project.toml](https://github.com/fesmc/kryos/blob/main/tools/levante/Project.toml), `Manifest.toml` | The Julia project of all tools: CairoMakie, NCDatasets, DelimitedFiles, Julia 1.13 |
| [jl.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/jl.sh) | Runs a Julia script of the tools as a small job on the shared partition (1 h, 4 cpus, 16 GB; `JL_TIME`, `JL_MEM`) |

## Setup of the inputs

| File | What |
|---|---|
| [pm_prep.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_prep.jl) | Builds the inputs of the standalone `vilma.x` from the FesmData PaleoMIST product (reference, forcing variants A, B, S and the pm-diag variants). Subcommands: `budget` (ice volume against Gowan's), `margins` (registration against Gowan's polygons), `rslref` (the next reference correction), `taboo`, `selen` (PaleoMIST's SELEN sea level from its Stokes coefficients), `visc1d` |
| [pm_visc1d.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_visc1d.sh), [pm_visc1d.f90](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_visc1d.f90) | Computes the 1-D radial viscosity profile that VILMA2's `pre_spinup_1d` uses for a 3-D run, for the case `pm_bm.nml` (build against `libvilma.a`, run on the namelist of a run of `pm_3d.nml`, then `pm_prep.jl visc1d`) |
| [gia_refcorr.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/gia_refcorr.sh), [gia_refcorr.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/gia_refcorr.jl) | Computes the present-day reference correction r of a glacial-cycle run (standalone VILMA2 with the Earth of a kryos par file, iterated with `i_eq = 3` until its 0 ka bed is today's; the case folder is a PaleoMIST runme project, the iterations `it<n>` runme runs of its `vilma_rc.nml`), then writes the kryos reference files (`_REFCORR.nc`, the domain topography with `z_bed + r`) |
| [pm_iterref.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_iterref.sh) | pm-ref: iterates the present-day reference of the standalone B1 run until its modelled 0 ka bed equals today's (each iteration a runme run of a `paleomist/` case in `ROOT`; `BASE`, `EARTH`, `FORCE`, `LMAX`, `TMAX`, `RHOI`, `PM_P`, `NTHR`, `TOL`, `NMAX`, `ROOT`) |
| [pm_selen.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_selen.sh) | Computes PaleoMIST's own SELEN sea level at points with Gowan's `calc_sl`, as a check of `pm_prep.jl selen` |

## Runs

| File | What |
|---|---|
| [pm_project.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_project.sh) | `pm_project.sh [<folder>] [<data folder>]`: a folder of standalone VILMA2 runs (default `PM_DIR`) as a runme project for `vilma.x` of the VILMA2 checkout, with the links `input` (VILMA2's `input/`), `fesmdata`, `pm_data` (default `$PM_DIR/data`, else `$PM_REF/data`) |
| [paleomist/](https://github.com/fesmc/kryos/tree/main/tools/levante/paleomist) | The runme cases of the standalone VILMA2 under the PaleoMIST ice history, one per Earth: `pm_1d.nml`, `pm_l80.nml`, `pm_bm.nml`, `pm_3d.nml`, `pm_t3.nml` (complete `&ctl` and `&vilma`; the `-p` options in the header) |
| [members.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/members.sh) | Sourced by `plastic_esmda.sh` and `plastic_ensemble.sh`: the members of an iteration or wave as one job array, staged and submitted by runme (`--array THROTTLE`) from the calibration folder; prints the job id |
| [set_nml.py](https://github.com/fesmc/kryos/blob/main/tools/levante/set_nml.py) | Sets existing keys of a fully expanded par file in place (`group.key=value`) |
| [tests/plastic_slices.f90](https://github.com/fesmc/kryos/blob/main/tests/plastic_slices.f90) | Offline check against PaleoMIST: the plastic thickness at every slice on PaleoMIST's own bed, margin and tau_b (optionally with the scale and domain-factor files), against its `H_ice`; `plastic=1` |

## Calibration

| File | What |
|---|---|
| [tests/plastic_taub_fit.f90](https://github.com/fesmc/kryos/blob/main/tests/plastic_taub_fit.f90) | Stage 1: one tau_b factor per Gowan domain, fitted to the present-day thickness by Levenberg-Marquardt with a smoothness penalty, for a list of lambda; `plastic=1` |
| [taub_fit/](https://github.com/fesmc/kryos/tree/main/tools/levante/taub_fit) | The Stage 1 fit namelists of the released factor files (`fit_NH.nml`, `fit_ANT.nml`, `fit_PAT.nml`; `fesmdata/...` relative to the run folder), each one runme run of `taub_fit` |
| [plastic_taub_fit_report.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_taub_fit_report.jl) | Stage 1 report: L-curve and choice of lambda per grid, the `taub_domain_factor.nc` of all grids, present-day thickness per region before and after, maps |
| [plastic_taub_scale.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_taub_scale.jl) | Builds the scale groups of Gowan's domains from the `[scale]` table of a config, and the multiplier files (`taub_scale_file`) |
| [plastic_esmda.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_esmda.jl) | ES-MDA calibration of tau_b against the GAPSLIP 2.1 index points, VILMA2 in the loop: `init`, `update`, `selftest` |
| [plastic_esmda.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_esmda.sh) | The jobs of an ES-MDA calibration: `start` (init, then iteration 1), `submit` (an iteration's member array and its analysis) |
| [plastic_esmda.toml](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_esmda.toml) | Config of calibration v1: NH-32KM, 3 groups, independent model error |
| [plastic_esmda_v2.toml](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_esmda_v2.toml) | Config of calibration v2: NH-32KM and ANT-16KM, Stage 1 factors, present day pinned, correlated observation error |
| [plastic_ensemble.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_ensemble.jl) | Weighted ensemble: a Latin hypercube sample of the prior (`init`), and per-group correlated likelihood weights and ESS (`weights`; RSL index points and GNSS uplift rates); the Earth parameters (`earth-regions`, `earth-shift`) |
| [plastic_ensemble.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_ensemble.sh) | The jobs of a weighted ensemble: `start` (init, then wave 1), `submit` (a wave's members and the weights after them) |
| [plastic_ensemble_v2.toml](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_ensemble_v2.toml) | Config of ensemble v2: v1 plus the Earth ([earth]: global and regional log10-viscosity shifts of the 3-D Earth), the error model of ens-v1, all groups judging |
| [plastic_earth_shift.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_earth_shift.jl) | Included by `plastic_esmda.jl`: the Earth parameters of an ensemble, the region masks (`plastic_ensemble.jl earth-regions`) and the viscosity shift file of a member (VILMA2 `visc_shift_file`; `earth-shift`) |
| [plastic_ensemble_v1.toml](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_ensemble_v1.toml) | Config of ensemble v1, wave 1: NH-32KM, ANT-16KM and PAT-PAL-8KM on the 3-D Earth, 8 groups, per-group model error |
| [tests/plastic_closure_fit.f90](https://github.com/fesmc/kryos/blob/main/tests/plastic_closure_fit.f90) | Stage 1b of the implied mass balance: the flow closure (`yield` or `coulomb_n` parameter lists) against the observed present-day speed; with `prior_file`, the joint inversion on one set; `plastic=1` |

## Scoring

| File | What |
|---|---|
| [pm_compare.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_compare.jl) | The scorer: VILMA2 runs (standalone, or a kryos run's `geo_glob_vilma.nc`) against PaleoMIST and the GAPSLIP 2.1 indicators, by groups and marine-bed areas; included by the calibration and plot tools |
| [plastic_score.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_score.sh) | Scores one kryos run's VILMA2 stream against GAPSLIP 2.1 with `pm_compare.jl`, beside B1 and I910_it0 of `PM_REF`: `plastic_score.sh <run> <label> [<score folder>]` (default `$KRYOS_WORK/score/<label>/out`) |
| [plastic_pd_regions.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_pd_regions.jl) | Present-day ice thickness of the plastic backend against the FesmData topography hubs, per region |
| [plastic_residual_corr.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_residual_corr.jl) | Correlation of the index-point misfit in distance and age, and the model error sigma_m, per group, region and sector |
| [plastic_earth_compare.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_earth_compare.jl) | Plastic runs on different Earths: GAPSLIP scores by group, region, marine-bed area and site, site maps, RSL curves, volume and ESL |

## Plots and reports

| File | What |
|---|---|
| [plastic_lgm_map.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_lgm_map.jl) | Global maps of a plastic run after Gowan et al. (2021, Fig. 2): `mode=age`, `margins`, `smb` |
| [plastic_mapkit.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_mapkit.jl) | The drawing module of the global maps (Equal Earth projection, VILMA2 and instance fields, coastlines, margins, scored sites) |
| [plastic_ensemble_sd_map.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_ensemble_sd_map.jl) | sd of RSL and ice thickness across ensemble members at one age, with that age's index points |
| [plastic_ensemble_misfit.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_ensemble_misfit.jl) | Misfit of a wave by data group against PaleoMIST's ice (PaleoMIST, B1 1-D, B3 3-D): implausibility, bias, sites by share of chi2; d against age and per site |
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

## Archive

[tools/levante/archive/](https://github.com/fesmc/kryos/tree/main/tools/levante/archive): the gate, regression and
Yelmo-calibration scripts of the development phases (`compare_*`, `gate_S*`, `g2*`-`g7*`, `calib*`, `ant_friction*`,
`setup_diag*`, `run_kryos*`, `run_vilma2.sh`, `recon_t0.jl`, `vilma_link_test.sh`, `build_stack.sh`, ...), kept as they
were run so that the reports in [Developer notes](developer-notes.md) stay traceable; not maintained.
