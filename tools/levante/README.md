# Levante tools of kryos

The scripts that set up, run, calibrate, score and plot the plastic ice / VILMA2 experiments on DKRZ Levante
(documentation: https://fesmc.github.io/kryos/, "Experiments" and "Reference: Tools").

## Once

```bash
tools/levante/setup.sh
```

checks the paths, the builds and Julia; links `fesmdata` and `kryos_data` into the checkout and copies the shared grid
files and remapping weights (`$KRYOS_DATA/maps`) into its `maps/`; checks runme (>= 0.5.27) and creates
`.runme/config.toml` from `.runme/config.default.toml`; and instantiates the Julia project of the tools (`Project.toml`,
`Manifest.toml`: CairoMakie, NCDatasets, DelimitedFiles on Julia 1.13) in your own depot. You need your own juliaup with
the 1.13 channel and runme (`pip install --user git+https://github.com/fesmc/runme`; `setup.sh` prints how).

## The environment

Every script sources `env.sh`: the kryos checkout (`KRYOS`, the one the script is in), your writable work folder
(`KRYOS_WORK`, default `/work/ba1442/$USER/kryos-work`), the shared stores (`FESMDATA`, `KRYOS_DATA`, read-only), the
VILMA2 checkout of the standalone runs (`VILMA2`), the reference standalone runs (`PM_REF`), the Slurm account of the
small jobs of the tools (`ACCOUNT`, default that of `.runme/config.toml`), the runme command (`KRYOS_RUNME`) and the
julia command (`KRYOS_JULIA`). Override any of them in the environment or in `~/.kryos.env`, e.g.

```bash
KRYOS_WORK=/work/ba1442/$USER/paleo
ACCOUNT=ba1442
```

Model runs (`kryos.x`, `vilma.x`, the Stage 1 fits) are staged and submitted with runme
(https://github.com/fesmc/runme), from the checkout (its `.runme/`), e.g.

```bash
runme -r -s -q shared -w 08:00:00 --omp 64 -m 60G -o $KRYOS_WORK/runs/sci-3d -n par/kryos_NH-ANT-PAT_plastic_vilma2.nml
```

The calibration and ensemble members are one runme job array per iteration (`members.sh`, from the calibration
folder), the standalone VILMA2 runs runme runs of the cases in `paleomist/` from a folder made a project by
`pm_project.sh`. Scripts that run as Slurm jobs submit themselves: run them directly
(`tools/levante/plastic_score.sh ...`), not with `sbatch`. Julia tools run with
`julia --project=tools/levante tools/levante/<tool>.jl ...` (or `kryos_jl` after sourcing `env.sh`, or as a job with
`tools/levante/jl.sh <out file> tools/levante/<tool>.jl ...`).

## What is here

| purpose | scripts |
|---------|---------|
| one run | runme from the checkout (`.runme/`); `set_nml.py` (set keys of a par file in place) |
| scoring | `plastic_score.sh` (a run against GAPSLIP 2.1 beside the standalone B1), `pm_compare.jl` |
| Gowan et al. (2021) with standalone VILMA2 | `pm_prep.jl` (inputs), `pm_project.sh` (a folder of runs as a runme project), `paleomist/` (the runme cases, one per Earth), `pm_iterref.sh` (the iterated reference), `pm_visc1d.sh`, `pm_selen.sh`, `pm_maps.jl` |
| reference correction of the 3-D Earth | `gia_refcorr.sh`, `gia_refcorr.jl`, `gia_refcorr_plot.jl` |
| SMB tuning on imposed ice (data backend) | `smb_basin_budget.jl`, `smb_itm_fit.jl` |
| calibration | `taub_fit/` (Stage 1 namelists, runme runs of `taub_fit`), `plastic_taub_fit_report.jl`, `plastic_taub_scale.jl`, `plastic_esmda.sh` / `.jl` / `*.toml`, `plastic_ensemble.sh` / `.jl` / `_v1.toml`, `members.sh` (the member job array) |
| analysis and figures | `plastic_lgm_map.jl`, `plastic_mapkit.jl`, `plastic_earth_compare.jl`, `plastic_residual_corr.jl`, `plastic_ensemble_sd_map.jl`, `plastic_ensemble_misfit.jl`, `plastic_ensemble_spread.jl`, `plastic_esmda_report.jl`, `plastic_esmda_diag.jl`, `plastic_pd_regions.jl`, `plastic_maps.jl`, `plastic_implied_pd.jl`, `plastic_closure_maps.jl` |

The released calibration products (Stage 1 factor files, the ES-MDA v2 posterior) are in `$KRYOS_DATA/plastic`
(README there). `archive/` holds the gate, regression and Yelmo-calibration scripts of the development phases, kept as
they were run (their paths are those of their time; not maintained).
