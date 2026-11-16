---
title: "Ensembles"
---

How to run, monitor and score a weighted ensemble of plastic ice histories, and how to make the present-day factors it
starts from. The method (the prior, the likelihood, the weights, the constraints on the ice) and what the ensembles found:
[Our ice: the τ_b ensembles](../science/ice-ensembles.md). The data groups and the scores: [Data and scores](../science/data-scores.md).
How τ_b is composed from Gowan's values, the present-day factors and the multipliers:
[PlasticIceSheet, basal shear stress](../backend-plastic.md#basal-shear-stress). Status as of 2026-11-16.

## The ensembles

| | ens-v1 | ens-v2 |
|---|---|---|
| config | `tools/levante/plastic_ensemble_v1.toml` | `tools/levante/plastic_ensemble_v2.toml` |
| launch | `experiments/ens.sh v1` | `experiments/ens.sh v2` |
| par file | `par/kryos_NH-ANT-PAT_plastic_vilma2.nml` (NH-32KM, ANT-16KM, PAT-PAL-8KM; Bagge 6) | the same |
| members | 25 (m000 the prior mean) | 50 |
| parameters | τ_b multipliers per scale group and node | the same, plus shifts of log10 η of the 3-D Earth (`[earth]`: upper and lower mantle, and per region in the upper mantle) |
| data that judge | the RSL groups (Fennoscandia/Barents as one) | the RSL groups (Fennoscandia and Barents split), GNSS Fennoscandia and Hudson Bay, ESL |
| error model | `plastic_residual_corr.jl` on the ES-MDA v2 posterior (1-D Earth) | refitted from ens-v1 wave 1 (`plastic_ensemble.jl errmodel`) |
| where it ran | Levante (copied to albedo) | albedo |

```bash
experiments/ens.sh v2 [<folder>]        # default $KRYOS_WORK/ens-v2; the folder must not exist
```

One init job draws the members (seeded) and submits them as a runme job array, then the weights job after them. An ens-v2
member takes about 3 h on 64 threads on albedo. The scores of every member and every reference: `experiments/scores.sh <folder>`.

## Start, monitor, score

`experiments/ens.sh` calls `tools/levante/plastic_ensemble.sh`, which can be used directly:

```bash
tools/levante/plastic_ensemble.sh start <config.toml> <ensemble folder>   # init (wave 1) as a job, then its members and the weights
tools/levante/plastic_ensemble.sh submit <ensemble folder> <k>            # the members of wave k and the weights after them
```

- **Environment:** the scripts source `tools/levante/env.sh` and run the Julia tools with `KRYOS_JULIA --project=tools/levante`
  ([Getting started](../getting-started.md#the-hpc-tools)).
- **Config placeholders:** `{kryos}`, `{fesmdata}`, `{kryos_data}`, `{maps}`, `{work}` (`KRYOS`, `FESMDATA`, `KRYOS_DATA`,
  the checkout's `maps/`, `KRYOS_WORK`). `[run] kryos_dir`, `set_nml`, `maps`, `partition` (`KRYOS_QUEUE`) and
  `[score] fesmdata` have defaults. `[run] partition` is a runme queue alias of the machine (Levante `shared`, `compute`;
  albedo `12h`, `48h`). The init writes the expanded config into the folder, which then does not depend on the environment.
- **Binary:** a `kryos.x` built with `plastic=1 vilma=1` in `[run] kryos_dir`, and runme >= 0.5.27 (`tools/levante/setup.sh`).

The members of a wave are staged and submitted by runme from the ensemble folder (`tools/levante/members.sh`) as one job
array, throttled by `THROTTLE` of the folder's `slurm.env` (read again at every submission). The weights job runs after the
array (`afterany`).

To monitor:

- `jobs.txt`: the job ids of each submission; `init.out`: the init log.
- `iter_01/runme.out`: runme's log of the staging and submission.
- `iter_01/<m>/status`: a member's exit code and run time; `out.out`, `out.err`: its job and kryos output.
- `iter_01/weights.out`, `iter_01/weights.md`: the weights job and its report.
- `squeue -u $USER`.

Scores and weights by hand, for example after a member was rerun or the error model changed:

```bash
experiments/scores.sh <ensemble folder> [<wave>]                                     # a job: iter_<wave>/scores.csv
julia --project=tools/levante tools/levante/plastic_ensemble.jl weights <ensemble folder> [<wave>]
```

`scores.csv` has one row per series (member or reference), data group and window: `n`, RMSE, RMS z, χ²/n, I and log L. The
references are the `[[references]]` of the folder's config. A reference whose run is not finished is left out and named.

## The error model

The model error per data group (`[likelihood.<group>]`: `sigma_m`, `corr_dist_km`, `corr_time_kyr`) is refitted from the
residuals of a wave:

```bash
julia --project=tools/levante tools/levante/plastic_ensemble.jl errmodel <ensemble folder> [<wave>] [<top>] [refs=<regex>] [groups=<g1;g2>]
```

It writes `errmodel.toml` and `errmodel.md` in the folder; `weights` and `scores` read it from then on in place of the config's values. `top` fits on the best members;
`refs` fits the point groups on the references whose name matches; `groups` refits only those groups (the Earth search refits
the GNSS groups on its own cases this way). ens-v2's error model is ens-v1 wave 1's (top 6 members), with Fennoscandia and Barents
refitted after their split on 2026-11-12.

## The ensemble folder

| Path | Content |
|---|---|
| `config.toml`, `config.source.toml` | the config expanded with defaults; the config as given |
| `kryos.x`, `input/`, `par.nml` | frozen copies: the binary and `input/` of `[run] kryos_dir`, and the par file with `[run] set` applied (with `taub_domain_factor.nc`) |
| `.runme/`, `fesmdata`, `kryos_data`, `maps` | the folder is a frozen runme project, and the links runme puts into every member |
| `slurm.env` | `THREADS`, `MEM`, `TIME`, `PART`, `THROTTLE` from `[run]` |
| `groups.csv`, `obs.csv` | the scale groups; the data per group |
| `errmodel.toml` | the refitted error model, if any |
| `iter_01/` | `theta.csv` (the parameters of every member), `members.txt`, `par.nml`, `scale/mNNN.nc` (the τ_b multipliers of each member), `scale_files.txt` (runme's `-i` table), `array.json`, `weights.csv`, `weights.md`, `scores.csv` |
| `iter_01/<m>/` | member m, staged by runme: `par.nml` with its scale file, the run, `out.out`, `out.err`, `status` |

Every member writes the product: the 2-D ice (`H_ice`, `z_srf`, `z_bed`, `z_sl`, `f_ice`, `M_ice`) of each scaled instance
every 2.5 kyr, its 0-D series, and the VILMA2 stream. Folders made before 2026-11-11 have the members in `iter_kk/mNNN/`; the
tools still read them.

## Scale groups and the multiplier file

`tools/levante/plastic_taub_scale.jl` (included by `plastic_ensemble.jl`) reads the `[scale]` table of a config and writes a
multiplier file:

```bash
julia --project=tools/levante tools/levante/plastic_taub_scale.jl <config.toml> <out.nc> [<log factor>]
```

| `[scale]` key | Meaning |
|---|---|
| `taub_files`, `regions_files` | the τ_b products and REGIONS maps of the grids the scaled instances run on |
| `nodes_ka` | the time nodes [ka BP] |
| `pinned_ka` | nodes held at log factor 0 in every group (the present day, fitted in Stage 1) |
| `[[scale.groups]]` | `name`; `region` (Gowan's region: 1 North America with Greenland, 2 Eurasia, 3 Antarctica, 4 Patagonia); optionally `codes`, `codes_var`, `codes_in` or `greenland`; `tie` (node ages that share one parameter); `corr_time_kyr` |

A domain belongs to a group when more than half of its cells lie in the group's codes (FesmData Regions v2). A domain that
matches no group is unscaled; one that matches two stops the tool. The file has `domain_id(domain)`, `group(domain)`,
`time(node)` [years ago] and `log_factor(scale_group, node)`, with the attribute `group_names`. The groups of v1 and v2:
Laurentide, Greenland, Fennoscandia, Barents-Kara, British-Irish, West Antarctica, East Antarctica, Patagonia (no RSL data:
prior only); nodes 80, 60, 40, 25, 20, 15, 10, 5 and 0 ka, 0 ka pinned. The `ensemble.md` written at init lists the parameters
with their prior sd.

## Stage 1: the present-day factors

The present-day factors, one per τ_b domain, are fitted offline to today's ice thickness by
[tests/plastic_taub_fit.f90](https://github.com/fesmc/kryos/blob/main/tests/plastic_taub_fit.f90) (the objective:
[PlasticIceSheet](../backend-plastic.md#factors-and-multipliers)). It needs a `plastic=1` build:

```bash
make libkryos/bin/plastic_taub_fit.x plastic=1
cd $KRYOS
runme -r -s -e taub_fit -n tools/levante/taub_fit/fit_NH.nml  -o $KRYOS_WORK/taub-fit/NH  -q $KRYOS_QUEUE -w 02:00:00 --omp 16 -m 24G
runme -r -s -e taub_fit -n tools/levante/taub_fit/fit_ANT.nml -o $KRYOS_WORK/taub-fit/ANT -q $KRYOS_QUEUE -w 02:00:00 --omp 16 -m 24G
runme -r -s -e taub_fit -n tools/levante/taub_fit/fit_PAT.nml -o $KRYOS_WORK/taub-fit/PAT -q $KRYOS_QUEUE -w 02:00:00 --omp 16 -m 24G
```

The namelists (`&fit`: `hist`, `margin_var`, `taub_file`, `topo_file`, `lambda` (a list), `sigma_min`, `sigma_rel`, `du`,
`u_max`, `out`, ...) are in [tools/levante/taub_fit/](https://github.com/fesmc/kryos/tree/main/tools/levante/taub_fit). Each
run writes its `out` file (`fit_NH-32KM.nc`, `fit_ANT-16KM.nc`, `fit_PAT-PAL-8KM.nc`). The report picks λ per grid (the
discrepancy principle, else the corner of the L-curve) and writes the factor file the backend reads:

```bash
julia --project=tools/levante tools/levante/plastic_taub_fit_report.jl <out folder> $KRYOS_WORK/taub-fit/NH/fit_NH-32KM.nc \
      $KRYOS_WORK/taub-fit/ANT/fit_ANT-16KM.nc $KRYOS_WORK/taub-fit/PAT/fit_PAT-PAL-8KM.nc [<grid>=<lambda> ...]
```

The released factor files are in `$KRYOS_DATA/plastic`: `taub_domain_factor_v2.nc` (NH-32KM, ANT-16KM) and
`taub_domain_factor_v3.nc` (with PAT-PAL-8KM; the one the par files and both ensembles use).

## Job limits and cost

| | albedo | Levante |
|---|---|---|
| jobs at a time | at most 50 of your own (queued and running) | 30 |
| member queue | `12h` (`[run] time = "11:00:00"`) | `shared` |
| a 3-D member | about 3 h on 64 threads | 3.0-4.3 h on 64 threads |

The member arrays are throttled (runme `--array THROTTLE`; v1 and v2: 25 members plus the weights job plus one other run).
Lower `THROTTLE` in the folder's `slurm.env` before the next submission when other jobs run. On Levante the shared partition
also caps the CPUs one user can hold (`QOSMaxCpuPerUserLimit`, about 512 in November 2026: 8 members of 64 threads at a time).
A wave of 50 3-D members costs about 12,000 core-h.

## ES-MDA (retired)

Calibrations v1 and v2 (2026-11, Levante) moved the τ_b multipliers by an ensemble smoother with multiple data assimilation
(`tools/levante/plastic_esmda.jl`, `plastic_esmda.sh`, configs `plastic_esmda.toml`, `plastic_esmda_v2.toml`). Its posterior
collapsed (sd 0.27 to 0.11) and moved Antarctica to an implausible LGM volume; the weighted ensemble replaced it. The tools stay
in the repository: the ensemble shares their folder layout, members and scale groups. The record of both calibrations is in
[docs/dev/plastic-backend.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-backend.md) ("Calibration v1",
"Calibration v2"); their posterior multiplier file is `$KRYOS_DATA/plastic/taub_scale_esmda_v2_posterior.nc`.

The reports and figures of an ensemble (`plastic_ensemble_spread.jl`, `plastic_ensemble_sd_map.jl`, `plastic_ensemble_misfit.jl`,
`plastic_residual_corr.jl`): [Analysis and figures](analysis.md#ensembles).
