---
title: "Calibration and ensembles"
---

The plastic ice sheet (see [PlasticIceSheet](backend-plastic.md)) takes its thickness from the margin chronology and
the basal shear stress tau_b of Gowan et al.'s (2021) domains. Gowan tuned tau_b by hand. kryos instead calibrates it
in two steps and then samples it:

1. **Stage 1**: one constant factor per tau_b domain, fitted offline to the present-day ice thickness.
2. **ES-MDA**: time-dependent multipliers per *scale group* of domains, calibrated online (VILMA2 in the loop)
   against the GAPSLIP 2.1 sea-level index points (calibrations v1 and v2).
3. **Weighted ensemble**: a space-filling sample of the same multipliers. Each member gets a likelihood weight under
   a correlated error model per data group, and no member is discarded.

All of it runs on Levante: the runs, the scoring and the reports. The detailed record of each calibration, with its
checks and numbers, is
[docs/dev/plastic-backend.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-backend.md) ("Calibration
of tau_b", "Calibration v1", "Calibration v2"). The ensemble design is in
[docs/dev/plastic-ensemble-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-ensemble-design.md).

## How tau_b is composed

With `taub_method = "domains"`, the backend builds the tau_b of a domain at time t from three factors:

    tau_b(d, t) = tau_Gowan(d, t) * exp(u_d) * exp(theta_g(d)(t))

| Factor | Source | Key (`[ice_plastic]`) |
|---|---|---|
| `tau_Gowan(d, t)` | the FesmData product `{grid}_IceHistory-PaleoMIST_taub.nc`: Gowan's series per domain, linear in time | `taub_domain_file` |
| `exp(u_d)` | Stage 1: one log factor per domain, constant in time | `taub_domain_factor_file` |
| `exp(theta_g(t))` | the calibration or ensemble member: one log factor per scale group and time node, linear between the nodes and held outside them | `taub_scale_file` |

`taub_default` (the cells without a domain, 5000 Pa in the par files) is never scaled. With `taub_smooth_km`
(default 50 km), tau_b is then smoothed after these factors are applied. One file of each kind can serve every
instance: an id of another grid has no cells. A file of zeros gives Gowan's tau_b bit for bit (tested).

## Stage 1: present-day factors per domain

**Program:** [tests/plastic_taub_fit.f90](https://github.com/fesmc/kryos/blob/main/tests/plastic_taub_fit.f90),
offline and deterministic. It uses the backend's own modules and the library, so it needs a `plastic=1` build.

Stage 1 runs the backend at t = 0. The margin is the history's 0 ka slice (`margin_var > margin_min`, at the
backend's `margin_position`). The bed is that of the topography file, with z_sl = 0. tau_b is Gowan's at 0 ka times
exp(u_d) on the cells of domain d, followed by a cold solve. The objective is:

- **Fit cells:** the ice cells of the margin where the observation is not floating (`mask /= 3`).
- **Misfit:** chi2(u) = sum over the fit cells of ((H(u) - H_obs) / sigma)^2, with sigma = max(`sigma_min`,
  `sigma_rel` H_obs) (defaults: 50 m and 0.1).
- **Free domains:** those that have fit cells. All other domains keep u = 0.
- **Penalty:** lambda sum over adjacent free domains of w (u_d - u_d')^2, where w is the length [km] of their shared
  boundary.

The objective is minimised by Levenberg-Marquardt. The Jacobian comes from forward differences (`du`, one solve per
free domain, OpenMP over the domains), with |u| <= `u_max` (default 2.3: factors 0.1-10). Each lambda in the list is
solved in turn, starting from the solution of the previous one.

```bash
make libkryos/bin/plastic_taub_fit.x plastic=1
plastic_taub_fit.x <fit.nml>
```

The namelist `&fit` has these keys: `hist`, `margin_var`, `margin_min`, `margin_position`, `taub_file`,
`taub_default`, `topo_file`, `bed_var`, `H_var`, `mask_var`, `mode`, `lambda` (the list), `sigma_min`, `sigma_rel`,
`du`, `u_max`, `tol`, `n_iter`, `out`. The namelists of the released fits and their job are in
[tools/levante/taub_fit/](https://github.com/fesmc/kryos/tree/main/tools/levante/taub_fit) (`fit_NH.nml`, `fit_ANT.nml`,
`fit_PAT.nml`; lambda = 100, 30, 10, 3, 1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 0):

```bash
cd <fit folder>
K=<kryos checkout> sbatch $K/tools/levante/taub_fit/fit.sh $K/tools/levante/taub_fit/fit_NH.nml
```

The report picks lambda for each grid and writes the factor file that the backend reads:

```bash
julia --project=<julia-plot> tools/levante/plastic_taub_fit_report.jl <out folder> <fit_<grid>.nc> ... [<grid>=<lambda> ...]
```

The report chooses lambda by the discrepancy principle: the largest lambda with chi2/n <= 1. Where no lambda
reaches that, it takes the corner of the L-curve instead (the largest Menger curvature). It writes
`taub_domain_factor.nc` (`domain_id`, `log_factor` of the free domains of every grid), `report.md`, `lcurve.png` and
`maps_<grid>.png`.

The v2 result is recorded in plastic-backend.md, "Stage 1":

- NH-32KM: lambda 0.03 at the L-curve corner, chi2/n 12.54 -> 1.80.
- ANT-16KM: lambda 3 by the discrepancy principle, chi2/n 3.76 -> 0.94.
- Greenland: mean H - H_obs +132 -> -20 m.

The factors of small ice caps (Svalbard, the Canadian Arctic) go to 0.1-0.2. Part of that is a representation
error.

The ensemble config points at a later Stage 1 file that includes Patagonia,
`runs/taub-fit-v3/report/taub_domain_factor.nc`: the v2 fits of NH-32KM and ANT-16KM plus the fit of PAT-PAL-8KM
(`fit_PAT.nml`; section "Stage 1 v3" of [docs/dev/plastic-backend.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-backend.md)).

## Scale groups and nodes

**Tool:** [tools/levante/plastic_taub_scale.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_taub_scale.jl).
It is included by `plastic_esmda.jl` and reads the `[scale]` table of a calibration or ensemble config.

```bash
julia --project=tools tools/levante/plastic_taub_scale.jl <config.toml> <out.nc> [<log factor>]
```

This writes a multiplier file with that log factor for every group and node (default 0: Gowan's tau_b), and
`<out>.groups.csv`.

| `[scale]` key | Meaning |
|---|---|
| `taub_files`, `regions_files` | the tau_b products and REGIONS maps of the grids the calibrated instances run on (or `taub_file`, `regions_file`); cells are counted over all of them |
| `nodes_ka` | the time nodes [ka BP], any order |
| `pinned_ka` | nodes held at log factor 0 in every group; they are not parameters (the present day, fitted in Stage 1) |
| `[[scale.groups]]` | in order: `name`; `region` (Gowan's region of the domain: 1 North America including Greenland, 2 Eurasia, 3 Antarctica, 4 Patagonia); optionally `codes`, `codes_var` (default `"region_3"`) and `codes_in` (default true), or `greenland` (true/false, short for code 103 of `region_2`); `tie` (lists of node ages that share one parameter); `corr_time_kyr` (the group's own prior correlation time) |

A domain belongs to a group when more than half of its cells on the grids lie in the group's codes. With
`codes_in = false`, it belongs when not more than half do. The codes come from FesmData Regions v2: `region_2` 103
Greenland; `region_3` 10201 Russian Arctic, 10202 British Isles, 10203 Svalbard, 10204 Fennoscandia, 20101 West
Antarctica, 20102 East Antarctica, 20103 Antarctic Peninsula. A domain that matches no group is unscaled (group 0).
A domain that matches two groups stops the tool.

The file has the dimensions `domain`, `scale_group` and `node`, and the variables `domain_id(domain)`,
`group(domain)`, `time(node)` [years ago] and `log_factor(scale_group, node)`, with the attribute `group_names`.

The groups of the two configs in use (both with nodes 80, 60, 40, 25, 20, 15, 10, 5 and 0 ka, 0 ka pinned):

| Config | Group | Domains | Tied nodes |
|---|---|---|---|
| `plastic_esmda_v2.toml` | Laurentide | region 1, not Greenland | none |
| | Greenland | region 1, Greenland | 80=60 |
| | Eurasia | region 2 | 80=60=40 |
| | Antarctica | region 3 | 80=60=40=25=20 |
| `plastic_ensemble_v1.toml` | Laurentide | region 1, not Greenland; `corr_time_kyr = 10` | 25=20, 15=10 |
| | Greenland | region 1, Greenland | 80=60 |
| | Fennoscandia | region 2 without 10201, 10202, 10203 | 80=60=40 |
| | Barents-Kara | region 2, 10201, 10203 | 80=60=40 |
| | British-Irish | region 2, 10202 | 80=60=40=25 |
| | West Antarctica | region 3, 20101, 20103 | 80=60=40=25=20 |
| | East Antarctica | region 3, not 20101, 20103 | 80=60=40=25=20 |
| | Patagonia | region 4 (no RSL data: prior only) | 80=60=40=25=20 |

v2 has 25 parameters (plastic-backend.md). By its nodes and ties, ensemble v1 has 42. The `ensemble.md` written at
init lists them with their prior sd.

## Observations and the misfit

Both methods score a member the way
[tools/levante/pm_compare.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/pm_compare.jl) does:

- **Data:** the GAPSLIP 2.1 primary index points (0-80 ka) of pm_compare.jl's `GROUPS`: Laurentide, Greenland,
  Fennoscandia/Barents, Antarctica, far field.
- **Model RSL:** read from the run's VILMA2 stream `geo_glob_vilma.nc` (`read_run`, `series_run`). It is relative to
  0 ka, bilinear in space and linear in time between the 2.5 kyr records.
- **Misfit:** d = model - obs.
- **Indicator error:** sigma_i is the scorer's `ip_sigma`. It combines half the 2-sigma elevation range on the side
  of d with the age uncertainty times the model slope.
- **Limiting data:** not assimilated or weighted. Their violation fractions are reported as diagnostics.

The observation error is

$$R = \mathrm{diag}(\sigma_i^2) + \sigma_m^2 C,\qquad C_{ij} = \exp(-\mathrm{dist}_{ij}/L_x)\,\exp(-|t_i - t_j|/L_t)$$

with dist the great-circle distance and t the indicator age. sigma_m is the model error that the indicator's sigma
lacks: the GIA and Earth model, the resolution and the RSL interpolation.

## ES-MDA

The ensemble smoother with multiple data assimilation (Emerick and Reynolds 2013) runs an ensemble of parameter
sets, then moves every member by a Kalman-type update. The update is built from the ensemble covariance between the
parameters and the predicted data. "Multiple data assimilation" means the same data are assimilated N_a times, each
time with the observation error inflated by alpha_k (sum of 1/alpha_k = 1). This approximates an iterative solution
of a nonlinear problem.

**Tools:**

- [plastic_esmda.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_esmda.jl): init, scoring, update,
  and a self-test.
- [plastic_esmda.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_esmda.sh): the jobs.
- [plastic_esmda.toml](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_esmda.toml) (v1) and
  [plastic_esmda_v2.toml](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_esmda_v2.toml) (v2): the
  configs.

```bash
julia --project=<tools> plastic_esmda.jl init <config.toml> <calibration folder>   # prior ensemble = iteration 1
julia --project=<tools> plastic_esmda.jl update <calibration folder> <k>           # score iteration k, update, write k + 1
julia --project=<tools> plastic_esmda.jl selftest                                  # the update against a linear-Gaussian posterior
```

What one calibration does:

- **Parameters:** theta, the log multipliers per group and node (group-major), with tied nodes sharing one value and
  pinned nodes left out.
- **Prior:** N(0, C) within a group, C_pq = s(t_p) s(t_q) exp(-|t_p - t_q| / L). Here t_p is the mean age of a
  parameter's nodes, and s(t) = sigma sin^2(pi t / (2 `ramp_ka`)) below `ramp_ka`, sigma above. The groups are
  independent.
- **Clipping:** after an update, every member is clipped to |theta| <= `theta_max`.
- **Fixed error:** sigma_i is evaluated at the ensemble-mean RSL of the iteration and held fixed within it.
- **Gross-error cut:** fixed in iteration 1 at the prior mean member (m000). Index points with
  |d| / sqrt(sigma_i^2 + sigma_m^2) > `z_cut` are not assimilated (`obs_kept.csv`). Every score still covers all the
  points.
- **Update:** R = L L' (Cholesky, once per iteration). The predictions and observations are whitened by L^-1, the
  observations are perturbed with N(0, alpha I) in the whitened space, and the update is done in the ensemble
  subspace (an Ne x Ne solve). `corr_dist_km = 0` gives v1's independent `sigma_add`.
- **Mean member:** member 0 is the prior mean in iteration 1 and the ensemble mean after each update. After the last
  update, `final = "mean"` runs the posterior mean.
- **Failures:** a member that fails is dropped. If fewer than `min_members` succeed, the calibration stops.

| `[esmda]` / `[prior]` | v1 (`plastic_esmda.toml`) | v2 (`plastic_esmda_v2.toml`) |
|---|---|---|
| Instances | NH-32KM (`par/kryos_North_plastic_vilma2_b1.nml`) | NH-32KM, ANT-16KM (`par/kryos_NH-ANT_plastic_vilma2_b1.nml`) |
| Stage 1 factors | none | `runs/taub-fit/taub_domain_factor.nc`, frozen in the folder |
| Parameters | 24 (3 groups) | 25 (4 groups, 0 ka pinned) |
| Prior | sigma 0.3, correlation 10 kyr | sigma 0.3 ramped from 0 at 0 ka (`ramp_ka = 10`), correlation 20 kyr |
| Ne, alpha | 40; 9.333, 7, 4, 2 | same |
| Model error | `sigma_add` 5 m, independent | `sigma_m` 5 m, `corr_dist_km` 150, `corr_time_kyr` 2 |
| Gross-error cut | `z_cut` 5 | same |
| Members | 8 threads, 10 GB, throttle 29 | same |

The results are in plastic-backend.md. v1 thinned today's Greenland by about a quarter, and its factors alternate
between neighbouring deglacial nodes. v2 holds the present day. Its RMS z over all index points goes from 26.7 at the
prior mean to 21.6 at the posterior mean. But Antarctica's 80-20 ka node goes to a factor of 1.63 (+24 m SLE at the
LGM) while its own index points fit worse: "the posterior is usable for the NH groups, not for Antarctica or the
global ESL".

### The calibration folder

| Path | Content |
|---|---|
| `config.toml`, `kryos.x`, `input/`, `par.nml` | frozen copies: the config, the binary and input of `[run] kryos_dir`, and the par file with `[run] set` applied (and the Stage 1 file, `taub_domain_factor.nc`) |
| `slurm.env` | `THREADS`, `MEM`, `TIME`, `PART`, `ACCT`, `THROTTLE`, `MAPS` from `[run]`; read again at every submission |
| `groups.csv`, `obs.csv`, `obs_kept.csv` | the scale groups; the data per group; the index points kept by the cut |
| `esmda.md` | the report: per iteration, the mean member's scores and the ensemble's, chi2, ESL at 20 ka, theta before and after |
| `next` | the iteration to run, or `done` |
| `jobs.txt` | one line per submission: job ids of the members and the analysis |
| `iter_kk/` | `theta.csv`, `members.txt`, `scores.csv` (per member and set), `members.csv` (ok, chi2, ESL at 20 ka, run time), `slurm-<m>.out`, `update.out` |
| `iter_kk/mNNN/` | `taub_scale.nc`, `par.nml`, the run, `status` (exit code and seconds) |

Every member writes the product: the 2-D ice (`H_ice`, `z_srf`, `z_bed`, `z_sl`, `f_ice`, `M_ice`) of each scaled
instance on its `out_2D` schedule (2.5 kyr), its 0-D series, and the VILMA2 stream of the par file.

## Weighted ensemble

**Tools:**

- [plastic_ensemble.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_ensemble.jl): the sample and
  the weights.
- [plastic_ensemble.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_ensemble.sh): the jobs.
- [plastic_ensemble_v1.toml](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_ensemble_v1.toml): the
  config.

The design is [docs/dev/plastic-ensemble-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-ensemble-design.md).

Why weights instead of ES-MDA: ES-MDA assumes Gaussian errors and returns one compromise; its spread collapsed from
sd 0.27 to 0.11 in v2. The ensemble instead gives a range of ice histories (80-0 kyr) and their RSL, each with its
plausibility. The aim is that low, mid and high ice volumes and regional trade-offs all come out.

**Basis:**

- **Par file:** `par/kryos_NH-ANT-PAT_plastic_vilma2.nml`, with three plastic instances: NH-32KM, ANT-16KM and
  PAT-PAL-8KM.
- **Earth:** the VILMA2 3-D science setting (Bagge et al. 2021 viscosity, lmax 170, refcorr reference).
- **tau_b:** Gowan's domains times the Stage 1 factors (with Patagonia).
- **Members:** the same folder layout as a calibration, with each member writing the product.

**Sample (wave 1).** A Latin hypercube in the standard normal space of the parameters: n strata per parameter, one
point in each, in random order. It is correlated by the Cholesky factor of the prior covariance and clipped to
|theta| <= `theta_max`. Member m000 is the prior mean. `[ensemble] n = 24` members besides m000, `seed = 1`.

**Likelihood.** For each data group g, the likelihood uses the group's index points (all of them: no gross-error
cut), d = model - obs, and the error

$$R_g = \mathrm{diag}(\sigma_i^2) + \sigma_{m,g}^2 C_g,\qquad \log L = \sum_g -\tfrac12\left(d^\top R_g^{-1} d + \log\det R_g\right)$$

Here sigma_i is the member's own `ip_sigma` (its slope and the side of d). The groups are independent.

**Weights.** w_m ∝ L_m: the sample follows the prior, so the prior and the sampling density cancel. The weights are
normalised over the members that ran. The effective sample size is ESS = 1 / sum w^2, computed overall and per data
group (from that group's likelihood alone). No member is discarded: a member with a tiny weight stays in the set.

**Output:** `iter_kk/weights.csv` (per member: ok, log L per group and total, chi2/n per group, ESL at 20 ka, weight)
and `weights.md` (the ESS, the members by weight, the weighted parameter means and sd).

The code differs from the design in four ways:

- **Sampling:** the design uses Sobol points and 100 members per wave. v1 is a Latin hypercube with 24 members,
  marked "proof of concept".
- **Refinement waves:** `plastic_ensemble.jl` writes wave 1 only (`init`). The proposal fitted to the weighted members,
  and the importance weights over several waves, are not implemented. `plastic_ensemble.sh submit <folder> <k>`
  submits the members of a wave folder that already exists.
- **Products:** the products of the design (weighted quantiles, scenario members, a held-out check) and
  `plastic_ensemble_report.jl` are not in the repository.
- **Earth:** the Earth as an extra dimension of the ensemble is planned for later.

### Model-error parameters per data group

The ensemble's sigma_m and correlation lengths come from
[tools/levante/plastic_residual_corr.jl](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_residual_corr.jl),
applied to the v2 posterior on the 1-D Earth. That tool computes three things for each group:

- the non-centred correlation of d over pairs of index points, against distance (pairs with |dt| <= 1 kyr) and
  against age difference (pairs within 25 km);
- an exponential c exp(-h / L) fitted to the bins;
- the error beyond the observations, sigma_m^2 = mean(d^2) - mean(sigma_i^2).

| `[likelihood.<group>]` | `sigma_m` [m] | `corr_dist_km` | `corr_time_kyr` |
|---|---|---|---|
| Laurentide | 12 | 350 | 10 |
| Greenland | 30 | 500 | 5 |
| Fennoscandia/Barents | 10 | 250 | 5 |
| Antarctica | 8 | 500 | 10 |
| far field | 13 | 50 | 5 |
| (v2, all groups) | 5 | 150 | 2 |

The design note says these values are to be re-estimated on the first 3-D wave.

## Levante job limits and cost

- **At most 30 jobs** may be queued or running at a time. The member arrays are throttled (`--array=<ids>%THROTTLE`):
  - ES-MDA uses `throttle = 29`: 29 members plus the queued analysis.
  - Ensemble v1 uses `throttle = 25`: 25 members plus the weights job plus one other run.
  - `THROTTLE` is read from the folder's `slurm.env` at every submission. Lowering it there before the next iteration
    or wave is submitted (when other jobs run) takes effect from that submission.
- **Billing:** the shared partition bills max(cpus, mem / 0.9 GB). An ES-MDA member at 8 threads and 10 GB is billed
  12 cpus (about 2.5 cpu-h). VILMA2 barely scales: 12.5 min on 8 threads, 8 min on 64.
- **CPU cap:** the shared partition also caps the CPUs one user can hold at once (`QOSMaxCpuPerUserLimit` in `squeue`):
  in November 2026 about 512 CPUs, i.e. 8 members of 64 threads running at a time; the rest of the throttled array
  waits. This, not the throttle, bounds a 3-D wave (25 x 64 cpus: about 14 h for a wave of 25). Short single jobs fit
  better on the compute partition (whole nodes, 128 cpus).
- **Recorded cost:** v1 took 165 runs, 3.0 h of wall time and 464 cpu-h; v2 took 522 cpu-h. A 3-D coupled run takes
  3.5-4.2 h on 64 cores (`threads = 64`, `mem = "30G"`, `time = "08:00:00"` in ensemble v1). The design estimates
  about 24,000 core-h for a wave of 100 members.

## Start, monitor, score

Both job scripts:

- **Tools and Julia setup:** they find the tools next to themselves (or in `ESMDA_TOOLS`) and run
  `julia +1.13 --project=/work/ba1442/robinson/kryos-gate/p4w7/kryos/tools` with
  `JULIA_DEPOT_PATH=/work/ba1442/robinson/kryos-gate/julia-depot`.
- **Binary:** they need a `kryos.x` built with `plastic=1 vilma=1` in `[run] kryos_dir`
  (`libkryos/bin/kryos.x`), plus the maps folder of `[run] maps`.
- **Folder:** the target folder must not exist yet.

ES-MDA calibration:

```bash
bash tools/levante/plastic_esmda.sh start tools/levante/plastic_esmda_v2.toml /work/ba1442/robinson/kryos-ice-plastic/runs/esmda-<name>
```

`start` submits the init as a job. The init writes the prior ensemble and submits iteration 1: the members as an
array, then the analysis (`afterany`). The analysis scores, updates and submits the next iteration until `next` says
`done`. To resubmit an iteration by hand:

```bash
bash tools/levante/plastic_esmda.sh submit <calibration folder> <k>
```

Weighted ensemble:

```bash
tools/levante/plastic_ensemble.sh start tools/levante/plastic_ensemble_v1.toml <ensemble folder>
tools/levante/plastic_ensemble.sh submit <ensemble folder> <k>     # the members of wave k and the weights after them
```

The members run as calibration members (`plastic_esmda.sh member`). The weights job runs after the members
(`afterany`).

To monitor:

- `jobs.txt` holds the job ids of each submission.
- `init.out` is the init log.
- `iter_kk/slurm-<m>.out` is a member's job log. `iter_kk/mNNN/status` holds its exit code and run time, and
  `out.out` and `out.err` hold kryos's output.
- `iter_kk/update.out` (ES-MDA) and `iter_kk/weights.out` (ensemble) are the analysis logs.
- `next` and `esmda.md` give the state of a calibration, and `iter_kk/weights.md` that of an ensemble wave.

The SLURM queue shows the running jobs (`squeue -u $USER`).

Scoring is part of the jobs: an ES-MDA iteration is scored by its analysis (`plastic_esmda.jl update`, which also
writes iteration k + 1), and an ensemble wave by its weights job. To recompute the weights of a wave by hand, for
example after a member was rerun:

```bash
julia --project=<tools> plastic_ensemble.jl weights <ensemble folder> [<wave>]
```

The reports and figures of a calibration or ensemble are covered in [Analysis and figures](experiment-analysis.md):
`plastic_esmda_report.jl`, `plastic_esmda_diag.jl`, `plastic_ensemble_spread.jl`, `plastic_ensemble_sd_map.jl` and
`plastic_residual_corr.jl`.
