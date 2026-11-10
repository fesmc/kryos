---
title: "Simulating our own ice"
---

The ice sheets of PaleoMIST rebuilt online, 80 ka to present: a perfectly-plastic ice sheet
(backend `plastic`, see [PlasticIceSheet](backend-plastic.md)) on PaleoMIST's margins and
Gowan's basal shear stress (tau_b) domains, solved at every 100 yr step on the bed and sea
level of the global VILMA2 Earth (see [VILMA2](backend-vilma2.md)) that its load deforms. The
ice is imposed, not integrated (no surface, basal or calving flux); what the run changes is
tau_b and the Earth. With Gowan's tau_b and his 1-D Earth this is the bridge of
[Reproducing Gowan et al. (2021)](experiment-gowan.md); with the present-day factors, the
calibrated multipliers and the 3-D Earth it is our own ice history.

- **Program:** `kryos.x` built with `plastic=1 vilma=1` (`yelmo=0` is enough; `openmp=1`),
  linked against the PlasticIceSheet Fortran library (`libplastic.a`, a sibling checkout
  `PlasticIceSheet`, built with the same `openmp=`; `config/common.mk`). See
  [Getting started](getting-started.md).
- **Config:** the par files below (see also [Par files](reference-par-files.md)).
- **Scripts:** no single-run script; the run folder and job follow the members of
  [plastic_esmda.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_esmda.sh)
  (Running, below). Scoring and comparison: `tools/levante/plastic_score.sh`, `pm_compare.jl`,
  `plastic_earth_compare.jl`.
- **Design notes:** [docs/dev/plastic-backend.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-backend.md).

| Par file | Domains | Earth | Reference |
|---|---|---|---|
| [kryos_North_plastic.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_North_plastic.nml) | North | none (today's bed, `z_sl = 0`): smoke test | |
| [kryos_North_plastic_vilma2_b1.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_North_plastic_vilma2_b1.nml) | North | 1-D Gowan (B1) | relaxed, no correction |
| [kryos_North_plastic_vilma2.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_North_plastic_vilma2.nml) | North | 3-D (Bagge 2021) | today's + r |
| [kryos_NH-ANT_plastic_vilma2_b1.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT_plastic_vilma2_b1.nml) | North, Antarctica | 1-D Gowan (B1) | relaxed, no correction |
| [kryos_NH-ANT_plastic_vilma2.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT_plastic_vilma2.nml) | North, Antarctica | 3-D (Bagge 2021) | today's + r |
| [kryos_NH-ANT-PAT_plastic_vilma2.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT-PAT_plastic_vilma2.nml) | North, Antarctica, PatagoniaPaleo | 3-D (Bagge 2021) | today's + r |

`kryos_NH-ANT_plastic_vilma2_b1.nml` is the par file of the calibration v2 members,
`kryos_NH-ANT-PAT_plastic_vilma2.nml` that of the weighted ensemble v1
([Calibration and ensembles](experiment-calibration.md)).
`kryos_NH-ANT-PAT_plastic_vilma2_inv.nml` adds the mass-balance inversion:
[Implied mass balance](experiment-implied-mb.md). The 3-D North-only file has not been run with
`mode = "surface"` (its earlier `"flat"` run was cancelled in VILMA2's spin-up).

## Domains and instances

| Domain | Hub and ice grid | Ice | Geo port | Topography (`[topo_*] path`) |
|---|---|---|---|---|
| `North` | NH-32KM | `ice_nor` | `geo_nor` | `TOPO-GEBCO2026` (1-D) or `TOPO-GEBCO2026-refcorr` (3-D) |
| `Antarctica` | ANT-16KM | `ice_ant` | `geo_ant` | `TOPO-BedMachine-v4` or `TOPO-BedMachine-v4-refcorr` |
| `PatagoniaPaleo` | PAT-PAL-8KM | `ice_pat` | `geo_pat` | `TOPO-GEBCO2026-refcorr` (IceBoost v2 ice) |

- Each ice instance is `backend = "plastic"` on its hub grid (`grid = ""`), every step
  (`dt = 0`), `drive = "taub"`, `mode = "surface"` (`"full"` does not converge on 32 km beds),
  densities 910 and 1028 (the run's `phys_const = "Earth"`). A domain has no climate, surface
  or shelf: the plastic ice needs none.
- Each geo port is `backend = "vilma2"` on the ice grid, `dt = 100` (it must equal
  `[geo_glob] dt`, else the run stops).
- The global Earth `geo_glob` (`backend = "vilma2"`, Gauss grid of `lmax` 170, `dt = 100`,
  `ice_bg = "globe_ice"`). The ice outside the footprints of the domains is the PaleoMIST
  history: `[globe_ice]` on GLOBAL-0.5DEG, `ext = "file"`, the background X (`[X_z_bed]`,
  `[X_H_ice]`) plus the anomaly of `[X_dH_ice]` (`fesmdata/Global/GLOBAL-0.5DEG/GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc`,
  `H_ice`, relative to 0 ka). See [Global hub](global-hub.qmd).

## Ice inputs

The margin and tau_b come from FesmData IceHistory (v2.3.0; PAT-PAL-8KM v2.4.0) on each grid:

- Margin: `margin_file = "fesmdata/{domain}/{grid_name}/{grid_name}_IceHistory-PaleoMIST.nc"`
  (33 slices, 80 ka to present, "years ago"), ice where `f_ice` > 0.5 (the ice fraction of
  Gowan's unsmoothed 5 km grids: his margin), interpolated in time by the signed distance
  (`margin_interp = "sdist"`), placed within the cells (`margin_position` default `"sdist"`);
  outside the slices' range the run stops (`margin_extrapolate = "error"`).
- tau_b: `taub_method = "domains"`, `taub_domain_file =
  "fesmdata/{domain}/{grid_name}/{grid_name}_IceHistory-PaleoMIST_taub.nc"`: the domain map
  `taub_domain`, the ids `domain_id`, their values `taub` and the series `taub_t` (linear in time
  between the records; tau_b is time dependent in PaleoMIST 1.0). Cells without a domain take
  `taub_default = 5000` Pa. `taub_smooth_km` (default 50) smooths tau_b over the ice.

tau_b of a domain d at time t:

    tau_b(d, t) = taub_t(d, t) * exp(log_factor_d) * exp(log_factor_g(d)(t))

| Key (per instance) | File | Content |
|---|---|---|
| `taub_domain_file` | FesmData `..._IceHistory-PaleoMIST_taub.nc` | Gowan's tau_b (unscaled) |
| `taub_domain_factor_file` | Stage 1: `/work/ba1442/robinson/kryos-ice-plastic/runs/taub-fit-v3/report/taub_domain_factor.nc` (NH-32KM, ANT-16KM, PAT-PAL-8KM); v2: `.../runs/taub-fit/taub_domain_factor.nc` (NH-32KM, ANT-16KM) | `domain_id(domain)`, `log_factor(domain)`: one factor per domain, constant in time, from the present-day thickness |
| `taub_scale_file` | `taub_scale.nc` of a calibration member, `<calibration>/iter_kk/mNNN/taub_scale.nc` | `group(domain)`, `time(node)`, `log_factor(scale_group, node)`: one multiplier per scale group, linear in time between the nodes |

- All par files ship with `taub_domain_factor_file = ""` and `taub_scale_file = ""`: Gowan's
  tau_b, unscaled. `""` is no file; a file of zeros gives tau_b bit for bit.
- One file of each kind holds the domains of all grids: set the same path on every instance
  (an id of another grid has no cells).
- Stage 1 (`tests/plastic_taub_fit.f90`, report and file by
  `tools/levante/plastic_taub_fit_report.jl`) fits the factors to today's thickness of the
  topography hub; `taub-fit-v3` is the fit with PAT-PAL-8KM, as used by
  `tools/levante/plastic_ensemble_v1.toml` ("Stage 1 with PAT"; section "Stage 1 v3" of
  [docs/dev/plastic-backend.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-backend.md)).
- The scale file of a member is written by `plastic_esmda.jl` (or `plastic_ensemble.jl`) into
  its folder. The posterior mean of an ES-MDA calibration is member `m000` of the iteration
  after the last update (e.g. `runs/esmda-v2/iter_05/m000/taub_scale.nc`).
  <!-- TODO(check): the posterior-mean folder of runs/esmda-v2 (iter_05/m000, from plastic_esmda.jl write_iteration) -->

The 3-D North-only file (`kryos_North_plastic_vilma2.nml`) does not list the two keys; they take
the default `""` of `input/kryos_defaults.nml`.

## 1-D or 3-D Earth

The `_b1` files and their 3-D counterparts differ only in the `[vilma]` group, the background X
and the topography paths (`diff par/kryos_NH-ANT_plastic_vilma2_b1.nml
par/kryos_NH-ANT_plastic_vilma2.nml`; the North pair likewise):

| | `_b1`: Gowan's 1-D Earth (B1) | 3-D: the science Earth |
|---|---|---|
| `[vilma]` layers | `n_layer = 12`: 0-80-120 km elastic, 4e20 to 670 km, 4e22 to the CMB | `n_layer = 11`: 0-80 km elastic, 4e20 to 670 km, 1e22 below (the radial profile of the reference) |
| lateral viscosity | `l_visc_3d = .false.` | `l_visc_3d = .true.`, `visc_3d_file = "fesmdata/Global/GLOBAL-T170/GLOBAL-T170_SolidEarth-Bagge2021.nc"`, `visc_3d_structure = 6` |
| spin-up | `pre_spinup_1d = .false.`, `equil_time_max = 2.0e6` | `pre_spinup_1d = .true.`, `equil_time_max = 2.0e4` (after the 1-D pre-spin-up) |
| tolerance | `equil_rate_tol = 1.0e-5` | same |
| rotation | `rotation = .true.` (set) | `.true.` (the default of `input/vilma_defaults.nml`) |
| `[X_z_bed]`, `[X_H_ice]` | PaleoMIST's 0 ka bed and ice, `fesmdata/Global/{grid_name}/{grid_name}_IceHistory-PaleoMIST.nc`, `mode = "transient"` (read at `ext_tref` = 0) | `kryos_data/refcorr/PaleoMIST-sci/{grid_name}/{grid_name}_REFCORR.nc` (PaleoMIST's 0 ka bed + r, its 0 ka ice), `mode = "static"` |
| `[topo_*] path` | `fesmdata/.../{grid_name}_TOPO-GEBCO2026.nc`, `..._TOPO-BedMachine-v4.nc` (today's, no correction) | `kryos_data/refcorr/PaleoMIST-sci/{domain}/{grid_name}/{grid_name}_TOPO-...-refcorr.nc` (today's + r) |

Both use `lmax = 170`, `earth = "PREM"`. r is the present-day reference correction of the
science case (`tools/levante/gia_refcorr.sh`, `CASE=PaleoMIST-sci`): the rebound still to come
at 0 ka, so that the run ends on today's bed (up to 396 m over Quebec-Labrador). A different
Earth or ice history needs its own case; see
[docs/dev/vilma2-paleomist.md](https://github.com/fesmc/kryos/blob/main/docs/dev/vilma2-paleomist.md),
"Present-day reference correction for kryos".

## Timeline

- `[clock]`: `time_init = -80e3` (the first PaleoMIST slice), `time_end = 0`, `dtt = 100` yr
  (one plastic solve and one VILMA2 step per 100 yr), `method = "rel"`, `time_ref = 1950`.
- Cold start (`[kryos] restart = "None"`): the first plastic solve at -80 ka on the reference
  bed with `z_sl = 0`; then VILMA2's relaxation to the 80 ka load (`allow_no_spinup = False`:
  the run needs it) and the solve again on the spun-up Earth.
- **Output:** 1-D every 100 yr (`out_1D`), 2-D every 2.5 kyr (`out_2D`: the PaleoMIST slice
  times), restart bundles every 10 kyr (`out_rst`) and at the end.

## Running

On Levante, a run folder holds the executable, the par file and four links (the layout of a
calibration member, `plastic_esmda.sh member`). `K` is the kryos checkout the executable was
built in:

```bash
K=/work/ba1442/robinson/kryos-ice-plastic
R=$K/runs/<name>
mkdir -p $R && cd $R
cp $K/libkryos/bin/kryos.x .
cp $K/par/kryos_NH-ANT-PAT_plastic_vilma2.nml par.nml
ln -s $K/input input
ln -s /work/ba1442/fesmdata fesmdata
ln -s /work/ba1442/kryos_data kryos_data
ln -s /work/ba1442/robinson/kryos-ice-plastic-maps maps
```

- `fesmdata` and `kryos_data` are shared stores: read them, never write through the links.
- `maps` holds the grid descriptions and the remapping weights (the first run writes the
  weights). It needs `grid_ANT-16KM.txt` (copied from the checkout's `maps/`) and, with
  Patagonia, `grid_PAT-PAL-8KM.txt`. The folder of our runs,
  `/work/ba1442/robinson/kryos-ice-plastic-maps`, has the grid files of NH-32KM, ANT-16KM, PAT-PAL-8KM,
  GLOBAL-0.5DEG and GAUSS-L170 and all their weights: link it (or copy it) rather than starting a new one.
- `kryos.x` takes the par file as its only argument (no `-p`). Change keys in `par.nml` with
  [tools/levante/set_nml.py](https://github.com/fesmc/kryos/blob/main/tools/levante/set_nml.py) (as `plastic_esmda.jl`
  does):

```bash
F=/work/ba1442/robinson/kryos-ice-plastic/runs/taub-fit-v3/report/taub_domain_factor.nc
python3 tools/levante/set_nml.py par.nml \
  "ice_nor.taub_domain_factor_file=\"$F\"" "ice_ant.taub_domain_factor_file=\"$F\"" \
  "ice_pat.taub_domain_factor_file=\"$F\""
```

[tools/levante/plastic_run.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_run.sh) does all of this
in one step: it stages the run folder (the driver, the links, the par file with the factor and scale files set on every
instance, further `group.key=value` overrides) and submits the job:

```bash
tools/levante/plastic_run.sh runs/sci-3d par/kryos_NH-ANT-PAT_plastic_vilma2.nml \
    /work/ba1442/robinson/kryos-ice-plastic/runs/taub-fit-v3/report/taub_domain_factor.nc
tools/levante/plastic_run.sh runs/b1-short par/kryos_NH-ANT_plastic_vilma2_b1.nml "" "" 'clock.time_end=-60e3'
```

Its environment variables (`PARTITION`, `THREADS`, `MEM`, `TIME`, `MAPS`, `SUBMIT=0` to stage only) are in its header.

The job `plastic_run.sh` writes (shared partition, 64 threads, 60 GB, 8 h by default), by hand:

```bash
#!/bin/bash
#SBATCH --account=ba1442
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=30G
#SBATCH --time=08:00:00
export OMP_PROC_BIND=true OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OMP_STACKSIZE=512M OMP_MAX_ACTIVE_LEVELS=1
ulimit -s unlimited
./kryos.x par.nml > out.out 2> out.err
```

submitted from the run folder with `sbatch <job script>`.
At most 30 Levante jobs at a time (queued and running), also for ensembles.

Cost:

| Configuration | Wall time |
|---|---|
| North, 1-D (`_b1`) | 7 min 25 s on 64 threads (old margin); 13.6 min on 16 threads |
| North + Antarctica, 1-D (`_b1`) | 13.4 min on 16 threads; calibration members 11-25 min on 8 threads |
| North + Antarctica (+ Patagonia), 3-D | 3.5-4.2 h on 64 threads (`plastic_ensemble_v1.toml`) |

VILMA2 dominates and barely scales with threads on the 1-D Earth (8 threads with 10 GB is the
cheapest per run on the shared partition). A run continues from a bundle with
`kryos.restart = "restart-<time>"` and `clock.time_init = <time>` (see Restarts in
[Getting started](getting-started.md)).

## Outputs

In the run folder (`[kryos] out_dir = "./"`):

| File | Schedule | Content |
|---|---|---|
| `ice_nor.nc`, `ice_ant.nc`, `ice_pat.nc` | 2.5 kyr | all exports of the plastic instance (`out_vars = ""`): `H_ice`, `z_srf`, `z_bed`, `z_sl`, `f_ice`, `f_grnd`, `dHidt`, `mb_net`, `H_grnd`, ... |
| `ice_*_plastic.nc` | 2.5 kyr (`out_plastic = "out_2D"`) | the backend's stream: `mask`, `sdist`, `tau_b`, `taub_domain`, convergence (`converged`, `n_outer_used`, stress residual) |
| `ice_*_ts.nc` | 100 yr | 0-D exports, e.g. `M_ice` [Gt] |
| `geo_glob_vilma.nc` | 2.5 kyr (`out_vilma = "out_2D"`) | VILMA2's own diagnostics on the Gauss grid: `rsl`, `dz_bed`, `dz_ss`, `z_bed`, `bsl`, `bsl_ive`, ...; the file the GAPSLIP scores read |
| `geo_glob.nc`, `geo_glob_ts.nc` | 2.5 kyr, 100 yr | `dz_bed`, `dz_ss`, `C_ocn` on the Gauss grid; the series |
| `geo_nor.nc`, `geo_ant.nc`, `geo_pat.nc` (+ `_ts`) | 2.5 kyr, 100 yr | the Earth on each ice grid |
| `restart-<time>/` | 10 kyr and the end | restart bundles |
| `timing.txt` | | timing table (`[kryos] timing = True`) |
| `out.out`, `out.err` | | log (from the job above) |

The variable tables: `input/kryos-variables-ice.md` (the exports), `input/kryos-variables-ice-plastic.md`
(the `plastic` stream), `input/vilma-variables.md` (the VILMA2 stream).

## Scoring and plots

GAPSLIP 2.1 scores of one run beside B1 (`B1_l170`) and `I910_it0`, the same scoring as
[Reproducing Gowan et al. (2021)](experiment-gowan.md):

```bash
sbatch tools/levante/plastic_score.sh <kryos run folder> <label>
```

Output `/work/ba1442/robinson/kryos-ice-plastic/runs/score/<label>/out` (`pm_indicators.md`, ...;
the run under its folder name). `pm_compare.jl` reads the run's `geo_glob_vilma.nc` in place.

Several runs against each other (e.g. 1-D against 3-D Earth, Gowan's tau_b against the Stage 1
factors):

```bash
julia --project=<julia-plot> tools/levante/plastic_earth_compare.jl <out folder> <label>=<run folder> ...
```

At least two runs; the first is the reference of the maps. It reads `geo_glob_vilma.nc` and
`<instance>_ts.nc` of `ice_nor`, `ice_ant`, `ice_pat` and writes `compare.md` (ice volume and
ESL; scores by group, by GAPSLIP region and by marine-bed area, with RMS z and RMS z_m =
d / sqrt(sigma^2 + sigma_m^2), sigma_m = 5 m), `sites.csv`, `sites_map.png` (per site, RMS z_m of
each run minus that of the first) and `curves.png` (RSL at the best-sampled sites).

Maps and present-day checks of one run:

```bash
julia --project=<julia-plot> tools/levante/plastic_lgm_map.jl <run folder> <out.png> [mode=age|margins|smb] [age=20000]
julia --project=<tools> tools/levante/plastic_pd_regions.jl <run folder>
```

`plastic_lgm_map.jl` (after Gowan et al. 2021, Fig. 2): the ice surface above the sea level of
the age, RSL against today, the coastlines and the GAPSLIP sites scored; needs the 2-D output of
the plastic instances. `plastic_pd_regions.jl`: today's thickness per region against the
observations of the topography hubs. The Julia environments are described in [Analysis and figures](experiment-analysis.md).

### Results so far

From [docs/dev/plastic-backend.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-backend.md),
1-D Earth (`_b1`), ice volume in m SLE (ratio to PaleoMIST) and ESL at 20 ka:

| | NH 20 ka | ANT 20 ka | ESL 20 ka [m] | RMS z, all index points |
|---|---|---|---|---|
| PaleoMIST | 112.5 | 76.63 | -115.6 | 29.5 |
| North, Gowan's tau_b (baseline) | | | -131.7 | 25.9 |
| North + Antarctica, Stage 1 factors (v2 baseline) | 117.7 (1.05) | 76.80 (1.00) | -124.8 | 26.7 |
| v2 posterior mean (scale file of calibration v2) | 95.9 (0.85) | 100.9 (1.32) | -129.7 | 21.6 |

The 3-D runs (`sci-v2-prior`, `sci-v2-post`, `sci-v3-base` under `/work/ba1442/robinson/kryos-ice-plastic/runs`) are
summarised in the section "The 3-D Earth" of [docs/dev/plastic-backend.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-backend.md): at the same tau_b the 3-D Earth fits
the index points worse globally (RMS z_m 2.88 -> 3.47 with the prior tau_b), better in Hudson Bay (3.70 -> 2.51), the
North American Arctic (3.89 -> 2.49) and Greenland (7.51 -> 6.71), much worse in Fennoscandia (1.8 -> 4.8); ESL at 20 ka
4.4 m lower. The tau_b were calibrated on the 1-D Earth; the weighted ensemble recalibrates them on the 3-D Earth
([Calibration and ensembles](experiment-calibration.md)).
