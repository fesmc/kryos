---
title: "Analysis and figures"
---

The analysis of the plastic and VILMA2 runs is written in Julia, with CairoMakie for the figures, in
`tools/levante/`. Everything runs on Levante, next to the runs. This page lists, for each tool, what it shows, its
command line (from the file header) and the run output it reads.

## Environment

All tools share one Julia project, `tools/levante` (`Project.toml`, `Manifest.toml`: CairoMakie, NCDatasets,
DelimitedFiles on Julia 1.13), instantiated once by `tools/levante/setup.sh` (see
[Getting started](getting-started.md#the-levante-tools)). From the kryos root:

```bash
julia --project=tools/levante tools/levante/<tool>.jl ...      # or kryos_jl tools/levante/<tool>.jl ... after sourcing env.sh
tools/levante/jl.sh <out file> tools/levante/<tool>.jl ...      # the same as a small shared-partition job (1 h, 4 cpus, 16 GB)
```

`julia` is juliaup's, with the 1.13 channel (`KRYOS_JULIA` of `env.sh`: `$HOME/.juliaup/bin/julia +1.13`; Levante's
python module puts an old julia first on `PATH`).

Some conventions are shared by the tools:

- **Data paths:** the tools read FesmData at `FESMDATA` (`env.sh`, default `/work/ba1442/fesmdata`). The GAPSLIP 2.1
  indicators are read through `pm_compare.jl` (`read_indicators`).
- **Includes:** several tools include `pm_compare.jl` or `plastic_mapkit.jl` from the same folder, so they must stay
  together.
- **z_m:** the tools that report z_m use z_m = d / sqrt(sigma^2 + sigma_m^2) with sigma_m = 5 m, as
  `plastic_esmda_report.jl` does.

The run output the tools read:

| File | Written by | Needed for |
|---|---|---|
| `geo_glob_vilma.nc` | the VILMA2 stream of `geo_glob` | RSL, the coastlines, the load; every scoring |
| `<instance>.nc` (`ice_nor.nc`, `ice_ant.nc`, `ice_pat.nc`) | the 2-D output of a plastic instance (`out_2D`) | ice maps, ice spread, implied SMB |
| `<instance>_plastic.nc` | the backend's stream `plastic` (`out_plastic`) | `q_land`, `q_marine` of the implied mass balance |
| `<instance>_ts.nc` | the 1-D series | ice volume per instance |

Calibration and ensemble members write the 2-D ice product every 2.5 kyr (see
[Calibration and ensembles](experiment-calibration.md)).

## Maps of a run

### plastic_lgm_map.jl

Global maps of one kryos plastic run, after Gowan et al. (2021, Fig. 2). The drawing is in `plastic_mapkit.jl`.

```bash
julia --project=tools/levante tools/levante/plastic_lgm_map.jl <run folder> <out.png> [mode=age|margins|smb] [age=20000] [window=2500] \
      [ages=25000,20000,15000,12500,10000,7500] [label=<text>] [lon0=0]
```

| Mode | What it shows |
|---|---|
| `age` (default) | The run at one age. The ice surface elevation above that age's sea level (the plastic instances on their grids; the rest of the world's ice from VILMA2's load). The RSL of that age against today over that age's ocean, filled and contoured. The coastlines of that age and of today, and the target (PaleoMIST) margin. The GAPSLIP 2.1 sites with index points within age +- window, coloured by the RMS of z_m of those points. |
| `margins` | The target margins at `ages` as isochrones, today's coastline, and every site with index points, coloured by the RMS of z_m over all its points (0-80 kyr). |
| `smb` | The implied mass balance at the age (`mass_balance = "implied"`): `smb` on the ice [m w.e. yr-1] (including the land-margin ablation in the margin cells), the totals of `smb` and `cmb` per instance [Gt yr-1], the target margin and the coastlines. |

Triangles mark the ice-sheet data groups (Laurentide, Greenland, Fennoscandia/Barents, Antarctica); circles mark the
far field.

**Inputs:** `geo_glob_vilma.nc`, and the instances' 2-D output at the age (`H_ice`, `z_srf`; `smb`, `cmb` for
`mode=smb`).

### plastic_mapkit.jl

Not a script. This is the drawing module that `plastic_lgm_map.jl` and `plastic_ensemble_sd_map.jl` include. It
provides:

- the Equal Earth projection (Savric et al. 2018), written out in the file, with no GeoMakie;
- fields on structured meshes drawn as triangles, and contours;
- the VILMA2 fields of a run at an age (`geo_glob_vilma.nc`);
- the 2-D fields of the plastic instances (`<instance>.nc` with the grid's `lon2D`, `lat2D`);
- coastlines, ice margins, and the GAPSLIP 2.1 index points scored by `pm_compare.jl`.

The regional fields it reads are `H_ice`, `z_srf - z_sl` and, when written, `smb` and `cmb`.

## Ensembles

### plastic_ensemble_sd_map.jl

The spread of an ensemble at one age on a global map. It shows two standard deviations across the members:

- **RSL** against today (VILMA2), over the ensemble-mean ocean of the age;
- **Ice thickness** of the plastic instances on their grids, where every member has ice (the margin is prescribed).

The map also shows the GAPSLIP 2.1 sites with index points within age +- window, the target margin, and the
coastlines of the ensemble-mean ocean.

```bash
julia --project=tools/levante tools/levante/plastic_ensemble_sd_map.jl <out.png> <member folder> ... [age=20000] [window=2500] \
      [label=<text>] [lon0=0]
```

**Inputs:** member folders with `status` 0, `geo_glob_vilma.nc` and, for the ice, `<instance>.nc` with 2-D output at
the age.

### plastic_ensemble_spread.jl

Where the sea-level data constrain the ice and where they do not. It compares the spread of the surface elevation
z_srf - z_sl across the members of one iteration of a calibration (`plastic_esmda.jl`) with the spread across the
`ngood` members of lowest chi2/n (the members that fit about equally well). The comparison is on the cells with ice
in every member.

```bash
julia --project=tools/levante tools/levante/plastic_ensemble_spread.jl <calibration folder> <iteration> <2-D folder> <out folder> \
      [age=20000] [ngood=10]
```

**Inputs:**

- `<2-D folder>/mNNN`: the members of the iteration, rerun with `ice_nor.nc` and `ice_ant.nc` (2-D, 2.5 kyr records).
- `iter_kk/members.csv` of the calibration: chi2/n = d' R^-1 d / n over the assimilated index points.

**Output:**

| File | Content |
|---|---|
| `spread_<age>ka.png` | per grid: the mean of the good members, the sd of all members, the sd of the good members, and the ratio of the two sds (1: the data do not narrow the spread) |
| `members.png` | chi2/n against the ice volume per scale group at the age and against the ESL |
| `spread_<age>ka.md` | the members by chi2/n; the effective ensemble size of likelihood weights exp(-n chi2 / 2) and of tempered ones; the median sd per grid |

## Comparing runs against GAPSLIP

### pm_compare.jl

The scorer. It was written for the PaleoMIST test of the standalone VILMA2 (the runme cases of
`tools/levante/paleomist/`, [Reproducing Gowan et al. (2021)](experiment-gowan.md)), and its scoring functions
are used by every tool below and by `plastic_esmda.jl`. A `key=folder` argument points a run key at any folder, and
`read_run` also reads a kryos run's VILMA2 stream `geo_glob_vilma.nc` in place.

```bash
julia --project=tools/levante tools/levante/pm_compare.jl <fesmdata> <runs folder> <out folder> [<key>=<folder> ...]
```

How it scores:

- **Grid:** everything is put on the 0.5 degree PaleoMIST grid; the VILMA2 fields are interpolated bilinear in
  lon-lat.
- **Model RSL:** relative to present at the GAPSLIP 2.1 samples (primary data, age <= 80 ka), bilinear in space and
  linear in time between the 2.5 kyr records.
- **Index points:** d = model - rsl and sigma^2 = sigma_z^2 + (dRSL/dt sigma_t)^2. It reports the RMS of d/sigma and
  the fraction with |d/sigma| <= 2.
- **Limiting data:** a point counts as violated if no model RSL over its 2-sigma age range reaches its 2-sigma range.
- **Grouping:** results are given by the groups of Gowan's regions (`GROUPS`) and by the marine-bed areas (`FOCUS`).

**Output:**

- `pm_metrics.md`, `pm_metrics_B1_PM.csv`, `pm_sites.csv`, `pm_gmsl.csv`;
- `pm_fields.nc` (read by `pm_maps.jl`);
- `pm_indicators.md` (the scores), `pm_ind_samples.csv`, `pm_ind_sites.csv`, `pm_ind_curves.csv`.

`tools/levante/plastic_score.sh` scores one kryos run this way, beside the references B1 (`B1_l170`) and `I910_it0`
of `PM_REF`, as a small job it submits itself; the report goes to `<score folder>/out` (default
`$KRYOS_WORK/score/<label>/out`):

```bash
tools/levante/plastic_score.sh <kryos run folder> <label> [<score folder>]
```

### plastic_earth_compare.jl

The plastic ice on different Earths. It compares several kryos runs:

- GAPSLIP 2.1 scores by data group, by GAPSLIP region, by marine-bed area and by site;
- the change per site against the first run, on maps;
- RSL curves at the best-sampled sites;
- the ice volume per instance and the ESL.

```bash
julia --project=tools/levante tools/levante/plastic_earth_compare.jl <out folder> <label>=<run folder> ...
```

**Inputs:** kryos run folders with `geo_glob_vilma.nc` and `<instance>_ts.nc`.

**Output:**

| File | Content |
|---|---|
| `compare.md` | tables of volume and ESL, by group, by GAPSLIP region and by marine-bed area |
| `sites.csv` | per site and run: number of index points, RMS z, RMS z_m, mean d |
| `sites_map.png` | per site, RMS z_m of each run minus that of the first (blue: better); Antarctica in polar view |
| `curves.png` | RSL at the sites with the most index points per group |

The 3-D against 1-D comparison it produced is in the section "The 3-D Earth" of [docs/dev/plastic-backend.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-backend.md).

### plastic_residual_corr.jl

The structure of the misfit of the index points. It is used to choose the model-error part of the observation error:
R = diag(sigma_i^2) + sigma_m^2 C, with C_ij = exp(-dist / L_x) exp(-|dt| / L_t).

```bash
julia --project=tools/levante tools/levante/plastic_residual_corr.jl <out folder> <label>=<run folder> ...
```

For each run, it takes d = model - obs and the scorer's 1 sigma per index point. Then, within each group of
`pm_compare.jl`, each GAPSLIP region and each sector, it computes:

- **Correlation:** rho = mean(d_i d_j) / sqrt(mean(d_i^2) mean(d_j^2)) over pairs of index points. It is not centred,
  so a common bias also counts as a correlated error. rho is binned against distance (pairs with |dt| <= 1 kyr) and
  against age difference (pairs within 25 km).
- **Fit:** an exponential c exp(-h / L), fitted to the bins and weighted by the number of pairs.
- **Model error:** sigma_m^2 = mean(d^2) - mean(sigma_i^2), the error beyond the observations.

**Output:**

- `corr_distance.png`, `corr_time.png`: rho by group, with v2's assumption (150 km, 2 kyr);
- `corr_map.png`: the fitted L_x and sigma_m per GAPSLIP sector;
- `residual_corr.md`: the tables.

Applied to the v2 posterior on the 1-D Earth, this gave the `[likelihood]` values of the weighted ensemble.

## Implied mass balance

See [Implied mass balance](experiment-implied-mb.md) for the method.

### plastic_implied_pd.jl

The implied mass balance of a kryos plastic run with `mass_balance = "implied"` at the present day, against
observations:

- **SMB:** MAR v3.14 for Greenland and RACMO2.4p1 for Antarctica, both 1981-2010.
- **Surface speed:** Joughin et al. (2018) and Mouginot et al. (2019).

All of these are FesmData SMB and Velocity products.

```bash
julia --project=tools/levante tools/levante/plastic_implied_pd.jl <run folder> <out folder> [age=0]
```

**Inputs:**

- `ice_nor.nc` and `ice_ant.nc` with the 2-D exports `smb`, `cmb`, `ux_s`, `uy_s`, `H_ice` (`out_vars` all);
- `<instance>_plastic.nc` with `q_land`, `q_marine` (`out_plastic`).

**Grids:** for Greenland, the observations are on GRL-32KM, and the model (NH-32KM) is taken at the nearest cell
centre with ice, within 30 km. For Antarctica, both are on ANT-16KM. The comparison covers the model's ice cells,
and also the cells where the observation has ice.

**Output:** `implied_pd.md` (integrals [Gt yr-1] and speed statistics) and `implied_pd.png` (maps and a speed
scatter).

### plastic_closure_maps.jl

Maps of the present-day implied mass balance of `tests/plastic_closure_fit.f90` parameter sets, one NetCDF file per
set. For each set it maps:

- the implied SMB rho_i (div q - q_land) on the ice, smoothed conservatively as the backend's `mb_smooth_km` does
  (scatter form, sums unchanged);
- log10(model / observed) of the surface speed.

Optionally, an observed SMB on the same grid is shown as the first column.

```bash
julia --project=tools/levante tools/levante/plastic_closure_maps.jl <out.png> <sd_km> [ref=<smb.nc>] <label>=<set.nc>[:inv] ...
```

The suffix `<set.nc>:inv` maps the set's joint inversion instead (the `prior_file` option of `plastic_closure_fit`:
`smb_inv`, `speed_inv`, `q_marine_inv`), unsmoothed.

For each set, the tool prints the integrals [Gt yr-1] (net, accumulation, ablation, outflow to the ocean) and the
ablation area (the fraction of ice cells with smoothed SMB < 0). It prints the same for `ref` over the model's ice.

## Other report tools

| Tool | Command | What |
|---|---|---|
| `plastic_esmda_report.jl` | `julia --project=tools/levante tools/levante/plastic_esmda_report.jl <calibration folder> <out folder> [<prior slices.nc> <posterior slices.nc>]` | a calibration: per-iteration scores, factor curves (`theta_factors.png`), ice volume of the mean members, present-day thickness per region, an offline Greenland check from `plastic_slices` output; `report.md` |
| `plastic_esmda_diag.jl` | `julia --project=tools/levante tools/levante/plastic_esmda_diag.jl <calibration folder> <out folder> <prior run> <posterior run>` | where the posterior moves the ice: volumes by scale group and ESL, the prior-ensemble sensitivity, residuals by age, the Antarctic sites, ANT-16KM maps at 20 ka |
| `plastic_taub_fit_report.jl` | see [Calibration and ensembles](experiment-calibration.md) | Stage 1: L-curves, lambda, the factor file, maps |
| `plastic_pd_regions.jl` | `julia --project=tools/levante tools/levante/plastic_pd_regions.jl <run dir> [<fesmdata>]` | present-day thickness per region against the FesmData topography hubs |
| `plastic_maps.jl` | `julia --project=tools/levante tools/levante/plastic_maps.jl <out.png> <ages ka, comma separated> <label>=<slices.nc> ...` | the offline check (`tests/plastic_slices.f90`): H - PaleoMIST H_ice at given ages |
| `pm_maps.jl` | `julia --project=tools/levante tools/levante/pm_maps.jl <compare folder>` | figures from `pm_compare.jl` output (the standalone VILMA2 PaleoMIST test) |
