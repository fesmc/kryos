---
title: "Online ice"
---

How to run our own ice: the perfectly-plastic ice sheet ([PlasticIceSheet](../backend-plastic.md)) solved at every 100 yr
step on PaleoMIST's margins and Gowan's basal shear stress (τ_b) domains, on the bed and sea level of the global VILMA2
Earth ([VILMA2](../backend-vilma2.md)) that its load deforms, 80 ka to present. The ice is not integrated (no surface,
basal or calving flux): a run is set by τ_b, the margins and the Earth. These runs are the T, D and E rows of the
[results](../science/results.md); the concept: [Methods](../science/methods.md). Status as of 2026-11-16.

- **Program:** `kryos.x` built with `plastic=1 vilma=1` (`yelmo=0` is enough), linked against the PlasticIceSheet library
  ([Getting started](../getting-started.md)).
- **Config:** the par files below ([Par files](../reference-par-files.md)).
- **Running:** one run is staged and submitted with [runme](https://github.com/fesmc/runme) from the checkout; the ladder
  rows with `experiments/refs-tau.sh`, the ensembles with `experiments/ens.sh` ([Ensembles](ensembles.md)).
- **Design notes:** [docs/dev/plastic-backend.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-backend.md).

| Par file | Domains | Earth | Reference |
|---|---|---|---|
| [kryos_North_plastic.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_North_plastic.nml) | North | none (today's bed, `z_sl = 0`): smoke test | |
| [kryos_North_plastic_vilma2_b1.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_North_plastic_vilma2_b1.nml) | North | Gowan's 1-D | relaxed, no correction |
| [kryos_North_plastic_vilma2.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_North_plastic_vilma2.nml) | North | 3-D, Bagge structure 6 | today's + r |
| [kryos_NH-ANT_plastic_vilma2_b1.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT_plastic_vilma2_b1.nml) | North, Antarctica | Gowan's 1-D | relaxed, no correction |
| [kryos_NH-ANT_plastic_vilma2.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT_plastic_vilma2.nml) | North, Antarctica | 3-D, Bagge 6 | today's + r |
| [kryos_NH-ANT-PAT_plastic_vilma2.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT-PAT_plastic_vilma2.nml) | North, Antarctica, PatagoniaPaleo | 3-D, Bagge 6 | today's + r |
| [kryos_NH-ANT-PAT_plastic_vilma2_nh16.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT-PAT_plastic_vilma2_nh16.nml) | the same, North on NH-16KM | 3-D, Bagge 6 | today's + r |

`kryos_NH-ANT-PAT_plastic_vilma2.nml` is the par file of the ensembles ens-v1 and ens-v2. Its Earth, Bagge's structure 6,
is one Earth among those the [Earth search](../science/earth-search.md) tests, not a preferred one: another Earth is set in
the `[vilma]` group (below). `kryos_NH-ANT-PAT_plastic_vilma2_inv.nml` adds the mass-balance inversion
([Implied mass balance](../implied-mb.md)). The `_nh16` file runs the Northern Hemisphere on 16 km; 32 km stays the default
(16 km changes the scores little and costs 2 m of ESL).

## Domains and instances

| Domain | Hub and ice grid | Ice | Geo port | Topography (`[topo_*] path`) |
|---|---|---|---|---|
| `North` | NH-32KM | `ice_nor` | `geo_nor` | `TOPO-GEBCO2026` (1-D) or `TOPO-GEBCO2026-refcorr` (3-D) |
| `Antarctica` | ANT-16KM | `ice_ant` | `geo_ant` | `TOPO-BedMachine-v4` or `TOPO-BedMachine-v4-refcorr` |
| `PatagoniaPaleo` | PAT-PAL-8KM | `ice_pat` | `geo_pat` | `TOPO-GEBCO2026-refcorr` (IceBoost v2 ice) |

- Each ice instance is `backend = "plastic"` on its hub grid (`grid = ""`), every step (`dt = 0`), `drive = "taub"`,
  `mode = "surface"` (`"full"` does not converge on 32 km beds), densities 910 and 1028 (`phys_const = "Earth"`). A domain
  has no climate, surface or shelf: the plastic ice needs none.
- Each geo port is `backend = "vilma2"` on the ice grid, `dt = 100` (it must equal `[geo_glob] dt`, else the run stops).
- The global Earth `geo_glob` (`backend = "vilma2"`, Gauss grid of `lmax` 170, `dt = 100`, `ice_bg = "globe_ice"`). The
  ice outside the domains is the PaleoMIST history: `[globe_ice]` on GLOBAL-0.5DEG, `ext = "file"`, the background X
  (`[X_z_bed]`, `[X_H_ice]`) plus the anomaly of `[X_dH_ice]`. See [Global hub](../global-hub.qmd).

## Ice inputs

The margins and τ_b come from FesmData IceHistory (v2.3.0; PAT-PAL-8KM v2.4.0) on each grid:

- Margins: `margin_file = "fesmdata/{domain}/{grid_name}/{grid_name}_IceHistory-PaleoMIST.nc"` (33 slices, 80 ka to
  present), ice where `f_ice` > 0.5, interpolated in time by the signed distance (`margin_interp = "sdist"`); outside the
  slices' range the run stops (`margin_extrapolate = "error"`).
- τ_b: `taub_method = "domains"`, `taub_domain_file = "fesmdata/{domain}/{grid_name}/{grid_name}_IceHistory-PaleoMIST_taub.nc"`
  (Gowan's domains and their values through time), times two optional factor files per instance:

| Key (per instance) | Shipped value | Content |
|---|---|---|
| `taub_domain_factor_file` | `kryos_data/plastic/taub_domain_factor_v3.nc` | the present-day factors (Stage 1): one per domain, constant in time |
| `taub_scale_file` | `""` | the multipliers of an ensemble member: one per scale group, linear in time between nodes |

How τ_b is composed from these, and how the factors are fitted:
[PlasticIceSheet, basal shear stress](../backend-plastic.md#basal-shear-stress).

- `""` is no file: with both keys `""` the τ_b is Gowan's, unscaled. With runme an empty value is
  `-p ice_nor.taub_domain_factor_file=` (nothing after `=`), one per instance.
- One file of each kind holds the domains of all grids: set the same path on every instance.
- The scale file of an ensemble member is written by `plastic_ensemble.jl` into `<ensemble>/iter_01/scale/mNNN.nc`.
- The fixed-bed smoke file and the 3-D North-only file do not list the two keys; they take the default `""` of
  `input/kryos_defaults.nml`.

## The Earth

The `_b1` files and their 3-D counterparts differ only in the `[vilma]` group, the background X and the topography paths
(`diff par/kryos_NH-ANT_plastic_vilma2_b1.nml par/kryos_NH-ANT_plastic_vilma2.nml`):

| | `_b1`: Gowan's 1-D Earth | 3-D par files: Bagge structure 6 |
|---|---|---|
| `[vilma]` layers | `n_layer = 12`: 120 km elastic, 4e20 to 670 km, 4e22 to the CMB | `n_layer = 11`: 80 km elastic, 4e20 to 670 km, 1e22 below |
| lateral viscosity | `l_visc_3d = .false.` | `l_visc_3d = .true.`, `visc_3d_file = "fesmdata/Global/GLOBAL-T170/GLOBAL-T170_SolidEarth-Bagge2021.nc"`, `visc_3d_structure = 6` |
| spin-up | `pre_spinup_1d = .false.`, `equil_time_max = 2.0e6` | `pre_spinup_1d = .true.`, `equil_time_max = 2.0e4` |
| `[X_z_bed]`, `[X_H_ice]` | PaleoMIST's 0 ka bed and ice (`mode = "transient"`, read at 0 ka) | `kryos_data/refcorr/PaleoMIST-sci/{grid_name}/{grid_name}_REFCORR.nc` (PaleoMIST's 0 ka bed + r), `mode = "static"` |
| `[topo_*] path` | today's topography, no correction | `kryos_data/refcorr/PaleoMIST-sci/.../{grid_name}_TOPO-...-refcorr.nc` (today's + r) |

Both use `lmax = 170`, `earth = "PREM"`. r is the present-day correction of the case PaleoMIST-sci (Bagge 6 under
PaleoMIST's ice); its definition: [Global hub, the reference topography](../global-hub.qmd#the-reference-topography-refcorr).

**Another Earth.** Change the `[vilma]` group: a 1-D profile (`n_layer`, `r_top`, `r_bot`, `eta`, `rheology`), another Bagge
structure (`visc_3d_structure = k`), or a sampled field (`visc_3d_file = <field.nc>`, `visc_3d_structure = 0`). The Earth
search writes such par files (`plastic_ensemble.jl earth-search-sample3d`): [VILMA2, Earth structure](../backend-vilma2.md#earth-structure-1-d-and-3-d).
The reference topography must then be the r of that Earth: run its imposed-ice case first ([Imposed ice](imposed-ice.md))
and point `[topo_*]` and `[X_*]` at its files.

## Timeline

- `[clock]`: `time_init = -80e3` (the first PaleoMIST slice), `time_end = 0`, `dtt = 100` yr (one plastic solve and one
  VILMA2 step per 100 yr), `method = "rel"`, `time_ref = 1950`.
- Cold start (`[kryos] restart = "None"`): the first plastic solve at -80 ka on the reference bed with `z_sl = 0`; then
  VILMA2's relaxation to the 80 ka load and the solve again on the spun-up Earth.
- Output: 1-D every 100 yr (`out_1D`), 2-D every 2.5 kyr (`out_2D`), restart bundles every 10 kyr (`out_rst`) and at the end.

## Running

A run folder holds the executable, the par file and the links `input`, `fesmdata`, `kryos_data`, `maps`, staged and
submitted with runme from the checkout ([Getting started](../getting-started.md#running-with-runme)). On either machine:

```bash
cd $KRYOS
runme -r -s -q $KRYOS_QUEUE -w 08:00:00 --omp 64 -m 60G -o $KRYOS_WORK/runs/sci-3d -n par/kryos_NH-ANT-PAT_plastic_vilma2.nml
runme -r -s -q $KRYOS_QUEUE -w 01:00:00 --omp 16 -m 30G -o $KRYOS_WORK/runs/b1-short -n par/kryos_NH-ANT_plastic_vilma2_b1.nml \
    -p clock.time_end=-60e3
```

- runme copies `libkryos/bin/kryos.x` (build first: `make kryos plastic=1 vilma=1`) and the par file into the run folder,
  and writes the `-p` keys into the copy (a key must exist in the par file or in `input/kryos_defaults.nml`).
- A scale file on every instance: `-p ice_nor.taub_scale_file=<file> ice_ant.taub_scale_file=<file> ice_pat.taub_scale_file=<file>`.
- `fesmdata` and `kryos_data` are the shared stores: read them, never write through the links.
- The job runs `./kryos.x <par file>` with `OMP_PROC_BIND=true OMP_STACKSIZE=512M` and `ulimit -s unlimited`, logs to
  `out.out` and `out.err`, and writes `status` (exit code and seconds) when it ends.
- Without `-r` runme only stages the folder.

The ladder rows T1, D1 and D3 are made from ens-v2's member m000 by `experiments/refs-tau.sh <ens-v2 folder> <id>`
(T1: Gowan's τ_b on Gowan's 1-D Earth; D1: with the present-day factors; D3: the Northern Hemisphere on 16 km).

Cost:

| Configuration | Wall time |
|---|---|
| North + Antarctica, 1-D (`_b1`) | 13.4 min on 16 threads |
| North + Antarctica + Patagonia, 3-D | about 3 h on 64 threads (albedo, ens-v2), 3.5-4.2 h (Levante, ens-v1); NH-16KM about 3.3 h |

VILMA2 dominates the cost (the plastic ice about 0.7 s per kyr). A run continues from a bundle with
`kryos.restart = "restart-<time>"` and `clock.time_init = <time>` ([Getting started, restarts](../getting-started.md#restarts)).
Job limits per machine: [Getting started](../getting-started.md#the-hpc-tools).

## Outputs

In the run folder:

| File | Schedule | Content |
|---|---|---|
| `ice_nor.nc`, `ice_ant.nc`, `ice_pat.nc` | 2.5 kyr | the exports of the plastic instance: `H_ice`, `z_srf`, `z_bed`, `z_sl`, `f_ice`, `f_grnd`, ... |
| `ice_*_plastic.nc` | 2.5 kyr | the backend's stream: `mask`, `sdist`, `tau_b`, `taub_domain`, convergence |
| `ice_*_ts.nc` | 100 yr | 0-D exports, e.g. `M_ice` [Gt] |
| `geo_glob_vilma.nc` | 2.5 kyr | VILMA2's diagnostics on the Gauss grid: `rsl`, `dz_bed`, `dz_bed_rate`, `z_bed`, `bsl`, ...; the file the RSL, ESL and GNSS scores read |
| `geo_glob.nc`, `geo_glob_ts.nc` | 2.5 kyr, 100 yr | `dz_bed`, `dz_ss`, `C_ocn` on the Gauss grid; the series |
| `geo_nor.nc`, `geo_ant.nc`, `geo_pat.nc` (+ `_ts`) | 2.5 kyr, 100 yr | the Earth on each ice grid |
| `restart-<time>/` | 10 kyr and the end | restart bundles |
| `timing.txt` | | timing table |

The variable tables: `input/kryos-variables-ice.md`, `input/kryos-variables-ice-plastic.md`, `input/vilma-variables.md`.

## Scoring and plots

A run is scored as a reference or a member of an ensemble config (`experiments/scores.sh`;
[Ensembles](ensembles.md#start-monitor-score)). For a quick look at one run, the older site scores against B1:

```bash
tools/levante/plastic_score.sh <kryos run folder> <label> [<score folder>]
```

Runs against each other (`plastic_earth_compare.jl`), maps (`plastic_lgm_map.jl`) and present-day checks
(`plastic_pd_regions.jl`): [Analysis and figures](analysis.md).
