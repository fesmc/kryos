---
title: "Imposed ice"
---

How to run a published ice history forced onto VILMA2 with a given Earth. There is no ice model: the ice is read from
the history, so a change of score is the Earth's or the solver's. These runs are the I rows of the
[results](../science/results.md) and every run of the [Earth search](earth-search.md). The concept and its role in the
reconstruction: [Methods](../science/methods.md). Status as of 2026-11-16.

Every imposed-ice run is a **case** of `tools/levante/gia_refcorr.sh`: standalone VILMA2 (`$VILMA2/bin/vilma.x`, not
`kryos.x`) with the `&vilma` group of a kryos par file, under a FesmData ice history, with its own present-day correction
r iterated until the run ends on today's bed. What r is and why a run needs it: [Global hub, the reference topography](../global-hub.qmd#the-reference-topography-refcorr).

## One case

```bash
CASE=ice-PMIST_E1D-PMIST PAR=par/kryos_NH-ANT_plastic_vilma2_b1.nml NCPU=16 tools/levante/gia_refcorr.sh
```

The script submits itself (`KRYOS_PART_LONG`, 2 days, `NCPU` cpus, default 128, 60 GB) and runs the iterations within that
job. The environment of a case:

| variable | meaning | default |
|---|---|---|
| `CASE` | the case name: everything goes to `$KD/refcorr/$CASE` | required |
| `PAR` | the kryos par file whose `&vilma` group is the Earth | required |
| `HIST` | the ice history, a FesmData IceHistory product | PaleoMIST on GLOBAL-0.25DEG |
| `ZREF`, `XHIST`, `XZREF` | for a history without a bed (GLAC-1D 80 ka): today's bed from the PaleoMIST product; the history on the 0.5° background grid | none |
| `INTERP` | `linear` (records linear in time) or `steps` (SELEN's Heaviside steps) | `linear` |
| `TOL`, `NMAX` | stop when the 0 ka bed misfit is below `TOL` m, or after `NMAX` runs | 1.5, 8 |
| `R0` | run once with this r, without iterating (a borrowed r); no reference files | none |
| `KD` | the parent of `refcorr/` | `$KRYOS_WORK/kryos_data` |
| `KEEP` | `last`: once converged, keep only the converged iteration's `vilma_out.nc` | `all` |
| `DOMAINS` | the domain topographies to correct, `<domain>:<hub grid>:<topo file under fesmdata> ...` | North (NH-16KM, NH-32KM), Antarctica |
| `NCPU`, `MEM`, `NTHR` | job size and OpenMP threads | 128, 60G, `NCPU` (at most 64) |

A 1-D Earth needs 16-32 cpus and about 3 GB at degree 170. A 3-D Earth needs a node or half of one.

Run times on albedo: 1-D 3-5 min per iteration, VILMA2's 1Db profile about 17 min, 3-D 2-30 h depending on the Earth (soft
upper mantles are slow); 3-5 iterations to converge.

## A case folder

`$KD/refcorr/<CASE>/`:

| file | content |
|---|---|
| `iter.log` | the iterations, the 0 ka bed misfit of each, and the last line `== done` (or `not converged`, or `failed`) |
| `it<n>/` | iteration n, a runme run of `vilma_rc.nml`: `vilma_out.nc`, `out.out`, `status` |
| `r_<n>.nc`, `r.nc` | the corrections; `r.nc` links to the converged one |
| `<grid>/<grid>_REFCORR.nc`, `<domain>/<hub>/<topo>-refcorr` | the reference files a kryos run reads (today's bed + r) |
| `data/` | `force.nc` (the ice history as VILMA2 forcing) and `ref.nc` (the history's 0 ka bed and ice) |
| `provenance.txt` | appended at every start: the environment, the git commits of kryos and VILMA2, md5 of the par file and of `vilma.x`, the inputs |

The case is restartable: an iteration whose `status` is 0 is not rerun. A case that ended in `not converged` or `failed` is
not retried by the pools below: delete its folder to rerun it.

## The ladder rows

The I rows are made by `experiments/refs-ice.sh`, one call per row:

```bash
experiments/refs-ice.sh I1            # or I1 ... I8, or all
```

| id | case(s) | par file (Earth) | ice |
|---|---|---|---|
| I1 | `ice-PMIST_E1D-PMIST` | `kryos_NH-ANT_plastic_vilma2_b1.nml` (Gowan's 1-D) | PaleoMIST |
| I2 | `ice-PMIST_E1D-B21-6`, `…-layers` | `…_bagge1d.nml`, `…_bagge1d_layers.nml` (Bagge 6's radial mean) | PaleoMIST |
| I3 | `ice-PMIST_E3D-B21-6` | `kryos_NH-ANT-PAT_plastic_vilma2.nml` (Bagge 6) | PaleoMIST |
| I4 | `ice-PMIST_E3D-B21-<k>`, k = 1-5, 7, 8 | the 3-D par file with `visc_3d_structure = k` | PaleoMIST |
| I5 | `ice-GLAC1D_E1D-PMIST` | Gowan's 1-D | GLAC-1D 80 ka (nn9894) |
| I6 | `ice-GLAC1D_E1D-B21-6-layers` | Bagge 6's radial mean as layers | GLAC-1D |
| I7 | `ice-GLAC1D_E3D-B21-6` | Bagge 6 | GLAC-1D |
| I8 | `ice-GLAC1D-9927_E1D-PMIST`, `…_E3D-B21-6` | Gowan's 1-D, Bagge 6 | GLAC-1D, North America nn9927 |

The cases go to `$KRYOS_WORK/kryos_data/refcorr/`. The GLAC-1D histories (FesmData `IceHistory-GLAC1D_80k`) have no bed:
today's comes from PaleoMIST (`ZREF`). Before that FesmData release, point `GLAC1D_GLOBAL` at a local build of the product.

## Many cases: pools

Many small cases (the Earth search: hundreds of 1-D Earths) run in few jobs with `tools/levante/gia_refcorr_pool.sh`. A
pool job has `NWORK` workers, each on its own share of the cores. Each worker runs the next case of a list that is neither
done nor claimed, until the list is done or the job's time is nearly up:

```bash
LIST=cases.txt NJOB=4 NWORK=8 NCPU=128 MEM=64G TIME=12:00:00 tools/levante/gia_refcorr_pool.sh
```

The list has one case per line, `CASE PAR [VAR=value ...]`, the variables being those of `gia_refcorr.sh` for that case
(`HIST`, `ZREF`, `R0`, ...). A worker claims a case with a link `<case>/claim -> <job id>`. A claim of a job that no longer
runs is taken over, so resubmitting the same list resumes the cases of a job that was cut off. A worker claims no case after
`TIME - MARGIN`, nor when less than 1.2 times the longest finished case of the list is left. `PART` sets the partition and
QOS (default `KRYOS_PART`; `PART="$KRYOS_PART_LONG"` with a `TIME` over 12 h on albedo). Each case writes `pool.out`; the
job writes `$KD/refcorr/pool-<job>.out`.

## A borrowed r

A run with `R0=<r.nc>` is one pass with that r instead of 3-5 iterations. The 3-D runs of the Earth search take the r of
the 1-D run of their own radial profile, and the Earth-only ensemble `earth-pca` (`experiments/earth-pca.sh`) took the r
of PaleoMIST-sci. Such a run ends near, not on, today's bed: its score includes that offset.

## Scoring

An imposed-ice case is scored as a reference of an ensemble config: `[[references]]` in
`tools/levante/plastic_ensemble_v2.toml` (its `run` is the case's converged iteration), then

```bash
experiments/scores.sh <ensemble folder>
```

writes `<folder>/iter_01/scores.csv` (one row per series, data group and window). The data, the groups and the scores:
[Data and scores](../science/data-scores.md). The mechanics of the scoring folder: [Ensembles](ensembles.md#start-monitor-score).

## History: the reproduction of Gowan et al. (2021) {#history}

Before `gia_refcorr.sh`, PaleoMIST's ice was imposed on VILMA2 as run **B1** (2026-10-31): Gowan's 1-D Earth, 931 kg m⁻³,
2.5 kyr steps and no present-day correction. It reproduced PaleoMIST to about 10 m RMS globally and scored as PaleoMIST
against GAPSLIP. Its cases (`tools/levante/paleomist/pm_*.nml`), the tools `pm_prep.jl`, `pm_project.sh`,
`pm_iterref.sh`, `pm_visc1d.sh`, `pm_selen.sh`, `pm_compare.jl` and `pm_maps.jl`, and the analysis of the remaining
differences (sea-level equation, reference state, ice density) are in
[docs/dev/vilma2-paleomist.md](https://github.com/fesmc/kryos/blob/main/docs/dev/vilma2-paleomist.md). The published runs
`B1_l170`, `ref/I910_it0` and `B3_l170` are in `$PM_REF` (default `$KRYOS_DATA/paleomist`), read-only. I1 replaces B1 as
the reference: our solver with PaleoMIST's ice and Earth, 910 kg m⁻³, linear in time, with its own r.

`pm_compare.jl` still reads any standalone or kryos run and is the site-by-site comparison tool:
[Analysis and figures](analysis.md#comparing-runs-against-gapslip).
