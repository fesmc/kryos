---
title: "The Earth search"
---

How to run the Earth search: PaleoMIST's and GLAC-1D's ice imposed on sampled Earths, scored, emulated and cut, in waves. The
question, the parameters, the waves and what they found: [The Earth search](../science/earth-search.md). The design record:
[docs/dev/earth-search-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/earth-search-design.md). Status as of
2026-11-16.

Everything goes through one script, `experiments/earth-search.sh`, on albedo or Levante:

| step | command | does |
|---|---|---|
| wave A | `sample-a` | 200 1-D Earths of 5 radial parameters (seed 1) and Gowan's Earth in those parameters: `$W/earths`, `$W/cases.txt` (402 cases) |
| | `run-a [njob]` | the cases in pools of 16 workers × 8 cores (default 9 jobs; resubmit to go on) |
| scores | `score [b\|b2]` | a job: the scoring folder, the GNSS error model refitted on the cases, `scores.csv`, the emulators; with `b` the wave-B emulator, with `b2` the wave-B2 one |
| wave B | `sample-b` | 50 3-D Earths of 13 parameters, the radial ones from the space wave A left under either ice: `$W/earthsB` |
| | `run-b 1d [njob]` | the 1-D runs of their radial profiles (pools of 4 × 8 cores, default 25 jobs) |
| | `run-b 3d [njob]` | the 3-D runs, one pass each with the r of the profile's 1-D run (pools of 2 × 32 cores for 36 h); after `run-b 1d` |
| wave B2 | `sample-b2` | 30 3-D Earths E101-E130 in the boxes of `experiments/earth-search-b2-boxes.toml`: `$W2/earths` |
| | `run-b2 1d\|3d [njob]` | as `run-b`; the 3-D pools run 44 h (`KRYOS_PART_LONG`): the soft Earths of box S are slow |

Run `experiments/earth-search.sh` without arguments for the full usage.

## Environment

| variable | meaning | default |
|---|---|---|
| `W` | the search folder | `$KRYOS_WORK/earth-search` |
| `W2` | the wave-B2 folder (its cases in `$W2/kryos_data/refcorr`) | `$KRYOS_WORK/earth-search/waveB2` |
| `KD` | the parent of the cases' `refcorr/` | `$W/kryos_data` |
| `ENS` | the ens-v2 folder whose members and config the scoring folder links (read-only) | the ens-v2 run on albedo |
| `ERRMODEL` | the error model of the RSL groups (the refit of refs-test, Fennoscandia and Barents split) | the refs-test `errmodel.toml` on albedo |
| `GLAC1D` | the folder of the GLAC1D_80k IceHistory products (until FesmData has them) | a local build on albedo |
| `BASIS`, `REGIONS` | the Earth basis v2 and the regional masks | `$KRYOS_WORK/kryos_data/plastic/earth_basis_v2.nc`, `$KRYOS_DATA/plastic/earth_regions_v2.nc` |

The script sets `KEEP=last` (only the converged iteration's output is kept) and `DOMAINS` to the North on NH-32KM. The 1-D cases do
not count against the job limit (many workers share a job); the 3-D ones do, and `run-b 3d` and `run-b2 3d` wait before each submission while you have 50 jobs or more.

## What a step makes

- **Earths** (`plastic_ensemble.jl earth-search-sample`, `earth-search-sample3d`, `earth-search-sample3d-box`): per member
  `E<nnn>.nml` (the par file: a 1-D profile, or the 3-D template with `visc_3d_file = E<nnn>.nc`, `visc_3d_structure = 0`),
  `E<nnn>.nc` (the 3-D field), `E<nnn>_1d.nml` (the 1-D par file of its radial profile), and `members.csv` (the parameters) and
  `origin.csv`. How the field is built: [VILMA2, Earth structure](../backend-vilma2.md#earth-structure-1-d-and-3-d).
- **Cases** (`gia_refcorr.sh` through `gia_refcorr_pool.sh`): `$KD/refcorr/ice-<ICE>_E1D-ES1-E<nnn>` (wave A),
  `ice-<ICE>_E1D-ES2-E<nnn>` and `ice-<ICE>_E3D-ES2-E<nnn>` (wave B, B2), ICE = PMIST or GLAC1D. A 1-D case iterates its own r;
  a 3-D case runs once with `R0` = its profile's r. Pools, claims and resubmission: [Imposed ice, pools](imposed-ice.md#many-cases-pools).
- **Scores** (`score`): `$W/score/ens/`, a linked copy of ens-v2 with its config (the split RSL groups, the GNSS groups, ESL with
  the bar) and every case as a reference; `iter_01/scores.csv`. `config_nognss.toml` is the same without GNSS in the cut.
- **Emulators** (`plastic_earth_emulate.jl`): `$W/emulate` (wave A), `$W/emulate-nognss`, `$W/emulateB` (wave B), `$W2/emulateB2`
  (wave B2: the radial GPs on every 1-D run, over the union of the boxes): `summary.md`, `gp.csv`, `loo.csv`, `cases.csv`,
  `sample.csv`, `nroi_*.png` ([Analysis and figures](analysis.md#the-earth-search)).

## Watching the runs

A pool job logs to `$KD/refcorr/pool-<job>.out` (`worker <w>: <case> <start>`, then `exit <code>` and the last line of the case's
`iter.log`). A case is finished when its `iter.log` ends in `== done`. To count the finished 3-D cases of wave B:

```bash
grep -ho "ice-[A-Z0-9]*_E3D-ES2-E[0-9]* exit 0" $KD/refcorr/pool-*.out | sort -u | wc -l
```

A case cut off by the job's time limit keeps its claim link; the next pool on the same list takes it over and starts that case again.

## Where the runs are

| wave | folder (albedo) |
|---|---|
| A, B (2026-11-13 to 11-16) | `/albedo/work/user/alrobi001/kryos-search` (`W`), to move into `$KRYOS_WORK/earth-search` |
| B emulator | `$KRYOS_WORK/earth-search/emulateB` |
| B2 (from 2026-11-16) | `$KRYOS_WORK/earth-search/waveB2` |

`$KRYOS_WORK` on albedo is `/albedo/work/projects/p_forclima/kryos-work`.
