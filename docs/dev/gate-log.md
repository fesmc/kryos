# Phase 0 sign-off on Levante: gate log (R0, R1)

Status 2026-10-24 18:25 CEST: **Gate 0 passes (R1 == R0 bit for bit, both
gates)**, and `cont2 == cont` holds on both stacks. `chunk` vs `cont` is equal
in every prognostic field and in the 1 kyr restart. It differs only in
diagnostics at the restart time t = 500. Decision (Alex): the restart-time
record, and the bundle a chunk writes at its own start, are not compared
(regression-design §1.2, `Compare.compare_chunks`). Under that rule
`chunk == cont` passes. `$G/ref` was made read-only at 18:25.

## Setup

- Remote: `ssh levante` (user a271116). Successive ssh calls land on different
  login nodes, so `pgrep` over ssh is unreliable. Check the logs, not the processes.
- `$G = /work/ba1442/robinson/kryos-gate` (created): `stacks/`, `tools/`, `reports/`.
- configme 0.15.7 and runme 0.5.26 are already in `~/.local/bin`.
- Modules (login profile): intel-oneapi-compilers/2023.2.1, netcdf-c/4.9.3 and
  netcdf-fortran/4.6.2 (openmpi-5.0.8, intel-2021.5.0), python3/2023.01.

## Stack v24 (R0)

```
configme install yelmox:v2.4 -m dkrz_levante -c ifx --dir $G/stacks/v24/yelmox   # log: stacks/v24/install.log
```

| package | ref | commit |
|---|---|---|
| yelmox | v2.4 | 918d996 |
| fesm-utils | v1.4 | d994f08 |
| yelmo | v2.4 | 7bd4316 |
| FastHydrology | v0.4 | a0d5e8e |
| elsa | v3.1.0 | 74d426d |
| tracer | v0.2 | 208ac3a |
| FastIsostasy | v0.3-fesmc | ce0d6ff |
| chion | v0.2.0 | 0de32ba |
| insol | v1.2.0 | 42b254e |
| rembo1 | v1.3 | 9c24c42 (not linked: yelmox default `rembo=0`, stub) |

`$G/stacks/<stack>/stack.lock` (read-only, written 18:05) holds the full commit of
each checkout, all of them clean, plus the sha256 of `yelmox.x`. Over ssh, run
git as `bash -lc` (it comes from a module; the non-login `.bashrc` fails on `module`).

## Stack S (R1)

- The `configme install yelmox:v2.4+fesm-utils:phase0-kryos` form is a literal
  list (only yelmox and fesm-utils, no expansion; see `configme install --help`).
  So the stack was installed like v24 (log: `stacks/S/install.log`), and then
  inside `$G/stacks/S/yelmox/fesm-utils` the following ran:
  `git checkout -- SHTns; git fetch origin; git checkout -b phase0-kryos --track origin/phase0-kryos; configme config fesm-utils -m dkrz_levante -c ifx`.
- fesm-utils commit **5a31a7c** (`5a31a7cc893d67da49da8a4a42f7db13b963ec31`) is
  the tip at 17:1x. It is newer than 78f79fe by 38e0067 (public `coupler_method`)
  and 5a31a7c (multi-template `nml_validate`; its message says "17 quoted units
  keys differ vs v1.4, metadata only"). It is 29 commits ahead of v1.4.
- The vendored fftw/lis/SHTns are unchanged between v1.4 and phase0-kryos.
- All other packages are at the v24 commits above.

## Gate build

- Variable check: `DFLAGS_NODEBUG` is right for every Fortran build. The machine
  fragment sets `-Ofast -march=znver3 -traceback`, and the command-line value
  overrides it and passes through `$(MAKE) -C` to fesm-utils, yelmo, FastIsostasy,
  chion and insol.
- The fesm-utils external libraries do **not** use `DFLAGS_NODEBUG`. Their flags
  are `CFLAGS_LIB`/`FFLAGS_LIB` (`config/libs.mk`, hard-coded `-Ofast -march=znver3
  -traceback -fPIC`). They were rebuilt with those variables overridden.
  - LIS prepends `-O3 -fomit-frame-pointer` (last `-O2` wins).
  - SHTns keeps its own scheme and is not linked by yelmox.
  - Only the omp variant was rebuilt. The serial externals are still the -Ofast
    builds from the install.
- Script: `$G/tools/build_stack.sh <stackdir>` (`SKIP_LIBS=1` skips the
  externals). It runs:
  1. `make -j16 fftw|lis openmp=1 CFLAGS_LIB=... FFLAGS_LIB=...`
  2. `make clean` in every package
  3. `make yelmox openmp=1 DFLAGS_NODEBUG="-O2 -fp-model precise -march=znver3 -traceback"`
- v24:
  - The first build finished at 17:19 (`stacks/v24/build-gate.log`).
  - Because of the node issue above, a second `make yelmox` was started in the
    same tree by mistake and ran at the same time (`build-gate-2.log`).
  - So a clean Fortran rebuild (`SKIP_LIBS=1`, `build-gate-fortran.log`) was
    started and was still compiling.
  - Externals: fftw and lis built OK (17:14), before the overlap.
- S: `stacks/S/build-gate.log` was still compiling (past the externals and the cleans).
- Every checked compile line carries `-O2 -fp-model precise -march=znver3 -traceback`.
- Outcome:
  - v24: `build-gate-fortran.log` finished at 17:24 with no errors, and all
    168 compile lines carry the gate flags. `yelmox.x` sha256 c0b49bf….
  - S: `make yelmox` finished at 17:22 ("yelmox.x is ready"). After that the
    log ends in `line 29: lean: command not found` / syntax errors. Cause:
    `build_stack.sh` was edited at 17:19:38 while this run was in `make yelmox`,
    and bash then resumed reading the file at the old byte offset. Everything up
    to and including `make yelmox` ran from the old script, so the binary is
    valid; only the final `ls`/`done` lines were lost. 245 of the 247 compile
    lines carry the gate flags; the other 2 are LIS test programs (`test2f`,
    `test6f`, not linked). `yelmox.x` sha256 9350f7f….
  - Lesson: never edit a script on Levante while a run of it is still going.

## Data links and runme

- In each stack root: `ice_data -> /work/ba1442/ice_data`,
  `fesmdata -> /work/ba1442/fesmdata`, `isostasy_data -> /work/ba1442/isostasy_data`.
- `.runme/config.toml` (gitignored) had `account = "CHANGEME"`. It was set to `ba1442`.

## Maps

- yelmox reads `maps/` in the run folder, which runme links to the stack's `maps/`.
  Only `grid_*.txt` are tracked there.
- Plan: copy `v24/yelmox/maps` to `$G/maps/S`, then run one short R0 `mapgen` run
  per gate (`submit_ref.sh v24 <gate> mapgen`, time_end 10, maps -> `$G/maps/S`)
  so that fesm-utils v1.4 generates the maps. Then `chmod -R a-w $G/maps/S`, then
  submit everything with maps -> the frozen folder.
- Done: `$G/maps/S` was seeded from `v24/yelmox/maps` (grid descriptions,
  `genmap.sh`, `remap.sh`). mapgen jobs 28398532 (grl) and 28398535 (ant) ran
  on v24 and COMPLETED with empty `out.err`. They generated `scrip-bilin_GLOBAL-1DEG_GRL-16KM`
  and the ANT-16KM<->ANT-32KM `bilin`/`con`/`nn` maps. `chmod -R a-w` was applied at 18:12.
  If any later run needs a map that does not exist yet, it fails on write instead
  of generating that map silently.
- Existing Levante caches with ANT-16KM<->ANT-32KM maps: `yelmox-v24-test/maps`
  and `yelmox-chion-*`, both of unknown generator version.

## Tools on Levante

- `$G/tools/submit_ref.sh <stack> <gate> cont|cont2|chunk|mapgen [MAPS]`
  - Stages with `runme -s -q compute -w 00:30:00 --omp 16 -e yelmox -n <base> -o <dir> -p ...`.
  - Repoints `<dir>/maps`, then runs `sbatch --parsable`.
  - chunk: c00 runs 0–500; c01 runs 500–1000 with `ctrl.restart=c00/restart-0.500-kyr` and afterok, as in `scripts/run_chunks.sh`.
  - Output goes to `$G/ref/<gate>/yx-918d996/<stack>/<run>`.
- `$G/tools/ncbitcmp.py A B [--tmin T] [--tmax T] [--files GLOB]`
  - Compares each `.nc` recursively, bit for bit (uint view).
  - Labels value-equal but bit-different results as `bits`.
  - Aligns records on common times; reports the first differing time, nrec and max abs diff.
  - Tested on synthetic files.

## Runs (submitted 18:12)

| gate | stack | cont | cont2 | chunk c00, c01 |
|---|---|---|---|---|
| grl | v24 | 28398583 | 28398586 | 28398587, 28398588 |
| grl | S | 28398590 | 28398591 | 28398592, 28398593 |
| ant | v24 | 28398595 | 28398596 | 28398597, 28398598 |
| ant | S | 28398599 | 28398600 | 28398601, 28398602 |

All 16 jobs COMPLETED with exit code 0 (cont about 2.5 min, chunk about 2 min each).

## Results

`$G/tools/compare_R0R1.sh` writes to `$G/reports/R0R1/` (one `.txt` per pair, plus `summary.txt`).

| pair | grl | ant |
|---|---|---|
| R1 vs R0, cont | IDENTICAL (27 files, 1301 vars) | IDENTICAL (27 files, 1395 vars) |
| R1 vs R0, chunk (c00 + c01) | IDENTICAL (42 files) | IDENTICAL (44 files) |
| cont2 vs cont, v24 and S | IDENTICAL | IDENTICAL |
| c00 vs cont, t ≤ 500 | all vars equal | all vars equal |
| c01 vs cont, t ≥ 500 | 35 vars differ, all at t = 500 | 34 vars differ, all at t = 500 |

- The "files missing" in the chunk pairs are restart folders that exist only
  on one side (`restart-0.000-kyr` and `restart-1.000-kyr`). This is expected.
- c01 vs cont: the results are identical in v24 and S, so this is not caused by fesm-utils. The
  differences are:
  - `restart-0.500-kyr/yelmo_restart.nc`: 33 diagnostic variables (`pc_tau*`,
    `visc_bar`, `z_bed_filt`, `bmb_ref`, `cmb_flt_*`, `calv_rate_*`, `uz_*`,
    `T_prime`, `bc_b`, …; 9 of them sign-of-zero only). This is the restart that c01 writes
    at its own start, before the first step, so those fields are not yet recomputed.
  - `yelmo.nc`: `visc_bar` at t = 500 only (1 of 6 records).
  - `yelmo_ts*.nc`: `mb_clip_tot` at t = 500 only (1 of 51 records).
  - Every record after t = 500, and `restart-1.000-kyr`, are bit-identical.
- So yelmox v2.4 is restart-exact in its prognostic state. The record at the
  restart time differs in diagnostics that are not restored.
- `compare_chunks(cont, chunk; gate)` run locally on copies of the v24 runs:
  grl PASS (1104 vars), ant PASS (1364 vars). Negative checks on grl: a change
  at t = 500 in c01 still passes, a change at t = 700 fails with first
  divergence 700 `c01/yelmo.nc:H_ice`, and a missing c01 bundle at 1000 is an error.
  The S chunk runs are bit-identical to the v24 ones (R1 vs R0, chunk).

## Stack S2 and baseline R2 (2026-10-25)

Decision (Alex, 2026-10-24): Yelmo is fed a prepared topography from the hub
(regression-design §1.3 item 4). Commits, brought to Levante as git bundles
(`$G/bundles`, nothing pushed):

| repo | branch | commit | base |
|---|---|---|---|
| fesm-utils | `topo-prep` | 15746a9 | 5a31a7c (`phase0-kryos`) |
| yelmo | `topo-prepared` | b5ab5d2c | v2.4 (7bd4316) |
| yelmox | `topo-prep` | 74f3a2b | v2.4 (918d996) |

- `$G/stacks/S2` is a copy of S with these checked out (fesm-utils reconfigured),
  built with `SKIP_LIBS=1` (S's gate-built fftw/lis). Scripts: `tools/setup_S2.sh`
  (its submissions failed: `$G/ref` was read-only; parent folders reopened),
  `tools/submit_S2.sh`, `tools/compare_S2.sh` (reports in `$G/reports/S2/`).
- Over ssh, start long scripts as `setsid nohup bash -l <absolute path> ... < /dev/null &`:
  a plain `nohup ... &` dies with the session, and `bash -l` starts in `$HOME`.

Results (all jobs COMPLETED, exit 0):

| pair | grl | ant |
|---|---|---|
| R1(S2) vs R0, cont (yelmox v2.4 on S2: new code off) | IDENTICAL | IDENTICAL |
| R2 cont2 vs cont | IDENTICAL | IDENTICAL |
| R2 c00 vs cont (t ≤ 500) | equal | equal |
| R2 c01 vs cont (t ≥ 500) | differs only at t = 500 (same diagnostics as R0) | same |
| R2 vs R0, cont | 788 of 1301 vars differ | 943 of 1395 vars differ |

- So the new fesm-utils and yelmo code is neutral when off; R2 is deterministic
  and restart-exact under the chunk rule: R2 is a valid baseline. Refs frozen.
- The hub preparation in R2 (`out.out`): lakes removed; ice < 1 m removed (grl
  225 cells, ant 293 on the 16 km hubs); `grad_lim_zb = 0.5` needed 0 passes;
  Yelmo logs "prepared ... cleaning skipped" on both paths.
- R2 vs R0 at t = 0 (`yelmo.nc`, after the cold start):
  - grl: `H_ice_ref` identical (the hub's preparation equals Yelmo's
    present-day path on the same grid); initial `H_ice` differs in 636 cells,
    21 by > 10 m, max 356 m: R0 started from Yelmo's other (initial-state)
    path (thin ice 0.1–10 m set to 10 m; lakes with `bnd%z_sl`), R2 from the
    reference.
  - ant: `H_ice_ref` differs in 2456 cells, 657 by > 10 m, max 272 m: the
    preparation runs on the 16 km hub and is then remapped (`con`) to 32 km,
    where R0's Yelmo prepared on 32 km (D3, accepted).
  - Through 1 kyr: `V_ice` differs by at most 0.011 (grl) and 0.084 (ant)
    [Yelmo units, 1e6 km3]; `z_bed` from t = 10 via isostasy.
- Gate 1 is now `K == R2` (refs `$G/ref/<gate>/yx-74f3a2b/S2/`).

## kryos on Levante; C6 check against R2 (2026-10-25)

- `$G/kryos`: clone of kryos main (bundle), dependency slots linked to the S2
  checkouts, configured with configme 0.16 from source
  (`PYTHONPATH=$G/tools/configme-src/src python3 -m configme config kryos --only
  -m dkrz_levante -c ifx`, so S2's packages are not reconfigured), built with
  the gate flags (`tools/setup_kryos.sh`, log `stacks/kryos-build.log`). It
  compiles cleanly against S2's yelmo (v2.4 + `prepared`).
- C6 dump (`test_domain.x dump`, par files `bundles/dump_{grl,ant}.nml`, run dirs
  `$G/kryos-dump/{grl,ant}`; run dirs need links `input fesmdata ice_data
  isostasy_data fesm-utils maps`), compared with R2 (`tools/cmp_c6.py`):

| field | grl | ant |
|---|---|---|
| refs(ice) `z_bed` vs Yelmo `z_bed_ref` (R2 restart t0) | identical | identical |
| refs(ice) `z_bed_sd` vs Yelmo `z_bed_sd` | identical | identical |
| refs(ice) `H_ice` vs Yelmo `H_ice_ref` | 1 cell (1.3 m) | 2 cells (0.5 m) |
| hub `z_bed` at t0 vs R2 `htopo.nc` t0 | identical | identical |

- The `H_ice` cells are where Yelmo's `mask_ice` removes reference ice (expected).
  Hub `H_ice`, `z_srf` at t0 differ from R2's `htopo.nc` t0, which is written
  after the cold start (the hub then follows the ice; Antarctica after its 10 yr
  equilibration): that comparison belongs to C9.

## Gate 1, Greenland: K == R2 (2026-10-25)

- Run: kryos `a8666f2` on stack S2 (`tools/run_kryos.sh grl cont`, job 28411222,
  `$G/runs/grl/kx-a8666f2/S2/cont`), 1000 yr in 2:51 (yelmox R2 2:34), par file
  `par/gate/grl.nml` (converted, `tools/bin/yelmox2kryos.jl`).
- Compare (`tools/bin/kryos-compare.jl R2 K --map tools/gate/map_phase1.toml
  --gate grl`, run locally on copies): **PASS**, 202 variables, 0 failing, first
  divergence none, 12 kryos-only variables (N). Report:
  `docs/dev/gate-reports/grl-K-a8666f2-vs-R2.md`.
- On the way:
  - wiring rule changed (`9fadf2f`): rows never round; sp kernels narrow
    their operands; consumers narrow once (the pr / bmb_shlf / Q_geo round trips
    were not exact before);
  - ifx: five compiler bugs worked around (`40f54fa`..`a8666f2`): the output
    table cache, host-associated arrays in tests duplicated under `-qopenmp`,
    assumed-rank in a nested SELECT RANK, and `mold=` from a non-polymorphic
    actual (the `*_adopt` wrappers now take an allocatable prototype); the whole
    test suite passes under ifx with the gate flags;
  - gate map: the surface, shelf and climate adapters write yelmox's own
    fields and units (class E), the default writer writes exports as float
    (C rows wrapped in `Float32`).
- Not yet: restart bundles (C10), chunk vs cont, Antarctica (C11).

## Gate 1, Antarctica: K == R2 (2026-10-25)

- Run: kryos `2461747` (= `a8666f2` + the converter fix below) on stack S2,
  job 28411639, `$G/runs/ant/kx-2461747/S2/cont`, 1000 yr in 3:23.
- Compare: **PASS**, 350 variables, 0 failing, first divergence none
  (`docs/dev/gate-reports/ant-K-2461747-vs-R2.md`). This is the first run with
  hub (ANT-16KM) and ice (ANT-32KM) on different grids, the 10 yr
  `init_equil_time`, the region files (APIS, EAIS, WAIS), the ocean profiles to
  the shelf (bilin 32 -> 16 km) and `bmb_shlf` back to the ice (con 16 -> 32
  km): the predicted rare tie differences of `bmb_shlf` (adapters 1.2) did not
  occur.
- Fixes on the way: the converter forced `remap = "bilin"` on the geothermal
  and sediment fields, so datasrc asked for an ANT-32KM -> ANT-32KM map (into the
  read-only maps folder); now `auto` (a copy on the instance grid, as yelmox's
  short-circuited bilin). `Compare.jl`: a row naming variables pairs the files
  only if one of them exists (the Antarctica gate has no forcing coefficients).
- Gate 1 (continuous runs) passes for both gates. Open: restart bundles and
  chunk vs cont (C10, in progress).

## Phase 1 gate complete: restarts and chunks (2026-10-25)

- kryos `a2cf248` (C10: restart bundles, manifest, restore) on stack S2; runs
  staged by `tools/run_kryos.sh <gate> cont|chunk` (chunk: c00 0->500, c01
  500->1000 from `c00/restart-500`, afterok; overrides set in the par file by
  `tools/set_nml.py`). All six jobs COMPLETED (cont 2:55 / 3:33, chunks
  2:14+2:26 / 2:22+2:15).

| comparison | grl | ant |
|---|---|---|
| K cont vs R2, outputs + restart bundles 500, 1000 (`--restarts`) | PASS, 924 vars | PASS, 1036 vars |
| K chunk vs K cont (`--chunks`, the §1.2 rule) | PASS, 1470 vars | PASS, 1710 vars |

- Yelmo's `pc_tau_max` (running maximum over 50 internal steps, not
  restored) is refilled within the 500 yr chunk: the bundle at 1000 is
  identical. Reports in `docs/dev/gate-reports/*-a2cf248-*`.
- With Gate 1 (continuous runs, above), kryos phase 1 reproduces the yelmox
  baseline R2 bit for bit for both gates, including restart exactness.
- Found on the way (C10, local test): with Yelmo's default `vm-m16` calving a
  Yelmo restart is not exact (diverges in the first step after a restore);
  the gate configs do not use it. Yelmo issue, not kryos.

## Base length: 15 kyr (2026-10-25)

- R2 (yelmox `74f3a2b`) and K (kryos `a2cf248`) on stack S2, the gate `-p`
  with `ctrl.time_end=15000` (`submit_ref.sh S2 <gate> base`, `run_kryos.sh
  <gate> base`, 2 h limit): grl R2 19:10, K 22:34; ant R2 34:24, K 38:08.
- Compare with `--restarts` (bundles every 500 yr): **grl PASS, 924 variables;
  ant PASS, 1036 variables**; first divergence none in both
  (`docs/dev/gate-reports/*-15kyr.md`). 1500 coupling steps bit for bit.
- The Antarctica comparison ran on Levante: KryosTools with juliaup Julia 1.13
  (`JULIA_DEPOT_PATH=$G/julia-depot`, FesmUtils.jl cloned to `$G/FesmUtils.jl`
  so the tools' relative source path resolves); Greenland locally on copies.
- Phase 1 regression gate: complete (1 kyr cont, restarts, chunks, 15 kyr).

## Performance: stack S3 (2026-10-25)

- fesm-utils `phase0-kryos` `bb39f2f`: ncio attribute routines take an
  `ncid`; `fields_restart_write`/`read` and the clock/cadence restart writers
  open each file once; `put_copy` tests finiteness with `abs(x) <= huge(x)`.
  kryos `af8f1b9`: one-grid rows put without a copy, lapse skips `exp` where the
  exponent is zero, backends declare their restart files (no `ls`).
- Stack S3 = S2 with that fesm-utils, every Fortran package rebuilt (an ncio
  interface change needs all dependents recompiled: a stale yelmo/fastisostasy
  build segfaulted in `nc_exists_var` locally). `tools/setup_S3.sh`; kryos's
  dependency links and `run_kryos.sh` now point at S3.
- Identity checks (`--restarts`, all PASS): R1(S3) vs R0 (grl 1274, ant 1368
  vars), R2(S3) vs R2 (same), K(S3) vs K(S2) 1 kyr (grl 1249, ant 1341) and
  15 kyr grl (15641). Reports in `$G/reports/S3/`.
- 15 kyr grl wall time (timers, s): kryos 1125 (was 1256), yelmox 1066:
  output 220 (178), geo 310 (309), ice 469 (463), climate 13.5 (9.5),
  surface 96 (93), shelf 15 (13). The remaining output difference is the extra
  kryos files (`_ts` every 10 yr, default-writer 2D every 100 yr): their opens.

## WP0: stack S4, NOR baselines, data check (2026-10-29)

Scope narrowed to the science case (VILMA2 global + NH-32KM + ANT,
`phase2-plan.md` top): S4, its identity checks, NOR only. Reports:
`gate-reports/wp0-baselines.md`, `gate-reports/wp0-data-check.md`.

- yelmox `topo-prep` **5ac89c8**: `[topo]` keys in the other 8 bases
  (Q1, per-base values; PYR `grad_lim_zb = 1.0`, others 0.5). Not pushed
  (bundles in `$G/bundles`).
- Stack **S4** = S3 + yelmo `topo-prepared` ea90fee7 (vm-m16 restart fix)
  + rembo1 v1.3 linked (`rembo=1`); `tools/setup_S4.sh`, binaries
  `stacks/S4/bin/{R1,R2}`. Identity (`--restarts`, all PASS): R1(S4) vs R0,
  R2(S4) vs R2, K(S4) vs K(S3), grl and ant. kryos's links and
  `run_kryos.sh` take `STACK=S4`.
- Maps: `$G/maps/S4` = `maps/S` + 8 pairs (mapgen R0 runs on v24), frozen.
  NOR needs `bilin_GLOBAL-1DEG_NH-32KM`; the others were made before the cut.
- NOR (1 kyr, gate `-p`, ~11 min per run): R1 == R0 (cont, c00), cont2 ==
  cont for R0, R1, R2'; R1 and R2' chunk == cont (vm-m16 fix). R0 chunk vs
  cont and R1 vs R0 c01 differ from t = 510 (v2.4 has no vm-m16 fix;
  expected). Gate `nor` is K == R2' (`$G/ref/nor/yx-5ac89c8/S4`).
- R2' vs R0, NOR t = 0: `H_ice_ref` differs (5436 cells, ≤ 10 m in 5425,
  max 50 m) because recon sets `H_ice_ref = H_ice` after its recipe, which
  starts from Yelmo's initial ice (R0 raises 0.1-10 m ice to 10 m, R2'
  takes the prepared field); initial `H_ice` differs after the recipe's
  200 yr of dynamics. grl/ant (none/equil) keep the prepared reference.
  WP4 reproduces the recipe on the hub's prepared field.
- Open: the hub log reports 8385 cells removed by `H_ice_min` on NH-32KM,
  the raw TOPO-GEBCO2026 file has 1313 cells with 0 < H < 1 m; the hub field
  before preparation is not the raw file (to trace in WP4).
- Data check (Q17): every file of the 8 bases resolves through the gate
  links; no links made. Three yelmox series files are only in the stack's
  `input/` (PYR, GRP, rGRL: deferred bases).
- Julia compares on Levante need `julia +1.13`.

## WP1: experiment plugin mechanism (2026-10-29)

Branch `wp1`, commits as listed in `phase2-plan.md` WP1 (hooks, `prepare`,
`PASS_SPINUP`, `step_time_t%tick`, `init_time_thrm`, experiment lists, restart).
No experiment is active in grl or ant: every list is empty.

- Local, `make check` at every commit; `test_run_sequence` `plain` (cont, c00,
  c01) compared with `kryos-compare --restarts` (outputs and bundles, 1151
  variables) against the previous commit's build: PASS after each of
  PASS_SPINUP, `init_time_thrm` (off), `prepare`, `tick`, the config
  templates, the step hooks, the cold-start hooks and the restart commit. A
  cont run with `init_equil_time = 5` (the equilibration, which the test does
  not exercise) is identical to the build before the refactor.
- Guard: a host experiment that resets Yelmo's recorded tick at `ICE_BND`
  makes `update` stop with `kryos_ice_yelmo::update ... was not prepared`
  (`test_run_sequence`, child process).
- Levante, kryos `f81e004` (wp1 as a bundle, never pushed) on stack S4,
  `par/gate/{grl,ant}.nml` regenerated from yelmox `topo-prep` 5ac89c8 (new keys
  `experiments = ""`, `init_time_thrm = 0.0`, `phys_const_file`; same values
  otherwise). Jobs (`run_kryos_wp1.sh`, COMPLETED, 2-3.5 min): grl cont 28537662,
  cont2 28537678, chunk 28537684 + 28537686; ant cont 28537688, cont2 28537689,
  chunk 28537691 + 28537692. Compares with `compare_wp1.sh` (`julia +1.13`),
  reports in `$G/reports/WP1`, `gate-reports/*-K-f81e004-*`:

| comparison | grl | ant |
|---|---|---|
| K cont vs R2 (yx-74f3a2b S2; mapped, `--restarts`) | PASS, 924 vars | PASS, 1036 vars |
| K cont2 vs R2 (same) | PASS, 924 vars | PASS, 1036 vars |
| K cont2 vs K cont (`--restarts`) | PASS, 1249 vars | PASS, 1341 vars |
| K chunk vs K cont (`--chunks`) | PASS, 1470 vars | PASS, 1710 vars |

  First divergence none in all eight; 12 kryos-only variables (N) in the
  mapped compares, as in phase 1. K chunk == R2 follows from the last two rows.
- The kryos checkout on Levante had an untracked `input/yelmo_phys_const.nml`
  (now a tracked synchronised copy); moved to
  `$G/yelmo_phys_const.nml.untracked-in-K`.

## Phase 4 WP0: stack S5, G0 (2026-10-29)

Report `gate-reports/p4-wp0.md`. S5 = S4 + fesm-utils phase0-kryos 5662899 + SHTns + vilma2 main 5e5a9ee, all
Fortran packages rebuilt (`tools/setup_S5.sh`, `tools/gate/stack-S5.toml`). libvilma.a links on Levante (vilma2 host tests
PASS, job 28542741). G0 (`--restarts`): K(S5, vilma=0) == K(S4) and R2(S5) == R2(S4), grl and ant, all PASS
(jobs 28542775/77, 28542778/79). `$G/maps/S5` = S4 copy; Gauss and GLOBAL-0.5DEG maps pending (WP2/WP3 code).
Data check: PaleoMIST 0.5DEG (0 ka = last record) and both Bagge files present; `fesmdata-prepared` is not a gate link.

## Phase 4 WP1: contracts (2026-10-29)

Report `gate-reports/p4-wp1.md`. Branch `p4-wp1` on main 79461c5: `dphi` (fastisostasy `dphi = bsl`), ice `footprint`,
`H_ice_pd`, the `flotation` binding (Yelmo's own rule; == `tpo%f_grnd` bit for bit on grl), geo `host`/`port` options with
the rank-2 global declaration, `T_GIA` rows (extra sources in `wire_t%srcs`) and the `wire_resolve` rule for core-filled
imports. Mac, `make check` at every commit; `test_run_sequence` `plain` (cont, c00, c01) compared with
`kryos-compare --restarts` against the base build: PASS at every commit (1151 / 620 / 620 variables, 0 failing, 0 absent),
head against a fresh 79461c5 build too. Levante G1 (ant, nor, `K3-sci`) is WP7.

## Phase 4 WP2: the vilma2 geo backend (2026-10-29)

Report `gate-reports/p4-wp2.md`. Branch `p4-wp2` on main 21320db: `vilma=1` build (libvilma, SHTns, `vilma-static`, sync copies of
`vilma_defaults.nml` and `vilma-variables.md`, `vilma2:main` pin), `vilma_earth_t` (global, Gauss grid `GAUSS-L<lmax>`) and `vilma_geo_t` (port),
stub, factory, `[geo_vilma2]` template, `gia_check` (run-level contract checks). Mac, `make check` with `vilma=1` and `vilma=0` at every commit;
the adapter equals `solid_earth_*` called directly, bit for bit (reference, seed, spin-up, steps, restart round trip). The `test_run_sequence`
plain runs of both builds compared with `kryos-compare --restarts` against 21320db: PASS at every commit (1159 / 625 / 625 variables, 0 failing,
0 absent). Open: planet radius 6371221 vs VILMA2's 6371000 m (tolerance 1e-4 chosen), fixed group name `vilma`, port cadence in WP4.

## Phase 2 WP2-WP4: experiments relax, opt cb_ref, opt tf, recon (2026-10-29)

Branch `p2-exp` (design in `phase2-plan.md`, "WP2-WP4 approved design"). Mac: `make check` (`yelmo=1`; `yelmo=0` too, stubs) and
the `plain` run of `test_run_sequence` (cont, c00, c01) compared with `kryos-compare --restarts` against a fresh build of the
base commit at every step: PASS (0 failing, 0 absent); the plugin sets (`xr/`) check chunk == cont with the windows straddling the restart.

### `ant-cf` (W2d, kryos `aeee5dd`, stack S5, kryos built with `vilma=0`)

R2: yelmox topo-prep (S5 checkout 5ac89c8, code-identical to 74f3a2b), `submit_ref.sh S5 ant-cf cont`, shipped ANT with
`opt.opt_tf=none` (relax + opt_cf on), `ref-S5/ant-cf/yx-5ac89c8/S5/cont`, job 28545775. K: `run_kryos_p2.sh ant-cf cont|cont2|chunk`
(jobs 28545844, 28545861, 28545886/89), compared by `compare_p2.sh` (job 28547483, `julia +1.13`, outputs and restart bundles):

| comparison | result |
|---|---|
| K cont vs R2 (mapped, `--restarts`) | PASS, 1036 variables, 0 failing, 0 absent, 14 new (N) |
| K cont2 vs R2 | PASS, 1036 variables |
| K cont2 vs K cont | PASS, 1350 variables |
| K chunk vs K cont (`--chunks`) | PASS, 1722 variables |

R2 chunk is not compared (yelmox restarts the window time per chunk). Relax and the cb_ref optimisation are bit-identical to yelmox.

### Hub count 8385 vs 1313 (W4a, design Q8; closes the open item of the WP0 NOR notes)

`tools/levante/wp4_hub_count.py` (offline, numpy, float32 in the order of `htopo_init`) on `NH-32KM_TOPO-GEBCO2026.nc`:

| step | 0 < H < 1 | H < 0 | htopo count (H < 1, H /= 0) |
|---|---:|---:|---:|
| raw file | 1313 | 0 | 1313 |
| after `remove_lakes` (the hub's order) | 1009 | 7376 | **8385** |

Both hypotheses hold: `remove_lakes` runs before `H_ice_min` and rewrites 11423 cells (2600 lakes filled, 8823 erroneously
thick reduced to `H = z_srf - z_bed`; 7376 of those are negative where `z_srf < z_bed`), and the logged count includes the
negatives (`H < H_ice_min .and. H /= 0`). The hub field is the raw file on the same grid (NH-32KM, `remap = "con"`
identity). Nothing to change: the 8385 is the log's counting rule, not extra ice removed.

### `ant-full` (W3d, kryos `5a663a3`, S5)

Shipped ANT: relax, opt_cf and opt_tf = L21-points (the shelf corrections `tf_corr`, `tf_corr_basin` now the shelf component's, restart
file `shelf_ant_corr_restart.nc`; `corr_method = tf`, four basin constants). R2: `ref-S5/ant-full/yx-5ac89c8/S5/cont`, job 28548078.
K cont, cont2, chunk: jobs 28548453, 28548455, 28548456/57; `compare_p2.sh` job 28548977 (map rows send yelmox's `marine_shelf.nc:tf_corr`
and `tf_corr_basin` to `<shelf>_corr_restart.nc`):

| comparison | result |
|---|---|
| K cont vs R2 (mapped, `--restarts`) | PASS, 1042 variables, 0 failing, 0 absent, 16 new (N); `tf_corr` and `tf_corr_basin` of both bundles pass (145161 cells) |
| K cont2 vs R2 | PASS, 1042 variables |
| K cont2 vs K cont | PASS, 1362 variables |
| K chunk vs K cont | PASS, 1736 variables |

The thermal-forcing optimisation is bit-identical to yelmox: the ice exports (`H_ice`, `H_grnd`, `dHidt`, `H_ice_pd`, `H_grnd_pd`)
remapped ANT-32KM -> ANT-16KM (bilinear) through the domain coupler give yelmox's `optimize_tf_corr` inputs bit for bit.

### `grl` and `ant` after the shelf move (kryos `5a663a3`, S5)

`par/gate/{grl,ant}.nml` regenerated (the marine_shelf keys moved to `&shelf_corr_<sfx>`, `corr` on the shelf instance), K against
the S5 R2 runs of G0 (`ref-S5/{grl,ant}/yx-5ac89c8/S5/cont`), same comparisons (K jobs 28549571-77 and 28549580-88, `compare_p2.sh` jobs
28550099 and 28550100):

| comparison | grl | ant |
|---|---|---|
| K cont vs R2 (`--restarts`) | PASS, 930 variables | PASS, 1042 variables |
| K cont2 vs R2 | PASS, 930 | PASS, 1042 |
| K cont2 vs K cont | PASS, 1270 | PASS, 1362 |
| K chunk vs K cont | PASS, 1496 | PASS, 1736 |

Bit-identical in all outputs and bundles: the move of the corrections to the shelf component (and the new `H_grnd`, `H_grnd_pd`
exports) changes no value.

### `nor` (W4e, kryos `888884d`, S5)

NOR (NH-32KM): cold start from the ICE-6G_C reconstruction (`recon_regions = region:North_America`, `recon_time = -26e3`), R2' =
yelmox topo-prep: `ref-S5/nor/yx-5ac89c8/S5/cont`, job 28549546. K cont, cont2, chunk: jobs 28549631, 28549632, 28549633/34;
`compare_p2.sh` job 28550508:

| comparison | result |
|---|---|
| K cont vs R2' (mapped, `--restarts`) | PASS, 911 variables, 0 failing, 0 absent, 16 new (N) |
| K cont2 vs R2' | PASS, 911 variables |
| K cont2 vs K cont | PASS, 1251 variables |
| K chunk vs K cont | PASS, 1472 variables |

The recon recipe (reference ice from the file, regions, smoothing, `yelmo_init_topo`, the thin-floating cleaning, `H_ice_ref = H_ice`,
`yelmo_ice_publish`, `domain_refresh_boundary`, the 200 yr stabilisation) reproduces yelmox bit for bit, including the hub and the
surface refreshed inside the recipe. R2' vs R0 differences are those documented for WP0.

## Phase 4 WP3: kryos_globe_ice (2026-10-29)

Report `gate-reports/p4-wp3.md`. Branch `p4-wp3` on main a6e91bf: the gather (§4) and the reference (§5.1-5.3) with one reference definition (D18),
owner weights for overlapping domains, the hub output `globe_ice.nc` (rounded, lossless under the ncio switch) with its conservation check, and the
offline merge in Julia (`KryosTools.GlobeIce`). Mac, `make check` at every commit; `test_run_sequence` plain runs compared with
`kryos-compare --restarts` against a build of the base: PASS at every commit; head against a fresh a6e91bf build: 1171 / 632 / 632 variables, 0 failing, 0 absent. G1b locally: hub closure 1e-15,
order swap identical, offline merge equal to 1e-12. The run with and without `[globe_ice]` and the real-run offline merge are WP4/WP7.

## Phase 3 WP1, WP5: two domains, local (2026-10-29)

Branch `p3-dom` on main e6b44ed. Merge step `tools/bin/kryos-merge.jl` (`Convert.merge_par`): `par/kryos_NH-ANT.nml` (from `par/kryos_North.nml`, new, and
`par/kryos_Antarctica.nml`) and `par/gate/sci.nml` (from the `nor` and `ant` gate files); 46 groups renamed and 46 keys repointed in each; the run-wide keys
that differ are resolved explicitly (recorded in the file header): `clock.time_const` (-21e3 | 0.0) = 0.0, `barysealevel.sl_path` (read only by bsl `file`) =
"none", for the base file also `clock.time_end`. Gate tooling: `tools/bin/kryos-gate-kk.jl`, `tools/gate/sci.toml`, `sci-const.toml`.
`tests/test_run_domains.f90` (Mac, GRL-16KM + ANT-32KM, synthetic snapshots, 80 yr, bundles at 40 and 80, files from the merge step), all checks through
`kryos-gate-kk.jl`:

| check | result |
|---|---|
| bsl const: each domain of the two-domain run vs its single-domain run (every file, bundles at 40 and 80, globe files; skipped: `forcing_imports.nc:M_ice`, `geo_glob_imports.nc:dV_bsl` of the bundles) | PASS, Greenland 1144 variables, Antarctica 1097 |
| bsl fastiso: domains swapped in `[kryos] domains` vs the original order (files and bundles) | PASS |
| bsl fastiso: chunk 0 -> 40 -> 80 vs continuous (files after 40, bundles; `pc_tau_max` skipped as in `test_run_sequence`) | PASS |
| bsl budget `bsl(t) - bsl(t-dt) = -sum_d dV_bsl_d / A_ocn` recomputed from the exports | closes: max residual / scale 2.9e-8 (0.25 eps of sp) |
| recon (Greenland) + relax (Antarctica), const: each domain vs its single-domain run; chunk vs continuous | PASS |
| derived `dTo` produced by both climates | stops at init: `run_check_producers` (new; `forcing_set`'s check was never reached) names the coefficient |

`make check` passes; the `test_run_sequence` runs (cont, chunk c01) compared with `kryos-compare --restarts` against a build of main e6b44ed: PASS (1171 / 632
variables, 0 failing). Levante gates: below.

### Levante: `sci-const`, `sci` (kryos 99589d0, stack S5b)

Code: branch `p3-dom` at 99589d0 (bundle `$G/bundles/p3dom.bundle`; checkout and build `$G/p3dom/kryos`, linked to S5b, gate flags, `vilma=0`), maps `$G/maps/S5`.
Runs staged by `tools/levante/run_kryos_kk.sh` (par `par/gate/{sci,nor,ant}.nml` with the `kp` overrides of the gate files), checks by
`tools/levante/compare_kk.sh` (`kryos-gate-kk.jl check`, `julia +1.13`), reports `$G/reports/KK/{sci,sci-const}-99589d0-S5b.md`. 1 kyr, gate `-p`
(`write_lossless`, 2D every 100, 1D and 2Dsm every 10, bundles at 500 and 1000). Two-domain runs 13.7 min (16 threads), NOR alone 11.8, ANT alone 2.6.

`sci-const` (bsl const; the single-domain references are the `nor` and `ant` gate files on the same binary, `nor` with the merged clock `time_const = 0.0`);
jobs cont 28554086, nor 28554093, ant 28554097:

| check | result |
|---|---|
| North: two-domain run vs `nor` (every instance file, bundles at 500 and 1000, globe files) | PASS, 1249 variables, 0 failing, 2 skipped (aggregates: `forcing_imports.nc:M_ice`, `geo_glob_imports.nc:dV_bsl`) |
| Antarctica: two-domain run vs `ant` | PASS, 1360 variables, 0 failing, 2 skipped |

`sci` (bsl fastiso); jobs cont 28554103, cont2 28554104, chunk 28554105/06, swap 28554108:

| check | result |
|---|---|
| cont2 vs cont (`--restarts`) | PASS, 2575 variables |
| chunk vs cont (c00 0-500, c01 500-1000 from `restart-500`, bundles) | PASS, 3166 variables (no skips: `pc_tau_max` is refilled in a 500-yr chunk) |
| `[kryos] domains` swapped vs original order (files and bundles) | PASS, 2575 variables |
| `bsl` budget (`geo_glob_ts.nc` bsl, A_ocn; `geo_nor_ts.nc`, `geo_ant_ts.nc` dV_bsl) | closes: 100 records, max residual / scale 2.2e-8 (0.19 eps of sp), max residual 3.8e-6 m on a bsl that moves to -90.8 m; dV_bsl of both domains about 1.0e13 to 1.5e13 m3 per step at the end |

The budget bound stays 4 eps of sp (4.8e-7); the two runs (local 2.9e-8, Levante 2.2e-8) sit more than a factor 10 below it.

Clock choice (`clock.time_const = 0.0` in the merged file, NH has -21e3): NH alone, bsl fastiso, `time_const = 0.0` (`sci-const` run `nor-t0`, job 28554099) against the
`nor` gate's K cont (kryos 888884d, S5, -21e3), identity `--restarts` (`$G/reports/KK/nor-t0-vs-kx-888884d_S5_cont.md`): **1247 of 1251 variables equal bit for bit,
4 differ**: `time_init` and `time_prev` of `forcing_alpha_restart.nc` in both bundles (-21000 vs 0, the forcing time the tsgen coefficient stores; nothing
else, no output file and no field of any instance differs). Not bit for bit by the letter: reported to the coordinator.

After the rebase onto main 1f8a7b8 (p3-out, p3-smb merged): `make check` and `make check-long` pass; the long test found a latent crash in the climate init of
a domain without terms (`size` of an unallocated coefficient table; the Antarctica single run segfaulted), fixed in f4866ae (the table is empty, no value
changes). The Levante gates above ran on 99589d0 (before the rebase and the fix; the code paths they exercise are unchanged).

Accepted (Alex, 2026-10-29, option a): `clock.time_const = 0.0` in the merged file. The four differing values of the clock check (`time_init`, `time_prev` of
`forcing_alpha_restart.nc` in the bundles at 500 and 1000) are the stored clock constant, which enters no physics: an accepted difference
(`regression-design.md` §1.6). The science configuration will set the clock explicitly ("in the future it will be right": the real transient clock).

### Rerun on the rebased head and registration of `K3-sci` (kryos 2d80b43, stack S5b)

Code: p3-dom rebased on main 1f8a7b8 with the snapshots fix (2d80b43; the later commits of the branch are docs only), built in `$G/p3dom/kryos` on S5b as before.
Jobs: `sci-const` cont 28555723, nor 28555737, ant 28555738; `sci` cont 28555739, cont2 28555740, chunk 28555743/44, swap 28555746 (two-domain runs 13.5 min).
Reports `$G/reports/KK/{sci-const,sci}-2d80b43-S5b.md` (jobs 28556423/24):

| gate | check | result |
|---|---|---|
| `sci-const` | North: two-domain run vs `nor` | PASS, 1251 variables, 0 failing, 2 skipped (aggregates) |
| `sci-const` | Antarctica: two-domain run vs `ant` | PASS, 1362 variables, 0 failing, 2 skipped |
| `sci` | cont2 vs cont (`--restarts`) | PASS, 2579 variables |
| `sci` | chunk vs cont (bundles) | PASS, 3170 variables |
| `sci` | domains swapped vs original order | PASS, 2579 variables |
| `sci` | bsl budget | closes: 100 records, max residual / scale 2.2e-8 (bound 4.8e-7), max residual 3.8e-6 m |

(The variable counts are above the 99589d0 run: main's output-precision changes add variables.)

**`K3-sci` registered**: `$G/ref/K3-sci/kx-2d80b43/S5b/{sci,sci-const}/cont` (`register_k3.sh 2d80b43 S5b`; 122 nc files each, files and bundles at 500 and 1000,
`par.nml`, `kryos_commit.txt`; read-only, 1.4 GB; `registered.txt`). Phase 4 G1 compares against it (`kryos-compare.jl --map identity --restarts`).

## Phase 4 WP7, part 1: stack S6, G0-G4 (2026-10-30)

Report `gate-reports/p4-wp7a.md`. Branch `p4-wp7` on main cd85c14. S6 = S5b + fesm-utils 0a626e9 + vilma2 2676879 (`tools/setup_S6.sh`; kryos cd85c14 built `vilma=1` and `vilma=0`);
identity (K v0 vs K(S5b), R2(S6) vs R2(S5b), grl and ant, `--restarts`): PASS. G1 (grl, ant, ant-full, nor against the earlier K; sci, sci-const against K3-sci): PASS, bit identical.
`maps/S6` frozen (S5 + 15 files: GAUSS-L170 pairs, GLOBAL-0.5DEG -> GAUSS/NH/ANT, ANT-16KM pairs).
G0: `make check` fails on Levante at cd85c14 (ifx): `test_forcing_increment` (name clash), `test_geo_vilma2` (main-program arrays; module fix), `test_globe_ice` and every new con map
(fesm-utils `conservative.f90`: allocatable listed `private` after a serial alloc/free; patch `gate-reports/p4-wp7a-fesm-utils-conservative-ifx.patch`); with the fixes both builds pass.
G1b: with `[globe_ice]` bit identical to K3-sci (2579 variables), cont2 and domain swap identical, `globe_ice.nc` equals the offline merge; the hub closure is 2.6e-10 (reference) and 1.5e-10 (ice), not 1e-12,
and VILMA2's host check rejects `host_w` = 1 + 2.8e-9 at the polar row: runs with VILMA2 need a decision (variant diff `gate-reports/p4-wp7a-kryos-hub-variant.diff`).
G2a (lmax 170): adapter == `vilma.x` to <= 6e-10 of max |field| (step 1) and <= 1e-12 after, with equal spin-up; the default spin-up tolerance stops the two at different passes (21 % at the cap centre).
G3 (400 yr): merge closure 7.1e-13, `mass_resid` 2.1e-15, `bsl_ive` against the gathered ice 4.3e-7 as volumes, `A_ocn_ref` +0.23 %. G4: cont2 == cont and rounding (max error / tolerance 0.977, bundles identical) pass; chunk != cont at round-off (1e-14 `bsl`, 7e-12 m `dz_bed`) in the Earth's state.
Science par `par/kryos_NH-ANT_vilma2.nml` (+ conversion note, generator `tools/bin/kryos-vilma2-par.jl`); `par/kryos_Antarctica.nml`, `par/kryos_NH-ANT.nml` regenerated (they predated the shelf move).

## Phase 4 WP7, part 1b: S7 (2026-10-30)

Report `gate-reports/p4-wp7b.md`. fesm-utils `phase0-kryos` f30d952 (pushed): conservative_planar pre-pass arrays renamed (ifx), `test_conservative_project` (fails on the old code with ifx). kryos `p4-wp7` 7e631a2 on main ffe745c:
producer-side clip of `host_w` and `host_C` in the gather (`HOST_CLIP_TOL` 1e-6, excess logged), `HUB_TOL` 1e-9, test fixes, plan 4.5 and 13.2, vilma-variables.md synchronised with vilma2 d2006c8, `fill_int` in the `vilma` stream,
science par from PaleoMIST at -80 ka in both domains. Stack S7 = S6 + fesm-utils f30d952 + vilma2 d2006c8 + kryos 7e631a2 (`tools/setup_S7.sh`).
G0: `make check` passes on stock code, both builds. G1: grl, ant, ant-full, nor, sci, sci-const bit identical to the previous K and K3-sci. G1b: sci-hub bit identical to K3-sci, hub closure 2.6e-10 / 1.5e-10 (< 1e-9), `globe_ice.nc` equals the offline merge.
Stock science run (400 yr, PaleoMIST start): no stop, clip excess 2.8e-9; t0: recon takes record 1, ice vs PaleoMIST -80 ka ESL difference -0.17 m (North) and -0.48 m (Antarctica). G3: closure 1.2e-12, `mass_resid` 2.2e-15,
(iii) full balance with `fill_int` 1.6e-13, (iv) as volumes 1.7e-6, `A_ocn_ref` +0.23 %. G4: cont2 == cont, rounding passes (0.977, bundles identical); chunk vs cont differs only in `pc_tau_max` (known) and the Maxwell memory `tau_*`, at 1e-28 Pa.

## Phase 4, S8: ncio returns the values as stored (2026-10-30)

Report `gate-reports/p4-s8.md`. fesm-utils `phase0-kryos` a13cd38 (pushed): `nc_read` no longer sets values below 1e-30 to zero (the flush was in the reader only), `test_ncio_tiny`; vilma2 main 4419636 (pushed): `test_restart_3d` with tiny values in the six memories (fails on the old ncio),
the `read_exact` bypass dropped. Stack S8 = S7 + a13cd38 + 4419636 + kryos e4453dc (`tools/setup_S8.sh`).
G1: grl, ant, ant-full, nor, sci, sci-const bit identical to the previous K and K3-sci. G4: cont2 == cont; chunk == cont except `pc_tau_max` (the known Yelmo diagnostic); the VILMA2 memories `tau_*` bit identical (S7: 1e-28 Pa); rounding 0.977, bundles identical.

## Phase 4, S9: exact `pc_tau_max` and `-0.0` (2026-10-30)

Report `gate-reports/p4-s9.md`. yelmo `topo-prepared` 30361f96 and `dev` 18b41e3d (pushed): the restart holds the `pc_taus` window and `pc_tau`, `pc_tau_masked`, `pc_tau_max`, read back; `tests/bench/check_restart.sh` (fails on the old code).
fesm-utils `phase0-kryos` 45e09c2 (pushed): `nc_read` keeps a stored `-0.0`, `test_ncio_tiny` extended. Stack S9 = S8 + 45e09c2 + 30361f96 + kryos 624f686 (`tools/setup_S9.sh`).
G1: grl, ant, ant-full, nor, sci, sci-const bit identical to K(S8) and K3-sci (0 failing, 0 absent, 0 new). G4: cont2 == cont; chunk == cont with no skip (3206 variables): `pc_tau_max` exact, VILMA2 memories exact; rounding 0.977, bundles identical.

## Phase 4, S10: state-only restore (2026-10-30)

Report `gate-reports/p4-s10.md`. fesm-utils `fields-diag` afae7a0 and kryos `rst-state` (not pushed): diagnostics (`fields_add(..., restart = .false.)`) are not in the bundle and are put again on
`PASS_RESTART`; a diagnostic an older bundle holds is ignored and logged. Stack S10 = S9 + afae7a0 + vilma2 e449d58 + kryos 3074043 (`tools/setup_S10.sh`; every package rebuilt, yelmo included).
G0: `make check` and `make check-long` (vilma=1), `make check` (vilma=0) pass. G1: grl, ant, ant-full, nor against K(S9), sci, sci-const against K3-sci: outputs and every state field of the bundles
bit identical; the only differences are the diagnostics the bundles no longer hold (regional `bsl`, `tf_shlf`). G4: sci-vilma2 cont2 == cont, chunk == cont (3200 variables, bundles included),
rounding 0.977, bundles identical; sci (fastiso) chunk == cont (3170). Backward compatibility: an S9 bundle restored by S10 + one extra diagnostic, chunk == cont (3200, bundles included, the S9
diagnostics ignored and logged); an S10 bundle restored by S10 without `bsl_ive_host`: identical except that variable, absent.

## Phase 4 WP7c: G2c, G5, G6, G7 on S8 (2026-10-30)

Report `gate-reports/p4-wp7c.md`. Branch `p4-wp7c` rebased on main b99fff4 (scripts, `par/gate/{g5-vilma2,g5-fastiso,g6-smb}.nml`, gates `g6`, `g6-smb`); S8 binaries as built. About 9.7 node-hours.
G2c (10 runs, `equil_time_max` 0 to 2e4 yr x `equil_rate_tol` 1e-3, 1e-5): at 1e-3 the cap is irrelevant (the full phase converges after 2 passes, results bit-identical) and the seed holds about 15 % of the converged deflection; 1e-5 brings it to 9 % (1-D phase alone, +3.7 min), 4 % (cap 2e3, +9 min), 2 % (5e3, +18), 0.9 % (1e4, +32), reference 2e4 (converged at 14 kyr, +43); all jobs below 1 h. Alex chose `equil_rate_tol` 1e-5, `equil_time_max` 2e4 (science par and conversion note), bundles every 5 kyr (Q12), cadence 10 yr (Q6).
G5 (ANT only, 20 to 0 ka, 2000 steps, no stop, 1.7 node-h; FastIsostasy twice for the differences): the ANT base has no time-dependent forcing (paleo terms added as the yelmox ant-paleo `-p`); the ice grows (+28 m of sea level), not a deglaciation; gained cells subside (98-99 %), lost cells uplift from -12 ka (46 to 100 %); G3 (i) 1.3e-12, (ii) 2.4e-15, (iii) 1.6e-14, (v) +0.23 %; (iv) as volumes 1.7e-2 at -8 ka (bound 1e-5 met only after -5 ka).
G6 on S8 (placeholder spin-up) and rerun on S9 with the chosen spin-up: 2 kyr from -80 ka, no stop, G3 each step passes ((iv) on the host part 5.8e-13), chunk == cont bit identical (3206 variables); smb_simple 1 kyr likewise (3014 variables) with `par/gate/g6-smb.nml` (the conversion note's `-p` list stops at `nml_validate`; note now points at the generator). Drift per kyr with the chosen spin-up: NH -0.5 m, ANT +2.05 m of sea level (K3-sci +3.7, +5.0). G2b not run (no harness yet).
G7: VILMA2 cost 0.67 s + 0.175 s/yr x dt per call (proportional to simulated time), step 4.8 s at -80 ka, bundle 340 MB written in 21 s, rounded output 1.7 x smaller; proposals: geo cadence stays 10 yr, restart every 5 kyr; 80 kyr run about 14 node-h, at least three 8-h jobs, under 25 GB.

## Phase 4, S11: VILMA2 row order (2026-10-31)

Reports `gate-reports/p4-s11.md`, `gate-reports/p4-g2b.md`. Bug: the passthrough of kryos's Gauss grid (latitude ascending) into VILMA2's north-first state copied
the rows unreversed, so the kryos side was mirrored consistently except the 3-D viscosity read in real geography (G2b diagnosis on S9: host NH `dz_bed` 51 m rms, flat with lmax;
4 m with a mirrored viscosity file). vilma2 `host-orient` c17d513 (reversing passthrough both ways, `se%C`, the legacy reader by the lat axis; `test_host_orient`, `test_drive`),
kryos `vilma-orient` (adapter on VILMA2's host-grid outputs; G2a lon-lat variant; S11 scripts), neither pushed. Stack S11 = S10 + c17d513 + kryos 3607918 (`tools/setup_S11.sh`).
G0 PASS (check, check-long, vilma=0). G2a: equal spin-up <= 2.5e-14 after the seed, no spin-up 2.9e-10 (bound 1e-9); lon-lat variant within the load's remap error (1.8e-3 against 3.0e-3; mirrored 1.4).
G2b: the difference falls with lmax, `bsl` 0.36-0.42 % of the range, `dz_bed` host NH 0.44-0.47 %, ANT 0.64-1.0 % of max; thresholds proposed. G3 (sci-vilma2, g6, g6-smb) all PASS ((iii) 1.6e-13 after aligning
`g3_budget.jl` by lat). G4: cont2 == cont, chunk == cont (3200, bundles identical), rounding 0.977. G6: no stop, G3 each step, chunk == cont (3200; smb 3008); drift as S9 (NH -0.53, ANT +2.02 m per kyr).
G2c: 1e-5 / 2e4 still converges (full phase 7 passes, 8.6 kyr, 0.87 tol; 39 min job against 55 on S8). VILMA2 restarts from before S11 are invalid. G5 and G7 predate the fix (G7's cost holds; G5 to be rerun).

## Phase 4, S11: G5 rerun (2026-10-31)

Report `gate-reports/p4-s11.md` (G5 rerun). ANT only, -20 to 0 ka, VILMA2 global on S11 with the chosen spin-up; FastIsostasy runs of S8 reused. 2000 steps, no stop, 2.85 node-h.
G3: (i) 1.3e-12, (ii) 1.8e-15, (iii) 1.7e-14, (iv) host part 9.8e-8 (series) / 6.7e-13 (streams), (v) +0.233 %: **PASS**. ANT still grows (+24.6 m, S8 +27.8): not a deglaciation.
Against S8: lost cells uplift from -15 ka (S8 -12 ka), 87-100 % per interval from -12 ka, +1 to +8 m per 500 yr; gained cells subside 98.7-100 %; min `dz_bed` -470 m (S8 -430).
The full spin-up from the LGM load hit `equil_time_max` (2e4 yr, v 2.2 tol).
