# WP0 baselines (2026-10-28/29)

Scope after the cut (coordinator, Alex): stack S4, S4 identity checks (grl, ant), NOR maps and NOR R0/R1/R2'.
Gate root `$G = /work/ba1442/robinson/kryos-gate`. Scripts: `tools/setup_S4.sh`, `tools/levante/*` (copies of
`$G/tools/*_wp0.sh`, `build_stack_wp0.sh`, `compare_wp0.sh`, `wp0_r2diff.py`; the older `$G/tools` scripts are unchanged).

## 1. yelmox topo-prep commit 2 (step 1)

- `topo-prep` was clean at 74f3a2b. New commit **5ac89c8** (+4 lines in each of the 8 base files, 32 in all).
  Values read from each base's own Yelmo keys: `grad_lim_zb` = ytopo value (PYR 1.0, others 0.5), `z_bed_f_sd = 0.0`
  (yelmo_init_topo, every base), `remove_lakes = True`, `H_ice_min = 1.0`. All as expected, no deviation.
- Moved by git bundles (`$G/bundles/yx2.bundle` 918d996..topo-prep, `ym4.bundle` b5ab5d2c..ea90fee7). Nothing pushed.

## 2. Stack S4 (step 2)

`$G/stacks/S4` = copy of S3 (not modified), `tools/setup_S4.sh` (log `stacks/s4-setup.log`, 12 min). Gate flags
`-O2 -fp-model precise -march=znver3 -traceback` on every compile line (R1 186 of 186, R2' 31 of 31).

| package | commit |
|---|---|
| yelmox | R1: v2.4 918d996; R2': topo-prep 5ac89c8 (74f3a2b + [topo] keys of 8 bases; code identical to 74f3a2b) |
| yelmo | topo-prepared **ea90fee7** (b5ab5d2c + vm-m16 restart fix; S3: b5ab5d2c) |
| fesm-utils | bb39f2f (as S3) |
| FastHydrology, elsa, tracer, FastIsostasy, chion, insol | as S3 (a0d5e8e, 74d426d, 208ac3a, ce0d6ff, 0de32ba, 42b254e) |
| rembo1 | v1.3 9c24c42, **linked** (`make yelmox rembo=1`; `climate_rembo.o`, `-lrembo`) |

- Binaries kept (runme copies the exe at staging): `stacks/S4/bin/R1/yelmox.x`, `bin/R2/yelmox.x`
  (sha256 7cec96f1..., 2cbbc616...). Lock: `stacks/S4/stack.lock`.
- Every Fortran package was cleaned and rebuilt (`build_stack.sh`, new `MAKE_EXTRA="rembo=1"`).
- kryos (`$G/kryos`, commit af8f1b9, unchanged) relinked to S4, `make clean`, rebuilt by `run_kryos.sh` (new `STACK=S4`).

Identity checks (`kryos-compare` identity, `--restarts`, cont 1 kyr, reports `$G/reports/WP0/`):

| pair | grl | ant |
|---|---|---|
| R1(S4) vs R0 (v24) | PASS, 1274 vars | PASS, 1368 vars |
| R2(S4, 5ac89c8) vs R2 (74f3a2b, S3) | PASS, 1274 vars | PASS, 1368 vars |
| K(S4) vs K(S3), kryos af8f1b9 | PASS, 1249 vars | PASS, 1341 vars |

0 failing in all six. Jobs: R1(S4) grl 28527985, ant 28527986; R2(S4) grl 28527987, ant 28527988;
K(S4) grl 28528005, ant 28528006 (`runs/{grl,ant}/kx-af8f1b9/S4/cont`). The vm-m16 fix is neutral for grl/ant.

## 3. Maps (step 3)

mapgen R0 runs on v24 (one per base, private copy of `maps/S` each, time_end = time_init + 10), then merged into
**`$G/maps/S4`** (= `maps/S` + new pairs, merge script checks duplicates; frozen `chmod -R a-w`; `maps/S` untouched).
Jobs and new pairs (after the cut only NOR is used; the others were made before the cut arrived):

| base | job | wall | new map files |
|---|---|---|---|
| NOR | 28527693 | 6:44 | `bilin_GLOBAL-1DEG_NH-32KM` (2.2 MB) |
| LIS | 28527695 | 2:34 | `bilin_GLOBAL-1DEG_LIS-32KM` |
| PYR | 28527697 | 6:18 | `bilin_PYR-2KM_PYR-250M`, `bilin`/`con`/`nn_PYR-250M_PYR-2KM` |
| SRG | 28527699 | 1:35 | `con_SRG-250M_SRG-16KM` (`bilin_SRG-16KM_SRG-250M` was not generated) |
| eGRL | 28527704 | 1:48 | `bilin_GRL-16KM_GRL-8KM` |
| GRP | 28527692 | 1:24 | none |
| eANT | 28528010 | 5:06 | none (ANT-8KM esm data on its own grid) |

eANT needed 60 GB (OOM at 15 GB; 28527703, 28527746, 28527967 failed: OOM, OOM, `--cpus-per-task` conflict). On `shared`
(940 MB per CPU) more memory means more CPUs: `submit_ref.sh` `NCPU`, threads stay 16 (`runme --cpus` does not change the template).
rGRL was not run on v24 (v24 has the rembo stub; it needs S4).

## 4. Baselines: NOR (step 4)

`gate nor`: `-p ctrl.write_lossless=True tm_2D.dt=100 tm_2Dsm.dt=10 tm_1D.dt=10 tm_rst.dt=500`, 1 kyr, maps `maps/S4`,
partition `shared`, 16 threads, 1 h limit (`tools/levante/submit_ref.sh`, `submit_wp0.sh`). Runs in
`$G/ref/nor/yx-918d996/v24` (R0), `.../yx-918d996/S4` (R1), `.../yx-5ac89c8/S4` (R2').
All 12 jobs COMPLETED, exit 0, empty `out.err`.

| run | R0 (v24) | R1 (S4) | R2' (S4) |
|---|---|---|---|
| cont | 28528164, 11:08 | 28528181, 10:35 | 28528191, 10:28 |
| cont2 | 28528421, 10:35 | 28528426, 10:50 | 28528429, 10:39 |
| chunk c00 / c01 | 28528423 / 28528424, 8:45 / 7:29 | 28528427 / 28528428, 9:07 / 7:30 | 28528430 / 28528431, 8:51 / 7:33 |

Comparisons (`tools/levante/compare_wp0.sh base nor`, KryosTools on Julia 1.13; reports `$G/reports/WP0/nor-*.md`, `summary.txt`):

| comparison | result |
|---|---|
| R1 vs R0, cont, `--restarts` | **PASS**, 1248 vars |
| cont2 vs cont: R0, R1, R2' | PASS, PASS, PASS (1248 vars each) |
| R1 chunk vs R1 cont (`--chunks`) | **PASS**, 1080 vars |
| R2' chunk vs R2' cont | **PASS**, 1080 vars |
| R0 chunk vs R0 cont | FAIL, 346 of 1080 vars, first divergence t = 510 (`c01/yelmo_sm.nc:H_ice`, ...) |
| R1 vs R0, chunk c00 | PASS, 894 vars |
| R1 vs R0, chunk c01 | FAIL, 346 vars (the same R0 failure) |
| R2' vs R0, cont | 112 of 186 vars differ (documented below) |

- Gate rule "R1 == R0 bit for bit" holds for the continuous run and for c00. c01 differs because R0 (Yelmo v2.4)
  restarts inexactly with `vm-m16` calving (`eps_eig_1/2` not in the restart, gate-log) and R1 has the fix: the
  chunk run is restart-exact only on S4 (R1 and R2' chunk == cont). This is the Q9 expectation.
- R2' vs R0 (`tools/levante/wp0_r2diff.py nor`): hub preparation log of R2': lakes removed; ice < 1 m removed in
  8385 cells; `grad_lim_zb = 0.5` needed 0 passes; Yelmo "prepared ... cleaning skipped" on both paths.
  At t = 0 (`yelmo.nc`): `H_ice` vs `H_ice_ref` within one run differs in 32308 cells (R0) and 31435 (R2'),
  as NOR starts from its recon/initial state; R2' vs R0: `H_ice` 21419 cells (118 by > 10 m, max 234.7 m),
  `H_ice_ref` 5436 cells (11 by > 10 m, max 50.3 m), `z_bed` identical. Through 1 kyr `V_ice` differs by at most
  0.0116 (1e6 km3), 0.0111 at the end. First divergence at time 0 (`isos.nc:z_ss`, `bsl`, `H_ice`, `H_ice_ref`).
  Differing files: chion, htopo, isos, isos_ts, mshlf, snap, yelmo, yelmo_sm, yelmo_ts.

## 5. Dropped bases (scope cut)

- No mapgen or baseline job was submitted after the cut; no queued job existed to cancel (all submitted jobs had
  started). The cut arrived after R0/R1/R2' `cont` of GRP, LIS, PYR, SRG, eANT, eGRL had been submitted (21 jobs,
  rGRL never submitted). Not compared. State at the time of this report: GRP, LIS, SRG, eGRL cont completed
  (2:22 to 15:29 wall); the PYR (R0 28528167, R1 28528183, R2' 28528193) and eANT (28528171, 28528185, 28528195)
  cont jobs were still running (4 h limit) and finish on their own. Outputs: `$G/ref/{grp,lis,pyr,srg,eant,egrl}/`.
- rGRL: held (coordinator: tm_1D.method = "file", R0 on v24 impossible: stub). No rGRL job.

## 6. Open points

- `SRG-16KM -> SRG-250M` bilin map not generated (probably unused).
- The Julia compare needs `julia +1.13` (the depot is for 1.13; the default channel is 1.12).
- Phase 1 refs (`$G/ref/{grl,ant}`) and S/S2/S3/v24 untouched; new directories only.

## 7. NOR: why `H_ice_ref` and initial `H_ice` differ at t = 0 (R2' vs R0)

Read-only analysis of the R0/R2' `out.out` and outputs, and of the code (yelmox `topo-prep`, Yelmo `topo-prepared`
ea90fee7). Scripts `tools/levante/wp0_nor_why.py` (run as `$G/tools/nor_why.py`; a second small script printed the
largest `H_ice_ref` differences). No job was submitted.

**Ruled out** (identical in R0 and R2', from the par file and the logs):
- the Yelmo cleaning keys of the North base equal the `[topo]` values of 5ac89c8: `grad_lim_zb` 0.5 (0 smoothing
  passes), `z_bed_f_sd` 0, `smooth_H_ice` = `smooth_z_bed` = 0; `z_bed` is identical at t = 0 (0 cells);
- `remove_lakes`: both remove lakes with sea level 0 (R0 in `ydata_load` / `yelmo_init_topo`, R2' in the hub);
- the 1 m threshold on the present-day reference: Yelmo's PD path does `where (H_ice < 1) H_ice = 0`
  (`yelmo_data.f90`, `ydata_load`), equal to the hub's `H_ice_min = 1`. So for the PD reference
  (`dta%pd%H_ice`) the two paths agree, as for grl.

**Cause for `H_ice_ref` (5436 cells, max 50.3 m).** NOR has `init_method = "recon"`. `domain_init_recon`
(`libs/kryos_startup.f90`) overwrites `bnd%H_ice_ref` first with the ICE-6G_C slice, then, after the 1 yr fixed-topography
equilibration and the removal of floating ice thinner than 50 m, with `tpo%now%H_ice`
("Store the clean thickness as the reference state"). The `H_ice_ref` written at t = 0 is therefore the initial
ice after the recon recipe, not the prepared PD reference, and it inherits the initial-state path of Yelmo
(`yelmo_init_topo`), which is the one that differs:
- R0 (cleaning applied): `0.1 <= H_ice < 10 m` is raised to **10 m** (`yelmo_ice.f90`, `yelmo_init_topo`; ice < 0.1 m removed).
- R2' (`prepared`, cleaning skipped): the hub's field as prepared: ice < 1 m removed, ice in 1 to 10 m kept as is.
  (Same split as the grl "initial H_ice" difference, gate-log: 636 cells.)
- Outside the recon region (or where the bed is below -500 m) the initial ice is this field; inside it the recon
  ice is the same in both, but `smooth_gauss_2D(H_ice, f_sigma = 2)` over the whole field spreads the 10 m steps
  over a few cells. Numbers: of the 5436 differing cells 5425 differ by at most 10 m; 5300 have R0 > R2'
  (thin ice raised); 3560 lie within 6 cells of a source cell with 0.1 <= H < 1 m (the hub field before
  preparation has 8385 cells with 0 < H < 1 m according to the hub log).
- The two largest differences (50.3 and 47.8 m) are cells where R2' has 0 and R0 has 50.3 / 47.8 m: the 10 m shift
  moves the smoothed thickness across the 50 m threshold of the floating-ice removal (`mask_bed == 5 and H < 50 -> 0`).
  The next largest are 13 to 15 m (smoothing of neighbouring raised cells). Only 11 cells differ by more than 10 m.
- grl/ant differ because there `init_method` is `none`/`equil`: `H_ice_ref` stays the prepared PD reference and is
  identical when hub and ice grid agree.

**Cause for the initial `H_ice` (21419 cells, 118 by > 10 m, max 234.7 m).** The `H_ice` written at t = 0
is after the whole recon recipe, which ends with 200 yr of dynamics with free topography
(`yelmo_update_equil`, `beta_min = 100`). It starts from the initial ice that differs as above (R0 thin ice at
10 m, R2' removed or thinner; total volume 3.5263e7 vs 3.5261e7 in the output units), and the hub/surface/climate are
refreshed from it, so the differences grow and spread:
- 4242 of the 21419 cells are within 6 cells of a differing `H_ice_ref` cell, 7010 within 30 cells; the rest is
  dynamic spread over the 200 yr (and the isostasy reference follows: `isos.nc:z_ss` differs by 0.76 m at t = 0,
  `bsl` by 0.017 m);
- 21340 of the cells have ice in both runs (35 / 44 only in R0 / R2'), i.e. margin and thickness changes, not new ice
  patches. The few large differences (118 cells > 10 m, max 234.7 m) are individual dynamic cells.
This chain was not isolated cell by cell; doing so needs a run that dumps `H_ice` before the smoothing in both
variants. It is not needed for the baseline (R2' is deterministic, restart-exact, and documented).

**Consequence for the NOR gate:** K == R2' (WP4, recon) must reproduce the recon recipe exactly, including the
Yelmo initial-state path being skipped (`prepared`): the hub's prepared field is the initial ice outside the recon
region and the smoothing runs on it; the R0-style 10 m raise does not apply in R2'.
