# Phase 4 WP7, part 1: stack S6, G0, G1, G1b, G2a, G3, G4 (2026-10-30)

`$G = /work/ba1442/robinson/kryos-gate`. Branch `p4-wp7` (main cd85c14), scripts `tools/setup_S6.sh`, `tools/levante/{gate_S6,compare_S6,check_kryos,run_vilma2,g34_vilma2,g2a,cmp_g1b,cmp_g4}.sh`,
`tools/levante/{g2a_compare,g3_budget,hub_offline,round_check}.jl`, `tests/{g2a_cap,hub_areas,mapgen_pairs}.f90`, gates `tools/gate/{sci-hub,sci-vilma2}.toml`.
Reports on Levante: `$G/reports/S6/`, `$G/reports/cmp-*.log`. Runs: `$G/runs/<gate>/kx-cd85c14/S6*/`.

## Stack S6

S5b + fesm-utils `phase0-kryos` 0a626e9 (every package rebuilt, gate flags), vilma2 2676879 (`vilma-static`), yelmox topo-prep 5ac89c8 (R2', `bin/R2`).
kryos cd85c14 in two checkouts, both built (setup log `stacks/s6-setup.log`, `stack.lock`): `stacks/S6/kryos` (`vilma=1`), `stacks/S6/kryos-v0` (`vilma=0`).
`maps/S6` frozen (read-only): maps/S5 (44 files, unchanged, `cmp`) + 15 new (below).

## Identity (`--restarts`, cont, 1 kyr)

| comparison | grl | ant |
|---|---|---|
| K(S6, vilma=0) vs K(S5b) | PASS 1270 variables | PASS 1362 |
| R2(S6) vs R2(S5b) | PASS 1274 | PASS 1368 |

0 failing, 0 absent, 0 new in all four (jobs 28561241/44, ref-S6 28561242/46, compare 28561331).

## G0: builds, `make check`

Both builds succeed. `make check` at **cd85c14 fails on Levante (ifx 2023.2, `-O2 -fp-model precise`, 16 threads)**, in three tests; the rest pass:

| test | build | failure | cause | fix (on `p4-wp7` or in fesm-utils) |
|---|---|---|---|---|
| `test_forcing_increment` | both | does not compile: internal function `series` conflicts with the module name `series` | ifx only | rename (`7899b72`) |
| `test_geo_vilma2` | vilma=1 | adapter and direct calls differ at the reference; VILMA2 stops "host fraction outside its bounds" | with ifx the main program's body did not see the whole-array assignments of its internal procedures (`reference_terms` set `w`, `z_ref`: zero in the body, right inside). The same pattern in my G2a harness gave all-zero exports. Not a kryos or VILMA2 defect | the arrays in a module (`1d408d5`) |
| `test_globe_ice` | both | `forrtl 151: allocatable array is already allocated` in `conservative_planar` (fesm-utils `conservative.f90:313`) | ifx: an allocatable that is allocated and freed in the serial part (`cidx`, `cd2` of the `do_project` pre-pass) and then listed `private` in the OpenMP region is "already allocated" on entry. **It also stops every new conservative map (mapgen)**; the old maps (made by v24 code) hid it | fesm-utils: rename the serial temporaries (`cidx0`, `cd20`); patch `p4-wp7a-fesm-utils-conservative-ifx.patch` (the two `if (allocated)` guards in it are not needed) |

With these fixes (checkouts `p4w7/kryos-fx1` vilma=1, `kryos-fx0` vilma=0 against a copy of the stack's fesm-utils with the patch): **`make check` passes, both builds** (jobs 28562835/36, 19.5 and 20.4 min). The S6 kryos binaries are those of cd85c14 (the failing tests are test code and the map generation). **The fesm-utils fix needs a commit in fesm-utils (and a stack S7) to make G0 pass as specified.**

## Maps (GAUSS-L170, `$G/maps/S6`)

The new con maps could only be made with the fesm-utils patch above (a copy `p4w7/fu-patched`, kryos `p4w7/kryos-mg` linked to it; the maps are data, no numerical change).
New (15): `grid_GAUSS-L170.txt` (written by kryos, radius = r_earth), `grid_GLOBAL-0.5DEG.txt` (FesmData), con ice -> GAUSS-L170 for NH-32KM and ANT-32KM, con ice -> GLOBAL-0.5DEG (NH-32KM, ANT-32KM),
con GLOBAL-0.5DEG -> GAUSS-L170, NH-32KM, ANT-32KM, bilin GAUSS-L170 -> NH-32KM, ANT-32KM (the run's own coupler made all of these); GAUSS-L170 <-> ANT-16KM con and bilin (`tests/mapgen_pairs.f90`, from the grid descriptions; not used by the science case, whose ports are on the ice grids).

## G1: no regression (`--restarts`, maps S5 files, cont, S6 vilma=1)

| gate | against | result |
|---|---|---|
| grl | K(S5b, 5a663a3) | PASS 1270 variables |
| ant | K(S5b, 5a663a3) | PASS 1362 |
| ant-full | K(S5, 5a663a3) | PASS 1362 |
| nor | K(S5, 888884d) | PASS 1251 |
| sci | K3-sci | PASS 2579 |
| sci-const | K3-sci | PASS 2579 |

0 failing, 0 absent, 0 new (jobs 28561552-58, compare 28561734).

## G1b: the hub as an observer (gate `sci-hub`: `par/gate/sci-hub.nml` = sci.nml + `[globe_ice]`, target GAUSS-L170; binary = the S6 kryos with the changes below)

| check | result |
|---|---|
| (1) run with `[globe_ice]` vs K3-sci (sci cont) | **PASS, bit identical**, 2579 variables (jobs 28562816, compare 28563174) |
| (3) cont2 vs cont; domains swapped vs cont | **PASS, bit identical**, 2586 variables each (bsl budget: 100 records, 2.2e-8) |
| (4) `globe_ice.nc` (last record) vs the offline merge (`hub_offline.jl`, areas from `hub_areas.x`) | equal: coverage 2.2e-16, `H_ice` 1.4e-3 m = float32 rounding of the inputs (ulp 5e-4 m at 4.5 km); tolerance 3e-3 m. The ice of the last gather is that of t_out - 10 yr |
| (2) hub closure | **FAIL at 1e-12 as specified**: reference 2.55e-10, merged ice 1.48e-10 to 1.51e-10 at every output. Passes at 1e-9. The run stops at the first hub closure (HUB_TOL = 1e-12) |

**Two stops with the real maps, both at the first VILMA2-path gate (before any run reached its step)**:
1. `HUB_TOL = 1e-12` (kryos_globe_ice.f90:77): 2.55e-10 at the reference (the merge weights of the con maps are exact only to about 1e-10).
2. VILMA2's `solid_earth_host_check` (tolerance 1e-9) rejects `host_w = 1 + 2.8e-9` at the polar row of the Gauss grid (cell 85,1; Antarctica); kryos's own `check_host` admits it (`overlap_tol` 0.01).

To go on, a **variant** (`$G/p4w7/kryos-mg` = cd85c14 + diff `p4-wp7a-kryos-hub-variant.diff`: `HUB_TOL = 1e-9`, `host_w` clipped to 1 and `host_C` to `host_w` in the reference and the gather) built on the patched fesm-utils was used for sci-hub and all vilma2 runs below.
These two are decisions for you (tolerance, or clip in the merge, or exact maps); the variant is not proposed code.

## G2a: VILMA2 standalone vs the adapter (`tests/g2a_cap.f90`, lmax 170, 3D Bagge structure 6, PREM, 64 threads)

Land cap (radius 12 degrees at 90E 58N, 2 km, +25 % per 10-yr step, 6 steps) inside a land box (host region, w = 1, no coast, no external ice), same field to `vilma.x` (Gauss grid, `remap_input = .false.`, `i_eq = 2`) and to `geo_glob`. Differences relative to max |field| (`g2a_compare.jl`):

| spin-up | record 0 (seed) dz_bed / dz_ss | records 1-6 dz_bed | records 1-6 dz_ss | bsl, bsl_ive, C_ocn |
|---|---|---|---|---|
| none (both) | not comparable: the standalone seed solves with the reference ice (zero, F7), the host seed with the new load | 3.2e-10 (step 1), then <= 8e-13 | 5.8e-10 (step 1), then <= 1.4e-13 | bsl 1.5e-15, bsl_ive and C_ocn bit identical |
| `pre_spinup_1d`, `equil_rate_tol` 1e-6 (both, 12 passes) | 4.1e-6 / 2.6e-5 | <= 3.4e-14 | <= 4.0e-13 | as above |
| `pre_spinup_1d`, default tol 1e-3 (science case) | 0.21 / 0.027 | 0.11 to 0.22 | 0.012 to 0.027 | bsl 1.5e-15, bsl_ive and C_ocn bit identical |

The adapter equals VILMA2 to round-off when both do the same spin-up. With the default tolerance they stop the 1-D spin-up at different passes: the host seed already holds the instantaneous response, so pass 1 (66.7 yr) moves the bed by 0.07 m (max 9.7 m) against 7.1 m (max 82 m) in the standalone run, whose pass 1 contains the jump; kryos reads "converged" after pass 1 (mean rate 3.4e-4 m/yr < 1e-3), the standalone run after pass 2. The centre of the cap is then 13 m (21 %) apart, a constant offset to the end: an incompletely relaxed state, not an adapter error. **Proposed G2a threshold: 1e-9 of max |field| for records after the first (observed <= 6e-10), spin-up held equal; the seed record and the spin-up are G2c.** The adapter was compared to `vilma.x` only on a Gauss-grid load: the gather (merge) is tested by G3 and `test_globe_ice`; G2b is the realistic comparison.

Cost (64 threads, Levante AMD, `-O2`): one `solid_earth_update` at dt = 10 yr, T170 3D: **0.65 to 0.93 s** (standalone driver mean 0.74 s, plus 0.24 s output); setup + reference ~15 s; `geo_init` 0.4 s; peak RSS 8.5 GB for the harness.

## G3: budgets (gate `sci-vilma2`, 400 yr from -80 ka, `geo = vilma2`, lmax 170, 3D Bagge 6, 1-D pre-spin-up, dt 10; `g3_budget.jl`, report `$G/reports/S6/g3-budget.md`)

| budget | result | bound |
|---|---|---|
| (i) merge closure, 40 gathers | max **7.1e-13** | 1e-9: PASS |
| (ii) VILMA2 `mass_resid`, 41 steps | max **2.1e-15** | 1e-6: PASS |
| (iii) `int C rsl dA` against `-(rho_i/rho_w) int dI_g dA + int F dA`, `dI_g` recomputed from the gathered terms (lossless run) | the recomputed `-(rho_i/rho_w) int dI_g dA` equals VILMA2's `bsl_ive A_ocn` to **8e-8** (float32 storage); the balance leaves a residual of 2.47e14 to 2.49e14 m3 (1.05 % of `int C rsl dA`), nearly constant in time | 1e-6 not demonstrable: the residual is the own part's subgrid fill, `zeta0 (C_own - C0_own)`, which needs `zeta0` (not an output) |
| (iv) `bsl_ive` against the sea-level equivalent of the gathered grounded-ice change (domains + external) | as sea level: **3.4 %** (-69.20 m vs -66.96 m); as volumes (`bsl_ive A_ocn` against ESL times the sea surface area of the constants): **4.3e-7** | 1 %: the difference is the ocean area (VILMA2 divides by the area of its solve, 3.50e14 m2 at the end, against 3.618e14 m2 of the constants in the log's ESL) |
| (v) `A_ocn_ref` | 3.6264e14 m2, +0.23 % of 3.618e14 | 1 %: PASS |

The NH domain holds -67.3 m ESL of ice against the reference at t0 (the recon's LGM ice, see the conversion note), Antarctica -0.26 m, external +3.6e-5 m; `bsl` -67.4 m at t0, -67.06 m at the end. Proposals: (iv) as volumes at 1e-5, or `ESL` normalised with VILMA2's ocean area; (iii) a VILMA2 diagnostic of `int fill dA`, or the check `bsl_ive A_ocn` against the recomputed term (8e-8) plus (ii).

## G4: chunk == cont, rounding (same runs, 200-yr chunks at -79800)

| check | result |
|---|---|
| cont2 vs cont (`--restarts`) | **PASS**, 2588 variables |
| chunk vs cont, bundles included | **FAIL**: 3188 variables, 41 failing, all in c01 and all in the Earth's state and its copies: `geo_glob_restart.nc` (`rsl`, `z_bed`, `dz_bed`, `dz_ss`, `bsl`, `dphi`, `tau_*`, `rot_*mem`), `geo_{glob,nor,ant}_{imports,exports}.nc`, `geo_glob_vilma.nc` (bsl, dphi), `geo_glob_ts.nc:mass_resid`. Size of the differences: `bsl` 1.4e-14 m (rel 1.5e-15), `dphi` 9e-14 m, `dz_bed` 7e-12 m (max rel 3.8e-9), `mass_resid` 2e-16. c00 and the bundle of c00 at -79800 equal cont's; the first difference is within five steps after the restore. Every Yelmo, climate, surface and shelf file passes, except `ice_ant_restart.nc:pc_tau_max` (the known diagnostic that a 500-yr chunk refills and a 200-yr chunk does not, skip in `test_run_sequence`) |
| rounded cont vs lossless: bundles | **identical**, 2 bundles, 242 variables |
| rounded cont vs lossless: files | **max error / tolerance 0.977 <= 1**, no variable outside, 30 files (`round_check.jl`) |

The restart is exact in the harness at lmax 16 (`test_geo_vilma2`, `test_gia_sequence`); at T170 with the 3D memory, rotation, sub-stepping and the 1-D spin-up it is exact to round-off only. Not investigated further; candidates are the iteration and stepper state that is not in the file (warm start of the sea-level fixed point, `sigma_lm`, the stepper seed), which would enter at the tolerance of the iteration. A test of `test_geo_vilma2` at lmax 170 with 3D viscosity would isolate it.

## Measurements for G7 (science case, 64 threads, exclusive node)

| quantity | value |
|---|---|
| wall of the 400-yr run | 15.4 min (init, with the Yelmo cold starts, about 12 min) |
| per 10-yr step (timers, cumulative at 400 yr / 40 steps) | total 5.3 s: `globe` (gather + VILMA2 + forcing, column `forcing`) 2.1 s, ice (two domains) 1.6 s, output 1.0 s, climate 0.25 s, surface 0.22 s, ports 0.006 s, hub 0.007 s |
| VILMA2 call | 0.65 to 0.93 s (G2a, same node) |
| spin-up at t0 (1-D, tol 1e-3) | 6 passes of 64 to 70 yr, 1030 cpu-s summed over threads; with `equil_rate_tol` 1e-6, 12+ passes (pass 12: 490 sub-steps, 3557 cpu-s; the run took 21.5 min) |
| peak RSS | 11.6 GB (kryos.x, 64 threads, T170 3D + NH + ANT); sci-hub (no VILMA2, 16 threads) 3.5 GB |
| restart bundle | 341 MB: `geo_glob_restart.nc` 211.4 MB, geo_glob imports 3.3 MB, exports 3.1 MB, the domains the rest |
| output of 400 yr, 2D every 50 yr | rounded: 181 MB all files, `geo_glob_vilma.nc` 19.0 MB, `geo_glob.nc` 6.2 MB; lossless: `geo_glob_vilma.nc` 57.7 MB, `geo_glob.nc` 11.7 MB, `geo_glob_gia.nc` 9.2 MB |

## G2c: what I would test (needs Alex's `equil_time_max`)

`pre_spinup_1d` alone relaxes about 400 yr and stops at a mean bed rate below 1e-3 m/yr (the cap-centre is then still moving at ~0.2 m/yr, G2a). Proposed grid: `equil_time_max` = 0, 2e3, 5e3, 1e4, 2e4 yr, each with `equil_rate_tol` = 1e-3 and 1e-5 m/yr; the report per value: bed rate and `dz_bed` at the seed, drift of `dz_bed` over the first 5 intervals, wall. Cost grows steeply (sub-steps per pass 3 to 490 within 12 passes), so the long values need a node-hour estimate first. The criterion (mean rate) is dominated by the globe; a local one (max rate in the host region) is worth a look.

## Open for you

1. G0: commit the fesm-utils fix (and S7), and the test fixes (`p4-wp7`).
2. G1b (2) and the host check: tolerance or clip (above).
3. G4: chunk != cont at round-off in the Earth.
4. The NH start at -80 ka (recon at -26 ka) and alpha (conversion note).
5. `par/kryos_Antarctica.nml` and `par/kryos_NH-ANT.nml` were stale (no `[shelf_corr_ant]`, `shelf_ant` stopped at init): regenerated with the converter (only that group moves).
