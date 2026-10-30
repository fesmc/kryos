# Phase 4 WP7, part 1b: S7, the fixes, reruns of G0, G1, G1b (2), a stock science run, G3 and G4 (2026-10-30)

Decisions of Alex on `p4-wp7a.md`, done. `$G = /work/ba1442/robinson/kryos-gate`. Stack `$G/stacks/S7` (`tools/setup_S7.sh`, log `stacks/s7-setup.log`); kryos at the head of `p4-wp7` (**7e631a2**, rebased on main ffe745c)
in two checkouts (`kryos` vilma=1, `kryos-v0`), built stock: no variant, no local patch. Runs `$G/runs/*/kx-7e631a2/S7/`, reports `$G/reports/S7/`, `$G/reports/cmp-*-2856*.log`.

## Commits

- fesm-utils `phase0-kryos` **f30d952** (pushed): `conservative_planar` pre-pass has its own arrays (`cidx0`, `cd20`), no guards; `test_conservative_project` (no cdo; two lat-lon maps onto a
  polar-stereographic grid; constant preserved) **fails on the old code with ifx** (forrtl 151, checked on Levante with the old and the new source) and passes on the new; CHANGELOG entry;
  `make check` passes on the Mac (gfortran, `openmp=0` and `openmp=1`). Shared checkout fast-forwarded, `include-serial` and `include-omp` rebuilt.
- kryos `p4-wp7` (not pushed), on main ffe745c: the test fixes (`series` rename, `test_geo_vilma2` arrays in a module); the producer-side clip in `kryos_globe_ice` (`host_clip`, `HOST_CLIP_TOL = 1e-6`, largest
  excess logged at the reference and every gather, a stop above the bound; `test_host_clip`); `HUB_TOL = 1e-9`; plan 4.5 and 13.2 (G1b (2) 1e-9, G3 (iii) full balance with VILMA2's fill integral, (iv) as volumes at 1e-5);
  `input/vilma-variables.md` synchronised with vilma2 d2006c8 (`dt_try_s`, `fill_int`; `make check` fails on a stale copy, and the restart write stops without `dt_try_s`); the `vilma` stream also writes `fill_int` and
  `fill_int_own`; `par/kryos_NH-ANT_vilma2.nml`: alpha series kept, both domains from PaleoMIST at -80 ka (`[recon_nor]`, `[recon_ant]`: `recon_regions = "all"`, `recon_zb_min = -1e5`, `recon_smooth_sigma = 0`,
  `recon_H_float_min = 0`; `[Antarctica] experiments = "recon_ant"`), conversion note choice 1; scripts and reports.
- Mac: builds with `make -j4`, the tests run with `OMP_NUM_THREADS=2`: all pass except `test_run_sequence`, which I stopped at 110 s (the rule: nothing long on the Mac); `make check` for both builds is the Levante run below.

## S7

S6 + fesm-utils f30d952 + vilma2 **d2006c8** (first built at 2676879 in the setup run, then moved with a clean rebuild; `setup_S7.sh` pins d2006c8) + kryos 7e631a2; yelmox R2' rebuilt; every package cleaned and rebuilt.

## G0: `make check` on stock code, both builds: PASS

`kryos` (vilma=1) and `kryos-v0` (vilma=0), gate flags, jobs 28569163 and 28569164 (22 and 21.6 min): all tests passed, `test_globe_ice` and `test_geo_vilma2` included.

## G1: PASS, bit identical (`--restarts`, maps S5 files, S7 vilma=1; jobs 28569167-76, compare 28569668)

| gate | against | variables |
|---|---|---|
| grl | K(S5b) | 1270 |
| ant | K(S5b) | 1362 |
| ant-full | K(S5) | 1362 |
| nor | K(S5) | 1251 |
| sci | K3-sci | 2579 |
| sci-const | K3-sci | 2579 |

0 failing, 0 absent, 0 new. (The same six passed at the previous S7 head 64317fc, before the stream change.)

## G1b, hub as observer, stock code (sci-hub cont, job 28569179, compare 28569669)

| check | result |
|---|---|
| (1) vs K3-sci | PASS, bit identical, 2579 variables |
| (2) hub closure | reference **2.55e-10**, merged ice 1.48e-10 to 1.51e-10 at every output: **PASS at 1e-9** |
| (4) `globe_ice.nc` vs the offline merge | equal (coverage 2.2e-16, `H_ice` 1.4e-3 m, the float32 rounding of the inputs) |

cont2 and the domain swap were not rerun (they passed at S6; the hub code of the observer is unchanged by the clip, which is a no-op there except for the 2.8e-9 excess).

## Science run on stock code (gate `sci-vilma2`, 400 yr from -80 ka, lmax 170, 3D, 1-D pre-spin-up): no stop

Jobs 28569181 (cont), 28569182 (cont2), 28569185/86 (chunk), 28569188 (lossless); 13.7 min, peak RSS 11.6 GB as before.
- The clip works at the polar row: excess of `host_w` **2.83e-9** at the reference and every gather (bound 1e-6), of `host_C` 0; the 1-D spin-up of the PaleoMIST start is one pass (66.7 yr, mean bed rate 4.1e-4 m/yr).
- **t0 check, PaleoMIST at -80 ka in both domains** (`recon_t0.jl`): recon reads record 1 (80 000 years ago) in both. The run's first ice record against the file's record 1 on the ice grid:
  North: max |dH| 1102 m (rms 43.6 m), volume 5.102e15 vs 5.032e15 m3, ESL -12.482 vs -12.312 m, **difference -0.170 m**;
  Antarctica: max |dH| 1671 m (rms 80.7 m), 2.798e16 vs 2.778e16 m3, ESL -68.448 vs -67.967 m, **difference -0.481 m** (ESL = -(rho_i/rho_sw) V / 3.618e14 m2, rho_i 910, rho_sw 1028).
  The loads differ by 0.17 m and 0.48 m of sea level (0.4 to 0.7 % of the domain's ESL); the cells of the largest differences (1.1 to 1.7 km) were not located. The gathered t0 load against the reference (present day):
  Antarctica -3.81 m, North -4.88 m, external +3.6e-5 m (the base's LGM recon gave North -67 m).

## G3 (the stock cont and lossless runs; `g3_budget.jl`, `$G/reports/S7/g3-budget.md`)

| budget | result | bound |
|---|---|---|
| (i) merge closure, 40 gathers | max 1.2e-12 | 1e-9: PASS |
| (ii) `mass_resid` | max 2.2e-15 | 1e-6: PASS |
| (iii) full balance `int C rsl dA = -(rho_i/rho_w) int dI_g dA + fill_int`, 9 records (lossless run), `dI_g` recomputed from the gathered terms | residual **1.6e-13** of `abs(int C rsl dA)`; the recomputed term equals VILMA2's `bsl_ive A_ocn` to 9.5e-8 (float32); `fill_int - fill_int_own` equals the host `int F dA` to 8e-14 | 1e-6: PASS |
| (iv) `bsl_ive` vs the sea-level equivalent of the gathered grounded ice, as volumes | max **1.7e-6** (as sea level 3.6e-4) | 1e-5: PASS |
| (v) `A_ocn_ref` | 3.6264e14 m2, +0.23 % | 1 %: PASS |

## G4 (400 yr, chunks at -79800): cont2 and rounding pass; chunk == cont passes except 6 variables at 1e-28

| check | result |
|---|---|
| cont2 vs cont | **PASS**, 2596 variables |
| rounded cont vs lossless | bundles **identical** (248 variables); max error / tolerance **0.977 <= 1**, no variable outside |
| chunk vs cont, bundles included | 3198 variables, **7 fail**: `ice_ant_restart.nc:pc_tau_max` (the known diagnostic that a 200-yr chunk does not refill; skip in `test_run_sequence`) and `geo_glob_restart.nc:tau_{a,b,c}_{re,im}`, the Maxwell memory, with max |diff| **1.0e-28 to 1.8e-28 Pa** in 75 000 to 127 000 of 12.9e6 cells. |

Everything else is bit identical now, the Earth's state, `rsl`, `dz_bed`, `dz_ss`, `bsl`, `dphi`, `rot_*` and every file of the other components (at S6 41 variables differed, up to 7e-12 m). The remaining differences are at 1e-28 Pa (the memory
stress is of the order of 1e4 to 1e6 Pa where it matters): numbers far below any physical effect, but not bit identical; which cells they are in was not looked at. I did not add a skip: pass or fail by the letter is for you and the VILMA2 agent.

## Open

1. G4: the 1e-28 Pa residue in `tau_*` (VILMA2) and the `pc_tau_max` skip for 200-yr chunks.
2. The 1.1 to 1.7 km maximum differences of the recon start against the file (rms 44 and 81 m; the ESL difference is small).
3. cont2 and the domain swap of the hub gate, and G2a with d2006c8, were not rerun.
