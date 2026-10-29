# Phase 4 WP0: stack S5, gate G0, maps, data check (2026-10-29)

`$G = /work/ba1442/robinson/kryos-gate`. Scripts: `tools/setup_S5.sh`, `tools/gate/stack-S5.toml` (id 95a3007c),
`tools/levante/{run_kryos,vilma_link_test,compare_g0}.sh` (on Levante `run_kryos_p4.sh`; log `stacks/s5-setup.log`).

## Stack S5 (`$G/stacks/S5`, `stack.lock`)

S4 plus: fesm-utils `phase0-kryos` **5662899**; SHTns 3.7.5 built in fesm-utils (`make shtns`, serial and omp);
vilma2 `main` **5e5a9ee** (`make vilma-static openmp=1`, `$G/stacks/S5/vilma2`, `fesm-utils` linked to the stack's).
Every Fortran package cleaned and rebuilt with the gate flags (all 186 yelmox-side compile lines carry them);
fesm-utils serial library rebuilt. yelmox topo-prep 5ac89c8 (R2', rembo=1) -> `bin/R2/yelmox.x` (sha256 c6e15552...); no R1.
Unchanged from S4: yelmo ea90fee7, FastHydrology a0d5e8e, elsa 74d426d, tracer 208ac3a, FastIsostasy ce0d6ff,
chion 0de32ba, insol 42b254e, rembo1 9c24c42; fftw/lis externals not rebuilt.
vilma2 was first built at 40a94db (setup run), then moved to 5e5a9ee (output precision) with `make clean` + rebuild;
`setup_S5.sh` pins 5e5a9ee. kryos `main` ee1f91d (pulled from origin) built `vilma=0`: `kryos/vilma2` links the stack's.

## libvilma.a link test (job 28542741, COMPLETED, 26 s)

vilma2's own tests, built against `libvilma.a` with the host link line (gate flags, omp): `test_params`, `test_host_terms`,
`test_host_bitident`, `test_restart_host`, `test_remap` (Gauss grid), `test_coupling`, `test_io_prec`: all PASS.
`test_visc_fesmdata` (FesmData Bagge file, structure 6, lmax 170, read only) reads the file but FAILS its 1e-6 dex
check: max |FesmData - reference| = 3.9e-3 dex (`Bagge2021_visc.nc`) and 7.8e-3 dex (`GLOBAL-T170_SolidEarth-Bagge2021.nc`).
The reference is vilma2's vendored `input/bagge2021.nc`: `isostasy_data` has no `bagge2021.nc` (the file the test names).
Not a link problem; to settle which file the coupled run uses (below).

## G0: identity (`--restarts`, `julia +1.13`, `compare_g0.sh`)

| comparison | jobs | grl | ant |
|---|---|---|---|
| K(S5, vilma=0, ee1f91d) cont vs K(S4, f81e004) cont | 28542775 / 28542777 | PASS, 1249 vars | PASS, 1341 vars |
| R2(S5) cont vs R2(S4) cont (yelmox 5ac89c8) | 28542778 / 28542779 | PASS, 1274 vars | PASS, 1368 vars |

0 failing, 0 absent, 0 new, first divergence none in all four. main ee1f91d differs from f81e004 only in `docs/`.
Runs: `runs/{grl,ant}/kx-ee1f91d/S5/cont` (maps `$G/maps/S5`), `ref-S5/{grl,ant}/yx-5ac89c8/S5/cont`
(a new root; `ref/` untouched). Reports: `$G/reports/P4WP0/`.

## Maps

`$G/maps/S5` = copy of `maps/S4` (44 files, read-only). The G0 runs used it. **No new pair made**; pending:

| pair (plan §13.1) | status | needs |
|---|---|---|
| ANT-16KM <-> ANT-32KM (bilin both ways, con 16->32, nn) | exists (from S4) | - |
| GLOBAL-1DEG -> NH-32KM (bilin) | exists (from S4); not a phase 4 pair | - |
| GAUSS-L170 <-> NH-32KM, ANT-32KM, ANT-16KM (con + bilin) | pending | `grid_GAUSS-L170.txt` is written by kryos at `geo_glob` init (§3.5, WP2/WP4) and the pairs primed by `wire_prime`; fesm-utils reads a Gaussian description (explicit `yvals`) and its coupler builds the map, but no kryos run asks for it yet |
| GLOBAL-0.5DEG -> GAUSS-L170, NH-32KM, ANT-32KM (con) | pending | `[globe_ice]` datasrc on GLOBAL-0.5DEG (WP3); no existing config reads that grid into these |

Adding them: a new folder (S5 copy + new maps), or reopen S5 (`chmod -R u+w`) before the first run that needs them.

## Data check (read only, no copies, no links)

Gate links are `fesmdata`, `ice_data`, `isostasy_data` (-> `/work/ba1442/<name>`); `fesmdata-prepared` is **not** linked.

| file | exists | notes |
|---|---|---|
| `fesmdata/Global/GLOBAL-0.5DEG/GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc` (59 MB) | yes, via gate link | `grid_name = GLOBAL-0.5DEG`, 720 x 360, lon -179.75..179.75, lat -89.75..89.75 (ascending); `time` = 80000..0, step 2500, 33 records, units "years ago", **0 ka is the last record** (index 33); `z_bed`, `H_ice` (also `z_topo`, `dz_sl`, `f_ice`, `f_valid`) float, `_FillValue -9999`, **no missing cells**; 0 ka: `z_bed` -8614..5636 m (mean -2120), `H_ice` 0..4404 m (mean 219) |
| `/work/ba1442/fesmdata-prepared/Bagge2021_visc/Bagge2021_visc.nc` (341 MB) | yes, absolute path only (not under a gate link) | `log10_visc(structure, depth, lat, lon)` float, NaN fill; lon 512 x lat 256 (Gaussian lats), depth 164, structure 8, structure 6 = `v_1.0_s16`; readable by vilma2 at lmax 170 (link test) |
| `fesmdata/Global/GLOBAL-T170/GLOBAL-T170_SolidEarth-Bagge2021.nc` (63 MB) | yes, via gate link | the path in plan §3.4; same layout, bit-rounded (0.05%), `grid_name = GLOBAL-T170`; also readable; differs from the prepared file by ~4e-3 dex at structure 6 |
| `GLOBAL-T170` grid description `fesmdata/Global/GLOBAL-T170/grid_GLOBAL-T170.txt` | yes | 512 x 256, not VILMA2's 680 x 342 (plan F1) |
| VILMA2 `input/vilma_defaults.nml` (12 KB) | yes (`stacks/S5/vilma2/input`) | not yet in kryos `input/` (WP2: `sync.txt` line) |
| VILMA2 `input/vilma-variables.md` (5.6 KB, `prec` column at 5e5a9ee) | yes (vilma2) | read at run time by `vilma_io` (default `input/vilma-variables.md`): **missing from the planned `sync.txt` lines** (§3.3 lists only `vilma_defaults.nml`); WP2 needs it or the table path |
| VILMA2 `input/bagge2021.nc` (default `visc_3d_file`, vendored) | yes (vilma2) | not needed when `visc_3d_file` is the FesmData file with `visc_3d_structure = 6`; used by `test_coupling` |
| `isostasy_data/earth_structure/viscosity/bagge2021.nc` | **missing** | only `pan2022.nc`, `wiens2022.nc` there; `test_visc_fesmdata` names it |
| earth structure (PREM, `n_layer = 11`) | built in | no file |

Not checked: G5's "21 -> 0 ka window of the ANT base forcing" (phase4-plan §13.2 G5 says WP0's data check settles it; not in this task).
