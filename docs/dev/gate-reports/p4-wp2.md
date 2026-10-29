# Phase 4 WP2: the `vilma2` geo backend (2026-10-29)

Branch `p4-wp2` on main 21320db (WP0 and WP1), never pushed. Mac, gfortran. Each commit: `make check vilma=1` (own build dir via `objdir=`/`bindir=`,
against `~/models/vilma2` main 5e5a9ee, `libvilma.a` current, fesm-utils `phase0-kryos` 5662899 with its SHTns) and
`make check vilma=0`; then the `test_run_sequence` plain runs (cont 0 -> 80, chunk 0 -> 40, 40 -> 80, grl, `yelmo=1 fastisostasy=1`) of
the `vilma=1` and of the `vilma=0` build compared with `kryos-compare --restarts` against the runs of a build of 21320db.

| commit | what |
|---|---|
| cdcd0e3 | build: `vilma=1` replaces the stop: `VILMAROOT`, `INC/LIB_VILMA`, `LIB_SHTNS` (serial or omp, before FFTW), `libvilma.a` and SHTns as file dependencies, `vilma-static`, `build_info_libs`; `vilma_defaults.nml` and `vilma-variables.md` in `input/sync.txt`; `vilma2:main` in the manifest |
| 2830b75 | `kryos_geo_vilma2.f90` (`vilma_earth_t`, `vilma_geo_t`, constraints), stub, factory (`vilma2` by scope), template `[geo_vilma2]` (`vilma`, `ice_bg`), the Gauss grid, `init_ctx_t%map_fldr`; `test_geo_vilma2`, `test_geo_vilma2_stub`, `test_build` |
| 24835d0 | `kryos_gia_check` (run-level contract checks, called by `run_init`), `test_gia_check`; the adapter logs VILMA2's gravity beside the run's |
| 44069fd | `test_geo_vilma2`: the 1D stream of the global part |

## Bit-neutral against 21320db

Baseline: the three `plain` runs of a build of main 21320db (same `src/`, `tests/`, `input/` as WP1's head), copied before the first
commit. After each commit, both builds. `kryos-compare --restarts` (outputs and bundles), identity map; 0 failing, 0 absent, 0 new.

| commit | build | cont | chunk c00 | chunk c01 |
|---|---|---|---|---|
| cdcd0e3 | vilma=1 | PASS, 1159 vars | PASS, 625 | PASS, 625 |
| cdcd0e3 | vilma=0 | PASS, 1159 vars | PASS, 625 | PASS, 625 |
| 2830b75 | vilma=1 | PASS, 1159 vars | PASS, 625 | PASS, 625 |
| 2830b75 | vilma=0 | PASS, 1159 vars | PASS, 625 | PASS, 625 |
| 24835d0 | vilma=1 | PASS, 1159 vars | PASS, 625 | PASS, 625 |
| 24835d0 | vilma=0 | PASS, 1159 vars | PASS, 625 | PASS, 625 |
| 44069fd | vilma=1 | PASS, 1159 vars | PASS, 625 | PASS, 625 |
| 44069fd | vilma=0 | PASS, 1159 vars | PASS, 625 | PASS, 625 |

Only grl runs locally (as WP1); ant, nor and `K3-sci` are Levante (WP7, G1).

## Tests

- `test_geo_vilma2` (`vilma=1`, lmax 16, 64 x 34, prescribed host terms, no ice instance): factory admission; the Gauss grid (`GAUSS-L16`, 64 x 34,
  Gaussian, longitudes from 0, cell areas sum to `4 pi R^2` to 1e-12, one cell area = `R^2 dlon w_j`, description written to the map folder);
  densities of VILMA2 = the run's `rho_ice`, `rho_sw`; reference (exports zero, `C_ocn` the reference's), the seed with and without a spin-up
  (`pre_spinup_1d`), three updates: `dz_bed`, `dz_ss` (= VILMA2's minus `dphi`), `C_ocn`, `dphi`, `bsl`, `bsl_ive`, `mass_resid`, `A_ocn` bitwise equal to
  `solid_earth_init/spinup/update` called with the same arrays (also two models in one process: OpenMP reproducible at this size);
  restart round trip (write, rebuild at the reference, read, fields, `PASS_RESTART`, continue: bitwise equal to the continued run and to the direct
  calls); the port copies imports exactly with their validity, `PASS_GEOMETRY` puts zeros and `z_bed_ref`; the 1D stream; child runs: wrong grid
  (`se%remap`), a grid file of another size, an inadmissible host term (`C > w`), an update before the build, an import not published, a radius
  mismatch, `vilma` key naming another group.
- `test_geo_vilma2_stub` (`vilma=0`): constraints unlinked; the factory stops with "rebuild kryos with `vilma=1`", global and regional.
- `test_gia_check`: no vilma2 passes untouched; a consistent run passes; stops for FastIsostasy + VILMA2 (both ways), a port without the global part,
  a domain without geo or ice, missing `has_flotation`, `footprint`, `H_ice_pd`, a port `dt` that differs.
- `test_build`: the vilma2 factory admission when linked.

## Decisions to confirm

1. **Radius.** fesm-utils' "Spherical Earth" has R = 6371221 m, VILMA2's `r_earth` is 6371000 m (3.5e-5; areas 6.9e-5). The plan says the adapter stops on a
   mismatch; it stops above 1e-4 relative and logs the difference. Alternatives: a planet of the grid at VILMA2's radius (fesm-utils change), the
   Gauss area set in the coupler (`coupler_set_area`, WP4) to VILMA2's sphere, or the tolerance as is. All budgets are homogeneous in `R^2`,
   so ratios (bsl, ESL) are unaffected; only Gt are off by 7e-5 against VILMA2's own sphere.
2. **Group name.** `vilma_par_load` looks the defaults up under the group name it reads, so key `vilma` can only be `vilma` (the adapter stops
   otherwise, as the FastIsostasy adapter pins `bsl`). A free name needs a `defaults_group` argument in VILMA2.
3. **`phys_const_compare`** (plan 3.1) does not apply: the record has no `r_earth`, and VILMA2's gravity (9.8172, from its Earth model) has no key.
   Both are logged instead (radius above, `g` at `PASS_GEOMETRY`).
4. **Port cadence.** The plan has the port take its cadence from the global geo. `gia_check` stops when the port's `dt` is set and differs; an unset
   `dt` (0, every step) is not an error, so WP4's scatter must update the port with the global's cadence whatever its own.
5. **`make check` needs a `vilma2` checkout** whatever the switch: two of the synchronised copies come from it.
6. Not in WP2: the 2D output (the default writer stops for a global grid: lon-lat writer, rows of `kryos-variables-geo.md`, WP5); footprint
   overlap (WP3, merge); `globe%cpl` registration of the Gauss grid (WP4; the adapter holds the grid in `gauss`).
