# Phase 4 WP1: contracts (2026-10-29)

Branch `p4-wp1` (from main ee1f91d, never pushed). Mac, gfortran, `yelmo=1 fastisostasy=1`, `make check` at every commit.

| commit | what |
|---|---|
| 32f755b | geo `dphi` export (fastisostasy `dphi = bsl`: init and `geo_scatter`); row `geo.dz_ss -> ice.z_sl` adds `dphi` |
| 43fd58b | ice static exports `footprint`, `H_ice_pd` (yelmo: `mask_ice /= MASK_ICE_NONE`, copy of `dta%pd%H_ice`; put once at `PASS_INIT`) |
| 03b840b | ice `flotation` binding, `ice_flotation` wrapper, `caps%has_flotation`; yelmo: `calc_H_grnd` + `gl_sep` subgrid rule |
| 434533a | `geo_opt_t%host`, `%port`; rank-2 declaration of a host global part; `component_core_filled` |
| fd8a417 | rows: globe -> port (`dz_bed`, `dz_ss` bilin, `dphi`), six `T_GIA` rows; `wire_lag_status` over every `T_GIA` source |
| 59cfb6e | `wire_resolve`: an import the core fills and no row names is no gap, logged; `wire_node%core` |

## Bit-neutral: `test_run_sequence` `plain` against main's build

Baseline: the `plain` runs (cont 0 -> 80, chunk 0 -> 40, 40 -> 80) of this tree at ee1f91d, copied before the first commit.
After each commit the same runs of that commit's build, `kryos-compare --restarts` (outputs and bundles), identity map.

| commit | cont | chunk c00 | chunk c01 | new variables (N) |
|---|---|---|---|---|
| 32f755b | PASS, 1151 vars | PASS, 620 | PASS, 620 | 0 |
| 43fd58b | PASS, 1151 | PASS, 620 | PASS, 620 | 2 (`footprint`, `H_ice_pd`) |
| 03b840b | PASS, 1151 | PASS, 620 | PASS, 620 | 2 |
| 434533a | PASS, 1151 | PASS, 620 | PASS, 620 | 2 |
| fd8a417 | PASS, 1151 | PASS, 620 | PASS, 620 | 2 |
| 59cfb6e | PASS, 1151 | PASS, 620 | PASS, 620 | 2 |

0 failing, 0 absent in all 18. The chunk == cont check of `test_run_sequence` passes at every commit. At 59cfb6e also
`make check yelmo=0 fastisostasy=0` (stubs; the new tests need no backend).

## Tests

- `tests/test_ice_contract.f90` (new): `footprint`/`H_ice_pd` declared by option only; toy ice with a flotation binding through
  the wrappers; stops without `has_flotation` and off the instance grid.
- `tests/test_geo_contract.f90` (new): regional, port and host-global declarations; `geo_scatter` (`dphi` for a barystatic geo,
  not for a port); a host part without a grid stops; `wire_node` carries the core-filled names.
- `tests/test_wiring.f90`: table (59 rows, `T_GIA`, globe -> port rows), `tc_gia` (rows, handles, no map and no workspace, never
  an exchange row, `wire_exchange` refuses, lag over every source, inactive row naming a missing source, the host-import rule and
  its negatives), abort cases `gia_exchange`, `host_unmarked`.
- `tests/test_run_sequence.f90` (Yelmo, grl): `dphi == bsl` and `z_sl = dphi + dz_ss`; `footprint == (mask_ice /= NONE)`,
  `H_ice_pd == dta%pd%H_ice`, put once (t1 = 0 after the steps; present after a restore); `ice_flotation` == `tpo%now%f_grnd` bit for
  bit (0 cells differ) and == the exported `f_grnd`, at t0 and t = 20.
- `test_geo_fastisostasy`: `dphi == bsl` at the cold start and every step; `test_output`: the tables hold `dphi`, `footprint`, `H_ice_pd`.

Not run here: the flotation test on ant and nor (Levante, WP7/G1); `grl` only locally.
