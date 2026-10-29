# Phase 4 WP1: contracts (2026-10-29)

Branch `p4-wp1`, rebased onto main 79461c5 (WP0), never pushed. Mac, gfortran, `yelmo=1 fastisostasy=1`, `make check` at every commit.

| commit | what |
|---|---|
| a71a921 | geo `dphi` export (fastisostasy `dphi = bsl`: init and `geo_scatter`); row `geo.dz_ss -> ice.z_sl` adds `dphi` |
| 8287e16 | ice static exports `footprint`, `H_ice_pd` (yelmo: `mask_ice /= MASK_ICE_NONE`, copy of `dta%pd%H_ice`; put once at `PASS_INIT`) |
| 0581f7d | ice `flotation` binding, `ice_flotation` wrapper, `caps%has_flotation`; yelmo: `calc_H_grnd` + `gl_sep` subgrid rule |
| bd844f4 | `geo_opt_t%host`, `%port`; rank-2 declaration of a host global part; `component_core_filled` |
| c0ad7e5 | rows: globe -> port (`dz_bed`, `dz_ss` bilin, `dphi`), six `T_GIA` rows (extra sources in `wire_t%srcs`); `wire_lag_status` over every `T_GIA` source |
| bb6d753 | `wire_resolve`: an import the core fills and no row names is no gap, logged; `wire_node%core` |
| d6236b4 | `T_GIA` extra sources as an explicit list `wire_t%srcs`, not `aux` (review change) |

## Bit-neutral: `test_run_sequence` `plain` against main's build

Baseline: the `plain` runs (cont 0 -> 80, chunk 0 -> 40, 40 -> 80) of a build of ee1f91d, copied before the first commit;
the table is from the first six commits on that base. After the review change and the rebase, head d6236b4 is compared with
a fresh build of main 79461c5 (same `src/`, `tests/`, `input/` as ee1f91d): see the last row.
After each commit the same runs of that commit's build, `kryos-compare --restarts` (outputs and bundles), identity map.

| commit | cont | chunk c00 | chunk c01 | new variables (N) |
|---|---|---|---|---|
| a71a921 | PASS, 1151 vars | PASS, 620 | PASS, 620 | 0 |
| 8287e16 | PASS, 1151 | PASS, 620 | PASS, 620 | 2 (`footprint`, `H_ice_pd`) |
| 0581f7d | PASS, 1151 | PASS, 620 | PASS, 620 | 2 |
| bd844f4 | PASS, 1151 | PASS, 620 | PASS, 620 | 2 |
| c0ad7e5 | PASS, 1151 | PASS, 620 | PASS, 620 | 2 |
| bb6d753 | PASS, 1151 | PASS, 620 | PASS, 620 | 2 |

0 failing, 0 absent in all 18. (Hashes are those after the rebase; the trees are unchanged by it.)


Head d6236b4 (after the review change and the rebase) against a fresh build of main 79461c5: cont PASS, 1151 vars;
chunk c00 PASS, 620; c01 PASS, 620; 0 failing, 0 absent, 2 new. `make check` passes. The chunk == cont check of `test_run_sequence` passes at every commit. At bb6d753 also
`make check yelmo=0 fastisostasy=0` (stubs; the new tests need no backend).

## Tests

- `tests/test_ice_contract.f90` (new): `footprint`/`H_ice_pd` declared by option only; toy ice with a flotation binding through
  the wrappers; stops without `has_flotation` and off the instance grid.
- `tests/test_geo_contract.f90` (new): regional, port and host-global declarations; `geo_scatter` (`dphi` for a barystatic geo,
  not for a port); a host part without a grid stops; `wire_node` carries the core-filled names.
- `tests/test_wiring.f90`: table (59 rows, `T_GIA` with `srcs`, globe -> port rows), `tc_gia` (rows, handles, no map and no workspace, never
  an exchange row, `wire_exchange` refuses, lag over every source, inactive row naming a missing source, the host-import rule and
  its negatives), abort cases `gia_exchange`, `host_unmarked`.
- `tests/test_run_sequence.f90` (Yelmo, grl): `dphi == bsl` and `z_sl = dphi + dz_ss`; `footprint == (mask_ice /= NONE)`,
  `H_ice_pd == dta%pd%H_ice`, put once (t1 = 0 after the steps; present after a restore); `ice_flotation` == `tpo%now%f_grnd` bit for
  bit (0 cells differ) and == the exported `f_grnd`, at t0 and t = 20.
- `test_geo_fastisostasy`: `dphi == bsl` at the cold start and every step; `test_output`: the tables hold `dphi`, `footprint`, `H_ice_pd`.

Not run here: the flotation test on ant and nor (Levante, WP7/G1); `grl` only locally.
