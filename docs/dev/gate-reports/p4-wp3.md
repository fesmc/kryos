# Phase 4 WP3: `kryos_globe_ice` (2026-10-29)

Branch `p4-wp3`, rebased onto main a6e91bf, never pushed. Mac, gfortran, `make check` at every commit
(during development `make -o input-check check` was used for a while: the shared `vilma2` checkout had moved ahead of the
synchronised `vilma_defaults.nml`; main a6e91bf re-synced it).

| commit | what |
|---|---|
| aefc529 | `[globe_ice]` template; `globe_ice_init` (grids, datasrc of the globe); `globe_ice_add_domain` (true areas, footprint, coverage) |
| 04e0c42 | reference (§5.1): X, host terms, spliced reference, ocean area |
| 357af89 | gather of a step (§4.1): external ice, host terms, checks and diagnostics (§4.5) |
| 398b144 | one reference definition (H_ice_pd in the footprint, X elsewhere) for VILMA2 and the hub; rim metrics documented (code, plan 4.4, 5.1) |
| 1695fb2 | owner weights for overlapping domains (priority, `coupler_owner_reduce`) |
| f082dec | hub output `globe_ice.nc` (background grid, table `kryos-variables-globe_ice.md` with `prec`), hub conservation check |
| 923af09 | `KryosTools.GlobeIce` offline merge, `tools/bin/kryos-globe-ice.jl`; the test compares `globe_ice.nc` with it |

## Final state (rebased onto a6e91bf)

The history's time units are read as they are (fesm-utils 60aea73 decodes "years ago"); the test file has `years ago`, the par groups
no `time_units`. Full `make check` (with input-check) passes on the head against shared fesm-utils 3ba56a9; the `plain` runs compared with
`kryos-compare --restarts` against a fresh a6e91bf build: cont PASS 1171 variables, c00 and c01 PASS 632, 0 failing, 0 absent, 0 new.

## Bit-neutral

Baseline: the `plain` runs of `test_run_sequence` (cont, c00, c01) of a build of aeee5dd. After each commit the same runs,
`kryos-compare --restarts`: PASS, 1159 / 625 / 625 variables, 0 failing, 0 absent, 0 new (the last three commits also after the stack
rebuild below). Not run locally: ant (Levante, WP7).

## G1b pieces tested locally (`tests/test_globe_ice.f90`)

- Hub conservation: the spliced reference and the merged ice close to 1e-15 (stop above 1e-12), logged at the reference and at every write.
- Rounded file within its `prec` and not exact; lossless under `nc_set_lossless`: exact.
- Domain order swapped (registration and add order, with owner weights): `globe_ice.nc` identical (`cmp`), reference and gather bit identical.
- Offline merge (Julia, from the coupler's SCRIP maps and the model areas): lossless file equal to 1e-12 (9e-13 m on ~1e4), rounded file within its rounding.
- Not local: (1) the run with and without `[globe_ice]` (WP4 wiring), (3) `cont2` determinism, the offline merge on a real run's output (WP7).

## Environment

fesm-utils phase0-kryos moved to 3ba56a9 (sphere radius in `grid_init`, "years ago", merge log threshold) while the yelmo, FastIsostasy and
chion archives were still built against 5662899: `test_run_sequence` failed with a corrupt planet name until `make clean` and
`make deps-static` rebuilt them in the shared trees. Any worktree built before that needs the same.
