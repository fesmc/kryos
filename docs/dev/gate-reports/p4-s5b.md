# Stack S5b (phase 4): identity with S5

S5b = S5 with fesm-utils `phase0-kryos` 3ba56a9 (was 5662899) and vilma2 main 2676879 (was 5e5a9ee), every dependent package
rebuilt (`tools/setup_S5b.sh`; yelmox topo-prep 5ac89c8 and the other packages as in S5, `stack.lock` in the stack folder).

fesm-utils 5662899..3ba56a9: optional sphere radius for `grid_init` / `grid_init_gaussian` (areas, signature; default unchanged), time
units `years ago` / `years before present` equivalent to `years BP`, `coupler_merge_end` logs coverage only above 1 + 1e-9.

## Identity, `cont`, `kryos-compare --restarts` (`tools/levante/gate_S5b.sh`, `compare_S5b.sh`)

K = kryos 5a663a3 (`vilma=0`), built on S5b in its own checkout `stacks/S5b/kryos`; maps S5 (no map file added). R2 = yelmox topo-prep, S5b build.
Runs: K grl 28552184, K ant 28552189, R2 grl 28552188, R2 ant 28552190; compare job 28552538.

| comparison | grl | ant |
|---|---|---|
| K(S5b) vs K(S5) | PASS, 1270 variables | PASS, 1362 |
| R2(S5b) vs R2(S5) | PASS, 1274 | PASS, 1368 |

0 failing, 0 absent, 0 new, no first divergence: the radius, time-unit and log changes change no existing result.

vilma2 2676879 host-API link test (`vilma_link_test.sh`, job 28552215): test_params, host_terms, host_bitident, restart_host, remap, coupling,
io_prec pass; `test_visc_fesmdata` fails only on the script's reference file `input/bagge2021.nc`, which vilma2 main no longer vendors
(the script's `REFV` default needs updating).
