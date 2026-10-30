# Phase 4 G2b: PaleoMIST replay, kryos against the standalone VILMA2 (diagnosis on S9, gate on S11, 2026-10-30)

`$G = /work/ba1442/robinson/kryos-gate`. Harness `tests/g2b_replay.f90`, scripts `tools/levante/{g2b.sh,g2b_compare.jl,g2b_diag.jl,g2b_maps.jl,g2b_prep.jl,g2b_visc_mirror.sh}`
(branch `p4-g2b`, merged into `vilma-orient`). Runs `$G/runs/g2b/kx-3e263ed/S9/` (diagnosis), `$G/runs/g2b/kx-3607918/S11/` (gate); reports `$G/p4g2b/`, `$G/reports/S11/g2b*.md`.

## Set-up

- kryos: a mock ice (Yelmo's flotation rule, `gl_sep = 1`) replays the PaleoMIST thickness on NH-32KM and ANT-32KM through the globe of the
  run (gather, host-resolved terms, `geo_glob` VILMA2, bilinear scatter into the ports), lag 1; -21 ka to 0, 2100 steps of 10 yr.
- Standalone: `vilma.x` on the global 0.5 degree PaleoMIST history (`g2b_prep.jl`, lagged by one step), `remap_input = .true.`,
  `i_eq = 1` (present-day reference of the same product), the same `&vilma` and densities.
- Both: lmax 85, 128, 170 (GAUSS-L<lmax>), Bagge2021 3-D viscosity (structure 6), spin-up of the science par (`equil_rate_tol` 1e-5,
  `equil_time_max` 2e4). Differences kryos - standalone on the Gauss grid; regions: host (kryos `host_w` > 0.5) NH and ANT, rim
  (0 < w <= 0.5), outside (w = 0). The two differ by design in the host region: kryos's load is the hub's (Yelmo flotation on the ice
  grids), the standalone's is VILMA2's own remap and flotation of the 0.5 degree field.

## Diagnosis on S9 (before the row-order fix)

| | lmax 85 | 128 | 170 |
|---|---|---|---|
| `dz_bed` host NH, max rms [m] (at 0 ka) | 51.30 | 51.29 | 51.31 |
| `dz_bed` host ANT, max rms [m] | 12.33 | 12.37 | 12.42 |
| `bsl` max abs difference [m] (/ range) | 0.95 (8.6e-3) | 0.92 (8.3e-3) | 0.91 (8.3e-3) |

- The difference did not fall with lmax (ratios 1.00) and grew through the deglaciation (host NH rms 9 m at -21 ka, 34 m at -10 ka, 51 m at 0 ka).
- The load difference did not explain it (`g2b_diag.jl`, l085): the local isostatic estimate of the load difference, -(rho_w/rho_m) dL, is
  5 to 11 m rms in the host NH; the time-integrated load difference correlates with the `dz_bed` difference at 0 ka with -0.74 (NH), 0.06 (ANT).
- Without spin-up (both runs): the same growth after the seed (host NH 52 m rms at -15 ka).
- With the 3-D viscosity mirrored north-south (`g2b_visc_mirror.sh`, the data reversed along lat in a fesmdata view used by kryos only):
  host NH max rms **4.2 m** (1.5 m at 0 ka), host ANT 2.8 m, `bsl` max 0.46 m, end 0.03 m; without spin-up 0.2 to 0.8 m after the seed.
- Cause: kryos's Gauss grid ascends in latitude, VILMA2's state is north first; the passthrough (`grid_is_gauss`) copied the rows
  unreversed both ways, so the whole kryos side was mirrored consistently, except what VILMA2 reads from files in real geography: the
  Laurentide and Fennoscandian loads sat on the ~60 S viscosity. G2a passed because its standalone (`remap_input = .false.`) read the
  ascending Gauss file without an orientation check, mirrored the same way.

Fixed in vilma2 `host-orient` (c17d513: the passthrough reverses the rows of an ascending Gauss host, the legacy reader goes by the lat
axis) and kryos `vilma-orient` (the adapter exports VILMA2's host-grid outputs); stack S11.

## Gate on S11 (vilma2 c17d513, kryos 3607918)

l085 on S11 equals the S9 mirrored-viscosity diagnosis run to every printed digit (4.170 m, 2.752 m, 0.4624 m): the fix does what the
mirrored file did.

| lmax | `bsl`: max abs d [m] (/ range) | end d [m] | `bsl_ive` max abs d [m] | `dz_bed` host NH, max rms [m] (/ max abs) | host ANT | rim | outside | `rsl` host NH | host ANT | C host NH: max int dC / ocean area |
|---|---|---|---|---|---|---|---|---|---|---|
| 85 | 0.462 (4.2e-3) | 0.032 | 1.60 | 4.17 (4.7e-3) | 2.75 (1.0e-2) | 0.70 | 0.31 | 4.24 | 2.83 | 0.142 |
| 128 | 0.394 (3.6e-3) | 0.003 | 1.34 | 3.95 (4.5e-3) | 2.09 (7.8e-3) | 0.62 | 0.31 | 4.02 | 2.18 | 0.145 |
| 170 | 0.398 (3.6e-3) | 0.017 | 1.16 | 3.89 (4.4e-3) | 1.71 (6.4e-3) | 0.62 | 0.31 | 3.96 | 1.80 | 0.148 |

Ratio to the next finer lmax (> 1: the difference falls): `dz_bed` host NH 1.05, 1.02; host ANT 1.32, 1.22; `bsl` 1.17, 0.99.

- The largest difference is at the seed (-21 ka) and falls through the run (host NH rms 4.2 m at -21 ka, 2.0 at -10 ka, 1.5 at 0 ka; l085).
  At the seed the two ocean functions differ most (host NH 3349e9 m2, l085): the kryos seed gather sees z_sl = 0 and the present-day bed
  (land on the present-day shelf in the standalone run, 2764e9 m2 of it). The local isostatic estimate of the load difference (9.3 m rms at
  -21 ka, 2.0 m at 0 ka) bounds the `dz_bed` difference; the correlation with the time-integrated load difference is now -0.29 (NH), -0.35 (ANT).
- `bsl` differs by at most 0.4 m (0.36 % of the 110 m range), at the seed; 0.003 to 0.03 m at the end.
- Sign test (ice grids, 500-yr intervals, loss > 100 m): uplift in median 100 % (NH) and 93-94 % (ANT) of the lost cells in kryos, as the
  standalone (100 %, 91-93 %); the mean bed change over them within 3 % (NH 10.3 m, ANT 2.0 m). Gained cells subside (100 %).
- Cost (64 threads, exclusive node): kryos 16 / 30 / 60 min, standalone 14 / 30 / 61 min (lmax 85 / 128 / 170), with the spin-up
  (6 min at l085).

## Thresholds (plan 13.2; set 2026-10-31)

| quantity | bound | observed (S11, lmax 85-170) |
|---|---|---|
| `bsl` max abs difference / range of the standalone `bsl` | 1 % | 0.36-0.42 % |
| `dz_bed` and `rsl`, host regions: max over records of the rms / max abs `dz_bed` of the standalone | 2 % | NH 0.44-0.47 %, ANT 0.64-1.0 % |
| `dz_bed`, rim and outside: the same | 1 % | 0.17-0.58 % |
| convergence: max rms of `dz_bed` (host NH, ANT) at lmax 170 against lmax 85 | not larger | 0.93, 0.62 |
| sign test, lost cells: median uplift fraction in kryos | >= 0.9 and within 0.05 of the standalone | 0.93-1.00, within 0.03 |

The S9 numbers (5.8 % host NH, flat with lmax) would have failed the second and the fourth bound. S11 passes all five.
The seed (z_sl = 0 on the present-day bed in the kryos seed gather) stays: its difference decays and is inside the bounds.
