# Phase 4 WP7c: G2c, G5, G6, G7 on stack S8 (2026-10-30); G2b not run

`$G = /work/ba1442/robinson/kryos-gate`. Stack `$G/stacks/S8` as built (kryos e4453dc, vilma2 4419636, fesm-utils a13cd38; S9 not touched). Branch `p4-wp7c` (main 624f686, not pushed).
Scripts `tools/levante/{g2c,g5,g6,g7,run_vilma2,jl,cmp_g6}.sh`, `tools/levante/{g2c_report,g5_report,g5_iv,g5_spike,g6_report,g7_report,g7_cadence}.jl`, `tools/bin/{kryos-ant-par,kryos-smb-par}.jl`,
par files `par/gate/{g5-vilma2,g5-fastiso,g6-smb}.nml`, gates `tools/gate/{g6,g6-smb}.toml`. Runs `$G/runs/{g2c,g5,g5mg,g6,g6-smb,g7}/kx-e4453dc/S8/`, full reports `$G/reports/WP7c/`.
All runs: 64 threads, exclusive node (G5 fastiso 32 threads), `-O2 -fp-model precise`, S8 binaries. Used about 9.7 node-hours in all (G2c 4.1, G5 3.0, G6 1.7, G7 0.9).
`run_vilma2.sh`: only change is `SINGLE=1` (any run name is one run like `cont`).

## G2c: spin-up (for Alex: `equil_time_max`, `equil_rate_tol`)

Science cold start from -80 ka (gate `g2c`: the science par, `kryos.write_lossless`, 2D and the `gia`/`vilma` streams every 10 yr, 100 yr), `equil_time_max` = 0, 2e3, 5e3, 1e4, 2e4 yr, each with `equil_rate_tol` 1e-3 and 1e-5 m/yr. `equil_rate_tol` is not in the par file: added to `[vilma]` in a copy.
Reference for the differences: the most converged run (2e4, 1e-5). Host region: `w_ref` > 0.5 (the ice footprints, 20.8 % of the area); means are cos(lat)-weighted. Wall: job elapsed (init 10.5 min is in all), node-h = wall/60.

| tol [m/yr] | cap [yr] | 1-D phase | full phase (relaxed yr, passes, sub-steps) | stop | seed `dz_bed`, host mean (max) [m] | seed diff to reference, host [m] (share) | rate of interval 1 diff, host [mm/yr] | drift 5 intervals, host mean [m] | wall [min] (node-h) |
|---|---|---|---|---|---|---|---|---|---|
| 1e-3 | 0 | 1 pass, 200 yr, v 4.1e-4 | - | 1-D converged | 1.62 (57) | 10.4 (94 %) | 0.465 | 0.213 | 11.6 (0.19) |
| 1e-3 | 2e3, 5e3, 1e4, 2e4 | same | 519, 2, 46; v 5.5e-4 (0.55 tol) | converged | 2.39 (135) | 9.3 (84 %) | 0.461 | 0.095 | 12.8-13.1 (0.22), results bit-identical |
| 1e-5 | 0 | 12 passes, 93.5 kyr, v 6.7e-6 | - | 1-D converged | 10.84 (528) | 1.03 (9 %) | 0.133 | 0.096 | 15.3 (0.26) |
| 1e-5 | 2e3 | same | 1994, 5, 174; v 3.6 tol | cap | 10.99 (564) | 0.49 (4 %) | 0.031 | 0.043 | 20.3 (0.34) |
| 1e-5 | 5e3 | same | 4990, 6, 432; v 1.8 tol | cap | 11.05 (575) | 0.28 (2 %) | 0.007 | 0.043 | 29.7 (0.50) |
| 1e-5 | 1e4 | same | 9989, 8, 863; v 0.86 tol | converged | 11.10 (582) | 0.10 (0.9 %) | 0.002 | 0.043 | 43.2 (0.72) |
| 1e-5 | 2e4 | same | 13968, 8, 1206; v 0.76 tol | converged | 11.12 (584) | reference | reference | 0.043 | 54.8 (0.91) |

Global means (same runs): seed 0.93 m (stock), 1.09 (1e-3 with the full phase), 3.58-3.69 (1e-5). The global mean rate that the criterion uses (log, last pass): 4.1e-4 m/yr (stock 1-D), 5.5e-4 (full phase), 6.7e-6 to 7.6e-6 (1e-5).
Host region next to global: the mean |`dz_bed`| and the mean rate in the host region are 1.7 to 3.0 times the global values (stock: rates of intervals 2-5 4.8-6.0 mm/yr host against 1.4-1.7 global); the criterion at 1e-3 m/yr is 1 mm/yr globally, 2.5 to 6 mm/yr over the ice.
Rates of interval 1 are those of the relaxation under the held t0 load (lag 1: no load change in the first call), rates 2-5 contain the response to the changing ice; the drift (`dz_bed` at the 5th interval against the seed) has the load response in it, which is why it plateaus at 0.043 m (host) from the 2e3 run on.

Findings:
1. At 1e-3 the cap does nothing: the 1-D phase stops after its first pass (200 yr relaxed, mean rate 4.1e-4) and the full phase after two (519 yr, 5.5e-4). Results for 2e3 to 2e4 are bit-identical. The stock state holds about 15 % of the host mean deflection of the converged state (1.6 m of 11.1 m); the mean absolute difference of the seed fields is 94 % of the converged mean (84 % with the full phase). The column "share" below is this difference over the mean |`dz_bed`| of the reference seed. The rest relaxes during the run (not measured beyond 100 yr).
2. The tolerance is what matters. 1e-5 in the 1-D phase alone (93.5 kyr of relaxation, 12 passes, +3.7 min) brings the seed to 9 % of the converged one; the full phase adds the 3-D structure: 4 % at 2e3 (+8.7 min), 2 % at 5e3 (+18 min), 0.9 % at 1e4 (+32 min, converged by the criterion just below the cap), reference at 2e4 (converged after 14.0 kyr, +43 min).
3. Cost of the full phase is about 115 cpu-s per sub-step (1.8 s wall at 64 threads) with a mean sub-step of 11.5 yr: 1200 sub-steps for 14 kyr. All jobs are far below 8 h (largest 55 min, 0.91 node-h); nothing was skipped.
4. The 20 ka LGM start of G5 (below) behaves differently at the stock tolerance: the 1-D phase needs 7 passes (8.6 kyr relaxed; pass 1 has mean rate 0.47 m/yr) before the 1e-3 criterion is met. The rate criterion is a good stop for a load far from equilibrium and too loose for the mild -80 ka load.
5. Equilibrium under the t0 load is the assumption of decision 7; nothing here says it is realistic at -80 ka (the NH was not in isostatic equilibrium then). The numbers say what each setting does to the seed.

Pick proposed for review: `equil_rate_tol = 1e-5`; `equil_time_max` 1e4 (seed within 1 %, +32 min once per cold start) or 2e3 (within 4 %, +9 min). The 80 kyr run costs about 14 node-h (G7), so either is small. G6 below ran with the placeholder (0, 1e-3); it takes 25 min to rerun with the choice.

## G5: ANT deglaciation (ANT only, VILMA2 global, NH as external PaleoMIST)

Par `par/gate/g5-vilma2.nml` (`tools/bin/kryos-ant-par.jl` from the science par): `domains = "Antarctica"`, clock -20e3 to 0 (`method = "rel"`), alpha series, `[recon_ant]` PaleoMIST, `[Antarctica] experiments = "recon_ant"`; stock spin-up (placeholders). NH ice is the external `[globe_ice]` history (ESL -97.3 m against the reference at the start).
Window: 20 to 0 ka. The ice history has records every 2.5 kyr and the recon takes the nearest one, so a start at 21 ka is the 20 ka slice 1 kyr early. The ANT base (`yelmox_Antarctica`, and `kryos_Antarctica.nml`) has **no time-dependent forcing**: constant RACMO and ISMIP7 snapshots, `coefs = ""`, `terms = ""`. The paleo forcing is the documented `-p` of yelmox `scripts/ant-paleo/params.sh`: coefs alpha, dTo and terms lgm, ocean (alpha series covers 125 kyr, 0.9 at 21 ka), with DJF as the second PMIP3 temperature; I added the terms, `[pmip3_mean_tas_djf]` and the DJF group in `[snapshot_picontrol]`, `[snapshot_lgm]` (North is not active). This needed a new map GLOBAL-1DEG -> ANT-32KM (bilin), made by a mapgen run (`$G/maps/gen/S6-g5mg-e4453dc`, frozen; the other maps are S6's).
The yelmox ant-paleo setup also changes the physics (32 km tuning, optimised `cb_ref` and `tf_corr` from a spin-up bundle): not used; the Yelmo configuration is that of the science case.

Run `vilma2-full`: 2000 steps, no stop, 101 min (1.7 node-h), 9.5 GB, 2.95 s per step (VILMA2 and the gather 2.15, ice 0.34, output 0.19). The same with FastIsostasy: `fastiso-full` (`[barysealevel] method = "fastiso"`: sea level from the ANT ice only), `fastiso-sealevel-file` (method "file", `input/sl_deglaciation.dat`, the global deglaciation curve), 36-38 min.

**The run is not a deglaciation.** Antarctica gains ice in all three runs: V_sle 67.9 m at the start to 95.7 m (VILMA2), 84.3 m (FI, own sea level), 92.0 m (FI, file). The gathered ANT grounded ice goes from -14.0 m to -45.5 m ESL against the present-day reference (+31.4 m), the sea level rises by 67 m (VILMA2: `bsl` -108.2 m to -41.3 m). The unoptimised ANT base (relax, opt off) is not in equilibrium under present-day or glacial forcing (K3-sci: +5 m per kyr, G6: +2.6 m per kyr). So "uplift under loss" is exercised only by local retreat (about 650-700 cells lost more than 100 m at the end); the gate as specified has no large ice loss to test against. A supplementary test with prescribed ice loss (a mock history, as G2a for growth) would settle the sign; it is part of what G2b needs.

(a) Sign and size of `dz_bed` (ANT port, change since the start; `g5.md`): cells that gained more than 100 m subsided in 97.8 to 99.6 % of cases (cumulative) (mean -7.5 m at 500 yr, -144 m at 0 ka for +848 m of ice); on the start ice mask the mean change is -4.7 m (500 yr) to -134 m (0 ka), minimum -430 m, maximum +71 m; `dz_ss` below 43 m everywhere. Cumulatively, lost cells mostly subside (uplift in 0 to 7 %: neighbouring ice gain and the ocean load of the rising sea level are in it, not separated). Per 500-yr interval the regression of the bed change on the ice change has the expected sign in 38 of 40 intervals (slope -0.01 to -0.17 m/m, R -0.3 to -0.7; positive and weak (0.005, R 0.1) only in the last two intervals, 3 and 13 lost cells); lost cells subside before -12 ka and **uplift in 46 to 74 % of cells from -12 to -9 ka and 66 to 100 % from -8.5 ka on** (+0.4 to +2.5 m per 500 yr, while few cells are left after -2.5 ka); before -12 ka the sea level rises by 3 to 6 m per 500 yr.
(b) `bsl` against the ESL (`g5.md` (b)): ESL ANT + external (gather at t - dt) -113.9 m at -18 ka, `bsl_ive` -120.6 m (the ratio of the ocean areas, 3.405e14 against 3.618e14 m2, explains it), `bsl` -111.1 m: `bsl` is the ocean mean of the relative sea level and exceeds the sum by 2.9 to 5.0 m (the fill of the flooded area). ANT part: -(change of the gathered ANT ESL) = 31.4 m against Yelmo's `V_sle` change 27.8 m: the gather counts the grounded ice below flotation (3.6 m), `V_sle` only the ice above flotation.
G3 (all 2000 steps; (iii) at the 200 records of the lossless `gia` and `vilma` streams, every 100 yr; `g5-g3.md`):

| budget | result | bound |
|---|---|---|
| (i) merge closure | max 1.3e-12 | 1e-9: PASS |
| (ii) `mass_resid` | max 2.4e-15 | 1e-6: PASS |
| (iii) full balance with `fill_int` | residual / abs(int C rsl dA) max **1.6e-14**; the recomputed term equals `bsl_ive A_ocn` to 7e-8 | 1e-6: PASS |
| (iv) `bsl_ive` against the gathered ESL, as volumes | max **1.7e-2** (-7.95 ka); 1e-3 to 5e-3 while the NH is melting (-20 to -9 ka); at the last step 7.8e-6 (as sea level: 6.0 %, the ocean area) | 1e-5: **exceeded** during the deglaciation, met for the last 5 kyr |
| (v) `A_ocn_ref` | +0.233 % | 1 %: PASS |

(iv): the series (`g5-iv.md`, `g5-spike.log`) is smooth, not a glitch: an excess of `bsl_ive A_ocn` over ESL times the sea surface of 2.6e14 to 2.8e14 m3 (0.75 m of sea level) between about -8.6 and -7.4 ka, 1e13 and below before and after. (iii) holds to round-off, so VILMA2 is consistent with the terms it got: the difference is between the kryos ESL of the gathered grounded ice and VILMA2's grounded-ice increment `dIg` while the NH marine basins flood (VILMA2 documents that `bsl_ive` reads high over cells that flooded since the reference). Not investigated further. The bound of 1e-5 came from the -80 ka case (1.7e-6); for an LGM-scale deglaciation 2e-2 as volumes is what is observed.

(c) Documented differences against FastIsostasy (`g5.md` (c), same ice, climate, forcing, start; columns on the start ice mask, change since the start):

| t [ka] | V_sle VILMA2 / FI-own sea level / FI-file [m] | bsl VILMA2 / FI-own / FI-file [m] | mean `dz_bed` VILMA2 / FI-own / FI-file [m] | max `dz_bed` | min `dz_bed` |
|---|---|---|---|---|---|
| -15 | 83.9 / 76.5 / 83.9 | -93.1 / -18.2 / -98.3 | -47.6 / -98.4 / -134.8 | 12.4 / 2.6 / 23.4 | -137.9 / -403.7 / -483.6 |
| -10 | 89.6 / 80.9 / 87.9 | -54.1 / -22.6 / -40.5 | -89.2 / -129.0 / -171.5 | 21.7 / 8.8 / 21.6 | -245.2 / -452.8 / -506.0 |
| -5 | 93.0 / 83.1 / 89.6 | -37.3 / -24.7 / -2.3 | -116.1 / -146.6 / -187.5 | 36.1 / 17.0 / 22.2 | -348.4 / -483.4 / -548.3 |
| 0 | 95.7 / 84.3 / 92.0 | -41.3 / -25.9 / 0.0 | -134.3 / -157.1 / -204.4 | 53.6 / 26.0 / 15.0 | -430.2 / -520.7 / -604.0 |

- Sea level: FI with its own method sees only the ANT ice (falls by 26 m as ANT grows); VILMA2 includes the NH melt (far field) and the self-consistent rise; FI with the file curve rises by 135 m (VILMA2 108 m to 41 m, since ANT is 45 m above the reference at the end). That is the largest difference in the ice (V_sle 84.3 m against 95.7 m at 0 ka). With the same sea level history (file) the ice agrees to 0.006 m at -15 ka and differs by 1.6 m at -10 ka and 3.7 m (4 %) at 0 ka.
- Bed: FI subsides faster and deeper (mean -204 m against -134 m at 0 ka, minimum -604 m against -430 m) and uplifts less (maximum +15 m against +54 m): a regional LV-ELVA Earth of the ANT 1e21 Pa s layers against the global 3-D viscosity, and no far-field load (NH ice, ocean load) in the FI run beyond the imposed sea level; VILMA2 also has rotation. Lost-cell uplift fraction (FI-file, > 100 m loss): 1 to 9 %, as VILMA2's (0.3 to 6.6 %).
- A VILMA2 standalone replay of the recorded ANT load (plan 13.2 G5) was not run (no harness; see G2b).

## G6: science case smoke (`par/kryos_NH-ANT_vilma2.nml`, both domains, from -80 ka)

2 kyr (`cont`, 200 steps, 25 min, 0.42 node-h, 11.5 GB), restart at -79000 (chunks of 100 + 100 steps, 17.3 and 16.7 min). No stop, no warning. Streams `gia`, `vilma` lossless every 50 yr.

| check | result |
|---|---|
| G3 every step: (i) / (ii) / (iii) / (iv) as volumes / (v) | 1.2e-12 / 2.6e-15 / 1.6e-13 (41 records) / 1.7e-6 (as sea level 2.2e-3) / +0.23 %: all PASS |
| chunk (restart at the middle) vs cont, bundles included (`kryos-gate-kk.jl`) | **PASS, bit identical**, 3198 variables, 0 failing (`pc_tau_max` passes at 1 kyr chunks) |

Drift against K3-sci (the fastiso `sci` run over its 1 kyr; change since the start, `g6-drift.md`):

| | t [yr] | 100 | 250 | 500 | 750 | 1000 | 2000 |
|---|---|---|---|---|---|---|---|
| North, vilma2 | d V_sle [m] | 0.12 | 0.21 | 0.12 | -0.04 | -0.05 | 0.07 |
| North, K3-sci | d V_sle [m] | -0.85 | 0.21 | 1.66 | 2.78 | 3.74 | |
| Antarctica, vilma2 | d V_sle [m] | 0.36 | 0.83 | 1.50 | 2.09 | 2.61 | 4.23 |
| Antarctica, K3-sci | d V_sle [m] | 1.73 | 2.40 | 3.37 | 4.22 | 4.96 | |

North is stable (|d V_sle| < 0.25 m in 2 kyr, V_ice 5.1e6 at the start against 32.6e6 for the ICE-6G LGM start of K3-sci, which grows 3.7 m per kyr under constant PD forcing). Antarctica grows 2.6 m per kyr (K3-sci 5.0): the same unoptimised ANT base as in G5. `bsl` -7.75 m to -12.08 m.

`smb_simple` variant (surface of both domains from `smb_simple`): **the `-p` list of the conversion note stops** at `nml_validate`: `surface_chion`, `chion`, `orbit` are unknown keys of `[surface_nor]`, `[surface_ant]` once the backend is `smb_simple`. Fixed in `par/gate/g6-smb.nml` (`tools/bin/kryos-smb-par.jl`: the three keys out, `smb_simple = "smb_simple"`; the module defaults of `[smb_simple]` need no group). 1 kyr (`cont`, 100 steps, 16 min), restart at -79500: no stop; G3: 1.3e-12 / 1.5e-15 / 3.1e-13 / 3.0e-7 volumes (2.1e-3 sea level) / +0.23 %; chunk == cont **PASS**, 3006 variables. Drift per kyr: NH +1.07 m, ANT +0.84 m (untuned `smb_simple`; chion: -0.05 and +2.6).

## G7: cost

Science case, 400 yr from -80 ka at geo cadence dt = 10, 20, 50, 100 yr (`g7`; production-like output: 2D every 100 yr, streams off, bundles at -79800 and -79600); 64 threads, exclusive node; `g7.md`, `g7-cadence.md`.

| dt geo [yr] | calls | gather + VILMA2 + forcing [s/step] | per call [s] | ice, two domains [s/step] | output [s/step] | climate + surface + shelf | total [s/step] | wall [min] | node-h | peak RSS [GB] |
|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 40 | 2.14 | 2.14 | 1.30 | 0.86 | 0.50 | 4.82 | 13.7 | 0.23 | 11.6 |
| 20 | 20 | 2.06 | 4.12 | 1.29 | 0.78 | 0.50 | 4.64 | 13.3 | 0.22 | 11.5 |
| 50 | 8 | 1.99 | 9.97 | 1.25 | 0.80 | 0.50 | 4.55 | 13.6 | 0.23 | 11.6 |
| 100 | 4 | 1.79 | 17.88 | 1.15 | 0.77 | 0.50 | 4.22 | 13.0 | 0.22 | 11.5 |

- **The VILMA2 cost is proportional to the simulated time**: per call 0.67 s + 0.175 s/yr x dt. In the 3-D model the explicit sub-step is about 11.5 yr and costs about 1.8 s (115 cpu-s, the same as in the spin-up), so a call at dt = 10 is one sub-step plus 0.7 s of gather and solve. A coarser cadence saves only the fixed part: dt = 100 saves 17 % of the Earth's cost, 12 % of the step. The cost is VILMA2's sub-step ceiling, not the cadence.
- Error of a coarser cadence (coupled, 400 yr at -80 ka, slow loads; `g7-cadence.md`): `dz_bed` against dt = 10: max 0.16 m (dt 20), 0.26 m (50), 0.32 m (100) with rms 2.6, 4.0, 7.5 mm after 400 yr; `bsl` within 0.002 m.
- Peak memory 11.6 GB; init (two Yelmo cold starts, 12 min) and spin-up (stock 38 s cpu) are in the wall. Output of the step: 0.12 s without output, 5.8 s for a 2D output of all components.
- Restart: bundle 340 MB (`geo_glob_restart.nc` 212 MB), written in 21.3 s (the output column jumps by that at the step after the restart; the 30-step average hides it).
- Output of 400 yr with 2D and streams every 50 yr (G4 runs): rounded 175 MB for 30 files, lossless 298 MB (1.71 x); per file 1.2 to 1.7 x, the Earth's: `geo_glob_vilma.nc` 16.2 against 57.9 MB (3.6 x), `geo_glob.nc` 4.6 against 11.9 (2.6 x), `geo_glob_gia.nc` 2.3 against 8.6 (3.7 x). A 2D record of all components is about 20 MB rounded.
- ANT only (G5, 2000 steps): 2.95 s per step (ice 0.34 s), 0.85 node-h per 10 kyr.

Recommendations:
- **Q6 (geo cadence): keep 10 yr.** A longer cadence does not save time (12 % at 100 yr) and adds a lag of up to 0.3 m in `dz_bed` at -80 ka; the lever is VILMA2's explicit step ceiling, not the cadence (a VILMA2 question; with a 50-yr ceiling the Earth's share would fall by a factor of about 4: an estimate, not measured).
- **Q12 (restart schedule): a bundle every 5 kyr** (16 bundles, 5.4 GB, 340 MB and 21 s each, 6 min of writing in 15 h), plus the one that ends each job. No thinner schedule for `geo_glob` is needed (212 MB). Restart exactness is the G4/G6 result (bit identical at S8).
- Estimate for the 80 kyr science run (NH + ANT, `par/kryos_NH-ANT_vilma2.nml`, dt = 10, 8000 steps, 64 threads on an exclusive node): per step 4.8 s at the -80 ka ice (measured) up to about 7.5 s with LGM-size NH and ANT ice (K3-sci: ice 3.4 s, surface 0.8 s per step); wall **11 to 17 h**, taken as **about 14 node-h**, plus init 0.2 h and the spin-up (G2c: +0.15 to +0.7 h). The compute partition limits a job to 8 h: **at least three jobs** of about 25 to 30 kyr each, restarting from the bundle. Storage with the science par output (2D every 200 yr, small 2D every 50 yr, 1-D every 10 yr, streams off): about 9 GB rounded (15 GB lossless) + 5.4 GB of bundles (5 kyr), under 25 GB; `/work/ba1442` has 692 GB free (96 % used). Untested: 128 threads on the node (half of its cores idle at 64).

## G2b: not run

The compute budget is not the limit (a 21-0 ka replay at three lmax and the standalone: about 7 node-h; 80 kyr about 35). What is missing is a harness and three decisions: there is no mock-ice component (the `ice` backends are `yelmo` and its stub), and `vilma.x` needs the same history as a global field. Proposal: `tests/g2b_replay.f90` like `g2a_cap.f90`, which puts a PaleoMIST history (NH-32KM, ANT-32KM files) into the ice exports at every step and drives the real gather and `geo_glob` (merge, host-resolved flotation, ports), against `vilma.x` on the same history from the global 0.5 degree PaleoMIST file (VILMA2's own remap and flotation at lmax), at lmax 85, 128, 170, difference fields of `dz_bed`, `dz_ss`, `bsl`.
To decide: (1) window: 21-0 ka (or 20-0) or the whole 80 ka; (2) the standalone's input: the global PaleoMIST (what VILMA2 would do alone; the NH/ANT footprint part then differs by the grid) or the hub's merged field; (3) the same stock spin-up (placeholders) or the G2c choice. This also gives the ice-loss sign test missing from G5.

## Open for you

1. G2c: `equil_rate_tol` and `equil_time_max` (above), then rerun G6 (and G5 if wanted) with them.
2. G5: the ANT base does not deglaciate (grows 28 m of sea level). Options: accept G5 as run (budgets, signs for growth, FastIsostasy differences) and test loss with the G2b mock history; or tune the ANT configuration (the yelmox ant-paleo `opt` spin-up is not in kryos).
3. G3 (iv) as volumes for a large deglaciation: bound 1e-5 is exceeded (1.7e-2 at -8 ka; the budget (iii) holds at 1.6e-14): relax the bound, or find the definition difference (the flooding of marine cells) with VILMA2.
4. The conversion note's `smb_simple` `-p` list is wrong (above); `par/gate/g6-smb.nml` replaces it.
5. G2b design (above), Q6 and Q12 (above).
6. The ANT paleo terms (lgm, ocean, DJF) are in the G5 par only; the science par's ANT has constant present-day climate (alpha acts on the NH only).
