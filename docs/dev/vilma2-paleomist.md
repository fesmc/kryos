# VILMA2 standalone under the PaleoMIST ice history (Gowan et al. 2021)

Branch `vilma-paleomist`. Question: does the standalone VILMA2 (`vilma.x`, vilma2 523b660), driven by the PaleoMIST 1.0 ice history (minimal MIS 3 scenario) on Gowan's 1-D Earth, reproduce the sea level and bed of PaleoMIST; what do resolution and the 3-D Earth change; and how far is the LGM from equilibrium with the solid Earth.

- Scripts `tools/levante/pm_prep.jl` (inputs), `pm_run.sh` (sbatch, one run), `pm_compare.jl` (metrics, fields), `pm_maps.jl` (figures). Runs, inputs and comparison under `$G/runs/paleomist/` (`G=/work/ba1442/robinson/kryos-gate`).
- vilma.x of the S11 stack (`$G/stacks/S11/vilma2`, c17d513, tree identical to 523b660). Runs on the shared partition, 64 OpenMP threads.

## Setup

| | Gowan et al. (2021), SELEN | here, VILMA2 |
|---|---|---|
| sea-level equation | shoreline migration, grounding-line adjustment, rotation | same (`rotation = .true.`, flotation, sloping coast) |
| resolution | degree 256, 1 degree output | lmax 256 (and 170), Gauss grid 1024 x 514; compared on the 0.5 degree PaleoMIST grid |
| Earth | 1-D: 120 km elastic lithosphere, UM 4e20, LM 4e22 Pa s | `earth = "PREM"`, 12 layers: elastic 0-80-120 km, Maxwell 4e20 to 670 km (splits 220, 400, 600), 4e22 to the CMB (splits 771, 1428, 2084, 2741 km), fluid core; rho, mu the r^2 shell means of incompressible PREM |
| load history | 2.5 kyr steps, 80-0 ka | the 33 records as steps (`pm_force_B.nc`): record t_k enters over [t_k, t_k + 1 yr] and is held to t_{k+1} |
| before 80 ka | not stated | equilibrium with the 80 ka load: spin-up to `equil_rate_tol` 1e-5 m/yr |
| densities | not stated | VILMA2 defaults rho_ice 931, rho_water 1000 kg/m3 |
| reference | | present day = PaleoMIST 0 ka `z_bed` and `H_ice`, relaxed (`i_eq = 1`) |

- Ice: FesmData `GLOBAL-0.25DEG_IceHistory-PaleoMIST.nc` (the 0.25 degree product on the cell-centred FesmData grid), remapped conservatively to the Gauss grid by vilma.x (`remap_input`).
- Step convention: vilma.x ramps the ice linearly between two forcing slices, so an exact step needs a zero-length interval. A dt = 0 update applies the new load, but its diagnostics are those of the previous ice (KNOWN ISSUE in `vilma_timestep`), so the step is a 1-yr ramp instead; 1 yr is far below the shortest Maxwell time (100 yr in the 1-D Earth, 15 yr at the Bagge floor). Record t_k is compared with the state at t_k + 1 yr.
- Definitions (checked on the product): PaleoMIST `dz_sl` is the sea-level change from present, positive where the sea is higher than today against the land, `z_bed` is relative to the sea level of the time and z_bed(t) - z_bed(0) = -dz_sl(t) to 0.01 m. VILMA2's rsl has the same sign and z_bed = z_bed_eq - rsl, so Gowan's bed change relative to present is -RSL in both and has the RSL metrics with the bias negated. VILMA2's rsl is relative to the relaxed reference and is not zero at 0 ka (the Earth is still rebounding), while `dz_sl` is S(t) - S(0) by construction: the B runs are compared as rsl(t) - rsl(0 ka), as Gowan's. A - B compares both states against the same relaxed reference.
- Runs (`$G/runs/paleomist/`):

| run | load | Earth | lmax | spin-up cap |
|---|---|---|---|---|
| A1_l256, A1_l170 | 22.5 ka ice held, spin-up | 1-D Gowan | 256, 170 | 2e6 yr |
| B1_l256, B1_l170 | equilibrium with 80 ka, then 80 -> 0 ka | 1-D Gowan | 256, 170 | 2e6 yr |
| A3_l170 | as A1 | 3-D: science setting (PREM, 80 km lithosphere, Bagge 2021 structure 6, `pre_spinup_1d`) | 170 | 2e5 yr (3-D phase) |
| B3_l170 | as B1 | as A3 | 170 | 2e5 yr (3-D phase) |

Smoke: S_1d, S_3d (lmax 64, 25 -> 20 ka, spin-up cap 2 kyr).

## Convergence and cost

Spin-up (`solid_earth_spinup`): passes of 200 yr growing by 1.6, converged when the mean bed speed over a pass is below 1e-5 m/yr. All converged, none hit the cap.

| run | phase | passes | relaxed time | final mean speed [m/yr] | last pass: length, max bed change | wall time | memory |
|---|---|---|---|---|---|---|---|
| A1_l256, A1_l170 | 1-D | 15 | 384 kyr | 6.8e-6 | 144 kyr, 9.8 m | 10 min, 5 min | 4.3, 2.3 GB |
| B1_l256, B1_l170 | 1-D | 12 | 93 kyr | 6.6e-6 | 35 kyr, 3.6 m | 6 min, 3 min | 4.3, 2.7 GB |
| A3_l170 | 1-D pre-spin / 3-D | 15 / 12 | 384 / 93 kyr | 6.7e-6 / 7.8e-6 | 144 kyr, 19 m / 35 kyr, 2.8 m | 2 h 03 | 9.4 GB |
| B3_l170 | 1-D pre-spin / 3-D | 12 / 7 | 93 / 8.6 kyr | 7.4e-6 / 9.3e-6 | 35 kyr, 19 m / 3.4 kyr, 3.2 m | 1 h 44 | 9.5 GB |

Wall times include the spin-up and, for B, the 80 kyr transient (1-D: sub-steps of 100 yr; 3-D: 11.6 yr, 1.5 h for the 80 kyr). The LGM equilibrium needs 384 kyr of relaxation in 1-D (Gowan's lower mantle 4e22, and the lateral mean of Bagge): a 2e4 yr cap on a full-model spin-up is far from it. After the 1-D pre-spin the 3-D phase needs 93 kyr (A3) or 8.6 kyr (B3). The criterion is a mean over a pass: the last 1-D pass of A still moves the bed by up to 10 m.

## Metrics

RSL relative to present [m] on the 0.5 degree grid, area weighted: RMS, max |d| and bias = mean(d), d = first - second. Regions: glaciated = ice at some record and ice-free today (Laurentide, Fennoscandia, Barents-Kara, Patagonia, the West Antarctic embayments); ocean = present ocean outside it; global. Gowan's bed change relative to present is -RSL: same RMS and max, bias negated. Full tables: `figures/vilma2-paleomist/pm_metrics.md`, all 33 records of B1 - PaleoMIST: `pm_metrics_B1_PM.csv`.

| comparison | ka | glac RMS | glac max | glac bias | ocean RMS | ocean max | ocean bias | global RMS | global max | global bias |
|---|---|---|---|---|---|---|---|---|---|---|
| B1 l256 - PaleoMIST | 22.5 | 23.7 | 91 | -0.9 | 8.7 | 33 | 8.5 | 10.2 | 91 | 7.5 |
| | 15 | 29.0 | 106 | -6.7 | 7.5 | 22 | 7.3 | 10.2 | 106 | 6.0 |
| | 10 | 28.8 | 110 | -8.2 | 5.2 | 14 | 5.1 | 8.8 | 110 | 4.0 |
| B1 l170 - B1 l256 | 22.5 | 0.4 | 1.7 | 0.0 | 0.1 | 1.4 | 0.0 | 0.1 | 1.9 | 0.0 |
| | 15 | 0.4 | 2.0 | 0.0 | 0.1 | 1.7 | 0.0 | 0.2 | 2.0 | 0.0 |
| | 10 | 0.3 | 2.1 | 0.0 | 0.1 | 1.2 | 0.0 | 0.1 | 2.2 | 0.0 |
| B3 - B1 (l170) | 22.5 | 52.6 | 196 | -7.2 | 6.9 | 112 | -0.2 | 14.6 | 196 | -0.3 |
| | 15 | 57.1 | 230 | -13.0 | 7.5 | 118 | -0.2 | 16.1 | 230 | -0.4 |
| | 10 | 45.8 | 180 | -8.8 | 5.3 | 107 | -0.1 | 12.4 | 180 | -0.4 |
| B3 - PaleoMIST | 22.5 | 60.1 | 228 | -8.1 | 10.8 | 107 | 8.3 | 18.1 | 228 | 7.2 |
| A1 - B1 (l256), RSL | 22.5 | 146.9 | 470 | 112.9 | 11.8 | 168 | -0.2 | 37.8 | 470 | 9.0 |
| A1 - B1 (l256), u | 22.5 | 154.9 | 497 | -117.5 | 26.5 | 176 | 11.9 | 45.0 | 497 | 0.0 |
| A3 - B3 (l170), RSL | 22.5 | 166.4 | 476 | 119.9 | 12.2 | 125 | 0.0 | 42.4 | 476 | 9.6 |
| A3 - B3 (l170), u | 22.5 | 178.1 | 501 | -127.4 | 25.6 | 136 | 12.4 | 49.3 | 501 | 0.0 |
| B1 l256 at 0 ka vs relaxed reference, RSL | 0 | 96.2 | 294 | 59.6 | 8.8 | 74 | -2.3 | 25.2 | 294 | 3.2 |
| B3 at 0 ka vs relaxed reference, RSL | 0 | 102.1 | 429 | 59.0 | 8.3 | 56 | -2.2 | 26.4 | 429 | 3.3 |

A - B is against the relaxed reference (u = the radial bed displacement, dz_bed). At 22.5 ka the barystatic sea level (ocean mean of RSL relative to present) is -110.3 m in B1 and B3 against -116.3 m for the ocean mean of PaleoMIST's `dz_sl` (and -101.7 m for Gowan's ice-volume ESL at 910 kg/m3, 361e6 km2); at Barbados B1 gives -95.9 m against -105.0 m.

## Figures

![B1 against PaleoMIST](figures/vilma2-paleomist/pm_B1_vs_PM.png)

RSL relative to present of PaleoMIST and B1 (lmax 256) and their difference; black: the 22.5 ka ice margin.

![resolution and 3-D effects](figures/vilma2-paleomist/pm_effects.png)

Top: lmax 170 - 256 (B1), +-5 m. Bottom: 3-D - 1-D (B3 - B1 at lmax 170), +-30 m.

![sites](figures/vilma2-paleomist/pm_sites.png)

RSL relative to present at the sites: every output step of the B runs (the step at each record, then the relaxation over the 2.5 kyr) and the PaleoMIST records. B1 l170 lies on B1 l256.

![GMSL](figures/vilma2-paleomist/pm_gmsl.png)

VILMA2 barystatic sea level (bsl - bsl(0)) and ice-volume equivalent (bsl_ive) against PaleoMIST's ocean mean of `dz_sl` and its ESL series.

![A - B](figures/vilma2-paleomist/pm_A_minus_B.png)

Equilibrium under the 22.5 ka ice (A) minus the transient (B) at 22.5 ka, RSL and bed displacement u, +-300 m.

![present-day disequilibrium](figures/vilma2-paleomist/pm_rsl0.png)

RSL of B at 0 ka against the relaxed present-day reference.

## Findings

1. **VILMA2 reproduces PaleoMIST to about 10 m globally and 25-30 m RMS over the formerly glaciated areas** (B1, 22.5-10 ka). Fennoscandia matches within about 15 m through both glacial cycles (Gulf of Bothnia). The misfit is largest where the ice sat on deep marine beds: in the Ross and Weddell embayments VILMA2's RSL rise is about half of PaleoMIST's (up to 100 m less at 10 ka), over Hudson Bay 60 m less at 20-15 ka, and the Barents peaks are lower. VILMA2 grounds 99 % of the Antarctic LGM ice volume, so its flotation is not the cause. A marine-ice load in SELEN without the displaced water would have this pattern (the net load is the ice minus the water column it replaces); this cannot be checked without Gowan's load files. The other differences from SELEN (elastic structure, normal modes against VILMA2's explicit memory stepping) would not single out the marine beds.
2. **The far field is 5-8 m less deep in VILMA2 at the LGM** (ocean bias +8.5 m at 22.5 ka, Barbados -96 against -105 m; barystatic -110 against -116 m), while the ice-volume equivalent is -121 m. Part of it is the definition of the ocean mean (VILMA2's ocean function against PaleoMIST's z_bed < 0 and ice-free), the rest is not resolved.
3. **lmax 170 is enough.** lmax 170 differs from 256 by 0.4 m RMS and at most 2.4 m, for half the cost.
4. **The 3-D science Earth changes the deglacial RSL by 50 m RMS over the glaciated areas, up to 230 m**: much smaller over Fennoscandia and the Barents Sea (stiff cratonic mantle; Gulf of Bothnia peaks of 60-100 m against 300 m), larger over the Laurentide margins and West Antarctica. The global and ocean means are unchanged (bias -0.2 m). It moves the fit to PaleoMIST from 24 to 60 m RMS, as expected: PaleoMIST is built on the 1-D Earth. B3 - B1 is the science setting against Gowan's Earth, so it includes the 80 against 120 km lithosphere and the radial profile, not the lateral structure alone.
5. **The LGM is far from equilibrium with the solid Earth.** At 22.5 ka the equilibrium bed under the same ice is on average 120 m lower over the glaciated areas than the transient one, up to 500 m west of Hudson Bay (RSL 113 m higher, up to 470 m), with both Earths (1-D 147 m RMS, 3-D 166 m). The far field differs by about 12 m RMS. An LGM start from an equilibrium spin-up therefore overstates the LGM bed depression by hundreds of metres under the large ice sheets.
6. **The present day is not relaxed after this history** (up to 294 m of RSL left over Hudson Bay in 1-D, 429 m in 3-D). With `i_eq = 1` the present-day state of a deglaciation run is not the reference; comparisons with PaleoMIST must be relative to the 0 ka state, as Gowan's.

## Not as in Gowan et al. (2021)

- Steps as 1-yr ramps (see Setup); records compared 1 yr after their step.
- Earth: VILMA2's incompressible PREM shell means in 12 layers; the 670 km boundary, PREM and incompressibility assumed. VILMA2's explicit memory scheme (`fe`, cfl 1, sub-steps of 100 yr) against SELEN's normal modes.
- Equilibrium with the 80 ka load before the history (assumed); densities 931/1000 (assumed); VILMA2's degree-1 frame (centre of mass) and its flotation and sloping-coast ocean load.
- Inputs from FesmData: the 0.25 degree product remapped to the cell-centred FesmData grid (forcing), the 0.5 degree product (comparison), not Gowan's gridline-registered 0.25 degree file.
- The 3-D runs use the science setting (80 km lithosphere, Bagge 2021 on the CLIMBER-X radial profile), not Gowan's layers.
