---
title: "Simulations and scores"
---

A working table of the glacial-cycle simulations: what each one changes, and how well it fits the GAPSLIP 2.1 RSL index
points. The page is not linked from the navigation. Status as of 2026-11-12.

The bar is not a perfect fit. It is the published reconstructions, PaleoMIST 1.0 (Gowan et al. 2021) and GLAC-1D
(Tarasov et al. 2012, Briggs et al. 2014): their misfits were accepted when they were published. A simulation is an
improvement when it fits better than they do, and the steps below are there to show which change made the difference.

## Naming

A simulation is written `ICE | EARTH`, with `@16` when the Northern Hemisphere ice runs on the 16 km grid (32 km by
default). Except for `pub`, the GIA solver is always VILMA2: toroidal flow on, ice density 910 kg m⁻³, and the present-day
correction r that belongs to the same ice and Earth (`gia_refcorr.sh`).

**Ice history**

| code | ice history |
|---|---|
| `PM` | published PaleoMIST 1.0 as a whole: its ice and its own RSL (SELEN, Gowan's 1-D Earth). Only as `PM \| pub` |
| `G1` | published GLAC-1D (PMIP4) as a whole: its ice and its own RSL (26-0 ka). Only as `G1 \| pub` |
| `PMi` | PaleoMIST ice **imposed** (FesmData `IceHistory-PaleoMIST`, minimal MIS 3 scenario) |
| `G1i` | GLAC-1D ice **imposed**, North America nn9894 (FesmData `IceHistory-GLAC1D_80k`); `G1i-9927` with nn9927 |
| `PI·τPM` | **PlasticIceSheet online** (the kryos plastic backend, every step, PaleoMIST's margins) with PaleoMIST's **published τ_b** (Gowan's domains and their values through time) |
| `PI·τPD` | PlasticIceSheet with τ_PM times the **present-day factors** (Stage 1: fitted to today's Greenland and Antarctica) |
| `PI·τE` | PlasticIceSheet with **sampled τ_b**, an ensemble member: `PI·τE036` is member m036 |
| `PI·τC` | PlasticIceSheet with **calibrated τ_b** (the ES-MDA calibrations of 2026-11-09, retired) |

In file and folder names: `PI-tPM`, `PI-tPD`, `PI-tE036`.

**Solid Earth**

| code | Earth |
|---|---|
| `pub` | the reconstruction's own GIA (SELEN for PaleoMIST, Tarasov's for GLAC-1D) |
| `E1G` | Gowan's 1-D Earth, PaleoMIST's: 120 km elastic lithosphere, upper mantle 4×10²⁰ Pa s, lower mantle 4×10²² Pa s |
| `E1B` | the radial mean of a Bagge structure (VILMA2's case 1Db: the 3-D setup with `visc_3d_tol = 99`) |
| `E3B` | Bagge et al. (2021) 3-D, structure 6 (v_1.0_s16); `E3B-k` for structure k of the eight in FesmData (`GLOBAL-T170_SolidEarth-Bagge2021`: reduction factor 0.2857, 0.3, 0.4, 0.5, 0.75, 1.0 with profile s16, 0.4 with sc06 and sc06b) |
| `E3Bx` | E3B with sampled log10 viscosity shifts (upper and lower mantle, five regions): a member's draw |
| `E1x` | a 1-D Earth with sampled upper- and lower-mantle viscosity and lithosphere (not built) |

## Scores

On the index points of age ≤ 25 ka (GLAC-1D's own RSL ends at 26 ka), under the error model refitted from ens-v2
(`runs/ens-v2/errmodel.toml`), from `plastic_ensemble.jl scores` (`runs/ens-v2/iter_01/scores.csv`):

- **log L**, the log-likelihood over all groups: higher is better. Differences of a few tens are large.
- **I** per data group, (χ² − n)/√(2n): ≤ 3 passes. The bar raises the cut to the best published score in a group:
  Greenland 7.9 (PaleoMIST), every other group 3. **passed** counts the groups within the bar.
- **RMSE**, the root mean square of model − obs [m] over all points.
- **RMS z**, the root mean square of (model − obs)/σ_i, the classic score of the first comparisons (2026-11-09). A few
  points with tiny σ_i dominate it, Greenland above all.

Groups: L Laurentide (892 points), G Greenland (276), F Fennoscandia/Barents (1088), A Antarctica (47), FF far field
(2248).

## The ladder

Each step changes one thing against its parent, so the change in score is that thing's.

| id | simulation | changes against | status | log L | passed | RMSE [m] | RMS z | I: L / G / F / A / FF |
|---|---|---|---|---|---|---|---|---|
| **P1** | `PM \| pub` | published | done | −9318 | 2/5 | 20.0 | 30.1 | 6.7 / 7.9 / 2.9 / 3.4 / 6.4 |
| **P2** | `G1 \| pub` | published | done | −9096 | 3/5 | 19.0 | 18.6 | −0.3 / 9.8 / 4.6 / −0.2 / −0.6 |
| **Imposed ice, our solver** | | | | | | | | |
| I1 | `PMi \| E1G` | P1: the solver (SELEN → VILMA2, 931 → 910 kg m⁻³, steps → linear, r) | done | **−8866** | **5/5** | **16.2** | 28.5 | −3.5 / 6.8 / 2.2 / 0.4 / −0.6 |
| I2 | `PMi \| E1B` | I1: the radial profile (Gowan → Bagge's mean) | running | | | | | |
| I3 | `PMi \| E3B` | I2: the lateral structure | running | | | | | |
| I4 | `PMi \| E3B-k` | I3: the other Bagge structures | proposed | | | | | |
| I5 | `G1i \| E1G` | I1: the ice (PaleoMIST → GLAC-1D) | waits for FesmData | | | | | |
| I6 | `G1i \| E1B` | I5: the radial profile | waits | | | | | |
| I7 | `G1i \| E3B` | I6: the lateral structure | waits | | | | | |
| I8 | `G1i-9927 \| E1G`, `\| E1B`, `\| E3B` | I5-I7: North America nn9927 | waits | | | | | |
| **Our ice, published τ_b** | | | | | | | | |
| O1 | `PI·τPM \| E1G` | I1: the ice model (ICESHEET 5 km offline → PlasticIceSheet 32 km online) | running | | | | | |
| O2 | `PI·τPM \| E3B` | O1: the 3-D Earth | proposed | | | | | |
| **Our ice, present-day factors** | | | | | | | | |
| D1 | `PI·τPD \| E1G` | O1: the present-day factors | running | | | | | |
| D2 | `PI·τPD \| E3B` | D1: the 3-D Earth (m000 of ens-v2) | done | −9243 | 3/5 | 20.3 | 34.7 | 3.2 / 7.0 / 13.7 / 0.5 / 1.4 |
| D3 | `PI·τPD \| E3B @16` | D2: Northern Hemisphere resolution | running | | | | | |
| **Our τ_b: ensembles** | | | | | | | | |
| E1 | `PI·τE \| E3B` | D2: sampled τ_b (ens-v1, 25 members, Levante) | done, not rescored | | | | | |
| E2 | `PI·τE \| E3Bx` | E1: the Earth sampled too (ens-v2, 50 members) | done | | | | | |
| | best log L: `PI·τE048 \| E3Bx` | | | −9017 | 3/5 | 19.7 | 32.2 | −0.5 / 13.0 / 5.1 / 0.9 / −0.6 |
| | closest to the bar: `PI·τE036 \| E3Bx` | | | −9063 | 4/5 | 18.5 | 27.3 | 1.3 / 6.2 / 7.1 / 2.9 / 0.4 |
| | median of the 50 | | | −9426 | | 23.0 | 37.3 | 5.9 / 7.8 / 12.2 / 1.0 / 3.0 |
| E3 | `PI·τE \| E1G` | E2: our τ_b on PaleoMIST's Earth | proposed | | | | | |
| E4 | `PI·τE \| E1x` | E3: the 1-D Earth sampled | proposed | | | | | |
| E5 | `PI·τE \| E3B-k` | E1: the Bagge structure as a parameter | proposed | | | | | |

The scored runs, on Albedo under `/albedo/work/user/alrobi001`:

| id | run |
|---|---|
| P1, P2 | FesmData `IceHistory-PaleoMIST` (`dz_sl`), `IceHistory-GLAC1D` (topography) |
| I1-I8 | `kryos-work/kryos_data/refcorr/<case>`: PaleoMIST-1D, PaleoMIST-Bagge1D, PaleoMIST-3D, GLAC1D-1D, ... |
| O1, D1 | `kryos-plastic/kryos/runs/pm-online/gowan`, `.../pdfac` |
| D2, E2 | `kryos-plastic/kryos/runs/ens-v2/iter_01/<m>` (m000; all members) |
| D3 | `kryos-plastic/kryos/runs/res16/m000-NH16` |
| E1 | Levante `/work/ba1442/robinson/kryos-ice-plastic/runs/ens-v1` |

## So far

- **Our solver with PaleoMIST's ice beats the publication** (I1 against P1): every group passes, RMSE 16 m against 20 m.
- **No ensemble member reaches I1 yet.** The best members match GLAC-1D (P2) on log L but fail Fennoscandia/Barents
  or Greenland. ens-v2 samples the Earth only around Bagge's 3-D structure. I2 and I3 show whether that Earth is what
  holds the members back, and O1, D1 and E3 show what our ice does on PaleoMIST's Earth.
- **Shared misfit**, common to every good series: west Greenland (Sisimiut −100 m), the St. Lawrence (+16 to +56 m),
  and +10 to +20 m in the far field (Barbados, the Great Barrier Reef, the Bonaparte Gulf).
  `plastic_ensemble_misfit.jl` maps it.

## Earlier runs, not on this scale

| run | in this language | note |
|---|---|---|
| B1, B3 (2026-10-31, `kryos_data/paleomist`) | `PMi \| E1G`, `PMi \| E3B` | VILMA2 523b660, 931 kg m⁻³, 2.5 kyr steps, no r; RMS z of B1 30.9 |
| plastic online test (2026-11-09) | about `PI·τPM \| E1G` | RMS z 25.9 against PaleoMIST 29.5 (RMS z of P1 on all points today: 29.5) |
| esmda-v1, -v2 (Levante) | `PI·τC \| E1G` | ES-MDA calibrations, retired |
| ens-v1 (Levante) | `PI·τE \| E3B` | E1 above, scored under the ens-v1 error model |
