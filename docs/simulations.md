---
title: "Simulations and scores"
toc: false
page-layout: full
---

A working table of the glacial-cycle simulations: what each one changes, and how well it fits the GAPSLIP 2.1 RSL index
points. The page is not linked from the navigation. Status as of 2026-11-12.

The bar is not a perfect fit. It is the published reconstructions, PaleoMIST 1.0 (Gowan et al. 2021) and GLAC-1D
(Tarasov et al. 2012, Briggs et al. 2014): their misfits were accepted when they were published. A simulation is an
improvement when it fits better than they do, and the ladder below is there to show which change made the difference.

## Naming

Every component is named `<thing>-<source>-<variant>`: what it is, where it comes from, and which version of it, the
variant left out when there is only one. A simulation is its ice and its Earth, written `<ice> | <earth>`.

**Things**

| thing | what it is |
|---|---|
| `pub` | a published reconstruction as a whole: its ice, its Earth and its own GIA solver, and the RSL it published |
| `ice` | an ice history **imposed** on VILMA2 (the ice does not respond to the Earth) |
| `tau` | basal shear stress for **PlasticIceSheet run online** (the kryos plastic backend: the ice is computed every step from τ_b and the margins, on the bed that VILMA2 deforms) |
| `E1D` | a 1-D (radially symmetric) Earth |
| `E3D` | a 3-D Earth (lateral viscosity) |

**Sources**

| source | what it is |
|---|---|
| `PMIST` | PaleoMIST 1.0 (Gowan et al. 2021): ice of the minimal MIS 3 scenario (FesmData `IceHistory-PaleoMIST`), τ_b of its domains through time (`IceHistory-PaleoMIST_taub`), its 1-D Earth (120 km elastic lithosphere, upper mantle 4×10²⁰ Pa s, lower mantle 4×10²² Pa s) |
| `GLAC1D` | GLAC-1D (Tarasov et al. 2012, Briggs et al. 2014): from 80 ka, North America nn9894 (FesmData `IceHistory-GLAC1D_80k`) |
| `B21` | Bagge et al. (2021): the 3-D viscosity structures of FesmData `SolidEarth-Bagge2021` |
| `L24`, `P22` | Lloyd et al. (2024), Pan et al. (2022): other 3-D viscosity fields in FesmData (not used yet) |
| `ENS1`, `ENS2` | our weighted ensembles (ens-v1, ens-v2) |
| `ES1`, `ES2` | the broad Earth search (docs/dev/earth-search-design.md): wave A, 200 1-D Earths of 5 radial parameters; wave B, 50 3-D Earths of 13 parameters (`E1D-ES2-E<nnn>` the radial profile of `E3D-ES2-E<nnn>`) |

**Variants**

| variant | of | meaning |
|---|---|---|
| `9927` | `ice-GLAC1D` | North America nn9927 instead of nn9894 |
| `PD` | `tau-PMIST` | times the present-day factors (Stage 1, fitted to today's Greenland and Antarctica) |
| `E036` | `tau-ENS2`, `E3D-B21-6` | ensemble member m036: its sampled τ_b, or its sampled viscosity shifts |
| `E<nnn>`, `PMIST` | `E1D-ES1`, `E1D-ES2`, `E3D-ES2` | Earth E<nnn> of the search; `E1D-ES1-PMIST` is Gowan's Earth in the search's parameters (the check) |
| `1` ... `8` | `E3D-B21`, `E1D-B21` | Bagge structure k (below); `E1D-B21-k` is its radial mean (VILMA2's case 1Db, `visc_3d_tol = 99`) |

The Bagge structures (`structure` in FesmData; lateral reduction factor and radial viscosity profile):

| k | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| Bagge name | v_0.2857_s16 | v_0.3_s16 | v_0.4_s16 | v_0.5_s16 | v_0.75_s16 | v_1.0_s16 | v_0.4_sc06 | v_0.4_sc06b |

Bagge's reference structure is 3 (v_0.4_s16). Ours so far is 6 (v_1.0_s16), the one with the largest lateral contrast.

**Defaults.** The GIA solver is VILMA2 except for `pub`: toroidal flow on, ice density 910 kg m⁻³, degree 170, and the
present-day correction r that belongs to the same ice and Earth (`gia_refcorr.sh`). PlasticIceSheet runs the Northern
Hemisphere on 32 km, Antarctica on 16 km and Patagonia on 8 km. A suffix `@NH16` marks the Northern Hemisphere on 16 km.

Examples: `ice-PMIST | E1D-PMIST` is PaleoMIST's ice and Earth with our solver; `tau-PMIST-PD | E3D-B21-6` is our ice with
PaleoMIST's τ_b and the present-day factors on Bagge's structure 6 (m000 of ens-v2); `tau-ENS2-E036 | E3D-B21-6-E036` is
member m036 of ens-v2.

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

GNSS groups (`plastic_ensemble.jl`, `[[gnss]]`; columns to come): the present-day uplift rate at GNSS stations,
elastic-corrected, against the run's `dz_bed_rate` at its last VILMA2 record, with the same likelihood (sigma_m in mm/yr,
correlated over distance only): GNSS Fennoscandia (BIFROST2015) and GNSS Hudson Bay (Schumacher et al. 2018 and
MIDAS) judge; GNSS Greenland (GNET) and GNSS Antarctica (Martín-Español et al. 2016) are reported only. Their RMSE is
in mm/yr. A simulation without
`dz_bed_rate` (P1, P2 and every run of VILMA2 before 2026-11-12) has no GNSS score. The model rate is in VILMA2's
centre-of-mass frame, the data in ITRF/IGS: a few tenths of a mm/yr apart, not corrected.

## The ladder

Each step changes one thing against its parent, so the change in score is that thing's. Steps:

- **P**: the published reconstructions.
- **I**: imposed ice on our solver, with the Earth varied.
- **T**: our ice online with PaleoMIST's τ_b.
- **D**: our ice with the present-day factors.
- **E**: ensembles of τ_b (and Earth).

Columns (see [Scores](#scores)): **log L** the log-likelihood over all groups, higher is better; **passed** the groups
within the bar; **RMSE** over all points [m]; **RMS z** of the normalized residuals; **I** the implausibility
(χ² − n)/√(2n) of one data group, ≤ 3 passes (Greenland ≤ 7.9, the bar), lower is better. The groups: **Laur.** Laurentide
(892 points), **Greenl.** Greenland (276), **Fenn./Bar.** Fennoscandia and the Barents Sea (1088), **Antarc.** Antarctica
(47), **far field** (2248).

| id | simulation | changes against | status | log L | passed | RMSE [m] | RMS z | I Laur. | I Greenl. | I Fenn./Bar. | I Antarc. | I far field |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **P1** | `pub-PMIST` | published | done | −9318 | 2/5 | 20.0 | 30.1 | 6.7 | 7.9 | 2.9 | 3.4 | 6.4 |
| **P2** | `pub-GLAC1D` | published | done | −9096 | 3/5 | 19.0 | 18.6 | −0.3 | 9.8 | 4.6 | −0.2 | −0.6 |
| **Imposed ice, our solver** | | | | | | | |  |  |  |  |  |
| I1 | `ice-PMIST \| E1D-PMIST` | P1: the solver (SELEN → VILMA2, 931 → 910 kg m⁻³, 2.5 kyr steps → linear, r) | done | **−8866** | **5/5** | **16.2** | 28.5 | −3.5 | 6.8 | 2.2 | 0.4 | −0.6 |
| I2 | `ice-PMIST \| E1D-B21-6` | I1: the radial profile (PaleoMIST's → Bagge's mean) | running | | | | |  |  |  |  |  |
| I3 | `ice-PMIST \| E3D-B21-6` | I2: the lateral structure | running | | | | |  |  |  |  |  |
| I4 | `ice-PMIST \| E3D-B21-k`, k ≠ 6 | I3: Bagge's other structures, 3 his reference | proposed | | | | |  |  |  |  |  |
| I5 | `ice-GLAC1D \| E1D-PMIST` | I1: the ice (PaleoMIST → GLAC-1D) | waits for FesmData | | | | |  |  |  |  |  |
| I6 | `ice-GLAC1D \| E1D-B21-6` | I5: the radial profile | waits | | | | |  |  |  |  |  |
| I7 | `ice-GLAC1D \| E3D-B21-6` | I6: the lateral structure | waits | | | | |  |  |  |  |  |
| I8 | `ice-GLAC1D-9927 \| E1D-PMIST`, `\| E1D-B21-6`, `\| E3D-B21-6` | I5-I7: North America nn9927 | waits | | | | |  |  |  |  |  |
| **Our ice, PaleoMIST's τ_b** | | | | | | | |  |  |  |  |  |
| T1 | `tau-PMIST \| E1D-PMIST` | I1: the ice model (ICESHEET 5 km offline → PlasticIceSheet 32 km online) | running | | | | |  |  |  |  |  |
| T2 | `tau-PMIST \| E3D-B21-6` | T1: the 3-D Earth | proposed | | | | |  |  |  |  |  |
| **Our ice, present-day factors** | | | | | | | |  |  |  |  |  |
| D1 | `tau-PMIST-PD \| E1D-PMIST` | T1: the present-day factors | running | | | | |  |  |  |  |  |
| D2 | `tau-PMIST-PD \| E3D-B21-6` | D1: the 3-D Earth (m000 of ens-v2) | done | −9243 | 3/5 | 20.3 | 34.7 | 3.2 | 7.0 | 13.7 | 0.5 | 1.4 |
| D3 | `tau-PMIST-PD \| E3D-B21-6 @NH16` | D2: Northern Hemisphere resolution, 32 → 16 km (ESL at 20 ka −128.4 against −130.4 m; par `..._nh16.nml`) | done | −9256 | 3/5 | 20.6 | 35.4 | 3.1 | 7.8 | 13.8 | 0.5 | 1.5 |
| **Ensembles** | | | | | | | |  |  |  |  |  |
| E1 | `tau-ENS1 \| E3D-B21-6` | D2: sampled τ_b (ens-v1, 25 members, Levante) | done, not rescored | | | | |  |  |  |  |  |
| E2 | `tau-ENS2 \| E3D-B21-6-ENS2` | E1: the Earth sampled too (ens-v2, 50 members) | done | | | | |  |  |  |  |  |
| | best log L: `tau-ENS2-E048 \| E3D-B21-6-E048` | | | −9017 | 3/5 | 19.7 | 32.2 | −0.5 | 13.0 | 5.1 | 0.9 | −0.6 |
| | closest to the bar: `tau-ENS2-E036 \| E3D-B21-6-E036` | | | −9063 | 4/5 | 18.5 | 27.3 | 1.3 | 6.2 | 7.1 | 2.9 | 0.4 |
| | median of the 50 | | | −9426 | | 23.0 | 37.3 | 5.9 | 7.8 | 12.2 | 1.0 | 3.0 |
| E3 | ensemble τ_b on `E1D-PMIST` | E2: our τ_b on PaleoMIST's Earth | proposed | | | | |  |  |  |  |  |
| E4 | ensemble τ_b on a sampled 1-D Earth | E3: the 1-D Earth sampled | proposed | | | | |  |  |  |  |  |
| E5 | ensemble τ_b on `E3D-B21-k` | E1: the Bagge structure as a parameter | proposed | | | | |  |  |  |  |  |

The scored runs, on Albedo under `/albedo/work/user/alrobi001`:

| id | run |
|---|---|
| P1, P2 | FesmData `IceHistory-PaleoMIST` (`dz_sl`), `IceHistory-GLAC1D` (topography) |
| I1-I8 | `kryos-work/kryos_data/refcorr/<case>`: PaleoMIST-1D (I1), PaleoMIST-Bagge1D (I2), PaleoMIST-3D (I3), GLAC1D-1D, ... |
| T1, D1 | `kryos-plastic/kryos/runs/pm-online/gowan`, `.../pdfac` |
| D2, E2 | `kryos-plastic/kryos/runs/ens-v2/iter_01/<m>` (m000; all members) |
| D3 | `kryos-plastic/kryos/runs/res16/m000-NH16` |
| E1 | Levante `/work/ba1442/robinson/kryos-ice-plastic/runs/ens-v1` |

## The Earth search

Imposed ice on sampled Earths (docs/dev/earth-search-design.md), scored as the rows above but on the split groups, the GNSS
groups and ESL: the error model of refs-test (Fennoscandia, Barents refitted after the split), the GNSS groups refitted on the
wave-A runs (Fennoscandia 0.7 mm/yr and 750 km, Hudson Bay 1.5 mm/yr and 500 km), ESL at its starting values. log L over the RSL
groups (on this scale I1 −8860, I5 −8835); passed of 9 groups (6 RSL, GNSS Fennoscandia and Hudson Bay, ESL; Greenland ≤ 7.9 and
ESL ≤ 7.0 by the bar). Earth: lithosphere [km], log10 η UM / TZ / LM1 / LM2.

| id | simulation | Earth | log L | passed | I L / G / F / B / A / FF | I GNSS F / HB | I ESL |
|---|---|---|---|---|---|---|---|
| S0 | `ice-PMIST \| E1D-ES1-PMIST` | Gowan's, in the parameters (the check: I1 again) | −8860.3 | 7/9 | −3.5 / 6.8 / 1.0 / 0.7 / 0.4 / −0.6 | 7.4 / 3.5 | 5.9 |
| S1 | `ice-PMIST \| E1D-ES1-E<nnn>` | wave A, 200 Earths: not ruled out 20.2 % of the space (30.3 % without GNSS) | | | | | |
| | best log L: `ice-PMIST \| E1D-ES1-E117` | 111, 20.13 / 21.75 / 23.15 / 21.69 | −8800.7 | 7/9 | −3.2 / 2.0 / 1.6 / 1.4 / 0.2 / −0.7 | 5.7 / 6.4 | 5.9 |
| | most groups: `ice-PMIST \| E1D-ES1-E017` | 87, 20.56 / 20.61 / 21.68 / 21.02 | −8905.5 | 8/9 | −0.7 / 2.5 / 0.6 / 1.0 / 0.8 / −1.6 | 4.9 / 2.1 | 6.0 |
| S2 | `ice-GLAC1D \| E1D-ES1-E<nnn>` | wave A: not ruled out 13.5 % (26.3 % without GNSS); under either ice 28.8 %, both 4.9 % | | | | | |
| | best log L: `ice-GLAC1D \| E1D-ES1-E165` | 67, 20.77 / 20.08 / 22.54 / 21.88 | −8708.9 | 7/9 | −3.5 / 2.6 / 0.8 / 3.3 / −1.4 / −4.4 | 3.6 / 2.9 | 3.3 |
| | all groups: `ice-GLAC1D \| E1D-ES1-E033` (14 Earths pass all 9) | 106, 20.84 / 21.43 / 21.08 / 22.08 | −8745.6 | 9/9 | −2.4 / 1.7 / 2.2 / 0.7 / −1.1 / −4.4 | 0.5 / 1.1 | 3.4 |

## So far

- **Our solver with PaleoMIST's ice beats the publication** (I1 against P1): every group passes, RMSE 16 m against 20 m.
- **No ensemble member reaches I1 yet.** The best members match GLAC-1D (P2) on log L but fail Fennoscandia/Barents
  or Greenland. ens-v2 samples the Earth only around Bagge's structure 6. I2 and I3 show whether that Earth is what
  holds the members back, and T1, D1 and E3 show what our ice does on PaleoMIST's Earth.
- **Shared misfit**, common to every good series: west Greenland (Sisimiut −100 m), the St. Lawrence (+16 to +56 m),
  and +10 to +20 m in the far field (Barbados, the Great Barrier Reef, the Bonaparte Gulf).
  `plastic_ensemble_misfit.jl` maps it.

## Earlier runs, not on this scale

| run | in this language | note |
|---|---|---|
| B1, B3 (2026-10-31, `kryos_data/paleomist`) | `ice-PMIST \| E1D-PMIST`, `ice-PMIST \| E3D-B21-6` | VILMA2 523b660, 931 kg m⁻³, 2.5 kyr steps, no r; RMS z of B1 30.9 |
| plastic online test (2026-11-09) | about `tau-PMIST \| E1D-PMIST` | RMS z 25.9 against PaleoMIST 29.5 (RMS z of P1 on all points today: 29.5) |
| esmda-v1, -v2 (Levante) | calibrated τ_b on `E1D-PMIST` | ES-MDA calibrations, retired |
| ens-v1 (Levante) | `tau-ENS1 \| E3D-B21-6` | E1 above, scored under the ens-v1 error model |
