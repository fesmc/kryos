---
title: "Simulations and scores"
toc: false
page-layout: full
---

A working table of the glacial-cycle simulations: what each one changes, and how well it fits the GAPSLIP 2.1 RSL index
points. Status as of 2026-11-13. How each simulation is made and scored, and
the concepts behind the ladder: [the workflow](plastic-workflow.md).

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
(`errmodel.toml`), from `plastic_ensemble.jl scores` (`experiments/scores.sh`):

- **log L**, the log-likelihood over all groups: higher is better. Differences of a few tens are large.
- **I** per data group, (χ² − n)/√(2n): ≤ 3 passes. The bar raises the cut to the best published score in a group
  where both publications fail: Greenland 7.9 (PaleoMIST), every other group 3. ESL keeps the plain cut. **passed**
  counts the groups within the bar.
- **RMSE**, the root mean square of model − obs [m] over all RSL points.
- **RMS z**, the root mean square of (model − obs)/σ_i, the classic score of the first comparisons (2026-11-09). A few
  points with tiny σ_i dominate it, Greenland above all.

Groups: Laurentide (892 points), Greenland (276), Fennoscandia (718), Barents (370), Antarctica (47), far field (2248),
and ESL, the eustatic sea level of the Spratt and Lisiecki (2016) stack from 25 to 10 ka against the run's ice-volume
sea level.

GNSS groups (`plastic_ensemble.jl`, `[[gnss]]`; columns to come): the present-day uplift rate at GNSS stations,
elastic-corrected, against the run's `dz_bed_rate` at its last VILMA2 record, with the same likelihood (sigma_m in mm/yr,
correlated over distance only): GNSS Fennoscandia (BIFROST2015) and GNSS Hudson Bay (Schumacher et al. 2018 and
MIDAS) judge; GNSS Greenland (GNET) and GNSS Antarctica (Martín-Español et al. 2016) are reported only. Their RMSE is
in mm/yr. A simulation without `dz_bed_rate` (P1, P2 and every run of VILMA2 before 2026-11-12) has no GNSS score. The
model rate is in VILMA2's centre-of-mass frame, the data in ITRF/IGS: a few tenths of a mm/yr apart, not corrected.

## The ladder

Each step changes one thing against its parent, so the change in score is that thing's. Steps:

- **P**: the published reconstructions.
- **I**: imposed ice on our solver, with the Earth varied.
- **T**: our ice online with PaleoMIST's τ_b.
- **D**: our ice with the present-day factors.
- **E**: ensembles of τ_b (and Earth).

Every row is made by one script of [experiments/](https://github.com/fesmc/kryos/tree/main/experiments), named in the
block's header row.

Columns (see [Scores](#scores)): **log L** the log-likelihood over all groups, higher is better; **passed** the groups
within the bar, of 7; **RMSE** over all RSL points [m]; **RMS z** of the normalized residuals; **I** the implausibility
(χ² − n)/√(2n) of one data group, ≤ 3 passes (Greenland ≤ 7.9, the bar), lower is better. The groups: **Laur.**
Laurentide, **Greenl.** Greenland, **Fenn.** Fennoscandia, **Bar.** the Barents Sea, **Antarc.** Antarctica, **FF**
far field, **ESL** the eustatic sea level.

![Every reference simulation of the ladder: log L (left) and I per data group against the bar (right; blue within it,
red beyond it, the colour the distance from the cut).](figures/simulations/ladder.png)

### Published reconstructions

The bar. Their own ice, Earth and GIA solver, and the RSL they published (FesmData `IceHistory-PaleoMIST` `dz_sl`,
`IceHistory-GLAC1D` topography); nothing is run. P2 has no ESL score: the published topography has no ice volume, so it
passes 4 of 6 groups.

| id | simulation | changes against | status | log L | passed | RMSE [m] | RMS z | I Laur. | I Greenl. | I Fenn. | I Bar. | I Antarc. | I FF | I ESL |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **P1** | `pub-PMIST` | published | done | −9311 | 3/7 | 20.0 | 30.1 | 6.7 | 7.9 | 0.4 | 2.6 | 3.4 | 6.4 | 7.0 |
| **P2** | `pub-GLAC1D` | published | done | −9117 | 4/6 | 19.0 | 18.6 | −0.3 | 9.8 | 1.8 | 5.6 | −0.2 | −0.6 | – |

### Imposed ice on our solver

A published ice history forced onto VILMA2 (`experiments/refs-ice.sh <id>`): the ice is fixed, so a change of score is
the solver's or the Earth's. I1 and I5 put each publication's ice on our solver with PaleoMIST's Earth; I2-I4 and I6-I8
change the Earth: Bagge's radial mean, then his lateral structure, then the strength of that structure.

| id | simulation | changes against | status | log L | passed | RMSE [m] | RMS z | I Laur. | I Greenl. | I Fenn. | I Bar. | I Antarc. | I FF | I ESL |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| I1 | `ice-PMIST \| E1D-PMIST` | P1: the solver (SELEN → VILMA2, 931 → 910 kg m⁻³, 2.5 kyr steps → linear, r) | done | −8860 | 6/7 | 16.2 | 28.5 | −3.5 | 6.8 | 1.0 | 0.7 | 0.4 | −0.6 | 5.9 |
| I2 | `ice-PMIST \| E1D-B21-6` | I1: the radial profile (PaleoMIST's → Bagge's mean, VILMA2's 1Db) | done | −8833 | 5/7 | **14.5** | 21.2 | −3.9 | 3.5 | 3.3 | 1.5 | 0.6 | −0.7 | 5.9 |
| | `ice-PMIST \| E1D-B21-6-layers` | the same radial mean as plain 1-D layer means | done | −8884 | 5/7 | **14.5** | 21.2 | −3.4 | 3.0 | 5.0 | 1.9 | 0.7 | −0.6 | 5.9 |
| I3 | `ice-PMIST \| E3D-B21-6` | I2: the lateral structure | done | −9243 | 3/7 | 19.6 | 35.3 | 3.1 | 5.3 | 15.2 | 3.5 | 0.9 | −1.0 | 6.0 |
| I4 | `ice-PMIST \| E3D-B21-1` | I3: Bagge's structure 1 (weakest lateral contrast) | done | −8902 | 5/7 | 15.0 | 20.1 | −3.3 | 2.7 | 5.7 | 1.8 | 0.5 | −0.6 | 5.9 |
| | `ice-PMIST \| E3D-B21-2` | structure 2 | done | −8913 | 5/7 | 15.1 | 20.1 | −3.0 | 2.7 | 6.0 | 1.9 | 0.5 | −0.6 | 5.9 |
| | `ice-PMIST \| E3D-B21-3` | structure 3, Bagge's reference | done | −8992 | 5/7 | 16.1 | 21.1 | −1.5 | 3.3 | 8.1 | 2.2 | 0.6 | −0.6 | 6.0 |
| | `ice-PMIST \| E3D-B21-4` | structure 4 | done | −9059 | 5/7 | 17.0 | 22.6 | −0.4 | 3.8 | 9.8 | 2.5 | 0.6 | −0.6 | 6.0 |
| | `ice-PMIST \| E3D-B21-5` | structure 5 | done | −9174 | 4/7 | 18.5 | 29.0 | 1.6 | 4.7 | 13.0 | 3.1 | 0.8 | −0.6 | 6.0 |
| | `ice-PMIST \| E3D-B21-7` | structure 7 (sc06 profile) | done | −9128 | 4/7 | 18.9 | 33.5 | 1.8 | 8.9 | 7.0 | 1.2 | 0.6 | −1.1 | 6.0 |
| | `ice-PMIST \| E3D-B21-8` | structure 8 (sc06b profile) | done | −9097 | 5/7 | 17.8 | 26.5 | 0.8 | 6.2 | 8.1 | 1.7 | 0.6 | −0.6 | 6.0 |
| I5 | `ice-GLAC1D \| E1D-PMIST` | I1: the ice (PaleoMIST → GLAC-1D) | done | −8835 | 6/7 | 17.8 | 28.2 | −3.0 | 5.7 | 0.5 | 1.2 | −0.2 | −0.8 | 3.3 |
| I6 | `ice-GLAC1D \| E1D-B21-6-layers` | I5: the radial profile (Bagge's mean as layers) | done | −8871 | 5/7 | 15.4 | **18.8** | −1.2 | 1.7 | 3.7 | 1.2 | −1.2 | −0.9 | 3.3 |
| I7 | `ice-GLAC1D \| E3D-B21-6` | I6: the lateral structure | done | −9226 | 3/7 | 20.4 | 26.1 | 4.0 | 6.4 | 14.6 | 3.0 | −0.4 | −0.7 | 3.4 |
| I8 | `ice-GLAC1D-9927 \| E1D-PMIST` | I5: North America nn9927 | done | **−8825** | 6/7 | 17.9 | 26.3 | −3.2 | 5.3 | 0.6 | 1.3 | −0.2 | −1.0 | 5.0 |
| | `ice-GLAC1D-9927 \| E3D-B21-6` | I7: North America nn9927 | done | −9196 | 4/7 | 20.5 | 27.1 | 4.6 | 6.1 | 14.7 | 3.0 | −0.5 | −2.0 | 5.0 |
| | `ice-PMIST \| E3D-PCA-E001` ... `E030` | I3: 30 Earths sampled from the expert basis (`earth-pca.sh`) | running | | | | | | | | | | | |

![PaleoMIST's ice on Bagge's structures 1-6, ordered by their lateral factor, and on the radial mean of structure 6 as a
1-D Earth (left). Dashed: the cut, I = 3.](figures/simulations/bagge.png)

Fennoscandia fails as soon as there is lateral structure, and the more the stronger it is (5.7 to 15.2); the radial mean
alone passes (3.3). Laurentide and Greenland degrade the same way, more slowly.

### Our ice, PaleoMIST's τ_b

PlasticIceSheet online on VILMA2, from PaleoMIST's basal shear stress (`experiments/refs-tau.sh <ens-v2> <id>`). T1
against I1 is the ice model: the same τ_b, but our ice is computed every step on the deforming bed. D adds the
present-day factors fitted to today's Greenland and Antarctica (Stage 1).

| id | simulation | changes against | status | log L | passed | RMSE [m] | RMS z | I Laur. | I Greenl. | I Fenn. | I Bar. | I Antarc. | I FF | I ESL |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T1 | `tau-PMIST \| E1D-PMIST` | I1: the ice model (ICESHEET 5 km offline → PlasticIceSheet 32 km online) | done | −9309 | 5/7 | 17.9 | 26.9 | −1.1 | 7.9 | −0.1 | 1.4 | 0.1 | 17.1 | 1.7 |
| T2 | `tau-PMIST \| E3D-B21-6` | T1: the 3-D Earth | proposed | | | | | | | | | | | |
| D1 | `tau-PMIST-PD \| E1D-PMIST` | T1: the present-day factors | done | −9030 | 5/7 | 17.8 | 27.7 | −0.9 | 8.8 | −0.1 | 1.1 | 0.7 | 5.0 | 2.3 |
| D2 | `tau-PMIST-PD \| E3D-B21-6` | D1: the 3-D Earth (m000 of ens-v2) | done | −9261 | 4/7 | 20.3 | 34.7 | 3.2 | 7.0 | 14.3 | 3.3 | 0.5 | 1.4 | 2.7 |
| D3 | `tau-PMIST-PD \| E3D-B21-6 @NH16` | D2: Northern Hemisphere resolution, 32 → 16 km (par `..._nh16.nml`) | done | −9275 | 3/7 | 20.6 | 35.4 | 3.1 | 7.8 | 14.4 | 3.4 | 0.5 | 1.4 | 3.1 |

The present-day factors fix most of the far-field misfit (T1 17.1 → D1 5.0). On Bagge's structure 6 (D2) Fennoscandia fails as it does for imposed ice (I3); 16 km in the north (D3) changes little.

### Ensembles

τ_b sampled per region and domain, and in ens-v2 the Earth too (shifts of structure 6); the members are weighted by their
likelihood (`experiments/ens.sh v1|v2`). The rows report the best members and the median.

| id | simulation | changes against | status | log L | passed | RMSE [m] | RMS z | I Laur. | I Greenl. | I Fenn. | I Bar. | I Antarc. | I FF | I ESL |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E1 | `tau-ENS1 \| E3D-B21-6` | D2: sampled τ_b (ens-v1, 25 members, Levante) | done, not rescored | | | | | | | | | | | |
| E2 | `tau-ENS2 \| E3D-B21-6-ENS2` | E1: the Earth sampled too (ens-v2, 50 members) | done | | | | | | | | | | | |
| | best log L: `tau-ENS2-E048 \| E3D-B21-6-E048` | | | −9017 | 4/7 | 19.7 | 32.2 | −0.5 | 13.0 | 4.5 | 1.1 | 0.9 | −0.6 | 3.0 |
| | most groups passed: `tau-ENS2-E036 \| E3D-B21-6-E036` | | | −9063 | 6/7 | 18.5 | 27.3 | 1.3 | 6.2 | 6.8 | 1.2 | 2.9 | 0.4 | 0.2 |
| | median of the 50 | | | −9474 | | | | | | | | | | |
| E3 | ensemble τ_b on `E1D-PMIST` | E2: our τ_b on PaleoMIST's Earth | proposed | | | | | | | | | | | |
| E4 | ensemble τ_b on a sampled 1-D Earth | E3: the 1-D Earth sampled | proposed | | | | | | | | | | | |
| E5 | ensemble τ_b on `E3D-B21-k` | E1: the Bagge structure as a parameter | proposed | | | | | | | | | | | |
| P1, P2 | FesmData `IceHistory-PaleoMIST` (`dz_sl`), `IceHistory-GLAC1D` (topography) |
| I1-I8, Earth-only | `kryos-work/kryos_data/refcorr/<ice>_<earth>`, e.g. `ice-GLAC1D_E1D-PMIST` |
| T1, D1 | `kryos-work/runs/pm-online/gowan`, `.../pdfac` |
| D2, E2 | `kryos-plastic/kryos/runs/ens-v2/iter_01/<m>` (m000; all members) |
| D3 | `kryos-work/runs/res16/m000-NH16` |
| E1 | Levante `/work/ba1442/robinson/kryos-ice-plastic/runs/ens-v1` |

The best members beat P2 on log L but not the imposed references (I1, I5), and Fennoscandia limits most of them: 4 of
the 50 pass it on structure 6 (best m013, I 0.2). E3-E5 move the ensemble onto other Earths; the Earth search below
decides which.

### Where the runs are

The scored runs, on albedo under `/albedo/work/user/alrobi001` (the scores: `kryos-work/refs-test/ens/iter_01/scores.csv`):

| id | run |
|---|---|
| P1, P2 | FesmData `IceHistory-PaleoMIST` (`dz_sl`), `IceHistory-GLAC1D` (topography) |
| I1-I8, Earth-only | `kryos-work/kryos_data/refcorr/<ice>_<earth>`, e.g. `ice-GLAC1D_E1D-PMIST` |
| T1, D1 | `kryos-work/runs/pm-online/gowan`, `.../pdfac` |
| D2, E2 | `kryos-plastic/kryos/runs/ens-v2/iter_01/<m>` (m000; all members) |
| D3 | `kryos-work/runs/res16/m000-NH16` |
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

How the search works, step by step, and what each piece is for: [the workflow](plastic-workflow.md#the-earth-search-step-by-step).

![Wave A: log L (RSL groups and ESL) of every 1-D Earth against each parameter, by ice history; dashed: I1 and
I5.](figures/earth-search/wave-a-logL.png)

The upper mantle orders the scores; a lithosphere thinner than about 150 km helps; the transition zone and the lower mantle
hardly matter in 1-D. On log L over the RSL groups, 33 Earths beat Gowan's under GLAC-1D's ice (I5) and 13 under
PaleoMIST's (I1).

![Wave A: the fraction of the 1-D Earth space the emulator does not rule out under either ice history, per parameter
pair.](figures/earth-search/wave-a-nroi-either.png)

**Wave B** (running, 2026-11-13): 50 3-D Earths, their radial part drawn from the space above and 8 lateral parameters
added (the experts' lateral pattern, its amplitude, regional upper-mantle shifts). Each Earth runs as its radial profile
alone (`E1D-ES2-E<nnn>`, done: 100 runs) and as the 3-D Earth from the profile's correction (`E3D-ES2-E<nnn>`, running).
The rows and the emulator follow when the runs are in.

| id | simulation | Earth | status |
|---|---|---|---|
| S3 | `ice-PMIST \| E1D-ES2-E<nnn>`, `ice-GLAC1D \| E1D-ES2-E<nnn>` | wave B radial profiles, 50 under each ice | done, not scored |
| S4 | `ice-PMIST \| E3D-ES2-E<nnn>`, `ice-GLAC1D \| E3D-ES2-E<nnn>` | wave B 3-D Earths, one pass from the profile's r | running |


## So far

- **Our solver beats both publications with their own ice.** I1 (PaleoMIST's ice) −8860 against −9311 (P1), and I5
  (GLAC-1D's ice) −8835 against −9117 (P2); both pass every RSL group. GLAC-1D's ice with nn9927 (I8, −8825) is the best
  simulation so far.
- **ESL separates the ice histories.** PaleoMIST's ice fails ESL on every Earth (I 5.9-6.0, about 7 m short at 20 ka),
  GLAC-1D's ice is close (3.3), and nn9927 is further off (5.0).
- **Bagge's radial profile is as good as Gowan's; his lateral structure is not.** The radial mean of structure 6 (I2)
  gives the lowest RMSE (14.5 m). With lateral structure, Fennoscandia fails, the more the larger the lateral contrast
  (structure 1: 5.7, structure 6: 15.2), on either ice (I3, I7). ens-v2 and the D rows sit on structure 6.
- **No online member reaches the imposed references yet.** The best ens-v2 members (m048 −9017, m036 −9063) beat GLAC-1D
  (P2) on log L but fail Fennoscandia (on structure 6) or Greenland. Our ice on Gowan's Earth (D1, −9030) passes
  Fennoscandia but fails the far field (5.0).
- **Shared misfit**, common to every good series: west Greenland (Sisimiut −100 m), the St. Lawrence (+16 to +56 m),
  and +10 to +20 m in the far field (Barbados, the Great Barrier Reef, the Bonaparte Gulf).
  `plastic_ensemble_misfit.jl` maps it.

## The figures

The ladder and Bagge figures are made from the score table by `tools/levante/plastic_sim_figures.jl` (on albedo or
Levante, then copied to `docs/figures/simulations/`):

```bash
julia --project=tools/levante tools/levante/plastic_sim_figures.jl $KRYOS_WORK/refs-test/ens/iter_01/scores.csv docs/figures/simulations
```

The figures of the Earth search (`docs/figures/earth-search/`) are those of `plastic_earth_emulate.jl`.

## Earlier runs, not on this scale

| run | in this language | note |
|---|---|---|
| B1, B3 (2026-10-31, `kryos_data/paleomist`) | `ice-PMIST \| E1D-PMIST`, `ice-PMIST \| E3D-B21-6` | VILMA2 523b660, 931 kg m⁻³, 2.5 kyr steps, no r; RMS z of B1 30.9 |
| plastic online test (2026-11-09) | about `tau-PMIST \| E1D-PMIST` | RMS z 25.9 against PaleoMIST 29.5 (RMS z of P1 on all points today: 29.5) |
| esmda-v1, -v2 (Levante) | calibrated τ_b on `E1D-PMIST` | ES-MDA calibrations, retired |
| ens-v1 (Levante) | `tau-ENS1 \| E3D-B21-6` | E1 above, scored under the ens-v1 error model |
