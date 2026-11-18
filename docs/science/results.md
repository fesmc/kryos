---
title: "Results"
toc: false
page-layout: full
---

Every scored simulation, on one scale: what each changes against its parent ([Methods](methods.md#the-ladder)) and how
well it fits ([Data and scores](data-scores.md)). Status 2026-11-18 (ESL rescored: the barystatic sea level, linear in time;
[Data and scores](data-scores.md#esl)).

**Columns.** **log L**, the log-likelihood over the six RSL groups, higher is better (differences of a few tens are
large); **passed**, the groups within the bar, of 7 (RSL and ESL) or of 9 with the GNSS rates (of 8 in the Earth search,
which reports ESL but does not judge it); **RMSE** over all RSL
points [m]; **RMS z** of the normalized residuals; **I** per group, ≤ 3 passes (Greenland ≤ 7.9; ESL ≤ 3 for all:
the Spratt and Lisiecki stack is its target; PaleoMIST no longer sets the ESL bar): **Laur.** Laurentide, **Greenl.**
Greenland, **Fenn.** Fennoscandia, **Bar.** Barents,
**Antarc.** Antarctica, **FF** far field, **ESL** the barystatic sea level against the Spratt and Lisiecki (2016) stack.

![Every reference simulation of the ladder: log L (left) and I per data group against the bar (right; blue within it,
red beyond it, the colour the distance from the cut).](figures/simulations/ladder.png)

## Published reconstructions

The bar: their own ice, Earth and GIA solver, and the RSL they published; nothing is run. P2 has no ESL score (the
published topography has no ice volume).

| id | simulation | changes against | status | log L | passed | RMSE [m] | RMS z | I Laur. | I Greenl. | I Fenn. | I Bar. | I Antarc. | I FF | I ESL |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **P1** | `pub-PMIST` | published | done | −9311 | 3/7 | 20.0 | 30.1 | 6.7 | 7.9 | 0.4 | 2.6 | 3.4 | 6.4 | 7.2 |
| **P2** | `pub-GLAC1D` | published | done | −9117 | 4/6 | 19.0 | 18.6 | −0.3 | 9.8 | 1.8 | 5.6 | −0.2 | −0.6 | – |

## Imposed ice on our solver

A published ice history forced onto VILMA2: the ice is fixed, so a change of score is the solver's or the Earth's. I1
and I5 put each publication's ice on our solver with PaleoMIST's Earth; I2-I4 and I6-I8 change the Earth: Bagge's
radial mean, then his lateral structure, then the strength of that structure.

| id | simulation | changes against | status | log L | passed | RMSE [m] | RMS z | I Laur. | I Greenl. | I Fenn. | I Bar. | I Antarc. | I FF | I ESL |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| I1 | `ice-PMIST \| E1D-PMIST` | P1: the solver (SELEN → VILMA2, 931 → 910 kg m⁻³, 2.5 kyr steps → linear, r) | done | −8860 | 6/7 | 16.2 | 28.5 | −3.5 | 6.8 | 1.0 | 0.7 | 0.4 | −0.6 | 10.2 |
| I2 | `ice-PMIST \| E1D-B21-6` | I1: the radial profile (PaleoMIST's → Bagge's mean) | done | −8833 | 5/7 | **14.5** | 21.2 | −3.9 | 3.5 | 3.3 | 1.5 | 0.6 | −0.7 | 10.1 |
| | `ice-PMIST \| E1D-B21-6-layers` | the same radial mean as plain 1-D layer means | done | −8884 | 5/7 | **14.5** | 21.2 | −3.4 | 3.0 | 5.0 | 1.9 | 0.7 | −0.6 | 10.1 |
| I3 | `ice-PMIST \| E3D-B21-6` | I2: the lateral structure | done | −9243 | 3/7 | 19.6 | 35.3 | 3.1 | 5.3 | 15.2 | 3.5 | 0.9 | −1.0 | 10.3 |
| I4 | `ice-PMIST \| E3D-B21-1` | I3: Bagge's structure 1 (weakest lateral contrast) | done | −8902 | 5/7 | 15.0 | 20.1 | −3.3 | 2.7 | 5.7 | 1.8 | 0.5 | −0.6 | 10.2 |
| | `ice-PMIST \| E3D-B21-2` | structure 2 | done | −8913 | 5/7 | 15.1 | 20.1 | −3.0 | 2.7 | 6.0 | 1.9 | 0.5 | −0.6 | 10.2 |
| | `ice-PMIST \| E3D-B21-3` | structure 3, Bagge's reference | done | −8992 | 5/7 | 16.1 | 21.1 | −1.5 | 3.3 | 8.1 | 2.2 | 0.6 | −0.6 | 10.2 |
| | `ice-PMIST \| E3D-B21-4` | structure 4 | done | −9059 | 5/7 | 17.0 | 22.6 | −0.4 | 3.8 | 9.8 | 2.5 | 0.6 | −0.6 | 10.2 |
| | `ice-PMIST \| E3D-B21-5` | structure 5 | done | −9174 | 4/7 | 18.5 | 29.0 | 1.6 | 4.7 | 13.0 | 3.1 | 0.8 | −0.6 | 10.2 |
| | `ice-PMIST \| E3D-B21-7` | structure 7 (sc06 profile) | done | −9128 | 4/7 | 18.9 | 33.5 | 1.8 | 8.9 | 7.0 | 1.2 | 0.6 | −1.1 | 10.4 |
| | `ice-PMIST \| E3D-B21-8` | structure 8 (sc06b profile) | done | −9097 | 5/7 | 17.8 | 26.5 | 0.8 | 6.2 | 8.1 | 1.7 | 0.6 | −0.6 | 10.3 |
| I5 | `ice-GLAC1D \| E1D-PMIST` | I1: the ice (PaleoMIST → GLAC-1D) | done | −8835 | 6/7 | 17.8 | 28.2 | −3.0 | 5.7 | 0.5 | 1.2 | −0.2 | −0.8 | 7.6 |
| I6 | `ice-GLAC1D \| E1D-B21-6-layers` | I5: the radial profile (Bagge's mean as layers) | done | −8871 | 5/7 | 15.4 | **18.8** | −1.2 | 1.7 | 3.7 | 1.2 | −1.2 | −0.9 | 7.5 |
| I7 | `ice-GLAC1D \| E3D-B21-6` | I6: the lateral structure | done | −9226 | 3/7 | 20.4 | 26.1 | 4.0 | 6.4 | 14.6 | 3.0 | −0.4 | −0.7 | 7.5 |
| I8 | `ice-GLAC1D-9927 \| E1D-PMIST` | I5: North America nn9927 | done | **−8825** | 6/7 | 17.9 | 26.3 | −3.2 | 5.3 | 0.6 | 1.3 | −0.2 | −1.0 | 9.8 |
| | `ice-GLAC1D-9927 \| E3D-B21-6` | I7: North America nn9927 | done | −9196 | 4/7 | 20.5 | 27.1 | 4.6 | 6.1 | 14.7 | 3.0 | −0.5 | −2.0 | 9.7 |
| | `ice-PMIST \| E3D-PCA-E001` ... `E030` | I3: 30 Earths sampled from the expert basis | run, not scored | | | | | | | | | | | |

![PaleoMIST's ice on Bagge's structures 1-6, ordered by their lateral factor, and on the radial mean of structure 6 as a
1-D Earth (left). Dashed: the cut, I = 3.](figures/simulations/bagge.png)

- **Our solver beats both publications with their own ice:** I1 −8860 against P1 −9311, I5 −8835 against P2 −9117,
  both passing every RSL group. GLAC-1D's ice with nn9927 (I8) is the best of the ladder.
- **Both ice histories fail ESL.** PaleoMIST's ice on every Earth (I 10.1-10.4, about 16 m short at 20 ka: 110 against
  126 m), GLAC-1D's ice less (7.4-7.6), nn9927 further off (9.7-9.8).
- **Bagge's radial profile is as good as Gowan's, his lateral structure is not.** The radial mean of structure 6 (I2)
  gives the lowest RMSE. With lateral structure Fennoscandia fails, the more the larger the lateral contrast (structure
  1: 5.7, structure 6: 15.2), under either ice (I3, I7). Laurentide and Greenland degrade the same way, more slowly.

## Our ice, PaleoMIST's τ_b

PlasticIceSheet online on VILMA2 from PaleoMIST's τ_b. T1 against I1 is the ice model: the same τ_b, but our ice is
computed every step on the deforming bed. D adds the present-day factors.

| id | simulation | changes against | status | log L | passed | RMSE [m] | RMS z | I Laur. | I Greenl. | I Fenn. | I Bar. | I Antarc. | I FF | I ESL |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T1 | `tau-PMIST \| E1D-PMIST` | I1: the ice model (ICESHEET 5 km offline → PlasticIceSheet 32 km online) | done | −9309 | 5/7 | 17.9 | 26.9 | −1.1 | 7.9 | −0.1 | 1.4 | 0.1 | 17.1 | 3.0 |
| D1 | `tau-PMIST-PD \| E1D-PMIST` | T1: the present-day factors | done | −9030 | 4/7 | 17.8 | 27.7 | −0.9 | 8.8 | −0.1 | 1.1 | 0.7 | 5.0 | 5.0 |
| D2 | `tau-PMIST-PD \| E3D-B21-6` | D1: the 3-D Earth (m000 of ens-v2) | done | −9261 | 3/7 | 20.3 | 34.7 | 3.2 | 7.0 | 14.3 | 3.3 | 0.5 | 1.4 | 5.9 |
| D3 | `tau-PMIST-PD \| E3D-B21-6 @NH16` | D2: Northern Hemisphere 32 → 16 km | done | −9275 | 3/7 | 20.6 | 35.4 | 3.1 | 7.8 | 14.4 | 3.4 | 0.5 | 1.4 | 6.6 |

The present-day factors remove most of the far-field misfit (T1 17.1 → D1 5.0). On structure 6 (D2) Fennoscandia
fails as it does under imposed ice (I3). 16 km in the north (D3) changes the scores little and costs 2 m of ESL; 32 km
stays the default.

## Our ice, the Coulomb law

PlasticIceSheet online on VILMA2 with the Coulomb yield law instead of PaleoMIST's τ_b (`taub_method = "coulomb"`, λ on
the relaxed bed, one c_f per ice sheet): four (H_N, λ_min) of the band that meets the 20 ka targets (H_N 30, 100, 300,
1000 m with λ_min 0.2, 0.3, 0.5, 0.75; c_f of the calibration, Patagonia Eurasia's), each on Gowan's Earth and on the
Earth search's E117, with the r of PaleoMIST's ice on that Earth. Scored on the Earth search's scale, but with ESL judged
(9 groups; the bars Greenland 7.9, ESL 3 against the Spratt stack: every Coulomb run but C5 at H_N 1000 fails ESL).

| id | simulation | changes against | status | log L | passed | RMSE [m] | I L / G / F / B / A / FF | I GNSS F / HB | I ESL |
|---|---|---|---|---|---|---|---|---|---|
| C1 | `tau-CB-H030 \| E1D-PMIST` | D1: the Coulomb law, H_N 30 m | done | −9016 | 5/9 | 22.3 | −2.3 / 4.9 / 2.4 / 6.8 / 2.8 / −0.9 | 7.4 / 11.1 | 10.3 |
| | `tau-CB-H100 \| E1D-PMIST` | H_N 100 m | done | −8977 | 5/9 | 21.5 | −2.5 / 4.5 / 1.8 / 6.0 / 2.8 / −1.4 | 7.0 / 10.9 | 11.1 |
| | `tau-CB-H300 \| E1D-PMIST` | H_N 300 m | done | −8918 | 4/9 | 20.1 | −3.5 / 3.9 / 1.1 / 5.0 / 3.4 / −0.9 | 6.2 / 9.8 | 11.0 |
| | `tau-CB-H1000 \| E1D-PMIST` | H_N 1000 m | done | −8885 | 4/9 | 18.5 | −2.5 / 3.7 / −0.2 / 4.5 / 4.0 / −0.8 | 6.6 / 11.8 | 9.8 |
| C2 | `tau-CB-H030 \| E1D-ES1-E117` | C1: the Earth E117 | done | −8925 | 5/9 | 18.0 | −0.9 / 1.5 / 1.0 / 5.8 / 2.6 / −0.8 | 8.1 / 6.3 | 10.5 |
| | `tau-CB-H100 \| E1D-ES1-E117` | | done | −8910 | 5/9 | 17.5 | −1.1 / 1.4 / 0.5 / 5.0 / 2.6 / −1.4 | 7.4 / 6.4 | 11.2 |
| | `tau-CB-H300 \| E1D-ES1-E117` | | done | −8840 | 5/9 | 16.7 | −1.9 / 1.2 / −0.3 / 4.1 / 3.0 / −1.7 | 6.9 / 6.6 | 9.1 |
| | `tau-CB-H1000 \| E1D-ES1-E117` | | done | −8829 | 4/9 | 15.8 | −1.1 / 0.9 / −1.0 / 3.9 / 4.0 / −1.1 | 7.0 / 8.6 | 10.0 |
| C3 | `tau-CBS-H1000 \| E1D-ES1-E117` | C2 H1000: Eurasia split into Fennoscandia and BarentsKara | done | −8811 | 5/9 | 14.0 | −1.2 / 0.9 / 0.6 / 2.2 / 4.0 / −2.1 | 6.4 / 8.6 | 10.5 |
| | `tau-CBS-H2000 \| E1D-ES1-E117` | H_N 2000 m, λ_min 1.0 | done | −8850 | 5/9 | 14.0 | −0.2 / 0.8 / 1.0 / 2.3 / 4.7 / −0.9 | 6.8 / 11.0 | 9.9 |
| | `tau-CBS-H3000 \| E1D-ES1-E117` | H_N 3000 m, λ_min 1.0 | done | −8871 | 5/9 | 14.3 | 0.5 / 0.8 / 1.2 / 2.3 / 5.0 / −1.3 | 7.0 / 12.8 | 10.9 |
| C4 | `tau-CBS-H1000-pPG \| E1D-ES1-E117` | C3 H1000: Greenland's margin from PaleoGrIS 1.0, 15-5 ka (a margin patch) | done | −8811 | 5/9 | 13.9 | −1.1 / 0.6 / 0.6 / 2.2 / 4.1 / −2.0 | 6.3 / 8.6 | 10.5 |
| C5 | `tau-CBD-H1000 \| E1D-ES1-E117` | C3 H1000: the floor of N, δ = 0.05, c_f recalibrated (lgm-calib-v3) | done | −9194 | 5/9 | 14.2 | −0.9 / 1.2 / −0.1 / 2.5 / 4.7 / 16.6 | 5.8 / 9.7 | 2.1 |
| | `tau-CBD-H300 \| E1D-ES1-E117` | C2 H300 (λ_min 0.5) on the maps v2, with the floor, lgm-calib-v3 | done | −8918 | 4/9 | 13.6 | −1.8 / 1.3 / −0.1 / 2.6 / 3.7 / 6.0 | 5.5 / 6.8 | 5.2 |
| C5b | `tau-CBDC-H1000 \| E1D-ES1-E117` | C5 H1000 on the band's centre line: λ_min 0.34 | done | −8894 | 5/9 | 14.5 | 1.9 / 1.1 / −0.3 / 2.2 / 4.2 / −0.9 | 6.1 / 12.7 | 10.7 |
| | `tau-CBDC-H300 \| E1D-ES1-E117` | C5 H300 on the centre line: λ_min 0.24 | done | −8872 | 5/9 | 14.3 | 0.4 / 1.4 / −0.3 / 2.6 / 3.6 / −0.9 | 5.7 / 8.9 | 10.9 |
| C6 | `tau-CBDA-H1000 \| E1D-ES1-E117` | C5 H1000 with Antarctica of ant-calib: λ_min 0.34, zb −600 / 0 m, c_f 0.0194 (λ on BedMachine v4) | done | −8815 | 5/9 | 14.1 | −0.6 / 0.7 / −0.2 / 1.8 / 4.0 / −1.6 | 6.2 / 9.7 | 12.3 |
| | `tau-CBDA-H300 \| E1D-ES1-E117` | C5 H300, Antarctica λ_min 0.24, zb −600 / 0 m, c_f 0.0406 | done | −8827 | 5/9 | 14.1 | −1.4 / 1.0 / −0.2 / 2.2 / 3.2 / −1.2 | 5.7 / 6.8 | 12.2 |
| | `tau-CBDAR-H1000 \| E1D-ES1-E117` | CBDA H1000, Antarctica's c_f 0.01896 calibrated with λ on the runs' refcorr bed | done | −8819 | 5/9 | 14.2 | −0.6 / 0.7 / −0.2 / 1.8 / 4.0 / −1.6 | 6.3 / 9.7 | 12.8 |
| | `tau-CBDAR-H300 \| E1D-ES1-E117` | CBDA H300 with λ_min 0.22, c_f 0.04154 (refcorr bed) | done | −8822 | 5/9 | 14.2 | −1.4 / 1.0 / −0.2 / 2.2 / 3.2 / −1.4 | 5.7 / 6.8 | 12.7 |
| C7 | `tau-CBF-H1000 \| E1D-ES1-E117` | CBDAR H1000 with λ fixed on the file (`taub_lambda_bed = "file"`) in every instance, Antarctica's c_f 0.01943 recalibrated on it | done | −8810 | 5/9 | 14.3 | −0.5 / 0.6 / −0.1 / 1.7 / 4.0 / −2.0 | 6.3 / 8.2 | 13.7 |
| | `tau-CBF-H300 \| E1D-ES1-E117` | CBDAR H300 with λ on the file, Antarctica's c_f 0.04212 | done | −8788 | 6/9 | 14.5 | −1.3 / 0.9 / 0.0 / 2.0 / 3.0 / −2.7 | 5.7 / 4.5 | 15.2 |
| C8 | `tau-CBF-H300-pDA \| E1D-ES1-E117` | C7 H300: NAIS's margin from Dalton et al. (2020), 25-20 ka (a margin patch) | done | −8788 | 5/9 | 14.5 | −1.4 / 0.9 / −0.0 / 2.0 / 3.0 / −2.6 | 5.7 / 4.5 | 15.2 |
| | `tau-CBF-H300-pDE \| E1D-ES1-E117` | C7 H300: Eurasia's margin from DATED-1 (most-credible), 27.5-7.5 ka | done | −8769 | 6/9 | 14.4 | −1.3 / 0.9 / −0.2 / 1.7 / 3.0 / −3.0 | 5.3 / 4.5 | 15.2 |
| | `tau-CBF-H300-pDAE \| E1D-ES1-E117` | both | done | −8779 | 5/9 | 14.4 | −1.4 / 0.9 / −0.2 / 1.8 / 3.0 / −2.7 | 5.3 / 4.5 | 15.2 |

The far field passes on every run (D1: 5.0) and log L rises with H_N on both Earths; Barents fails on all (too much
Eurasian ice at 20 ka, 20-21 m against PaleoMIST's 18), Antarctica worsens with H_N. A larger H_N lowers τ_b where the
ice is thin: less MIS 3 ice (NAIS 26 instead of 31 m SLE at H_N 1000) and a τ_b that drops with the deglacial ice.

C3 splits Eurasia's c_f (region maps v2): BarentsKara (the Barents and Kara seas with their islands) and Fennoscandia
(with the British-Irish ice), each calibrated to PaleoMIST's 20 ka volume (6.6 and 11.2 m; lgm-calib-v2, NAIS,
Greenland and Antarctica as before). With one c_f, Fennoscandia carried 15 m and BarentsKara 5-6 m at 20 ka: too high
RSL in the White Sea and on Kola, too low on Franz Josef Land. The split raises BarentsKara's c_f and lowers
Fennoscandia's (H_N 1000: 0.0216 and 0.0106 against 0.0151); online they carry 7.2-7.5 and 12.2-12.5 m. Barents passes
(I 3.9 → 2.2) and C3 at H_N 1000 has the best log L of the Coulomb runs (−8811, against −8801 for PaleoMIST's ice on
the same Earth). H_N 2000 and 3000 (λ_min 1.0, the most the band allows) lose log L again: Antarctica and Hudson Bay's
GNSS worsen, and the domes reach 5.0-5.3 km (NAIS) and 7.2-7.4 km (Antarctica).

C4 takes Greenland's margin from PaleoGrIS 1.0 (Leger et al. 2024) as a margin patch (`margin_patches`,
[backend-plastic.md](../backend-plastic.md#margin-patches)): in region 2 of the maps v2 from 15 to 5 ka, the slices of
`tools/levante/margin_patch.jl` (PaleoGrIS's 13 isochrones, 13.5-6.75 ka, between PaleoMIST's 15 and 5 ka slices);
PaleoMIST's margin elsewhere and at other times, bit for bit (every record: no mask cell differs outside Greenland, and
none inside it before 15 or after 5 ka). C3's run is A; the same set-up with a patch equal to the base (PaleoMIST's own
file in region 2) gives it back bit for bit. PaleoGrIS keeps more shelf ice until 12.5 ka and puts the margin on land
about 1 kyr earlier: Greenland holds 0.18-0.30 m SLE more at 14-12.5 ka, 0.72-0.89 m less at 10-9 ka and 0.14-0.16 m
more at 7.5-6 ka (A 10.52, 9.31, 7.08 m at 14, 10, 6 ka; C4 10.70, 8.59, 7.22 m). Greenland's RSL improves (log L
−858.7 → −854.5, RMSE 21.7 → 21.0 m, I 0.9 → 0.6) and its GNSS do not change (I −2.0); the far field (−2.1) and
Laurentide (−1.8, NAIS's ice the same) lose about as much, so the total stays at −8811 and ESL, the groups passed and the
other ice sheets do not change (NAIS, Fennoscandia, BarentsKara within 0.003 m SLE at every record). The runs part
before 15 ka by up to 0.08 m SLE in Greenland (19 ka) and by 0.125 m of thickness elsewhere from 26 ka: the static
footprint (the base's, or the patch's in its region) gains 37 cells of PaleoGrIS's shelf ice, which changes what the
global gather passes to VILMA2. A run with the patch's slices moved after 0 ka (never used: its footprint only, `pFP`)
equals C4 bit for bit up to 15 ka.

C5 adds the floor of the effective pressure, N = ρ_i g H max((1 − H_f/H)^p, δ) with δ = 0.05 (`taub_delta`;
[backend-plastic.md](../backend-plastic.md)). Without it (C1-C3) the Antarctic marine troughs sat at flotation or
grounded, the solve had several solutions, and the Antarctic mass jumped by up to 2.0 × 10⁶ Gt per 100-yr step, in
alternating pairs. With the floor the jumps are gone (Antarctica, per 100-yr step):

| run | max \|ΔM\| [Gt] | alternating pairs > 5 × 10⁴ Gt | North: max \|ΔM\|, pairs |
|---|---|---|---|
| C3 `tau-CBS-H1000` (no floor) | 2.03 × 10⁶ | 36 | 3.13 × 10⁶, 73 |
| C5 `tau-CBD-H1000` | 3.9 × 10⁵ | 0 | 1.13 × 10⁶, 6 |
| C2 `tau-CB-H300 \| E1D-ES1-E117` (no floor) | 1.56 × 10⁶ | 41 | 1.77 × 10⁶, 37 |
| C5 `tau-CBD-H300` | 2.5 × 10⁵ | 0 | 9.2 × 10⁵, 3 |

The northern maxima of C5 are at 10.3 ka, the deglacial retreat of the margin, not an alternation; Patagonia has no
step above 5 × 10⁴ Gt in any run. The floor keeps the troughs grounded above flotation, so it grounds more ice:
recalibrated at 20 ka (lgm-calib-v3; NAIS, Fennoscandia, BarentsKara to their targets, Greenland and Antarctica from
the present-day fit) the c_f hardly change (H_N 1000, λ_min 0.75: NAIS 0.0085, Fennoscandia 0.0104, BarentsKara 0.0194,
Greenland 0.0221, Antarctica 0.0192; H_N 300, λ_min 0.5: 0.0197, 0.0217, 0.0480, 0.0443, 0.0430), but the Antarctic 20 ka
excess rises from 12.7 to 35.2 m (H_N 1000) and from 11.8 to 26.9 m (H_N 300): both points leave the feasible band
(excess ≤ 20 m), which now spans λ_min 0.1-0.5. Online, Antarctica carries 31.4 and 23.9 m above present at 20 ka
(C3: 14.6, C2 H300: 12.8; PaleoMIST 11.2); the northern sheets are as before (NAIS 76.3, Fennoscandia 11.9-12.2,
BarentsKara 6.9-7.2 m). The far field, passed by every earlier Coulomb run, fails (I 16.6 at H_N 1000, 6.0 at 300: too
much Antarctic ice), Antarctica still fails (4.7, 3.7); ESL drops to 2.1 and 5.2, H_N 1000 the only Coulomb run under the
ESL bar (3). log L drops against C3 (−9194, −8918), 5/9 and 4/9 groups. The floor itself works;
these two points need a smaller λ_min (or a recalibrated Antarctic c_f) to stay in the band.

C5b reruns both H_N on the centre line of the new band, by the rule of C1 (the Antarctic 20 ka excess near 12.5 m),
with λ_min refined in steps of 0.02 (lgm-calib-v3/fine): λ_min 0.34 at H_N 1000 (excess 12.3 m; c_f NAIS 0.0111,
Fennoscandia 0.0121, BarentsKara 0.0360, Greenland 0.0222, Antarctica 0.0210) and 0.24 at H_N 300 (12.4 m; 0.0247,
0.0246, 0.0821, 0.0445, 0.0460); Patagonia takes Fennoscandia's c_f. Online, Antarctica carries 11.9 and 11.0 m above
present at 20 ka (C3: 14.6) and no alternating jumps (max |ΔM| 3.0 and 2.9 × 10⁵ Gt; North 6 and 3 pairs, at the
deglacial retreat). The far field passes again (−0.9) and ESL is back at 10.7-10.9, but the smaller λ_min grounds more ice in
the North: NAIS 81.6 and 80.7 m at 20 ka (C3: 77.3, domes 4.9 and 4.5 km), Hudson Bay's GNSS worsens (12.7 at H_N 1000,
8.9 at 300) and so does Laurentide at H_N 1000 (1.9). Both pass 5/9, as C3, with a lower log L (−8894, −8872 against
−8811): the floor removes the jumps at a cost in log L of 60-80 on the centre line.

Up to C6 the runs took λ on their own reference topography relaxed by local isostasy (`taub_lambda_bed = "relaxed"`),
which on E1D-ES1-E117 is today's bed + r (refcorr), while the calibrations took it on GEBCO / BedMachine. Over Hudson Bay
r is +150 to +264 m, so the runs saw Hudson Bay harder than the calibration (λ up to 0.14 higher at λ_min 0.75), and the
north's c_f calibrated on the softer bed grew too much ice there (coulomb-taub/nais-diag; CBDR, the north recalibrated on
the refcorr bed: −8886, −8868). λ stands for the geology of the bed and must not depend on the Earth: since lambda-file λ
is a fixed field per grid, `z_lambda` of `$KRYOS_DATA/plastic/{grid}_lambda_bed.nc` (today's GEBCO / BedMachine bed +
(ρ_i/3300) H_ice, `tools/levante/plastic_lambda_bed.jl`), read by the backend (`taub_lambda_bed = "file"`, the default)
and by the calibration programs alike. The file reproduces the calibrations' relaxed bed bit for bit, so the north's c_f
of lgm-calib-v4 equal lgm-calib-v3's (H_N 1000, λ_min 0.75: NAIS 0.008536, Fennoscandia 0.010375, BarentsKara 0.019428,
Greenland 0.022060; H_N 300, λ_min 0.5: 0.019699, 0.021662, 0.047985, 0.044263).

C6 (coulomb-taub/ant-calib, ant-runs) recalibrates only Antarctica on C5, by the least present-day χ²/n with the 20 ka
excess 12.5 ± 1 m (and at most a quarter of the Ross, Weddell, Amundsen expansion on the floor of N): λ_min 0.34 at H_N
1000 and 0.22-0.24 at 300, with zb −600 / 0 m. CBDA took λ on BedMachine, CBDAR on the runs' refcorr bed; both lift
log L to −8815 to −8827 (Antarctica 11.6-13.3 m at 20 ka, far field −1.2 to −1.6, ESL 12.2-12.8).

C7 is CBDAR with λ on the file in every instance (`experiments/coulomb-band.sh` set `file`) and Antarctica's c_f
recalibrated on it by C6's rule (lgm-calib-v4/ant: 0.01943 at H_N 1000, 0.04212 at 300, excess 13.5 and 12.6 m; the north
and Patagonia keep their c_f). The run's μ is c_f λ of the file (within the 10⁻⁶ of the output). The north's interior
beds are now as soft as in the calibration: NAIS carries 74.6 and 73.5 m at 20 ka (CBDAR 75.6, 76.0; target 72.5; dome
4.3 and 4.0 km), Antarctica 12.7 and 11.2 m. Hudson Bay's GNSS improves most (I 9.7 → 8.2 at H_N 1000, 6.8 → 4.5 at 300),
Greenland, Barents and the far field gain 0.1-1.3, Laurentide and Fennoscandia lose 0.1-0.2, Antarctica holds (4.0) or
gains (3.2 → 3.0), and ESL is worse (13.7, 15.2 against 12.8, 12.7: less northern ice at the LGM). log L −8810 and −8788: H_N 300 is the best Coulomb run so far (6/9; better than PaleoMIST's ice on the same
Earth, −8801), H_N 1000 passes 5/9. No alternating jumps in Antarctica (max |ΔM| 3.6 and 2.0 × 10⁵ Gt); the North keeps
5 and 3 pairs at the deglacial retreat (10.3 ka), as CBDAR.

C8 (coulomb-taub/dated-margins) takes the dated margin chronologies of FesmData's `IceHistory-Dalton2020` (North
America, 36 isochrones 22.1-0.91 cal ka) and `IceHistory-DATED1` (Eurasia, 25-10 ka every 1 kyr, most-credible) as margin
patches on C7 H300 (`tools/levante/margin_patch.jl`, `experiments/coulomb-dated-margins.sh`): Dalton in region 1 from 25
to 20 ka (its slices 22.1, 21.1, 20.5 ka between PaleoMIST's 25 and 20 ka), DATED-1 in regions 3 and 6 from 27.5 to
7.5 ka (all its slices), and both. PaleoMIST 1.0 is close to both at its own slices (intersection over union of the
masks 0.97-1.0 for DATED-1 at 25, 20, 15, 10 ka, 0.91-0.94 for Dalton at 22-17 ka), so the patches change the ice only
between them: NAIS +0.2 m at 22.5 ka, −0.2 m at 20 ka, Eurasia within 0.05 m; ESL is C7's, I 15.2 (at 25 ka far above
the stack's −135.4 m). DATED-1 improves Fennoscandia (I 0.0 → −0.2), Barents (2.0 → 1.7), the far field and
Fennoscandia's GNSS (5.7 → 5.3): log L −8769 (C7 −8788). Dalton gives nothing (−8788) and its hand-over to PaleoMIST at
20 ka is a step of 1.6 × 10⁶ Gt (20.4 ka; C7's largest 8.6 × 10⁵ at 10.3 ka); Antarctica's I 3.008 against 2.995 is the
one group fewer passed. The early LGM ice that ESL wants (esl-diag's ED-pNE25: −8774, ESL I 5.2) is not in these
chronologies: Dalton 2020 starts at 22.1 ka, DATED-1 equals PaleoMIST at 25 ka.

## Ensembles

τ_b sampled per scale group and node, in ens-v2 the Earth too (shifts of structure 6); members weighted by their
likelihood ([The ice: τ_b ensembles](ice-ensembles.md)).

| id | simulation | changes against | status | log L | passed | RMSE [m] | RMS z | I Laur. | I Greenl. | I Fenn. | I Bar. | I Antarc. | I FF | I ESL |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E1 | `tau-ENS1 \| E3D-B21-6` | D2: sampled τ_b (ens-v1, 25 members) | done, not rescored | | | | | | | | | | | |
| E2 | `tau-ENS2 \| E3D-B21-6-ENS2` | E1: the Earth sampled too (ens-v2, 50 members) | done | | | | | | | | | | | |
| | best log L: `tau-ENS2-E048 \| E3D-B21-6-E048` | | | −9017 | 4/7 | 19.7 | 32.2 | −0.5 | 13.0 | 4.5 | 1.1 | 0.9 | −0.6 | 7.0 |
| | most groups passed: `tau-ENS2-E036 \| E3D-B21-6-E036` | | | −9063 | 6/7 | 18.5 | 27.3 | 1.3 | 6.2 | 6.8 | 1.2 | 2.9 | 0.4 | 2.2 |
| | median of the 50 | | | −9474 | | | | | | | | | | |
| E4 | `tau-JE-<m> \| E1D-ES1-<Earth>` | C8 pDE with PaleoGrIS: H_N, λ_min, c_f, the LGM timing and 53 of wave A's 1-D Earths sampled together ([The joint ice + Earth ensemble](joint-ensemble.md)) | smoke test done (3 members), wave 1 (200) set up | | | | | | | | | | | |
| E3 | the final ensemble: our τ_b on the Earths of the Earth search | E2: the Earth chosen without our ice, a wide τ_b prior | proposed | | | | | | | | | | | |

The best members beat P2 but not the imposed references (I1, I5); Fennoscandia limits most of them (4 of the 50 pass
it on structure 6, best m013, I 0.2).

## The Earth search

Imposed ice on sampled Earths ([The Earth search](earth-search.md)), judged on 8 groups: the six RSL groups, GNSS
Fennoscandia and Hudson Bay; ESL is reported, not judged ([ESL in the Earth search](earth-search.md#esl-in-the-earth-search)). Earth: lithosphere [km], log10 η UM / TZ / LM1 / LM2 [Pa s].

| id | simulation | Earth | log L | passed | I L / G / F / B / A / FF | I GNSS F / HB | I ESL |
|---|---|---|---|---|---|---|---|
| S0 | `ice-PMIST \| E1D-ES1-PMIST` | Gowan's, in the search's parameters (the check: I1 again) | −8860 | 6/8 | −3.5 / 6.8 / 1.0 / 0.7 / 0.4 / −0.6 | 7.4 / 3.5 | 10.2 |
| | `ice-GLAC1D \| E1D-PMIST` | I5 on this scale | −8835 | 7/8 | −3.0 / 5.7 / 0.5 / 1.2 / −0.2 / −0.8 | 7.1 / 2.7 | 7.6 |
| S1 | `ice-PMIST \| E1D-ES1-E<nnn>` | wave A, 200 1-D Earths: 20.2 % of the space not ruled out | | | | | |
| | best log L: `ice-PMIST \| E1D-ES1-E117` | 111, 20.13 / 21.75 / 23.15 / 21.69 | −8801 | 6/8 | −3.2 / 2.0 / 1.6 / 1.4 / 0.2 / −0.7 | 5.7 / 6.4 | 10.2 |
| | most groups: `ice-PMIST \| E1D-ES1-E017` | 87, 20.56 / 20.61 / 21.68 / 21.02 | −8906 | 7/8 | −0.7 / 2.5 / 0.6 / 1.0 / 0.8 / −1.6 | 4.9 / 2.1 | 10.7 |
| S2 | `ice-GLAC1D \| E1D-ES1-E<nnn>` | wave A: 13.5 % not ruled out; under either ice 28.8 %, both 4.9 % | | | | | |
| | best log L: `ice-GLAC1D \| E1D-ES1-E165` | 67, 20.77 / 20.08 / 22.54 / 21.88 | −8709 | 6/8 | −3.5 / 2.6 / 0.8 / 3.3 / −1.4 / −4.4 | 3.6 / 2.9 | 7.9 |
| | all groups: `ice-GLAC1D \| E1D-ES1-E033` (14 Earths pass all 8) | 106, 20.84 / 21.43 / 21.08 / 22.08 | −8746 | 8/8 | −2.4 / 1.7 / 2.2 / 0.7 / −1.1 / −4.4 | 0.5 / 1.1 | 8.0 |
| S3 | `ice-<ICE> \| E1D-ES2-E<nnn>` | wave B, the radial profiles of the 50 3-D Earths, under each ice | done, scored (in the emulator of S4) | | | | | |
| S4 | `ice-PMIST \| E3D-ES2-E<nnn>` | wave B, 50 3-D Earths: 39 % of the 13-D space not ruled out (all 100 runs; 42 % on 96) | | | | | |
| | best: `ice-PMIST \| E3D-ES2-E014` | 55, 19.68 / 22.47 / 21.13 / 21.09; α 0.74 (Lloyd-like) | −8836 | 7/8 | −2.5 / 2.6 / 1.4 / 0.8 / 2.6 / −1.3 | 5.6 / 1.0 | 10.7 |
| | `ice-PMIST \| E3D-ES2-E043` | 157, 20.98 / 20.60 / 22.09 / 21.02; α 0.59 | −8966 | 7/8 | 2.9 / 1.4 / 4.9 / 1.3 / 0.6 / −1.7 | 2.3 / 0.9 | 10.5 |
| S5 | `ice-GLAC1D \| E3D-ES2-E<nnn>` | wave B: 45 % not ruled out; under both 33 % (36 % on 96 runs) | | | | | |
| | best: `ice-GLAC1D \| E3D-ES2-E014` | as above | **−8717** | 7/8 | −3.0 / 3.2 / 0.7 / 2.0 / −0.3 / −3.2 | 7.8 / 1.2 | 8.1 |
| | `ice-GLAC1D \| E3D-ES2-E036` | 124, 21.39 / 20.64 / 21.24 / 21.75; α 0.08 (nearly 1-D) | −8893 | 7/8 | 2.1 / 6.7 / 3.8 / 0.7 / −0.8 / −4.5 | 0.5 / 2.1 | 8.1 |
| S6 | `ice-<ICE> \| E3D-ES2-E1<nn>` | wave B2, 30 3-D Earths in two widened boxes | running | | | | | |

No 3-D Earth passes all 8 groups. E014 is the only one that beats Gowan's Earth under both ice histories. ESL, reported
only: PaleoMIST's ice scores 10.0-11.0 on every Earth, GLAC-1D's 7.3-8.4 ([Data and scores](data-scores.md#the-bar)).

## Misfit common to every simulation

West Greenland (Sisimiut, −100 m), the St. Lawrence (+16 to +56 m) and +10 to +20 m in the far field (Barbados, the Great
Barrier Reef, the Bonaparte Gulf) fail for every good series, the published ones included. Lateral Earth structure
cannot remove them (Milne et al. 2018 for west Greenland; Pan et al. 2022 for the far field, where its effect at
Barbados has the wrong sign), so they point to the ice: too little LGM ice in the far field, too thick local ice in the
St. Lawrence.

## Reproducing the tables

Each row is made by a script of [experiments/](https://github.com/fesmc/kryos/tree/main/experiments) (refs-ice.sh,
refs-tau.sh, ens.sh, earth-search.sh) and scored by `experiments/scores.sh`. Every imposed-ice case writes
`provenance.txt` (its inputs, the commits of kryos and VILMA2, the md5 of the par file and the binary). The figures come
from `tools/levante/plastic_sim_figures.jl` (the ladder and Bagge figures) and `tools/levante/plastic_earth_emulate.jl`
(the Earth search).

The runs are on AWI's albedo:

| id | where |
|---|---|
| P1, P2 | FesmData `IceHistory-PaleoMIST` (`dz_sl`), `IceHistory-GLAC1D` (topography) |
| I1-I8 | `/albedo/work/user/alrobi001/kryos-work/kryos_data/refcorr/<ice>_<earth>`; scores `kryos-work/refs-test/ens/iter_01/scores.csv` |
| T1, D1 | `/albedo/work/user/alrobi001/kryos-work/runs/pm-online/gowan`, `.../pdfac` |
| C1, C2 | `/albedo/work/projects/p_forclima/kryos-work/coulomb-taub/coupled/runs/tau-CB-H<H_N>_<Earth>`; scores `.../coupled/score/iter_01/scores.csv` |
| C3 | `/albedo/work/projects/p_forclima/kryos-work/coulomb-taub/split/runs/tau-CBS-H<H_N>_E1D-ES1-E117`; scores `.../split/score/iter_01/scores.csv`; calibration `.../coulomb-taub/lgm-calib-v2` |
| C4 | `/albedo/work/projects/p_forclima/kryos-work/coulomb-taub/margin-patch/runs/tau-CBS-H1000-{pPG,pPM,pFP}_E1D-ES1-E117`; scores `.../margin-patch/score/iter_01/scores.csv`; the patch `.../kryos-work/margins/NH-32KM/NH-32KM_IceHistory-PaleoGrIS-patch.nc` |
| C5 | `/albedo/work/projects/p_forclima/kryos-work/coulomb-taub/delta/runs/tau-CBD-H<H_N>_E1D-ES1-E117`; scores `.../delta/score/iter_01/scores.csv`; calibration `.../coulomb-taub/lgm-calib-v3`; jumps `.../delta/jumps/jumps.md`; volumes `.../delta/coupled/coupled.md` |
| C5b | `.../coulomb-taub/delta/runs/tau-CBDC-H<H_N>_E1D-ES1-E117`; scores as C5; calibration `.../lgm-calib-v3/fine`; jumps `.../delta/jumps-centre/jumps.md`; volumes `.../delta/coupled-centre/coupled.md` |
| C6 | `.../coulomb-taub/ant-runs/runs/tau-{CBDA,CBDAR}-H<H_N>_E1D-ES1-E117`; scores `.../ant-runs/score/iter_01/scores.csv`; calibration `.../coulomb-taub/ant-calib/fine` (CBDA), `.../ant-calib/fine-ref/report.md` (CBDAR); jumps `.../ant-runs/jumps/jumps.md`; volumes `.../ant-runs/coupled/coupled.md` |
| C7 | `.../coulomb-taub/lambda-file/runs/tau-CBF-H<H_N>_E1D-ES1-E117`; scores `.../lambda-file/score/iter_01/scores.csv`; calibration `.../coulomb-taub/lgm-calib-v4` (north), `.../lgm-calib-v4/ant/pick.txt` (Antarctica); jumps `.../lambda-file/jumps/jumps.md`; volumes `.../lambda-file/coupled/coupled.md`; λ files `/albedo/work/projects/p_forclima/kryos_data/plastic/{grid}_lambda_bed.nc` |
| C8 | `.../coulomb-taub/dated-margins/runs/tau-CBF-H300-{pDA,pDE,pDAE}_E1D-ES1-E117`; scores `.../dated-margins/score-esl/iter_01/scores.csv` (ESL rescored; `score` has the old comparator); patches `.../dated-margins/margins`; jumps `.../dated-margins/jumps/jumps.md`; volumes `.../dated-margins/coupled/coupled.md`; ESL by age `.../dated-margins/diag/esl.csv`; products (staging) `/albedo/work/projects/p_forclima/kryos-work/fesmdata-staging/fesmdata` |
| E4 | `/albedo/work/projects/p_forclima/kryos-work/joint-ens/<set>/runs/tau-JE-<m>_E1D-ES1-<Earth>` (sets `smoke`, `wave1`); design `.../<set>/design.csv`; scores `.../<set>/score/iter_01/scores.csv`; emulator `.../<set>/emulate`; Earths `.../joint-ens/earths.csv`, `refcorr/`; calibration `.../joint-ens/calib` |
| D2, E2 | `/albedo/work/user/alrobi001/kryos-plastic/kryos/runs/ens-v2/iter_01/<m>` |
| D3 | `/albedo/work/user/alrobi001/kryos-work/runs/res16/m000-NH16` |
| E1 | `/albedo/work/projects/p_forclima/kryos-work/levante/ens-v1` (copied from Levante) |
| S0-S5 | `/albedo/work/user/alrobi001/kryos-search` (`kryos_data/refcorr/`, `score/ens/`); the emulators in `/albedo/work/projects/p_forclima/kryos-work/earth-search/emulateB` |
| S6 | `/albedo/work/projects/p_forclima/kryos-work/earth-search/waveB2` |

The ESL numbers since 2026-11-18 (all other scores unchanged) are from the rescore of every score folder above with the
same config and error model: `/albedo/work/projects/p_forclima/kryos-work/rescore-esl/<folder>/score/iter_01/scores.csv`
(folders `refs-test` for P, I, T, D, E; `earth-search` for S0-S5; `coupled`, `split`, `margin-patch`, `delta`,
`ant-runs`, `lambda-file` for C1-C7; `nais-diag`, `delta-north`, `esl-diag`, `paleogris`, `rtest` for the tests), the
emulators in `.../rescore-esl/earth-search/emulate*`, the comparison old / new in `.../rescore-esl/summary.md`.

The runs under `/albedo/work/user/alrobi001` move to the project's shared folder
(`/albedo/work/projects/p_forclima/kryos-work`) once their last jobs end.

## Earlier runs, not on this scale

| run | in this language | note |
|---|---|---|
| B1, B3 (2026-10-31) | `ice-PMIST \| E1D-PMIST`, `ice-PMIST \| E3D-B21-6` | VILMA2 523b660, 931 kg m⁻³, 2.5 kyr steps, no r; RMS z of B1 30.9 |
| plastic online test (2026-11-09) | about `tau-PMIST \| E1D-PMIST` | RMS z 25.9 against PaleoMIST 29.5 |
| esmda-v1, -v2 (Levante) | calibrated τ_b on `E1D-PMIST` | ES-MDA calibrations, retired |
