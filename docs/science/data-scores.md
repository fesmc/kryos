---
title: "Data and scores"
---

What every simulation is scored against, and how. All data come from
[FesmData](https://github.com/fesmc/FesmData) products; the scores from `plastic_ensemble.jl scores`, the same code and
error model for every simulation, published reconstructions included.

## Data

| group | data | n | judges |
|---|---|---|---|
| Laurentide | GAPSLIP 2.1 RSL index points | 892 | yes |
| Greenland | GAPSLIP 2.1 | 276 | yes |
| Fennoscandia | GAPSLIP 2.1 (mainland and the British Isles) | 718 | yes |
| Barents | GAPSLIP 2.1 (Svalbard, Franz Josef Land, Novaya Zemlya, the Kara Sea) | 370 | yes |
| Antarctica | GAPSLIP 2.1 | 47 | yes |
| far field | GAPSLIP 2.1 | 2248 | yes |
| ESL | the eustatic sea level of the Spratt and Lisiecki (2016) stack, 25-10 ka, against the run's ice-volume sea level | 7 | yes |
| GNSS Fennoscandia | present-day uplift rates, BIFROST2015 (Kierulf et al. 2021) | 164 | yes |
| GNSS Hudson Bay | Schumacher et al. (2018), elastic-corrected, and MIDAS stations (Blewitt et al. 2016) | 61 | yes |
| GNSS Greenland | GNET (Berg et al. 2024) | 55 | reported only |
| GNSS Antarctica | Martín-Español et al. (2016), elastic-corrected | 67 | reported only |

The RSL points are those of age ≤ 25 ka: GLAC-1D's own RSL ends at 26 ka, so the published reconstructions can be scored
on the same points. GNSS rates are compared with the run's rate of bed change at its last record (elastic-corrected
data; VILMA2's centre-of-mass frame against ITRF/IGS, a few tenths of a mm/yr apart, not corrected). In Greenland and
Antarctica present-day uplift is dominated by the elastic response, Little Ice Age loading and transient rheology, so
those groups are reported, not judged.

## The error model

Per data group g, the residuals d = model − observation have the covariance

R = diag(σ_i²) + σ_m,g² C,  C_ij = exp(−dist_ij / L_x,g) · exp(−|Δt_ij| / L_t,g)

with σ_i the data errors and a model error σ_m,g correlated over distance L_x,g and time L_t,g. The model error covers
what no simulation can fit (the resolution, the ice outside the reconstruction, the solver). Its parameters are fitted
from the residuals of the best members of an ensemble (ens-v2 for the RSL groups; wave A of the Earth search for GNSS:
Fennoscandia 0.7 mm/yr over 750 km, Hudson Bay 1.5 mm/yr over 500 km) and then held fixed for every simulation.

## The scores

- **log L**, the Gaussian log-likelihood −½ d′R⁻¹d − ½ log|R| summed over the six RSL groups: higher is better.
  Differences of a few tens are large. This is the one scale of all the tables.
- **I**, the implausibility of a group, (χ² − n)/√(2n) with χ² = d′R⁻¹d: about 0 when the residuals are as large as the
  error model expects, ≤ 3 passes.
- **passed**, the number of judging groups within the bar: of 9 (six RSL groups, ESL, GNSS Fennoscandia and Hudson
  Bay), of 7 for a simulation without GNSS rates (the published reconstructions and runs before 2026-11-12).
- **RMSE** [m] and **RMS z** (the root mean square of the residuals normalized by σ_i, the score of the first
  comparisons) are reported for reference only. A few points with tiny σ_i dominate RMS z, Greenland above all.

## The bar

A group passes at I ≤ 3. Where both published reconstructions fail a group, the cut is raised to the better of their
two scores: Greenland 7.9 (PaleoMIST). Every other RSL group and the GNSS groups keep 3. ESL depends on what is judged:

- **Our ice** (the ensembles): the plain cut 3. ESL is an independent constraint on the ice, which the publications were
  not built against.
- **An Earth under imposed ice** (the Earth search): PaleoMIST's own score, 7.04. Under imposed ice ESL is a property of
  the ice, not of the Earth (5.9-6.1 under PaleoMIST's ice on every Earth), so a plain cut would rule Earths out for the
  ice's sake.

A simulation is *as good as the publications* when it passes every group, and *better* when it does so with a higher
log L. The Earth search uses the same bar to rule Earths out ([The Earth search](earth-search.md)).
