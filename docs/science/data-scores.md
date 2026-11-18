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
| ESL | the Spratt and Lisiecki (2016) sea-level stack, 25-10 ka every 2.5 kyr, against the run's barystatic sea level ([below](#esl)) | 7 | yes |
| GNSS Fennoscandia | present-day uplift rates, BIFROST2015 (Kierulf et al. 2021) | 164 | yes |
| GNSS Hudson Bay | Schumacher et al. (2018), elastic-corrected, and MIDAS stations (Blewitt et al. 2016) | 61 | yes |
| GNSS Greenland | GNET (Berg et al. 2024) | 55 | reported only |
| GNSS Antarctica | Martín-Español et al. (2016), elastic-corrected | 67 | reported only |

The RSL points are those of age ≤ 25 ka: GLAC-1D's own RSL ends at 26 ka, so the published reconstructions can be scored
on the same points. GNSS rates are compared with the run's rate of bed change at its last record (elastic-corrected
data; VILMA2's centre-of-mass frame against ITRF/IGS, a few tenths of a mm/yr apart, not corrected). In Greenland and
Antarctica present-day uplift is dominated by the elastic response, Little Ice Age loading and transient rheology, so
those groups are reported, not judged.

### ESL

The data: `sl_short` of the stack (FesmData `SeaLevel-Spratt2016`, 1 sigma `sl_short_sd`) relative to its own present
(+8.5 m at 0 ka), linear in time between its 1-kyr records at 25, 22.5, ..., 10 ka (22.5 ka is the mean of 22 and 23 ka).
The model: VILMA2's barystatic sea level `bsl` of the run's global instance (`geo_glob_vilma.nc`; `vilma_out.nc` of a
standalone run) relative to its last record (0 ka), linear between its records. `bsl` is the ocean mean of C·rsl,
∫C rsl dA / ∫C dA over the ocean of the time (VILMA2 `update_bsl`); by the mass balance of the sea-level equation it is the
melt water of the grounded-ice change plus the fill of newly flooded area. The ocean that fills the cavity of lost marine
ice does not raise it, so in effect the ice above flotation counts, with the change of the ocean area: the
ice-volume-equivalent sea level the stack stands for. The error model of the group has only a time term (σ_m 5 m, L_t 10
kyr). `[esl] model_var` of the scoring config selects the model's series (`tools/levante/plastic_ensemble.jl`).

Until 2026-11-18 the model's series was `bsl_ive`, −(ρ_i/ρ_w) of the whole grounded-column change over the ocean area, the
ice below flotation included (about 10 m more ice at the LGM than `bsl`), and the data were taken at their nearest record
(22, 17, 12 ka for 22.5, 17.5, 12.5 ka). `model_var = "bsl_ive"` gives that comparator back (with the data linear in
time). With `bsl` every simulation reads 10-11 m less ice at 20 ka and its ESL I rises by about 4-5 (where the run is
short of the stack); the interpolation alone moves I by up to ±0.8.

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
- **An Earth under imposed ice** (the Earth search): PaleoMIST's own score, 7.22. Under imposed ice ESL is a property of
  the ice, not of the Earth (10.0-11.0 under PaleoMIST's ice on every Earth, 7.3-8.4 under GLAC-1D's), so a plain cut
  would rule Earths out for the ice's sake. PaleoMIST's own ESL (`SeaLevel-PaleoMIST`) is its whole ice volume, below
  flotation included, while the runs are scored on `bsl`: since the rescore of 2026-11-18 both ice histories on our
  solver score worse than the bar, so ESL rules out every Earth
  ([The Earth search](earth-search.md#esl-after-the-rescore)). Open: a bar on the same comparator.

A simulation is *as good as the publications* when it passes every group, and *better* when it does so with a higher
log L. The Earth search uses the same bar to rule Earths out ([The Earth search](earth-search.md)).
