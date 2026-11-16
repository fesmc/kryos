---
title: "Reconstructing the last glacial cycle"
---

Ice sheets and sea level of the last glacial cycle (80-0 ka) reconstructed with
[kryos](https://fesmc.github.io/kryos/): a perfectly-plastic ice sheet computed from basal shear stress on a deforming,
laterally varying solid Earth, scored against relative sea-level (RSL) index points, the eustatic sea-level stack and
present-day GNSS uplift rates. These pages are the material of a paper in preparation: the question, the method, the
results so far and what is still open. How to run things is in kryos's
[how-to pages](https://fesmc.github.io/kryos/howto/imposed-ice.html).

## The question

Can an ice sheet computed from basal shear stress (τ_b) on a deforming Earth reproduce the sea-level record of the last
glacial cycle better than the published reconstructions? And which ingredient makes the difference: the ice, the Earth,
or the solver?

The bar is not a perfect fit. It is the two published reconstructions whose misfits were accepted when they were
published: PaleoMIST 1.0 (Gowan et al. 2021; a 1-D Earth and a τ_b tuned by hand) and GLAC-1D (Tarasov et al. 2012,
Briggs et al. 2014; a calibrated ice-sheet model). A reconstruction is an improvement when it is as good as they are in
every data group, and better overall ([Data and scores](data-scores.md)).

## The approach

Three steps, each removing one unknown before the next:

1. **The solver.** Put the published ice histories, unchanged, on our solid-Earth model (VILMA2) with the published
   Earth. Any change of score is the solver's. Then vary the Earth under the same ice: the change is the Earth's
   ([Methods](methods.md), the ladder).
2. **The Earth.** Search the Earth structure broadly under both published ice histories, imposed, with an emulator to
   cover what was not run. Two ice histories, so that an Earth effect is not an ice error in disguise
   ([The Earth search](earth-search.md)).
3. **The ice.** Compute our own ice from τ_b on a handful of Earths the search leaves plausible, sample τ_b widely and
   weight the histories by their fit ([The ice: τ_b ensembles](ice-ensembles.md)). Because the Earths were chosen
   without our ice, the gain over the published reconstructions splits into the ice's and the Earth's.

## Status (2026-11-16)

| step | state | next |
|---|---|---|
| 1. Solver and published Earths | done: both published ice histories on VILMA2, Gowan's and Bagge's Earths | - |
| 2. Earth search, wave A (200 1-D Earths) | done, emulated | - |
| 2. Earth search, wave B (50 3-D Earths) | 96 of 100 runs done and emulated; the last 4 queued | the final pass at 100 runs |
| 2. Earth search, wave B2 (30 3-D Earths, a widened box) | sampled, stability test running | the runs (2-4 days), then the emulator over A, B and B2 |
| 3. τ_b ensembles on Bagge's structure 6 (ens-v1, ens-v2) | done | - |
| 3. The final ensemble on the chosen Earths | designed (proposed) | choose 4-6 Earths after wave B2 |

## Findings so far

- **Our solver with the published ice beats both publications.** PaleoMIST's ice on VILMA2 with Gowan's Earth scores
  log L −8860 against −9311 for the publication, GLAC-1D's ice −8835 against −9117 ([Results](results.md)).
- **Bagge's radial profile fits as well as Gowan's; his lateral structure does not.** With Bagge's lateral structure
  Fennoscandia fails, the more the stronger the lateral contrast, under either ice history.
- **The data constrain the upper mantle, not the lower mantle.** Wave A rules out 71 % of the 1-D Earths under either ice history, 95 % under both. The
  lithosphere and the upper mantle are constrained, the transition zone and the lower mantle are not, and the upper
  mantle and the transition zone trade off along a ridge ([The Earth search](earth-search.md)).
- **Lateral structure usually costs fit, with one clear exception.** Most 3-D Earths of wave B score worse than their
  own radial profile. One, E014 (thin lithosphere, soft upper mantle, stiff transition zone, Lloyd-like lateral
  pattern), fits better than every published Earth under both ice histories and passes 8 of 9 groups under each. It sits
  at the edge of the prior, which wave B2 tests.
- **Our own ice does not yet reach the imposed references.** The best ens-v2 members beat GLAC-1D overall but fail
  Fennoscandia, on Bagge's structure 6, which the Earth search has since shown to be a poor Earth.
- **Some misfit is common to every reconstruction.** West Greenland, the St. Lawrence and parts of the far field fail
  for every simulation, the published ones included. Earth structure cannot explain them (Milne et al. 2018; Pan et al.
  2022), so they point to the ice.

## Open questions

- **Which Earth?** Whether the best region lies outside the prior (wave B2), and how many Earths carry into the final
  ensemble.
- **Far field against ESL.** Ice histories that fit the ESL stack overshoot the far-field RSL, and the reverse.
- **The final ensemble.** The width of the τ_b prior, and whether the RSL of a candidate ice history can be predicted
  from step responses of each Earth ([The ice: τ_b ensembles](ice-ensembles.md#the-final-ensemble-proposed)).

## The pages

| page | in a paper |
|---|---|
| [Methods](methods.md) | the model, the simulations, the naming and the ladder |
| [Data and scores](data-scores.md) | the data, the likelihood and the bar |
| [The Earth search](earth-search.md) | method and results of the Earth search |
| [The ice: τ_b ensembles](ice-ensembles.md) | method and results of our ice, and the final ensemble |
| [Results](results.md) | every scored simulation in one table |
