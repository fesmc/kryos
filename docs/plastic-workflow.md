---
title: "Glacial-cycle reconstructions: the workflow"
toc: true
---

The working procedure behind [Simulations and scores](simulations.md): what we are trying to show, how each simulation
is made and scored, and where to continue. Status as of 2026-11-13.

## The question

Can an ice sheet computed from basal shear stress (PlasticIceSheet, online in kryos) on a deforming Earth (VILMA2)
reproduce the relative sea-level (RSL) record of the last glacial cycle better than the published reconstructions? And
which ingredient makes the difference: the ice, the Earth, or the solver?

The bar is not a perfect fit. It is the two published reconstructions whose misfits were accepted when they were
published: PaleoMIST 1.0 (Gowan et al. 2021; 1-D Earth, a hand-tuned τ_b, ICESHEET offline at 5 km) and GLAC-1D
(Tarasov et al. 2012, Briggs et al. 2014; a calibrated ice-sheet model). A simulation is an improvement when it scores
better than they do, group by group.

## The concepts

- **Imposed ice** (`ice-<src>`): a published ice history forced onto VILMA2. The ice does not respond to the Earth. This
  isolates the Earth and the solver.
- **Online ice** (`tau-<src>`): PlasticIceSheet computes the ice every step from τ_b and the margins, on the bed that
  VILMA2 deforms. This is our reconstruction method.
- **Present-day correction r**: a glacial-cycle run must end with today's bed. `gia_refcorr.sh` iterates the reference
  topography (today's + r) until the 0 ka bed of the run matches today's (to 1.5 m). r belongs to one ice history and
  one Earth. A run with a borrowed r (`R0`, the Earth-only ensemble) is one pass instead of 3-5.
- **The ladder**: each simulation changes one thing against its parent, so the change of score is that thing's
  (P → I → T → D → E; [Simulations and scores](simulations.md#the-ladder)).
- **Scores**: the log-likelihood log L over all data, and per data group the implausibility
  I = (χ² − n)/√(2n) under an error model with model error and spatial correlation, refitted from ens-v2
  (`errmodel.toml`). A group passes at I ≤ 3, or at the best published score where both publications fail (the bar). The
  ESL group (the Spratt et al. 2016 stack, 25-10 ka) keeps the plain cut.
- **Data**: GAPSLIP 2.1 RSL index points of age ≤ 25 ka in six groups (Laurentide, Greenland, Fennoscandia, Barents,
  Antarctica, far field), the ESL stack, and GNSS uplift rates (Fennoscandia, Hudson Bay; Greenland and Antarctica
  reported only). All come from FesmData products.
- **Ensembles**: τ_b is sampled (scales per region and domain factors), the members are weighted by their likelihood
  with tempering, and the Earth is either fixed (Bagge structure 6) or sampled.

## The Earth search, step by step

Which Earth do the data want, independent of our ice? The expert 3-D viscosity fields (Bagge, Lloyd, Pan) disagree by about
a factor of 10 in the lower mantle and the lithosphere, and Bagge's lateral structure fails Fennoscandia (the I rows of
[Simulations and scores](simulations.md#imposed-ice-on-our-solver)). So the search samples Earths broadly, scores them under
imposed ice, and uses an emulator to cover the space that was not run. The design and the decisions:
`docs/dev/earth-search-design.md`. Results: [Simulations and scores](simulations.md#the-earth-search).

**0. Held fixed in every wave**

| piece | what it is | what it is for |
|---|---|---|
| two imposed ice histories, PaleoMIST and GLAC-1D | every Earth runs under both | an Earth effect must not be an ice error in disguise |
| present-day correction r | per Earth and ice, by `gia_refcorr.sh` | the run's bed today matches the observed bed, so a score measures the Earth, not a topography offset |
| the scores | log L and I per group: 6 RSL groups, GNSS uplift rates in Fennoscandia and Hudson Bay, ESL | GNSS constrains today's upper mantle directly; ESL mainly checks the ice |
| the bar | I ≤ 3, Greenland ≤ 7.9, ESL ≤ 7.04 (the best published score where both publications fail) | an Earth is plausible when it is as good as the published reconstructions in every group |

**1. Wave A: 1-D Earths (done 2026-11-13).** 200 radial profiles of 5 parameters (lithosphere thickness; log10 viscosity of
the upper mantle, the transition zone, the upper and the lower part of the lower mantle), a Latin hypercube over wide
ranges. 400 runs of 3-5 passes of about 5 min. *For*: a cheap, broad map of the radial profiles the data allow.

**2. Emulator and cut of wave A (done).** Per data group and ice history, a Gaussian process predicts the score
(log χ²/n) and its uncertainty anywhere in the 5-D space (`plastic_earth_emulate.jl`). A point is ruled out when the emulator
is confident (3 sd) that it fails the bar in a group. *For*: the space wave B samples from (29 % of the 1-D space, under
either ice), and what 1-D Earths alone can say (the lithosphere and the upper mantle are constrained; the transition zone
and the lower mantle are not).

**3. Wave B: 3-D Earths (running).** 50 Earths: the radial part from the space left by wave A, plus 8 lateral parameters,
the amplitude α of the experts' mean lateral pattern, two coefficients c₂, c₃ of the expert basis (Lloyd against Pan,
Bagge's reduction factor) and five regional upper-mantle shifts (Fennoscandia, Barents, Laurentide, Greenland, West
Antarctica). Per Earth and ice two runs: its radial profile alone (1-D, iterated to its r), then the 3-D Earth in one pass
from that r (`R0`). *For*: does lateral structure fit better than the best 1-D Earths, and where? The pair 3-D minus 1-D
is the effect of the lateral structure on the same radial profile.

**4. Emulator and cut of wave B (next).** As step 2, in all 13 parameters. Each Earth also gets its distance from the
expert fields (radial and lateral), a second axis next to the fit, not a weight in it. *For*: the plausible 3-D Earth
space, and the front of the fit against the closeness to the experts. An Earth far from the experts that fits much better
is a finding; one that fits a little better is not.

**5. Wave C: our own ice (after B).** A τ_b ensemble (as ens-v2) on a handful of plausible Earths from wave B, spread
over the space left. *For*: the target, beating the publications with our ice, and splitting the gain between the ice
and the Earth, since the Earths were chosen without our ice.

## Reproducing a simulation

1. Install and build kryos ([Getting started](getting-started.md)): `configme install kryos`, then
   `make kryos plastic=1 vilma=1 yelmo=0`. Run `tools/levante/setup.sh` for the data links, maps, runme and the Julia
   project of the tools.
2. Run the row's script in [experiments/](https://github.com/fesmc/kryos/tree/main/experiments) (its README maps ids to
   scripts). Example, PaleoMIST's ice on Gowan's 1-D Earth:

   ```bash
   experiments/refs-ice.sh I1
   ```

3. Score: `experiments/scores.sh <ensemble folder>`. It scores every member and every reference of the folder's config,
   and writes `iter_01/scores.csv`.
4. Update the table in [Simulations and scores](simulations.md) from `scores.csv` (the window ≤ 25 ka rows).

Every `gia_refcorr.sh` case writes `provenance.txt` (inputs, the commits of kryos and VILMA2, md5 of the par file and
the binary). A kryos run records its build. The results so far were made on albedo (AWI); the paths of the table are
under `/albedo/work/user/alrobi001`.

**Reproducibility check (2026-11-13).** I1 made on the older VILMA2 (fe17b6d, old namelist keys) and again on VILMA2
main (1f01ca7) converge in the same iterations, and their r agree to 4 µm. The references made on the older stack (I1-I4,
T1, D1) therefore stand for main.

## Tools

| Tool | Does |
|---|---|
| `tools/levante/gia_refcorr.sh` | an imposed-ice case: VILMA2 standalone, r iterated, reference topographies written |
| `tools/levante/plastic_ensemble.sh`, `.jl` | ensembles: sampling, the runs, weights, scores, error model, the Earth basis and Earth samples |
| `tools/levante/plastic_ensemble_misfit.jl` | maps and tables of the misfits, the shared misfit of the best members |
| `tools/levante/plastic_earth_basis.jl` | the PCA basis of the expert viscosity fields, fields from it, the 1-D par files of the Earth search |
| `tools/levante/gia_refcorr_pool.sh` | many small `gia_refcorr.sh` cases in few jobs (the Earth search) |
| `tools/levante/plastic_earth_emulate.jl` | the emulators and the cut of the Earth search |
| `tools/levante/plastic_sim_figures.jl` | the figures of [Simulations and scores](simulations.md) |
| `tools/levante/pm_compare.jl` | the RSL, GNSS and ESL data and the model series at the sites |

## Findings so far

- **Our solver with PaleoMIST's ice beats the publication.** I1 (`ice-PMIST | E1D-PMIST`) has log L −8860 against −9311
  for P1. The changes are the solver (SELEN → VILMA2), the ice density (931 → 910 kg m⁻³), linear time interpolation
  instead of 2.5 kyr steps, and the run's own r.
- **GLAC-1D's ice does better still on our solver.** I5 −8835 and I8 (nn9927) −8825, both passing everything but ESL
  (I 3.3 and 5.0).
- **Bagge's lateral structure hurts in Fennoscandia, his radial profile does not.** On PaleoMIST's ice, the
  Fennoscandia I grows with the lateral contrast of Bagge's structures (1: 5.7 → 6: 15.2), while the radial mean of
  structure 6 as a 1-D Earth scores 3.3. The online ensemble ens-v2 sits on structure 6, which explains part of its
  Fennoscandia misfit.
- **No online member reaches the imposed references yet.** The best ens-v2 members (m048 −9017, m036 −9063) match
  GLAC-1D (P2 −9117) but not I1 or I5.
- **ESL**: PaleoMIST's ice is about 7 m short of the stack at 20 ka (I 5.9); GLAC-1D's ice fits it (P2 0.0).

## Open questions

- **Which Earth?** The Earth search ([above](#the-earth-search-step-by-step)): wave A done, wave B running.
- **Structural misfit.** West Greenland, the St. Lawrence and parts of the far field fail for every simulation,
  published ones included (the shared-misfit maps of `plastic_ensemble_misfit.jl`).
- **Far field against ESL.** Members that fit the ESL stack overshoot the far-field RSL and the reverse.
- **Resolution.** 16 km in the north changes the scores little (D3 against D2) but costs 2 m of ESL. 32 km stays the
  default for ensembles, and 16 km is ready (par `..._nh16.nml`, its maps and reference topographies).

## Where things are

| What | Where |
|---|---|
| The launch scripts | [experiments/](https://github.com/fesmc/kryos/tree/main/experiments) |
| Configs of the ensembles and their references | `tools/levante/plastic_ensemble_v1.toml`, `_v2.toml` |
| Par files | `par/kryos_NH-ANT-PAT_plastic_vilma2*.nml`, `par/kryos_NH-ANT_plastic_vilma2_b1.nml` |
| Designs | `docs/dev/plastic-ensemble-design.md`, `docs/dev/earth-search-design.md`, `docs/dev/solid-earth-literature.md` |
| Results (albedo) | `kryos-work/kryos_data/refcorr/<case>`, `kryos-plastic/kryos/runs/ens-v2`, `kryos-work/refs-test/ens` (scores) |
