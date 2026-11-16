---
title: "The ice: τ_b ensembles"
---

Our own ice: PlasticIceSheet online on VILMA2, its thickness computed from the margin chronology and the basal shear
stress τ_b. What we sample, how members are weighted, what the ensembles so far show, and the final ensemble that
follows the Earth search. Status 2026-11-16: ens-v1 and ens-v2 done on Bagge's structure 6; the final ensemble designed,
not run.

## What is sampled

τ_b on Gowan's domains (drawn by Gowan et al. 2021 from topography and surficial geology), as the product of three
factors ([how τ_b is composed](https://fesmc.github.io/kryos/backend-plastic.html)):

τ_b(d, t) = τ_Gowan(d, t) · exp(u_d) · exp(θ_g(t))

- τ_Gowan(d, t), Gowan's series per domain;
- u_d, the **present-day factors**, one per domain, constant in time, fitted offline to today's ice thickness of
  Greenland and Antarctica and then held fixed;
- θ_g(t), the sampled **log multipliers**, per scale group of domains (Laurentide, Greenland, Fennoscandia,
  Barents-Kara, British-Irish, West and East Antarctica, Patagonia) and time node (5-25 kyr apart), the present day
  pinned at 0. About 45 parameters.

A first approach, ES-MDA calibrations of the same multipliers, was retired: it assumes Gaussian errors and returns one
compromise, and its spread collapsed (sd 0.27 to 0.11 in calibration v2).

## The weighted ensemble

- **Members.** A space-filling sample (Sobol points) of the prior, each member a full coupled run of 80 kyr
  (about 240 core-h on a 3-D Earth).
- **Weights.** w ∝ L(θ) p(θ)/q(θ), with the likelihood of [Data and scores](data-scores.md) and q the sampling density.
  No member is discarded; the effective sample size and the weight at the edges of the range are the diagnostics.
- **The Earth as a dimension** (ens-v2). Shifts of log10 η of Bagge's structure 6, sampled with τ_b: one above and one
  below 670 km, and regional upper-mantle shifts (Laurentide, Greenland, Fennoscandia/Barents, West and East
  Antarctica; prior sd 0.5 dex).

## Results: ens-v1 and ens-v2

Both on Bagge's structure 6 (the Earth with the largest lateral contrast of Bagge's), before the Earth search showed it
to be a poor Earth ([Results](results.md#ensembles)).

- **The best members beat GLAC-1D overall, not the imposed references.** ens-v2's best, m048 (log L −9017) and m036
  (−9063, 6 of 7 groups), beat the published GLAC-1D (−9117) but not the published ice on our solver (−8860, −8835).
- **Fennoscandia limits them.** 4 of the 50 ens-v2 members pass Fennoscandia; on structure 6 it fails as it does under
  imposed ice. The Fennoscandian misfit is the Earth's: with PaleoMIST's ice, Gowan's 1-D Earth is 20-40 m high around
  the Gulf of Bothnia in the early Holocene, structure 6 is 20-85 m low, and every member follows structure 6.
- **The prior allowed ice no ice sheet can make.** m036 loses about 12×10⁶ Gt of Antarctic ice between 20 and 15 ka,
  the interior thinning by 1500 m in 5 kyr, because its East Antarctic multiplier drops from ×1.82 to ×0.71. Three
  constraints follow for the ensembles to come:
  1. a **random-walk prior** in time from the pinned present (step sd 0.15 per 5 kyr: sd 0.3 at 20 ka, 0.6 at 80 ka);
  2. a **bound on the rate of ice change** per scale group: at most 1.5 times the fastest of PaleoMIST, GLAC-1D, ICE-6G_C
     and ICE-5G over 2.5 kyr (it rules out 21 of the 50 ens-v2 members, 20 for East Antarctica);
  3. the **Antarctic excess volume** at 20 ka within 9.6-15.4 m sea-level equivalent (Blasco et al. 2021; ens-v2 spans
     −6 to +45 m, 9 of 50 inside).

## The final ensemble (proposed)

Designed 2026-11-13, to run on the Earths of the Earth search (its wave C). ens-v2 is not redone.

**Why a new design.** The ensembles so far multiply Gowan's τ_b, which Gowan tuned by hand against RSL with his 1-D
Earth. Centring the prior on his values biases our ice toward that Earth, which the Earth search shows is not the one the
data want. A wide prior is needed, but most coupled runs would then be spent on ice no ice sheet could make.

**The pieces:**

1. **A τ_b prior not tied to Gowan's values.** Per domain, between Gowan's minimum and base value, widened; Gowan's own
   history is one point inside it, not its centre. The present-day factors stay pinned.
2. **A screen on the ice alone.** The plastic ice costs about 1 min of a 4 h coupled run. Per Earth, thousands of τ_b
   histories run on a fixed bed and sea level and are ruled out by what needs no RSL: the rate of ice change, the
   Antarctic volume, ESL from the ice volume, today's ice.
3. **RSL predicted from step responses.** VILMA2's Earth is linear and constant in time, so the RSL of any load history
   is a convolution of the load with that Earth's response to a load step. About 24 step-load runs per Earth (3 load
   patterns per scale group) then predict the RSL and the GNSS rates of every screened history in milliseconds. The
   nonlinear parts (the moving ocean function, marine grounding, the bed feedback on the ice) are measured against
   coupled runs and enter as an emulator variance.
4. **History matching per Earth.** Score the predictions, rule out, resample what is left, repeat.
5. **Confirmation.** About 10 coupled runs per Earth over what is left: the prediction error, the bed feedback and the
   final scores. These are the members of the final ensemble.

On 4-8 Earths, the same prior and screen on each, so the ice-Earth trade-off shows directly: which histories fit on
which Earth. About 7000 core-h per 3-D Earth, against 12,000 for one coupled ensemble of 50 members, with 100 times more
histories judged.

**To decide:** the widening of Gowan's per-domain range; the number of load patterns, and whether the step responses
are per scale group or per region; whether ESL judges in the screen or only at the end; and the Earths, after wave B2.
