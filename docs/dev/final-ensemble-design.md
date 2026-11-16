# Our ice on sampled Earths: the final ensemble (design, 2026-11-13)

Status: proposed. Follows the Earth search (`earth-search-design.md`, wave C) and the constraints on the ice
(`plastic-ensemble-design.md`, "Constraints on the ice"). ens-v2 is not redone.

## Why

The ensembles so far multiply Gowan's tau_b. Gowan drew the domains from topography and surficial geology and tuned their
values by hand against RSL and GIA data with SELEN on his 1-D Earth: raised to maximise the LGM volume, lowered during the
deglaciation, and for Greenland and Antarctica first fitted to today's ice surface (`backend-plastic.md`). His values are a fit
made with one Earth. Centring the prior on them biases our ice towards that Earth, which the Earth search showed is not the
one the data want under either published ice history. A narrow prior around his values also limits what the ensemble can
find, while a wide one wastes most coupled runs on ice no ice sheet could make (ens-v2: m036, `plastic-ensemble-design.md`).

The aim: a wide prior on tau_b, made affordable by ruling out cheaply what can be ruled out without a coupled run, and by
predicting the RSL of a candidate ice history on a given Earth without running it.

## The pieces

**1. A tau_b prior not tied to Gowan's values.** Keep his domains (geology). Per domain, tau_b between his `taub_min` and
`taub_max` (FesmData `IceHistory-PaleoMIST_taub`: the base value and the minimum his deglacial adjustment goes to), widened by a
factor (to be set, e.g. 2 each side), sampled per scale group and node in log space with the random walk in time (step 0.15 per
5 kyr). Gowan's own history is one point inside this prior, not its centre. The present-day factors (Stage 1) stay pinned.

**2. Screen on the ice alone.** The plastic ice needs no Earth beyond a bed and a sea level: in a coupled run it costs about
1 min of 80 kyr out of 4 h (`timing.txt` of ens-v2: the ice 0.7 s per kyr, VILMA2 the rest). Per Earth, run thousands of
tau_b histories on a fixed bed and sea level (those of the Earth's reference run: the published ice on that Earth,
`gia_refcorr.sh`), and rule out with the checks that need no RSL:

- the rate of ice change per scale group (`[rate]`),
- the excess volume windows (`[[volume]]`, Antarctica at 20 ka: Blasco et al. 2021),
- ESL against the Spratt and Lisiecki (2016) stack (from the ice volume, no GIA needed to first order),
- today's Greenland and Antarctica (Stage 1).

What survives is the ice space the physics allows, without narrowing the prior.

**3. Predict RSL on each Earth from a few runs.** VILMA2's Earth is linear (Maxwell) and does not change in time, so the RSL
change at a site from a load is a convolution of the load history with that Earth's response to a load step. Per Earth and
scale group, K load patterns (the leading EOFs of the screened ice loads of the group, K about 3) each imposed as a step in
VILMA2 standalone (gia_refcorr.sh, one pass): about 8 groups × 3 = 24 runs per Earth. The RSL and GNSS rates of any screened
history then follow in milliseconds: project its loads on the patterns, convolve with the step responses, add up. The ocean
load is part of each run's response (sea-level equation with that run's ocean function).

What this misses, and where it is checked (piece 5): the ocean function and the grounding of marine ice change with the
history (the sea-level equation is then not linear), the plastic ice responds to the bed it sits on (the screen uses a fixed
bed), and loads outside the K patterns. The misfit of the prediction against coupled runs enters the likelihood as an
emulator variance, as in the Earth search.

**4. History matching in waves, per Earth.** Score the screened histories with the predicted RSL (the error model of the
references, the bar), cut, refocus the sampling on what is not ruled out, repeat. The range narrows only by the data.

**5. Confirm with coupled runs.** Per Earth, about 10 coupled runs spread over the histories not ruled out (and a few just
outside): the prediction error (piece 3), the bed feedback (if large: screen again on the bed of the median coupled run), and
the final scores. These are the members of the final ensemble.

**Earths.** A handful (5-8) of plausible Earths from wave B, spread over the space it leaves and including the closest to the
expert fields, and Gowan's for reference. The same tau_b screen and prior on all of them, so the ice-Earth trade-off shows
directly: which histories fit on which Earth.

## Cost (per Earth)

| piece | runs | cost |
|---|---|---|
| 2 screen | 5,000-10,000 plastic histories | about 1 min each (to be measured standalone): 80-170 core-h serial |
| 3 response | about 24 VILMA2 standalone, one pass | 1-D Earth: minutes each; 3-D: about 190 core-h each (wave B), 4,500 core-h |
| 4 waves | the predictions | negligible |
| 5 confirm | about 10 coupled runs | about 240 core-h each, 2,400 core-h |

On a 3-D Earth about 7,000 core-h, against about 12,000 for one coupled ensemble of 50 members (ens-v2), with a far wider
prior and 100 times more histories judged.

## Tools

| what | where |
|---|---|
| the tau_b prior from `taub_min`/`taub_max` | plastic_esmda.jl `params`, `prior_cov` (a per-domain range instead of a per-group multiplier of Gowan's values) |
| the ice alone on a prescribed bed and sea level, many histories per job | a kryos par without VILMA2 (geo from the Earth's reference run), batched like `gia_refcorr_pool.sh` |
| the screen | plastic_ensemble.jl: `[rate]`, `[[volume]]` (done), ESL and Stage 1 checks on the screen's output |
| load patterns and step responses | new: EOFs of the screened loads per group, step-load ice histories for gia_refcorr.sh |
| the RSL prediction and the waves | new, beside plastic_earth_emulate.jl: convolution, the likelihood of plastic_ensemble.jl, the cut |
| the confirm runs | plastic_ensemble.sh as now |

## To decide

1. The widening of Gowan's per-domain range (piece 1).
2. K, and whether the step responses are per scale group or per Gowan region (fewer runs, coarser).
3. Whether the ESL stack judges in the screen or only in the final scores (the bar for ESL: see `earth-search-design.md`).
4. The Earths: how many, and which.
