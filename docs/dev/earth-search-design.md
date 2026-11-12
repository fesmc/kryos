# A broad search of the solid Earth with emulation (design, 2026-11-12)

Status: for discussion, nothing implemented. Follows the literature review (`solid-earth-literature.md`), the Earth basis from
expert fields (`plastic_ensemble.jl earth-basis`) and the first scores against the published reconstructions
(`docs/simulations.md`).

## Why

The expert 3-D viscosity fields (Bagge et al. 2021, Lloyd et al. 2024, Pan et al. 2022) disagree by about 1 dex in the lower
mantle and in the lithosphere, and no study ranks them against RSL data. Inversions in the literature constrain a few regions
(Fennoscandia, Hudson Bay) and leave the rest open. Searching only around one expert field (ens-v2: shifts of Bagge's structure 6)
is too narrow: PaleoMIST's ice fits all groups on Gowan's 1-D Earth (`ice-PMIST | E1D-PMIST`) and on none of the 3-D Earths tried.

So: search broadly for the Earth patterns the data want, with an emulator to make it affordable, and report how far each one is
from the expert fields, instead of building that distance into the score.

## Parameters

The field is a radial profile, plus the lateral pattern of the experts scaled, plus regional shifts:

log10 η(x, z) = P(z; θ_r) + α A(x, z) + Σ_k c_k φ_k(x, z) + Σ_j s_j w_j(x) [z < 670 km]

| part | parameters | range | note |
|---|---|---|---|
| radial profile P | lithosphere (elastic) thickness | 50-200 km | `E1D-PMIST` 120 km |
| | upper mantle (lithosphere to 410 km) | 10^19.5-10^21.5 Pa s | `E1D-PMIST` 4×10^20 |
| | transition zone (410-670 km) | 10^20-10^22.5 | |
| | lower mantle, upper (670-1500 km) | 10^21-10^23.5 | `E1D-PMIST` 4×10^22 |
| | lower mantle, lower (1500-2891 km) | 10^21-10^23.5 | |
| lateral | α, amplitude of the experts' mean lateral anomaly A (basis v2 mean minus its radial mean) | 0-1.5 | α = 0: 1-D |
| | c_2, c_3: the basis v2 components (Lloyd against Pan; Bagge's reduction factor) | the experts' range | c_1 is mostly radial: in P |
| regional | s_j, upper-mantle shifts in Fennoscandia, Barents, Laurentide, Greenland, West Antarctica (masks of ens-v2) | ±1 dex | |

13 parameters, of which the first 5 make a 1-D Earth. Fields are written as VILMA2 viscosity files (`earth-field`); a 1-D
member is a 1-D par file (`earth-radial` style), which runs in about 3 min instead of 1-2 h.

## Waves

1. **Wave A, 1-D (5 radial parameters).** About 200 Earths, Latin hypercube, each under two imposed ice histories (PaleoMIST
   and GLAC-1D 80 ka), with its own present-day correction (gia_refcorr iterates in 4 passes of about 3 min). About 400 runs,
   1-2 days on a few nodes.
2. **Emulator and cut.** Per data group, a Gaussian process of the score (log of χ²/n, or I_g) in the standardized parameters
   (Matérn 5/2, fitted by maximum likelihood), with its own variance as part of the implausibility. Implausible: worse than the
   bar (the published reconstructions) in a group that judges, for both ice histories. What remains is the 1-D Earth space the
   data allow.
3. **Wave B, 3-D.** About 40-60 Earths in the remaining radial space with the lateral and regional parameters (13 in all), under
   both ice histories, one pass each with the present-day correction of their radial profile (R0). About 100 runs of 1-2 h.
   Emulator and cut as in 2.
4. **Wave C, our ice.** A τ_b ensemble (as ens-v2) on a handful of Earths from wave B: the plausible ones, spread over the
   remaining space.

## Scores

The data groups of the v2 config: the RSL groups (Fennoscandia and Barents split), GNSS Fennoscandia and Hudson Bay (present-day
rates: the most direct constraint on the upper mantle there), ESL (Spratt and Lisiecki 2016: mainly the ice, little the Earth).
The bar: as good as the best published reconstruction in every group that judges.

## Agreement with the experts, a second axis

For each Earth: its distance from the expert fields, in two parts:

- radial: the RMS difference of log10 η of its profile from the nearest expert's radial mean (and from Gowan's 1-D), over depth
  weighted by volume;
- lateral: the distance of (α, c_2, c_3) from the experts' coefficients (α = 1, c from `coef` of the basis).

The result is a front: the fit to the data against the distance from the experts. An Earth far from the experts that fits much
better is a finding; one that fits only a little better is not worth it. This stays a choice for us to read, not a weight in the
score.

## The ice-Earth trade-off

An Earth search with fixed ice blames the ice's errors on the Earth. Two ice histories (PaleoMIST, GLAC-1D) that differ where
they differ: an Earth pattern counts when it improves both. Our own ice comes in at wave C, on Earths that are already plausible
for both.

## To decide

1. The ranges and the 13 parameters above.
2. The emulator: our own small Gaussian process (LinearAlgebra only, about 100 lines) or a package (AbstractGPs.jl: a new
   dependency of tools/levante).
3. The size of waves A and B.
4. Whether GNSS judges from wave A on, or is reported only until its error model is fitted.
