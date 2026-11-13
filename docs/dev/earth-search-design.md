# A broad search of the solid Earth with emulation (design, 2026-11-12)

Status: decided, being implemented. Follows the literature review (`solid-earth-literature.md`), the Earth basis from
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

## Decided (2026-11-12)

1. The 13 parameters and ranges above.
2. The emulator: AbstractGPs.jl (with KernelFunctions.jl), new dependencies of tools/levante; hyperparameters by maximum
   likelihood on a grid (no optimizer package).
3. Wave A: 200 1-D Earths, wave B: 50 3-D Earths, each under both ice histories.
4. GNSS judges from wave A on (the runs need vilma2 with `dz_bed_rate`: they run on a stack built on vilma2 main).

## Implementation and wave A (2026-11-12/13)

**Tools** (tools/levante):

- `plastic_ensemble.jl earth-search-sample <1-D par> <outdir> [n] [seed]`: wave A, a Latin hypercube of the 5 radial parameters,
  uniform over the ranges above, as 1-D par files (`plastic_earth_basis.jl` `earth1d_par`); `earth-search-par` writes one.
  VILMA2 takes density and rigidity per layer as PREM's shell average (`earth = "PREM"`), so a 1-D Earth keeps the layering of its
  template (Gowan's 12 layers of `par/kryos_NH-ANT_plastic_vilma2_b1.nml`: the PREM discontinuities and the lower mantle in parts)
  and adds the boundaries of the viscosity (the lithosphere base, 410, 670 and 1500 km); a template boundary within 5 km of one of
  those is dropped. Elastic (rheology 0, η = 10⁴⁰) above the lithosphere base, Maxwell below, the fluid core below 2891 km.
- `gia_refcorr.sh`: `NCPU`, `MEM` (a 1-D Earth needs a fraction of a node), `KEEP=last` (only the converged iteration's
  `vilma_out.nc`), `ZREF`/`XZREF` (today's bed from PaleoMIST for GLAC1D_80k, which has no `z_bed`), `provenance.txt` per case.
- `gia_refcorr_pool.sh`: many small cases in few jobs, a pool of workers per job (each on its own cores), each claiming the next
  case of a list (an atomic link to the job id).
- `plastic_ensemble.jl errmodel ... refs=<regex> groups=<g;...>`: the error model of the point groups (GNSS, ESL) from run
  references, here the wave-A cases (no ensemble member has `dz_bed_rate`).
- `plastic_earth_emulate.jl`: the emulators and the cut. Emulated per ice history and group: y = log(χ²/n), not I_g. Both are
  monotone in χ², so the cut is the same (y_cut = log(1 + cut √(2/n))), but I_g spans orders of magnitude over the space (GNSS
  Hudson Bay from 3 to 3000) and a stationary GP of it follows the worst Earths. Cut: I_emu = (E[y] − y_cut)/sd_emu > 3 in a group
  that judges rules a point out.

**Check of the parameterization.** Gowan's Earth in the parameters (120 km; 4×10²⁰, 4×10²⁰, 4×10²², 4×10²² Pa s) against the b1
par: the same η and rheology at every km of depth, all 12 boundaries of b1 kept, two added (410 and 1500 km, the same η on both
sides), the rest of the par identical. Run as `ice-PMIST | E1D-ES1-PMIST`, it repeats `ice-PMIST | E1D-PMIST` (same stack):
the bed misfit per iteration 269.93/37.37/7.16 m against 269.74/37.33/7.15 m, log L −8891.5 against −8891.4 (7 groups), every
I_g within 0.1.

**Runs.** 200 Earths × 2 ice histories, 402 cases with the two of Gowan's Earth, on the kryos-gnss stack (vilma2 main, with
`dz_bed_rate`), in pools of 8-16 cores per case. A case converges in 2-5 passes (most in 3-4); a pass is 4-5 min (PaleoMIST,
34 records) to 7-8 min (GLAC-1D, 351 records) on 8-16 cores. 8 cores per case cost half the core-hours of 16 (2200 against 4400
core-s per PaleoMIST pass); 32 cores cost 7000. In all about 12,000 core-hours, a third of it lost to a bug of the pool's first
claim (`ln -s` into an existing claim folder succeeds), which ran 141 cases twice at once; those were deleted and rerun. 400 of
402 cases scored (two converged late, scored with wave B). The GLAC-1D output is 2.6 GB per case (all 351 records).

**Scores.** As the references of `docs/simulations.md`: the window ≤ 25 ka, the error model of refs-test
(`kryos-work/refs-test/ens/errmodel.toml`, Fennoscandia and Barents refitted after the split), the GNSS groups refitted on the
wave-A runs (top 6 per group: GNSS Fennoscandia σ_m 0.7 mm/yr, L_x 750 km; GNSS Hudson Bay 1.5 mm/yr, 500 km; the starting
value 1 mm/yr ruled out every Earth under PaleoMIST), ESL at its starting values. Groups that judge: the 6 RSL groups, GNSS
Fennoscandia and Hudson Bay, ESL. Cut: Greenland 7.88, ESL 7.04 (PaleoMIST), every other group 3. On this scale I1 =
`ice-PMIST | E1D-PMIST` −8860.3 and I5 = `ice-GLAC1D | E1D-PMIST` −8835.5 (log L over the RSL groups).

**Emulator.** Leave-one-out R² 0.96-0.998 for the RSL groups except the far field (0.78-0.89, a narrow range of y and a
nugget), 0.85-0.98 for GNSS, 0.97-0.98 for ESL; 91-98 % of the cases within 2 LOO sd. The LOO prediction never rules out a case
that passes the cut (2 of 400 for GNSS Fennoscandia under PaleoMIST); it keeps many that fail, which is the conservative side.

**Not ruled out** (100,000 points; 5-95 % ranges):

| set | fraction | without GNSS | lithosphere [km] | UM | TZ | LM1 | LM2 |
|---|---|---|---|---|---|---|---|
| PaleoMIST | 20.2 % | 30.3 % | 54-167 | 19.7-21.2 | 20.1-22.4 | 21.1-23.4 | 21.1-23.3 |
| GLAC-1D | 13.5 % | 26.3 % | 54-145 | 20.4-21.4 | 20.1-22.2 | 21.1-23.4 | 21.1-23.4 |
| both | 4.9 % | 19.8 % | 54-140 | 20.3-21.4 | 20.0-22.3 | 21.1-23.4 | 21.1-23.4 |
| either | 28.8 % | 36.8 % | 54-162 | 19.8-21.3 | 20.1-22.4 | 21.1-23.4 | 21.1-23.4 |

The data constrain the upper mantle and the lithosphere (< 170 km) and leave the transition zone and the lower mantle open in 1-D;
the pairs show trade-offs (a stiffer UM with a soft TZ). By group, Fennoscandia rules out half the space (a thin lithosphere and a
soft UM and TZ fit it), Greenland a quarter (it wants a stiffer UM), GNSS Fennoscandia 38-62 %, GNSS Hudson Bay 23-36 %; the far
field, Antarctica and ESL rule out nothing. The two ice histories want different Earths: PaleoMIST a soft UM (20.1-20.4) with a
stiffer TZ, GLAC-1D a UM of 20.6-21. Requiring both would rule Earths out for the ice, not the Earth; wave B samples the radial
part from "either".

**Best runs** (log L over the RSL groups; parameters: lithosphere, log10 η UM/TZ/LM1/LM2):

| ice | Earth | parameters | log L | against the reference Earth | groups passed (of 9) |
|---|---|---|---|---|---|
| PaleoMIST | E117 | 111 km, 20.13/21.75/23.15/21.69 | −8800.7 | +59.6 on I1 | 7 (fails GNSS F 5.7, HB 6.4) |
| PaleoMIST | E185 | 97 km, 20.21/21.80/22.11/22.94 | −8802.5 | +57.8 | 7 |
| GLAC-1D | E165 | 67 km, 20.77/20.08/22.54/21.88 | −8708.9 | +126.6 on I5 | 7 (Barents 3.3, GNSS F 3.6) |
| GLAC-1D | E017 | 87 km, 20.56/20.61/21.68/21.02 | −8712.3 | +123.2 | 8 (GNSS F 6.7) |
| GLAC-1D, all 9 groups | E033 | 106 km, 20.84/21.43/21.08/22.08 | −8745.6 | +89.9 | 9 (14 Earths pass all 9) |
| PaleoMIST, most groups | E017 | 87 km, 20.56/20.61/21.68/21.02 | −8905.5 | −45.2 | 8 (fails GNSS F 4.9; none passes all 9) |
| both | E120 | 98 km, 21.04/20.64/22.02/22.27 | −8784.7 (GLAC-1D), −8913.2 (PaleoMIST) | +51, −53 | 9, 8 |

Gowan's Earth is not among the best under either ice: 1-D Earths gain 60 (PaleoMIST) and 127 (GLAC-1D) in log L, mostly
in Greenland and, under GLAC-1D, the far field. No 1-D Earth fits the GNSS rates of Fennoscandia and Hudson Bay under PaleoMIST's
ice. Files (Albedo, `/albedo/work/user/alrobi001/kryos-search`): `earths/` (members.csv, the par files), `kryos_data/refcorr/`
(the cases), `score/ens/` (config, errmodel, `iter_01/scores.csv`), `emulate/` and `emulate-nognss/` (summary.md, gp.csv,
loo.csv, sample.csv, the figures nroi_*.png).
