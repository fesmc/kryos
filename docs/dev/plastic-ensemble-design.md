# Weighted ensemble of plastic ice histories on the 3-D Earth (design draft, 2026-11-10)

Status: draft for discussion, nothing implemented. Follows calibration v2 (ES-MDA, docs/dev/plastic-backend.md) and the
residual and spread diagnostics of 2026-11-09/10.

## Goal

A range of ice histories (80-0 kyr) and their RSL, each consistent with GAPSLIP 2.1 within an honest error model, on the
3-D Earth (Bagge 2021, science setting): weights instead of rejection, so that low, mid and high ice volumes and
regional trade-offs ("thicker here, thinner there") come out with their plausibility. ES-MDA is not used: it assumes
Gaussian errors and returns one compromise (its spread collapsed from sd 0.27 to 0.11 in v2).

## Basis (fixed in every member)

- Earth: VILMA2 3-D science setting (Bagge 2021 structure 6, lmax 170, refcorr reference PaleoMIST-sci). Earth
  permutations come later as an extra dimension of the same ensemble (below).
- Instances: NH-32KM, ANT-16KM, PAT-PAL-8KM (`par/kryos_NH-ANT-PAT_plastic_vilma2.nml`), PaleoMIST margins, Gowan's tau_b
  domains times the present-day factors (Stage 1, `runs/taub-fit-v3`).
- Output: the product (2-D ice every 2.5 kyr, the VILMA2 stream), as plastic_esmda.jl now writes.

## Parameters

Log multipliers of tau_b per scale group and time node, the present day pinned (0):

| group | nodes [kyr] | change from v2 |
|---|---|---|
| Laurentide | 80, 60, 40, 25, 15, 5 | fewer deglacial nodes (zigzag), prior correlation 10 kyr |
| Greenland | 80=60, 40, 25, 20, 15, 10, 5 | - |
| Fennoscandia | 80=60=40, 25, 20, 15, 10 | split from Eurasia |
| Barents-Kara | 80=60=40, 25, 20, 15, 10 | split from Eurasia |
| British-Irish | 25, 20, 15 | split from Eurasia |
| West Antarctica | 80=60=40=25=20, 15, 10, 5 | split; prior only beyond the Holocene data |
| East Antarctica | 80=60=40=25=20, 15, 10, 5 | split |
| Patagonia | 25=20, 15 | new; no RSL data: prior only (its load reaches the far field) |

About 45 parameters. The splits give the regional freedom the v2 spread showed missing (one multiplier per group scales
a whole sheet). Prior: N(0, 0.3^2) per node (factor 0.74-1.35 at 1 sd), exponential correlation in time within a group,
independent groups.

## Likelihood

d = model - obs at the index points (pm_compare.jl), R = diag(sigma_i^2) + sigma_m,g^2 C, C_ij = exp(-dist / L_x,g)
exp(-|dt| / L_t,g), per data group g, from the residual diagnostics (plastic_residual_corr.jl, v2 posterior on the 1-D Earth;
to be redone on the first 3-D wave):

| group | sigma_m [m] | L_x [km] | L_t [kyr] |
|---|---|---|---|
| Laurentide | 12 | 350 | 10 |
| Greenland | 30 | 500 | 5 |
| Fennoscandia/Barents | 10 | 250 | 5 |
| Antarctica | 8 | 500 | 10 |
| far field | 13 | 50 | 5 |

(v2: 5 m, 150 km, 2 kyr for all.) Limiting data: one-sided terms (probability of the model RSL being on the right side,
with sigma_i and sigma_m), optional, off in the first wave. log L = -1/2 d' R^-1 d (+ the limiting terms).

## Sampling and weights

- Wave 1: N = 100 members, Sobol points in the standardized prior space (inverse normal), so the prior is covered evenly
  to 2.5 sd.
- Weights w_m ∝ L(theta_m) p(theta_m) / q(theta_m), q the sampling density of the member's wave (mixture over the waves:
  deterministic-mixture importance sampling). No member is ever discarded; a member with a tiny weight stays in the set.
- Diagnostics per wave: effective sample size ESS = (sum w)^2 / sum w^2 (overall and per data group), the weight of
  the members on the edge of the sampled range, posterior means and sd of the parameters.
- Refinement: when ESS < 20 or weight sits at the edge of the range, the next wave samples a proposal fitted to the
  weighted members (Gaussian, covariance inflated 1.5^2, extending into the low-weight region), N = 50-100. The error
  model is not tempered: if the weights collapse, the error model is revisited (a collapse with the R above would mean
  the model misses something systematic).

## Products

Weighted mean and quantiles (5, 50, 95 %) of H_ice, z_srf, RSL at every 2.5 kyr; ice volume per group and ESL
through time; scenario members picked from the high-weight set (low, mid, high LGM ESL; the largest regional
trade-offs, e.g. Laurentide against Antarctica); posterior parameter correlations (which groups trade off through the
far field); a held-out check (10 % of the sites by region, scored by the weighted ensemble).

## Earth as a dimension

Wave 1 (ens-v1) showed that the misfit in Fennoscandia/Barents is the Earth's: with PaleoMIST's ice, the 1-D Earth (B1)
is 20-40 m high around the Gulf of Bothnia in the early Holocene, the 3-D Earth (B3) 20-85 m low, and every member follows
B3 (plastic_ensemble_misfit.jl). PaleoMIST itself misfits west Greenland. So the Earth is sampled with tau_b (v2,
`plastic_ensemble_v2.toml`), and all data groups judge plausibility:

- Parameters (prior N(0, 0.5^2) dex each, independent, |shift| <= 1.5): shifts of log10(eta) of the Bagge field, a global
  one above 670 km (the upper mantle with the lithosphere) and one below (the lower mantle), and one per region above 670
  km: Laurentide, Greenland, Fennoscandia/Barents, West and East Antarctica. Above the lid threshold (log10 eta > 22) the
  field is elastic, so an upper-mantle shift also moves the lithosphere's base.
- Regions: masks from the tau_b domains of their scale groups (300 km buffer, 500 km cos^2 taper, shared where they
  overlap), `kryos_data/plastic/earth_regions_v2.nc` (`plastic_ensemble.jl earth-regions`). Iceland's domains belong to
  the Laurentide tau_b group and Iceland lies within Greenland's buffer, so no regional shift applies there (`exclude_codes`,
  REGIONS region_2 105, tapered over 200 km).
- Each member: a shift file `iter_kk/earth/mNNN.nc` (VILMA2 `visc_shift_file`, `plastic_earth_shift.jl`).
- The present-day reference correction r (refcorr) stays that of the unperturbed Earth; its error is measured on
  PaleoMIST-sci with all shifts at +0.5 and -0.5 dex against a rerun without shift (`gia_refcorr.sh`, 2026-11-11).

## Cost

A coupled 3-D run took 3.5-4.2 h on 64 cores (about 240 core-h): N = 100 is about 24,000 core-h per wave, 30 at a time
(about 13 h of wall time per wave). Before the first wave: measure the 3-D scaling (16 and 32 cores per member may be
cheaper per run), and whether a member can start from a shared spun-up 80 kyr state.

## Tools

`plastic_ensemble.jl` (sample, write members as plastic_esmda.jl does, score, weights, refine), `plastic_ensemble_report.jl`
(the products above, from the existing map and comparison tools), config `plastic_ensemble_v1.toml`.

## Open questions

1. The group splits and nodes above.
2. sigma_m and correlation lengths: from the 1-D v2 residuals now, re-estimated after wave 1?
3. N per wave and the core-hours.
4. Limiting data in the likelihood from the start, or later?
