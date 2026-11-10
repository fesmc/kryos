# Implied mass balance by joint inversion (design draft)

Status (2026-11-10): implemented and tested at present day, **paused**. The priority is the ice configurations and the
sea-level estimates; the inversion needs more development before it is effective (section "State and directions" at the
end). It replaces the forward closure of `mass_balance = "implied"` as the way to get the mass balance of a reconstructed
plastic ice sheet; the forward closure stays as the flow-law prior.

## Why

The forward closure (docs/dev/plastic-backend.md, "Implied mass balance") computes smb = div q from a flow law on the
plastic geometry. Every error of the flow law (the sliding is unknown) appears in smb, with either sign: at present day
the best fits leave ablation on 20-40 % of the Antarctic ice cells (RACMO: 0.4 %), and fits that match the observed
speed give several times the observed discharge. What we know about the climate, e.g. that Antarctica has no surface
melt today, is not used.

Ice shelves do not change the mass budget of the grounded ice: in a time step the flux across the grounding line is the
basin's SMB less its volume change. The grounding-line flux stays the sink cmb (calving and shelf basal melt together);
the shelves' buttressing sets how fast the ice crosses it, which we need not simulate because the geometry is imposed.

## Problem

Per solve (time t, geometry of the plastic solve: H, s, mask, tau_b; dH/dt over the step):

- unknowns: on each ice cell the depth-averaged speed U_i along the plastic flow direction -grad s/|grad s| (fixed by the
  geometry), and the mass balance a_i;
- hard constraint, mass conservation with the steady solver's face fluxes: div(H U n) = a - dH/dt (face flux: face
  thickness and face speed from its two cells; a margin face: outflow at the ice cell's speed, to land or to the ocean
  as now). For fixed geometry this is linear in U: D U = a - dH/dt, so a can be eliminated;
- objective (Gaussian, independent errors):
  J(U) = sum_i ((ln U_i - ln U_law,i) / sigma_u)^2 + sum_i ((a_i(U) - a_prior,i) / sigma_a,i)^2 (+ a smoothness term on
  ln U - ln U_law if needed);
  U_law the forward closure's speed at tau = tau_b (coulomb_n or yield, capped at tau_b), a_prior the prior SMB;
- bounds: U_i >= U_def,i (the deformation speed at tau_b: sliding >= 0); a_i >= 0 where the prior is cold (a_prior >= 0
  with no ablation term), free elsewhere.

The forward closure is the limit sigma_u -> 0; pure balance flux of the prior is sigma_a -> 0; a sign rule alone is
sigma_u -> 0 with only the bound.

## Solver

Gauss-Newton in x = ln U: each step a sparse SPD system (W_u + J_a' W_a J_a) dx = r (5-point divergence, so a ~13-point
normal matrix), solved with Lis CG; the bounds by an active set (projected Newton), as the steady solver's H >= 0. Size:
NH-32KM ~4e4 ice cells, ANT-16KM ~5e4, PAT-PAL-8KM ~2e4: seconds per solve. Warm start from the previous step's U.
Posterior: the inverse Hessian's diagonal by random probes (later), for sd maps of a and U.

## Prior SMB

From a surface instance of the domain (the plastic backend's optional `smb` import), so the prior is whatever surface
backend runs:
- `smb_simple` (snow line in latitude, maritime accumulation decaying inland, lapse-rate ablation on the synthetic or
  actual elevation) everywhere: one scheme for every domain and time; a rough first estimate, weak in accumulation;
- later, for the known domains at present day and as a pattern for the past: RACMO2.4p1 / MAR v3.14 with a lapse rate
  and a temperature-scaled accumulation (Clausius-Clapeyron), as a second surface backend.

sigma_a large to match: e.g. max(0.5 |a_prior|, 0.1 m ice/yr) for smb_simple.

## Placement

- PlasticIceSheet: `plastic_smb_inversion(ps, z_bed, z_sl, mask, rheo, slide, smb_prior, sigma_smb, dHdt, sigma_u, smb,
  speed [, cold, tau_cap, margin_frac, q_land, q_marine, stat])`, sharing the face geometry of the steady fluxes; test
  against the forward closure (sigma_u -> 0) and the balance flux (sigma_a -> 0).
- kryos plastic backend: `mass_balance = "inverted"`, keys `inv_sigma_u` (ln units), `inv_sigma_smb_rel`,
  `inv_sigma_smb_min`, `inv_cold` ("prior": a >= 0 where a_prior >= 0 | "none"); needs the smb import (a surface
  instance in the domain). Exports smb, cmb (= -rho_i q_marine), smb_tot, cmb_tot, ux_s, uy_s as now, plus the implied
  sliding U - U_def and the two misfit terms in the plastic stream.

## Checks

- Present day: Antarctica against RACMO (ablation area, net, discharge per basin) and Mouginot speed; Greenland against
  MAR and Joughin.
- 20 kyr: the NH and Antarctic totals, the ablation zones (southern Laurentide, Eurasian margins), the grounding-line
  discharge of the Antarctic shelf sectors.
- Sensitivity to sigma_u, sigma_a (the two limits above).

## Open questions

1. Speed misfit in ln U (Gauss-Newton) or linear in U (one SPD solve, quadratic program)?
2. Cold bound from the prior's sign, or from a temperature (surface temperature below a threshold)?
3. A smoothness term on the speed correction, or independent cells?

Decisions (2026-11-10): ln U with Gauss-Newton; the cold bound from the prior's sign; smoothness on the deviation of the
mass balance from the prior (`smooth_smb`), not on the speed correction (which made it worse, below).

## State and directions (2026-11-10)

### What is implemented

- PlasticIceSheet `src/plastic_inversion.f90` (`plastic_smb_inversion`, tests `tests/test_inversion.f90`, 9 checks):
  upwind face speeds (the mean of two cells decoupled odd and even cells along a flowline); the no-sliding bound soft
  (`sigma_def`, the uncertainty of the rate factor); Gauss-Newton in ln U with backtracking; each step by conjugate
  gradients applied matrix-free, preconditioned with the exact inverse of the mass-balance part (two triangular sweeps
  in descending-surface order, the upwind operator being lower triangular there) and, with `smooth_smb`, an inner CG for
  the correlated weight.
- kryos: `mass_balance = "inverted"` (the update step, with dH/dt over the step and the imported smb as the prior),
  `inv_*` keys, `u_inv`, `u_law` and the misfit terms in the plastic stream; `smb_simple` with `mask_ice` (the current
  ice) and `tsl_const` (no climate) as the prior; `par/kryos_NH-ANT-PAT_plastic_vilma2_inv.nml`; the offline check
  `tests/plastic_closure_fit.f90` with `prior_file`.
- Never run: a glacial run of `"inverted"` (the LGM run of the `_inv` par file was cancelled when the work was paused).

### Present-day Antarctica, RACMO2.4p1 as the prior (`runs/closure-cn/inv`)

Closure coulomb_n (A 3e-17, Hc 50 m, lambda_min 0.3, m 2, p 0.5, u_max 1000 m/yr), tau_b smoothed 50 km and capped,
sigma_u 0.7, sigma_a = max(0.5 |prior|, floor); Gt/yr, ablation area the fraction of the ice cells with a < 0:

| run | L_a (smooth_smb) | floor [kg m-2 yr-1] | net = discharge | ablation | ablation area | RMS log10 speed | wall time |
|-----|------------------|---------------------|-----------------|----------|---------------|-----------------|-----------|
| forward closure | - | - | 1905 | -106 | 25 % (100 km smoothed) | 0.557 | seconds |
| a0 | 0 | 100 | 1239 | -17.5 | 35 % | 0.482 | 6.5 min |
| a50 | 50 km | 100 | 1589 | -11.0 | 17 % | 0.490 | 67 min |
| a100 | 100 km | 100 | 1924 | -9.0 | 9.5 % | 0.511 | 84 min |
| a50f30 | 50 km | 30 | 1633 | -5.8 | 3.8 % | 0.499 | 61 min |
| RACMO2.4p1 | | | 2049 | -1 | 0.4 % | | |

- The inversion removes the forward closure's errors: positive SMB nearly everywhere, RACMO's coastal maximum, the
  speed closer to the observed (0.48-0.51 against 0.56). Smoothing the SMB deviation removes the streaks along converging
  flowlines and brings the total towards RACMO (1924 at 100 km); a tighter floor removes most of the ablation.
- The total stays below the prior: the slow flow-law speeds of the interior pull the flux down, and the misfit of a is
  cheap where sigma_a is large.
- 3000-7000 cells end below the deformation speed by more than sigma_def: there the geometry with A = 3e-17 deforms
  faster than the prior's accumulation can feed (A is uncertain in space; the soft bound absorbs it).
- Smoothing the speed correction instead (32, 100 km) made it worse (ablation -54, -157 Gt/yr): a smooth correction
  cannot give the grid-scale speed changes the flow routing on the plastic surface needs.

### Why it is paused

- Cost: with `smooth_smb` each solve takes about an hour on ANT-16KM (32 threads), almost all in the inner solve of the
  correlated weight inside the preconditioner. A glacial run solves at every step (800 steps over 80 kyr).
- The prior: `smb_simple` (latitude snow line, maritime accumulation, no time dependence beyond the elevation) is the
  only prior available for the past; the inversion is only as good as the prior where the flow law is weak.
- The 3-D-Earth ice configurations and their sea level are the priority (the weighted ensemble).

### Directions

1. **Cost of the correlated prior error.**
   - Replace the inner CG by a fast solver of the screened Laplacian (w + S_a): multigrid, or an FFT/DCT solve on the
     rectangular grid with w approximated by its mean (as a preconditioner of the inner solve or of the outer one).
   - Warm-start the inner solves from the previous preconditioner application.
   - Or parameterise the correlated error by a coarse basis (e.g. the deviation on a 4x coarser grid, interpolated):
     fewer unknowns, the correlation built in, no inner solve.
2. **Transient runs.** Solve at the output times (2.5 kyr slices) rather than every step, warm-started from the previous
   slice; use dH/dt over the slice interval.
3. **The flow-law prior.**
   - A spatially varying rate factor (a field of A with its own prior, or A from a temperature estimate) instead of the
     soft bound absorbing the misfit.
   - sigma_u per region or per flow regime; refit the closure (A, Hc, lambda_min) against the observed speed with the
     inversion in the loop rather than on its own.
4. **The SMB prior.**
   - A second surface backend: RACMO2.4p1 / MAR v3.14 at present day with a lapse rate, scaled to the past by a
     temperature anomaly (Clausius-Clapeyron for the accumulation, degree days for the ablation).
   - sigma_a from the spread of models, not a fixed fraction.
5. **The geometry.** The flow routing on the plastic surface carries the kinks of the tau_b domains: tau_b smoothing
   (`taub_smooth_km`) and, for the flow directions only, a smoothed surface.
6. **Checks still to do.** Greenland against MAR; the LGM run of `par/kryos_NH-ANT-PAT_plastic_vilma2_inv.nml` with
   the `smb_simple` prior (20 kyr totals, ablation zones of the southern Laurentide and Eurasian margins, the Antarctic
   grounding-line discharge); the sensitivity to sigma_u and sigma_a.
7. **Posterior uncertainty.** sd maps of a and U from the diagonal of the inverse Hessian (random probes).

