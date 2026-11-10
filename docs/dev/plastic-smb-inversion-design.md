# Implied mass balance by joint inversion (design draft)

Status: draft for discussion (2026-11-10). Replaces the forward closure of `mass_balance = "implied"` as the way to get the
mass balance of a reconstructed plastic ice sheet; the forward closure stays as the flow-law prior.

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
