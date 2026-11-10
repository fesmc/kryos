---
title: "Implied mass balance (in development)"
---

::: {.callout-warning}
## In development, paused

This page describes code that exists and has been tested at present day, but not in a glacial run. Work on it is paused
(2026-11-10): the state and the possible directions are under [Status and open issues](#status-and-open-issues) and in
[docs/dev/plastic-smb-inversion-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-smb-inversion-design.md).
The keys, defaults and par file below may change.
:::

The plastic ice sheet imposes its geometry: the margin chronology and tau_b set the thickness, and nothing is
integrated (see [PlasticIceSheet](backend-plastic.md)). By default, the backend exports `smb = cmb = bmb = 0`, and
the mass change that the margin and the bed impose is `mb_net`. The *implied* mass balance answers a different
question: what surface mass balance and discharge does this geometry need, under a flow law? That gives a climate
consistency check of a reconstruction, and the mass fluxes for the other components.

The key `mass_balance` of `[ice_plastic]`
([input/kryos_defaults.nml](https://github.com/fesmc/kryos/blob/main/input/kryos_defaults.nml)) selects how it is
computed:

| `mass_balance` | What | Needs |
|---|---|---|
| `"none"` (default) | smb = cmb = 0 | |
| `"implied"` | forward flow closure: smb and cmb from the divergence of the steady face fluxes on the geometry | `rheo_A > 0`, a `slide_law` |
| `"inverted"` | joint inversion of speed and mass balance, with the closure's speed and an imported smb as priors | as `"implied"`, plus `drive = "taub"`, `mb_smooth_km = 0`, and an `smb` import (a surface instance in the domain) |

With either `"implied"` or `"inverted"`, the backend exports `smb`, `cmb` and their sums `smb_tot` and `cmb_tot`
[Gt yr-1], plus the surface velocity `ux_s` and `uy_s`. Its stream `plastic` adds `div_q`, `q_land`, `q_marine`,
`abl_margin` and `discharge`.

## Forward closure: `mass_balance = "implied"`

The flux is the membrane-free DIVA flux of the library's steady SMB-driven solver, on staggered faces:

$$q = -H\left(u_b(\tau) + \frac{2A}{n+2}\,\tau^n H\right)\frac{\nabla s}{|\nabla s|},\qquad \tau = \rho_i g H |\nabla s|$$

The flux from the last ice cell to the margin is one-sided (`plastic_smb_flux_divergence`). After each solve, the
backend computes div q per ice cell and the outflow through its margin faces. That outflow is `q_land` where the
ice-free neighbour's bed is above sea level, and `q_marine` where it is below. The exports, as means over the step
[kg m-2 yr-1], are:

- **On the ice:** smb = rho_i (dHidt + div q - q_land). The outflow to land is ablated in the margin cell, a sub-grid
  ablation zone. cmb = -rho_i q_marine (calving).
- **On the cells the ice left over the step:** their loss rho_i dHidt, as cmb where the bed lies below sea level,
  else as smb.
- **Budget:** smb + cmb summed over the grid is `mb_net`.

The closure has three controls:

- **Stress cap:** with `drive = "taub"`, the driving stress of the face fluxes is capped at tau_b, the mean of its two
  cells on a face (or the ice cell's value on a margin face). Without the cap, the steep step from the ice to the
  margin surface drives the outflow of the margin cells far above the stress of the ice. In the first `coulomb_n`
  fit, the marine margin cells had 10 times the implied SMB of the interior.
- **tau_b smoothing:** `taub_smooth_km` (default 50 km) applies a Gaussian mean of tau_b over the cells with
  tau_b > 0. It is used for the solve and the closure alike, so the domain steps become ramps.
- **SMB smoothing:** `mb_smooth_km` (default 0) applies a conservative Gaussian smoothing of div q - q_land over the
  ice. The sums are unchanged. The plastic surface is piecewise (steps at the domain boundaries, kinks along the
  flowlines), and its grid-scale flux divergence is mostly that.

### Sliding laws

A plastic sheet sits at a driving stress of tau_b, so its deformation flux falls from the divide outwards. A steady
flowline needs u = q/H to rise 10-20 times from the divide to the equilibrium line, so the sliding has to rise towards
the margin. Two laws use tau_b:

| `slide_law` | Law | Keys |
|---|---|---|
| `"yield"` | u_b = u_y (tau / tau_b)^m, u_y = `slide_uy` (tau_b / 100 kPa)^`slide_uy_b`; `slide_uy_b` < 0 makes a weak bed fast | `slide_uy`, `slide_uy_b`, `slide_m` |
| `"coulomb_n"` | regularised Coulomb, tau = C (u_b / (u_b + u_0))^(1/m), C = tau_b lambda N / (rho_i g `slide_Hc`) | `slide_Hc`, `slide_lambda_min`, `slide_zb0`, `slide_zb1`, `slide_p`, `slide_u0`, `slide_m`, `slide_umax` |

The terms of `"coulomb_n"`:

- **Effective pressure:** N = rho_i g H (1 - H_f/H)^`slide_p` (Leguy et al. 2014), with H_f the flotation thickness.
- **Bed weakening:** lambda runs linearly from `slide_lambda_min` at a bed of `slide_zb0` (relative to sea level) or
  lower, to 1 at `slide_zb1`.
- **Speed at tau = tau_b:** the sliding speed is u_0 x / (1 - x), with x = (`slide_Hc` rho_i g / (lambda N))^m, capped
  at `slide_umax`.

So with `"coulomb_n"` the tau_b domains enter through the geometry, not the speed. A uniform Weertman beta does the
opposite of what is needed: in the first test, the LGM Antarctic discharge was 200 Gt/yr against 5000 Gt/yr today.

### Known limitation

The forward closure computes smb = div q from a flow law whose sliding is unknown. Every error of the flow law
therefore appears in smb, with either sign. The design note records two symptoms at present day:

- the best fits leave ablation on 20-40 % of the Antarctic ice cells (RACMO: 0.4 %);
- fits that match the observed speed give several times the observed discharge.

The closure also does not use what is known about the climate, for example that Antarctica has no surface melt today.

In short, the forward closure puts the discharge in the wrong places: with the stress capped at tau_b, the yield
law's Antarctic discharge falls to 86 Gt/yr (RACMO: 2049), and `coulomb_n` keeps the total (1905) only by fast sliding
near flotation, so the flux rises within a smoothing length of the Ross and Ronne grounding lines while the interior
flux stays too small (present-day tests, `runs/closure-cn`: section "Implied mass balance: present-day tests" of
[docs/dev/plastic-backend.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-backend.md)).

This is why `"inverted"` replaces it as the way to get the mass balance. The forward closure stays as the flow-law
prior.

## Joint inversion: `mass_balance = "inverted"`

The method is documented in three places:

- **Design:** [docs/dev/plastic-smb-inversion-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-smb-inversion-design.md).
- **Library:** `plastic_smb_inversion` in PlasticIceSheet,
  [src/plastic_inversion.f90](https://github.com/fesmc/PlasticIceSheet/blob/main/src/plastic_inversion.f90), with its
  API in PlasticIceSheet `docs/fortran-port.md`.
- **Backend:** [src/backends/ice/kryos_ice_plastic.f90](https://github.com/fesmc/kryos/blob/main/src/backends/ice/kryos_ice_plastic.f90).

At each update, the inversion works on the geometry at t1:

- **Unknowns:** the depth-averaged speed U of each ice cell along the plastic flow direction -grad s/|grad s|, as
  x = ln U. The flow direction is fixed by the geometry.
- **Mass conservation:** a = dH/dt + div(H U n), with the face geometry of the steady fluxes. A face takes the speed
  of its upstream cell. The outflow to land is ablated in its margin cell and left out of a; the outflow to the ocean
  is in a and in `q_marine`.
- **Objective:** (1) ln U close to the closure's speed U_law at tau = min(rho_i g H |grad s|, tau_b), with sd
  `inv_sigma_u`; (2) a smoothness of that correction over `inv_smooth_km`; (3) a close to the prior a_prior, with sd
  sigma_a = max(`inv_sigma_smb_rel` |a_prior|, `inv_sigma_smb_min`); (4) a smoothness of the deviation from the prior
  over `inv_smooth_smb_km`, which makes the prior error spatially correlated; (5) a penalty on a < 0 on the "cold"
  cells (no margin face, a_prior >= 0), with scale `inv_sigma_cold` (`inv_cold`); (6) a soft bound keeping U at least
  the deformation speed (sliding >= 0).
- **Solver:** Gauss-Newton on x, at most `inv_n_iter` steps. Each step solves the normal equations matrix-free by
  conjugate gradients, preconditioned with the inverse of the mass-balance part: two triangular sweeps in order of
  descending surface. Each step ends with a backtracking line search on J.
- **Limits:** `inv_sigma_u -> 0` gives the flow law's speed (the forward closure). sigma_a -> 0 gives the balance flux
  of the prior.

On the ice, smb = rho_i a and cmb = -rho_i q_marine, as for `"implied"`, but without `mb_smooth_km`. `ux_s` and `uy_s`
are the speed plus u_def/(n+1) along the flow direction. The stream `plastic` adds `u_inv` and `u_law` (the inverted
speed and the flow law's speed, the inversion's prior), and the 0-D `inv_iters`, `inv_J_law`, `inv_J_smb`,
`inv_J_cold` and `inv_n_bound`.

### Keys

| Key | Default | Meaning |
|---|---|---|
| `inv_sigma_u` | 0.7 | sd of ln(speed) about the closure's (0.7: a factor 2) |
| `inv_smooth_km` | 0 | [km] length of the smoothness of ln(speed / closure's); 0 = none |
| `inv_smooth_smb_km` | 50 | [km] length of the smoothness of the deviation of a from the prior; 0 = none |
| `inv_sigma_smb_rel` | 0.5 | sd of a about the prior, relative to \|prior\| |
| `inv_sigma_smb_min` | 100 | [kg m-2 yr-1] its floor |
| `inv_cold` | .true. | penalise a < 0 where the prior is >= 0, off the margin |
| `inv_sigma_cold` | 10 | [kg m-2 yr-1] the scale of that penalty |
| `inv_n_iter` | 50 | maximum Gauss-Newton steps |

### The prior: smb_simple

The prior is the smb that the domain's surface instance exports. In the par file below, it comes from
[smb_simple](backend-smb-simple.md) with its defaults, which are not tuned (`&smb_simple_plastic`), and two settings:

- `mask_ice = True`: its mask is the current ice (the hub's `H_ice`, at every call).
- `tsl_const = 263.15`: no `tsl_ann` import, for a domain without a climate. Only `T_srf` uses it, and the plastic
  ice does not read `T_srf`.

The design note calls this "a rough first estimate, weak in accumulation", so its sd is large. RACMO2.4p1 or MAR v3.14
with a lapse rate is proposed as a later second surface backend.

### Par file

[par/kryos_NH-ANT-PAT_plastic_vilma2_inv.nml](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT-PAT_plastic_vilma2_inv.nml)
is `par/kryos_NH-ANT-PAT_plastic_vilma2.nml` with the following changes:

- **Prior:** a surface instance `smb_simple` per domain (`surface_nor`, `surface_ant`, `surface_pat`), as above.
- **Inversion:** `mass_balance = "inverted"`, `taub_smooth_km = 50`, `mb_smooth_km = 0` and the `inv_*` defaults,
  written out per instance. The stream `plastic` is written on `out_2D`, as in the base file.
- **Flow law:** the `coulomb_n` closure of the present-day fits (`runs/closure-cn`, the stress capped at tau_b). Its
  values per instance:

| | `ice_nor` (NH-32KM) | `ice_ant` (ANT-16KM), `ice_pat` (PAT-PAL-8KM) |
|---|---|---|
| `rheo_A` [Pa-3 yr-1], `rheo_n` | 1.0e-16, 3 | 3.0e-17, 3 |
| `slide_Hc` [m] | 100 | 50 |
| `slide_lambda_min` | 0.1 | 0.3 |
| `slide_zb0`, `slide_zb1` [m] | -500, 0 | -500, 0 |
| `slide_p` | 1.0 | 0.5 |
| `slide_u0`, `slide_m`, `slide_umax` | 100 m yr-1, 2, 1000 m yr-1 | same |

Build with `plastic=1 vilma=1`. The run setup is that of [Simulating our own ice](experiment-plastic.md). The closure
values come from the present-day scans of `runs/closure-cn` (the sets within the observed totals; section "Implied
mass balance: present-day tests" of [docs/dev/plastic-backend.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-backend.md)). The first LGM run of this par file
is queued on Levante (2026-11-10); no result yet.

## Present-day closure fit: tests/plastic_closure_fit.f90

[tests/plastic_closure_fit.f90](https://github.com/fesmc/kryos/blob/main/tests/plastic_closure_fit.f90) ("Stage 1b")
fits the flow closure against the observed present-day surface speed, offline, on the present-day plastic geometry
of Stage 1 (see [Calibration and ensembles](experiment-calibration.md)). That geometry is the 0 ka margin, with
Gowan's tau_b at 0 ka times the Stage 1 factors (`factor_file`), smoothed as `taub_smooth_km`, on the topography
file's bed with z_sl = 0.

```bash
make libkryos/bin/plastic_closure_fit.x plastic=1
plastic_closure_fit.x <fit.nml>
```

**Parameter sets:** for each set of the lists, the program computes the surface velocity, the flux divergence and the
margin outflow (stress capped at tau_b, as in the backend):

- `law = "yield"`: each (A, uy, b) of `A_list`, `uy_list`, `b_list`, with `m` fixed.
- `law = "coulomb_n"`: each (A, Hc, lambda_min, m) of `A_list`, `Hc_list`, `lam_list`, `m_list`, with `u0`, `umax`,
  `p`, `zb0` and `zb1` fixed.

**Scores:** over the ice cells with an observed speed above `speed_min` (default 0.1 m yr-1):

- the RMS of log10(model / observed) and the median ratio;
- the integrals [Gt yr-1] of the implied SMB, its positive part, and the outflow to land and to the ocean;
- the ablation area (the fraction of ice cells with div q - q_land < 0).

**Output:** a text table of the sets sorted by RMS, and a NetCDF file for the best set.

**`prior_file`** (one parameter set only) runs the joint inversion on that set, with `prior_var` (default `"smb"`,
[kg m-2 yr-1]) of `prior_file` as the prior, dH/dt = 0, and the `inv_*` settings named as the backend's keys (plus
`inv_lin_maxit`, `inv_sigma_def`, `inv_verbose`). The NetCDF file then also holds `speed_inv`, `smb_inv` and
`q_marine_inv`.

The namelist `&fit` has these keys: `hist`, `margin_var`, `margin_min`, `taub_file`, `taub_default`, `factor_file`,
`taub_smooth_km`, `topo_file`, `bed_var`, `mask_var`, `obs_file`, `obs_var` (default `"uxy_srf"`), `speed_min`,
`mode`, `law`, the lists, `n_glen`, `m`, `u0`, `umax`, `p`, `zb0`, `zb1`, `out`, `prior_file`, `prior_var` and the
`inv_*` settings. The repository holds no example `fit.nml`.

The figures:

```bash
julia --project=tools/levante tools/levante/plastic_closure_maps.jl <out.png> <sd_km> [ref=<smb.nc>] <label>=<set.nc>[:inv] ...
julia --project=tools/levante tools/levante/plastic_implied_pd.jl <run folder> <out folder> [age=0]
```

`plastic_closure_maps.jl` maps the fit's sets (`:inv` for the inversion). `plastic_implied_pd.jl` compares a coupled
run's implied mass balance with MAR, RACMO and the observed speeds. `plastic_lgm_map.jl mode=smb` maps it at any age.
See [Analysis and figures](experiment-analysis.md).

## Status and open issues

**Paused (2026-11-10).** The inversion is implemented and tested at present day; it needs more development before it is
effective, and the priority is the ice configurations and their sea level. The full state, the present-day results and
the possible directions are in the section "State and directions" of
[docs/dev/plastic-smb-inversion-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-smb-inversion-design.md).

Present-day Antarctica with RACMO2.4p1 as the prior (Gt/yr; ablation area: the fraction of the ice cells with a < 0):

| | net = discharge | ablation | ablation area | RMS log10 speed | wall time |
|---|---|---|---|---|---|
| forward closure (`coulomb_n`) | 1905 | -106 | 25 % | 0.557 | seconds |
| inversion, no smoothing | 1239 | -17.5 | 35 % | 0.482 | 6.5 min |
| inversion, `inv_smooth_smb_km` 100 | 1924 | -9.0 | 9.5 % | 0.511 | 84 min |
| inversion, 50 km, floor 30 kg m-2 yr-1 | 1633 | -5.8 | 3.8 % | 0.499 | 61 min |
| RACMO2.4p1 | 2049 | -1 | 0.4 % | | |

- **What works:** positive SMB nearly everywhere with RACMO's coastal pattern, a speed closer to the observed than the
  forward closure, and with the smoothing of the SMB deviation no streaks and a total near RACMO.
- **What does not yet:** the cost (about an hour per solve with the smoothing, against some 800 solves in a glacial
  run); the total below the prior where the flow law is slow; `smb_simple` is the only prior for the past; no glacial
  run yet.
- **Directions:** a fast solver for the correlated prior error (multigrid, a DCT solve, or a coarse basis); solving at
  the output slices only; a spatially varying rate factor; a RACMO/MAR-based prior scaled to the past; posterior sd maps.
