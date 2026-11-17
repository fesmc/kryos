---
title: "PlasticIceSheet"
---

The ice backend `plastic` is a perfectly-plastic ice sheet: at every update the ice
mask of a **margin chronology** and a **basal shear stress** field τ_b at the end of
the step give the thickness on the bed and sea level of the geo component. This is
how Gowan et al. (2021) built PaleoMIST, run online in kryos, so that the ice and the
solid Earth (VILMA2, see [VILMA2](backend-vilma2.md)) see each other at every coupling
step. The solver is the Fortran library
[PlasticIceSheet](https://github.com/fesmc/PlasticIceSheet) (module `plastic_ice`,
`libplastic.a`, double precision, no dependencies), a clean-room port of
PlasticIceSheet.jl.

- **Build:** `plastic=1`. `config/common.mk` expects a sibling checkout or symlink
  `PlasticIceSheet` with `libplastic/include` (`libplastic.a` and the `.mod` files),
  installed and configured by configme (`configme install kryos`); `make kryos` builds it with the same `openmp=` as
  kryos and with `lis=1` against the same LIS.
  Without it the stub stops the run ("rebuild kryos with plastic=1"). `yelmo=0` is
  enough; with VILMA2 add `vilma=1`. See [Getting started](getting-started.md).
- **Select:** `backend = "plastic"` in the ice instance group of a domain; the keys
  are those of the template `[ice_plastic]` in
  [`input/kryos_defaults.nml`](https://github.com/fesmc/kryos/blob/main/input/kryos_defaults.nml).
- **Source:** [`src/backends/ice/kryos_ice_plastic.f90`](https://github.com/fesmc/kryos/blob/main/src/backends/ice/kryos_ice_plastic.f90)
  (the backend), [`plastic/kryos_plastic_margin.f90`](https://github.com/fesmc/kryos/blob/main/src/backends/ice/plastic/kryos_plastic_margin.f90)
  (the margin chronology), [`plastic/kryos_plastic_taub.f90`](https://github.com/fesmc/kryos/blob/main/src/backends/ice/plastic/kryos_plastic_taub.f90)
  (τ_b), [`plastic/kryos_plastic_friction.f90`](https://github.com/fesmc/kryos/blob/main/src/backends/ice/plastic/kryos_plastic_friction.f90)
  (the sliding laws that follow τ_b), [`plastic/kryos_plastic_bed.f90`](https://github.com/fesmc/kryos/blob/main/src/backends/ice/plastic/kryos_plastic_bed.f90)
  (the bed factor λ of `coulomb_n` and the Coulomb law). The four helper modules use no library and are
  compiled whatever `plastic=`, so their tests run with the stub.
- **Library API:** PlasticIceSheet
  [`docs/fortran-port.md`](https://github.com/fesmc/PlasticIceSheet/blob/main/docs/fortran-port.md).
- **Par files:** `par/kryos_North_plastic.nml` (NH-32KM, fixed bed, no geo),
  `par/kryos_North_plastic_vilma2_b1.nml` and `par/kryos_NH-ANT_plastic_vilma2_b1.nml`
  (Gowan's 1-D Earth), `par/kryos_NH-ANT_plastic_vilma2.nml` and
  `par/kryos_NH-ANT-PAT_plastic_vilma2.nml` (the 3-D Earth), see
  [Par files](reference-par-files.md).

## A perfectly-plastic ice sheet

A perfectly-plastic ice sheet carries at its base exactly its yield stress τ_b: the
driving stress equals τ_b everywhere on the ice,

$$
\rho_i g H \, |\nabla z_s| = \tau_b ,
$$

so that, given the grounded-ice mask, the bed and τ_b, the surface z_s follows from a
static eikonal equation, solved from the margin inwards. No climate, flow law or time
integration enters: the margin and τ_b set the geometry. On a flat bed with uniform
τ_b this is the Nye profile, $H = \sqrt{2 \tau_b d / (\rho_i g)}$ at a distance $d$
from the margin.

Gowan et al. (2021) built PaleoMIST 1.0 (80 ka to present, every 2.5 kyr) this way
with ICESHEET 2.0: hand-drawn grounded-ice margins every 2.5 kyr; the plastic
equation along flowlines from each margin on 5 km regional grids (North America with
Greenland, Eurasia, Antarctica, Patagonia); τ_b set per **domain** (polygons drawn
from topography and surficial geology: higher in mountainous terrain, lower on
continuous sediment and for major ice streams), tuned by hand to sea-level and GIA
data, raised to maximise the LGM volume and lowered during the deglaciation, and for
Greenland and Antarctica first fitted to the present ice surface; the bed and sea
level from SELEN on a 1-D Earth. The FesmData products used here are
`IceHistory-PaleoMIST` (the ice, `f_ice`, `z_bed`, ... on the kryos grids) and
`IceHistory-PaleoMIST_taub` (the domains and their τ_b). PaleoMIST is a set of
2.5 kyr slices, each reconstructed on its own (`margin_interp = "nearest"` reproduces
that); online, kryos interpolates the margin in time and solves on the bed and sea
level of its own Earth at every step. See
[Imposed ice](howto/imposed-ice.md) and [Online ice](howto/online-ice.md).

## What a step does

At every ice update (`dt = 0`: every model step):

1. the bed `z_bed` and sea level `z_sl` from the imports (geo); until they are
   published, and always in a domain without geo, the reference bed and `z_sl = 0`;
2. the ice mask at t1 from the margin chronology;
3. τ_b at t1; every cell of the mask must have τ_b > 0, else the run stops with the
   count;
4. the solve (below); a solve that does not converge is a warning (the outer
   iteration count, the last change, the library's message) and the run goes on, bad
   input stops;
5. with `mass_balance = "implied"` or `"inverted"`: the mass balance the geometry
   needs (see [Flow closure](#flow-closure));
6. the exports.

The ice is imposed, not integrated: there is no surface, basal or calving flux, and
the mass change that the margin and the bed impose is `mb_net`, so that the budget of
the global geo and the sea level closes.

## Drive

| `drive` | Library call | Input | Needs |
|---|---|---|---|
| `"taub"` (default) | `plastic_solve(ps, z_bed, z_sl, tau_b, mask, warm=.true.)` | τ_b of the keys below | |
| `"smb"` | `plastic_solve_smb` | the imported `smb` / ρ_i (steady mass conservation) | `rheo_A > 0`, `slide_law` |

`drive = "smb"` does not converge on 32 km grids (library); the runs use `"taub"`.

## Margin chronology

`margin_file` holds the slices (x, y, time) of `margin_var` on the instance grid
(`{domain}`, `{grid_name}` are replaced); a cell is ice where the slice exceeds
`margin_H_min` (a fill value is no ice). The time axis is converted to model time from
its units: `years ago` and `BP` count into the past, a leading `k` means thousands of
years. The FesmData product `{grid_name}_IceHistory-PaleoMIST.nc` has 33 slices, 80 ka
to present, in `years ago`.

- **Variable.** `margin_var = "f_ice"`, `margin_H_min = 0.5` (default) is the ice
  fraction of Gowan's unsmoothed 5 km reconstruction on each cell: his margin at the
  grid's resolution. PaleoMIST's `H_ice` is smoothed with a 25 km Gaussian filter, so
  `H_ice > 1 m` reaches beyond the margin and makes the plastic ice too large
  (`margin_var = "H_ice"`, `margin_H_min = 1.0`).
- **Signed distance.** At init every slice becomes its signed distance to the margin
  (m, positive inside): the exact Euclidean distance between cell centres (separable
  transform, Felzenszwalb and Huttenlocher), less half a cell, so that the zero level
  lies midway between an ice cell and its ice-free neighbour and `sd > 0` gives the
  slice back exactly. A slice without ice (or without an ice-free cell) has the
  domain diagonal.
- **Interpolation** (`margin_interp`). `"sdist"` (default): the signed distance
  interpolated linearly in time between the two slices that bracket t, the mask is
  `sd > 0`. The margin moves smoothly and linearly in time, the slices come back
  exactly at their times, and a cell is ice only where one of the two slices has ice.
  `"nearest"`: the nearest slice (a tie keeps the earlier one), Gowan's offline
  behaviour.
- **Position within the cells** (`margin_position`, the library's `margin_frac`).
  `"sdist"` (default): on each ice-free cell next to the ice the margin lies at
  $-sd/\sqrt{dx\,dy}$ of the spacing from the cell centre (clamped to [0, 0.99]), with
  sd the interpolated signed distance; at a slice's time this puts the margin on the
  cell faces, between the slices it moves within the cells. `"face"` and `"centre"`
  (the ice-free cell centre, as PlasticIceSheet.jl) pass no fraction.
- **Outside the time range** of the slices (and of a τ_b series):
  `margin_extrapolate = "error"` stops, `"hold"` holds the first or last slice.

The static export `footprint` is 1 where some slice has ice: every interpolated mask
lies in it.

### Margin patches

`margin_patches` names groups (template `[margin_patch]`), each a regional chronology
that replaces `margin_file`'s in its region within its time range: `file`, `var`,
`H_min` as `margin_file`, `margin_var`, `margin_H_min`, and the region: the cells of a
region map (`region_file`, `region_var`, integer ids on the instance grid) whose id is
one of `region_ids` (each id must occur in the map). At t the mask and signed distance
are `margin_file`'s, and in the region of a patch whose slices span t (the ends
included), the patch's (`margin_set_t` of `kryos_plastic_margin`). Each chronology
interpolates between its own slices only (with the instance's `margin_interp`), so
outside the patch regions, and inside them before and after a patch's slices, the
margin is `margin_file`'s bit for bit; a patch whose first and last slices are the
base's slices at those times hands over continuously. Regions must not overlap. The
signed distance is discontinuous at a region's boundary by the two chronologies'
difference there; `margin_position = "sdist"` clamps the fraction it gives. The
footprint is `margin_file`'s, or a patch's in its region: a patch with ice where
`margin_file` has none in any slice grows the static footprint, which changes the global
gather (and so the run) from the start, before the patch's time range.

```
&ice_nor
    ...
    margin_patches = "margin_grl"
/
&margin_grl
    file        = "/albedo/work/projects/p_forclima/kryos-work/margins/NH-32KM/NH-32KM_IceHistory-PaleoGrIS-patch.nc"
    var         = "f_ice"
    H_min       = 0.5
    region_file = "kryos_data/plastic/{grid_name}_taub_regions_v2.nc"
    region_var  = "region"
    region_ids  = "2"
/
```

Greenland from PaleoGrIS 1.0 (Leger et al. 2024): `tools/levante/margin_patch.jl` writes
that patch on NH-32KM (PaleoGrIS's isochrones, 13.5 to 6.75 ka, and PaleoMIST's 15 and
5 ka slices; used in region 2 of the maps v2, Greenland).

## Basal shear stress

`taub_method` (with `drive = "taub"`):

- `"const"`: `taub_const` [Pa] everywhere;
- `"file"`: the static field `taub_var` (x, y) [Pa] of `taub_file`;
- `"domains"`: Gowan's domain map, below;
- `"coulomb"`: the Coulomb yield law, below (`mode = "surface"`).

### Gowan's domains

`taub_domain_file` is the product `{grid_name}_IceHistory-PaleoMIST_taub.nc` (FesmData
IceHistory): `taub_domain_var` (x, y), the domain id of each cell (≤ 0: none; the
product uses -9999), `taub_id_var` (n), the ids, `taub_value_var` (n), their τ_b [Pa],
and, when the file has it, `taub_time_var` (n, time; either dimension order), a series
per id, linear in time. The product has 949 domains (all regions in every file),
`domain_id` = 1000 × region + Gowan's domain number (region 1 North America with
Greenland, 2 Eurasia, 3 Antarctica, 4 Patagonia); `taub` is Gowan's base value
(`domains_max`) and `taub_t` (time, domain; years ago, 80 ka to 0, oldest first) the
value at each slice: τ_b is time dependent in PaleoMIST 1.0, so the series is used.
Ids 1000, 2000, 3000 are Gowan's "outside every polygon" (5000 Pa). In remapping,
`taub_domain` is taken as classes (the domain covering most of each cell); the 2-D
base field `taub_map` can be read with `taub_method = "file"`.

- `taub_override` names a group of a table (keys `codes`, `taub` [Pa]) that sets the
  value of some ids, constant in time, or gives an id of the map that the file has no
  value for; an id of the map without a value stops at init.
- `taub_default` [Pa]: the cells without a domain (0: none, and such a cell in the ice
  mask stops the run). The par files use 5000 Pa, Gowan's value outside every polygon.

### Factors and multipliers

Two optional NetCDF files scale the domain values (`"domains"` only; `""` = none; both
take `{domain}`, `{grid_name}`). One file can hold the domains of all grids, so one
file serves several instances. tau_b of a domain is

$$
\tau_b = \tau_{b,\mathrm{Gowan}}(t) \; e^{u_d} \; e^{v_{g(d)}(t)} ,
$$

the value after the override table, times the domain factor, times the multiplier of
its scale group. `taub_default` is scaled by neither.

| Key | File | What |
|---|---|---|
| `taub_domain_factor_file` | dimension `domain`: `domain_id(domain)`, `log_factor(domain)` (natural log $u_d$) | a factor per domain, constant in time: the present-day fit (Stage 1, `tests/plastic_taub_fit.f90`, written by `tools/levante/plastic_taub_fit_report.jl` as `taub_domain_factor.nc`; released as `$KRYOS_DATA/plastic/taub_domain_factor_v2.nc`, `_v3.nc`). An id of the product missing from the file keeps factor 1; an id of the file without τ_b stops unless its log_factor is 0. |
| `taub_scale_file` | dimensions `domain`, `scale_group`, `node`: `domain_id(domain)`, `group(domain)` (1 to the number of groups, 0 = unscaled), `time(node)` (any order, units as the margin's), `log_factor(scale_group, node)` ($v_g$) | multipliers per scale group, time dependent: the parameters of a weighted-ensemble member (`tools/levante/plastic_taub_scale.jl`, written by `plastic_ensemble.jl` as `<ensemble>/iter_01/scale/mNNN.nc`). $v_g$ is linear in time between the nodes and held outside them. An id of the product missing from the file is unscaled; an id of the file with a group but no τ_b stops. |

A file of zeros gives the unscaled τ_b bit for bit, and so does leaving the key out.
This page is the reference for how τ_b is composed; the other pages link here.

**The present-day factors (Stage 1)** are fitted offline by `tests/plastic_taub_fit.f90`
to today's ice thickness of the topography hub: the backend at t = 0 on the history's
0 ka margin and the hub's bed (z_sl = 0), τ_b = Gowan's 0 ka value × e^{u_d}, a cold
solve. The misfit is χ²(u) = Σ ((H(u) − H_obs)/σ)² over the grounded ice cells of the
margin, σ = max(`sigma_min`, `sigma_rel` H_obs) (50 m, 0.1), plus λ Σ w (u_d − u_d')²
over adjacent free domains (w the length of their shared boundary). It is minimised by
Levenberg-Marquardt (|u| ≤ `u_max` = 2.3: factors 0.1-10) for a list of λ; the report
`plastic_taub_fit_report.jl` picks λ per grid (the largest with χ²/n ≤ 1, else the
corner of the L-curve). Released: `taub_domain_factor_v3.nc` (NH-32KM, ANT-16KM,
PAT-PAL-8KM), the file of the par files and the ensembles. How to run it:
[Ensembles, Stage 1](howto/ensembles.md#stage-1-the-present-day-factors).

**The multipliers** are the parameters of the weighted ensembles: one log factor per
scale group (Laurentide, Greenland, Fennoscandia, Barents-Kara, British-Irish, West
and East Antarctica, Patagonia) and time node, sampled from a prior and weighted by
the sea-level data ([Our ice: the τ_b ensembles](science/ice-ensembles.md)). The
ES-MDA calibrations that first set them (v1, v2) are retired.

### Coulomb yield law

With `taub_method = "coulomb"` τ_b is not given but follows from the geometry (Robinson
et al. 2020, as in Yelmo):

$$
\tau_b = c_f \, \lambda(z_b - z_{sl}) \, N_\mathrm{eff}, \quad
N = \rho_i g H \max\!\left[\left(1 - \frac{H_f}{H}\right)^{p}, \delta\right], \quad
\frac{1}{N_\mathrm{eff}} = \frac{1}{N} + \frac{1}{\rho_i g H_N},
$$

N the effective pressure of Leguy et al. (2014) (H_f = (ρ_w/ρ_i) max(0, z_sl − z_b), the
bracket 0 at or below flotation) with the floor δ = `taub_delta` (default 0.05) of the
overburden, capped smoothly at ρ_i g `taub_HN`, and λ linear from
`taub_lam_min` at a bed `taub_zb0` or more below sea level to 1 at `taub_zb1` or higher
(`kryos_plastic_bed`, the same λ as `coulomb_n`). λ is geology, N (through H_f) the ocean
connection: with `taub_lambda_bed = "relaxed"` (default) λ is a fixed field of the
reference bed relaxed by local isostasy, z_b,ref + (ρ_i/ρ_m) H_ice,pd (ρ_m = 3300 kg m⁻³,
relative to present sea level), computed at init and held in time, so that the
isostatically depressed interior beds of today's ice sheets do not read as marine;
`"current"` takes the bed and sea level of each solve. In the plastic equation this is
a slope, |∇z_s| = c_f λ N_eff/(ρ_i g H): the constant c_f λ where H ≪ H_N on land (a wedge),
c_f λ δ H_N/(δ H + H_N) at flotation, and for H ≫ H_N the Nye profile of τ_b = c_f λ ρ_i g H_N.

The floor makes the solve well posed. Without it (N = 0 at flotation) the slope vanishes at
flotation, and a marine trough whose flotation surface rises inland (Weddell, Ross) can sit
wholly at flotation or grounded above it: the discrete problem has several solutions, and
which one the sweep reaches flips on decimetre changes of the bed. Coupled to VILMA2 (the
elastic and geoid response to the change) the Antarctic mass then jumped by up to
1.5 × 10⁶ Gt per 100-yr step, alternating between the two states. With δ > 0 every grounded
cell has a positive slope, the solution is unique and continuous in the bed and sea level
(the library's uniqueness argument, PlasticIceSheet docs/fortran-port.md); δ = 0.05 removed
the jumps (max |ΔM| 2.5 × 10⁵ Gt, no alternation), δ = 0.02 did not. Nothing depends on
time explicitly. The backend computes μ = c_f λ on the bed and sea level of each solve and
calls the library's `plastic_solve_coulomb` (mode `"surface"` only), which returns τ_b
(output `tau_b`; `taub_mu` is μ).

c_f is one value per **region** (ice sheet), constant in time: `taub_region_file` is a map
`region` (x, y) of region ids (0 = none; a cell of the ice mask without a region stops the
run), `taub_cf_file` has `region_id(region)`, `log_cf(region)` (natural log; an id not on the
grid is allowed, a region of the map without c_f stops). The maps
`$KRYOS_DATA/plastic/{grid_name}_taub_regions_v2.nc` (NH-32KM, ANT-16KM, PAT-PAL-8KM; written by
`tools/levante/plastic_cf_regions.jl`, which also writes c_f files) have 1 Laurentide (North
America with the Cordillera and Innuitian), 2 Greenland, 3 Fennoscandia (with the
British-Irish ice), 4 Antarctica, 5 Patagonia, 6 BarentsKara (the Barents and Kara seas with
their islands), from the geography of Gowan's domains (no τ_b of the product): whole domains,
Greenland and BarentsKara by the majority of their cells in REGIONS codes; the cells off the
domains take the nearest region, so that v2 covers the grid (any margin chronology). v1
(`{grid_name}_taub_regions.nc`, the runs C1, C2) has 1-5 with all of Eurasia as region 3. `taub_smooth_km`
smooths μ. The keys of `"domains"` (override, factors, scale groups) do not apply.

### Smoothing

`taub_smooth_km` (default 50; `drive = "taub"`, any method; `"coulomb"`: μ): τ_b at t, after the
factors, is replaced on the cells with τ_b > 0 by its Gaussian mean over those cells
(sd `taub_smooth_km`, within 3 sd; a normalised convolution), for the solve and the
flow closure alike. The steps between domains become ramps of a few sd (the ice
averages its basal stress over several thicknesses); the value of a uniform region is
kept. 0 = none. A par file that leaves the key out gets 50 km.

## Solver

The keys of the library's `plastic_par_t`, passed through:

- `mode = "surface"` (default, recommended): the plastic surface equation solved for
  z_s directly by one monotone fast-sweeping solve, with the bed relief and no outer
  iteration (ICESHEET's formulation; not in PlasticIceSheet.jl). `"flat"`: no bed
  slope. `"full"`: w = H² with the bed slope by a damped fixed point; it does not
  converge on 32 km beds.
- `relax`, `warm_frac_max`, `n_outer`, `outer_tol`: the outer iteration and the warm
  start of `"full"` (a warm start only if at most `warm_frac_max` of the ice cells
  are new); `max_sweeps`, `tol`: the sweeps.
- `rho_i`, `rho_w`, `g`: 910, 1028, 9.81, the kryos `Earth` constants (library
  defaults 917, 1025, 9.80665). A density other than the run's `rho_ice` / `rho_sw`
  is logged as a warning: the ice mass the globe gathers for VILMA2 is `H_ice` times
  the run's `rho_ice`.
- `H_min`: floor of H inside the mask.

## Exports

| Export | Value |
|---|---|
| `H_ice`, `z_srf` | the solve (`z_srf`: the library's on the ice, max(z_bed, z_sl) elsewhere) |
| `z_bed`, `z_sl` | the bed and sea level the solve used |
| `f_ice` | 1 where H > 0 |
| `f_grnd` | 1 on the ice (plastic ice is grounded by construction), the flotation rule elsewhere |
| `H_grnd` | H − (ρ_w/ρ_i)(z_sl − z_bed) |
| `dHidt` | (H − H_prev)/(t1 − t0), mean over the step |
| `mb_net` | ρ_i dHidt [kg m-2 yr-1], mean over the step |
| `smb`, `cmb` | 0, or the mass balance of the flow closure (`mass_balance`), mean over the step |
| `bmb` | 0 |
| `M_ice` | ρ_i Σ H dx dy [Gt] |
| `z_bed_ref`, `H_ice_pd` | static: the domain's reference bed and present-day ice |
| `H_ice_ref` | the present-day ice (`H_ice_pd`), put with the geometry at every update |
| `footprint` | static: 1 where some margin slice has ice (a patch's slice in its region) |
| `ux_s`, `uy_s`, `smb_tot`, `cmb_tot` | with `mass_balance = "implied"` or `"inverted"`: the surface velocity of the closure; smb and cmb summed over the grid [Gt yr-1] |

`H_ice_ref` is not the ice the backend produces but a fixed reference: today's ice of
the domain's hub topography. The solid Earth measures its load against it: for the
global VILMA2 the core builds the reference load from each domain's present-day ice
(`H_ice_pd`), bed (`z_bed_ref`) and `footprint` merged with the background ice
outside the footprints, and the rows ice → geo pass `z_bed_ref`, `H_ice_ref` at init
(see [Global hub](global-hub.qmd)). These exports are what a geo backend needs from any
ice. The backend also binds the flotation rule (grounded where ρ_i H ≥ ρ_w (z_sl −
z_bed), with its own densities).

**Imports:** `z_bed`, `z_sl` (optional: the reference bed and 0 until published);
`smb` only with `drive = "smb"` or `mass_balance = "inverted"`. `T_srf` and `bmb_shlf`
are not read, so a domain needs no climate, surface or shelf.

**Restart:** `<instance>_restart.nc` with H, w = H², the mask w belongs to and whether
a warm start is usable (the library's warm-start state), `H_prev` and `z_srf`,
lossless. `prepare` is not implemented, so no experiment hook works on this backend.

## Output

The exports go to `<instance>.nc` (`out_2D`) and `<instance>_ts.nc` (`out_1D`), table
`ice`. The backend's own stream `plastic` (`out_plastic`, `<instance>_plastic.nc`,
table [`input/kryos-variables-ice-plastic.md`](https://github.com/fesmc/kryos/blob/main/input/kryos-variables-ice-plastic.md))
holds the state of the last solve and its inputs:

| Variable | Units | What |
|---|---|---|
| `H_ice`, `z_srf`, `z_bed`, `z_sl` | m | the solve and the bed and sea level it used |
| `mask` | 1 | ice mask of the margin chronology |
| `sdist` | m | signed distance to the margin (positive inside) |
| `tau_b` | Pa | basal shear stress |
| `taub_domain` | 1 | τ_b domain (0 = none) |
| `taub_region`, `taub_mu` | 1 | `"coulomb"`: the region of c_f, μ = c_f λ |
| `div_q`, `q_land`, `q_marine` | m yr-1 | `"implied"`, `"inverted"`: flux divergence and the outflow through the margin faces to land and to the ocean (ice equivalent) |
| `u_inv`, `u_law` | m yr-1 | `"inverted"`: the speed of the inversion and of the flow law |
| `converged`, `n_outer_used`, `stress_res_median`, `stress_res_max` | | 0-D: convergence and the stress residual \|ρ_i g H \|∇z_s\| − τ_b\| / τ_b over the ice (1 at the first ice cell by construction) |
| `abl_margin`, `discharge` | Gt yr-1 | 0-D, `"implied"`, `"inverted"` |
| `inv_iters`, `inv_J_law`, `inv_J_smb`, `inv_J_cold`, `inv_n_bound` | | 0-D, `"inverted"` |

## Flow closure

A plastic sheet has no velocity. With a flow law, `rheo_A > 0` (Glen rate factor
[Pa-n yr-1], exponent `rheo_n`) and a sliding law `slide_law`, the backend computes the
flux the geometry carries and from it the mass balance it needs (`mass_balance`).
`drive = "smb"` needs the same closure.

| `slide_law` | Law | Keys |
|---|---|---|
| `"none"` (default) | no closure | |
| `"linear"` | library: linear, uniform coefficient | `slide_beta` |
| `"weertman"` | library: Weertman | `slide_beta`, `slide_m` |
| `"coulomb"` | library: regularised Coulomb | `slide_C`, `slide_u0`, `slide_m`, `slide_umax` |
| `"yield"` | library Weertman with β per cell, β = τ_b / u_y^(1/m), u_y = `slide_uy` (τ_b / 100 kPa)^`slide_uy_b` | `slide_uy`, `slide_uy_b`, `slide_m`; `drive = "taub"` |
| `"coulomb_n"` | library regularised Coulomb with C per cell, C = τ_b λ N / (ρ_i g `slide_Hc`) | `slide_Hc`, `slide_lambda_min`, `slide_zb0`, `slide_zb1`, `slide_p`, `slide_u0`, `slide_m`, `slide_umax`; `drive = "taub"` |

`"yield"` and `"coulomb_n"` (`kryos_plastic_friction.f90`) set their coefficient from
the imposed τ_b at every step. In `"coulomb_n"`, N = ρ_i g H (1 − H_f/H)^`slide_p` is
the effective pressure (Leguy et al. 2014; H_f = (ρ_w/ρ_i) max(0, z_sl − z_bed), N = 0
at or below flotation) and λ rises linearly from `slide_lambda_min` at a bed
`slide_zb0` or more below sea level to 1 at `slide_zb1` or higher (weak marine beds).
At τ = τ_b the sliding speed then depends on N and the bed alone: slow under thick
ice, fast where the ice thins below `slide_Hc`/λ or nears flotation, at most
`slide_umax` (`"coulomb_n"` needs one).

`mass_balance = "implied"` takes the divergence of the face fluxes of the library's
steady SMB-driven solver on the new geometry (with `drive = "taub"` the driving stress
of the fluxes capped at τ_b), the outflow through the margin to land ablated in the
margin cells and to the ocean calved; `mb_smooth_km` smooths it conservatively.
`mass_balance = "inverted"` (`drive = "taub"`) inverts jointly for the speed and the
mass balance against the closure's speed and the imported `smb` as the prior (e.g.
from [smb_simple](backend-smb-simple.md)), with the `inv_*` keys. Both are in
development: see [Implied mass balance](implied-mb.md).

## Keys

Template `[ice_plastic]` of `input/kryos_defaults.nml`; the instance group also has
the core keys of `[instance]` (`backend`, `grid`, `dt`, `out_2D`, ...).

| Key | Default | |
|---|---|---|
| `drive` | `"taub"` | `"taub"` or `"smb"` (needs `rheo_A`, `slide_law`) |
| `mode` | `"surface"` | `"surface"`, `"flat"`, `"full"` |
| `relax`, `warm_frac_max`, `n_outer` | 0.5, 0.1, 100 | outer iteration and warm start (`"full"`) |
| `max_sweeps`, `tol`, `outer_tol` | 200, 1e-6, 1e-3 | sweeps per eikonal solve; sweep tolerance (change of w = H² [m²]; of z_s [m] in `"surface"`); outer tolerance [m] |
| `rho_i`, `rho_w`, `g` | 910, 1028, 9.81 | [kg m-3], [kg m-3], [m s-2] |
| `H_min` | 1.0 | [m] floor of H inside the mask |
| `margin_file` | `""` | required: the margin slices |
| `margin_var`, `margin_H_min` | `"f_ice"`, 0.5 | ice where `margin_var` > `margin_H_min` |
| `margin_interp` | `"sdist"` | `"sdist"`, `"nearest"` |
| `margin_extrapolate` | `"error"` | `"error"`, `"hold"` |
| `margin_position` | `"sdist"` | `"sdist"`, `"face"`, `"centre"` |
| `margin_patches` | `""` | groups of regional chronologies (template `[margin_patch]`: `file`, `var`, `H_min`, `region_file`, `region_var`, `region_ids`) |
| `taub_method` | `"const"` | `"const"`, `"file"`, `"domains"`, `"coulomb"` |
| `taub_const` | 0.0 | [Pa] `"const"` (> 0) |
| `taub_file`, `taub_var` | `""`, `"taub"` | `"file"` |
| `taub_domain_file` | `""` | `"domains"`: the product |
| `taub_domain_var`, `taub_id_var`, `taub_value_var`, `taub_time_var` | `"taub_domain"`, `"domain_id"`, `"taub"`, `"taub_t"` | `"domains"`: its variables |
| `taub_override` | `""` | `"domains"`: group of a table (`codes`, `taub`) |
| `taub_default` | 0.0 | [Pa] `"domains"`: cells without a domain (0: none) |
| `taub_domain_factor_file` | `""` | `"domains"`: factor per domain (Stage 1) |
| `taub_scale_file` | `""` | `"domains"`: multipliers per scale group |
| `taub_smooth_km` | 50.0 | [km] Gaussian sd of the τ_b smoothing (`"coulomb"`: of μ); 0 = none |
| `taub_region_file`, `taub_cf_file` | `""`, `""` | `"coulomb"`: the region map, c_f per region |
| `taub_HN`, `taub_p` | 0.0, 1.0 | `"coulomb"`: [m] cap of N (required, > 0); exponent of N in [0, 1] |
| `taub_delta` | 0.05 | `"coulomb"`: floor of N as a fraction of the overburden ρ_i g H, in (0, 1] |
| `taub_lam_min`, `taub_zb0`, `taub_zb1` | 1.0, -300.0, 200.0 | `"coulomb"`: λ (1 = no bed dependence) |
| `taub_lambda_bed` | `"relaxed"` | `"coulomb"`: the bed of λ, `"relaxed"` (fixed) or `"current"` |
| `rheo_A`, `rheo_n` | 0.0, 3.0 | Glen rate factor [Pa-n yr-1] (0: no closure), exponent |
| `slide_law` | `"none"` | see [Flow closure](#flow-closure) |
| `slide_beta`, `slide_m` | 0.0, 1.0 | linear, weertman coefficient; stress exponent |
| `slide_C`, `slide_u0`, `slide_umax` | 0.0, 0.0, 0.0 | [Pa], [m yr-1], [m yr-1] (0 = no cap) |
| `slide_Hc`, `slide_lambda_min`, `slide_zb0`, `slide_zb1`, `slide_p` | 0.0, 1.0, -500.0, 0.0, 0.5 | `"coulomb_n"` |
| `slide_uy`, `slide_uy_b` | 0.0, 0.0 | `"yield"` |
| `mass_balance` | `"none"` | `"none"`, `"implied"`, `"inverted"` |
| `mb_smooth_km` | 0.0 | [km] `"implied"` |
| `inv_sigma_u`, `inv_smooth_km`, `inv_smooth_smb_km` | 0.7, 0.0, 50.0 | `"inverted"` |
| `inv_sigma_smb_rel`, `inv_sigma_smb_min` | 0.5, 100.0 | `"inverted"` (floor in kg m-2 yr-1) |
| `inv_cold`, `inv_sigma_cold`, `inv_n_iter` | `.true.`, 10.0, 50 | `"inverted"` |
| `out_plastic` | `"none"` | schedule of `<instance>_plastic.nc` |

## Configuration

Example: the NH-32KM instance of `par/kryos_NH-ANT_plastic_vilma2_b1.nml` (solver keys
at their defaults left out here):

```fortran
&ice_nor
    backend            = "plastic"
    grid               = ""             ! the hub grid, NH-32KM
    dt                 = 0.0            ! every step
    out_2D             = "out_2D"
    out_1D             = "out_1D"
    drive              = "taub"
    mode               = "surface"
    margin_file        = "fesmdata/{domain}/{grid_name}/{grid_name}_IceHistory-PaleoMIST.nc"
    margin_var         = "f_ice"
    margin_H_min       = 0.5
    margin_interp      = "sdist"
    margin_extrapolate = "error"
    taub_method        = "domains"
    taub_domain_file   = "fesmdata/{domain}/{grid_name}/{grid_name}_IceHistory-PaleoMIST_taub.nc"
    taub_default       = 5000.0         ! [Pa] Gowan's value outside every polygon
    taub_domain_factor_file = ""        ! Stage 1 factors; "" = none
    taub_scale_file    = ""             ! calibration multipliers; "" = none
    out_plastic        = "out_2D"       ! ice_nor_plastic.nc
/
```

The clock of these runs steps 100 yr from -80 ka to 0 (one plastic solve and one
VILMA2 step per 100 yr). On NH-32KM a `"surface"` solve takes about 10 ms.
