# ice backend `plastic`

A perfectly-plastic ice sheet as a kryos ice model, so that kryos reproduces Gowan et al. (2021) PaleoMIST
online: the margin chronology and the basal shear stress (tau_b) domains of PaleoMIST, and the bed and sea
level of the geo backend (VILMA2) at every coupling step. The solver is the Fortran library PlasticIceSheet
(module `plastic_ice`, `libplastic.a`, dp, no dependencies), the port of PlasticIceSheet.jl.

## Files

| File | What |
|---|---|
| `src/backends/ice/kryos_ice_plastic.f90` | the backend (`plastic_ice_t`), built with `plastic=1` |
| `src/backends/ice/kryos_ice_plastic_stub.f90` | the stub (`plastic=0`): the factory stops with "rebuild kryos with plastic=1" |
| `src/backends/ice/plastic/kryos_plastic_margin.f90` | the margin chronology: signed distance, interpolation (no library) |
| `src/backends/ice/plastic/kryos_plastic_taub.f90` | tau_b: const, gridded field, domains (no library) |
| `input/kryos_defaults.nml` `[ice_plastic]` | the keys of the instance group, with their defaults |
| `input/kryos-variables-ice-plastic.md` | the variables of the backend's own stream `plastic` |
| `par/kryos_North_plastic.nml` | NH-32KM, fixed bed (no geo) |
| `par/kryos_North_plastic_vilma2.nml` | NH-32KM with the global VILMA2 Earth of the science case (3-D) |
| `par/kryos_North_plastic_vilma2_b1.nml` | NH-32KM with VILMA2 on Gowan's 1-D Earth and B1's reference: the run scored against GAPSLIP beside B1 |
| `tests/test_plastic_margin.f90` | margin and tau_b (runs with `plastic=0`) |
| `tests/test_ice_plastic.f90` | the backend through the ice wrappers (`plastic=1` only) |
| `tests/plastic_slices.f90` | the offline check against PaleoMIST per slice (`plastic=1`) |
| `tools/levante/plastic_maps.jl`, `pm_kryos_out.jl`, `plastic_score.sh` | maps of the offline check; a kryos VILMA2 stream as a standalone `vilma_out.nc`; the GAPSLIP scoring with `pm_compare.jl` |
| `tools/levante/plastic_taub_scale.jl`, `plastic_esmda.jl`, `plastic_esmda.sh`, `plastic_esmda.toml` | the tau_b calibration: scale groups and multiplier files; ES-MDA; its jobs; its config |

Build: `plastic=1` adds `INC_PLASTIC` / `LIB_PLASTIC` (`config/common.mk`: the sibling `PlasticIceSheet`,
`libplastic/include` with `libplastic.a` and the `.mod` files; `make -C PlasticIceSheet static
config=dkrz_levante_ifx`). `yelmo=0` is enough for these runs. The library's API: PlasticIceSheet
`docs/fortran-port.md` (branch `fortran-port`).

## What a step does

`update(st)`, every ice step:

1. the bed `z_bed` and sea level `z_sl` from the imports (geo); until they are published (a domain without geo:
   always) the reference bed and `z_sl = 0`;
2. the ice mask at `t1` from the margin chronology (below);
3. tau_b at `t1` (below); every ice cell must have a positive tau_b, else the run stops naming the count;
4. `plastic_solve(ps, z_bed, z_sl, tau_b, mask, warm=.true., stat)` (`"surface"`: one sweep for z_s, warm
   ignored; `"full"`: warm-started from the last `w = H**2`) (`drive = "smb"`: `plastic_solve_smb` with the
   imported smb / rho_i, and tau_b the driving stress rho_i g H |grad z_srf| for the diagnostics);
   `stat = PLASTIC_NOT_CONVERGED` is a warning with the outer iteration count, the last change and the
   library's message, and the run goes on; `PLASTIC_BAD_INPUT` stops;
5. the velocity and the smb it implies (`plastic_velocity`, `plastic_smb_implied` with the speed, as
   Julia's `smb_from_velocity`) when the rheology and the sliding law are set;
6. the exports.

The ice is imposed, not integrated: no surface, basal or calving flux. The mass change the margin and the bed
impose is `mb_net`, so the budget of the global geo (and the sea level) closes.

## Exports (contracts-design 2.1)

| Export | Value |
|---|---|
| `H_ice`, `z_srf` | the solve (`z_srf`: the library's on the ice, `max(z_bed, z_sl)` elsewhere) |
| `z_bed`, `z_sl` | the bed and sea level the solve used |
| `f_ice` | 1 where H > 0 |
| `f_grnd` | 1 on the ice (plastic ice is grounded by construction), the flotation rule elsewhere |
| `dHidt` | (H - H_prev)/(t1 - t0), MEAN |
| `smb`, `bmb`, `cmb` | 0, MEAN |
| `mb_net` | rho_i dHidt [kg m-2 yr-1], MEAN |
| `M_ice` | rho_i sum(H) dx dy [Gt] |
| `H_grnd` | H - (rho_w/rho_i)(z_sl - z_bed) (as Yelmo's) |
| `footprint` | static: 1 where some margin slice has ice (every mask the interpolation gives lies in it) |
| `z_bed_ref`, `H_ice_pd` | static: the domain's reference bed and ice (present day) |
| `H_ice_ref` | the reference ice (present day), with the geometry (the row ice -> geo reads it at init) |
| `ux_s`, `uy_s` | only with rheology and sliding: the velocity of `plastic_velocity` at the cell centres (diagnostic) |

`z_bed_ref`, `H_ice_ref` (rows ice -> geo, init), `footprint`, `H_ice_pd` and the flotation binding (the gather
of a host-resolved geo, `gia_check`) are what a geo backend needs from any ice; they are not specific to the
plastic model.

Flotation binding (`caps%has_flotation`): grounded where rho_i H >= rho_w (z_sl - z_bed), with the backend's
densities.

Own stream `plastic` (`out_plastic`, `<instance>_plastic.nc`, table `ice-plastic`): `H_ice`, `z_srf`, `z_bed`,
`z_sl`, `mask`, `sdist`, `tau_b`, `taub_domain` (domains), `ux`, `uy`, `smb_implied` (velocity), and the 0-D
`converged`, `n_outer_used` and the stress residual `stress_res_median`, `stress_res_max`
(`plastic_stress_residual`: |rho_i g H |grad z_srf| - tau_b| / tau_b over the ice, cell-centred, so 1 at the first
ice cell by construction) of the last solve, computed when the stream is written.

## Margin chronology

The slices of `margin_var` (default `f_ice`) in `margin_file` (`{domain}`, `{grid_name}` replaced, as
`recon_path`), on the instance grid, with ice where the variable exceeds `margin_H_min` (a fill value is no
ice). The time axis is converted to model time from its units: `years ago` and `BP` count into the past, a
leading `k` is thousands of years (the FesmData PaleoMIST product `{grid}_IceHistory-PaleoMIST.nc`: 33
slices, 80 ka to present, `years ago`).

At init every slice becomes its signed distance to the margin, in m, positive inside: the exact Euclidean
distance between cell centres (separable transform, Felzenszwalb and Huttenlocher), from an ice cell to the
nearest ice-free cell and vice versa, less half a cell so that the zero level lies midway between them and
`sd > 0` is the slice again. A slice without ice (or without ice-free cells) has the domain diagonal.

At time t, `margin_interp`:

- `"sdist"` (default): the signed distance interpolated linearly between the two slices that bracket t; the
  mask is `sd > 0`. The margin migrates smoothly and linearly in time (two concentric disks: the radius is
  linear in t, tested), the slices come back exactly at their times, and a cell is ice only where one of the
  two slices has ice.
- `"nearest"`: the nearest slice (a tie keeps the earlier one): Gowan's offline behaviour, for comparison.

Position within the cells (`margin_position`, the library's `margin_frac`, PlasticIceSheet docs/fortran-port.md
"Margin position"): `"sdist"` (default) places the margin on each ice-free cell next to the ice at -sd / sqrt(dx dy) of
the spacing from its centre (`plastic_margin_frac_from_sdist`, clamped to [0, 0.99]), with sd the interpolated signed
distance (positive in the ice, m, cell centres: the library's convention). sd is the centre-to-centre distance to the
nearest opposite cell less half a cell, so at a slice's time it is the distance to the staircase of cell faces and the
margin lies at the faces (0.5); between the slices it moves within the cells, linearly in time. `"face"` and `"centre"`
(the ice-free cell centre, PlasticIceSheet.jl and this backend until 2026-11-09) set `par%margin` and pass no fraction.

Which variable: PaleoMIST's `H_ice` is Gowan's thickness smoothed with a 25 km Gaussian filter, so `H_ice > 1 m`
reaches beyond his margin, and the plastic ice there makes the NH volume too large (most of the 15-50 % excess
below). The product's `f_ice` is the ice fraction of the unsmoothed 5 km reconstruction on each cell, and
`f_ice > 0.5` is Gowan's margin at the grid's resolution: the default. `margin_var = "H_ice"`, `margin_H_min = 1` keeps
the earlier margin (the results below until 2026-11-09 used it).

Outside the slices' time range: stop (`margin_extrapolate = "error"`), or hold the first or last slice
(`"hold"`; also for a tau_b series).

## tau_b

`taub_method`:

- `"const"`: `taub_const` [Pa] everywhere;
- `"file"`: the static field `taub_var` (x, y) [Pa] of `taub_file`;
- `"domains"` (Gowan's): the product `{grid}_IceHistory-PaleoMIST_taub.nc` (FesmData, not yet merged):
  `taub_domain_var` (x, y) the domain id of each cell (<= 0 = none; the product: -9999), `taub_id_var` (n) the ids, `taub_value_var` (n) their tau_b
  [Pa], and, when the file has it, `taub_time_var` (n, time; either dimension order) a series per id, linear
  in time. The group named by `taub_override` (keys `codes`, `taub`, as the `till_cf` group of the yelmo
  backend) sets the value of some ids, constant in time, or gives an id of the map the file has no value for;
  an id of the map without a value stops at init. Cells without a domain take `taub_default` (0 = none: such
  a cell in the ice mask stops the run).

The product (949 domains, all regions in every file): `domain_id` = region*1000 + Gowan's id, `taub` (Gowan's
domains_max), `taub_t` (time, domain) in NetCDF order, years ago, 80 ka to 0, oldest first: tau_b is time
dependent in PaleoMIST 1.0, so the series is used, linear between its records. Ids 1000, 2000, 3000 are Gowan's
domains outside every polygon (5000 Pa); cells without a domain take `taub_default` (the par files: 5000 Pa).
Its 2-D base field `taub_map` can be read with `taub_method = "file"`.

`taub_scale_file` (`"domains"` only; `""` = none): multipliers of tau_b per scale group, time dependent (the calibration
below). NetCDF, dimensions `domain`, `scale_group`, `node`: `domain_id(domain)`, `group(domain)` (the scale group of each id,
1 to the number of groups, 0 = unscaled), `time(node)` (any order, its units as the margin's: `years ago`, ...),
`log_factor(scale_group, node)` (natural log of the multiplier). tau_b of a domain is its value (the series, after the
override table) times exp(log_factor of its group at t), log_factor linear in time between the nodes and held outside them;
the `taub_default` of the cells without a domain is not scaled. An id of the product missing from the file is unscaled; an id
of the file with a group but no tau_b stops. Init logs the number of groups, of scaled domains and nodes and the factor range.
Without the key the backend is unchanged (bit for bit, below).

## Passes (core-design 3.2)

| Pass | What the backend does |
|---|---|
| `PASS_INIT` | the first solve at `time_init` on the reference bed (`z_sl = 0`: no import is published yet), cold (`warm = .false.`); the static exports, the geometry and `M_ice` |
| `PASS_GEOMETRY` | the solve on the bed and sea level geo gave (cold start block a; with a host-resolved geo again in block c', on the spun-up Earth) and the geometry |
| `PASS_SPINUP` | the MEAN exports of the virtual previous step: no change (0), as Yelmo's equilibrium-free start |
| `PASS_RESTART` | after `restart_read` (H, w, H_prev, z_srf) and the core's restore of imports and exports: the inputs of the restored time (bed, sea level, mask, tau_b) and the diagnostics `ux_s`, `uy_s`; no state export |
| others | nothing |

`prepare` is not implemented (the default): the backend reads its imports in `update`, so an experiment hook
cannot change them (no experiment works on this backend).

Restart: `<instance>_restart.nc` with `H`, `w`, the mask `w` belongs to and whether it is usable (the
library's warm-start state, `ps%wk%mask_prev`, `ps%has_warm`, so that a restarted `"full"` solve starts as the
continuous one), `H_prev` and `z_srf`, lossless; the read checks the grid and the time.

## Keys (`[ice_plastic]`)

| Key | Default | |
|---|---|---|
| `drive` | `"taub"` | `"taub"` or `"smb"` (needs `rheo_A`, `slide_law`) |
| `mode` | `"surface"` | passed to the library as it is: `"surface"` (the plastic surface equation for z_s, one monotone sweep, with the bed), `"flat"` (no bed slope), `"full"` (w = H**2 with the bed slope by a fixed point; does not converge on 32 km beds) |
| `relax`, `warm_frac_max`, `n_outer`, `max_sweeps`, `tol`, `outer_tol` | 0.5, 0.1, 100, 200, 1e-6, 1e-3 | the library's iteration controls |
| `rho_i` | 910 | [kg m-3] the kryos `Earth` constant and the PaleoMIST runs (library default 917) |
| `rho_w`, `g` | 1028, 9.81 | the kryos `Earth` constants (library defaults 1025, 9.80665) |
| `H_min` | 1 | [m] floor of H inside the mask |
| `margin_file`, `margin_var`, `margin_H_min` | `""`, `"f_ice"`, 0.5 | the margin slices: ice where `margin_var` > `margin_H_min` (`"H_ice"`, 1 m: the smoothed thickness) |
| `margin_interp`, `margin_extrapolate` | `"sdist"`, `"error"` | see above |
| `margin_position` | `"sdist"` | the margin within the cells (below): `"sdist"`, `"face"`, `"centre"` |
| `taub_method`, `taub_const` | `"const"`, 0 | |
| `taub_file`, `taub_var` | `""`, `"taub"` | `"file"` |
| `taub_domain_file`, `taub_domain_var`, `taub_id_var`, `taub_value_var`, `taub_time_var` | `""`, `"taub_domain"`, `"domain_id"`, `"taub"`, `"taub_t"` | `"domains"` |
| `taub_override`, `taub_default` | `""`, 0 | `"domains"` |
| `taub_scale_file` | `""` | `"domains"`: the multipliers per scale group (above) |
| `rheo_A`, `rheo_n` | 0, 3 | Glen rate factor [Pa-n yr-1] (0: no velocity), exponent |
| `slide_law` | `"none"` | `"linear"`, `"weertman"` (`slide_beta`, `slide_m`), `"coulomb"` (`slide_C`, `slide_u0`, `slide_m`) |
| `out_plastic` | `"none"` | schedule of the stream `plastic` |

A density other than the run's `rho_ice` / `rho_sw` is logged as a warning: the ice mass the globe gathers for
VILMA2 is the export `H_ice` times the run's `rho_ice`.

## Imports

`z_bed`, `z_sl` (optional: the reference bed and 0 until published); `smb` only with `drive = "smb"` (then
required); `T_srf`, `bmb_shlf` are not read and may stay unwired, so a domain needs no climate, surface or
shelf.

## Tests

- `tests/test_plastic_margin.f90` (any build): the distance transform against brute force (dx /= dy); the
  signed distance of a slice; two disks: the slices at their times, the margin within one cell of the linear
  radius at nine times, the area halfway, the footprint, `nearest`; a margin file in `kyr BP` with a fill
  value; hold and the stop outside the range; tau_b domains (static, override, a series in either dimension
  order, an id added by the table, the stops); the tau_b scale (groups by id, an id added by the table, linear in
  time between unordered nodes, held outside, times the series, the default unscaled, factor 1 bit for bit, the
  stops).
- `tests/test_ice_plastic.f90` (`plastic=1`): the backend on a toy grid (flat bed, two disks, tau_b const,
  `mode = "flat"`): the plastic profile at the centre, the masks, the static exports, an update (dHidt,
  mb_net, M_ice), a restart round trip stepped on bit for bit, the flotation binding.
- `tests/test_build.f90`: the factory admits `plastic` when linked.

## Results (Levante, 2026-11-09)

Build: kryos `plastic=1 yelmo=0 vilma=1` (stack S11 dependencies) against PlasticIceSheet `fortran-port`
(`/work/ba1442/robinson/PlasticIceSheet-fortran/ifx`); runs in `/work/ba1442/robinson/kryos-ice-plastic/runs`.
tau_b: `NH-32KM_IceHistory-PaleoMIST_taub.nc` (FesmData IceHistory v2.3.0),
domains with the time series, `taub_default` 5000 Pa; the margin `H_ice` > 1 m (the default is now `f_ice` > 0.5).

**Fixed bed** (`par/kryos_North_plastic.nml`: today's bed, z_sl = 0, 80 ka -> 0, dtt = 100 yr, 800 solves): 25 s on
4 cores; `"surface"` 10 ms, `"flat"` 4 ms per solve; every solve converged; the ice mask is PaleoMIST's (H > 1 m) at
every slice. Restart (flat, tau_b const): the 40 ka bundle continued to 0 equals the continuous run bit for bit
(1284 variables and records).

The fixed-bed, offline and online results in this section up to "Baseline" used the margin `H_ice` > 1 m at the
ice-free cell centres (the library before `margin_frac`).

**Offline check** (`tests/plastic_slices.f90`): at each of the 33 slices PaleoMIST's own bed (relative to the sea
level of the time, so z_sl = 0), its margin, Gowan's tau_b of that time, a cold solve, against its `H_ice` over its
ice. Volume in m SLE (rho_i/rho_sw = 910/1028, 3.618e14 m2 ocean):

| mode | bias [m] | rms [m] | stress residual median | 80 ka | 65 ka | 45 ka | 20 ka | 0 ka |
|---|---|---|---|---|---|---|---|---|
| PaleoMIST | | | | 12.9 | 53.7 | 40.6 | 112.5 | 7.5 |
| surface | +257 to +502 | 309 to 558 | 0.025 to 0.14 | 18.7 | 65.6 | 49.2 | 130.0 | 10.5 |
| flat | +285 to +550 | 425 to 650 | 0.39 to 0.83 | 19.2 | 67.7 | 50.4 | 132.6 | 10.6 |
| full (none converged) | -196 to +306 | 355 to 580 | | 16.6 | 57.8 | 44.0 | 100.1 | 8.6 |

`"surface"` solves the plastic equation (median stress residual 3 %, the maximum 1 at the first ice cells), and
still makes the ice 260-500 m thicker than PaleoMIST, the NH volume 15-50 % larger. The misfit is broad and
positive over the interiors, larger towards the margins and on the Arctic islands, and Greenland is +300 m today
(figure: `"surface"` top, `"flat"` bottom, 20 ka and 0 ka). Not the solver then: tau_b (values, units or how
ICESHEET uses them), the resolution (Gowan: 5 km, PaleoMIST smoothed to 25 km) or the margin. The unconverged
`"full"` result being closer is not evidence (no fixed point).

![H - PaleoMIST](figures/plastic-backend/plastic_dH.png)

**Online with VILMA2** (`par/kryos_North_plastic_vilma2_b1.nml`, mode `"surface"`: Gowan's 1-D Earth, lmax 170, B1's
reference; NH-32KM plastic, the rest of the world PaleoMIST's history; Earth and ice step 100 yr): 7 min 25 s on 64
threads (shared partition), of which init and VILMA2's relaxation to the 80 ka load about 2 min and the 800 steps
333 s, 4.2 s per simulated kyr (VILMA2 221 s, output 101 s, plastic 8 s). H against PaleoMIST: bias +283 to +450 m,
volume 1.18 to 1.51 (as the offline check). ESL (`bsl_ive` - `bsl_ive`(0), m):

| ka | 80 | 65 | 60 | 45 | 30 | 22.5 | 20 | 15 | 10 | 5 |
|---|---|---|---|---|---|---|---|---|---|---|
| PaleoMIST | -10.0 | -51.1 | -78.6 | -38.4 | -61.7 | -110.9 | -115.6 | -81.0 | -29.3 | -5.0 |
| I910_it0 (B1 at 910) | -10.2 | -52.6 | -82.4 | -39.2 | -63.8 | -117.8 | -122.8 | -84.5 | -29.5 | -5.1 |
| online plastic | -12.5 | -61.4 | -94.4 | -44.7 | -72.4 | -131.7 | -137.7 | -96.2 | -34.0 | -5.2 |

GAPSLIP 2.1 indicators (`tools/levante/plastic_score.sh`: `pm_compare.jl` with the run's VILMA2 stream beside B1,
the same scoring as docs/dev/vilma2-paleomist.md; RMS z, median |z|, fraction |z| <= 2, mean misfit [m], marine /
terrestrial limiting data violated):

| set | model | RMS z | median abs(z) | abs(z) <= 2 | mean d | ML viol. | TL viol. |
|---|---|---|---|---|---|---|---|
| all | PaleoMIST | 29.5 | 4.42 | 0.30 | 1.5 | 0.45 | 0.30 |
| | B1_l170 | 30.4 | 4.41 | 0.30 | 2.9 | 0.42 | 0.30 |
| | I910_it0 | 30.3 | 4.36 | 0.31 | 3.2 | 0.42 | 0.30 |
| | online plastic | 25.5 | 4.03 | 0.32 | 1.1 | 0.41 | 0.36 |
| Laurentide | B1_l170 | 9.6 | 5.32 | 0.19 | 1.3 | 0.25 | 0.25 |
| | online plastic | 9.8 | 4.73 | 0.28 | 4.3 | 0.19 | 0.31 |
| Greenland | B1_l170 | 110.6 | 9.28 | 0.09 | -15.1 | 0.60 | 0.02 |
| | online plastic | 88.9 | 9.22 | 0.07 | -19.4 | 0.58 | 0.02 |
| Fennoscandia/Barents | B1_l170 | 16.5 | 5.82 | 0.27 | 3.8 | 0.26 | 0.44 |
| | online plastic | 14.6 | 4.82 | 0.28 | 2.9 | 0.23 | 0.52 |
| far field | B1_l170 | 16.0 | 3.23 | 0.39 | 5.1 | 0.44 | 0.18 |
| | online plastic | 15.3 | 3.07 | 0.38 | 1.4 | 0.47 | 0.17 |

The online run scores as well as B1 or better on most sets with 15-20 % more NH ice; the terrestrial limiting
data are violated more often (the sea level too high where the ice is too thick). Differences to B1 besides the
ice: the run's densities (910, 1028; B1 931, 1000), today's GEBCO bed on NH-32KM instead of PaleoMIST's 0 ka bed,
the 100-yr coupling.

PaleoMIST's ice area by threshold (NH-32KM, 1e6 km2; `H > 0`, `> 1 m`, `> 10 m`): 80 ka 5.55, 4.94, 4.60; 65 ka
16.31, 15.66, 15.30; 45 ka 12.05, 11.46, 11.08; 20 ka 28.33, 27.66, 27.23; 0 ka 3.13, 2.66, 2.42; the threshold moves
the plastic volume by about 10 % (flat, 80 ka: 21.6, 19.2, 17.9 m SLE).

**Baseline** (2026-11-09, restored PaleoMIST files; margin `f_ice` > 0.5, `margin_position = "sdist"`, the library
with `margin_frac`). Offline check, `"surface"`, volume [m SLE] by margin variable and position (at the slices `"sdist"` is
`"face"`):

| margin | bias [m] | rms [m] | volume / PaleoMIST | 80 ka | 65 ka | 45 ka | 20 ka | 0 ka |
|---|---|---|---|---|---|---|---|---|
| PaleoMIST | | | | 12.9 | 53.7 | 40.6 | 112.5 | 7.5 |
| `f_ice` > 0.5, sdist | +88 to +166 | 146 to 204 | 1.07 to 1.10 | 13.8 | 58.1 | 43.5 | 122.4 | 8.0 |
| `f_ice` > 0.5, centre | +171 to +244 | 211 to 292 | 1.10 to 1.17 | 14.9 | 60.2 | 45.0 | 124.3 | 8.6 |
| `H_ice` > 1 m, sdist | +225 to +386 | 273 to 440 | 1.14 to 1.39 | 17.4 | 63.4 | 47.6 | 128.3 | 9.8 |
| `H_ice` > 1 m, centre (before) | +257 to +502 | 309 to 558 | 1.16 to 1.50 | 18.7 | 65.6 | 49.2 | 130.0 | 10.5 |

The margin variable takes most of the excess, the position within the cells the rest of the half cell: 7-10 % more
volume than PaleoMIST remain (the native-grid check of the library: 5-7 % at 30 km is the resolution). The misfit is
positive and broad, larger near the margins; the patch at the head of Baffin Bay persists, smaller: at 20 ka 18 cells
over +500 m (max +936 m at 77.5 N 65.8 W, the Greenland side of Smith Sound) against 87 with the old margin
(max +1082 m); 6 cells over +500 m in the Canadian Arctic today (max +894 m, northern Ellesmere).

![H - PaleoMIST, new baseline](figures/plastic-backend/plastic_dH_fice.png)

Online with VILMA2 (`par/kryos_North_plastic_vilma2_b1.nml`, run `runs/vilma2-b1-fice`): 13.6 min on 16 threads (3 runs
on a node). ESL [m]:

| ka | 80 | 65 | 60 | 45 | 30 | 22.5 | 20 | 15 | 10 | 5 |
|---|---|---|---|---|---|---|---|---|---|---|
| PaleoMIST | -10.0 | -51.1 | -78.6 | -38.4 | -61.7 | -110.9 | -115.6 | -81.0 | -29.3 | -5.0 |
| I910_it0 (B1 at 910) | -10.2 | -52.6 | -82.4 | -39.2 | -63.8 | -117.8 | -122.8 | -84.5 | -29.5 | -5.1 |
| online plastic, old margin | -12.5 | -61.4 | -94.4 | -44.7 | -72.4 | -131.7 | -137.7 | -96.2 | -34.0 | -5.2 |
| online plastic, baseline | -10.1 | -55.7 | -88.8 | -41.4 | -67.2 | -125.4 | -131.7 | -92.3 | -32.6 | -5.0 |

GAPSLIP 2.1 (as above):

| set | model | RMS z | median abs(z) | abs(z) <= 2 | mean d | ML viol. | TL viol. |
|---|---|---|---|---|---|---|---|
| all | PaleoMIST | 29.5 | 4.42 | 0.30 | 1.5 | 0.45 | 0.30 |
| | B1_l170 | 30.4 | 4.41 | 0.30 | 2.9 | 0.42 | 0.30 |
| | I910_it0 | 30.3 | 4.36 | 0.31 | 3.2 | 0.42 | 0.30 |
| | online, old margin | 25.5 | 4.03 | 0.32 | 1.1 | 0.41 | 0.36 |
| | online, baseline | 25.9 | 3.95 | 0.34 | 1.8 | 0.39 | 0.33 |
| Laurentide | PaleoMIST | 9.0 | 4.13 | 0.26 | 2.9 | 0.20 | 0.28 |
| | B1_l170 | 9.6 | 5.32 | 0.19 | 1.3 | 0.25 | 0.25 |
| | I910_it0 | 9.4 | 5.12 | 0.20 | 1.3 | 0.25 | 0.24 |
| | online, old margin | 9.8 | 4.73 | 0.28 | 4.3 | 0.19 | 0.31 |
| | online, baseline | 9.4 | 4.72 | 0.28 | 3.0 | 0.21 | 0.28 |
| Greenland | PaleoMIST | 107.1 | 11.12 | 0.04 | -22.2 | 0.69 | 0.01 |
| | B1_l170 | 110.6 | 9.28 | 0.09 | -15.1 | 0.60 | 0.02 |
| | I910_it0 | 110.2 | 9.60 | 0.09 | -15.3 | 0.60 | 0.02 |
| | online, old margin | 88.9 | 9.22 | 0.07 | -19.4 | 0.58 | 0.02 |
| | online, baseline | 91.7 | 10.21 | 0.14 | -12.4 | 0.52 | 0.02 |
| Fennoscandia/Barents | PaleoMIST | 19.9 | 7.80 | 0.17 | 6.4 | 0.32 | 0.41 |
| | B1_l170 | 16.5 | 5.82 | 0.27 | 3.8 | 0.26 | 0.44 |
| | I910_it0 | 16.2 | 5.55 | 0.28 | 3.4 | 0.25 | 0.44 |
| | online, old margin | 14.6 | 4.82 | 0.28 | 2.9 | 0.23 | 0.52 |
| | online, baseline | 13.2 | 4.30 | 0.32 | 1.8 | 0.23 | 0.49 |
| far field | PaleoMIST | 13.6 | 2.85 | 0.41 | 1.3 | 0.46 | 0.16 |
| | B1_l170 | 16.0 | 3.23 | 0.39 | 5.1 | 0.44 | 0.18 |
| | I910_it0 | 16.1 | 3.20 | 0.39 | 6.0 | 0.44 | 0.20 |
| | online, old margin | 15.3 | 3.07 | 0.38 | 1.4 | 0.47 | 0.17 |
| | online, baseline | 15.4 | 3.14 | 0.39 | 3.1 | 0.45 | 0.17 |

The baseline ESL is 14 % above PaleoMIST's at 20 ka (old margin 19 %, I910_it0 6 %); it scores better than B1 on every NH
set but RMS z of Greenland (dominated by a few points, below) and violates fewer terrestrial limiting data than the old
margin.

Checks on the restored files: the baseline with `taub_scale_file = ""` and without the key give the same output, every
file bit for bit; `margin_position = "centre"` with the new library and the old margin (`H_ice` > 1 m) gives the old
binary's ice bit for bit, and VILMA2's sea-level diagnostics (`bsl`, `mass_resid`) and the stress residual to 2e-14.

## Calibration of tau_b (ES-MDA)

Gowan tuned tau_b by hand against RSL and GIA data; here per-group, time-dependent multipliers on his tau_b are calibrated
online, VILMA2 in the loop, against the GAPSLIP 2.1 sea-level index points, by an ensemble smoother with multiple data
assimilation (ES-MDA, Emerick and Reynolds 2013). Tools in `tools/levante`: `plastic_esmda.toml` (the config),
`plastic_taub_scale.jl` (groups, multiplier files), `plastic_esmda.jl` (init, scoring, update), `plastic_esmda.sh` (the jobs).

**Parameters**: theta = the log multipliers of 3 groups x 9 nodes (80, 60, 40, 25, 20, 15, 10, 5, 0 ka) = 27, linear in time
between the nodes (`taub_scale_file`). Groups (`[[scale.groups]]`, by Gowan's region of the domain): Laurentide (region 1 without
Greenland; 101 domains with cells on NH-32KM, 24771 cells), Greenland (109, 3830), Eurasia (region 2; 159, 10075); Antarctica is
not calibrated (no cells on NH-32KM). Greenland: a region-1 domain whose cells on the grid lie more than half in Greenland of the
FesmData REGIONS map (`region_2` = 103, a partition of land and shelf by territory). The rule is clear-cut except for four
domains at the Nares Strait and Baffin Bay (fraction in Greenland: 1056 0.49, 1063 0.27, 1000 0.16, 1055 0.04: Laurentide; the
next highest 0.86 and above: Greenland); `groups.csv` of a calibration lists them all. Nodes without direct data share one
parameter (`tie` of a group): Greenland 80/60 ka (no index points before 44 ka), Eurasia 80/60/40 ka (none before 25 ka in
Fennoscandia/Barents): 24 parameters. Prior: N(0, C), sigma 0.3 (a factor 0.74-1.35 at one sigma), exponential correlation in
time (10 kyr, between the mean ages of the parameters' nodes) within a group, groups independent; every member is clipped to
|theta| <= 1.5 (factor 0.22-4.5) after an update (the count is logged).

**Observations**: the index points of GAPSLIP 2.1 (primary, 0-80 ka) of pm_compare.jl's groups Laurentide (896),
Greenland (276), Fennoscandia/Barents (1088) and the far field (2445): 4705. The model RSL at a point is the scorer's
(`pm_kryos_out.jl` makes the run's VILMA2 stream a `vilma_out.nc`, `pm_compare.jl`'s `read_run` and `series_run`, RSL relative
to 0 ka, bilinear in space, linear in time between the 2.5-kyr records): the residual d = model - rsl is the d of the scores.
sigma is the scorer's `ip_sigma`: half the 2-sigma elevation range on the side of d, and the age uncertainty times the model
slope, evaluated at the ensemble-mean RSL of the iteration and fixed within it, and `sigma_add` = 5 m in quadrature: the model
error the indicator's sigma lacks (GIA and Earth model, resolution, RSL interpolation). A gross-error cut, fixed in iteration
1 at the prior mean member (Gowan's tau_b): index points with |d| / sqrt(sigma^2 + sigma_add^2) > 5 there are not assimilated
(`obs_kept.csv`); every score stays over all index points, comparable with B1. Index points without a model value in some
member are left out of that iteration. The limiting data
(marine 2436, terrestrial 1353 in those groups) are not assimilated; their violation fractions are reported.

**Update**: Ne members (40) and N_a iterations with inflation alpha (9.333, 7, 4, 2); perturbed observations
y + sqrt(alpha) sigma e; in units of sigma, theta_a = theta + dTheta (I + G S)^-1 G (y_pert - D) with S the scaled
anomalies of the predictions, dTheta those of theta and G = S'/alpha (Woodbury: the Kalman gain with an Ne x Ne solve; no
localisation). Each iteration also runs the mean member (m000: the prior mean, Gowan's tau_b, in iteration 1; the ensemble
mean after an update); after the last update the posterior mean is run (`final = "mean"`; or the ensemble, or nothing). A
member that fails (exit code, missing output) is dropped and logged; fewer than `min_members` (half) stops the calibration.

**Jobs**: `plastic_esmda.sh start <config> <folder>` runs the init as a job (the prior ensemble, the multiplier files and par
files of iteration 1, a frozen copy of `kryos.x`, `input/` and the par file with the `set` overrides) and submits iteration 1:
the members as a job array (one kryos run each, `-c threads`), then the analysis after them (`afterany`), which scores,
updates and submits the next iteration until `next` says done. The folder: `esmda.md` (the report: per iteration the mean
member's scores and the ensemble's, chi2, ESL at 20 ka, theta before and after as factors), `iter_kk/theta.csv`,
`scores.csv` (per member and set, the scorer's metrics), `members.csv` (ok, chi2, ESL, run time), `mNNN/` (the run).
Members write only `geo_glob_vilma.nc` (2.5 kyr) and the 1-D series (no restarts, no 2-D ice output).

**Checks** (Levante, 2026-11-09, the restored PaleoMIST files, the old margin and library): kryos before the key, with the key
unset, and with a multiplier file of zeros (members' outputs only) give the same `geo_glob_vilma.nc` and 1-D series bit for bit,
also at 4, 8, 16 and 64 threads, and the same as the earlier online run (all output files); `pm_compare.jl` after the
refactor scores that run as before. `plastic_esmda.jl selftest`: four ES-MDA steps (alpha 9.333, 7, 4, 2, N = 4000) on a
linear-Gaussian problem give the exact posterior mean to 0.04 sd and its sd to 1.3 %.

**Cost**: VILMA2 dominates and barely scales: a member takes 12.5 min on 8 threads, 15 on 4, 12-14 on 16, 8 on 64 (several
members on one node; 8 min on 16 threads alone), 2-7 GB. The shared partition bills max(cpus, mem / 0.9 GB): 8 threads with
10 GB (12 cpus billed) is the cheapest per member, 2.5 cpu-h, against 9 cpu-h on 64. The full calibration, 33 x 4 runs and the
posterior mean, is about 330 cpu-h (2.6 node-h), about 1-1.5 h of wall time if the array runs at once.

**Smoke test** (`runs/esmda-smoke`: Ne = 8, alpha = 1, then the posterior mean; margin `f_ice` > 0.5, the build before `margin_position`,
linked while the library was being rebuilt, so its margin position is not certain; its mean member scores within 2 % of
the baseline): the 9 members ran at once in 12-15 min each on 8 threads, all succeeded;
scored and updated in 2 min. Prior ensemble (factors 0.43-2.3 across the parameters): RMS z of all index points 18.0 to
31.4 (mean member 25.8), Greenland 48 to 115, ESL at 20 ka -107 to -148 m (mean member -132): tau_b moves the scores
well beyond the spread between B1 and the runs. With the indicators' own sigma, chi2/n is about 740 at the mean member:
1089 of the 4705 index points have |z| > 10 and carry 98 % of the sum of z^2 (104 have |z| > 50: southern Greenland,
Townsville and Bowen, the eastern Baltic). The one update with alpha = 1 and Ne = 8 < 27 parameters collapsed the ensemble
(sd ratio 0.03) and sent the mean to factors 0.0004-895; its run did not get through VILMA2's relaxation to the 80 ka load
in 26 min and was cancelled; the analysis then logged the failure and stopped (`next`: done). The same members, updated
offline with other settings (median |shift| of theta; factor range; sd after / before):

| alpha | sigma_add | index points | median shift | factors | sd ratio |
|---|---|---|---|---|---|
| 1 | 0 | all | 1.73 | 0.00-914 | 0.03 |
| 9.333 | 0 | all | 1.72 | 0.00-508 | 0.07 |
| 9.333 | 5 m | all | 0.39 | 0.42-2.7 | 0.23 |
| 9.333 | 10 m | all | 0.23 | 0.53-2.0 | 0.37 |
| 9.333 | 0 | abs(z) <= 10 at the mean | 0.11 | 0.66-1.34 | 0.13 |
| 9.333 | 5 m | abs(z) <= 10 at the mean | 0.14 | 0.68-1.51 | 0.35 |

So the indicators' sigma alone is far below the model error: the full calibration needs a model error (`sigma_add`,
5-10 m) and/or a gross-error check (|z| at the prior mean), else a few hundred points drive it.

Launch:

```
bash tools/levante/plastic_esmda.sh start tools/levante/plastic_esmda.toml /work/ba1442/robinson/kryos-ice-plastic/runs/esmda-<name>
```

## Calibration v1 (2026-11-09)

`runs/esmda-v1` (`plastic_esmda.toml` at fa734aa: the baseline above, margin `f_ice` > 0.5 at `"sdist"`, Gowan's 1-D Earth;
24 parameters, Ne = 40, alpha 9.333, 7, 4, 2, sigma_add 5 m, cut at |z| > 5, at most 29 members at once), reported by
`tools/levante/plastic_esmda_report.jl`. The cut kept 4385 of the 4705 index points: Laurentide 836 of 896, Greenland 150 of
276, Fennoscandia/Barents 1065 of 1088, far field 2334 of 2445. All 165 runs succeeded (11-24 min on 8 threads); 3.0 h of wall
time, 464 cpu-h with the analyses (12 cpus billed per member). Scores of the mean member over all index points (RMS z; marine
/ terrestrial limiting data violated), the ensemble's RMS z (median [min, max]) and ESL at 20 ka:

| iteration | all | Laurentide | Greenland | Fennoscandia/Barents | far field | ML / TL viol. | ensemble RMS z, all | ESL 20 ka [m] (ensemble) | sd(theta) |
|---|---|---|---|---|---|---|---|---|---|
| 1 (prior, Gowan) | 25.97 | 9.41 | 91.7 | 13.15 | 15.44 | 0.37 / 0.35 | 27.6 [14.8, 38.4] | -131.7 (-171 to -106) | 0.30 |
| 2 | 24.90 | 9.25 | 87.6 | 11.00 | 15.53 | 0.35 / 0.33 | 25.0 [22.2, 35.8] | -127.9 (-132 to -124) | 0.16 |
| 3 | 24.64 | 9.25 | 86.4 | 10.63 | 15.60 | 0.35 / 0.31 | 24.7 [22.3, 26.0] | -127.5 (-130 to -126) | 0.14 |
| 4 | 24.16 | 9.22 | 84.3 | 10.30 | 15.61 | 0.35 / 0.31 | 24.2 [22.3, 25.0] | -127.1 (-129 to -126) | 0.13 |
| posterior mean | 23.69 | 9.19 | 81.9 | 10.20 | 15.69 | 0.34 / 0.32 | | -127.4 | |

(PaleoMIST: all 29.5, B1_l170 30.4; ESL PaleoMIST -115.6, I910_it0 -122.8.) chi2/n of the assimilated points (their sigma)
fell from 2.06 to 1.67. No member failed and no update diverged; the spread shrank steadily. Clipping: none in the first two
updates, 6 values in the third, 40 in the last: Eurasia 10 ka at the bound in every member.

![posterior tau_b factors](figures/plastic-backend/esmda_v1_factors.png)

Posterior-mean factors (grey: prior +-1 sigma; red band: posterior ensemble +-1 sd; thin lines: the means after updates 1-3):

| group | 80 ka | 60 ka | 40 ka | 25 ka | 20 ka | 15 ka | 10 ka | 5 ka | 0 ka |
|---|---|---|---|---|---|---|---|---|---|
| Laurentide | 1.52 | 0.47 | 0.83 | 1.29 | 0.91 | 0.97 | 1.24 | 1.08 | 2.95 |
| Greenland | 1.27 | 1.27 | 0.98 | 0.94 | 0.79 | 2.02 | 0.36 | 0.27 | 0.44 |
| Eurasia | 1.41 | 1.41 | 1.41 | 1.30 | 0.84 | 0.39 | 4.48 | 2.89 | 0.91 |

NH ice volume (m SLE, NH-32KM; the coupled runs' `M_ice`): PaleoMIST 112.5 at 20 ka, 7.45 today; Gowan's tau_b 124.7 (1.11),
8.03 (1.08); the posterior mean 117.6 (1.05), 5.53 (0.74). Offline (`plastic_slices` with the posterior multipliers, PaleoMIST's
own bed): 115.2 (1.02) at 20 ka, 5.69 (0.76) today, 16.2 (1.26) at 80 ka, 75.9 (0.92) at 60 ka.

Greenland today, offline, over the Greenland region: PaleoMIST 7.24 m SLE, BedMachine v6 7.23; Gowan's tau_b 7.72 (H - BedMachine
+134 m mean, 203 m rms); the posterior mean 5.18 (-485 m mean, 576 m rms). The calibration thins today's Greenland ice sheet by
a quarter: it fails this sanity check.

Reading: the LGM and the deglaciation of the Laurentide are well constrained and barely move (0.9-1.3, sd 0.02-0.05); the
gain is in Fennoscandia/Barents (RMS z 13.2 -> 10.2) and Greenland (91.7 -> 81.9), little in the far field. The posterior
has features that are not physical: alternating factors between neighbouring deglacial nodes (Eurasia 0.39 at 15 ka, 4.48
at 10 ka; Greenland 2.0 at 15 ka, 0.27-0.44 from 10 ka), a factor 2.95 on the small Laurentide (Canadian Arctic) ice today,
and a Greenland ice sheet 28 % too thin today. tau_b likely absorbs misfits it should not explain: the 1-D Earth (one viscosity
profile for every region), the timing of the margin chronology, and the index points of one site counted as independent (the
cut removed 126 of Greenland's 276). The posterior is a fit, not yet a reconstruction.

## Open

- The NH ice is 7-10 % too large with Gowan's tau_b and the new margin (was 15-50 %): mostly resolution (above).
- Calibration v1 (above): oscillating deglacial factors, Eurasia 10 ka at the bound, Greenland too thin today. Options: a
  smoother prior in time (longer correlation, or fewer deglacial nodes), Greenland's present-day volume (BedMachine) as an
  observation or Greenland not calibrated, sites as units (correlated index points), a 3-D Earth (the science Earth).
- `drive = "smb"` does not converge on 32 km grids (library); `"full"` neither (use `"surface"`).
- The science-Earth variant (`par/kryos_North_plastic_vilma2.nml`, 3-D Bagge, lmax 170, refcorr reference) is not
  run with `"surface"` (its flat run was cancelled in VILMA2's spin-up).
