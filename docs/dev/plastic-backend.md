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
  order, an id added by the table, the stops).
- `tests/test_ice_plastic.f90` (`plastic=1`): the backend on a toy grid (flat bed, two disks, tau_b const,
  `mode = "flat"`): the plastic profile at the centre, the masks, the static exports, an update (dHidt,
  mb_net, M_ice), a restart round trip stepped on bit for bit, the flotation binding.
- `tests/test_build.f90`: the factory admits `plastic` when linked.

## Results (Levante, 2026-11-09)

Build: kryos `plastic=1 yelmo=0 vilma=1` (stack S11 dependencies) against PlasticIceSheet `fortran-port`
(`/work/ba1442/robinson/PlasticIceSheet-fortran/ifx`); runs in `/work/ba1442/robinson/kryos-ice-plastic/runs`.
tau_b: `NH-32KM_IceHistory-PaleoMIST_taub.nc` (`/work/ba1442/robinson/fesmdata-taub`, not yet in FesmData),
domains with the time series, `taub_default` 5000 Pa; the margin `H_ice` > 1 m (the default is now `f_ice` > 0.5).

**Fixed bed** (`par/kryos_North_plastic.nml`: today's bed, z_sl = 0, 80 ka -> 0, dtt = 100 yr, 800 solves): 25 s on
4 cores; `"surface"` 10 ms, `"flat"` 4 ms per solve; every solve converged; the ice mask is PaleoMIST's (H > 1 m) at
every slice. Restart (flat, tau_b const): the 40 ka bundle continued to 0 equals the continuous run bit for bit
(1284 variables and records).

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

## Open

- The NH ice is 15-50 % too large with Gowan's tau_b although `"surface"` solves the plastic equation: tau_b of
  the product, resolution or margin (above).
- `drive = "smb"` does not converge on 32 km grids (library); `"full"` neither (use `"surface"`).
- The science-Earth variant (`par/kryos_North_plastic_vilma2.nml`, 3-D Bagge, lmax 170, refcorr reference) is not
  run with `"surface"` (its flat run was cancelled in VILMA2's spin-up).
