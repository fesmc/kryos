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
| `par/kryos_North_plastic_vilma2.nml` | NH-32KM with the global VILMA2 Earth of the science case |
| `tests/test_plastic_margin.f90` | margin and tau_b (runs with `plastic=0`) |
| `tests/test_ice_plastic.f90` | the backend through the ice wrappers (`plastic=1` only) |

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
4. `plastic_solve(ps, z_bed, z_sl, tau_b, mask, warm=.true., stat)`, warm-started from the last
   `w = H**2` (the library's rule; no effect in `"flat"` mode) (`drive = "smb"`: `plastic_solve_smb` with the
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
`converged`, `n_outer_used` of the last solve.

## Margin chronology

The slices of `margin_var` (default `H_ice`) in `margin_file` (`{domain}`, `{grid_name}` replaced, as
`recon_path`), on the instance grid, with ice where the thickness exceeds `margin_H_min` (a fill value is no
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
| `mode` | `"flat"` | passed to the library as it is: `"flat"`, `"full"` (bed slope in the force balance; does not converge on 32 km beds, GRL-32KM, in Julia too), later `"surface"` |
| `relax`, `warm_frac_max`, `n_outer`, `max_sweeps`, `tol`, `outer_tol` | 0.5, 0.1, 100, 200, 1e-6, 1e-3 | the library's iteration controls |
| `rho_i` | 910 | [kg m-3] the kryos `Earth` constant and the PaleoMIST runs (library default 917) |
| `rho_w`, `g` | 1028, 9.81 | the kryos `Earth` constants (library defaults 1025, 9.80665) |
| `H_min` | 1 | [m] floor of H inside the mask |
| `margin_file`, `margin_var`, `margin_H_min` | `""`, `"H_ice"`, 1 | the margin slices |
| `margin_interp`, `margin_extrapolate` | `"sdist"`, `"error"` | see above |
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

## Smoke runs (Levante, 2026-11-09)

Build: kryos `plastic=1 yelmo=0 vilma=1` (stack S11 dependencies) against PlasticIceSheet `fortran-port`
(`/work/ba1442/robinson/PlasticIceSheet-fortran/ifx`); runs in `/work/ba1442/robinson/kryos-ice-plastic/runs`.
tau_b: the product `NH-32KM_IceHistory-PaleoMIST_taub.nc` (`/work/ba1442/robinson/fesmdata-taub`, not yet in
FesmData), domains with the time series, `taub_default` 5000 Pa. `mode = "flat"`.

`par/kryos_North_plastic.nml` (NH-32KM, today's bed, z_sl = 0, 80 ka -> 0, dtt = 100 yr, 800 solves): 22 s in
all on 4 cores, the solve about 0.004 s per step (`timing.txt`); every solve converged. The ice mask is
PaleoMIST's (H > 1 m) at every slice time. Restart: 0 -> 40 ka bundle, continued to 0 (tau_b const): every
record of `ice_nor_plastic.nc`, `ice_nor.nc`, `ice_nor_ts.nc` after the restart and every file of the final
bundle equal to the continuous run (1284 variables and records, bit for bit).

Against PaleoMIST's `H_ice` at its 33 slices, over its ice:

| Setup | bias [m] | rms [m] | volume / PaleoMIST |
|---|---|---|---|
| kryos run, today's bed, tau_b uniform 1e5 Pa | +640 to +1420 | 750 to 1550 | 1.59 to 2.10 |
| kryos run, today's bed, tau_b uniform 3.5e4 Pa | -68 to +218 | 475 to 660 | 0.96 to 1.26 |
| kryos run, today's bed, Gowan tau_b | +270 to +490 | 410 to 615 | 1.18 to 1.55 |
| `plastic_slices`, PaleoMIST's bed of the slice, Gowan tau_b | +285 to +550 | 425 to 650 | 1.18 to 1.55 (m SLE) |
| the same, `mode = "full"` (100 outer iterations, none converged) | -196 to +306 | 355 to 580 | 0.89 to 1.29 (m SLE) |

The offline check (`tests/plastic_slices.f90`: at each slice PaleoMIST's own bed, z_sl = 0, its margin, Gowan's
tau_b of that time, cold flat solve) does not reproduce PaleoMIST: the thickness is 300 to 500 m too high and the
NH volume 20 to 50 % too large (20 ka: 132.6 against 112.5 m SLE; 0 ka: 10.6 against 7.5 m). The bed of the time
changes little (the kryos run on today's bed gives the same misfit), so the difference is in the solve or in
tau_b (the flat mode ignores the bed slope; Gowan's ICESHEET integrates the surface along flowlines with the
bed). The `"full"` mode, not converged but usable, removes most of the bias (20 ka: 100.1 against 112.5 m SLE;
62.5 ka: 74.0 against 73.7; 0 ka: 8.6 against 7.5): the bed slope is the missing term. To be redone with the
library's direct `"surface"` mode.

PaleoMIST's ice area by threshold (NH-32KM, 1e6 km2; `H > 0`, `> 1 m`, `> 10 m`): 80 ka 5.55, 4.94, 4.60; 65 ka
16.31, 15.66, 15.30; 45 ka 12.05, 11.46, 11.08; 20 ka 28.33, 27.66, 27.23; 0 ka 3.13, 2.66, 2.42. The rims below
1 m are 2 to 15 % of the area; the margin threshold moves the plastic volume by about 10 % (80 ka: 21.6, 19.2,
17.9 m SLE for 0, 1, 10 m).
## Open

- `mode = "full"` and `drive = "smb"` do not converge on 32 km grids (library, Julia too): the runs use
  `"flat"` and `"taub"`; a direct `"surface"` mode of the library needs no backend change.
