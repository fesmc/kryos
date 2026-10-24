# kryos phase 1: backend adapters

Status: design note, 2026-10-24, no code. Scope: the eight adapters phase 1
needs to reproduce `yelmox_Greenland` and `yelmox_Antarctica`: `ice = yelmo`,
`geo = fastisostasy` (+ its global `bsl` part), `climate = snapshots`,
`surface = chion`, `shelf = marine_shelf` (incl. `bmb_method = "pico"`),
`geothermal = data`, `sediments = data`, `forcing = tsgen`.

Line references: yelmox `62bc7a6` (`libs/`, `yelmox/`), yelmo `f4f623f0`
(`src/`), FastIsostasy `22a458f` (`src/`), fesm-utils `phase0-kryos`
(`src/`). `design.md` decisions are cited as D1..D13.

## 0. What this note assumes from the parallel notes

`core-design.md` and `contracts-design.md` are being written now. This note
assumes the following; each line is a requirement on those notes, and a
mismatch is resolved there, not by changing an adapter silently.

| # | Assumption |
|---|---|
| A1 | Canonical names and units as in design §5 and `study/backends.md` §A; mass fluxes `kg m-2 yr-1`, + = ice gains (D3). Names used here that the table does not fix: `T_ocn`, `S_ocn`, `dT_ocn` (ocean profiles on axis `depth`), `lapse` (axis `month`), `m_ice` (scalar, Gt), `dbsl_vol` (scalar, m3). |
| A2 | `init(self, par_file, group, grid, time, cnst, ref)`: `ref` is a plain `kryos_ref_t` with the domain reference (`z_bed`, `H_ice`, `z_srf`, `z_bed_sd`, `dp`) and the `regions_class` on the instance grid, made by the core with the same conservative map for every grid (D4). Global instances get an empty `ref`. |
| A3 | `update(self, tm)`: `tm%t0`, `tm%t1` model time [yr, dp] from `cadence_interval`, `tm%t_rel` forcing time at `t1` (`clock_time_rel`). Backends never see ticks (D8). |
| A4 | `init_state(self, tm, stage)` is called by the core per stage, in the order of section 9; a backend ignores the stages it does not use. Stages: `INIT_GEOMETRY`, `INIT_FIELDS`, `INIT_STATE`, `INIT_FINAL`, `INIT_RESTART`. This generalises `domain_init_state` (`kryos_startup.f90:131-195`) and the tail of `domain_restart_read` (`:521-527`). |
| A5 | The core writes and restores every instance's export state with `fields_restart_write/read`, so a published export survives a restart without the adapter. Adapter files hold only backend state. |
| A6 | The re-elevation of the atmosphere (`remap_atm_to_surface`, `kryos_coupling.f90:508-539`) is a transform of the core's `climate -> surface` wiring row, not code in a surface adapter (every surface backend needs it). |
| A7 | Wiring rows phase 1 needs for bit-identity: `ice.H_ice -> geo.H_ice` from the ice exports, not from the hub (today `step_isostasy` reads Yelmo, `kryos_coupling.f90:146`); the hub mirror takes `ice.z_bed`, `ice.z_sl` as exported (widened from `bnd%`), not the core's own `dp` sums; `shelf.tf_shlf` is not wired to `ice` (Yelmo never receives it today). |
| A8 | A rank-3 export with `zdim = "depth"` carries its axis values through a core axis table filled at init by the producer. `fesm_fields` has the axis name but no coordinate values; this is a gap to close in contracts or fesm-utils. |

## 1. Cross-cutting rules

### 1.1 Shape of one adapter file

```fortran
type, extends(ice_t) :: yelmo_ice_t       ! src/backends/ice/kryos_ice_yelmo.f90
    type(yelmo_class) :: ylmo             ! the backend object
    type(yelmo_ini_t) :: ini              ! adapter parameters (instance group)
contains                                  ! the only bindings (D1 style rule)
    procedure :: init => init, init_state => init_state, update => update
    procedure :: restart_write => restart_write, restart_read => restart_read
    procedure :: write_output => write_output, end => finish
end type
! bound bodies stay short and call private free procedures:
! bnd_from_imports(self), exports_from_tpo(self, tm), ...
```

The stub `kryos_ice_yelmo_stub.f90` keeps the module and type names, holds no
backend object, and every binding calls `unavailable(self)` ("rebuild kryos
with `yelmo=1`"), as the skeleton already does. Data adapters
(`geothermal`, `sediments`) have no external library and no stub; they live
with their component (design §13).

Instance keys (backend, grid, dt, names of physics groups, adapter options)
are validated against one template group per backend, `[<component>_<backend>]`
in `input/kryos_defaults.nml`; backend physics groups against the library's
own defaults file (`yelmo_defaults.nml`, `chion_defaults.nml`, ...;
synchronised copies, design §9).

### 1.2 Precision: the bit-identity argument

Every phase 1 backend computes in `sp`; the exchange is `dp` (D2). Five
facts carry the bit-identity argument through this note:

- **P1 Widen and narrow.** `real(x_sp, dp)` is exact; `real(real(x_sp, dp), sp) == x_sp`.
  `fields_put`/`fields_get` have `sp` specifics that do exactly this
  (`fesm_fields.f90:411-505, 860-866`).
- **P2 One operation.** For `sp` operands `a`, `b` and `op` in `+ - * /`,
  `real(real(a,dp) op real(b,dp), sp) == a op b` evaluated in `sp`: double
  rounding is innocuous since 53 >= 2*24 + 2. It does not extend to two
  chained `dp` operations before the narrowing.
- **P3 Round trip.** For `sp` `x` and a `dp` factor `c`, `real((real(x,dp)*c)/c, sp) == x`:
  the `dp` error (a few 2^-53) is far below half an `sp` ulp (2^-25), and
  `x` is itself an `sp` number, never a tie. Holds when nothing but the
  identity remap lies between the two conversions.
- **P4 Remap.** The `sp` coupler remap widens, maps in `dp` and narrows
  (`coords/mapping.f90:1301-1354`, "accumulate in dp, store back in sp").
  So the `dp` remap of widened `sp` values, narrowed by the consumer, equals
  today's `sp` remap bit for bit. Identity grids are a copy in both.
- **P5 Interval mean.** A `MEAN` field with one `put` per window publishes
  `(x*dt)/dt`; narrowed, this is `x` (P3).

Consequences, used per field below: an adapter exports the backend's `sp`
value widened (no `dp` arithmetic first); a consumer adapter narrows right
after its one conversion; a unit conversion and its inverse use the same
`dp` factor. Where a non-identity remap lies between a conversion and its
inverse (Antarctica: `shelf` ANT-16KM -> `ice` ANT-32KM), P3 does not hold
strictly: the result differs from today's only when today's value lies within
a few `dp` ulps of an `sp` rounding midpoint (probability ~2^-27 per value).

### 1.3 Kinds, ordering, constants

| Library | `wp` | Array order | Constants | Year length |
|---|---|---|---|---|
| Yelmo | `sp` (`yelmo_defs.f90:42`) | (nx,ny), aa | `cnst` at `yelmo_init` -> `bnd%c`; `sec_year` from Yelmo's own group (`yelmo_boundaries.f90:117`) | its `sec_year`, logged |
| FastIsostasy | `sp` (own `precision` module, name clash hidden by `-zmuldefs`, design §12 row 9) | (nx,ny), padded inside (`out2in`) | `cnst`: rho_ice, rho_w, rho_sw, rho_asth, g (`fastisostasy.f90:118-125`); rho_litho, E, nu, r_earth, m_earth from `[isos]` | none |
| barysealevel | `sp` (`bsl_now`, `A_ocean_now`) | scalar | `A_ocean_pd` from `[barysealevel]` | none |
| chion + host | `sp` (accumulators `dp`) | columns `icol = i+(j-1)*nx` inside the host (`surface_chion.f90:188`) | `cnst` at `chion_init` | host: 360-day cycle of 86400 s days (`:72`); `days_trop = 365.2422` for the solar longitude (`:77`) |
| marine_shelf, pico | `sp` (own, `marine_shelf.f90:24`) | (nx,ny) | `cnst` at `marshelf_init` | none |
| snapshots, climate_terms/mods | fesm-utils `wp = sp` | (nx,ny,12), (nx,ny,nz) | none | none |
| tsgen | `sp` values, `dp` time | scalar | none | none |

No transposition in phase 1. Every phase 1 library takes the run's
`phys_const_class`, so `phys_const_compare` is not needed; adapters log at
init: Yelmo's `sec_year` against the run calendar, `A_ocean_pd` (3.625e14 m2)
against `cnst%area_seasurf` (3.618e8 km2), and chion's 360-day totals used as
per-year rates (1.4 % below a 365.2422-day year; today's behaviour). The year
length snapshots and chion need for `pr` is the run calendar's (the record
has none by design) and cancels (P3).

### 1.4 Restart and output conventions

- Restart: `restart-<time>/<domain>/<instance>_restart.nc` for regional,
  `restart-<time>/<instance>_restart.nc` for global instances (design §9).
  Each adapter writes its backend's native restart format under that name;
  the export state is the core's (A5). A backend that cannot restart from
  its file stops in `restart_read` with the instance and file named; no
  re-derivation by stepping (today's "rebuilt by an annual update" warning in
  `surface_chion.f90:462-466` becomes a stop).
- Output: `<instance>.nc` (2D, native grid), `<instance>_ts.nc` (1D). An
  adapter overrides `write_output` when the backend's diagnostics are richer
  than its exports, with today's variables, units and precision (the writer
  moves from `kryos_output.f90` into the adapter).

## 2. `ice = yelmo`

**Type.** `yelmo_ice_t` extends `ice_t`: `ylmo` (`yelmo_class`), `ini`
(`init_method`, `init_equil_time`, `init_time_thrm`, `vars_2D(:)`,
`vars_2D_sm(:)`), `ref_ok` (logical, section 2.5). Constraints at
registration: none (any grid, several instances, restart). Capabilities: none
in phase 1 (`own_bmb` deferred, D6).

**Groups.**

```
&ice_grl  backend = "yelmo"  grid = "GRL-16KM"  dt = 10.0       ! dt = today's dtt
          yelmo = "yelmo"                 ! Yelmo's top group; its nml_* keys name ytopo, ydyn, ...
          init_method = "none"  init_equil_time = 10.0  init_time_thrm = 0.0 /   ! Antarctica: "equil"
```

`[yelmo]`, `[ytopo]`, `[ycalv]`, `[ydyn]`, `[ytill]`, `[yhyd]`, `[ymat]`,
`[ytherm]`, `[yelmo_masks]`, `[yelmo_init_topo]`, `[yelmo_data]` unchanged
(schema `input/yelmo_defaults.nml`). `[yelmo] phys_const` is ignored: the run
passes `cnst` (the duplicate in `[domain]` and `[yelmo]` goes, study §4.2).

### 2.1 Imports and exports

| Import | Units | Into Yelmo | Conversion (exactly) |
|---|---|---|---|
| `z_bed` | m | `bnd%z_bed` | `real(z_bed, sp)` |
| `z_sl` | m | `bnd%z_sl` | `real(z_sl, sp)` |
| `smb` | kg m-2 yr-1 | `bnd%smb` | `real(smb, sp) * bnd%c%conv_we_ie * 1e-3` (same `sp` expression and order as `kryos_coupling.f90:321`; exact because kg m-2 is numerically mm w.e. with `rho_w = 1000`) |
| `T_srf` | K | `bnd%T_srf` | `real(T_srf, sp)` |
| `bmb_shlf` | kg m-2 yr-1 | `bnd%bmb_shlf` | `real(bmb_shlf / real(bnd%c%rho_ice, dp), sp)` (inverse of section 6, same `cnst%rho_ice`) |
| `T_shlf` | K | `bnd%T_shlf` | `real(T_shlf, sp)` |
| `Q_geo` (opt) | W m-2 | `bnd%Q_geo` | `real(Q_geo * 1000.0_dp, sp)` (mW m-2, + into ice) |
| `H_sed` (opt) | m | `bnd%H_sed` | `real(H_sed, sp)`; unwired: Yelmo's default 0 (`yelmo_boundaries.f90:724`), today's `method = 0` |
| `Qd`, `tf_shlf`, `fmb_shlf` (opt) | | | declared, unwired in phase 1 (snapshots has no `has_Qd`; Yelmo never got `tf_shlf`, A7) |

| Export | Semantics | From | Conversion |
|---|---|---|---|
| `H_ice`, `z_srf`, `f_ice`, `f_grnd` | INSTANT t1 | `tpo%now%...` | widen |
| `z_bed`, `z_sl` | INSTANT t1 | `bnd%z_bed`, `bnd%z_sl` as used | widen (the hub mirror needs these, A7) |
| `z_bed_ref`, `H_ice_ref` | INSTANT, at init | `bnd%z_bed_ref`, `bnd%H_ice_ref` | widen; for the hub anomaly (`kryos_coupling.f90:565`) and section 2.5 |
| `dHidt`, `smb`, `bmb`, `cmb`, `mb_net` | MEAN [t0,t1] | `tpo%now%...` (Yelmo's mean over the call) | `real(x, dp) * real(bnd%c%rho_ice, dp)` (Yelmo's own density) |
| `ux_s`, `uy_s` (opt) | INSTANT | `dyn%now%ux_s/uy_s` averaged to aa | not used in phase 1 |
| `m_ice` (scalar) | INSTANT t1 | `reg%V_ice * bnd%c%rho_ice * 1e-3` in `sp` | widen; the forcing feedback variable [Gt], today `yelmox.f90:185` |

The `bmb` export keeps Yelmo's sign (melt < 0 is "+ = gain" already).

### 2.2 `init`

1. `yelmo_init_grid(ylmo%grd, grid)`; `ylmo%outfldr` and the region file
   names `<instance>_ts[_<region>].nc` (from `domain_regions_init`,
   `kryos_regions.f90:20-46`, before the first `yelmo_update`).
2. `topo` = narrowed `ref%z_bed`, `H_ice`, `z_srf`, `z_bed_sd` (A2): today's
   `topo_y` (`kryos.f90:303-306`) by P4, provided the hub keeps the file's
   `float` values unmodified (FesmData fields are `float`; the core checks).
3. `yelmo_init(ylmo, filename=par_file, grid_def="none", time, domain,
   grid_name, group=ini%yelmo, reg=ref%reg, topo_init=topo, topo_pd=topo,
   cnst=cnst)` (`kryos.f90:310-313`).
4. Declare and lock imports/exports; publish `z_bed_ref`, `H_ice_ref` and the
   geometry at `t_init` (geo's init needs `H_ice`, `z_bed`).

### 2.3 `init_state` (two-step init, A4)

- `INIT_GEOMETRY` (after geo's init, today `domain_init_isostasy` tail,
  `kryos_startup.f90:243-249`): `bnd%z_bed`, `bnd%z_sl` from imports;
  `calc_ytopo_diagnostic(tpo, dyn, mat, thrm, bnd)`; publish geometry exports.
- `INIT_STATE` (after climate, surface, shelf): imports into `bnd` (2.1);
  `yelmo_print_bound`; `yelmo_init_state(ylmo, t, thrm_method="robin-cold")`
  (`kryos_startup.f90:178-179`); then `domain_init_ice` (`:252-287`) for the
  base configs only: `equil` -> `yelmo_update_equil(ylmo, t, init_equil_time,
  dt=1.0, topo_fixed=.false.)`; `init_time_thrm > 0` -> `yelmo_update_equil(...,
  dt=<ice cadence>, topo_fixed=.true.)` (today `dtt`, equal in both configs).
  Publish exports.
- `INIT_RESTART` (after `restart_read` and the geo import): as
  `INIT_GEOMETRY`, plus `yelmo_regions_update(ylmo)` (`kryos_startup.f90:521-527`).

### 2.4 `update(tm)`

1. `bnd_from_imports(self)`: every wired import, converted per 2.1. Static
   imports (`Q_geo`, `H_sed`) are rewritten only when their publication time
   changed (once).
2. `yelmo_update(ylmo, tm%t1)` (`kryos_coupling.f90:386`).
3. `exports_from_tpo(self, tm)`: INSTANT at `t1`, MEAN over `[t0,t1]` (one
   put, P5).

Bed and sea level: today `bnd%z_bed = bnd%z_bed_ref + dz_bed_y`, `bnd%z_sl =
isos%now%bsl + dz_ss_y` in `sp` (`kryos_coupling.f90:225-226`); the core's
`dp` sums of widened `sp` terms over the identity remap (`grid_isos ==
grid_ice` in both configs) narrow to the same value (P2), **provided** the
core's `z_bed_ref` on the ice grid is Yelmo's `bnd%z_bed_ref` bit for bit (2.5).

### 2.5 The reference check

Yelmo processes the reference it is given: `ydata_load` fills gaps, adds
`z_bed_f_sd*z_bed_sd`, removes englacial lakes, zeroes `H_ice < 1 m`, limits
bed gradients (`yelmo_data.f90:254-289`); `ydata_mask_pd_topo` removes ice
where `mask_ice == MASK_ICE_NONE` (`:421-424`); `bnd%z_bed_ref`,
`bnd%H_ice_ref` are the result (`yelmo_ice.f90:1243-1245`). In `init` the
adapter compares both with the narrowed `ref` and logs count and maximum of
differing cells. A bed difference stops init (the core's `z_bed` would not be
Yelmo's); an ice difference is Q1 (today's isostasy reference is Yelmo's
processed `H_ice_ref`).

### 2.6 Restart, output, end

- `restart_write`: `yelmo_restart_write(ylmo, file, time)` (today
  `yelmo_restart.nc`, `kryos_startup.f90:427`).
- `restart_read`: `yelmo_restart_read_topo_bnd`, `yelmo_restart_init`,
  `par%use_restart = .true.` (`:493-497`); publish the geometry (geo's
  restart needs `H_ice`).
- `write_output(stream, time)`: `2D` -> `yelmo_write_step` to `<instance>.nc`
  (`vars_2D`, default `YELMO_VARS_2D`, `kryos_output.f90:31-35`); `2D_sm` ->
  `<instance>_sm.nc` (`YELMO_VARS_2D_SM`); `1D` -> `yelmo_regions_write`.
  The coefficients appended to `yelmo_ts.nc` today (`:182`) move to the
  forcing's file (section 8).
- `end`: `yelmo_end(ylmo, time)`.

## 3. `geo = fastisostasy` and its global part

D5: one component, two levels. The file `kryos_geo_fastisostasy.f90` holds
two extensions of `geo_t`: `fastiso_geo_t` (regional, one per domain) and
`fastiso_bsl_t` (global, at most one per run, allocated by the factory when
any domain selects `fastisostasy`). Constraints: regional: square cells
(`dx == dy`, `isos_utils.f90:289`), init serial (FFTW plan creation is not
thread-safe); global: `max_instances = 1`.

```
&geo_grl  backend = "fastisostasy"  grid = "GRL-16KM"  dt = 10.0  isos = "isos" /   ! grid_isos, dt_isos
&bsl      backend = "fastisostasy"  group = "barysealevel" /   ! globe; fixed name (barysealevel.f90:59-61)
```

`[barysealevel] method` is the `bsl` source of D5: `fastiso` = own, `const`,
`file`; `mixed` is refused in phase 1 (its accumulator `dbsl_total` needs the
per-domain bookkeeping of phase 3).

**Types.** `fastiso_geo_t`: `isos` (`isos_class`), `bsl_loc` (`bsl_class`, a
copy of the global state, section 3.2), `z_bed_ref_g` (`sp`, the domain
reference on the geo grid, kept for the check), `dx` (`sp`).
`fastiso_bsl_t`: `bsl` (`bsl_class`).

### 3.1 Imports, exports

| Instance | Import | Export |
|---|---|---|
| regional | init: `z_bed_ref`, `H_ice_ref` (domain, A2), `z_bed`, `H_ice` (ice exports); update: `H_ice` | `dz_bed` = `real(out%w + out%we, dp)` (sum in `sp`, as `kryos_coupling.f90:222`); `dz_ss` = `real(out%z_ss - now%bsl, dp)` (`sp`, as `:223`); `dz_bed_dt` = `out%dwdt`; `dbsl_vol` = `now%deltaV_bsl` [m3] (scalar, INSTANT: the step's volume change) |
| global | `dbsl_vol` per domain (gather) | `bsl` (scalar) |

`dwdt_corr` (Yelmo's `dzbdt_corr`): FastIsostasy copies it into a local array
that is never read (`fastisostasy.f90:812-816`), so wiring it or not is
bit-identical. Recommendation: declare nothing in phase 1 and drop the edge of
D4 until FastIsostasy uses it.

### 3.2 Step sequence (one coupling step, N domains)

1. Global, before any domain: `bsl_update(bsl, tm%t_rel)` (today the driver,
   `yelmox.f90:180`): `const` sets `bsl_now = bsl_init`, `file` interpolates
   the series, `fastiso` interpolates `A_ocean_now` from `bsl_now`.
2. Scatter: each regional `bsl_loc = bsl` (method, `bsl_now`, `A_ocean_now`).
3. Regional `update(tm)`: `H_ice` import narrowed;
   `isos_update(isos, H_ice_sp, tm%t1, bsl_loc)` (`kryos_coupling.f90:149`,
   model time). Inside, `fastiso` mutates
   `bsl_loc%bsl_now = bsl_loc%bsl_now - deltaV_bsl / A_ocean_now` and sets
   `isos%now%bsl` (`fastisostasy.f90:835-850`). Exports as 3.1.
4. Gather and global update (method `fastiso` only, and only in `update`,
   never in an init stage): `bsl%bsl_now = bsl%bsl_now - s` with
   `s = sum_d (dbsl_vol_d / bsl%A_ocean_now)` computed in `sp` in domain
   order; export `bsl`.

For N = 1 step 4 is FastIsostasy's own two `sp` operations in the same
order, so `bsl == isos%now%bsl` bit for bit and the core's `z_sl = bsl +
remap(dz_ss)` is today's `isos%now%bsl + (z_ss - isos%now%bsl)` (P1, P2). For
N > 1 the sum is order-independent up to the fixed domain order (today's
sequential mutation lets a domain see the previous one's change; removed on
purpose), and `dz_ss` stays relative to the domain's own `bsl_loc`.

### 3.3 `init`, `init_state`

- `init`: `isos_init(isos, par_file, isos_group, grid%G%nx, grid%G%ny, dx, dy,
  cnst=cnst, domain=..., grid_name=grid%name)` (`kryos.f90:348-350`) with
  `dx = real(grid%G%dx*grid%cs%xy_conv, sp)`. Today's `dx_isos =
  yelmo%grd%G%dx*(grid_i%G%dx/grid_y%G%dx)` (`kryos.f90:344`) is the same
  number for 16 and 32 km cells (16000 and 32000 are exact in `sp` by either
  route); the adapter logs it. Global: `bsl_init(bsl, par_file, tm%t_rel)` and
  `bsl_update` (`yelmox.f90:92-93`).
- `INIT_GEOMETRY` (cold start): `isos_init_ref(isos, real(z_bed_ref,sp),
  real(H_ice_ref,sp))`; `isos_init_state(isos, real(z_bed,sp), real(H_ice,sp),
  t, bsl_loc)`; publish exports (`kryos_startup.f90:236-241`).

  Why this equals today: today's reference is Yelmo's processed
  `bnd%z_bed_ref`, `bnd%H_ice_ref` on `grid_isos == grid_ice` (`:236-238`).
  The domain reference equals it where that processing is the identity
  (2.5): for `z_bed` (no gaps, no slope above `grad_lim_zb = 0.5`) the ice
  check proves it at init; for `H_ice` not in general (Q1). The current state
  (`z_bed`, `H_ice`) comes from the ice exports, as today (`:239-241`).

### 3.4 Restart, output, end

- `restart_write`: `isos_restart_write` (today `isos_restart.nc`); global
  `bsl_restart_write` (today `bsl_restart.nc`, `kryos_startup.f90:127`).
- `restart_read`: global first, `bsl_restart_read`; then scatter; regional:
  `isos%par%use_restart = .true.`, `isos%par%restart = file`,
  `isos_init_state(isos, z_bed, H_ice, t, bsl_loc)` with the restored ice
  exports (`kryos_startup.f90:503-508`; not a bare `isos_restart_read`, see
  the comment at `:458-463`). The redo guard in `isos_update` (`:801-804`)
  keeps `bsl` and `V_old` untouched.
- Check (today `check_isostasy_reference`, `kryos_coupling.f90:229-281`), now
  inside the adapter after every restart read: `max|out%z_bed - out%w -
  out%we - z_bed_ref_g| <= 1 m`, else `fesm_error` naming instance and file.
  It compares geo with the domain, as D4 states.
- `write_output`: `2D` -> `<instance>.nc` with `z_bed`, `z_ss`, `w`, `we`
  (backend internals, `prec_abs = 1e-3`, `kryos_output.f90:260-277`) plus the
  exports `dz_bed`, `dz_ss`; `1D` -> `<instance>_ts.nc` with `bsl`,
  `z_bed_mean/min/max`, `w_mean`, `we_mean` (`:383-399`). The global part
  writes `bsl`, `A_ocean` to `bsl_ts.nc`.
- `end`: `isos_end`; global `bsl_end`.

## 4. `climate = snapshots`

The adapter absorbs `yelmox_climate.f90` (snapshots branch only; the
`select case` goes), `snapshots.f90`, `climate_terms.f90`, `climate_mods.f90`
and `climate_out.f90`. These move into kryos as private modules of the
backend (`src/backends/climate/snapshots/`), with module names prefixed
`kryos_snap_*` so yelmox and kryos can be built together (design §1).
`climate_out_class` stays as the adapter's internal now/ref record; it is no
longer a framework type.

```
&clim_grl  backend = "snapshots"  grid = "GRL-16KM"  dt = 10.0   ! grid_clim, dt_clim
           snapshots = "snapshots"  reference = "pd"  terms = "lgm" "ocean"   ! today [forcing<sfx>]
           taper = ""  bias = "" /                                ! transform groups, "" = off
```

`[snapshot_*]`, `[term_*]`, `[coef_*]`, varslice field groups: unchanged, from
the par file or the run's library files (`var_defs`, design §9 `par_library`).

**Type.** `snapshots_climate_t` extends `climate_t`: `snap`
(`snapshots_class`), `terms` (`terms_class`), `mods`
(`climate_mods_class`), `out` (`climate_out_class`), `dair(:)` (derived
air-anomaly producers), `coefs` (`snap_coef_table_t`, below), `south`.
Capabilities: `has_ocn_prof = snap%par%with_ocn`, `has_tsl_ann = .true.`,
`has_ocn_shelf = .false.`, `has_smb = .false.`, `has_Qd = .false.`;
`with_atm` required (`yelmox_climate.f90:168-173`).

### 4.1 Forcing scalars without the registry

`climate_terms` and `yelmox_climate` call `coefs_has`, `coefs_is_derived`,
`coefs_derived`, `coefs_value`, `coefs_set` on the driver's `coefs_class`.
The adapter keeps these calls on a private table `snap_coef_table_t` (name,
derived flag and spec, value), filled at `init` from the `[coef_<name>]`
groups its terms name (as `kryos_forcing.f90:113-133`), and at each update
from its imports (narrowed, P1; `coef_class%f` is `sp` today,
`kryos_forcing.f90:53`). It declares one scalar import per tsgen coefficient
(`alpha`, and `dTo` as last set, A5) and one scalar export per derived
coefficient it produces (`dTo`, `from = "air_anomaly"`); a name without a
producer stops the wiring at init (D7).

### 4.2 `init`

`south = sum(grid%lat)/size(grid%lat) < 0` (`kryos.f90:365`);
`terms_init`; `snapshots_init(snap, par_file, domain, grid, tm%t_rel, south,
terms, snapshots_group, reference)` (`yelmox_climate.f90:164-165`);
`terms_check(terms, caps)`; the air-anomaly producers with their region mask
from `ref%reg` (`:211-244`); `climate_mods_init` (`:247`). Declare imports:
`z_srf` (hub geometry on the climate grid; the other geometry fields are
esm/rembo's and come in phase 2), forcing scalars. Declare exports (4.4).

### 4.3 `update(tm)` = `climate_update` with `method = "snapshots"`

In the order of `yelmox_climate.f90:294-343`, at `time = tm%t_rel`:
`terms_eval(terms, coefs)`; `snapshots_update(snap, z_srf, time, terms)`;
copy `snap%now/ref` into `out`; `out%now%z_srf = z_srf`, `out%lapse =
snapshots_lapse_monthly(snap)`, `out%now%beta_p = snap%now%beta_p`;
`terms_save_ocean`; `terms_apply_atm`; derived coefficients
`f = factor * sum(now%tsl_ann - ref%tsl_ann, mask) / n` in `sp` into the
table and the export; `terms_eval(derived=.true.)`; `terms_apply_ocn`;
`climate_mods_apply(mods, out, time)`. Then fill exports. `INIT_FIELDS` and
`INIT_RESTART` run the same body at `t` (today `init=.true.` and
`force=.true.`, which snapshots ignores; `kryos_startup.f90:101,162`): the
backend is stateless, so the re-derivation equals the restored exports bit
for bit whenever the restart time lies on the climate cadence, which kryos
requires (design §7).

### 4.4 Exports

| Export | Units | From | Conversion |
|---|---|---|---|
| `tas` (month) | K | `out%now%tas` | widen |
| `pr` (month) | kg m-2 yr-1 | `out%now%pr` [mm/d] | `real(pr, dp) * (sec_year/sec_day)` (run calendar); chion divides by the same factor (P3) |
| `z_srf` | m | `out%now%z_srf` (the surface the atmosphere refers to) | widen |
| `beta_p` | K-1 | `out%now%beta_p` | widen |
| `lapse` (month) | K m-1 | `out%lapse(m)` spread over the grid | widen |
| `tsl_ann`, `ta_ann` | K | `out%now` | widen (consumers in phase 2: smb_simple, glacial smb) |
| `T_ocn`, `S_ocn` (depth) | K, psu | `out%now%to_ann`, `so_ann` | widen; axis `depth` = `out%now%depth` (A8) |
| `dT_ocn` (depth) | K | `out%now%to_ann - out%ref%to_ann` in `sp` | widen (today `kryos_coupling.f90:612`, computed on the climate grid before the remap; exporting `to_ref` and differencing on the shelf grid would not be bit-identical in Antarctica) |
| `dTo` (scalar, derived) | K | the table value | widen; the core's `climate -> forcing` row calls `forcing_set` (section 8) |

All INSTANT at `t1`. Phase 1 exports profiles only: with `has_ocn_shelf`
false the shelf interpolates (section 6). The esm backend's shelf-base path
(`marshelf_interp_shelf` inside `esm_update`, `yelmox_climate.f90:842-849`,
on climate-grid geometry with the shelf's parameters) moves into the shelf
component in phase 2 (D6); that is a documented difference for esm configs,
not a phase 1 concern.

### 4.5 Output, restart, end

`write_output(2D)`: `<instance>.nc` (today `snap.nc`): `climate_write_init`
then `climate_write_2D` for snapshots (`:427-501`), same fields and
precision; no 1D. No restart file (stateless; `climate_restart_write`,
`:701-717`, writes nothing for snapshots). `end`: `snapshots_end`.

## 5. `surface = chion`

The adapter absorbs the yelmox chion host `surface_chion.f90` and
`orbit_forcing` (the insolation table).

```
&srf_grl  backend = "chion"  grid = "GRL-16KM"  dt = 10.0       ! grid_surface; today every step
          host = "surface_chion"  chion = "chion"  orbit = "orbit" /   ! groups unchanged
```

**Type.** `chion_surface_t` extends `surface_t`: `sc`
(`surface_chion_class`: `chion_class`, daily controls, insolation table,
annual fields `ann`). Constraints: none. The orbit stays host-side in phase 1;
orbital elements as forcing scalars are phase 3.

### 5.1 Imports and conversions

`tas`, `pr` (month, on the surface grid at the surface's `z_srf`, A6), `z_srf`,
`z_bed`, `H_ice`, `z_sl` (hub geometry, `kryos_coupling.f90:487-490`).
`t2m = real(tas, sp)`; `pr_mmd = real(pr / (sec_year/sec_day), sp)` (P3);
geometry narrowed. The re-elevation of A6 is the identity in both configs:
`grid_clim == grid_surface`, so yelmox returns before applying it
(`kryos_coupling.f90:529`); the core must keep that shortcut (apply only when
the source and target grids differ), because with `dt_clim > dt_surface` the
climate's `z_srf` lags the hub and `dz` is not zero.

### 5.2 `init`, `init_state`, `update`

- `init`: `surface_chion_init(sc, par_file, real(grid%lat, sp), group=host,
  chion_group=chion, orbit_group=orbit, cnst=cnst)` (`kryos.f90:410-412`).
- `INIT_FIELDS`: `surface_chion_spinup` (`time_equil` annual cycles) then
  `surface_chion_update` at `t_rel` (`kryos_coupling.f90:492-495`).
- `INIT_FINAL`: one more `surface_chion_update` on the geometry after the ice
  init (`kryos_startup.f90:192`). Both stages step the snowpack at `t_init`
  on purpose; `update` itself stays idempotent at `t0 == t1`.
- `update(tm)`: `surface_chion_update(sc, t2m, pr_mmd, z_srf, z_bed, H_ice,
  z_sl, tm%t_rel)`: **one 360-day annual cycle per call, whatever `t1 - t0`**,
  as yelmox (study §4.1 item 8). The SMB is that cycle's total over the
  360-day year, exported as the rate per year for `[t0,t1]`.

### 5.3 Exports

| Export | Units | From | Note |
|---|---|---|---|
| `smb` | kg m-2 yr-1 | `ann%smb` (sf + rf - runoff - subl, `surface_chion.f90:372`) | widen; MEAN [t0,t1]. No `sec_year`: the host sums daily kg m-2 over its 360-day cycle; kg m-2 yr-1 is numerically today's mm w.e./yr, so Yelmo receives the same number (2.1) |
| `T_srf` | K | `ann%tsrf` | widen, MEAN |
| `melt`, `runoff`, `refrz` | kg m-2 yr-1 | `ann%...` | widen, MEAN |
| `albedo` | 1 | `ann%alb_s` | masked: invalid where `MV` (PDD) |

### 5.4 Restart, output, end

`restart_write/read`: `surface_chion_restart_write/read` unchanged (chion
state + `ann_*`, today `chion_restart.nc`, `kryos_startup.f90:431-432,
518-519`); a file without `ann_*` stops (1.4). `write_output(2D)`:
`<instance>.nc` with the 13 annual fields of `surface_chion_write_step`
(`kryos_output.f90:326-361`), same units and precision. `end`: `chion_end`.

## 6. `shelf = marine_shelf` (incl. `bmb_method = "pico"`)

D6: the shelf owns obtaining `bmb`. In phase 1 the ocean -> shelf-base step
is marine_shelf's own `marshelf_update_shelf`, called by this adapter; when
`pico` as a separate backend or `data` arrive (phase 2), that vertical
interpolation is lifted into a free procedure of the shelf component used by
every shelf backend.

```
&shlf_grl  backend = "marine_shelf"  grid = "GRL-16KM"  dt = 10.0   ! grid_shelf (Antarctica ANT-16KM)
           marine_shelf = "marine_shelf" /       ! [pico] via pico_group (marine_shelf.f90:806)
```

**Type.** `marshelf_shelf_t` extends `shelf_t`: `mshlf` (`marshelf_class`),
`dx` (`sp`). Constraints: none. Requires `has_ocn_prof` or `has_ocn_shelf` of
its ocean source (init stops naming the producer otherwise).

### 6.1 `init`

`marshelf_init(mshlf, par_file, marine_shelf_group, nx, ny, domain,
grid%name, ref%reg, xc=real(grid%G%x, sp), yc=real(grid%G%y, sp), cnst=cnst)`
(`kryos.f90:425-427`); `dx = real(grid%G%dx*grid%cs%xy_conv, sp)` (today's
`dx_shelf`, the same number for 16/32 km cells). Imports: `H_ice`, `z_bed`,
`f_grnd`, `z_sl`, `z_srf` (hub, bilinear, `kryos_coupling.f90:592-596`), and
either `T_ocn`, `S_ocn`, `dT_ocn` (depth) or, with `has_ocn_shelf` (phase 2),
`T_shlf_ocn`, `S_shlf_ocn`, `dT_shlf_ocn`, `dS_shlf_ocn`.

### 6.2 `update(tm)` (also `INIT_FIELDS`): the two-call update

Profiles (phase 1, `kryos_coupling.f90:610-619`):

```fortran
call marshelf_update_shelf(mshlf, H_ice, z_bed, f_grnd, z_sl, dx, depth, to, so, dto_ann=dto)
call marshelf_update(mshlf, H_ice, z_bed, f_grnd, z_sl, dx=dx, z_srf=z_srf)
```

with narrowed imports and `depth` from the axis table (A8). Shelf-base path
(phase 2): `mshlf%now%T_shlf/S_shlf/dT_shlf/dS_shlf` from the imports, then
`marshelf_update` only (`:598-606`).

### 6.3 Exports

| Export | Units | From | Conversion |
|---|---|---|---|
| `bmb_shlf` | kg m-2 yr-1 | `now%bmb_shlf` [m ice/yr, melt < 0, `marine_shelf.f90:710`] | `real(bmb, dp) * real(mshlf%par%rho_ice, dp)`: exact product (24 + 10 bits); same `cnst%rho_ice` as Yelmo's inverse, so Greenland is exact by P3 and Antarctica (ANT-16KM -> ANT-32KM `con`) by 1.2's tie argument |
| `T_shlf`, `tf_shlf` | K | `now%T_shlf`, `now%tf_shlf` | widen |

INSTANT at `t1`: a diagnostic of the geometry at `t1`, used by the ice one
step later (lag 1 in the wiring, as today).

### 6.4 Restart, output, end

`restart_write/read`: `marshelf_restart_write/read` (today `marine_shelf.nc`,
`kryos_startup.f90:428, 514`; it holds `now`, incl. `tf_corr`).
`write_output(2D)`: `<instance>.nc` with `bmb_shlf` (m/yr, `prec_abs 1e-4`),
`T_shlf`, `tf_shlf` (`kryos_output.f90:279-294`). `end`: `marshelf_end`.

## 7. `geothermal = data`, `sediments = data`

Both are `component_t + datasrc_t`: `init` calls `datasrc_init(ds, par_file,
group, exports, [name], grid, "input/datasrc_defaults.nml", domain=...)`,
`datasrc_update(ds, t_init)` and `datasrc_fill` once; `update` is a no-op; no
restart file (re-served at init); default `write_output`. Today both are read
on `grid_ice` (`kryos.f90:352-356`): the instances use the ice grid, so
`{grid_name}` in the paths is unchanged.

**Geothermal**, export `Q_geo` [W m-2, + into ice]:

```
&ghf_grl  backend = "data"  grid = "GRL-16KM"  Q_geo = "ghf_50" /   ! Greenland: use_obs = False
&ghf_50   mode = "const"  value = 0.05 /                           ! [W m-2] = ghf_const 50 mW m-2
&ghf_ant                       ! Antarctica: use_obs = True, Shapiro2004, f_stdev = 0
    backend = "data"
    grid    = "ANT-32KM"
    Q_geo   = "ghf_shapiro2004"
    Q_geo_min = 1e-4           ! [W m-2] today's floor ghf_min = 0.1 mW m-2
/
&ghf_shapiro2004
    filename = "fesmdata/{domain}/{grid_name}/{grid_name}_GHF-Shapiro2004.nc"
    name     = "ghf"
    mode     = "static"
    fill     = "nearest"
/
```

Bit-identity: `0.05_dp * 1000 -> sp = 50.0`; the `float` file field in
mW m-2 is converted to W m-2 by `fesm_units` and back by Yelmo (P3). Checked
at init: datasrc fills the file's missing values, `read_ghf` every negative
value (`geothermal.f90:155-162`); they agree when the negatives are the
declared missing value (FesmData sets it).
The floor `Q_geo_min` (today `ghf_min`, `geothermal.f90:65,98`) is the one
line the geothermal data backend adds beyond datasrc (clamp after fill, in
`dp`; `real(1e-4_dp*1000, sp) == 0.1_sp`).

**Sediments**, export `H_sed` [m]: today `method = 0`, `H_sed = 0`
(`sediments.f90` `case(0)`). Phase 1 configs declare no sediments instance;
the ice import `H_sed` stays unwired and Yelmo keeps its default 0
(`yelmo_boundaries.f90:724`), bit-identical. `method = 1` is a field group
like `ghf_shapiro2004` (`filename = ..._Sediment-CRUST1.nc`, `name = "z_sed"`,
`remap = "bilin"`, `fill = "nearest"`), exported as `H_sed`.

Dropped: `f_stdev` x `obs_err_name` (an ensemble perturbation; both configs
use `f_stdev = 0`), and the integer `method` / `use_obs` switches (replaced by
the field group).

## 8. `forcing = tsgen`

Global, `max_instances = 1`. The adapter absorbs `kryos_forcing.f90`
(coefficient registry) minus its domain coupling.

```
&forcing  backend = "tsgen"  coefs = "alpha" "dTo" /   ! today [forcing] coefs; [coef_<name>] unchanged
```

**Type.** `tsgen_forcing_t` extends `forcing_t`: `coef(:)` (`coef_class` as
today: name, units, `restart_every_df`, `f`, derived flag and spec,
`tsgen_class`), `counter_restart`, `t_set(:)`, `setter(:)` (who set a derived
value this step).

- `init`: `coefs_init` body (`kryos_forcing.f90:76-160`) at `tm%t_rel`. One
  scalar export per coefficient (units from `[coef_<name>] units`), published
  at `t_init` (derived: 0, as today `:147`). Import: scalar `m_ice` [Gt]
  (optional; required only if a coefficient uses a feedback method).
- `update(tm)`: tsgen coefficients `tsgen_update(tsg, tm%t_rel, var=m_ice)`,
  `f = tsg%f_now` (`:162-176`), exported widened at `t1`; derived ones
  untouched. `m_ice` is today's `V_ice*rho_ice*1e-3` in `sp`
  (`yelmox.f90:185`), exported by the ice at the end of the previous step,
  so the same number (phase 3: summed over domains in domain order).
- `forcing_set(frc, name, value, t, producer)`: free procedure of the forcing
  component module, called only by the core's `climate -> forcing` wiring row
  (a backend never sees the forcing instance); stops on an unknown or tsgen
  name (`coefs_set`, `:178-195`) and on a second producer of one name at one
  `t` (the bipolar "both domains set dTo" bug, study §4.3).
- Run control, free procedures: `forcing_kill`, `forcing_restart_due`,
  `forcing_restart_fldr` (kill switch, `restart-<n>` bundles), called by
  `kryos_run` where the driver calls them today (`yelmox.f90:213-224`).
- `restart_write/read`: `<instance>_<name>_restart.nc` per tsgen coefficient
  (`tsgen_restart_write` writes whole files; today `coef_<name>_restart.nc`,
  `:328-388`). Derived values: the export state (A5).
- `write_output(1D)`: `<instance>_ts.nc` with `coef_<name>`, `_df_dt`,
  `_dv_dt` (`coefs_write_step`, `:390-427`). The component module cannot be
  called `kryos_forcing` (Q3).

## 9. Init and restart order the adapters rely on

Cold start (`kryos_startup.f90:131-195`), one domain:

| Stage | Instance | Action | Today |
|---|---|---|---|
| init | all | `init` (globe first: forcing, bsl; then ice, geo, climate, surface, shelf, data) | `domain_init`, `coefs_init`, `bsl_init` |
| GEOMETRY | geo, ice | geo ref + state from domain ref and ice exports; exchange geo -> ice; ice `z_bed`/`z_sl` + diagnostics; core: hub refresh | `domain_init_isostasy`, `couple_yelmo_to_htopo` (`:146-156`) |
| FIELDS | climate, surface, shelf | snapshots at `t_rel`; chion spin-up + 1 cycle; shelf two calls; exchange to ice | `:162-169` |
| STATE | ice | `bnd` from imports, `yelmo_init_state`, `equil`/thrm; core: hub refresh | `:178-186` |
| FINAL | surface | 1 cycle on the final geometry; exchange surface -> ice | `:192-193` |

Restart (`domain_startup`, `:95-103`): restore the core's export states (A5);
`restart_read` of forcing, bsl, ice, geo (needs `H_ice`), shelf, surface;
exchange geo -> ice; `INIT_RESTART` of ice (diagnostics, regions) and
climate (re-derivation); core: hub refresh.

## 10. Bit-identity per crossing field

Greenland: every grid is GRL-16KM, every remap a copy, the hub a mirror of
Yelmo (`kryos_coupling.f90:557-562`). Antarctica: hub and shelf ANT-16KM;
ice, geo, climate, surface ANT-32KM.

| Field | Path | Greenland | Antarctica |
|---|---|---|---|
| `z_bed`, `z_sl` | geo -> ice (identity) | exact (P1, P2), given 2.5 | exact (`grid_isos == grid_ice`), given 2.5 |
| `H_ice` | ice -> geo | exact (P1) | exact (identity) |
| `bsl` | global reduction | exact (3.2) | exact |
| `smb`, `T_srf` | chion -> ice | exact (P1, P5) | exact (identity) |
| `bmb_shlf` | shelf -> ice | exact (P3) | `con` 16 -> 32 km after `* rho_ice`: exact except at `sp` ties (1.2) |
| `T_shlf` | shelf -> ice | exact | exact (P4) |
| `Q_geo` | data -> ice | exact (P3) | exact (P3), given the fill condition (7) |
| `tas`, `pr` | climate -> chion | exact (P1, P3) | exact (identity; `pr` P3) |
| `T_ocn`, `S_ocn`, `dT_ocn` | climate -> shelf | exact | exact by P4 (`bilin` 32 -> 16 km of widened `sp` values) |
| geometry | hub -> climate, surface, shelf | exact (mirror, copies) | exact only if the hub's refresh narrows to `sp` where `htopo_update` does today; a `dp` hub (design §12 row 8) differs at round-off (Q2) |
| forcing scalars | forcing -> climate | exact (P1) | exact |

## 11. yelmox features these adapters drop or defer

| Feature | Where today | Phase 1 | Returns |
|---|---|---|---|
| `relax` (ytopo ramp) | `step_relax`, `domain_relax_init` | dropped | 2, plugin on `yelmo_ice_t` (Antarctica gate, Q2) |
| `opt` (cb_ref L21, tf_corr) | `step_optimize`, `domain_opt_init(_cb_ref)` | dropped; `tf_corr` stays in the shelf restart file | 2, plugins on `yelmo_ice_t` + `marshelf_shelf_t` (`tf_corr` plugin, D6) |
| NEGIS cb_ref and its cold-start init (`use_negis`) | `negis_update_cb_ref`, `domain_negis_init` | dropped | 2, plugin |
| glacial SMB scaling | `couple_surface_to_yelmo:326-330` | dropped | 2, plugin (needs climate `ta_ann`, exported already) |
| `lim_pd_ice` (-4 m/yr) | `:334-337` | dropped | 2, plugin |
| `init_marine_H`, `init_kill_shelves`, `recon`, `recon_ref` | `domain_init_ice`, `domain_init_recon` | dropped (`none`, `equil`, `init_time_thrm` kept) | 2, cold-start plugins |
| `ice_bnd` (prescribed `H_ice_ref` in fixed cells) | `step_ice_bnd` | dropped (inactive in both configs) | 2, `data` instance + optional `H_ice_ref` import of yelmo |
| `with_ice_sheet = False` and the other `with_*` | `[comps]` | dropped | 2, `data` backends |
| CMIP output | `cmip_output.f90` | dropped | 2, post-processing over ice exports (design §9) |
| coefficients in `yelmo_ts.nc` | `kryos_output.f90:182` | moved to `forcing_ts.nc` | 1 |
| region 1D files `yelmo_ts_<region>.nc` | Yelmo-internal | kept (Yelmo writes them) | 1 |
| `isos.nc`, `mshlf.nc`, `chion.nc`, `snap.nc` names | framework writers | instance names | 1 |
| `bsl method = "mixed"` | barysealevel | refused | 3 (per-domain accumulator) |
| `dzbdt_corr -> isos` | `step_isostasy:147` | not wired (no effect, 3.1) | when FastIsostasy reads it |
| esm, rembo; smbpal, smb_simple, `surface_method = "climate"` | yelmox_climate, kryos.f90 | not in phase 1 | 2 |
| geothermal `f_stdev` | `geothermal.f90:86-94` | dropped | 2 if needed (datasrc perturbation) |
| `series_<name>_restart.nc` fallback | `kryos_forcing.f90:361-370` | dropped | never |

## 12. Open questions for Alex

1. **The isostasy reference (D4 versus today).** Today FastIsostasy's
   reference is Yelmo's *processed* present-day geometry
   (`kryos_startup.f90:236-238`): `H_ice < 1 m` set to 0, englacial lakes
   removed, ice outside `mask_ice_dynamic` removed, bed gradient-limited
   (`yelmo_data.f90:254-289, 421-424`). D4 initialises geo from the domain
   (hub) reference, which keeps that ice, so the reference load differs where
   Yelmo's processing acts (e.g. land ice outside `region:Greenland,
   Canadian_Arctic_Archipelago`) and phase 1 is not bit-identical. Options: (a) a yelmox preparatory commit that initialises
   isostasy from the hub reference, which becomes the baseline kryos must
   match (recommended: keeps D4 clean); (b) the ice exports its processed
   reference and geo uses it in phase 1 (bends D4); (c) move Yelmo's
   reference processing into the hub load (larger, changes Yelmo's inputs).
   The bed part is probably the identity on 16/32 km grids; the ice adapter
   checks it at init (2.5).
2. **The Antarctica gate.** `yelmox_Antarctica.nml` runs with `relax = True`
   and `opt = True`, which phase 1 drops (plugins, phase 2). Gate against a
   yelmox run with `-p sim.relax=False sim.opt=False`? And with its
   different hub/shelf grids, strict bit-identity also needs the hub refresh
   to narrow to `sp` as `htopo_update` does today; with the `dp` hub of
   design §12 row 8 Antarctica gives round-off differences (plus the rare
   `bmb_shlf` ties, 1.2). Accept "Greenland strict, Antarctica documented
   differences within a tolerance", or require an `sp`-faithful hub in phase
   1?
3. **Module names.** yelmox has `kryos_forcing` and `kryos_output`; design §1
   forbids reusing yelmox module names, but design §13 lists core
   `kryos_output`, and the forcing component module would be `kryos_forcing`.
   Rename (e.g. `kryos_forcings`, `kryos_outputs`), or drop the rule if the
   two are never linked into one executable?
