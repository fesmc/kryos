# fesm-utils: a study for the kryos design

Source studied: `/Users/alrobi001/models/fesm-utils` (branch `dev`, HEAD `d994f08`, clean tree), plus its
consumers `yelmox/libs/kryos*.f90`, `yelmox/yelmo/src`, `yelmox/FastIsostasy/src`, `climber-x/src`, and
`vilma2/src`. Paths below are relative to `fesm-utils/src/` unless prefixed. Line numbers are as of the
date of this study (2026-10-24).

## TL;DR

- fesm-utils is a solid library of **mechanisms**: grids/projections/remapping (`coords`), a grid
  registry + map cache (`coupler`), netCDF (`ncio`), namelists (`nml`), slicing of forcing files along
  time (`varslice`), physical constants (`phys_constants`), FesmData conventions (`fesmdata`,
  `topodata`, `regions`), and small time utilities. It contains **no framework concepts**: no
  component abstraction, no field registry, no scheduler, no restart manager, no logging.
- The real workhorses downstream are `ncio`, `nml`, `coords` (mapping), `varslice`, `phys_constants`
  and `regions`. `coupler`, `htopo` and `tsgen` have exactly one consumer (yelmox `kryos*`). Several
  packages (`multitopo`, `zremap`, `series` outside tsgen, `timeout`, `staggering`, `topodata` outside
  htopo) have no external consumer.
- **The prototype of kryos already exists** in `yelmox/libs/kryos*.f90` (a module literally named
  `kryos`): a domain = bag of concrete component objects + a coupler + an htopo hub, a
  `cadence_due` predicate, restart bundles (one folder per time, one file per component), a coefficient
  registry. It is hard-wired to Yelmo/FastIsostasy/etc. It is the best guide to what kryos must
  generalize.
- **Global<->regional is not a first-class concept.** The mapping core can map between any lon-lat /
  Gaussian / projected grids, conservatively, but (a) the conservative apply normalizes by covered
  area, so integrals are not conserved across partial coverage and no coverage fraction is returned;
  (b) Gaussian cell areas come from axis midpoints, not Gauss weights (VILMA2 works around this with a
  global rescale); (c) the map cache is keyed by grid *names* only; (d) there is no merge/split of N
  regional fields into one global field.
- **Module-name collisions are a live hazard**: climber-x defines its own modules `precision`
  (wp = dp), `constants`, `timer`, `coupler`; FastIsostasy defines its own `precision`. All Linux
  builds link with `-Wl,-zmuldefs` to silence duplicate symbols. A from-scratch framework should prefix every
  shared module.
- Best candidates to **add to fesm-utils**: a shared error/log module (10 private copies exist), a
  field/exchange-state type, coupler extensions for conservative global<->regional exchange
  (coverage fraction, integral normalization, merge/split, grid identity hashing, Gauss-weight cell
  edges), an integer-tick cadence clock, and a restart-bundle/manifest helper. The **abstract
  component type** and the scheduling *policy* are better kept in kryos until a second consumer needs
  them.

---

## 1. Inventory

### 1.1 Who uses what (real consumers)

Counted as files with a `use <module>` line. "YX" = `yelmox/libs` (incl. `kryos*`), "Y" = yelmo,
"FI" = FastIsostasy, "CX" = climber-x/src, "V2" = vilma2/src.

| Module | YX | Y | FI | CX | V2 | Notes |
|---|---|---|---|---|---|---|
| `ncio` | 23 | 10 | 4 | 72 | 6 | universal |
| `nml` | 23 | 9 | 4 | 23 | 2 | universal |
| `precision` | 16 | – | (own) | (own, wp=dp) | (own) | see §8.1 |
| `coords` (umbrella) | kryos, kryos_output | 5 | – | ~29 | 3 | CX uses only grid/map |
| `mapping`, `coordinates`, `interp2D` | 3 / 7 / 3 | 2 / – / 2 | – | – | – | legacy direct uses |
| `varslice` | 4 | – | – | – | – | yelmox forcing |
| `phys_constants` | 5 | 1 | 1 | 2 | – | |
| `regions`, `fesmdata` | 3 / 2 | 2 / 1 | – | – | – | |
| `timestepping` | 5 | – | – | – | – | |
| `coupler` | kryos only | – | – | **own module of same name** | – | |
| `htopo`, `tsgen` | kryos only | – | – | – | – | |
| `variable_io`, `subgrid`, `gaussian_quadrature`, `distances` | – | 3 / 2 / 7 / 1 | – | – | – | yelmo internals |
| `timer` | – | – | – | **own module of same name** | – | |
| `constants` | 1 | – | – | **own module of same name** | – | |
| `multitopo`, `zremap`, `series`, `timeout`, `topodata`, `staggering`, `derivatives`, `root_finder`, `ncio_transpose` | – | – | – | – | – | no external consumer (series via tsgen, topodata via htopo) |

climber-x's effective fesm-utils surface is therefore just `ncio`, `nml`, `coords` and
`phys_constants` (`climber-x/src/main/constants.f90:26`, `timer.f90:26`, `coupler.f90:27`, `climber-x/src/utils/precision.f90:37`
shadow the fesm-utils modules of the same names).

### 1.2 Package table

| Package (file) | Provides | Public API shape | Maturity | Known issues |
|---|---|---|---|---|
| `precision` (13 lines) | kinds `dp`, `sp`, `wp=sp` (`precision.f90:7-11`) | parameters | stable | `wp=sp` is a library-wide compile-time switch; several downstream codes define their own `precision` with different `wp` |
| `constants` | `mv=-9999._wp`, `mv_dp`, `pi`, deg/rad, `TOL`, `ERR_DIST/IND` (`constants.f90:16-30`) | parameters | stable | `mv` is `wp`; coords internals use `mv_dp` |
| `phys_constants` (860) | `phys_const_class` record of 12 primitives + 10 derived, all `real(dp)`, sentinel `PHYS_CONST_UNSET` (`phys_constants.f90:120-153`); `phys_const_load(c, file, group, defaults_file, defaults_group)` (204), `_set` (283), `_derive` (325), `_require` (387), `phys_const_get(c, name, out)` generic sp/dp (158, 409, 437), `_write/_log/_to_nc/_compare`; calendar constants `sec_year_*` (62-78) | record + by-name getter | mature, well designed (no defaults in code, loud failure) | not yet adopted by VILMA2 (`vilma_constants`) or climber-x core (own `constants`) |
| `nml` (1528) | `nml_read(file, group, name, value, comment, init, defaults_file, defaults_group)` for scalars and vectors of str/dp/sp/int/logical (18-24, 681-1000); `nml_validate` (634); `nml_has_param/_has_group` (94,109); `nml_replace` (63); flat `group.param = value` syntax | free functions, rescan file per key | mature, heavily used | see §5 for limits |
| `ncio` (4572, generated from `ncio.fypp` 2751) | `nc_create/open/close`, `nc_write`/`nc_read` for all kinds, ranks 0-6, `start/count`, `missing_value`, `prec_abs/prec_rel` bit rounding, `nc_write_dim`, `nc_write_map` (CF grid mapping), `nc_time_index`, `nc_set_lossless` (`ncio.f90:181-191`) | stateless by filename (optional `ncid`) | mature, universal | serial netCDF only; opens/closes per call unless `ncid` passed |
| `ncio_transpose` | `nc_write_t/nc_read_t` (transpose in memory) | generic | small | |
| `variable_io` (717) | `var_io_type` (name, dims, units, long_name, prec) from a markdown variable table; `load_var_io_table`, `find_var_io_in_table`, `var_io_write` (`variable_io.f90:8-37`) | table-driven metadata + writer | used by yelmo only | max 1000 rows (`:107`); 32-char names/units |
| `varslice` (2366) | `varslice_class`: one variable from one file (or glob), with time slicing (`exact/range/interp/extrap/range_mean/...`), sub-annual axis, climatology, unit scale/offset, missing fill, online remap via `map_class` (`varslice.f90:25-95`, `varslice_update` 430, `varslice_init_nml` 1249) | object per variable | mature, yelmox forcing workhorse | time axis generated from namelist `time=[active,x0,x1,dx,nsub]`, not read from the file; data stored `real(wp)` 4D; units only as scale/offset |
| `fesmdata` (692) | FesmData v2 conventions: flag tables (CF `flag_values/meanings`), file grid via `grid_name` attribute + `grid_<GRID>.txt`, `fesmdata_read_field` (remapped if needed), `{domain}/{grid_name}` path placeholders, unit helpers `fesmdata_flux_conv`, `fesmdata_temp_ref` (`fesmdata.f90:41-47, 264-410`) | free functions | stable | unit parsing is ad hoc (water flux + temperature only) |
| `topodata` (351) | read `<GRID>_TOPO-*.nc` (z_bed, H_ice, z_srf, fractions, mask), remap con/nn (`topodata.f90:1-22`) | object | stable | only via htopo |
| `regions` (1046) | categorical integer layers + name tables, hierarchical codes, selection expressions -> logical masks (`regions_select`, `regions_mask`, `regions_basin_ids`) (`regions.f90:88-95`) | object + expressions | mature | |
| `coords/*` (~13k) | geometry, projections, Gaussian latitudes, k-d tree, polygons, `weight_map_t`, `map_class`, conservative clipping, SCRIP-superset I/O, cdo interop, `nc_read_interp` | see §2 | mature, extensively tested (34 test programs in `test/coords`) | see §2 and `docs/dev/coords-followups.md` |
| `coupler` (357) | grid registry + map cache + `remap` (§2.4) | object | young, one consumer | fixed capacities, no masks/fractions pass-through, no identity shortcut |
| `htopo` (307) | hi-res geometry hub: reference topo + regions + current `z_bed/H_ice/z_srf/f_grnd/z_sl`, `htopo_update(dz_bed, dH_ice, z_sl)` (`htopo.f90:48-61, 136-157`) | object | young, kryos-only | fields `real(wp)=sp`; `H_ice` clipped at 0 (not mass-conserving) |
| `multitopo` (528) | two-grid consistent topography (coarse<->fine, fill-level ice reconstruction) (`multitopo.f90:1-26`) | object | v1, unused | "v1: two grids, sharing a projection" |
| `zremap` (686) | elevation remapping of SMB/temperature via local f(z) profiles (`zremap.f90:1-48`) | object init/update/map/end | recent, unused downstream | |
| `timestepping` (350) | `tstep_class` (§3) | object | stable | units hard-wired to "year" |
| `timeout` (492) | output-time lists (§3) | object | stable | `nmax=100000`; uses bare `stop` |
| `timer` (514) | CPU/wall profiling, `ncomp_max=100` | object | stable | `real(8)` style |
| `series` (529), `tsgen` (1357) | tabulated series; transient scalar forcing (ramps, PI controllers, noise, restart) | objects | tsgen mature | |
| `subgrid`, `derivatives`, `distances`, `staggering` | array kernels | free functions | `staggering` is a **non-callable stub** (no public list, undefined helper; `docs/packages/gridops.qmd:17-28`) | |
| `gaussian_quadrature` (1427) | FE quadrature on hexahedra (yelmo) | `gq2D_class` | vendored CISM style | |
| `root_finder` | Newton, secant | functions | small | |

Infrastructure: `docs/` Quarto site (published to gh-pages from `dev`), 10 test folders driven by
`make test-*` targets (`config/Makefile:294-490`) but **no CI that builds or runs Fortran tests**
(`.github/workflows/publish-docs.yml` only renders docs). Builds with a configme-generated Makefile
into `include-serial/` and `include-omp/` (`README.md`), and vendors FFTW, LIS and SHTns.

### 1.3 Ongoing modernization (docs/dev)

- `sustainability-plan.md`: dead-code removal (done), coupler precision note (done), error handling
  with rich context piloted in `varslice` and three other modules, each with a **private copy** of
  `*_error(proc,msg,detail)` + `to_str` ("if a fifth module needs it, reconsider a shared module").
  There are now 10 copies (`coupler_err`, `fesmdata_error`, `htopo_error`, `gaussian_quadrature_error`,
  `regions_error`, `topodata_error`, `varslice_error`, `coordinates_error`, `zremap_error`,
  `ncio_interp_error`) and ~70 bare `stop`s remain (e.g. `series.f90` 12, `variable_io.f90` 7,
  `coords/mapping_scrip.f90` 6, `coords/planet.f90` 6). Parked: fypp templating, `make check` + CI,
  `fpm.toml`/CMake, slimming vendored trees, F2008 submodules.
- `coords-followups.md`: precision/constants consolidation is the headline proposal (partly done: the
  `constants` module now exists); ll->stereo accuracy caveat (max|coords-cdo| ~0.2 vs ~3e-3 for
  ll->ll, conservation exact).
- `coords-adoption.md`: plan to migrate downstream from `map_scrip_*` to `map_class` + `map_field`,
  with `gen="cdo"` kept as generator until the analytic clip is trusted.
- `coords-performance.md`: analytic conservative weights are 90-100x faster than cdo for ll->ll,
  1.5-3x for ll->stereo; structured nn/bilinear 2-10x.

---

## 2. The grid / mapping / coupling stack

### 2.1 Grid description

`coords/coordinates.f90:30-69`:

```fortran
type coord_system                      ! shared metadata, embedded as %cs
    character(256) :: mtype, units     ! latlon|gaussian|cartesian|stereographic|polar_stereographic|
                                       ! lambert_azimuthal_equal_area|transverse_mercator|rotated_pole
    type(planet_class) :: planet ; type(projection_class) :: proj
    real(dp) :: xy_conv ; logical :: is_cartesian, is_projection, is_lon180
end type
type grid_axis_class ; integer :: nx, ny ; real(dp) :: x0, dx, y0, dy ; real(dp), allocatable :: x(:), y(:) ; end type
type grid_class
    character(256) :: name ; integer :: npts ; type(coord_system) :: cs ; type(grid_axis_class) :: G
    real(dp), allocatable, dimension(:,:) :: x, y, lon, lat, area ; integer, allocatable :: border(:,:)
end type
type points_class   ! scattered points: x, y, lon, lat, area, dx, dy (1D)
```

- Grids are always **logically rectangular with 1D axes** (regular or non-uniform, e.g. Gaussian
  latitudes). There are no stored cell bounds: cell edges are axis midpoints, extrapolated at the
  ends (`grid_area`, `coordinates.f90:135-206`). Curvilinear data are handled as `points_class` (e.g.
  `varslice` with 2D lon/lat, `varslice.f90:83`). No unstructured meshes.
- Projections: oblimap2 oblique stereographic/polar stereographic/LAEA, transverse Mercator, rotated
  pole (`coordinates.f90:845-871`). Planet ellipsoid via `planet_class` (WGS84 etc.).
- Global wrap/poles are inferred heuristically from distances for latlon/gaussian grids (`border`
  flags, `coordinates.f90:586-621`); there is no explicit `is_global` attribute.
- Gaussian grids: `gaussian_latitudes_calc(nlat, lat, wts)` (`coords/gaussian_latitudes.f90:20`)
  returns latitudes **and Gauss weights**, but `grid_class` ignores the weights; areas come from the
  midpoint rule.
- Persistence: `grid_<name>.txt` cdo-compatible griddes files with a `#` header carrying the full
  `coord_system` (`grid_cdo_write_desc_short`, `grid_cdo_read_desc`, `coords/grid_cdo.f90:18,169`);
  the reader falls back on CF keys (`grid_mapping_name`, `standard_parallel`, `proj_params`, ...).
  `fesmdata_grid_read` finds `grid_<GRID>.txt` next to a FesmData file or in `maps/`.

### 2.2 Map building and caching

- `weight_map_t` (`coords/weight_map.f90:27-43`) is the single CSR sparse store, `kind =
  MAP_DISTANCE` (per-link distance/position/quadrant; weights computed at apply time so one cached
  map serves nn/shepard/quadrant and adapts to a time-varying source mask) or `MAP_WEIGHT` (baked
  weights: conservative, bilinear, imported SCRIP).
- `map_class` (`coords/mapping.f90:36-47`) = `weight_map_t` + target geometry + names + method.
- `map_init` generic: grid->grid, grid->points, points->grid, points->points, and by grid names
  (`mapping.f90:49-53`). For grid->grid (`mapping.f90:80-193`):
  - cache file `fldr/scrip-<method>_<src>_<dst>.nc` (`gen_map_filename`,
    `coords/mapping_scrip.f90:970-984`) is loaded if present (format auto-detected: coords superset or
    plain cdo SCRIP), else generated and saved (atomic tmp + rename, `mapping.f90:295-331`);
  - `gen="coords"`: `method="con"` -> analytic conservative (`coords/conservative.f90:102-131`:
    separable analytic overlap for same-plane cartesian and for ll->ll/Gaussian, Sutherland-Hodgman
    planar clip for ll->projected, great-circle clip with `NSUB=8` subsampling for
    projected->latlon); nn/bilinear on a regular source -> direct index location
    (`map_init_structured`, `mapping.f90:728`); otherwise k-d tree neighbours (`map_init_internal`);
  - `gen="cdo"`: shells out `cdo gen<method>` (`map_init_cdo`, `mapping.f90:333`);
  - `precompute_weights` bakes shepard/bilinear distance maps into weights at init;
    `weight_cutoff=1e-4` prunes sliver links (`weight_map_prune`).
- **Cache key = grid names + method only.** Two different grids with the same name collide silently.
  VILMA2 already tags grid names with dimensions to work around this and documents the residual risk
  (`vilma2/src/vilma_remap.f90:14-19`).

### 2.3 Apply semantics (important for conservation)

`map_field(map, name, var1, var2, stat, method, missing_value, radius, mask2, reset, mask_pack,
fill_method, filt_method, filt_par, verbose)` is generic over grid/points x dp/sp/int (12 specifics,
`mapping.f90:55-62`). Internally everything is `real(dp)`. Source masking is done by setting source
values to `missing_value`; `mask_pack` restricts targets; `fill_method`/`filt_method` post-process.

`weight_map_apply` with `stat="mean"` (`coords/weight_map.f90:154-275`) computes, per target,
`sum(w*v)/sum(w)` **over valid sources only**. For conservative maps `w` are overlap areas, so this
is the area-weighted *mean over the covered part* of the target cell (SCRIP "fracarea"
normalization). Consequences:

- Mapping a constant gives the constant (this is the "cons" check in `coords-performance.md`).
- The area **integral is conserved only where a target cell is fully covered** by valid source
  cells. At domain edges or under masks, the partially covered target gets the covered mean, i.e.
  the integral is inflated by `1/frac`.
- The covered fraction is **not returned**; `mask2` is a logical "got a value" flag
  (`weight_map.f90:173, 221, 268`). There is no "destarea" (divide by full target area) option and no
  `stat="sum"`.
- For ll->projected the overlap areas are computed in the target projection plane
  (`conservative.f90:7-16`), while `grid%area` of a projected grid is the true ellipsoidal area of
  the inverse-projected cell (`coordinates.f90:179-190`). The two metrics agree only to the extent
  the scale factor is constant over a cell; any integral check should use one metric consistently.

For ice-load exchange with a global solid-earth model this is the central gap (§2.8).

### 2.4 The coupler registry

`coupler.f90` (the whole module is 357 lines):

```fortran
integer, parameter :: GRID_MAX = 16, MAP_MAX = 64          ! coupler.f90:35-36
type coupler_class                                          ! coupler.f90:43-51
    type(grid_entry) :: grids(GRID_MAX) ; integer :: ngrids
    type(map_class)  :: maps(MAP_MAX)   ; integer :: nmaps
    character(512)   :: map_fldr = "maps" ; character(32) :: gen = "coords"
end type
call coupler_init(cpl, map_fldr, gen)                        ! :73
call coupler_add_grid(cpl, name, grid)                       ! :89 (stores a deep copy)
call coupler_prime(cpl, src, dst, method, gen)               ! :109 build/load eagerly
call remap(cpl, var_src, src, var_dst, dst, method, stat)    ! :53-65 dp/sp/int x 2D/3D
```

- Grids are nodes, maps are (src, dst, method) edges; lookup is a linear scan with `trim` string
  compares (`get_map`, `:140-179`). Unregistered names are resolved from `grid_<name>.txt`
  (`resolve_grid`, `:181-195`).
- Fixed-size arrays are deliberate so that adding a map never reallocates/deep-copies existing maps
  (`:11-14`). Overflow stops with "raise it in coupler.f90". (`allocatable` arrays of a wrapper
  holding `type(map_class), pointer` or `move_alloc` growth would remove the cap without copies.)
- `remap` forwards only `stat`. It does **not** forward `missing_value`, `mask_pack`, `fill_method`,
  `filt_*`, nor return `mask2`. It has **no identity shortcut**: yelmox wraps it to copy when
  `src == dst` (`yelmox/libs/kryos.f90:714-751`). 3D fields loop the 2D map over the last index.
- The `dp`/`sp` overloads are intentionally concrete kinds, not `wp` (`:57-62`, and
  `sustainability-plan.md` Step 2).
- Tests: cartesian coarse<->fine and polar-stereographic disk maps only (`test/coords/test_coupler.f90`).
  No global/regional case.

In yelmox each `kryos_domain` owns its own `coupler_class` (`kryos.f90:176`), so in the bipolar
driver two domains build two independent registries (`yelmox_bipolar/yelmox_bipolar.f90:71`).

### 2.5 Precision handling

- Geometry, weights and accumulation: always `dp` (coords design decision, `coords-design.md`).
- `map_field`/`remap`/`nc_read_interp`: generic over `dp`, `sp`, `int` (and logical for
  `nc_read_interp` 2D).
- Library state that is "model-like" is `real(wp)` with `wp=sp`: `htopo` fields (`htopo.f90:55-59`),
  `varslice%var` (`varslice.f90:77`), `zremap`, `tsgen` values, `constants%mv`.
- `phys_constants` is all `dp` and narrows on `phys_const_get` into either kind (`:158-161`).
- Downstream working precisions differ: yelmo `wp=sp` (`yelmo/src/yelmo_defs.f90:42`),
  FastIsostasy `wp=sp` but its ice module `dp`, climber-x `wp=dp`, VILMA2 `wp=real64`.

### 2.6 Hub and helpers: htopo, regions, zremap, multitopo

- `htopo` is the domain's hi-res geometry source of truth: reference (`topodata`) + regions + current
  geometry refreshed as reference + model anomalies (`htopo_update`, `htopo.f90:136-157`). The
  anomaly form is the right pattern for receiving uplift from a coarse/global solid-earth model, but
  `H_ice` is clipped (`max(.,0)`), so the hub is not mass-conserving, and it is single precision.
- `regions` provides named categorical layers and expression masks; it is the natural source for
  "ownership" masks in a global merge (which domain owns which cell).
- `zremap` (elevation downscaling of SMB/T) and `multitopo` (coarse<->fine consistent topography with
  ice-volume-conserving fill-level reconstruction) are regional sub-grid tools; neither is used yet.

### 2.7 Global <-> regional today

What works:

- Any direction between lon-lat, Gaussian and projected grids, conservative or bilinear/nn, in one
  call, cached. Regional source onto global target leaves out-of-domain targets unmapped (fixed after
  a real incident where an NH grid clamped onto Antarctica and deleted its ice sheet in climber-x,
  `test/coords/test_regional_drop.f90:1-12`).
- VILMA2 already couples a **host lon-lat grid** to its **SHTns Gauss-Legendre grid** with coords:
  conservative ll->Gauss for mass-bearing fields, bilinear Gauss->ll for smooth fields, plus a global
  rescale so the SHTns quadrature integral equals the source integral, because coords Gauss cell
  areas are not the quadrature weights (`vilma2/src/vilma_remap.f90:1-29`, error O(1/nlat^2)).
- Global scalars are handled outside: barystatic sea level is a driver-owned `bsl_class` passed to
  each domain's isostasy step (`yelmox_bipolar.f90:64,126,207,220,228`).

What is missing for VILMA2 (global, Gauss grid) exchanging loads/uplift with N regional domains:

| Gap | Why it matters | Where it should go |
|---|---|---|
| Coverage fraction + integral ("destarea") normalization in apply | load of a global cell partly covered by a regional domain must be `sum(overlap*H)/A_cell`, plus a fraction to blend with the background | `weight_map_apply` + `map_field(frac=, norm="dst")` |
| Merge N regional fields into a global field with ownership/priority masks and a background field | overlapping domains (e.g. Greenland inside an NH domain) must not double count; outside all domains use data or a global ice model | new `coupler_merge` |
| Split/scatter global -> each regional domain, anomaly form | uplift/geoid as anomalies relative to a reference, like `htopo_update` | `coupler` + htopo pattern |
| Gaussian cell edges from Gauss weights (`sin(phi_{j+1/2}) = sin(phi_{j-1/2}) + w_j`) or explicit bounds in `grid_class` | removes VILMA2's global rescale; makes coords areas equal the spectral quadrature | `coordinates` (optional `y_bnds` / `area` override) |
| Grid identity (hash of cs + axes) in the cache key | grids with equal names but different axes collide | `gen_map_filename` / map file metadata |
| Time accumulation of fluxes/loads over the coupling interval | ice updates every 1-10 yr, VILMA2 every 10-100 yr; loads should be interval means | exchange-state accumulator (§6.2) |
| Global integral diagnostics (ice volume, SLE) per grid, consistent area metric | verify conservation per exchange | `coupler_integral(cpl, grid, field)` |
| `dp` throughout the exchange path | sp sums of ~1e6 cells drift | exchange buffers in dp |
| Pass-through of masks/missing/fill in `remap`, and identity shortcut | today only `stat` is forwarded | `coupler` |

---

## 3. Time

| Piece | What it does | Limits |
|---|---|---|
| `timestepping` `tstep_class` (`timestepping.f90:17-39`) | one global clock in `real(dp)` years: `time`, `time_cal`, `time_rel` (= cal - `time_ref`, default 1950), `time_elapsed`, step count `n`; methods `const` (calendar frozen, e.g. equilibrium spin-up), `cal`, `rel`; Kahan-compensated accumulation (`:328`); `tstep_init` from args or namelist (`tstep_method, tstep_const, time_init, time_end, dtt`, `:55-97`); `is_finished` | units fixed to "year" (`:88,91`); no sub-annual calendar; dt passed in each `tstep_update(ts, dt)` call, no per-component clocks |
| `timeout` `timeout_class` (`timeout.f90:15-24`) | list of output times from `method = none | const (dt) | times | file` (file lines `t0:dt:t1`); `timeout_check(tm, time)` matches within `1e-6` (`:12,33-67`); always includes `time_end` | linear scan; `nmax=100000` read through `nml_read` of a 100000-element vector (25.6 MB automatic string array on the stack, `nml.f90:857`) |
| `timer` (`timer.f90:8-26`) | CPU/wall profiling per numbered component, table output | profiling only |
| `tsgen` | scalar forcing `f(t)` incl. feedback controllers, with restart | forcing, not scheduling |
| `phys_constants` | named year lengths (`sec_year_360d ... _tropical`, `:62-78`) | definitions only |

Elsewhere:

- yelmox: `cadence_due(time, dt) = mod(nint(time*100), nint(dt*100)) == 0`
  (`yelmox/libs/kryos.f90:215-222`), used for `dt_clim`, `dt_isos`, `dt_cmip`. Cadences are measured
  from absolute time 0, have 0.01-yr resolution, and `nint` on default integers overflows past
  |t| ~ 21 Myr.
- climber-x: its own `timer` module is a full integer-step multi-cadence scheduler (per component
  `nstep_*`, `dt_*`, and `time_call_*`, `time_soy_*`, `time_eom_*` flags held as module globals,
  `climber-x/src/main/timer.f90:40-100, 263-575`), with a 360-day calendar and asynchronous
  acceleration (`n_accel`, `n_year_geo`, `n_year_ice`).

**There is no model of "component X updates every dt_X, last updated at t_X"** in fesm-utils. Gaps for
multi-cadence scheduling:

1. Integer time base. Floating-point modulo/tolerance matching (`timeout_check`, `cadence_due`) is
   fragile; climber-x's integer steps are the robust model. Use an integer tick (e.g. 1 tick = 1 day,
   or a user-defined `tick` in seconds) in `integer(int64)`.
2. Per-component state: `dt`, phase/offset, `t_last`, `t_next`, and whether the component is
   "leading" (integrates over `[t, t+dt]`) or "lagging" (`[t-dt, t]`); this determines when
   exchanged fields are valid.
3. Interval semantics for exchanged fields: instantaneous vs mean over the last interval vs
   accumulated (fluxes). Needed when a fast component (SMB, monthly) feeds a slow one (isostasy).
4. Calendar: sub-annual cadences (monthly SMB, daily climate) need a calendar choice (360/365/365.25 d),
   which `phys_constants` names but nothing uses.
5. Restart alignment: on restart the scheduler must restore `t_last` per component (yelmox comments
   on re-deriving climate depending on whether the restart time is on the `dt_clim` cadence,
   `kryos_startup.f90:69-80`).

---

## 4. I/O

### 4.1 Reading fields along time: varslice

- One `varslice_class` per (file pattern, variable). Namelist keys: `filename` (with `{domain}`,
  `{grid_name}`, `{key}` substitutions and globbing), `name`, `units = in out` (metadata only),
  `scaling = scale offset`, `time = active x0 x1 dx nsub`, optional `fill`, `fill_value`, `remap`,
  `grid_src` (`varslice.f90:1755-1870`, `docs/packages/varslice.qmd:160-205`).
- `varslice_update(vs, time, method, fill, rep)` (`:430`) returns slices by `exact`, `range`, `interp`,
  `extrap`, `range_mean/min/max/sd/sum`, with `rep` for sub-annual cycles; caches the last slice; reads
  across multiple files (`nc_read_multifile`, `:2161`); online remap to the target grid at init
  (`varslice_remap_init`, `:211`).
- Limits for a generic data backend: the file's time axis is **reconstructed from the namelist**
  (uniform `dx`), not read from the time variable; data are `real(wp)` only; one variable per object;
  units are not read from the file (no check against what the consumer expects); no calendar-aware
  time (years only); remap is fixed at init (fine).

### 4.2 Locating datasets: fesmdata/topodata/regions

- FesmData v2 layout `$FESMDATA/<Domain>/<GRID>/` with a global attribute `grid_name` and a
  `grid_<GRID>.txt` next to the file; `fesmdata_read_field` remaps if the file grid differs from the
  target (`fesmdata.f90:1-16, 414-630`). Path templating `{domain}`, `{grid_name}`, `{key}`
  (`fesmdata_parse_path`, `:634`).
- There is no catalogue/registry of datasets (by variable/standard name), only path templates.

### 4.3 Writing

- `ncio` writes everything; output per step is hand-written per component (`yelmox/libs/kryos_output.f90`
  has `isos_write_step`, `mshlf_write_step`, `smb_write_step`, ... and `YELMO_VARS_2D` lists).
- `variable_io` provides table-driven metadata (markdown table with id, name, dims, units, long_name,
  prec) and `var_io_write` generics (`variable_io.f90:8-37, 415-663`). Only yelmo uses it. This is
  the seed of a metadata-driven writer.
- Output precision: `prec_abs/prec_rel` bit rounding for compression, `nc_set_lossless` for restarts
  (recent, `git log` 4020c7f, ee23a19).

### 4.4 Restarts

No restart convention in fesm-utils except `tsgen_restart_write/read` (`tsgen.f90:731-846`, including
RNG state). yelmox defines the convention: a **bundle folder per time**, one restart file per stateful
component (`yelmo`, `isos`, `bsl_restart.nc`, `marshelf`, `smbpal`, `chion`, climate, `coef_<name>`),
written by `run_restart_write`/`domain_restart_write` (`yelmox/libs/kryos_startup.f90:1-4, 106-129`).

### 4.5 Would these support a generic "data backend"?

Partly. `varslice` already does the hard part (file globbing, time slicing, sub-annual cycles,
online remap, fill). A data backend for *any* component needs on top:

- a mapping from the component's **export field names** to varslice groups (one namelist group per
  field, or a table);
- time-axis reading from the file (and calendar/unit handling), not namelist reconstruction;
- units read from the file and converted to the field's declared units;
- serving into the exchange state (§6.2) with the right grid tag and time stamp;
- handling of "static" fields (read once), "climatology" (cycle), and "transient" uniformly (varslice
  already distinguishes these).

This is a thin layer (a few hundred lines) and fits well in fesm-utils as `datasrc` (§6.5).

---

## 5. Parameters: nml

Conventions (`nml.f90`, `docs/packages/nml.qmd`):

- Fortran-namelist-looking files parsed by hand (not `read(nml=)`), so groups can contain unknown keys
  and values can be read one key at a time: `call nml_read(file, group, name, value [, comment, init,
  defaults_file, defaults_group])`.
- Two-layer configuration: `defaults_file` is the schema; every key must exist there; user file
  overrides (`nml_read_internal`, `:173-240`). `defaults_group` maps a user group to a schema group
  (yelmo group aliasing, e.g. `ydyn_north -> ydyn`). `nml_validate` reports user keys not in the
  schema (`:634-679`).
- Flat syntax `group.param = value` (several per line) for command-line-like overrides (`:346-460`).
- `nml_replace(s, text, rep)` for string substitution (path templating).
- Echo of every value read when `VERBOSE` (module global, `:7`).

Limits:

| Limit | Where |
|---|---|
| file scanned from the top for **every key** (O(keys x lines)) | `nml_find_param`, `:244-341` |
| max 5000 lines per file, 1000 chars per line, no continuation lines | `:283-284, 325` |
| value string truncated to 256 chars in typed readers, so long vectors are silently cut | `:691-857` |
| vector reader allocates `size(value)` 256-char temporaries on the stack | `:857` (timeout: 100000 elements) |
| max 32 flat assignments per line, 2000 keys in `nml_validate` | `:122, 642` |
| no `n*value` repeat syntax, no arrays of derived types | parser |
| module-global switches (`VERBOSE`, `ERROR_NO_PARAM`) | `:7-8` |

Consumers build their own conventions on top: yelmox `yelmox_nml` has a variable-definition table
(`VAR_DEFS_MAX`, `var_defs_find`, defaults file); kryos groups are `[domain]`, `[comps]`, `[sim]`,
`[output]` with per-domain suffixes (`kryos.f90:56-133`).

---

## 6. What is missing: candidates to add

Guiding principle: fesm-utils should keep providing **mechanisms with no knowledge of specific
components**. Anything that encodes *which* components exist, or the coupling sequence, belongs in
kryos. A candidate goes into fesm-utils when it is (a) model-agnostic, (b) useful to more than one
host (kryos, climber-x, VILMA2 standalone), and (c) small enough to keep the library dependency-light.

| Candidate | Exists partially? | Verdict |
|---|---|---|
| 6.1 Abstract component base type | yelmox `kryos_domain` (concrete), climber-x `*_model` init/update/end conventions | **kryos first**; promote later if climber-x adopts it |
| 6.2 Field / exchange state | no; `coupler` gives grid tags, `variable_io` metadata | **add to fesm-utils** (generic container) |
| 6.3 Global-domain abstraction | `coupler`, VILMA2 `vilma_remap`, driver-owned `bsl` | **mechanisms in fesm-utils** (merge/split/fraction), topology (which domains) in kryos |
| 6.4 Scheduler | yelmox `cadence_due`, climber-x `timer` | **clock + cadence type in fesm-utils**, ordering policy in kryos |
| 6.5 Data backend helper | `varslice`, `fesmdata` | **add `datasrc` to fesm-utils** on top of varslice |
| 6.6 Restart bundle manager | yelmox `kryos_startup`, `tsgen_restart_*` | small **helper in fesm-utils** (paths + manifest); iteration over components in kryos |
| 6.7 Output writer from registry | `variable_io`, `timeout`, `ncio` | **add to fesm-utils** after 6.2 exists |
| 6.8 Physical constants | `phys_constants` (good) | **exists**; adopt everywhere, extend fields |
| 6.9 Logging / errors | 10 private `*_error` copies | **add now** (cheapest, highest leverage) |
| 6.10 Units | `fesmdata_flux_conv/temp_ref`, varslice scale/offset | **small canonical-units check only**; not a udunits port |

### 6.1 Abstract component base type

Partial today: yelmox components share a naming convention (`*_init`, `*_update`, `*_restart_write`,
`*_restart_read`, `*_write_init`, `*_write_step`) but are concrete types held side by side in
`kryos_domain` (`kryos.f90:163-184`), with backend selection via strings and `select case`
(`ctl%climate`, `ctl%surface_method`, `kryos.f90:70-71`).

Sketch (Fortran 2008):

```fortran
module kryos_component
    use fesm_precision, only: dp
    use fesm_fields,    only: field_state_t
    implicit none
    private

    type, abstract, public :: component_t
        character(len=64) :: name    = ""       ! instance name, e.g. "isos_north"
        character(len=64) :: kind    = ""       ! role: "ice_sheet", "isostasy", "smb", ...
        character(len=64) :: backend = ""       ! "yelmo", "pism", "data", ...
        character(len=64) :: grid    = ""       ! native grid tag (coupler name)
        type(field_state_t) :: imports, exports
    contains
        procedure(comp_init),    deferred :: init
        procedure(comp_update),  deferred :: update
        procedure(comp_rst_w),   deferred :: restart_write
        procedure(comp_rst_r),   deferred :: restart_read
        procedure(comp_out),     deferred :: write_output
        procedure                         :: finalize => comp_finalize_default
    end type

    abstract interface
        subroutine comp_init(self, par_file, group, time)
            import :: component_t, dp
            class(component_t), intent(inout) :: self
            character(len=*),   intent(in)    :: par_file, group
            real(dp),           intent(in)    :: time
        end subroutine
        subroutine comp_update(self, t0, t1)        ! integrate over [t0, t1]
            import :: component_t, dp
            class(component_t), intent(inout) :: self
            real(dp),           intent(in)    :: t0, t1
        end subroutine
        subroutine comp_rst_w(self, dir, time)
            import :: component_t, dp
            class(component_t), intent(in) :: self
            character(len=*),   intent(in) :: dir
            real(dp),           intent(in) :: time
        end subroutine
        subroutine comp_rst_r(self, dir, time)
            import :: component_t, dp
            class(component_t), intent(inout) :: self
            character(len=*),   intent(in)    :: dir
            real(dp),           intent(out)   :: time
        end subroutine
        subroutine comp_out(self, file, time, init)
            import :: component_t, dp
            class(component_t), intent(inout) :: self
            character(len=*),   intent(in)    :: file
            real(dp),           intent(in)    :: time
            logical,            intent(in)    :: init
        end subroutine
    end interface

    type, public :: component_ptr                 ! heterogeneous collection
        class(component_t), allocatable :: c
    end type
end module
```

The driver holds `type(component_ptr), allocatable :: comps(:)` and allocates concrete types via a
factory (`select case (backend) ; case ("yelmo") ; allocate(yelmo_component_t :: comps(i)%c)`).

Critical notes:

- Data exchange should **not** go through type-bound getters per field (that recreates hard-coded
  coupling); use the import/export states (6.2). The component only declares names and fills/reads
  them in `update`.
- Global vs regional is a property of the component's grid tag plus a `scope` attribute
  (`"regional"`/`"global"`), not a different base type.
- Compiler risk: abstract types, deferred bindings, `class(...), allocatable` components and
  `allocate(T :: x)` are solid in gfortran >= 9, ifort >= 19 and ifx; avoid parameterized derived
  types, `class(*)` in hot paths, polymorphic arrays (`class(component_t), allocatable :: c(:)`
  forces one dynamic type; use the wrapper), and recursive allocatable components with finalizers
  (`final` on types with allocatable polymorphic components has had gfortran bugs up to ~12).
  Submodules for large backends are fine on all three.
- Performance: one virtual call per component per step is negligible.
- Why kryos first: the concrete backends (Yelmo, FastIsostasy, VILMA2, data) are wrapped in kryos
  anyway; climber-x would need large refactoring to adopt it. Promote `component_t` to fesm-utils
  when a second host uses it.

### 6.2 Field registry / exchange state

Nothing exists. Today coupling is hard-coded `couple_*_to_yelmo` and `step_*` routines
(`yelmox/libs/kryos_coupling.f90:27-31`), each naming source/target grids and remap methods.

Sketch:

```fortran
module fesm_fields
    use fesm_precision, only: dp
    implicit none
    integer, parameter :: FIELD_INSTANT = 1, FIELD_MEAN = 2, FIELD_ACCUM = 3   ! time semantics
    integer, parameter :: CONS_NONE = 0, CONS_EXTENSIVE = 1                     ! conserve integral?

    type :: field_t
        character(len=64)  :: name  = ""     ! standard short name, e.g. "H_ice", "smb", "dz_bed"
        character(len=32)  :: units = ""     ! canonical units string, e.g. "m", "m yr-1"
        character(len=64)  :: grid  = ""     ! coupler grid tag
        integer            :: tsem  = FIELD_INSTANT
        integer            :: cons  = CONS_NONE
        integer            :: nz    = 1      ! 1 = 2D field
        real(dp)           :: t0 = -huge(1._dp), t1 = -huge(1._dp)   ! validity interval
        real(dp), allocatable :: v(:,:,:)     ! owned buffer, dp
        logical,  allocatable :: valid(:,:)   ! optional coverage
        real(dp)           :: wsum = 0._dp    ! accumulated time weight (FIELD_MEAN/ACCUM)
    end type

    type :: field_state_t
        type(field_t), allocatable :: f(:)
        integer :: n = 0
    contains
        procedure :: add       ! (name, units, grid, nx, ny, nz, tsem, cons) -> handle
        procedure :: handle    ! (name) -> integer index, error if absent (resolve at init)
        procedure :: put       ! (h, var, t0, t1) generic dp/sp 2D/3D, copies + accumulates
        procedure :: get       ! (h, var) generic dp/sp, copies out (and resets accumulators)
        procedure :: has       ! (name) -> logical
    end type
end module
```

Plus in kryos (or fesm-utils `coupler`), a connection list built at init from configuration:
`(src_comp, src_field) -> (dst_comp, dst_field, method, norm)`, executed by a generic
`exchange(conn, cpl)` that calls `remap`. This replaces `couple_*` code with data.

Risks and choices:

- **Copies vs pointers**: pointers into component-owned arrays (needs `target`, fragile when
  components reallocate, and mixes sp/dp). Owned `dp` buffers with copy-in/copy-out cost a memory
  copy per exchange, negligible versus a remap, and make sp/dp hosts interoperate. Recommend copies.
- **By-name lookups**: resolve names to integer handles at init; per-step cost is then an array index.
  Linear scans over tens of names are fine even unresolved.
- Avoid `class(*)` payloads; fix `real(dp)` 3D (`nz=1` for 2D) and provide `put/get` generics for
  `sp`/`dp` 2D/3D and `integer` masks (stored as dp or a separate integer buffer).
- Keep the vocabulary of names in kryos (a table like `variable_io`), not hard-coded in fesm-utils.

### 6.3 Global-domain abstraction

Mechanisms to add to fesm-utils (extending `coords`/`coupler`):

```fortran
! coverage fraction and normalization choice in the core apply
call map_field(map, name, var1, var2, ..., frac2=frac, norm="dst")   ! norm: "frac" (today) | "dst" | "none"

! merge N regional contributions into one global field, conservatively
subroutine coupler_merge(cpl, dst, var_dst, srcs, vars_src, owner, background, method)
    type(coupler_class), intent(inout) :: cpl
    character(len=*),    intent(in)    :: dst            ! global grid tag
    real(dp),            intent(inout) :: var_dst(:,:)   ! in: background; out: merged
    character(len=*),    intent(in)    :: srcs(:)        ! regional grid tags
    type(field_ref),     intent(in)    :: vars_src(:)    ! wrappers around regional arrays
    real(dp), optional,  intent(in)    :: owner(:,:,:)   ! (nx_src, ny_src, k) weights in [0,1], resolve overlaps
    character(len=*), optional, intent(in) :: method     ! default "con"
end subroutine
! var_dst = (1 - sum_k frac_k) * background + sum_k frac_k * mean_k,   frac_k from norm="dst"

! scatter a global field to each regional grid (anomaly form optional)
subroutine coupler_split(cpl, src, var_src, dsts, vars_dst, ref_src, ref_dst, method)

! integral of a field on a registered grid (consistent area metric), for conservation checks
function coupler_integral(cpl, grid, var, mask) result(total)
```

plus: Gaussian cell bounds from Gauss weights (or `grid_init(..., y_bnds=)`), a grid hash in the map
cache key, `GRID_MAX`/`MAP_MAX` removed, identity shortcut and mask/fill pass-through in `remap`.

What belongs in kryos: the **topology** (one global "earth" domain + N regional domains, which
global components exist, ownership priorities), and the order of merge/solve/split in the coupling
step. Do not create a `global_domain_t` in fesm-utils; a global component is just a component whose
grid tag is a global grid.

Risks: overlap resolution needs explicit owner weights (two regional domains covering the same global
cell); conservative maps from fine regional grids to a coarse Gauss grid are cheap with caching but
the spherical-target clip is the slowest path (`NSUB=8` subsampling, `conservative.f90:32`);
`sp` hubs lose precision in volume sums.

### 6.4 Scheduler for multi-cadence components

Add the *clock* to fesm-utils; keep the *order* of calls in kryos.

```fortran
module fesm_clock
    use, intrinsic :: iso_fortran_env, only: int64
    use fesm_precision, only: dp
    type :: clock_t
        integer(int64) :: tick_s = 86400_int64     ! seconds per tick (1 day default)
        integer(int64) :: n0, n, n_end             ! ticks
        character(len=16) :: calendar = "365_day"  ! 360_day | 365_day | julian | none
    contains
        procedure :: time_yr    ! current time in years (dp)
        procedure :: advance    ! (nticks)
    end type
    type :: cadence_t
        character(len=64) :: name = ""
        integer(int64) :: dt     = 0               ! ticks; 0 = never
        integer(int64) :: offset = 0               ! phase relative to n0
        integer(int64) :: n_last = -huge(1_int64)
    contains
        procedure :: due        ! (clock) -> logical : mod(n - n0 - offset, dt) == 0
        procedure :: mark       ! (clock) : n_last = n
        procedure :: interval   ! (clock) -> [t0, t1] years since last call
    end type
end module
```

`tstep_class` can stay as the user-facing "time view" (cal/rel/const methods) on top of the integer
clock. Risks: converting existing yr-based namelists (`dt_isos = 10.0`) to ticks requires rounding
rules; asynchronous acceleration (climber-x `n_accel`) is policy and should stay in the host.

### 6.5 Generic data-backend helper

```fortran
type :: datasrc_t                         ! serves named fields from files
    type(varslice_class), allocatable :: vs(:)
    character(len=64),    allocatable :: names(:)     ! export field names
    character(len=16),    allocatable :: mode(:)      ! static | clim | transient
contains
    procedure :: init     ! (par_file, group_prefix, names(:), grid) -> one varslice per name
    procedure :: update   ! (time or [t0,t1]) -> varslice_update with the right method
    procedure :: fill     ! (state: field_state_t) -> put each field with its validity interval
end type
```

A kryos "data" backend for any component is then `component_t` + `datasrc_t`. Prerequisites in
varslice: read the time axis from the file (keep the namelist form as an override), `dp` storage or a
`dp` accessor, units read from the file and checked against the field's declared units.

### 6.6 Restart bundle manager

Small, generic helper (the yelmox convention made reusable):

```fortran
function restart_bundle_dir(root, time, fmt) result(dir)      ! e.g. root/restart-<time>/
subroutine restart_manifest_write(dir, time, entries)         ! component, backend, file, grid, git hash
subroutine restart_manifest_read(dir, entries, time)
```

The loop "for each component call `restart_write(dir, time)`" is kryos. Use `nc_set_lossless` for
restart files (exists). Include scheduler state (`n_last` per cadence) and exchange-state
accumulators in the bundle; tsgen already includes RNG state, which is the standard to follow.

### 6.7 Standard output writer driven by the registry

Combine `variable_io` (metadata table), `timeout` (when), `ncio` (how) and the field state (what):
`call outfile_init(of, file, grid, names(:), table)`, `call outfile_write(of, state, time)`. Lets each
component's output be a configuration list instead of a hand-written `*_write_step` (compare
`kryos_output.f90:260-401`). Backend-internal diagnostics remain the backend's own business.

### 6.8 Physical constants

Already right: one `dp` record, values from a namelist (`par/phys_const_earth.nml`), loud failure on
unset, `phys_const_compare` to check two components agree. Missing: planet radius/ellipsoid (lives in
`coords/planet.f90`), mantle density is there (`rho_asth`), ocean area (`area_seasurf`); add
`R_earth`, `rho_mantle` aliases only if a consumer needs them. The work is adoption (VILMA2's
`vilma_constants`, climber-x `constants`), not new code.

### 6.9 Logging, errors, diagnostics

Add `fesm_log` (one module, replaces the 10 private helpers):

```fortran
call log_init(level=LOG_INFO, unit=output_unit, prefix="")      ! per-run
call log_msg(LOG_INFO, "isostasy", "updated", detail)            ! levels: DEBUG/INFO/WARN/ERROR
call fesm_error(proc, msg, detail [, stat])                       ! error stop 1 unless stat present
interface to_str ; ... int/sp/dp/logical, 1D ; end interface
```

Route errors to `error_unit`, honour a global verbosity (replacing `nml`'s `VERBOSE`), optionally tag
lines with component/domain. Diagnostics (global integrals, conservation residuals per exchange) can
log through it. Risk: low. This is the cheapest, highest-leverage addition, and the team's own rule
("a fifth module") is already exceeded.

### 6.10 Units

Partial: `fesmdata_flux_conv` (water flux mass/time units, 1 yr = 365.25 d), `fesmdata_temp_ref`
(K/degC), varslice `scale/offset`. Recommend **not** porting udunits. Instead: canonical unit strings
on every exchanged field, a normalizer (`units_key` already exists, `fesmdata.f90:334`), and a short
table of known conversions (m <-> km, yr <-> s via the chosen calendar, kg m-2 <-> m w.e. <-> m i.e.
via `phys_constants`). Exchange fails at init if units differ and no conversion is known. Keeps the
scope small; full dimensional analysis is not worth it in Fortran.

### Also worth adding

- **Grid identity**: `grid_hash(grid)` (cs + axes) for caches and restart checks.
- **Identity remap** in `coupler` (src == dst).
- **`make check` + CI** that builds both variants and runs the test programs (parked in the
  sustainability plan); required before downstream depends on new framework pieces.

---

## 7. Conventions worth standardizing

| Topic | Today | Recommendation |
|---|---|---|
| Kinds | `precision` has `dp, sp, wp=sp`; downstream redefine `wp` (sp, dp) | library exports `sp`, `dp` only in public interfaces; generic `sp`/`dp` overloads as in `map_field`/`coupler`; exchange buffers `dp`; each component keeps its own `wp` internally. Never let a library interface depend on `wp` |
| Module names | generic names (`precision`, `constants`, `timer`, `coupler`, `index`, `planet`, `mapping`, `regions`, `series`) collide downstream | prefix new modules (`fesm_*` or `fu_*`); add prefixed aliases for old ones, deprecate the bare names |
| Type names | `*_class` (fesm-utils, yelmo), `*_t` (`weight_map_t`), `*_type` (`var_io_type`) | pick one (`_t` or `_class`) for new code |
| Procedure names | `<pkg>_init/_update/_end`, `_restart_write/_read`, `_write_init/_write_step` (htopo, zremap, varslice, yelmox components) | make these the component contract (6.1); `_end` vs `_finalize` choose one |
| Errors | mix of `stop`, `error stop`, 10 private `*_error` helpers; 413 `write(*,*)` vs 137 `error_unit` | one `fesm_error(proc,msg,detail[,stat])`; errors to `error_unit`; `error stop 1` |
| Namelist | `&group` + `defaults_file` schema + `nml_validate`; per-domain suffix groups (`regions_north`) | defaults file per component backend (`<backend>_defaults.nml`), group = component instance name, validation mandatory at init |
| Time | `real(dp)` years; `time_cal` vs `time_rel` (ref 1950) | integer ticks internally (6.4); years `dp` at the interface; document "yr" = which calendar |
| Units | SI-ish, yr for time, model-specific (m ice equivalent vs w.e.) | declare canonical units per exchanged field (ice thickness m, SMB m i.e. yr-1, temperatures K, uplift m, mass Gt only in diagnostics) |
| Missing values | `mv=-9999` (`constants.f90:16`), `regions_undefined=-1`, netCDF fill reset above `1e10` in varslice | keep `mv`; exchanged fields carry a `valid` mask instead of in-band missing values |
| Masks | integer codes per model (`topo_mask_ocean=0, land=1, grnd=2, flt=3` in `topodata.f90:33-36`; yelmo `mask_bed_*` and `MASK_ICE_NONE/FIXED/DYNAMIC`); logical masks from `regions` | exchange **fractions in [0,1]** (`f_grnd`, `f_ice`, `f_ocn`) and logical validity; integer codes stay inside components or carry a CF flag table |
| Array layout | `(nx, ny[, nz])`, x fastest, latitude ascending (coords and VILMA2 Gauss grid both assume ascending) | keep; document; flip at backend boundaries |
| Optional args | widespread `present(x)` defaults | fine; prefer keyword calls in drivers |

---

## 8. Other things a from-scratch design should know

1. **Module-name collisions and `-zmuldefs`.** climber-x ships `module precision` with `wp = dp`
   (`climber-x/src/utils/precision.f90:37`), `module constants`, `module timer`, `module coupler`;
   FastIsostasy ships `module precision` (`yelmox/FastIsostasy/src/precision.f90:1`). Which `.mod`
   wins depends on include order, and all Linux builds link with `-Wl,-zmuldefs`, explicitly to work
   around "duplicate symbols in the static deps" (`climber-x/config/common.mk:139-142`,
   `yelmox/config/common.mk:91-95`, `fesm-utils/config/common.mk:20-22`), which hides such clashes. kryos and new
   fesm-utils modules should use unique prefixes. Note also that yelmox already contains modules named
   `kryos`, `kryos_coupling`, `kryos_output`, ... (`yelmox/libs/kryos*.f90`): the new repo must not
   reuse these module names if both are ever built together.
2. **yelmox `kryos*` is the prototype** (2875 lines). Reusable ideas: component grids as strings in
   configuration (`grid_ice`, `grid_isos`, `grid_clim`, ...), one coupler per domain, the hi-res
   hub as geometry source of truth, `cadence_due`, restart bundles, the driver-owned global object
   (`bsl`), the coefficient registry (`kryos_forcing.f90` `coefs_class`: named scalar forcings with
   derived specs and restart). What it lacks is abstraction: concrete types, backend `select case`,
   coupling written as code.
3. **Two parallel time models** exist (fesm-utils `timestepping` float years; climber-x integer-step
   `timer` with calendar). kryos needs to pick one; integer ticks are the safer base.
4. **Sp/dp mix is real.** Yelmo is sp, VILMA2 dp, climber-x dp, htopo sp. Design exchange in dp.
5. **The conservative path is now in-house**; `cdo` is optional (`gen="cdo"`) but still shelled out
   via `execute_command_line` for that option and for some grid-description utilities
   (`grid_cdo_write_desc_via_cdo`, `call_system_cdo`, `coords/grid_cdo.f90:624,688`).
6. **Deep copies.** `coupler_add_grid` and `resolve_grid` copy whole `grid_class` objects with five
   `nx*ny` dp arrays; assigning a `kryos_domain` deep-copies every component. Keep large objects in
   place (allocate once, pass by reference, avoid intrinsic assignment of containers).
7. **OpenMP variants.** The library is built twice (`include-serial`, `include-omp`);
   `weight_map_apply` and conservative weight generation are OpenMP-parallel. A host must link the
   variant matching its own build.
8. **No MPI anywhere.** Everything is single-process (OpenMP at most). If several regional domains
   should run concurrently, that is a new design dimension (MPI communicators per component); keep
   the component interface free of global state so it remains possible.
9. **Build and packaging.** configme-generated Makefile, no `fpm`/CMake, no Fortran CI. `ncio.f90` is
   generated from `ncio.fypp` (`make ncio-gen`); edit the `.fypp`.
10. **Global state in modules.** `nml` (`VERBOSE`, `ERROR_NO_PARAM`), `ncio` (`nc_lossless`), climber-x
    `timer`/`control` globals. A multi-domain framework needs per-instance state; avoid new module
    globals.
11. **Documentation is good and current** (Quarto site per package, dev notes with decisions). Follow
    the same pattern for kryos (design notes with numbered, settled decisions).
12. A Julia mirror exists (`/Users/alrobi001/models/FesmUtils.jl`); not studied here, but relevant if
    kryos tooling (analysis, grid generation) is written in Julia.
