# Backend study: what each candidate backend needs and provides

Status: study (2026-10-24). Input for the kryos component interfaces.
Coupling as done today in yelmox and CLIMBER-X is studied separately; this
document looks at the backends themselves.

States read: yelmo (local `src/`), CISM `cism_main_3.00.001` (local),
FastIsostasy (local), vilma2 (`climber-x/vilma2`), chion, rembo1, insol,
yelmox `libs/` (smbpal, smb_simple, marine_shelf, pico, obm, geothermal,
sediments), FastHydrology (`yelmo/FastHydrology`), fesm-utils (local).
PISM is not available locally; its section is from memory and is marked
as such. The two yelmox plans `docs/_cism-backend-plan.md` and
`docs/_vilma2-backend-plan.md` were read and are referenced, not repeated.

Paths are relative to each repository root unless written in full.

## Summary

- Every backend except rembo1 and CISM's GLAD layer is a derived-type
  object with `init / [init_state] / update / restart_write / end`. That is
  the shape kryos should standardise on.
- Clocks differ in three ways: **advance to a target time** (Yelmo,
  FastIsostasy, rembo1, FastHydrology, PISM `run_to`), **advance by dt**
  (VILMA2, chion, obm), **one fixed step per call** (CISM
  `glissade_tstep`, smbpal annual). `update(comp, time_target)` with an
  adapter loop covers all of them.
- Only VILMA2 and obm are global. Every other backend is regional.
  FastIsostasy is regional with a shared scalar (barystatic sea level)
  that acts as the only global quantity today.
- Precision is mixed: Yelmo, fesm-utils, FastIsostasy, chion, smbpal,
  marine_shelf, obm are `sp`; VILMA2, CISM, rembo1, PISM, insol (default)
  are `dp`. Time is `dp` everywhere except VILMA2 (relative `wp=dp`) and
  rembo1. The framework must convert at the adapter boundary.
- Units and signs disagree on SMB (m ice/yr, m ice/s, mm w.e./yr,
  kg m-2 s-1), surface temperature (K vs degC), sub-shelf melt sign (Yelmo
  and marine_shelf: melt < 0; CISM, PISM: melt > 0), geothermal flux (Yelmo
  mW m-2 up; CISM W m-2 positive down; PISM W m-2 up), the year length
  (360 d, 365 d, 365.2422 d) and ice density (910, 917, 931).
- Instantiability: Yelmo, FastIsostasy, chion, smbpal, marine_shelf,
  VILMA2 can run several instances. rembo1 cannot (module globals, grid
  fixed at compile time). CISM can, up to 8, sharing one communicator.
  PISM probably can, behind PETSc's global initialisation.

## 1. Ice-sheet models

### 1.1 Yelmo

| Item | Yelmo |
|---|---|
| Role | ice sheet (thermomechanics, calving, basal hydrology via FastHydrology, tracers) |
| Grid | regular projected x-y (`grid_class` from fesm-utils coords, meters; `yelmo_grid.f90:126,246,267`); scalars on aa (cell centres), velocities on ac (C-grid); regional |
| Clock | owns its clock; `yelmo_update(dom,time)` sub-steps adaptively (predictor-corrector) until `time` (`src/yelmo_ice.f90:42`). Idempotent at `time == current` (dt = 0). `yelmo_set_time` rewinds all sub-clocks (`:1822`) |
| Time units | years, `real(dp)`, absolute model time (may be negative) |
| Precision | `wp = sp`, hard-coded (`src/yelmo_defs.f90:42`; `docs/numerics-precision.md`). No build switch |
| Params | namelist groups `yelmo`, `ytopo`, `ydyn`, `ymat`, `ytherm`, `yelmo_init_topo`, `yelmo_masks`, ... with `input/yelmo_defaults.nml` as schema; group name is an argument of `yelmo_init` |
| Constants | optional `cnst` (`phys_const_class`) at `yelmo_init`; otherwise `input/yelmo_phys_const.nml` (`sec_year = 31556926`, `rho_ice = 910`) |
| Restart | own netCDF file, `yelmo_restart_write(dom,file,time)` (`src/yelmo_io.f90:301`); read via `yelmo_restart_init` (`yelmo_ice.f90:1602`) or `restart=` in namelist |
| OpenMP / MPI | OpenMP, no MPI |

**Entry points** (`src/yelmo_ice.f90`, re-exported by `src/yelmo.f90`):

| Procedure | Signature (abridged) |
|---|---|
| `yelmo_init` (`:685`) | `(dom, filename, grid_def, time, load_topo, domain, grid_name, group, outfldr, cnst, reg, mask_ice, topo_pd, topo_init)` — grid, parameters, allocation; driver may pass constants, regions, ice mask and topography instead of file reads |
| `yelmo_init_topo` (`:1047`) | `(dom, filename, group, time, load_topo, topo_init)` — step 1: topography only, so that boundary fields (T_srf, smb) can be computed from it |
| `yelmo_init_state` (`:1473`) | `(dom, time, thrm_method, T_ice)` — step 2: temperature, material, dynamics, after `bnd` is filled |
| `yelmo_update` (`:42`) | `(dom, time)` |
| `yelmo_update_equil` (`:603`) | `(dom, time, time_tot, dt, topo_fixed, tpo_solver, dyn_solver)` — spin-up helper |
| `yelmo_restart_write` / `yelmo_restart_init` | see above |
| `yelmo_end` (`:1775`) | `(dom, time)` |

**Inputs** (`dom%bnd`, `type ybound_class`, `src/yelmo_defs.f90:887`; units
from `input/yelmo-variables-ybound.md`). All 2D on the Yelmo grid; the
driver writes them before each `yelmo_update`.

| Field | Units | Notes |
|---|---|---|
| `z_bed` | m | bedrock elevation; rate `dz_bed_dt` computed internally between calls |
| `z_bed_sd` | m | subgrid roughness (static) |
| `z_sl` | m | **2D** sea level |
| `H_sed` | m | sediment thickness (static) |
| `smb` | m ice/yr | + = gain |
| `T_srf` | K | ice-surface temperature |
| `bmb_shlf` | m ice/yr | **melt < 0** (net `mbal = smb+bmb+fmb+calv`, `physics/mass_conservation.f90:247`) |
| `fmb_shlf` | m ice/yr | frontal melt |
| `T_shlf`, `tf_shlf` | K | shelf temperature, thermal forcing |
| `Q_geo` | mW m-2 | positive into the ice |
| `Qd` | (Yelmo's units) | subglacial discharge for frontal melt (`yelmo_topography.f90:1068`) |
| `mask_ice` | 0/1/2 | where ice may exist / is fixed / is dynamic (`yelmo_defs.f90:69-71`) |
| `H_ice_ref`, `z_bed_ref`, `tau_relax`, `calv_mask`, `basins`, `regions` | | config-like, mostly static |
| `dzbdt_corr`, `z_bed_corr` | | internal correction for low-res isostasy |

**Outputs** (`dom%tpo%now`, `type ytopo_state_class`, `yelmo_defs.f90:255`;
plus `dyn`, `thrm`):

| Field | Units | Notes |
|---|---|---|
| `H_ice`, `z_srf`, `z_base` | m | |
| `f_grnd`, `f_ice` | 0-1 | subgrid grounded and ice-covered fractions |
| `mask_bed` | int | 0 ocean, 1 land, 2 frozen, 3 stream, 4 grline, 5 float, 6 island, 7 partial (`physics/topography.f90:12-19`) |
| `dHidt`, `mb_net`, `smb`, `bmb`, `fmb`, `cmb`, `dmb` | m ice/yr | applied mass-balance terms over the last call |
| `dyn%now%ux_s, uy_s, uxy_b`, ... | m/yr | velocities (ac nodes) |
| `thrm%now%T_ice` | K | 3D |
| `hyd%now%N, W_til, Q_b` | | basal hydrology |
| `reg` / `regs` | | regional integrals (V_ice, ...) |

**Pain points for embedding**
- `wp = sp` fixed; a `dp` framework converts at the boundary.
- Two-step init (`init_topo` → driver fills `bnd` → `init_state`) is
  needed because `T_srf`/`smb` depend on initial topography. Any
  framework must allow a component to expose geometry before its full
  state exists.
- Yelmo carries its own high-resolution reference bed; isostasy only adds
  anomalies (`z_bed = z_bed_ref + dz`). Interfaces should pass anomalies
  or the full bed explicitly — not both implicitly.
- Many yelmox features reach into Yelmo internals (opt/`cb_ref`, NEGIS,
  equilibrium spin-up, regions output). See `_cism-backend-plan.md`.
- No `save` module state in `src/`; module variable `yelmo_use_omp`
  (`yelmo_defs.f90:93`) is harmless. Several `yelmo_class` objects can
  coexist (yelmox_bipolar does). The C API (`src/yelmo_c_api.f90:8-10`)
  hard-wires two global instances `ylmo1`/`ylmo2` selected by alias
  string.

**C API** (`src/yelmo_c_api.f90`, `docs/c-api.md`): `yelmo_init_grid`,
`yelmo_init`, `yelmo_init_state`, `yelmo_step(time)`, `yelmo_restart_write`,
`yelmo_get_grid_sizes/info`, `yelmo_get_var2D/3D`, `yelmo_set_var2D/3D`
by string name (`"bnd_z_bed"`, `"bnd_smb_ref"`, `"tpo_H_ice"`, ...;
`:618-690`), plus a callback for `N_eff`. Arrays cross as `c_double`
and are converted to `wp`. This string-keyed get/set is a useful
precedent for a language-neutral kryos interface.

### 1.2 CISM (glissade)

Verified against `yelmox/docs/_cism-backend-plan.md`; that plan's
analysis holds. Key facts for the interface:

| Item | CISM |
|---|---|
| Role | ice sheet |
| Grid | regular projected x-y (`[grid] ewn, nsn, dew, dns`); thickness at cell centres, velocities at vertices (B-grid); 2-cell halo even in serial for glissade (`libglimmer/parallel_slap.F90:53-58`); regional |
| Clock | **one fixed step per call**: driver loops `time = time + tinc; call glissade_tstep(model,time)` (`cism_driver/cism_front_end.F90:288-332`; `libglissade/glissade.F90:1172`). `model%numerics%time` in years, `real(dp)` |
| Precision | `dp` throughout |
| Params | INI-style `.config` file, read from file only; overridable in memory with `ConfigSetValue` before `glide_config` |
| Constants | `libglimmer/glimmer_physcon.F90`: `scyr = 31536000` (365 d, `:90`), `rhoi = 910`, `rhoo = 1028` (`:67-68`) — module variables, not per instance |
| Restart | own `[CF restart]` output section; `restart = 0/1/2` in config |
| MPI | serial (`parallel_slap`) or MPI (`parallel_mpi`); host may pass a communicator with `parallel_set_info` |

**Entry points**: `glissade_initialise(model, evolve_ice)` (`glissade.F90:89`),
`glissade_diagnostic_variable_solve(model)` (`:2109`),
`glissade_tstep(model, time)` (`:1172`), `glide_finalise(model)`
(`libglide/glide_stop.F90:70`). State is `type(glide_global_type)`
(`libglide/glide_types.F90`), fields are `pointer` arrays set directly.

**Inputs** (`glide_types.F90`):

| Field | Location | Units / convention |
|---|---|---|
| `acab` | `model%climate%acab` (`:1526`) | **m ice/s** in code (m/yr ice in I/O) |
| `smb` | `model%climate%smb` (`:1531`) | mm/yr w.e. (alternative input path) |
| `artm` | `model%climate%artm` (`:1534`) | **degC** |
| `eus` | `model%climate%eus` (`:1586`) | **scalar** sea level |
| `topg` | `model%geometry%topg` (`:1210`) | m |
| `bheatflx` | `model%temper%bheatflx` (`:1788`) | **W m-2, positive down** (geothermal is negative, `:1778-1784`) |
| `bmlt_float_external` | `model%basal_melt` (`:1907`) | m/s, **> 0 is melt** (`:1898`); `whichbmlt_float = 4` |
| `thermal_forcing` | `model%ocean_data` | K, depth levels; `whichbmlt_float = 6` |

**Outputs**: `thck`, `usrf`, `lsrf`, `topg`, `f_ground` (`:1201-1222`),
`calving_flux` kg m-2 s-1 (`:1341`), `bmlt_applied`, velocities.

**Pain points**
- **No 2D sea level.** `eus` is a scalar. FastIsostasy/VILMA2 produce a
  spatial sea surface. A CISM adapter can only pass `eus = bsl` and fold
  the sea-surface anomaly into an effective bed (`topg' = z_bed - dz_ss`),
  which is exact for flotation but shifts grounded `usrf` by `dz_ss` (m).
  Needs an upstream change for exactness.
- After an external `topg`/`eus` change, the halo update and
  `marine_connection_mask` must be redone (`glissade.F90:896`); they only
  run inside the isostasy branch.
- Fatal errors call `parallel_stop` (MPI abort) plus ~20 bare `stop`s.
- Module globals: communicator, `nhalo`, model registry (max 8 models,
  `libglide/glide_model_registry.F90:61`), physical constants. Several
  instances work (CESM runs GrIS and AIS) but share constants and comm.
- Build: CMake, one static library `libglimmercismfortran.a`; Python
  needed at build time (generated I/O); netCDF required.
- GLAD (`libglad/`) is CESM's coupling API: `glad_initialize`,
  `glad_initialize_instance`, `glad_get_grid_size`, `glad_gcm`, ... It
  only takes SMB, surface temperature and ocean T/S, uses integer hours
  (overflow after ~245 kyr) and a hard-coded 1-yr mass-balance step. It is
  a pattern to copy, not an API to use (decision 6 in the CISM plan).

### 1.3 PISM (from memory; not verified locally)

Treat this section as informed recollection of PISM 2.x; details should be
checked against the PISM source before design decisions rest on them.

| Item | PISM |
|---|---|
| Role | ice sheet (C++, PETSc) |
| Grid | regular projected x-y (`pism::Grid`, formerly `IceGrid`), PETSc DMDA-distributed, cell-centred; regional |
| Clock | owns `pism::Time` (seconds since a reference date, CF calendar: `365_day`, `gregorian`, ...). `IceModel::run_to(t)` advances to a target with adaptive steps; `IceModel::step()` does one step |
| Precision | `double` |
| Params | `pism::Config` (netCDF config file + command-line style options); one config per `Context` |
| Restart | `-i` input file / `-o` output; any output file is a restart |
| MPI | required (PETSc); runs with one rank |

**Coupling surface**: physics boundary conditions come from component
class hierarchies, chosen by option strings and stacked with modifiers:
- `surface::SurfaceModel` → `ice_surface_mass_flux` (kg m-2 s-1, + gain)
  and `ice_surface_temp` (K), plus liquid water fraction and layer mass.
  Usually on top of `atmosphere::AtmosphereModel` (air temperature,
  precipitation) through PDD or elevation schemes.
  `-surface given` reads from a file; `-surface given,anomaly,lapse_rate`
  stacks modifiers.
- `ocean::OceanModel` → `shelf_base_mass_flux` (kg m-2 s-1; positive =
  mass loss, from memory), `shelf_base_temperature`, melange back pressure.
  Options include `given`, `constant`, `th`, `pico`.
- `ocean::sea_level::SeaLevel` → 2D `sea_level_elevation` (`-sea_level
  constant|given|...`). Unlike CISM, PISM takes a spatial sea level.
- `frontalmelt::FrontalMelt`.
- `bed::BedDef` → `topg` evolution: `none`, `iso`, `lc` (Lingle-Clark),
  `given` (bed change from file).
- `energy::BedThermalUnit` takes `bheatflx` (W m-2, upward).
- Hydrology (`null`, `routing`, `distributed`) is internal.

**Embedding from Fortran**: PISM has no Fortran API. The likely route is a
C shim: a C++ library that creates `pism::Context` on a communicator
(`MPI_COMM_SELF` for single-process), builds a `Grid` from the kryos grid
definition, subclasses `IceModel` to override `allocate_couplers()` /
`allocate_bed_deformation()` so that the surface, ocean, sea-level and
bed models are "external" classes reading from buffers the shim fills,
and exposes `extern "C"` functions (`pism_init`, `pism_set_field(name,
ptr)`, `pism_run_to(t)`, `pism_get_field`, `pism_write_restart`,
`pism_end`). Fortran calls them through `bind(C)`. PISM's Python bindings
(SWIG) show that this access is possible from outside C++.
Uncertainties: whether the needed `allocate_*` hooks are virtual and
public in the current version; PETSc must be initialised once per process
(`PetscInitialize`), with a global options database; distributed arrays
must be gathered to the shim's plain arrays. OpenMP inside yelmox and
MPI inside PISM must not conflict.

## 2. Isostasy and sea level

### 2.1 FastIsostasy

| Item | FastIsostasy |
|---|---|
| Role | regional GIA (LLRA, ELRA, ELVA, LV-ELVA) + regional sea-surface perturbation + barystatic sea level bookkeeping |
| Grid | regular x-y, **square cells only** (`src/isos_utils.f90:289`), padded for FFT (`init_domain_size`, `:277`); takes only `nx, ny, dx, dy` — no `grid_class` |
| Region outside domain | padding is unloaded land: `ref%z_bed = 1e6`, `Hice = 0` before the domain is copied in (`src/fastisostasy.f90:586-590`); load anomalies outside the domain are zero |
| Clock | **target time**: `isos_update(isos, H_ice, time, bsl, dwdt_corr)` integrates the ODE `do while (ode%t < tf)` with Euler/RK4/BS32/Tsit5 (`src/integrators.f90`), adaptive dt for BS32/Tsit5. Diagnostics (elastic, geoid) every `dt_diagnostics` |
| Time units | years, `real(dp)` |
| Precision | `wp = sp` (`src/precision.f90`) |
| Params | `&isostasy` and `&barysealevel` groups (`par/test_isostasy_test1a.nml`); group is an argument |
| Constants | optional `cnst` at `isos_init` (`fastisostasy.f90:116-123`) |
| Restart | `isos_restart_write/read` (`src/isostasy_io.f90:56,170`); `bsl_restart_write/read` separate (`src/barysealevel.f90:430,457`) |
| Dependencies | FFTW3 (C library, via `fftw3.f03`), fesm-utils ncio/nml |

**Entry points** (`src/fastisostasy.f90`):

| Procedure | Signature |
|---|---|
| `isos_init` (`:67`) | `(isos, filename, group, nx_ice, ny_ice, dx, dy, cnst, K, domain, grid_name)` |
| `isos_init_ref` (`:575`) | `(isos, z_bed, H_ice, bsl, dz_ss)` — reference (relaxed) state |
| `isos_init_state` (`:621`) | `(isos, z_bed, H_ice, time, bsl)` — `bsl` is a `bsl_class` (cold or restart) |
| `isos_update` (`:773`) | `(isos, H_ice, time, bsl, dwdt_corr)` |
| `isos_end` (`:982`) | `(isos)` |
| `bsl_init / bsl_update(bsl, year_bp)` | `barysealevel.f90:42,197`; methods `const`, `file`, `fastiso`, `mixed` |

**Inputs**: `H_ice` [m] on the isostasy grid (total thickness; flotation is
applied internally via `Haf`); optional `dwdt_corr` [m/yr]; the shared
`bsl_class` (in/out).

**Outputs** (`isos%out`, `src/isostasy_defs.f90:195`; cropped to the
domain): `z_bed` = `ref%z_bed + w + we` [m], `w` (viscous), `we`
(elastic), `dwdt` [m/yr], `z_ss` (sea surface) and `dz_ss` [m],
`maskocean/maskgrounded/maskcontinent` (logical), `canom_full` [kg m-2].

**Global coupling through `bsl`**: each domain's `isos_update` subtracts
its volume change from the shared `bsl%bsl_now` (`fastisostasy.f90:835-850`,
`deltaV_bsl / A_ocean`). Several domains thus aggregate into one scalar by
sequential mutation. Order matters, and restart redo needs special
handling (`:812-834`). This is the only global quantity today.

**Pain points**: square cells and FFT padding; `sp`; the domain grid is
anonymous (`nx, ny, dx`) so the framework must map it; the restart-redo
logic around `time_prognostics` is subtle; FFTW plan creation is not
thread-safe (instances are otherwise independent derived types, no `save`).

### 2.2 VILMA2

Read with `yelmox/docs/_vilma2-backend-plan.md`, whose analysis holds.

| Item | VILMA2 |
|---|---|
| Role | global 3D visco-elastic solid Earth + self-consistent sea-level equation + rotation |
| Grid | **global** Gauss-Legendre grid owned by the model (`sht_grid`, built from `lmax`, default 64; `vilma.nml:61-62`). Host may supply any coords lon-lat grid; VILMA2 builds a conservative host→Gauss and bilinear Gauss→host map (`src/vilma_coupling.f90:143-232`, `vilma_remap.f90`) |
| Clock | **advance by dt**: `solid_earth_update(se, h_ice, dt_yr)` (`:354`); load ramped linearly over the interval; internal adaptive stepper in seconds. Owns `se%time`, starting at **0 at init** (relative), `real(wp)` |
| Precision | `wp = real64` enforced equal to `c_double` (`src/vilma_precision.f90`) |
| Params | `&vilma` group over a complete defaults file: `vilma_par_load(p, filename, defaults_file, group)` (`src/vilma_params.f90:212`); `&ctl` only for the standalone driver |
| Constants | **hard-coded** `rho_ice = 931`, `sec_per_year = 31556926` (`src/vilma_constants.f90:27,34`) — not overridable from a shared `phys_const_class` |
| Restart | `vilma_restart_write(self, time, filename, folder, init)` / `vilma_restart_read(self, filename, time)` (`src/vilma_io.f90:66,414`): Maxwell memory, rotation, state |
| Dependencies | fesm-utils coords, FFTW, SHTns (C), netCDF; OpenMP; no library target yet (gap 3 of the plan) |

**Entry points**: `vilma_par_load`, `solid_earth_init(se, z_bed_eq,
h_ice_eq, grid, h_ice_init)`, `solid_earth_spinup(se, h_ice_lgm)` (`:456`),
`solid_earth_update`, `vilma_restart_write/read`, `solid_earth_finalize`
(`:600`).

**Inputs**: global reference bed `z_bed_eq` and reference ice `h_ice_eq`
at init; at each update the **global** ice thickness `h_ice` [m] on the
host grid (flotation applied internally by `ocean_function`).

**Outputs** (host grid and Gauss grid `se%gg`): `rsl` [m] (= N − u + Δφ),
`z_bed = z_bed_eq − rsl` (`:442-453`), ocean function `C`, scalar `bsl`.
Gaps (plan): bed uplift `u` and sea-surface change `N` are not exported
separately; `se%bsl` is an ice-volume-equivalent diagnostic, not the
barystatic rise (`update_bsl`, `:573-597`).

**Pain points**: needs a global load, so ice outside the modelled domains
must come from somewhere (reference or ice history); its time is relative
and its constants are its own; `save` scratch in
`vilma_radial_fe.f90:591` and a module-level I/O table
(`vilma_io.f90:47-48`) — fine for one Earth per run, which is the
intended use anyway.

## 3. Surface mass balance / snowpack

### 3.1 chion

| Item | chion |
|---|---|
| Role | snowpack / SMB: BESSI (energy balance), PDD, ITM |
| Grid | **grid-agnostic columns** (`ncol`); optional `chion_set_grid(chn, x, y, js, is, mask)` for output; active-column mask (`chion_set_active_mask`) |
| Clock | **advance by one step of `dt_days`**; host owns the time loop (`src/chion_api.f90:28-33,349`). Day of year and solar longitude are forcing fields |
| Precision | `wp = sp` by default, `precision=dp` build option; accumulators `wp_acc = dp` (`src/chion_defs.F90:63-67`) |
| Params | groups `chion`, `bessi`, `pdd`, `itm`, `chion_const`, validated against `input/chion_defaults.nml`; group names are arguments, designed for several instances |
| Constants | optional `cnst` at `chion_init` (`chion_const_from_phys`) |
| Restart | `chion_restart_write/read(chn, file, time)` (`src/chion_io.f90`) |
| OpenMP | yes, over columns |

**Entry points**: `chion_init(chn, filename, ncol, group, cnst)` (`:208`),
`chion_init_state` (`:314`), `chion_update(chn, dt_days)` (`:349`),
`chion_get_smb` (`:505`), `chion_get_surface(chn, t_srf, albedo)`
(`:685`), `chion_get_surface_fluxes(chn, melt, runoff, refrz, subl)`
(`:752`), `chion_get_surface_flux_totals` (`:823`),
`chion_tsrf_annual` (`:883`), `chion_end` (`:375`).

**Inputs** (host writes `chn%forc%<field>(:)`, `chion_defs.F90:497-545`):
`air_temperature` [K], `snowfall_rate`, `rainfall_rate` [kg m-2 s-1],
`shortwave_down` [W m-2], `wind_speed` [m s-1]; optional with `has_*`
flags: net SW, downward LW, sensible/latent heat, relative humidity,
air pressure, `surface_height` [m], prescribed albedo, `coszm`, cloud,
dust deposition, `z_sur_std`, host ice albedo, `toa_shortwave`;
`latitude_deg`, `H_ice` [m] (required for BESSI), annual `PDDs`,
`day_of_year`, `solar_longitude_deg`.

**Outputs**: `smb` [kg m-2 s-1, **+ = ice gains**, mean over the last step];
`t_srf` [K], `albedo`; `melt, runoff, refrz, subl` [kg m-2 s-1]; annual
`tsrf` for PDD/ITM via `chion_tsrf_annual`. Missing values are `MV`.

**Pain points**: needs sub-annual (daily or sub-daily) forcing, so a
monthly climate must be disaggregated by someone (smbpal does it
internally; chion does not); SMB is ice-facing and needs firn spin-up;
module-cached variable tables are read-only and shared (fine).

### 3.2 smbpal (yelmox `libs/smbpal/`)

| Item | smbpal |
|---|---|
| Role | SMB from monthly climate: ITM (insolation-temperature melt) or PDD |
| Grid | 2D, `smbpal_init(smb, filename, x, y, lats, group, itm_group, orbit_group)` (`smbpal.f90:105`); regional |
| Clock | **one year per call**: `smbpal_update_monthly(smb, t2m, pr, z_srf, H_ice, time_bp)` (`:209`); internal 360-day year (`:12`), ITM every 2 days, PDD every 10 days. `smbpal_update_monthly_equil` for spin-up (`:175`). `time_bp` only sets the orbit |
| Precision | `prec = sp` (`smbpal_precision.f90`); uses insol (default `dp`) through `insol_wp` |
| Params | groups `smbpal`, `itm`, `orbit` |
| Restart | `smbpal_restart_write/read` (`:632,694`): annual fields + snowpack |

**Inputs**: monthly `t2m` [K], `pr` [mm w.e./d] (nx, ny, 12); `z_srf`,
`H_ice` [m]; `time_bp` [yr rel. 1950] for insolation.
**Outputs**: `smb%ann%smb` [mm w.e./yr], `tsrf` [K], `melt`, `runoff`,
`refrz`, `H_snow`, `alb_s`, `PDDs`.

### 3.3 smb_simple (yelmox `libs/smb_simple.f90`)

| Item | smb_simple |
|---|---|
| Role | parametric SMB (synthetic snow-line profile; TG24 scheme not yet wired) |
| Grid | 2D `x, y, lat` stored at init (`:190`); regional |
| Clock | stateless: `smb_simple_update(smbs, z_srf, t_sl)` (`:284`) |
| Precision | `sp` |
| Inputs | `z_srf` [m], sea-level air temperature `t_sl` [K]; constant CO2 and insolation parameters; target ice mask from file or `H_ice_ref` (`smb_simple_set_mask`, `:246`) |
| Outputs | `smb` [mm w.e./yr], `t_srf` [K] (lapse-rate) |
| Restart | none needed |

## 4. Climate / atmosphere

### 4.1 rembo1

| Item | rembo1 |
|---|---|
| Role | regional energy-moisture balance atmosphere + PDD/ITM SMB |
| Grid | **fixed at compile time**: `nxs = 106, nys = 181, dxs = 16 km` (`src/emb_global.f90:60-61`, Greenland on the Yelmo GRL-16KM grid); internal (j,i) ordering, transposed on output (`climate.f90:593-601`) |
| Clock | **target time with internal cadence**: `rembo_update(time, time_ins, dT_summer, z_srf, H_ice, z_sl, mask_relax, dT_mon, co2)` (`src/climate.f90:200`) decides itself whether the EMB (`dtime_emb`) and SMB (`dtime_smb`) are due (`:314-326`). 360-day year |
| Precision | `real(8)` |
| Params | `options_rembo` fixed-column text file (`options_rembo`) plus `rembo_Greenland.nml` |
| Restart | `rembo_restart_write(filename, time, z_srf, H_ice, z_sl)` (`:833`); read via options |
| State | **module globals**: `rembo_ann` (`climate.f90:38`), `ppcorr0`, ~175 variables in `emb_global` |

**Inputs**: `time`, `time_ins` (insolation time, yr BP), `dT_summer` [K]
or monthly `dT_mon(12)`, `co2` [ppm], `z_srf`, `H_ice`, `z_sl` [m],
`mask_relax`. Reads its own reference climatology and topography files.
**Outputs** (`rembo_ann`): `T_ann`, `T_jja`, `T_djf`, `T_srf` [K], `pr`,
`sf`, `smb`, `melt`, `runoff`, `refrz` (mm w.e./yr as consumed by yelmox
`climate_out_class`), `H_snow`, `pdds`.

**Pain points**: one instance per executable; grid must match the ice grid
exactly at compile time; anomaly forcing is driver-specific (summer dT);
legacy text options. Embedding requires either the existing stub pattern
(yelmox `climate_rembo.f90` / `climate_rembo_stub.f90`) or a refactor to
a derived type with a runtime grid.

### 4.2 insol

| Item | insol |
|---|---|
| Role | utility: daily-mean TOA insolation from Laskar (2004) orbits; orbital parameters |
| Grid | none, or lat array given at `insol_init(is, lats, ...)` (1D/2D) |
| Clock | stateless queries: `calc_insol_day(is, insol, day, time_bp)`, `calc_insol_days`, `calc_insol_ave`, `calc_insol_diurnal` (`src/insolation.f90:143-155`); `day_year` 360/365/366, default 360 |
| Precision | `wp = dp` default, `sp` with `-DINSOL_SP` (`src/insol_defs.F90:25-29`) |
| Notes | no dependencies; state read-only after init and thread-safe |

insol is not a component. It is a shared service: an "orbit/insolation
forcing" that several components (smbpal, chion, rembo, climate) need.
kryos should provide orbital elements as a run-level forcing and let each
backend call insol itself.

### 4.3 Existing neutral climate output

yelmox already has a backend-neutral climate record, `climate_out_class`
(`yelmox/libs/climate_out.f90`): monthly `tas` [K], `pr` [mm/d], annual
`ta_ann`, `tsl_ann`, `pr_ann`, the `z_srf` the fields are valid at, a
monthly lapse rate and `beta_p`; ocean as profiles (`to_ann`, `so_ann`,
`depth`) or at the shelf base (`T_shlf`, `S_shlf`, anomalies); optional
`smb`/`tsrf` and `Qd`, with `has_*` capability flags. It is a good
starting point for the kryos climate export.

## 5. Ocean and marine shelf

### 5.1 marine_shelf + pico (yelmox `libs/marine_shelf.f90`, `libs/pico/`)

| Item | marine_shelf |
|---|---|
| Role | sub-shelf melt: linear/quadratic/local/non-local thermal-forcing laws, PICO, ISMIP6 corrections |
| Grid | 2D regional (`marshelf_init(mshlf, filename, group, nx, ny, domain, grid_name, reg, xc, yc, cnst)`, `:628`) |
| Clock | stateless per call: `marshelf_update_shelf(...)` (`:161`) maps ocean profiles to the shelf base; `marshelf_update(mshlf, H_ice, z_bed, f_grnd, z_sl, dx, z_srf)` (`:393`) computes melt |
| Precision | `wp = sp` (`:20`) |
| Constants | `cnst` at init |
| Restart | `marshelf_restart_write/read` (`:1805,1922`) |

**Inputs**: geometry `H_ice, z_bed, f_grnd, z_sl, z_srf` [m], `dx` [m];
ocean `depth(:)`, `to_ann, so_ann, dto_ann` (nx, ny, nz) [K, psu], optional
`tf_ann`; basins/regions at init. PICO: `to, so, H_ice, f_grnd, basins,
is_ocn, dx` (`pico.f90:72`).
**Outputs**: `bmb_shlf` [m ice/yr, **melt < 0**, `:710`], `T_shlf`,
`T_fp_shlf` [K], `S_shlf` [psu], `tf_shlf`, `tf_basin`, `mask_ocn`
(0 land, 1 grounding line, 2 floating line, 3 ocean, 4 deep ocean,
5 lake; `:24-29`).

marine_shelf is really two things: an ocean→shelf-base interpolator and a
melt law. The melt law depends on ice geometry and could equally live in
the ice-sheet backend (CISM and PISM have their own).

### 5.2 obm (yelmox `libs/obm/`)

| Item | obm (nautilus) |
|---|---|
| Role | 4-box global ocean (N, tropics, tropical deep, S) for bipolar runs |
| Grid | none (scalars); **global** |
| Clock | **advance by dt**: `obm_update(obm, dt, name)` (`obm.f90:49`) |
| Precision | `sp` |
| Inputs | external freshwater fluxes `fn, ft, fs` [Sv]; `calc_fwf` (`ice2ocean.f90`) builds them from ice-sheet `mb_net/smb/bmb/cmb/dHidt` [m ice/yr] over a hydrographic mask |
| Outputs | box temperatures [degC], converted to K by `obm_temperature_kelvin` (`ocean2ice.f90`), salinities, overturning `m` |
| Restart | `write_obm_restart`, `obm_init(..., restart)` |

obm is the first example of a global component fed by several regional
ice sheets and feeding back a scalar or low-dimensional field.

## 6. Other components

| Backend | Role | Interface | Notes |
|---|---|---|---|
| geothermal (`libs/geothermal.f90`) | static `ghf` field | `geothermal_init(gthrm, filename, grid, domain, group)` (`:48`) | mW m-2; data from file or constant; `sp`; no update, no restart |
| sediments (`libs/sediments.f90`) | static `H` sediment thickness [m] | `sediments_init(sed, filename, grid, domain, group)` (`:42`) | remap on read; `sp`; no update, no restart |
| FastHydrology (`yelmo/FastHydrology/src/fast_hydrology.f90`) | basal hydrology (till bucket, K24 transport, N closures) | `hydro_init(hyd, filename, nx, ny, dx, dy, group, cnst, sec_year, ...)` (`:142`), `hydro_init_state` (`:197`), `hydro_update(hyd, H_ice, z_bed, z_sl, f_ice, f_grnd, mask, mdot, G, q_T, i_eb, uxy_b, A_glen, time, ...)` (`:300`) | target time in years (`dp`), SI rates (`:33-41`); `sp`; needs basal melt, velocities, `A_glen` — ice-model internals |

geothermal and sediments are already "data backends" in all but name.
FastHydrology is tightly coupled to ice dynamics (basal velocity, rate
factor, heat fluxes) and is called inside `yelmo_update`
(`yelmo_ice.f90:428`); CISM and PISM have their own hydrology.

## 7. fesm-utils facilities relevant to adapters

- `coords` `grid_class`: projected and lon-lat grids, including Gauss
  (`grid_cdo`, `gaussian_latitudes`). Yelmo, geothermal, sediments,
  VILMA2 take it; FastIsostasy, chion, smbpal, marine_shelf, rembo do not.
- `coupler` (`src/coupler.f90`): named grid registry (max 16) + cached maps
  (max 64), `remap(cpl, var_src, src, var_dst, dst, method, stat)` for
  `sp`/`dp`/int 2D/3D; methods `con`, `nn`, `shepard`, `quadrant`, `bilin`.
- `htopo` (`src/htopo.f90`): the high-resolution geometry hub —
  `z_bed, H_ice, z_srf, f_grnd, z_sl` refreshed by
  `htopo_update(htopo, dz_bed, dH_ice, z_sl)` (`:136`). It is already the
  minimal shared ice geometry.
- `timestepping` `tstep_class`: `dp` driver clock with `time`,
  `time_elapsed`, `time_cal`, `time_rel`.
- `phys_constants` `phys_const_class`: shared constants record, accepted by
  Yelmo, FastIsostasy, chion, marine_shelf, smb_simple, FastHydrology. Not
  by VILMA2, CISM, rembo1.
- `precision`: `wp = sp`.

## A. Minimum common interface per component

Conventions used below: fields are 2D on the component's own grid
(cell centres) unless noted; "rate" fields are means over the last
coupling interval. The proposed canonical units are SI-like but follow
the community's ice-sheet habit for rates per year; whichever is chosen,
the adapter owns the conversion (section D).

### A.1 Ice sheet (Yelmo, CISM, PISM, data)

| Direction | Field | Unit | Satisfied by |
|---|---|---|---|
| import | `z_bed` | m | all (CISM `topg`, PISM `topg` via bed model) |
| import | `z_sl` (2D) | m | Yelmo, PISM; CISM only scalar `eus` + bed folding |
| import | `smb` | kg m-2 s-1 or m ice/yr, + gain | all |
| import | `T_srf` | K | all (CISM in degC) |
| import | `bmb_shlf` | same as smb, + gain (melt < 0) | Yelmo, CISM (`bmlt_float_external`, sign flip), PISM (`shelf_base_mass_flux`, sign flip) |
| import (opt) | `T_shlf`, `tf_shlf` / ocean T,S profiles | K, psu | Yelmo; CISM `thermal_forcing`; PISM ocean model |
| import (opt) | `Q_geo` | W m-2, + into ice | all (Yelmo mW m-2; CISM negative-down) |
| import (opt) | `H_sed`, `Qd`, `fmb_shlf` | m, -, m ice/yr | Yelmo; not CISM |
| export | `H_ice`, `z_srf` | m | all |
| export | `z_bed` actually used | m | all |
| export | `f_grnd`, `f_ice` | 0-1 | Yelmo; CISM `f_ground`; PISM cell type → 0/1 |
| export | `z_base` | m | derivable (`z_srf - H_ice`) |
| export | `dHidt`, `smb_applied`, `bmb_applied`, `calv` (`cmb`), `mb_net` | rate | Yelmo; CISM (`calving_flux`, `bmlt_applied`, `acab_applied`); PISM (diagnostic fluxes) |
| export (opt) | `ux_s, uy_s` (cell centres) | m/yr | all, after averaging staggered velocities |
| export (opt) | `Q_b` / subglacial discharge | | Yelmo (FastHydrology), PISM; CISM partial |

Hidden behind the backend: 3D temperature, stresses, rate factor, basal
friction and its inversion (`cb_ref`), hydrology state, tracers, numerics,
regions output, spin-up helpers (`yelmo_update_equil`). Static
configuration (mask of allowed ice, basins, reference topography for
relaxation) is passed at init, not through the coupling fields.

Two choices to settle: (1) masks — export fractions (`f_grnd`, `f_ice`)
and let kryos derive one canonical integer mask, rather than mapping
Yelmo's 8-valued `mask_bed`, CISM's `thkmask` and PISM's cell type to each
other; (2) bed — pass the full `z_bed` (simplest, works for CISM/PISM) or
`dz_bed` relative to the backend's own reference (what Yelmo does now, to
keep its high-resolution bed). Passing the full bed with the reference
held in the framework's geometry hub (`htopo`) serves both.

A **data backend** for the ice sheet (prescribed `H_ice`, `z_srf`,
`f_grnd` from a file, e.g. ICE-6G) is needed for isostasy-only or
climate-only runs and fits this interface directly.

### A.2 Isostasy and sea level (FastIsostasy, VILMA2, data)

| Direction | Field | Unit | Notes |
|---|---|---|---|
| init | reference `z_bed_ref`, `H_ice_ref` | m | both backends: relaxed state |
| import | `H_ice` (total) | m | both apply flotation internally; one needs `z_sl` to do so consistently |
| import (opt) | `dzbdt_corr` | m/yr | FastIsostasy only |
| export | `dz_bed` (bed displacement vs. reference) | m | FastIsostasy `w + we`; VILMA2 needs `u` export (plan gap 1) |
| export | `dz_ss` (sea-surface change vs. reference) | m | FastIsostasy `z_ss`; VILMA2 needs `N` export |
| export | `bsl` (scalar barystatic) | m | FastIsostasy `bsl_class`; VILMA2 needs correct definition (gap 2) |
| export (opt) | `dz_bed_dt` | m/yr | FastIsostasy `dwdt` |

`z_sl = bsl + dz_ss` and `z_bed = z_bed_ref + dz_bed` are formed by the
framework on the ice grid. Hidden: viscosity/lithosphere fields, Maxwell
memory, SHT grid, ODE stepper, rotation. A data backend reads a
`z_bed`/`z_sl` (or `dz_bed`/`dz_ss`) history, or a `bsl` series only
(FastIsostasy `bsl` method `file` already does that).

### A.3 Surface / SMB (chion, smbpal, smb_simple, rembo's SMB, data)

| Direction | Field | Unit | Notes |
|---|---|---|---|
| import | `z_srf`, `H_ice` | m | all |
| import | near-surface air temperature | K | chion (sub-daily/daily), smbpal (monthly), smb_simple (sea-level annual) |
| import | precipitation (or snowfall + rainfall) | kg m-2 s-1 | chion split; smbpal total mm/d |
| import (opt) | radiation (SW down, LW down), wind, humidity, pressure, cloud | SI | chion BESSI only |
| import (opt) | orbit / `time_bp` | | smbpal ITM, chion via `toa_shortwave`/`day_of_year` |
| export | `smb` (interval mean, + gain) | kg m-2 s-1 | chion native; smbpal mm w.e./yr |
| export | `T_srf` (ice surface, ≤ 273.15 over ice) | K | all |
| export (opt) | `melt`, `runoff`, `refreeze`, `albedo` | kg m-2 s-1, 1 | chion, smbpal |

Hidden: snow layers, firn, albedo state, PDD statistics, daily
disaggregation. The import side splits into a **required core** (T, P,
geometry) and **capability-flagged extras** (chion's `has_*` pattern).
Whether an SMB backend needs daily forcing (chion) or monthly (smbpal) is
an adapter concern: the climate export stays monthly, and the chion
adapter disaggregates (as smbpal does internally).

### A.4 Climate / atmosphere (rembo1, snapshots/ESM forcing, data)

| Direction | Field | Unit | Notes |
|---|---|---|---|
| import | `z_srf`, `H_ice` (+ `z_sl`) | m | rembo |
| import | run forcing: time, `time_bp`/orbit, CO2, temperature anomaly | | rembo; snapshot/ESM backends |
| export | monthly `tas` | K | `climate_out_class` |
| export | monthly `pr` | kg m-2 s-1 (now mm/d) | |
| export | `z_srf_clim` (elevation the fields refer to), lapse rate | m, K/m | needed to downscale to another surface |
| export (opt) | radiation, wind, humidity | | for chion BESSI |
| export (opt) | `smb`, `T_srf` | | rembo provides an SMB directly (`has_smb`) |
| export (opt) | ocean `to, so` profiles + `depth` | K, psu, m | ESM forcing |
| export (opt) | `Qd` | | |

### A.5 Ocean / shelf (marine_shelf, pico, obm, data)

Split into two components:
- **ocean**: exports ocean `T, S` (profiles on `depth`, or at the shelf
  base) and imports freshwater flux from the ice (obm: Sv per box). Global
  (obm) or regional/data (ISMIP6/ESM fields).
- **shelf melt**: imports geometry (`H_ice, z_bed, z_srf, f_grnd, z_sl`)
  and ocean `T, S`; exports `bmb_shlf` (+ gain), `T_shlf`, `tf_shlf`.
  Backends: marine_shelf laws, PICO, or the ice-sheet backend's own
  (CISM `whichbmlt_float = 6`, PISM `-ocean pico`). The ice-sheet
  interface should allow "bmb computed internally from T/S" as a
  capability, so the shelf-melt component can be skipped.

### A.6 Geothermal, sediments

Export-only, static: `Q_geo` [W m-2, + into ice], `H_sed` [m]. The
current modules are already data backends. No update, no restart.

### A.7 Hydrology

No separate kryos component for now. All three ice models embed their
own hydrology, and FastHydrology needs ice internals (basal melt,
velocity, `A_glen`). The ice-sheet export may carry `Q_b` (basal
water/discharge) for the ocean/frontal-melt side (`Qd`).

## B. Time and clock conventions

| Backend | Native clock | Units | Adapter for `update(comp, t_target)` |
|---|---|---|---|
| Yelmo | target time, own clock | yr, dp | direct |
| CISM | one `tinc` step per call | yr, dp | loop `glissade_tstep` until `t_target`; require `(t_target - t)` a multiple of `tinc` |
| PISM | `run_to(t)`, own clock | s, CF calendar | convert years → seconds with PISM's calendar |
| FastIsostasy | target time | yr, dp | direct |
| VILMA2 | advance by `dt`, clock from 0 | yr, wp=dp | `dt = t_target - t_last`; adapter stores the absolute offset |
| chion | one step of `dt_days`, host owns time | days | loop sub-steps, pull forcing per sub-step (disaggregated), average fluxes |
| smbpal | one year per call | 360-d year | call once per year; `time_bp` from framework |
| smb_simple, marine_shelf | stateless | — | call at `t_target` |
| rembo1 | target time, internal cadence | yr | direct |
| obm | advance by `dt` | yr | `dt = t_target - t_last` |
| FastHydrology | target time | yr, dp | (inside ice backend) |
| insol | stateless query | day, `time_bp` | service |

Recommendation:
- One framework clock (`tstep_class`, `dp` years, absolute model time
  plus `time_bp` for orbital forcing). Every component is driven with
  `update(comp, t_target)` and keeps `t_last`. The adapter translates.
- Each component has its own coupling interval `dt_comp`; it is called
  when `t ≥ t_last + dt_comp`. Fields it exports are held between calls;
  rate fields it exports are means over its last interval.
- `update` must be idempotent for `t_target == t_last` (Yelmo already is;
  FastIsostasy has special logic for it). This makes restart redo trivial.
- Restarts only at times that are multiples of every component's
  interval (as GLAD requires and chion's accumulators imply).
- Year length is a framework constant passed to every backend that takes
  `phys_const_class`. Today it is 365.2422 d (Yelmo, FastHydrology,
  VILMA2), 365 d (CISM `scyr`), 360 d (smbpal, rembo1, insol default).
  The 360-day "year" of smbpal/rembo is a calendar for sub-annual cycles,
  not a conversion factor; keep it internal. CISM's and VILMA2's year
  cannot be overridden without patching.
- Multi-domain global solves (VILMA2, obm, bsl) need a three-phase step:
  regional `export` → global `solve` once → regional `import` (the
  load/solve/land split of the VILMA2 plan).

## C. Regional versus global

| Component | Backend | Scope | Grid | Global aggregation needed | Global → regional |
|---|---|---|---|---|---|
| ice sheet | Yelmo, CISM, PISM | regional | projected x-y | ice load for VILMA2; volume change for bsl; FWF for obm | — |
| isostasy | FastIsostasy | regional + shared scalar `bsl` | x-y padded | ΔV per domain → `bsl` | `bsl` |
| isostasy | VILMA2 | **global** | Gauss (T64) | `H_ice` of all domains spliced into a global load, plus ice outside the domains (reference or history) | `dz_bed`, `dz_ss`, `bsl` remapped to each domain |
| SMB | chion, smbpal, smb_simple | regional (chion: column set) | any | — | — |
| climate | rembo1 | regional (fixed grid) | x-y | — | — |
| climate | ESM/snapshot forcing | global data | lon-lat | — | regridded fields, anomalies |
| ocean | obm | **global** (boxes) | scalars | FWF from all domains (Sv) | box T per domain/hemisphere |
| shelf melt | marine_shelf, pico | regional | x-y | — | — |
| geothermal, sediments | data | regional | x-y | — | — |

Design consequence: a component instance is either per domain or per run.
Global instances (VILMA2, obm, global climate data, `bsl`) need an explicit
aggregation step that sums contributions from all domains on the global
grid (with a footprint fraction and an "outside the domains" fill), and a
distribution step back. The only fields that cross the boundary are ice
load (m ice or kg m-2), freshwater flux, and the returned `dz_bed`,
`dz_ss`, `bsl`, ocean T/S and climate anomalies. Load splicing must be
mass-conservative (`con` remap); returned fields are smooth (`bilin`).

## D. Precision, units, masks, signs

### Precision

| Backend | Fields | Time |
|---|---|---|
| Yelmo, FastHydrology, fesm-utils | `sp` (fixed) | `dp` |
| FastIsostasy | `sp` | `dp` |
| chion | `sp` (opt. `dp`); accumulators `dp` | host |
| smbpal, smb_simple, marine_shelf, obm, geothermal, sediments | `sp` | `dp` / none |
| insol | `dp` (opt. `sp`) | `dp` |
| VILMA2 | `dp` (`real64`, must equal `c_double`) | relative `dp` |
| CISM | `dp` | `dp` |
| rembo1 | `real(8)` | `real(8)` |
| PISM | `double` | `double` (seconds) |

The framework should store exchange fields in one kind (likely `dp`, or a
kryos `wp`) and convert inside each adapter. fesm-utils `coupler` already
remaps both `sp` and `dp`. Time stays `dp` everywhere.

### Units the framework must own

| Quantity | Yelmo | CISM | PISM | chion | smbpal / rembo / smb_simple | marine_shelf |
|---|---|---|---|---|---|---|
| SMB | m ice/yr | m ice/s (`acab`) or mm w.e./yr (`smb`) | kg m-2 s-1 | kg m-2 s-1 | mm w.e./yr | — |
| T_srf | K | degC (`artm`) | K | K | K | — |
| sub-shelf melt | m ice/yr | m ice/s | kg m-2 s-1 | — | — | m ice/yr |
| geothermal | mW m-2 | W m-2 | W m-2 | — | — | — |
| precipitation | — | — | kg m-2 s-1 | kg m-2 s-1 | mm w.e./d | — |
| ocean T | K | K | K (Kelvin, from memory) | — | — | K (obm degC) |
| time | yr | yr (s internal) | s | days | yr | — |

Conversions between ice and water equivalent need `rho_ice`, which
differs: Yelmo 910, CISM 910 (917 commented), VILMA2 931 (hard-coded),
FastHydrology default 917, marine_shelf 910. The framework should own one
`phys_const_class`, pass it to every backend that accepts it, and flag
those that do not (VILMA2, CISM, rembo1) so the mismatch is visible.
Recommendation: exchange mass fluxes as kg m-2 s-1 (no density ambiguity;
chion and PISM native), geometry in m, temperatures in K, heat flux in
W m-2.

### Signs

| Quantity | Yelmo | CISM | PISM (memory) | marine_shelf | chion |
|---|---|---|---|---|---|
| smb | + gain | + gain | + gain | — | + gain |
| basal / shelf melt | **melt < 0** (`bmb`) | **melt > 0** (`bmlt`) | **loss > 0** (`shelf_base_mass_flux`) | melt < 0 | — |
| geothermal | + into ice | **+ down** (geothermal < 0) | + up | — | — |
| FWF to ocean | — | — | — | — | obm: + for ice loss |

Proposal: every exchanged mass flux is "+ = ice gains", named `*mb`.
Adapters flip signs for CISM/PISM melt and CISM `bheatflx`.

### Masks

| Backend | Mask | Values |
|---|---|---|
| Yelmo | `mask_bed` | 0 ocean, 1 land, 2 frozen, 3 stream, 4 grline, 5 float, 6 island, 7 partial |
| Yelmo | `mask_ice` (input) | 0 none, 1 fixed, 2 dynamic |
| FastIsostasy | logical | `maskocean`, `maskgrounded`, `maskcontinent` |
| VILMA2 | ocean function `C` | 1 ocean, 0 land |
| marine_shelf | `mask_ocn` | 0 land, 1 grline, 2 floating line, 3 ocean, 4 deep ocean, 5 lake |
| CISM | `thkmask` bit mask; `f_ground` | |
| PISM (memory) | cell type | 0 ice-free land, 2 grounded, 3 floating, 4 ice-free ocean |

Exchange fractions (`f_ice`, `f_grnd`) and derive a small canonical mask
(ice-free land / grounded ice / floating ice / ice-free ocean) in the
framework. Backend masks stay internal.

## E. Embedding concerns

| Backend | Language | Build | Parallelism | Global/module state | Several instances? | Error handling |
|---|---|---|---|---|---|---|
| Yelmo | Fortran | `libyelmo.a` + deps (FastHydrology, elsa, tracer, Lis) | OpenMP | none of note | yes | `error stop`; kill file |
| CISM | Fortran 90 | CMake static lib, Python at build, netCDF | serial or MPI | comm, `nhalo`, registry (8), constants | yes (≤ 8, shared comm/constants) | `parallel_stop` → MPI abort |
| PISM | C++ | CMake, PETSc, MPI, netCDF, UDUNITS, GSL, FFTW | MPI (PETSc) | PETSc init + options database | probably, per `Context` | C++ exceptions (must be caught in the shim) |
| FastIsostasy | Fortran | `libisostasy.a`, FFTW | OpenMP variant | none | yes | `stop` |
| VILMA2 | Fortran 2008 | executables only (needs `libvilma.a`), FFTW, SHTns | OpenMP | `save` scratch, I/O table | yes in principle, one needed | `error stop` |
| chion | Fortran | `libchion.a` | OpenMP | read-only var tables | yes (by design) | `error stop` |
| smbpal, smb_simple, marine_shelf, pico, obm, geothermal, sediments | Fortran | inside yelmox | serial / OpenMP | none | yes | `stop` |
| rembo1 | Fortran 90 | `librembo.a` | serial | **all state in modules; grid compile-time** | **no** | `stop` |
| insol | Fortran | `libinsol.a`, no deps | thread-safe | read-only tables | yes | — |

Implications:
- **Optional linkage**: each external backend (CISM, PISM, VILMA2, rembo1,
  chion) behind a build flag with a stub module of the same interface,
  as yelmox does for rembo (`climate_rembo_stub.f90`). With an abstract
  type per component, the stub is simply an unregistered extension.
- **Language boundary**: a Fortran `bind(C)` layer is the common
  denominator. Yelmo already has one; PISM needs a C++ shim exposing
  `extern "C"` init/set/run_to/get/restart/end. A string-keyed
  get/set-field interface (Yelmo C API style) keeps the boundary small.
- **Threading**: yelmox is OpenMP-only. CISM serial or one-rank MPI, and
  PISM on `MPI_COMM_SELF`, must coexist with OpenMP regions in other
  components; MPI must be initialised by the framework, once.
- **Instances**: two Yelmo domains work today (yelmox_bipolar). rembo1
  prevents two climate domains and ties climate to one compile-time grid;
  it needs refactoring into a derived type or must stay a one-domain
  special case. CISM's registry and shared constants limit independent
  configurations.
- **Fatal errors**: no backend returns error codes; all stop. The
  framework cannot recover, only write diagnostics before calling them.
- **Constants and calendars** cannot be unified for CISM, VILMA2 and
  rembo1 without patches; the framework should check and report the
  differences at init.
- **Restart**: every stateful backend writes its own file (Yelmo, CISM,
  PISM, FastIsostasy + bsl, VILMA2, chion, smbpal, marine_shelf, obm,
  rembo1). The framework should own a restart bundle directory and give
  each component a file name, as yelmox's `restart-<kyr>-kyr/` does.

## Open questions raised by this study

1. Exchange precision: `dp` exchange fields with adapters converting, or
   a kryos `wp = sp` matching the current stack?
2. Exchange units for mass fluxes: kg m-2 s-1 (density-free) or m ice/yr
   (Yelmo's, needs a single `rho_ice`)?
3. Bed: full `z_bed` from a framework-held reference, or anomalies to
   each backend's own reference?
4. CISM's scalar `eus`: accept the bed-folding approximation, or require
   an upstream 2D sea level before CISM is supported with regional sea
   level?
5. Shelf melt: separate component, ice-sheet capability, or both?
6. rembo1: refactor to an instantiable type with a runtime grid, or keep
   it as a single-domain special case?
