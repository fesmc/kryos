# CLIMBER-X: lessons for kryos

Study of `/Users/alrobi001/models/climber-x` (state as of 2026-10). All paths below are
relative to that repo unless noted. Line numbers are file line numbers.

Scope: how CLIMBER-X couples regional ice-sheet domains (SICOPOLIS or Yelmo), the
SMB/BMB components, and the global `geo` component (topography, sea level, isostasy,
VILMA1/VILMA2) to a global climate model through a central `cmn` hub, and which patterns
kryos should copy or avoid.

---

## 0. Summary for the impatient

| Topic | What CLIMBER-X does | Verdict for kryos |
|---|---|---|
| Driver | One `program climber` with a hand-written fixed-order loop, every component guarded by `flag_X .and. time_call_X` | Keep the explicit, readable sequence; replace hand-written guards with a component registry + cadence objects |
| Hub | `cmn_class`: ~200 allocatable fields on the 5°×5° coupler grid, plus pairwise `X_to_Y` routines in a 4600-line `coupler.f90` | Avoid. Use per-component typed import/export states; the coupler only routes and remaps |
| Ice interface | `ice_class` (≈20 2D fields on the native ice grid) + `select case(model)` in every procedure of `ice_model.f90` | Keep the idea of a neutral exchange type with canonical units; drop the string dispatch, module-global backend arrays and backend output code in the interface module |
| Multi-domain | `ice(n)`, `smb(n)`, `bmb(n)` arrays indexed by domain; grids from a named-namelist registry (`ice_grids.nml`) | Keep: domain = named grid + bundle of regional components. Make the bundle an object |
| Regional→global | Overlay with `map_field(..., reset=.false.)` on top of a prescribed global background; coarse-cell override masks (`mask_ice`, `mask_smb`) decide regional vs global fallback | Keep the "background + regional overlay" idea; make overlap/priority and conservation explicit |
| Global geo/sea level | `geo` owns a global hi-res lat-lon topography; ice sees bed relative to current local sea surface, `z_sl=0` | Keep the separation (global iso/SL component owns bed+sea level); expose `z_bed` and `z_sl` explicitly in a fixed datum |
| VILMA2 coupling | `vilma_coupling`: `init(z_bed_eq,h_ice_eq,grid)`, `update(h_ice,dt)` → `rsl`, `z_bed`; model owns its grid + remap | Best interface in the codebase. Template for kryos components |
| Dummy/data | Three unrelated mechanisms: link-time stub modules (`.ice_model_dummy.f90`), `bnd/fake_*` data readers filling the same `cmn` slots, CPP `#ifdef VILMA1` stubs | Keep the concept (data backend writes the same export slots); unify into one runtime backend mechanism |

---

## 1. Overall architecture

### 1.1 Three grid tiers

| Tier | Grid | Defined in | Who lives there |
|---|---|---|---|
| Coupler / climate | 5°×5° lat-lon, `ni=72, nj=36`, hard-coded parameters | `src/main/climber_grid.f90:39-46` | `cmn`, atm, ocn, sic, lnd, bgc, `geo%grid` |
| Global hi-res geography | Regular lat-lon read from `geo_ref_file` (RTopo-2.0.1 at 0.125°), must be an integer subdivision of the 5° grid | `src/geo/geo_grid.f90:54-168` | `geo%hires` (bed, ice thickness, masks, lakes, runoff routing, q_geo, sediments) |
| Regional ice domains | Polar-stereographic / stereographic km grids, looked up by name in `nml/ice_grids.nml` | `src/main/climber_grid.f90` `ice_grid_init` (≈l.97-136) | `ice(n)`, `smb(n)`, `bmb(n)` — all on the same grid |

A VILMA Gauss-Legendre grid is a fourth tier, hidden inside `geo` (VILMA1) or inside the
VILMA2 library.

### 1.2 Driver and component inventory

`src/main/climber.f90` is the single driver program. Component state is held as local
variables of the program (`climber.f90:127-142`):

```fortran
type(cmn_class) :: cmn
type(atm_class) :: atm ; type(ocn_class) :: ocn ; ... ; type(geo_class) :: geo ; type(bnd_class) :: bnd
type(ice_class), allocatable :: ice(:)
type(smb_in_class) :: smb_in
type(smb_class), allocatable :: smb(:)
type(bmb_class), allocatable :: bmb(:)
type(grid_class), allocatable :: ice_grid(:)
```

Every component follows the same informal shape: `X_init`, `X_update`, `X_end`,
`X_write_restart`, `X_diag_init`, `X_diag`, plus a `X_def.f90` type module and an
`X_params.f90` namelist module. There is no abstract type or formal interface; the
uniformity is by convention only.

### 1.3 The `cmn` hub

`cmn_class` (`src/main/coupler.f90:128-367`) is a flat type with ~200 allocatable fields
on the coupler grid (fractions, elevations, radiation, fluxes, ocean T/S for BMB, lake
state, runoff routing indices, albedos from SEMI, coral hypsometry, ...).
All exchange routines live in `coupler.f90` (4613 lines) as pairwise copies:

- star pattern through `cmn`: `cmn_to_atm`/`atm_to_cmn`, `cmn_to_ocn`, `cmn_to_lnd`, `cmn_to_smb` (2209), `smb_to_cmn` (2429), `cmn_to_bmb` (2650), `ice_to_cmn` (2833), `cmn_to_geo` (3188), `geo_to_cmn` (3220), `bnd_to_cmn` (3362);
- direct component-to-component bypasses of `cmn`: `ice_to_smb` (2530), `smb_to_ice` (2621), `geo_to_smb` (2567), `ice_to_bmb` (2719), `bmb_to_ice` (2807), `geo_to_bmb` (2750), `ice_to_geo` (2906), `geo_to_ice` (3023), `bnd_to_geo` (3096).

So `cmn` is the hub only for the global climate grid; everything regional (ice/SMB/BMB)
and the global hi-res geo talk to each other directly through coupler routines that
know both types. The coupler also contains physics: runoff merging and routing
(`runoff_merge` 3765, `runoff_to_ocn` 3842), lake water balance (`lakes_update` 3620),
rain/snow partitioning of prescribed precipitation (in `bnd_to_cmn`, ≈3585-3600), sea-level
change from ice volume (`ice_to_geo`, 2973-3010), time relaxation of fluxes.

`coupler.f90` `use`s the `_def` module of every component (`coupler.f90:56-70`), so it is
the dependency sink of the whole model.

### 1.4 Coupling sequence per step

Initialisation order (`climber.f90:166-363`): control → timer → restart dirs →
`climber_grid_init` → `cmn_init` → `bnd_init` → `geo_init` (from bnd data) → ice grids →
`ice_init` per domain (from geo hi-res fields) → `ice_to_geo` → `geo_update` →
`geo_to_ice` → atm, ocn, bgc, sic, lnd → `smb_init`/`bmb_init` per domain → gases →
all `X_to_cmn`.

Main loop (`climber.f90:371-640`), one iteration per fastest component step:

| # | Block | Guard | Calls |
|---|---|---|---|
| 1 | bnd | always | `bnd_update`, `bnd_to_cmn` (414-415) |
| 2 | geo | `(flag_geo .or. ifake_geo==1) .and. time_call_geo` | `bnd_to_geo`, `ice_to_geo` (all domains), `cmn_to_geo`, `geo_update`, `geo_to_ice` (all domains), `geo_to_cmn` (424-440) |
| 3 | atm | `flag_atm .and. time_call_atm` | `cmn_to_atm`, `atm_update`, `atm_to_cmn` (444-455) |
| 4 | smb | `flag_smb .and. time_call_smb` | `cmn_to_smb` once; per domain: `ice_to_smb` or `geo_to_smb` at start of SMB year, `smb_update`, `smb_to_cmn`, `smb_to_ice` at end of SMB year (459-480) |
| 5 | bmb | `flag_bmb .and. time_call_bmb` | per domain: `cmn_to_bmb`, `ice_to_bmb`/`geo_to_bmb`, `bmb_update`, `bmb_to_cmn` (empty), `bmb_to_ice` (484-507) |
| 6 | ice | `flag_ice .and. time_call_ice` | per domain: `ice_update`, `ice_to_cmn`; exit on error (511-523) |
| 7 | lnd | `flag_lnd .and. time_call_lnd` | `cmn_to_lnd`, `lnd_update_wrapper`, `runoff_merge`, `lakes_update`, `lnd_to_cmn` (527-540) |
| 8-10 | sic, ocn, bgc | flags | analogous; `runoff_to_ocn` before ocean |
| 11 | ch4, n2o, co2 | `flag .and. time_eoy` | |
| 12 | restart | `time_write_restart` | `write_restart` (786-833) |

Coupling is sequential and explicit (lagged): geo uses the ice state from the previous
ice call; ice uses SMB from the last completed SMB year.

### 1.5 Cadences

All cadence logic is in `src/main/timer.f90` as global logical flags recomputed each step
(`time_call_X`, `time_soy_X`, `time_eoy_X`, `time_out_X`, ≈l.300-470). Control
parameters (`nml/control.nml:46-53`): `dt_day_*` per climate component,
`n_year_ice=1`, `n_year_smb=10`, `n_year_geo=10`, `n_accel=1`.

- Climate components run only in years where `mod(year,n_accel)==0` (`timer.f90:309`), i.e. asynchronous acceleration of the ice relative to climate.
- SMB runs one full year (daily) every `n_year_smb` years (`timer.f90:310`, 414-427); ice reuses the last annual SMB in between.
- Ice is called at end of every `n_year_ice` years (`timer.f90:448`). SICOPOLIS integrates `n_year_ice/dtime` internal steps (`ice_model.f90:158-162`); Yelmo integrates to the absolute `time` (`ice_model.f90:106`). Same argument, two different semantics.
- geo at start of year every `n_year_geo` (`timer.f90:457`); `dt_geo` is used inside GIA and VILMA.

### 1.6 On/off and backend switches

| Switch | Where | Mechanism |
|---|---|---|
| `flag_atm, flag_ocn, flag_lnd, flag_sic, flag_bgc, flag_smb, flag_bmb, flag_ice, flag_geo, flag_co2...` | `control.f90:39`, `nml/control.nml:6-19` | runtime `if` in driver and inside coupler routines |
| `ice_model_name = 'yelmo'|'sico'` | `control.f90:42` | runtime string `select case` in `ice_model.f90` |
| `ice_domain_name(1:10)` | `control.f90:44-46, 219-228` | number of non-empty names = `n_ice_domain` |
| `i_geo = 0..3` | `nml/geo_par.nml:7-10` | runtime `if` in `geo_update` (`geo.f90:603-648`): uniform SL, LLRA GIA, VILMA1, VILMA2 |
| `i_smb` (PDD/SEMI/simple), `i_bmb` | `smb_params.f90`, `bmb_params.f90` | runtime `if` inside the component |
| `ifake_atm/ocn/sic/ice/geo/dust` + files | `control.f90:171-184` | runtime choice of prescribed data, used when `flag_X=F` |
| `climber-clim`, `-clim-bgc`, `-clim-ice`, `-clim-bgc-ice` | `Makefile:142-202`, `config/Makefile_climber.mk:245-300` | build-time: link `.X_dummy.f90` stub modules with the same module name |
| `vilma1=0/1`, `vilma2=0/1` | `Makefile:15-27`, `config/common.mk:108-122` | build-time CPP `-DVILMA1/-DVILMA2`; stub `stop 1` if selected at runtime but not compiled (`geo/vilma1.F90:453-464`) |

Note `flag_X` is a mutable global: carbon-cycle spinup flips `flag_atm/lnd/sic` mid-run
(`climber.f90:398-403`).

---

## 2. The ice_model interface

### 2.1 Types

`src/ice/ice_def.f90:33-76`:

```fortran
type ice_class
  type(grid_class)           :: grid               ! native ice grid (coords)
  type(grid_ice_to_cmn_type) :: grid_ice_to_cmn    ! i_lowres(:,:), j_lowres(:,:), ncells(:,:)
  ! out
  logical :: error
  integer,  allocatable :: mask_extent(:,:)        ! "should be external in bnd? fixme"
  real(wp), allocatable :: z_sur, z_sur_std, z_bed_std, z_base, H_ice, calv, Q_b
  ! in
  integer,  allocatable :: mask_ocn_lake(:,:)
  real(wp), allocatable :: z_sl, z_bed, z_bed_fil, smb, accum, runoff, bmb_shlf,
                           temp_s, temp_g, q_geo, H_sed, t_ocn, s_ocn
end type
```

### 2.2 Fields across the boundary

| Field | Dir | Units (canonical) | Source → sink | Notes |
|---|---|---|---|---|
| `H_ice` | out | m | backend → geo, smb, bmb, cmn | |
| `z_sur` | out | m | backend → smb | |
| `z_base` | out | m | backend → bmb | |
| `z_bed` | in/out | m, rel. to current sea surface | geo → backend, and backend → `ice` | declared "in" but overwritten by `*_to_ice` (`ice_model.f90:517, 595`): ambiguous ownership |
| `calv` | out | m(ice)/s, positive loss | backend → cmn | Yelmo: `-(cmb+dmb)/sec_year` (526) |
| `Q_b` | out | m(ice)/s, positive melt | backend → cmn | Yelmo: `-bmb*f_ice/sec_year` (521) |
| `mask_extent` | out | 0/1 | backend | circular for Yelmo (504-505) |
| `error` | out | logical | SICO only | Yelmo never sets it |
| `smb, accum, runoff` | in | m(ice)/s | smb → ice (`coupler.f90:2621-2643`) | Yelmo uses only `smb` |
| `temp_s, temp_g` | in | °C | smb → ice | Yelmo maps `temp_g` to `T_shlf` (`ice_model.f90:557`): questionable |
| `bmb_shlf` | in | m(ice)/s | bmb → ice | sign flipped for SICO (640) |
| `t_ocn, s_ocn` | in | °C, psu | bmb → ice | SICO discharge param only |
| `z_sl` | in | m | geo (`z_sea_lake`, 0 over ocean) → ice | SICO uses it; Yelmo hard-codes `z_sl=0` (541) |
| `z_bed_fil, z_bed_std, z_sur_std` | in | m | geo/smb → ice | SICO-specific needs |
| `q_geo, H_sed` | in | W/m², m | geo → ice | |
| `mask_ocn_lake` | in | int | geo → ice → bmb | passes through ice to reach bmb |

Not exported: grounded fraction / floating mask / grounding-line info. Every consumer
re-derives floating ice from a flotation test on `H_ice`/`z_bed`
(`coupler.f90:2734-2738` bmb, `2870` cmn, `2981` geo volume, `geo.f90:674-700` with a
hysteresis variant). Four copies, not identical.

### 2.3 Procedures and backend selection

`src/ice/ice_model.f90` exposes five public routines (56-61), all taking the backend
name string and a domain index:

| Routine | Lines | Signature |
|---|---|---|
| `ice_init_domains(n_dom, model)` | 182-227 | allocates module-global `ylmo_doms(:)` or `sico_doms(:)` |
| `ice_init(ice, idx_dom, model, time, l_restart, grid, cmn_grid, geo_grid, z_bed_geo, z_bed_rel_geo, h_ice_geo, q_geo_geo, h_sed_geo)` | 229-452 | allocates `ice`, builds ice→cmn index map, remaps geo fields to ice grid, inits backend |
| `ice_update(ice, idx_dom, model, n_year_ice, time, time_out_ice)` | 65-180 | `ice_to_<backend>` → backend update → `<backend>_to_ice` + output |
| `ice_end(ice, idx_dom, model, time)` | 454-490 | |
| `ice_write_restart(model, ice, idx_dom, time, dir)` | 729-769 | file `ice_<backend>_<grid>_restart.nc` |

Backend state lives in module-level arrays (`ice_model.f90:46-54`):

```fortran
type(sico_class),  allocatable, target :: sico_doms(:)
type(yelmo_class), allocatable, target :: ylmo_doms(:)
type(yregions_class), allocatable :: ylmo_regs(:,:)   ! Yelmo-only 1D output buffer
```

`ice_class` carries no handle to its backend; the pair `(model, idx_dom)` is the handle.
Each public routine is a `select case(trim(model))` with one branch per backend and a
`stop` default. Each backend needs four adapter routines in the interface module:
`ice_to_<b>`, `<b>_to_ice`, plus the init and restart branches.

Build: `ice` target always compiles both SICOPOLIS (`ice-sico`) and Yelmo
(`yelmo-static`) (`Makefile:311`); choice is at runtime. One `ice_model_name` for all
domains: cannot run Yelmo for Antarctica and SICOPOLIS for NH in one run.

### 2.4 Backend leakage

| Leak | Where | Why it matters |
|---|---|---|
| Yelmo NetCDF writer (`yelmo_write_step_2D`, 130 lines) and 1D regional buffering (100-step buffer) inside the interface module | `ice_model.f90:111-145, 869-1003` | interface module depends on backend internals (`ylmo%tpo%now%...`) |
| SICOPOLIS (j,i) 0-based transposition and lithosphere smoothing filter inside the adapter | `ice_model.f90:584-601, 621-723` | numerics in the adapter; filter belongs in the backend or in geo |
| `ice_init` takes coupler grid, geo grid and five geo fields | `ice_model.f90:229-248` | ice interface knows host grid topology; it also builds the ice→cmn index map (279-305) which is coupler work |
| `set_ice_id` reads a global region-ID file with hard-coded IDs (EU=1, GRL=2, ...) | `ice_id.f90:37-48` | SICOPOLIS uses `id_mask==2` to switch Greenland-only discharge or freeze Greenland (`ice_sico/boundary_m.f90:191`, `calc_thk_m.f90:1380`) |
| Yelmo par file path fixed to `ice_yelmo_par.nml`, shared by all domains; restart path inside it | `ice_model.f90:346` | per-domain differences rely on Yelmo's `{domain}/{grid_name}` path templating (`yelmo/src/yelmo_defs.f90:1346-1383`) |
| `ice_restart` from control only reaches SICOPOLIS; Yelmo restart is configured in its own namelist | `ice_model.f90:423` vs 357 | two restart conventions |
| SICOPOLIS port `use`s CLIMBER-X `control`/`timer` modules | `ice_sico/sico_model.f90:9-10`, `sico_params.f90:7-8` | backend not standalone |
| Yelmo default boundary values set in the adapter (`T_srf=240`, `mask_ice=2`, ...) | `ice_model.f90:357-375` | |

### 2.5 Good and bad

Good:
- A neutral exchange type on the native ice grid with documented canonical units and sign conventions; conversion (yr↔s, K↔°C, W↔mW, sign of melt/calving) is done once in the adapter.
- Per-domain instancing is natural: `ice(n)` arrays, everything else loops `do n=1,n_ice_domain`.
- Adding a backend = one adapter pair + branches; consumers (`coupler.f90`) never see backend types.
- The physical-constants comparison at init (`ice_model.f90:383-389`, `phys_const_compare`) makes cross-model constant disagreements visible.

Bad:
- String dispatch duplicated in 5 places; module-global backend arrays; no object owns its backend.
- `ice_class` mixes import and export state in one type with one ambiguous field (`z_bed`) and some fields that are SICO-only (`accum`, `runoff`, `z_bed_fil`, `z_sur_std`, `t_ocn`, `s_ocn`) — the interface is the union of backend needs, not a minimal contract.
- No grounded/floating fraction, no mass-balance budget terms (dH/dt, SMB applied vs. requested) exported.
- Time semantics differ between backends (step count vs. target time).
- Error reporting only for one backend.

---

## 3. Regional-in-global mapping

### 3.1 Map machinery

All remapping uses fesm-utils `coords` (`map_init`, `map_field`) with methods
`"con"` (conservative), `"bil"`, `"nn"`, optional Gaussian filtering
(`filt_method="gaussian", filt_par=[sigma_km, dx]`), and weight caching:
`map_init(..., gen=map_gen, fldr="maps", load=.TRUE.)` writes/reads SCRIP files keyed
by method and grid names (`maps/scrip-bil_CMN-5x5_NH-32KM.nc`, ...). Grid names therefore
act as cache keys; some are derived from filename substrings with a hard-coded offset
(`geo_grid.f90:87`, `geo.f90:518`: `name=trim(file(7:ppos))` assumes an `input/` prefix).

Map objects are stored inconsistently: in the component (`smb%maps_cmn_to_ice`,
`bmb%maps_cmn_to_ice`), as `save` locals inside coupler routines
(`ice_to_geo` 2917-2935, `geo_to_ice` 3032-3049, `bnd_to_geo` 3105-3106), or rebuilt on
every call (`geo_to_smb` 2581, `geo_to_bmb` 2769, which reload from disk each time).

### 3.2 Directions

| From → To | Fields | Method | Where |
|---|---|---|---|
| cmn (5°) → smb/ice grid | ~35 atmospheric fields | `map_method` (nml) with Gaussian filter | `smb_model.f90:118, 213-460`; map built at 774 |
| cmn (5°) → bmb/ice grid | 3D ocean T/S at depth levels, lake T/S | bilinear | `bmb_model.f90:411` |
| ice → cmn (5°) | calving, basal melt (grd/flt) as kg/s per coarse cell | index aggregation via `grid_ice_to_cmn%i_lowres/j_lowres` (cell-centre-in-box test), area-weighted sum, time relaxation | `coupler.f90:2833-2899` |
| smb → cmn (5°) | monthly runoff, ice melt (kg/s); daily albedos (mean over ice cells) | same index aggregation (`grid_smb_to_cmn`) | `coupler.f90:2429-2523` |
| ice → geo hi-res | `H_ice` | bilinear for thickness + nearest-neighbour ice mask threshold 0.5, overlay with `reset=.false.` | `coupler.f90:2906-2964` |
| geo hi-res → ice | `z_bed`, `z_bed_fil`, `z_bed_std`, `z_sl` (sea/lake level), `q_geo`, `H_sed`, `mask_ocn_lake` | conservative (+Gaussian 50-100 km for some), nn for masks | `coupler.f90:3023-3089` |
| geo hi-res → geo coarse | fractions, mean elevations per surface type, runoff routing | block aggregation (`hires_to_lowres`) | `geo.f90:862`, `geo/hires_to_lowres.f90:46` |
| geo hi-res → smb/bmb (no ice model) | `z_topo`, `H_ice`, `z_bed_std`, ice-shelf mask | conservative | `coupler.f90:2567-2614, 2750-2800` |

The ice→cmn index map is computed by a brute-force O(N_ice × N_cmn) loop in
`ice_init` (`ice_model.f90:280-305`), duplicated verbatim in `smb_init`
(`smb_model.f90:745-772`) and in `bnd/fake_ice.f90` (own copy of
`grid_ice_to_cmn_type`, l.12-16). It assigns each fine cell to exactly one coarse cell by
its centre, so it is conservative only up to the projection-area approximation and
ignores fine cells straddling coarse-cell edges.

### 3.3 Masks and the "override" pattern

Two coarse masks record coverage (`cmn%mask_ice`, `cmn%mask_smb`, set in
`climber.f90:221-223, 303-305`). `runoff_merge` (`coupler.f90:3765-3835`) then picks,
per coarse cell:

```fortran
if (cmn%mask_smb(i,j).eq.0) then  ! not covered: use land-model ice runoff + prescribed dH/dt freshwater
   cmn%runoff_ice(i,j) = cmn%runoff_ice_l(i,j) - scale_dhdt_ice*cmn%dhdt_ice_i_mon(i,j,mon)
else                              ! covered: use SMB-model runoff
   cmn%runoff_ice(i,j) = scale_runoff_ice*cmn%runoff_ice_i_mon(i,j,mon)
endif
if (cmn%mask_ice(i,j).eq.0) then  ! not covered: land-model calving, no basal melt
...
```

Similarly `cmn_to_lnd` replaces the land model's ice albedo with SEMI albedo where
`mask_smb==1` (`coupler.f90:1471-1472`). This is the global/regional dichotomy in its
simplest form: a global fallback everywhere, replaced by the regional component where
the regional domain claims the coarse cell.

On the hi-res geo grid the same idea is a field overlay: `bnd_to_geo` writes the
prescribed global ice thickness (fake ice), then `ice_to_geo` maps each domain on top with
`reset=.false.` (cells outside the domain keep their previous value; fesm-utils
`coords/mapping.f90:927-993`). `nml/control.nml:363-365` documents it: interactive ice
sheets overwrite the fake ice thickness only over the ice-model domain(s).

### 3.4 Overlap and conservation

- Overlapping domains are not supported in any explicit way. In `ice_to_geo` the last domain wins in overlap areas, while the volume-above-floatation sum loops all domains (`coupler.f90:2977-2991`), so overlap double-counts sea-level contribution. In `ice_to_cmn`/`smb_to_cmn` each domain relaxes the coarse value toward its own contribution, so two domains sharing a coarse cell overwrite rather than add. Implicit rule: domains must not share coarse cells.
- Flux conservation ice→cmn: area-integrated sums in kg/s (good), but then time-relaxed with `relax_calv/relax_bmelt/relax_run` (2478-2479, 2886-2888), which delays rather than conserves on short times, and scaled by tunable `scale_calving_ice`, `scale_runoff_ice` (3812, 3825).
- Ice volume ↔ geo: thickness remapped bilinearly onto 0.125° and masked by nearest-neighbour; there is no check that remapped volume matches native volume (debug prints commented out at 2959-2960).
- `ice_init` takes initial `H_ice`, `z_bed` from the 0.125° geo grid (`ice_model.f90:318-328`), so a high-resolution domain (e.g. GRL-4KM) is initialised from ~14 km data unless Yelmo is restarted.

---

## 4. Geo / sea level / isostasy

### 4.1 What is global, what is regional

| Quantity | Owner | Grid | Notes |
|---|---|---|---|
| Bedrock `z_bed` (relative to current local sea surface) | geo | global hi-res | updated by GIA/VILMA or prescribed (fake_geo) |
| Global ice thickness `h_ice` | geo | global hi-res | background from bnd + overlay of all ice domains |
| Ocean/land/lake/floating masks, flood-fill ocean connectivity, lakes, runoff routing | geo | global hi-res → coarse | `geo_update` (583-970) |
| Scalar sea level `sea_level` | geo (or bnd) | scalar | bnd if `flag_geo=F`; `+= d_sea_level` for `i_geo=0,1`; ocean-area mean of `rsl` for `i_geo=2,3` (`geo.f90:928`) |
| `d_sea_level` from ice volume | coupler (`ice_to_geo`) | scalar | `-(ΔV_af)*(ρi/ρw)/A_ocn` (`coupler.f90:3007`), only from modelled domains |
| Relative sea level `rsl` | VILMA | global | `i_geo=2,3` |
| Sea/lake surface `z_sea_lake` | geo | hi-res | 0 over ocean, lake level over lakes |
| q_geo, H_sed | geo (read at init) | hi-res | passed to ice |
| Fractions, elevations for climate | geo coarse | 5° | `geo_to_cmn` |
| Ice dynamics, SMB, BMB | ice/smb/bmb | regional | per domain |

### 4.2 Datum convention

`geo%hires%z_bed` is always relative to the current local sea surface:
- `i_geo=0`: `z_bed = z_bed - d_sea_level` (`geo.f90:612`).
- `i_geo=1` (LLRA): shift to absolute, relax toward `z_bed_rel - ρi/ρa·H`, shift back (`gia.f90:39-75`).
- `i_geo=2,3`: `z_bed = z_bed_eq - rsl` (`vilma1.F90:135`; VILMA2 the same internally).

Consequently the ice model is told "sea level is 0" (Yelmo, `ice_model.f90:541`) or
gets only lake levels above 0 (SICO via `z_sl = z_sea_lake`). This is consistent and
avoids passing a spatial sea-level field to the ice, but it means the ice interface has
no explicit datum and the meaning of `z_bed` depends on geo internals.

### 4.3 VILMA1 and VILMA2 integration

Both are called only from `geo_update` (`geo.f90:626-648`) with the same contract:

```
in : h_ice on geo hi-res grid (global, merged)
out: rsl, z_bed on geo hi-res grid
```

- VILMA1 (`src/geo/vilma1.F90`): wrapper around a prebuilt static library (`src/vilma1/lib/vega_pism.a` + `.mod` files; source in `src/vilma_src/`). The wrapper owns the remap to the VILMA Gauss grid (`vilma1.F90:205-206`), writes ice thickness to a NetCDF file that VILMA1 then reads (`vilma1.F90:105-122`), calls `time_evolution`, and copies restart files with shell `cp` (`vilma1.F90:425-432, 626-632`). File-based coupling, module globals, `#ifdef VILMA1` around everything.
- VILMA2 (`src/geo/vilma2.F90`, 167 lines): thin. `vilma_par_load(se%par, geo_par.nml, defaults_file)`, `solid_earth_init(se, z_bed_eq, h_ice_eq, grid=geo_grid, h_ice_init)`, `solid_earth_update(se, h_ice, dt_yr)` → `se%rsl`, `se%z_bed`, `solid_earth_finalize`, `vilma_restart_write/read`. The library's own `vilma2/src/vilma_coupling.f90:1-30` states the contract: the model owns its Gauss grid and the host↔Gauss remap; the host reads outputs on the grid it supplied; only the smooth `rsl` is remapped back and `z_bed = z_bed_eq - rsl` is reconstructed on the host's high-resolution bed, preserving detail. VILMA2 also offers a `solver="v1"` backend behind the same API.
- Build: `-DVILMA1/-DVILMA2` with runtime fail-fast stubs (`vilma1.F90:453-464`, `vilma2.F90` analog).

### 4.4 Interface between geo and ice

geo → ice (`geo_to_ice`, every geo call): `z_bed`, `z_bed_fil`, `z_bed_std`, `z_sl`,
`mask_ocn_lake`, `q_geo`, `H_sed`. Ice → geo (`ice_to_geo`): `H_ice` only, plus `z_bed`
and `z_sl` read for the volume-above-floatation sum. Ice models run without their own
isostasy (Yelmo core has none; SICO gets `zl_neu` forced, `ice_model.f90:626`).

So in CLIMBER-X the solid Earth and sea level are a single global component sitting
beside the ice domains, receiving a global merged load and returning bed + sea level.
That is exactly the regional/global split kryos needs, but the sea-level budget
(`d_sea_level`) is computed in the coupler rather than in geo, and only for modelled
domains.

---

## 5. SMB and BMB: between global climate and regional ice

- SMB (`src/smb/`): one global input object `smb_in` on the coupler grid (`smb_def.f90:36-76`, filled once per step by `cmn_to_smb`, `coupler.f90:2209`), and one `smb(n)` per ice domain on the ice grid. Each `smb(n)` holds its own `maps_cmn_to_ice` and downscales ~35 fields (lapse-rate and topographic corrections, `smb/downscaling.f90`), then runs PDD, SEMI (energy balance + snowpack) or `smb_simple` (`i_smb`). Optional `fake_atm_hires` reads regional climate model data directly on the SMB grid (`smb_model.f90:472-487, 851-853`): a data backend inside the component.
- SMB → ice: annual means converted to m(ice)/s (`coupler.f90:2621-2643`). SMB → cmn: monthly runoff/ice melt and daily albedo back to the climate (`2429-2523`).
- Without an ice model, SMB/BMB take geometry from geo hi-res (`geo_to_smb`, `geo_to_bmb`): the regional component works with a "data" ice sheet.
- BMB (`src/bmb/`): per domain on the ice grid; 3D ocean T/S (coarse) and lake T/S are remapped bilinearly to the domain, extrapolated to ice-base depth, and turned into shelf melt (`i_bmb` 0/1/2) and discharge-param T/S (`t_disc`, `s_disc`). `bmb_to_cmn` is empty (`coupler.f90:2699-2711`): melt reaches the ocean via the ice model (`ice%Q_b` → `ice_to_cmn`).
- SMB and BMB are therefore regional components on the same grid as the ice, but separated from it. They are not inside the ice model (good), and could in principle run with any ice backend (good). They are, however, CLIMBER-specific in their inputs (`smb_in_class` mirrors CLIMBER atmosphere variables) and `use` global `timer`/`control` modules.

---

## 6. Dummy / data component patterns

Three distinct mechanisms coexist:

| Mechanism | Example | How it works | Assessment |
|---|---|---|---|
| Link-time stub (dummy) | `src/ice/.ice_model_dummy.f90`, `src/smb/.smb_model_dummy.f90`, `src/bmb/.bmb_model_dummy.f90`, `src/bgc-dummy/.bgc_*.f90` | Hidden file defines a module with the same name and same public signatures, bodies empty. Selected by Makefile target; `climber.f90` is compiled several times against different `.mod` sets (`config/Makefile_climber.mk:245-300`), as is `coupler.f90` (`coupler.o` vs `coupler_nobgc.o`) | Purpose: build without Yelmo/HAMOCC dependencies. The stub `_def` types must still contain every field the coupler touches, so they are full copies that drift: `.bgc_def.f90` already differs from `bgc/bgc_def.f90` (missing M4AGO fields, own `kt/ks` constants). `docs/TODO.md:8-9` already proposes replacing them. Not a model to copy |
| Prescribed data (fake) | `src/bnd/fake_atm.f90`, `fake_ocn`, `fake_sic`, `fake_lnd`, `fake_ice`, `fake_geo`, `fake_dust`; `smb/fake_atm_hires.f90` | Held in `bnd_class` (`bnd.f90:83-116`), updated in `bnd_update` when `flag_X=F` (247-285), with constant (`ifake_X=0`) or time-interpolated (`ifake_X=1`) NetCDF input. `bnd_to_cmn` writes them into the same `cmn` slots the real component's `X_to_cmn` would (e.g. `cmn%t2m`, `cmn%t_ocn_bmb`, `coupler.f90:3505-3560`) | Conceptually right: consumers are agnostic to whether a field was computed or read. In practice the fake readers live in one central `bnd` module, not next to the component; each reimplements time interpolation (`fake_ice.f90:125-150` etc.); some derived physics (rain/snow split) is done in the coupler; and "fake" inputs are also used as nudging/forcing data for active components (`ocn_restore_*`, `atm_fix_tau`, SMB bias correction), blurring the role |
| CPP backend gate | `src/geo/vilma1.F90`, `vilma2.F90` | `#ifdef VILMA1` around the whole body; runtime `i_geo` selects; stub stops with an explanatory message if not compiled | Clean enough for optional heavy dependencies; the fail-fast message is good practice |

Neither "dummy" nor "fake" is a backend in the sense kryos wants (a runtime-selectable
implementation of a component interface). The fake data readers are the closest, but
they are not organised per component.

---

## 7. What is good and worth keeping

1. **Named grid registry.** `nml/ice_grids.nml` holds one group per grid (`&GRL-16KM`, `&ANT-32KM`, ...) with projection parameters; `ice_grid_init(name)` builds a `coords` grid from it. Domains are chosen by listing names in `ice_domain_name`. Why: domains are data, not code; adding a grid needs no rebuild; grid names double as map-cache keys and output file suffixes (`ice_<grid>.nc`, `smb_<grid>_restart.nc`).
2. **Domain bundle on a shared grid.** ice, smb and bmb for domain `n` share `ice_grid(n)`; ice↔smb and ice↔bmb exchanges are plain array copies with unit conversion, no remapping (`coupler.f90:2530-2643, 2719-2826`). Why: the regional components that must be consistent with the ice geometry live on its grid; only regional↔global crossings remap.
3. **Neutral ice exchange type with canonical units.** All consumers use `ice_class` in SI-ish canonical units; backend conventions are converted in exactly one place. Why: consumers never change when a backend is added.
4. **Background + overlay for global fields.** Global hi-res ice thickness = prescribed background, overwritten by each domain where it has coverage (`reset=.false.`). Coarse-grid fluxes = global fallback model, replaced where `mask_ice/mask_smb` say a regional component covers the cell. Why: any subset of the globe can be modelled interactively while global components still see a complete field. This is the cleanest available answer to kryos's regional/global dichotomy.
5. **Global solid Earth with a narrow contract.** `h_ice(global) → rsl, z_bed`. VILMA2's `vilma_coupling` shows the best version: model owns its native grid and remap; host passes its grid at init; outputs come back on host grid; high-res bed detail preserved by reconstructing `z_bed = z_bed_eq - rsl` on the host grid. Why: the host never learns spectral/Gauss details and backends (LLRA, VILMA1, VILMA2, VILMA2-v1) are interchangeable behind it.
6. **Cached remapping weights** via `coords` (`maps/` SCRIP files keyed by grid names). Why: expensive conservative weights computed once across runs.
7. **Ice-free/ice-model-free modes of regional components.** SMB/BMB can run on geometry from geo when `flag_ice=F` (`geo_to_smb`, `geo_to_bmb`). Why: offline SMB experiments without the ice model.
8. **Per-component restart files per domain**, named by grid, in `restart_out/year_<YYYY>/`, with a crash-restart dump on error (`climber.f90:686-712`). Why: easy partial restarts (restart flags per component in `control.nml:388-399`).
9. **Physical-constant comparison at coupling time** (`ice_model.f90:383-389`). Why: silent constant mismatches between host and ice model are a classic source of drift.
10. **Fail-fast stubs for optional backends** with an actionable error (`vilma1.F90:453-464`).

---

## 8. What is legacy / awkward to avoid

1. **Giant flat hub (`cmn_class`, ~200 fields) + 4600-line pairwise coupler.** Every new exchange touches the hub type and the coupler; `coupler.f90` depends on all `_def` modules, so dummies must replicate full types. The hub has no notion of who produces a field, on which grid, at which cadence, or whether it is valid. kryos should give each component typed import/export states and let a thin coupler route fields.
2. **Physics in the coupler.** Sea-level budget (`ice_to_geo`), runoff merging and routing, lake balance, rain/snow split, flux time relaxation and tuning scalings live in `coupler.f90`. They cannot be swapped or tested in isolation.
3. **String dispatch + module-global backend arrays** (`ice_model.f90:46-47, 94, 192, 331, 468, 744`). Five `select case` blocks to edit per new backend; state not owned by an object; only one backend per run.
4. **Interface = union of backend needs.** SICO-only fields (`accum`, `runoff`, `z_bed_fil`, `z_sur_std`, `t_ocn`, `s_ocn`) sit in the shared type; Yelmo ignores several inputs (`z_sl`, ocean T/S) and misuses one (`temp_g`→`T_shlf`). Missing outputs (grounded fraction, masks) force consumers to re-derive flotation four times.
5. **Ambiguous field ownership.** `ice%z_bed` written by both geo and backend; `mask_extent` "should be external" (fixme at `ice_def.f90:51`).
6. **Host grid knowledge inside the ice interface.** `ice_init` takes `cmn_grid`, `geo_grid`, and geo fields and computes the ice→cmn index map (`ice_model.f90:229-328`); the same brute-force map is copied into `smb_init` and `fake_ice`.
7. **Build-time variants by module shadowing.** Same source compiled against different `.mod` sets (`climber_clim.o`, `climber_clim_ice.o`, ...), hidden dot-files, drifting dummy types (`.bgc_def.f90` vs `bgc_def.f90`). Three different selection mechanisms (Makefile targets, CPP, runtime strings/integers).
8. **Global mutable control state.** `flag_*`, `time_call_*`, `year_now` etc. are module globals `use`d everywhere, including inside backends (`ice_sico` uses `control`/`timer`; `bmb_model.f90:32-33`). Flags are flipped mid-run for spinup (`climber.f90:398-403`). Components are not reusable outside CLIMBER-X.
9. **Hard-coded grids and geography.** Coupler grid as parameters (`climber_grid.f90:39-40`), Panama closure by cell indices `f_ocn0(16,23)` (`geo.f90:491-494`), Greenland fix cell `(32,34)` (`geo.f90:718-725`), Hudson/Baltic closure by lon/lat boxes (`geo.f90:655-666`), region IDs `ice_id_grl=2` (`ice_id.f90:38-44`) driving SICO physics, grid names from `filename(7:ppos)`.
10. **Implicit no-overlap assumption** for multiple domains (§3.4), and sea-level contribution only from modelled domains.
11. **Mixed time semantics** (`n_year_ice` as step count for SICO, target time for Yelmo; `n_year_smb` must be a multiple of `n_accel`, `timer.f90:138-145`); cadence flags are 40+ globals rather than per-component schedules.
12. **Backend output in the interface module** (`yelmo_write_step_2D`, 1D buffer) and restart directories pre-created by running the whole timer loop once (`climber.f90:719-783`).
13. **File/shell-based coupling for VILMA1** (NetCDF hand-off, `execute_command_line('cp ...')`).
14. **SMB inputs tied to one climate model.** `smb_in_class` is literally CLIMBER's atmosphere variable list; swapping the climate source means rewriting `cmn_to_smb`.

---

## 9. Further lessons for a from-scratch kryos design

1. **Make "domain" a first-class object.** CLIMBER-X has the right ingredients scattered: a named grid, a list of regional components on it (ice, smb, bmb), per-domain restart/output naming, per-domain maps to global grids. In kryos a `domain` object should own its grid, its regional component instances (each with its chosen backend), its maps to every global component grid, and its coverage mask on those grids. kryos's primary use (one domain) is then the `n=1` case of the CLIMBER-X design, and multi-domain comes for free.
2. **Separate the two kinds of global component.** In CLIMBER-X the climate grid (`cmn`) and the global hi-res geography (`geo`) are two different "global" layers with different consumers. kryos will have the same: global forcing/climate (often data) and global solid Earth/sea level (VILMA, ELVA-like). Each global component should declare its grid; regional components import from it through maps owned by the domain.
3. **Regional→global aggregation needs explicit rules:** background provider, coverage/priority per domain, conservative accumulation (sum, not overwrite), and a volume check. CLIMBER-X shows the mechanism (`reset=.false.` overlay, coverage masks) but leaves the rules implicit.
4. **Sea-level budget belongs to the sea-level component.** Ice components should export grounded ice mass/volume above floatation (or `H_ice`, `f_grnd`) and the iso/SL component should compute barystatic and relative sea level. CLIMBER-X computes `d_sea_level` in the coupler from modelled domains only.
5. **Fix the datum in the interface.** CLIMBER-X works because everything is relative to the local sea surface and ice always receives `z_sl=0`. kryos should pass `z_bed` and `z_sl` both in an absolute datum (yelmox/isostasy convention) so that a domain-only isostasy model, a global SLE model and a data backend are interchangeable.
6. **Export masks and fractions from the ice model** (`f_grnd`, `f_ice`, `mask_bed`), so no consumer re-implements flotation.
7. **Model-owned grids and remaps (VILMA2 style)** for components with exotic native grids (spectral, Gauss, unstructured): the component takes a host grid at init and returns fields on it. For grid-native components, the domain-level map cache (fesm-utils `coords`) is enough.
8. **One backend mechanism**, runtime-selected per component per domain, with a `data` backend for every component reading the same export fields from files with shared time-interpolation code. Use CPP only to exclude heavy optional libraries, with a fail-fast stub as in `vilma1.F90`.
9. **Cadence per component**, expressed in the component's configuration (dt, phase), with the driver deriving call times; avoid 40 global logicals. Keep explicit target-time semantics (`update(t_end)`), as Yelmo does, not step counts.
10. **Keep components free of host globals.** `ice_sico`, `smb`, `bmb` all `use control`/`timer`; that is what prevents reuse. Pass time and paths in; let each component read its own namelist group, ideally with the Yelmo/VILMA2 `defaults_file` convention (complete defaults file + minimal overrides) and `{domain}/{grid_name}` path templating.
11. **Units/sign conversions at the adapter boundary only**, with a documented canonical table (CLIMBER-X `ice_def.f90` comments plus the conversions in `ice_to_yelmo`/`yelmo_to_ice` are a good start). Add the physical-constant comparison per backend.
12. **Errors and diagnostics per component instance**, reported uniformly (CLIMBER-X has `ice%error` only for SICO, crash restart dump for all).
13. **The virtual-cell design draft** (`docs/design/virtual-cells.md`) is a useful parallel: single-column physics, aggregation with conservation at the boundary, coupling contract unchanged at the coarse cell. The same principle applies to SMB/snowpack in kryos: keep column physics grid-agnostic and put all horizontal mapping in the domain/coupling layer.

---

## Appendix A: key file index

| File | Lines | Role |
|---|---|---|
| `src/main/climber.f90` | 836 | driver, init order, main loop, restart writing |
| `src/main/coupler.f90` | 4613 | `cmn_class` (128-367), all exchange routines, runoff/lakes physics |
| `src/main/control.f90` | 437 | `control.nml` flags, fake switches, restart flags, domain names |
| `src/main/timer.f90` | 648 | global cadence flags |
| `src/main/climber_grid.f90` | 140 | fixed 5° grid, `ice_grid_init` from `ice_grids.nml` |
| `src/ice/ice_def.f90` | 81 | `ice_class` |
| `src/ice/ice_model.f90` | 1005 | ice interface, Yelmo/SICO adapters, Yelmo output |
| `src/ice/ice_id.f90` | 91 | global region-ID mask |
| `src/ice/.ice_model_dummy.f90` | — | link-time no-op ice interface |
| `src/ice_sico/sico_model.f90` | 1928 | embedded SICOPOLIS port (`sico_init/update/end/write_restart`) |
| `src/geo/geo.f90` | 1166 | global hi-res geography, sea level, GIA/VILMA dispatch |
| `src/geo/geo_def.f90` | 157 | `geo_class`, `hires_type` |
| `src/geo/gia.f90` | 147 | LLRA |
| `src/geo/vilma1.F90` | 643 | VILMA1 wrapper (file-based) |
| `src/geo/vilma2.F90` | 167 | VILMA2 wrapper (API-based) |
| `vilma2/src/vilma_coupling.f90` | — | VILMA2 host coupling contract |
| `src/smb/smb_model.f90`, `smb_def.f90` | 1715, 431 | regional SMB, downscaling |
| `src/bmb/bmb_model.f90` | 597 | regional shelf/lake basal melt |
| `src/bnd/bnd.f90`, `src/bnd/fake_*.f90` | 659, ~2100 | boundary forcing and prescribed component data |
| `nml/control.nml`, `nml/ice_grids.nml`, `nml/geo_par.nml`, `nml/ice_yelmo_par.nml` | — | configuration |
| `config/Makefile_climber.mk`, `Makefile` | — | build variants, dummy linking, VILMA CPP gates |
| `docs/TODO.md`, `docs/design/virtual-cells.md` | — | owners' own notes on dummies and on subgrid design |
