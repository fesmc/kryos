# kryos core — phase 1 design note

Status: proposal, 2026-10-24, for review. No code. Scope: `src/core` for
phase 1 of `docs/design.md` §14 (one domain; ice = yelmo, geo = fastisostasy,
climate = snapshots, surface = chion, shelf = marine_shelf, geothermal and
sediments = data; a minimal globe with `forcing` (tsgen) and geo's global part
(bsl)). Gate: `yelmox_Greenland` and `yelmox_Antarctica` bit-identical.

Two notes written in parallel are referenced, not repeated: `contracts-design.md`
(vocabularies, abstract types + wrappers, capabilities, wiring tables with lags
and transforms, hub formulas, geo gather/update/scatter, forcing scalars) and
`adapters-design.md`. §10 lists what the core assumes from them. fesm-utils
names are those on branch `phase0-kryos` (`fesm_clock`, `fesm_fields`,
`fesm_log`, `fesm_datasrc`, `coupler`, `htopo`, `nml`, `phys_constants`).


**Review outcome (2026-10-24, Alex).** All recommendations accepted: Antarctica
gate with `relax = opt = False` (shipped base = phase 2 gate); geo initialised
from the ice instance's reference exports (design decision 4 refined); hub in
`sp` through the phase 1 gate, `dp` as a registered rebaseline after; wiring
transforms in the consumer's declared precision with yelmox's operand order;
`bmb_shlf` flips in Antarctica quantified on R1 output, registered as a
deviation if any; gate flags `-O2 -fp-model precise`, production `-Ofast` once
per phase as information; lapse-rate transform without the identical-grid
shortcut; yelmox `kryos_*` module names may be reused; `Namelist.jl`/`Resolve.jl`
move to FesmUtils.jl; multi-schema `nml_validate` and public `coupler_method`
added to fesm-utils; `dwdt_corr` edge dropped. Gate runs on Levante only (no
local data sync). Names and signatures: `contracts-design.md` wins where the
notes differ (`M_ice`, `dV_bsl`, `to_ann`/`so_ann`/`dto_ann`, `z_srf_clim`,
`init_ctx_t`, `step_time_t`, `PASS_*` + `PASS_RESTART`, `*_slot_t`, geo
`gather_into`/`geo_reduce`/`geo_scatter`); `pr` converts with `sec_year/86400`
via `fesm_units`, not 360.

## 1. Module map

| Module (`src/core/`) | Owns | Uses (kryos) |
|---|---|---|
| `kryos_config` | `config_t`: par file, `input/kryos_defaults.nml`, library list; read with template, `config_find` (`nml_find_group`), `config_validate` | — |
| `kryos_component` | `component_t`, `init_ctx_t`, `step_time_t`, pass constants (contracts-design) | — |
| `kryos_<comp>` (components/) | abstract `<comp>_t`, wrappers, factory submodule `kryos_<comp>_factory` | `kryos_component` |
| `kryos_wiring` | `wire_t`, the fixed tables (content: contracts-design), `wire_resolve`, `wire_exchange` on two field states, `wiring_describe` | `kryos_component` |
| `kryos_domain` | `domain_t`, grids + default chain, hub, coupler, per-grid references, instance slots, `exchange`, `domain_hub_update`, `domain_check_reference` | config, wiring, all `kryos_<comp>` |
| `kryos_globe` | `globe_t`: forcing slot, geo global slot (phase 3: global grids, coupler, global hub, ocean) | config, wiring, `kryos_forcing`, `kryos_geo` |
| `kryos_output` | `output_t`: streams (file × schedule), dispatch, hub and timing files | domain, globe |
| `kryos_restart` | bundle naming, manifest, write/restore order | domain, globe |
| `kryos_run` | `run_t`, `run_init`, `run_step`, `run_done`, `run_end` | everything above |
| `kryos_run_sequence` | submodule of `kryos_run`: the cold start, the restore, the step sequence (§3) | — |

No cycles: fesm-utils → config, component → components/backends (factories
as submodules) → wiring → domain, globe → output, restart → run. The core never
`use`s a backend; backends never `use` config, domain, globe or run. Wrappers
sit below the domain, so they take an `init_ctx_t` and a `step_time_t` (§2),
not `run_t`/`domain_t`. Against design §13: no `kryos_scheduler` (`fesm_clock`
schedules; cadences live in `component_t`), `kryos_config` and the submodule
added, `kryos_fesm_pending` deleted; logging is a convention (§8).

## 2. Types

```fortran
type :: run_t                                  ! kryos_run
    type(config_t)         :: cfg
    character(len=512)     :: out_dir = "./", restart = "None"
    type(phys_const_class) :: cnst              ! the run's one record (§9)
    type(clock_t)          :: clk
    type(schedule_t), allocatable :: sched(:)   ! output schedules by group name; one is the restart schedule
    integer                :: i_rst = 0
    type(globe_t)          :: globe
    type(domain_t), allocatable :: dom(:)       ! N = 1 in phase 1; no polymorphic arrays
    type(output_t)         :: out
    type(timer_class)      :: tmr               ! per-stage wall time (yelmox timer table)
end type

integer, parameter :: C_ICE = 1, C_GEO = 2, C_CLIMATE = 3, C_SURFACE = 4, C_SHELF = 5, &
                      C_GEOTHERMAL = 6, C_SEDIMENTS = 7, NCOMP = 7
integer, parameter :: GRID_DEFAULT(NCOMP) = [0, C_ICE, C_ICE, C_CLIMATE, 0, C_ICE, C_ICE]  ! 0 = hub

type :: grid_ref_t                             ! decision 4: reference topography on one grid
    character(len=64)     :: grid = ""
    real(dp), allocatable :: z_bed(:,:), H_ice(:,:), z_srf(:,:), z_bed_sd(:,:)
end type

type :: domain_t                               ! kryos_domain
    character(len=64)   :: name = ""            ! = its par-file group, e.g. "Greenland"
    character(len=64)   :: grid_hub = ""
    character(len=64)   :: inst(NCOMP) = ""     ! instance (group) name per component; "" = absent
    character(len=64)   :: grid(NCOMP) = ""     ! resolved native grid per component
    type(htopo_class)   :: hub                  ! sp in phase 1 (design §12 row 8 not done)
    type(field_state_t) :: hub_exports          ! z_srf, H_ice, z_bed, f_grnd, z_sl, z_srf_ref on grid_hub
    real(dp)            :: hub_t = 0.0_dp       ! time of the last hub update (lag checks)
    type(coupler_class) :: cpl                  ! grids by name, map cache (maps/)
    type(grid_ref_t), allocatable :: refs(:)    ! one per distinct grid, hub first
    class(ice_t),        allocatable :: ice
    class(geo_t),        allocatable :: geo
    class(climate_t),    allocatable :: climate
    class(surface_t),    allocatable :: surface
    class(shelf_t),      allocatable :: shelf
    class(geothermal_t), allocatable :: geothermal
    class(sediments_t),  allocatable :: sediments
    type(wire_t), allocatable :: wires(:)       ! rows resolved for this domain (active rows only)
end type

type :: globe_t                                ! kryos_globe
    class(forcing_t), allocatable :: forcing     ! tsgen: coefficients on the forcing time
    class(geo_t),     allocatable :: geo         ! global part: fastisostasy bsl (scope "global")
    type(wire_t), allocatable :: wires(:)        ! globe <-> domain scalar rows
end type
```

Assumed from contracts-design (fields the core fills or reads):

```fortran
type :: init_ctx_t                             ! one per instance, built by the core, freed after init
    character(len=512) :: par_file, library(8)
    character(len=64)  :: instance, component, backend, domain, grid_name, scope
    type(grid_class)   :: grid                  ! from maps/grid_<name>.txt
    type(regions_class):: reg                   ! domain regions on this grid (nn), as yelmox domain_regions_load
    type(grid_ref_t)   :: ref                   ! reference topography on this grid
    type(phys_const_class) :: cnst
    real(dp)           :: time, time_rel, dtt   ! t0
end type
type :: step_time_t
    real(dp) :: t0, t1      ! model-time interval (cadence_interval); t0 = t1 = time_init during init
    real(dp) :: t_eval      ! t1, or clock_time_rel at t1 when the instance has clock = "forcing"
    real(dp) :: dtt
end type
integer, parameter :: PASS_GEOMETRY = 1, PASS_BOUNDARY = 2, PASS_STATE = 3, PASS_FINAL = 4
```

## 3. Run lifecycle

```fortran
call run_init(run, par_file)            ! cold start or restore, then output at time_init
do while (.not. run_done(run))          ! run_done = clock_finished(run%clk)
    call run_step(run)                  ! advance, sequence, output + restart at the new time
end do
call run_end(run)
```

Writing output at the end of `run_step` (and of `run_init`) is the same order
as yelmox's top-of-loop writes; the final state and a bundle at `time_end`
follow from `schedule_due` being true once the clock is finished.

### 3.1 `run_init`

1. `config_init` ([kryos], §5); `log_set_level`; `phys_const_load` (§9); `nc_set_lossless`.
2. Clock: the `[kryos] clock` group read through `config_read` (template
   `[clock]`), then `clock_init(clk, time_init, time_end, dtt, calendar, year,
   method=, time_const=, time_ref=)` (argument form, so the group may have any name).
3. `globe_init` before the domains (yelmox: `bsl_init` + `bsl_update`, then
   `coefs_init`; the climate validates the coefficients it names): forcing,
   then geo global. Each instance: factory, `<comp>_init(inst, ctx)`,
   `cadence_init(inst%cadence, clk, inst%name, dt | every_step)`.
4. `domain_init` per domain: groups, grids, `htopo_init(dom%hub, par, name,
   topo, regions, name, grid_hub, cnst)`, `coupler_init(dom%cpl, "maps")`,
   `refs` (hub → grid by the coupler's sp path, §10 R4), instances in yelmox
   order ice, geo, geothermal, sediments, climate, surface, shelf;
   `wire_resolve` (stops on a gap: required import without producer, or a
   capability mismatch); `coupler_prime` for every row; static rows
   (geothermal, sediments → ice) exchanged once.
5. Schedules: one `schedule_init(s, clk, name, par, group)` per distinct
   group named by an instance `out_*` key or `[kryos] restart_out` (`restart=.true.`).
6. Restore (§3.3) if `[kryos] restart /= "None"`, else cold start (§3.2).
7. `clock_start` (checks every restart time and `time_end` against every
   cadence; logs `clock_describe`); describes (§8); output at `time_init`.

### 3.2 Cold start (yelmox `domain_init_state`, order kept exactly)

```fortran
st0 = step_time_t(t0, t0, t_eval, dtt)                       ! per instance, t0 = time_init
do d = 1, size(run%dom)
  associate (dom => run%dom(d))
    ! a) bed and sea level (yelmox domain_init_isostasy)
    call exchange(dom, "ice", "geo", t0)
    call exchange_globe(run%globe, "geo", dom, "geo", t0)     ! bsl after bsl_update(time_rel0)
    call geo_init_state(dom%geo, st0, PASS_GEOMETRY)          ! isos_init_ref + isos_init_state
    call domain_check_reference(dom)                          ! geo vs domain refs, 1 m (yelmox check)
    call exchange(dom, "geo", "ice", t0)
    call ice_init_state(dom%ice, st0, PASS_GEOMETRY)          ! yelmo: calc_ytopo_diagnostic
    ! b) boundary fields from the hub (steps with init = .true.)
    call domain_hub_update(dom, t0)
    call exchange(dom, "hub", "climate", t0); call exchange_globe(run%globe, "forcing", dom, "climate", t0)
    call climate_init_state(dom%climate, st0, PASS_BOUNDARY)
    call exchange(dom, "hub", "surface", t0); call exchange(dom, "climate", "surface", t0)
    call surface_init_state(dom%surface, st0, PASS_BOUNDARY)  ! chion: time_equil cycles + one update
    call exchange(dom, "hub", "shelf", t0);   call exchange(dom, "climate", "shelf", t0)
    call shelf_init_state(dom%shelf, st0, PASS_BOUNDARY)
    call exchange_into(dom, "ice", t0)                        ! every row into ice (yelmox couple_to_yelmo)
    ! c) ice state: yelmo_init_state(robin-cold), then the recipe of §5 (none | equil)
    call ice_init_state(dom%ice, st0, PASS_STATE)
    ! d) the surface on the final initial geometry: one more chion cycle, no spin-up
    call domain_hub_update(dom, t0)
    call exchange(dom, "hub", "surface", t0)
    call surface_init_state(dom%surface, st0, PASS_FINAL)
    call exchange(dom, "surface", "ice", t0)
  end associate
end do
```

Recipe: yelmox `[sim]` in both base configs has `init_marine_H = False`,
`init_kill_shelves = False`, `init_time_thrm = 0`; `init_method = "none"`
(Greenland) and `"equil"` with `init_equil_time = 10` (Antarctica:
`yelmo_update_equil(time_tot=10, dt=1, topo_fixed=.false.)`). Phase 1 keeps only
this path, as a yelmo-adapter key `init_equil_time` of the ice instance group
(0 = none), run inside `PASS_STATE`; recon, marine-ice start, kill-shelves and
the thermal equilibration become experiment plugins in phase 2 (hook
`HOOK_ICE_STATE` before `PASS_STATE`, where yelmox sets the opt `cb_ref` guess
and NEGIS).

### 3.3 Restore (restore, never step)

Order: `clock_restart_read(clk, dir/clock_restart.nc)` (continuation if the
bundle time = `time_init`, else a logged branch), `cadence_restart_read` for
every instance from the same file; `manifest` check (§7); globe: forcing, geo
global (bsl, `A_ocean`); per domain: ice → geo → climate → surface → shelf →
geothermal → sediments, each `<comp>_restart_read(inst, dir, t)` then
`fields_restart_read` of its imports and exports (`<inst>_fields.nc`). Ice
before geo because FastIsostasy's `isos_init_state(use_restart)` takes the
restored Yelmo geometry (adapters-design). Then `domain_check_reference`,
`domain_hub_update` (the hub is derived, never stored). No exchange runs:
imports are restored as received.

Differences from yelmox, both bit-neutral at aligned restart times: the
climate is restored from `climate_grl_fields.nc` instead of re-derived with
`force=.true.`, and the derived coefficient `dTo` is restored from
`forcing_restart.nc` instead of recomputed by that climate update.

### 3.4 `run_step` (yelmox.f90:176-208 order, N domains)

```fortran
call clock_advance(run%clk)
! 1. globe, start of step: bsl_update(time_rel) (A_ocean), then coefs_update(time_rel, V_ice)
call geo_global_begin(run%globe%geo, st)                            ! bsl_update: A_ocean(bsl)
call gather_into(run%globe, run%dom, "forcing", t)                  ! V_ice of every domain, lag 1
call forcing_update(run%globe%forcing, st);  mark
! 2. hooks: none in phase 1 (yelmox step_relax, step_optimize: see question 1)
! 3. geo: every domain, then the global reduction (gather -> update -> scatter, decision 5)
do d
  if (cadence_due(dom%geo%cadence, clk)) then
    call exchange(dom, "ice", "geo", t)                              ! H_ice (+dzbdt_corr) at t-dtt
    call exchange_globe(globe, "geo", dom, "geo", t)                 ! bsl at the start of the step
    call geo_update(dom%geo, st);  mark
  end if
end do
call geo_global_reduce(run%globe%geo, run%dom, t)                    ! bsl -= sum_d dV_d / A_ocean
! 4. per domain: ice, hub, climate, surface, shelf (yelmox couple_to_yelmo .. step_shelf)
do d
  call exchange_into(dom, "ice", t)            ! geo (lag 0), surface, shelf, climate Qd (lag 1)
  call ice_update(dom%ice, st);  mark
  call domain_hub_update(dom, t)               ! couple_yelmo_to_htopo
  if (cadence_due(dom%climate%cadence, clk)) then
    call exchange(dom, "hub", "climate", t); call exchange_globe(globe, "forcing", dom, "climate", t)
    call climate_update(dom%climate, st);  mark
  end if
  if (cadence_due(dom%surface%cadence, clk)) then   ! every step
    call exchange(dom, "hub", "surface", t); call exchange(dom, "climate", "surface", t)
    call surface_update(dom%surface, st);  mark
  end if
  if (cadence_due(dom%shelf%cadence, clk)) then     ! every step
    call exchange(dom, "hub", "shelf", t); call exchange(dom, "climate", "shelf", t)
    call shelf_update(dom%shelf, st);  mark
  end if
end do
! 5. end of step
if (forcing_kill(run%globe%forcing)) call clock_finish(run%clk)
call output_write_due(run%out, ...)             ! every due schedule evaluated once, then streams
if (schedule_due(run%sched(run%i_rst), run%clk)) call restart_write_bundle(...)
```

`mark` is `cadence_mark(inst%cadence, clk)` after the update; `st` comes from
`cadence_interval` and the instance's `clock`. Presence is tested with
`allocated(dom%climate)`, never a backend name; no experiment conditionals.
`exchange(dom, src, dst, t)` runs the domain's resolved rows `src → dst` in
table order; `"hub"` is a source. `exchange_into(dom, dst, t)` runs every
row into `dst`, in the order geo, surface, shelf, climate (yelmox
`couple_to_yelmo`); rows from static producers (cadence off) run once, at init. The global reduction between the two loops is what makes
`z_sl` at the ice see this step's `bsl`, as FastIsostasy's in-call update does
today (§10 R5).

### 3.5 `run_end`

`<comp>_end` for every instance (Yelmo last, after the last access),
`timer_print_totals`, final log line. No output (already written).

## 4. Timing, cadences, lags

| Instance | Cadence (phase 1 configs) | `clock` | yelmox source |
|---|---|---|---|
| `ice_*` | every step (`dt = 0`), dtt = 10 | model | `yelmo_update(ts%time)` |
| `geo_*` | `dt = 10` | model | `dt_isos`, `isos_update(ts%time)` |
| `geo_glob` (bsl) | every step | forcing | `bsl_update(time_rel)` in the driver |
| `climate_*` | `dt = 10` | forcing | `dt_clim`, snapshots at `time_rel` |
| `surface_*` | every step | forcing (orbit) | chion, one annual cycle per call |
| `shelf_*` | every step | model | marine_shelf |
| `forcing` | every step | forcing | `coefs_update(time_rel)` |
| `geothermal_*`, `sediments_*` | off (`dt < 0`), static | — | read once at init |

Cadences are on the absolute grid (`fesm_clock`), yelmox's on
`time_elapsed`: identical when `time_init` is a multiple of every `dt`, as in
both configs (`time_init = 0`). With `method = "const"` and `time_const = 0`
both configs freeze the forcing time at 0.

Lags (contracts-design owns the rows; the core sequence produces them):

| Row | Lag | Value the consumer sees in step n (t_{n-1} → t_n) |
|---|---|---|
| geo → ice `z_bed`, `z_sl` | 0 | geo integrated to t_n with load H(t_{n-1}) |
| ice → geo `H_ice` | 1 | H(t_{n-1}) |
| surface → ice `smb`, `T_srf` | 1 | surface at t_{n-1} on H(t_{n-1}), forcing at time_rel(t_{n-1}) |
| shelf → ice `bmb_shlf`, `T_shlf` | 1 | shelf at t_{n-1} |
| climate → ice `Qd` | 1 | (no Qd from snapshots) |
| hub → climate, surface, shelf; climate → surface, shelf | 0 | after the ice step |
| forcing → climate (coefficients) | 0 | `coefs_update(time_rel(t_n))` at the start of step n |
| ice → forcing `V_ice` | 1 | ice at t_{n-1} |
| globe bsl → geo | 0 | `bsl` of the last reduction, `A_ocean` from `bsl_update(time_rel(t_n))` |

Lag check in `wire_exchange`: the published `t1` must equal the producer's
`cadence_t_last` (hub: `hub_t`) and be ≤ t_n (lag 0) or ≤ t_{n-1} (lag 1);
otherwise stop naming the row (catches a producer that did not `put` and a
sequence edit that swapped the order).

## 5. Configuration

Par file: `nml`; every group is validated with `nml_validate` against its
template in `input/kryos_defaults.nml` (core keys) and, for an instance group,
against `input/<backend>_defaults.nml` group `[instance]` (backend keys) as
well; a key is valid if either declares it (question 3). Named library groups
(`snapshot_*`, `coef_*`, `term_*`, data field groups) are found with
`nml_find_group([par_file, library(:)], group)`; instance, domain and run groups
must be in the par file. Backend library groups (`[yelmo]`, `[isos]`,
`[marine_shelf]`, ...) keep their own schemas; their `restart` keys must be
`"None"` (restarts come from the bundle), checked by the adapters.

Templates in `input/kryos_defaults.nml` (one key per line in the file):

| Template | Keys (default) |
|---|---|
| `[kryos]` | `domains` (""), `globe` ("globe"), `clock` ("clock"), `restart` ("None"), `restart_out` ("none"), `out_dir` ("./"), `map_fldr` ("maps"), `phys_const` ("Earth"), `phys_const_file` ("input/phys_const.nml"), `library` (""), `write_lossless` (False), `log_level` ("info"), `timing` (True) |
| `[clock]` | `method` ("rel"), `time_ref` (1950), `calendar` ("365_day"), `year` ("tropical"); `time_init`, `time_end`, `dtt` required |
| `[schedule]` | `method` ("none"), `dt`, `times`, `file` (fesm_clock keys) |
| `[globe]` | `forcing` (""), `geo` ("") |
| `[domain]` | `grid_hub`, `topo` ("topo"), `regions` ("regions"), `basins`, `out_2D` ("none"); one key per component (`ice` … `sediments`, "" = absent) |
| `[instance]` | `backend`, `grid` ("" = chain), `dt` (0), `phase` (0), `clock` ("model"), `out_2D`, `out_1D` ("none"), `out_vars` ("" = all exports) |

`dt = 0`: every step; `dt > 0`: that cadence grid; `dt < 0`: never after
init. `grid = ""` follows `GRID_DEFAULT` (yelmox `kryos.f90:574-578`): ice ←
hub, geo ← ice, climate ← ice, surface ← climate, shelf ← hub, geothermal and
sediments ← ice; grids are `maps/grid_<name>.txt`, resolved and cached by the
domain's coupler. A domain is named by its group, its instances by the
component keys (one per component per domain), so instance names are unique
per run.

Minimal Greenland par file. Backend library groups (`[topo]`, `[regions]`,
`[coef_alpha]`, `[barysealevel]`, `[isos]`, `[yelmo]` … `[yelmo_data]`,
`[snapshots]`, `[snapshot_pd]`, `[orbit]`, `[surface_chion]`, `[chion]`,
`[itm_chion]`, `[chion_const]`, `[marine_shelf]`) are unchanged from
`yelmox_Greenland.nml` and not repeated:

```fortran
&kryos
    domains     = "Greenland"
    restart_out = "out_rst"
    library     = "input/clim_db.nml" "input/forcing_db.nml"
/
&clock
    method     = "const"
    time_const = 0.0
    time_init  = 0.0
    time_end   = 15e3
    dtt        = 10.0
/
&out_2D
    method = "const"
    dt     = 200.0
/
&out_2Dsm
    method = "const"
    dt     = 50.0
/
&out_1D
    method = "const"
    dt     = 10.0      ! yelmox tm_1D dt = 1: same times at dtt = 10; fesm_clock rejects dt < dtt
/
&out_rst
    method = "const"
    dt     = 20e3      ! + always time_end
/
&globe
    forcing = "forcing"
    geo     = "geo_glob"
/
&forcing
    backend = "tsgen"
    clock   = "forcing"
    coefs   = "alpha" "dTo"
    out_1D  = "out_1D"
/
&geo_glob
    backend = "fastisostasy"
    clock   = "forcing"
    bsl     = "barysealevel"
    out_1D  = "out_1D"
/
&Greenland
    grid_hub   = "GRL-16KM"
    basins     = "Zwally2012.group"
    out_2D     = "out_2D"
    ice        = "ice_grl"
    geo        = "geo_grl"
    climate    = "climate_grl"
    surface    = "surface_grl"
    shelf      = "shelf_grl"
    geothermal = "ghf_grl"
    sediments  = "sed_grl"
/
&ice_grl
    backend         = "yelmo"
    yelmo           = "yelmo"
    init_equil_time = 0.0     ! yelmox init_method = "none"; Antarctica: 10.0 ("equil")
    out_2D          = "out_2D"
    out_1D          = "out_1D"
    out_sm          = "out_2Dsm"
/
&geo_grl
    backend = "fastisostasy"
    dt      = 10.0            ! yelmox dt_isos
    isos    = "isos"
    out_2D  = "out_2D"
    out_1D  = "out_1D"
/
&climate_grl
    backend   = "snapshots"
    clock     = "forcing"
    dt        = 10.0          ! yelmox dt_clim
    snapshots = "snapshots"
    reference = "pd"
    terms     = "lgm" "ocean"
    out_2D    = "out_2D"
/
&surface_grl
    backend       = "chion"
    clock         = "forcing"
    surface_chion = "surface_chion"
    chion         = "chion"
    orbit         = "orbit"
    out_2D        = "out_2D"
/
&shelf_grl
    backend      = "marine_shelf"
    marine_shelf = "marine_shelf"
    out_2D       = "out_2D"
/
&ghf_grl
    backend = "data"
    dt      = -1.0
    Q_geo   = "ghf_const_50mW"   ! field group (adapters-design): 0.05 W m-2
/
&sed_grl
    backend = "data"
    dt      = -1.0
    H_sed   = "sed_none"         ! field group: 0 m
/
```

All grids default to the hub. Antarctica: `grid_hub = "ANT-16KM"`,
`[ice_ant] grid = "ANT-32KM"` (geo, climate, surface follow the chain to
`ANT-32KM`; shelf to the hub `ANT-16KM`), `init_equil_time = 10.0`,
`[forcing] coefs = ""`, `[climate_ant] terms = ""`, Q_geo from the Shapiro file.

Mapping from `yelmox_Greenland.nml`:

| yelmox | kryos |
|---|---|
| `ctrl.run_step` | `[kryos] clock` (names the clock group) |
| `ctrl.tstep_method`, `tstep_const`, `calendar`, `calendar_ref` | `[clock] method`, `time_const` (`time_const_cal` if calendar), `time_ref`; `calendar` now CF (`365_day`) |
| `ctrl.time_init`, `time_end`, `dtt`, `restart`, `write_lossless` | `[clock]` same keys; `[kryos] restart`, `write_lossless` |
| `comps.with_*` | dropped: a component runs iff the domain names an instance |
| `comps.climate`, `surface_method`, `dt_clim`, `dt_isos` | `[climate_grl] backend`, `[surface_grl] backend`, `[climate_grl] dt`, `[geo_grl] dt` |
| `sim.init_method = none/equil`, `init_equil_time` | `[ice_*] init_equil_time` (0 = none) |
| `sim.init_marine_H`, `init_kill_shelves`, `init_time_thrm`, `recon_*`, `use_negis`, `scale_glacial_smb`, `lim_pd_ice` | dropped (experiments, phase 2); all off in both configs |
| `sim.relax`, `sim.opt`, `[relax]`, `[opt]`, `[negis]` | question 1 |
| `output.write_*` | `out_2D`/`out_1D` = `"none"` per instance; `write_htopo` → `[Greenland] out_2D` |
| `output.write_cmip`, `dt_cmip` | dropped (False in both; CMIP becomes post-processing) |
| `domain.name`, `grid_hub`, `basins`; `grid_ice`… `grid_shelf` | group name, `[Greenland] grid_hub`, `basins`; `[<inst>] grid` |
| `domain.phys_const` | `[kryos] phys_const` (one per run) |
| `[topo]`, `[regions]` | unchanged, named by `[Greenland] topo`, `regions` |
| `tm_1D`, `tm_2D`, `tm_2Dsm`, `tm_rst` | `[out_1D]` (dt = 10), `[out_2D]`, `[out_2Dsm]`, `[out_rst]` |
| `forcing.coefs`; `reference`, `terms`; `var_defs` | `[forcing] coefs`; `[climate_grl] reference`, `terms`; `[kryos] library` |
| `[barysealevel]`, `[isos]` | unchanged, named by `[geo_glob] bsl`, `[geo_grl] isos` |
| `[sed] method = 0`, `[ghf] use_obs = False, ghf_const = 50` | data instances (H_sed = 0; Q_geo = 0.05 W m-2) |
| `[smbpal]`, `[smb_simple]`, `[itm]` | unused (not read, ignored) |

## 6. Output

A stream is (producer, file, schedule); producers are the instances, the
hub, the run (timing). Per instance `out_2D` → `<instance>.nc`, `out_1D` →
`<instance>_ts.nc`: the default `write_output` writes the exports in
`out_vars` with the `variable_io` table `input/kryos-variables-<component>.md`
(fixed `prec_abs`/`prec_rel`, `float` as yelmox). A backend may declare
extra streams (yelmo: `yelmo`, `yelmo_sm`, `yelmo_ts`), keys `out_<stream>`, files
`<instance>_<stream>.nc`. Instance names are unique in a run (`run_init`
stops on a duplicate), so file names need no domain or grid prefix. Each output time: `schedule_due` once per schedule
(it moves a cursor), then `<comp>_write_output(inst, file, t, stream, init)`
per due stream; files are created on a run's first write (anew per chunk).

| yelmox file | kryos phase 1 |
|---|---|
| `yelmo.nc` (`YELMO_VARS_2D`) | `ice_grl.nc` (exports) + `ice_grl_yelmo.nc` (backend set) |
| `yelmo_sm.nc`, `yelmo_ts.nc` + regions | `ice_grl_yelmo_sm.nc`, `ice_grl_yelmo_ts.nc`, `ice_grl_yelmo_ts_<region>.nc` (by the adapter) |
| `coef_*` in `yelmo_ts.nc` | `forcing_ts.nc` |
| `htopo.nc` | `Greenland_hub.nc` (`htopo_write_init/step`) |
| `isos.nc`, `isos_ts.nc` (`bsl`) | `geo_grl.nc` (+ backend `w`, `we`), `geo_grl_ts.nc`; `geo_glob_ts.nc` (`bsl`, `A_ocean`) |
| `chion.nc`, `snap.nc`, `mshlf.nc` | `surface_grl.nc`, `climate_grl.nc`, `shelf_grl.nc` |
| `timing.txt` | `timing.txt` (on `out_1D`) |
| `yelmo_cmip.nc`, `yelmo_ts_cmip.nc` | deferred: post-processing over ice exports (design §9) |

The regression tool compares by a variable map (yelmox file/name → kryos
file/name), with `write_lossless = True` on both sides.

## 7. Restart

```
restart-15000/                     <time> = clock time in years, shortest exact decimal
  manifest.nml                     [bundle] + one group per instance
  clock_restart.nc                 clock_restart_write + cadence_restart_write of every instance
  forcing_restart.nc               tsgen states and derived values ; forcing_fields.nc
  geo_glob_restart.nc              bsl_now, A_ocean_now ; geo_glob_fields.nc
  Greenland/
    ice_grl_restart.nc  (+ Yelmo sidecars)   ice_grl_fields.nc
    geo_grl_restart.nc  geo_grl_fields.nc    surface_grl_restart.nc  surface_grl_fields.nc
    shelf_grl_restart.nc  shelf_grl_fields.nc   climate_grl_fields.nc   ghf_grl_fields.nc  sed_grl_fields.nc
```

`manifest.nml` (read back with `nml_read`): group `[bundle]` (`time`,
`ticks`, `ticks_per_year`, `kryos_git`, `domains`) and one group per instance,
named as in the par file (`scope`, `component`, `backend`, `grid`, `files`,
`git`); hashes from a generated `kryos_build_info.f90`.

- Written when the restart schedule is due (always at `time_end`); write
  order = restore order (§3.3). `<instance>_restart.nc` only for stateful
  instances (no file is a valid choice of the adapter); `_fields.nc` for all.
  The hub is not stored. Folder creation as yelmox (`mkdir -p`).
- Validity: `clock_start` stops at init if a restart time or `time_end` is
  off any cadence grid; a killed run writes a bundle only if
  `clock_restart_allowed`, else logs why (fesm_clock §4.2).
- Restore checks the manifest against the par file: the same instances with
  the same component, backend and grid, or stop naming the instance (a
  changed `dt` is allowed, logged by fesm_clock).
- Chunks: chunk k sets `[clock] time_init` = previous `time_end` and
  `[kryos] restart = "restart-<time_end>"`; names are known at submission.
  Continuation keeps `n_origin`, so `time_elapsed` runs on across chunks
  (yelmox restarts it per chunk: relax/opt windows differ, documented).
- Deferred: forcing-increment bundles (`restart_every_df`, 0 in both
  configs) to phase 3; the kill switch stays.

## 8. Errors, logging, describe

- Core errors: `fesm_error("kryos_<module>::<proc>", msg, detail)`, no
  `stat` (all are configuration or wiring errors); messages name group, key,
  instance, row or file and the allowed values ("row `surface.smb -> ice.smb`:
  published at 20 yr, expected 10 yr (lag 1)").
- Prefix: wrappers set `log_set_prefix(<instance>)` around each backend call
  and restore it (process-global; contracts-design); core lines use `who =
  "kryos_run"` etc.; `write(*,*)` in backends stays unprefixed. One INFO line per step.
- Describe at init, in this order: run (par file, out dir, restart, N
  domains), `phys_const_log`, per domain grids (instance → grid, "(default:
  hub)") and primed maps, per instance backend + caps + cadence +
  `fields_describe` of imports and exports, `wiring_describe` (rows, method,
  lag, transform), the manifest on restore, `clock_describe` (from `clock_start`).

## 9. Physical constants

One `phys_const_class` per run: `phys_const_load(run%cnst, phys_const_file,
group=phys_const, defaults_file=<fesm-utils>/par/phys_const_earth.nml,
defaults_group="Earth")`, then `phys_const_require`. `input/phys_const.nml` is
yelmox's `yelmo_phys_const.nml` (incl. `sec_year = 31556926.0`); it goes to
every instance in `init_ctx_t` and to `htopo_init`. An adapter whose library
keeps its own record calls `phys_const_compare(cnst, own, "kryos", backend)`
and logs differences (none in phase 1). The clock's `year = "tropical"`
(31556926.08 s) is logged beside it; unifying them is a phase 2 item.

## 10. Requirements on the other two notes for bit-identity

- R1: the exchange remaps `real(sp)` values widened to `dp` (the `dp`
  specifics of `remap`); consumers narrow. Equal to yelmox's `sp` remap
  because `map_field_grid_grid_sp` widens, remaps in `dp` and narrows; a unit
  test proves it for every method used (missing value and mask arguments included).
- R2: no `dp` arithmetic between a remap and the narrowing where yelmox did
  `sp` arithmetic. Holds in both configs for geo → ice (identity grids) and for
  the Greenland hub (identity: the hub mirrors Yelmo incl. fractional `f_grnd`,
  not `htopo_update`). Not guaranteed for the Antarctica hub anomalies
  (ice 32 km → hub 16 km, yelmox subtracts in `sp` before the bilinear map) and
  for `bmb_shlf` (m ie/yr ↔ kg m-2 yr-1 around a `con` map shelf 16 km → ice 32 km): question 2.
- R3: wiring rows carry yelmox's explicit methods (`bilin` hub → shelf, `con`
  surface/shelf → ice, `bilin` climate → surface/shelf); `auto` only where
  yelmox used `remap_method_smooth` (same choice for these grids).
- R4: per-grid references are the hub's `sp` reference through the `sp`
  remap path, stored widened (Yelmo's `topo_init`/`topo_pd` bit-identical).
- R5: FastIsostasy updates `bsl` inside `isos_update`, in `sp`, before
  `calc_z_ss`. The regional adapter passes the globe's `bsl_now`/`A_ocean_now`
  into its own `bsl_class`; the global reduction computes `bsl_now -
  sum(dV)/A_ocean_now` in FastIsostasy's precision: identical for N = 1.
- R6: smb in kg m-2 yr-1 equals chion's mm w.e./yr; the yelmo adapter applies
  `smb*conv_we_ie*1e-3` in `sp`, in that order.
- R7: passes and order of §3.2 (the `PASS_FINAL` chion cycle changes chion's state).
- R8: a cold-start publication is at `[t0, t0]`; exports read at step 1 with
  lag 1 (`smb`, `T_srf`, `bmb_shlf`, `T_shlf`) must accept it (INSTANT, or a
  stated convention for MEAN).

## 11. Build and layout

- `.configme/manifest.toml`: `fesm-utils:phase0-kryos` until it merges to `dev`.
- Delete `src/core/kryos_fesm_pending.f90`; `kryos_component` uses `fesm_fields`, `fesm_clock`.
- New `src/core/`: `kryos_config`, `kryos_wiring`, `kryos_domain`, `kryos_globe`,
  `kryos_output`, `kryos_restart`, `kryos_run` + `kryos_run_sequence.f90`.
  Components `geo`, `climate`, `surface`, `shelf`, `geothermal`, `sediments`,
  `forcing` with factories; backends per adapters-design (`forcing/tsgen` has no
  stub: fesm-utils). Where yelmox `libs/` code (snapshots, climate_terms,
  marine_shelf, surface_chion) is compiled: adapters-design.
- Generated `kryos_build_info.f90` (git hashes of kryos and enabled backend roots).
- `input/`: `kryos_defaults.nml`, `phys_const.nml`, `kryos-variables-<component>.md`,
  backend defaults; `par/kryos_{Greenland,Antarctica}.nml`; `.runme/info.json`
  (links `input`, `maps`, `fesmdata`, `ice_data`, `isostasy_data`); `tools/`
  (Julia) yelmox ↔ kryos variable map and diff.

## 12. Implementation plan

Core commits (C*), with the contracts (K*) and adapters (A*) commits they
wait for. Each builds with `yelmo=0` and `yelmo=1` and passes `make check`.

| # | Commit | Verified by |
|---|---|---|
| C1 | fesm-utils pin to `phase0-kryos`; drop `kryos_fesm_pending`; `component_t` on `fesm_fields`/`fesm_clock` | build; existing driver runs |
| C2 | `kryos_config` + `input/kryos_defaults.nml` (+ multi-schema validation, question 3) | `test_config`: defaults, library lookup, unknown key and missing group abort |
| C3 | `kryos_run` clock, schedules, zero-domain loop, clock-only restart | `test_run_clock`: output and restart times = yelmox `timeout` times for both configs; continuation and branch |
| C4 | `test_remap_sp_dp` (R1) for `con`, `bilin`, `auto`, 2D and 3D, GRL-16KM, ANT-16KM ↔ ANT-32KM | bitwise equality |
| C5 | `kryos_wiring` engine: rows on two field states, transforms hook, lag check, describe | `test_wiring`: toy states, identity and remapped rows, every lag error aborts |
| C6 | `kryos_domain` without instances: grids, chain, hub, coupler, refs, hub output | hub and refs at t0 bit-identical to yelmox `htopo.nc` and Yelmo `topo_init` (both configs) |
| C7 | `kryos_globe` (after K: forcing, geo global; A: tsgen, fastisostasy bsl) | `forcing_ts.nc` `coef_alpha` vs yelmox `yelmo_ts.nc` (Greenland, `method = rel`, 21 kyr) |
| C8 | `kryos_output` streams and default writer | file layout test; variable tables load |
| C9 | sequence: cold start + `run_step` (after K and A for all phase 1 instances) | Greenland 1 step, then 1 kyr, lossless, bit-identical to yelmox |
| C10 | `kryos_restart`: bundle, manifest, restore | 0–1000 = 0–500 + 500–1000 bitwise; manifest mismatch aborts; Greenland 15 kyr on Levante |
| C11 | Antarctica: `init_equil_time`, non-identity grids (question 1, 2) | `yelmox_Antarctica` gate |
| C12 | driver arg handling, `par/` files, `.runme/info.json`, docs page | runme run of both par files |

## 13. Open questions

1. `yelmox_Antarctica` runs with `sim.relax = True` and `sim.opt = True`
   (`L21` cb_ref, `L21-points` tf_corr on the shelf), which design §14 puts in
   phase 2 as experiment plugins. Gate Antarctica in phase 1 against yelmox run
   with `sim.relax=False sim.opt=False` (plugins and their gate in phase 2,
   proposed), or bring the hooks and both plugins into phase 1?
2. With the `dp` exchange, the Antarctica gate cannot be guaranteed
   bit-identical where yelmox does `sp` arithmetic around a non-identity remap
   (hub anomalies ice → hub; `bmb_shlf` unit round trip shelf → ice). Accept
   documented ulp-level differences there (proposed), or require bit-identity
   by constraining the adapters (Yelmo exports its anomalies computed in `sp`;
   `bmb_shlf` crosses in m ie/yr, an exception to decision 3)?
3. Instance groups mix core keys and backend keys. Add a multi-schema form of
   `nml_validate` to fesm-utils (`phase0-kryos`; proposed), or repeat the core
   keys in every backend's `[instance]` template with a check in `make check`?
