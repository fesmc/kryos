# Component contracts, wiring and hub — design note (phase 1)

Status: proposal, 2026-10-24, for review. No code. Scope: the per-component
vocabularies, `component_t` and the per-component abstract types with their
wrappers and factories, the fixed wiring tables with lags and transforms, the
hub formulas, geo's two levels, the `forcing` component and masks, for phase 1
of `design.md` §14 (one domain; yelmo, fastisostasy, snapshots, chion,
marine_shelf; geothermal and sediments as data; a minimal globe with forcing and
geo's global part).

Referenced, not repeated: `core-design.md` (run, domain, globe, sequence,
configuration, output, restart; its names are adopted: `init_ctx_t`,
`step_time_t`, `grid_ref_t`, `PASS_*`, `exchange*`, `gather_into`,
`domain_hub_update`, `domain_check_reference`, the `sp` hub of phase 1),
`regression-design.md` (gate), and `adapters-design.md` (not yet written;
assumed to implement the "yelmox" columns of §2 as the contract). To fold into
core-design:

- **Slots** (decision 1 style rule): `domain_t`/`globe_t` hold `type(ice_slot_t)
  :: ice`, not `class(ice_t), allocatable :: ice`; the core never declares
  `class(...)`. Presence: `ice_present(dom%ice)`.
- **geo reduction**: wrappers cannot take `domain_t` (core §1), so
  `geo_reduce(globe%geo, run%dom, t)` becomes `gather_into` + `geo_reduce(globe%geo,
  st)` + `geo_scatter(dom%geo, st)` (§6).
- `V_ice` is `M_ice` [Gt] (a mass: `V_ice*rho_ice*1e-3`); `step_time_t` gains
  `dt` (§3.1); the hub has `hub_imports` and `hub_exports` gains `f_ice` (§5).


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

## 1. Conventions

- Units, signs, fractions: `design.md` §5. Names are fixed here; a backend
  never invents an exchange name.
- **tsem.** States are INSTANT; rates over a producer step are MEAN. Lagged
  rate imports are INSTANT (`MEAN → INSTANT`, fields-design §3.5); phase 1 has
  equal cadences on every such row (checked, §4.4); MEAN imports for slower
  consumers are phase 2. Cold start (core R8): MEAN exports put in
  `init_state` cover the virtual previous step `[t0 − dt, t0]` (the wrapper
  sets it), so lag-1 consumers find `t1 = t0` at step 1.
- **cons** is the field's nature; the remap method is the row's (§4.2), as in
  yelmox; `cons` decides only for method `""`.
- **Optional fields** exist only if the backend sets their flag in `init`
  (`self%opt%T_shlf`, `self%caps%has_ocn`); the component layer declares what
  the flags say; rows touching an undeclared field are inactive (§4.4).
- **Masked:** no phase 1 field carries `valid` (all producers fill everywhere;
  masked rows across grids wait for coupler pass-through), so the column is
  omitted.
- **Axes.** Axis 3 is `month` (12) or `depth` (the producer's ocean levels, m,
  positive down), declared by the producer as `axis_t` and passed to consumers
  in `ctx%axes` (climate is initialised before shelf).
- **Precision.** Each instance declares `prec` (its arrays' kind; `sp` for every
  phase 1 backend and the phase 1 hub). Exports are formed in `prec` and widened
  exactly; imports are narrowed once; transforms: §9.

## 2. Vocabulary

Columns: rank (zdim), tsem (I/M), cons (int/ext), req (R required, O
optional, C capability), yelmox counterpart and the adapter's conversion.

### 2.1 ice (regional; yelmo)

| Imports | units | rank | tsem | cons | req | yelmox (Yelmo `bnd%`) and conversion |
|---|---|---|---|---|---|---|
| `z_bed` | m | 2 | I | int | R | `z_bed`; narrow |
| `z_sl` | m | 2 | I | int | R | `z_sl`; narrow |
| `smb` | kg m-2 yr-1 | 2 | I | ext | R | `smb` [m ie/yr] = `smb_sp*conv_we_ie*1e-3`, in `sp`, that order (core R6) |
| `T_srf` | K | 2 | I | int | R | `T_srf` |
| `bmb_shlf` | kg m-2 yr-1 | 2 | I | ext | R | `bmb_shlf` [m ie/yr] = `/rho_ice` (domain `[Earth]`) |
| `T_shlf` | K | 2 | I | int | O | `T_shlf` |
| `tf_shlf` | K | 2 | I | int | O | `tf_shlf`; no row in phase 1 (yelmox does not pass it) |
| `Q_geo` | W m-2 | 2 | I | int | O | `Q_geo` [mW m-2] = `*1e3` |
| `H_sed` | m | 2 | I | int | O | `H_sed` |
| `Qd` | m3 s-1 | 2 | I | ext | O | `Qd` (per cell, as esm; no producer in phase 1) |

| Exports | units | rank | tsem | cons | req | yelmox (Yelmo) and conversion |
|---|---|---|---|---|---|---|
| `H_ice` | m | 2 | I | ext | R | `tpo%now%H_ice`; widen |
| `z_srf` | m | 2 | I | int | R | `tpo%now%z_srf` |
| `z_bed` | m | 2 | I | int | R | `bnd%z_bed` (the bed used; CISM un-folds, decision 11) |
| `z_sl` | m | 2 | I | int | R | `bnd%z_sl` (the sea level used) |
| `f_ice`, `f_grnd` | 1 | 2 | I | int | R | `tpo%now%f_ice`, `f_grnd` |
| `dHidt` | m yr-1 | 2 | M | ext | R | `tpo%now%dHidt` |
| `smb`, `bmb`, `cmb`, `mb_net` | kg m-2 yr-1 | 2 | M | ext | R | `tpo%now%smb`, ... [m ie/yr] `*rho_ice` (regression §3.2) |
| `M_ice` | Gt | 0 | I | ext | R | `reg%V_ice*bnd%c%rho_ice*1e-3`, in `sp` (yelmox.f90:185) |
| `ux_s`, `uy_s` | m yr-1 | 2 | I | int | O | `dyn%now%ux_s`, `uy_s` at cell centres |
| `dz_bed_dt_corr` | m yr-1 | 2 | I | int | O | `bnd%dzbdt_corr` |
| `z_bed_ref` | m | 2 | I | int | O | `bnd%z_bed_ref` (reference check, §5.3) |

### 2.2 geo (regional part; fastisostasy) and its global part

| Regional imports | units | rank | tsem | cons | req | source |
|---|---|---|---|---|---|---|
| `H_ice` | m | 2 | I | ext | R | ice (total, decision 4) |
| `dz_bed_dt_corr` | m yr-1 | 2 | I | int | O | ice; declared by fastisostasy only |
| `bsl` | m | 0 | I | int | R | global part (twice per step, §6) |
| `A_ocn` | m2 | 0 | I | ext | O | global part; fastisostasy only |

| Regional exports | units | rank | tsem | cons | req | yelmox (FastIsostasy) and conversion |
|---|---|---|---|---|---|---|
| `dz_bed` | m | 2 | I | int | R | `out%w + out%we`, sum in `sp` |
| `dz_ss` | m | 2 | I | int | R | `out%z_ss - now%bsl`, in `sp` |
| `bsl` | m | 0 | I | int | R | the global `bsl` this domain last solved with (set by `geo_scatter`, §6) |
| `dV_bsl` | m3 | 0 | I | ext | O | `now%deltaV_bsl` of this update; 0-D gather source |
| `z_bed_ref` | m | 2 | I | int | R | `out%z_bed - out%w - out%we` (reference check) |
| `dz_bed_dt` | m yr-1 | 2 | I | int | O | `now%dwdt`; no consumer in phase 1 |

Global part (0-D, INSTANT): exports `bsl` [m] (`bsl_class%bsl_now`), `A_ocn`
[m2] (`A_ocean_now`); import `dV_bsl` [m3], the gathered sum (§6).

### 2.3 climate (regional; snapshots)

| Imports | units | rank | tsem | req | source |
|---|---|---|---|---|---|
| `z_srf`, `H_ice`, `z_bed`, `f_grnd`, `z_sl` | m, m, m, 1, m | 2 | I | R | hub |
| `z_srf_ref` | m | 2 | I | R | hub (static reference, published once) |
| `coef_<name>` | coef units | 0 | I | C | forcing export `<name>`, one per coefficient its terms name (§7) |

| Exports | units | rank | tsem | cons | req | yelmox (`climate_out_class%now`) and conversion |
|---|---|---|---|---|---|---|
| `tas` | K | 3 (month) | I | int | R | `tas` |
| `pr` | kg m-2 yr-1 | 3 (month) | I | ext | R | `pr` [mm d-1] `*360` (regression §3.2, adapters) |
| `z_srf_clim` | m | 2 | I | int | R | `z_srf` (the surface `tas`, `pr` refer to) |
| `lapse` | K m-1 | 3 (month) | I | int | R | `lapse(12)`, spatially uniform for snapshots |
| `beta_p` | K-1 | 2 | I | int | R | `beta_p` |
| `to_ann`, `so_ann` | K, psu | 3 (depth) | I | int | C `has_ocn` | `to_ann`, `so_ann` |
| `dto_ann` | K | 3 (depth) | I | int | C `has_ocn` | `now%to_ann - ref%to_ann`, in `sp` (coupling:612) |
| `to_shlf`, `so_shlf`, `dto_shlf`, `dso_shlf` | K, psu, K, psu | 2 | I | int | C `has_ocn_shelf` | `T_shlf`, `S_shlf`, `dT_shlf`, `dS_shlf` (esm, phase 2) |
| `smb`, `T_srf` | kg m-2 yr-1, K | 2 | M | ext, int | C `has_smb` | `smb` [mm/a], `tsrf` (esm/rembo, phase 2) |
| `Qd` | m3 s-1 | 2 | I | ext | C `has_Qd` | `Qd` (esm, phase 2) |
| `<derived>` | K | 0 | I | int | C | `factor*mean(tsl_ann - ref%tsl_ann, region)` per derived coefficient `from = air_anomaly` (yelmox_climate.f90:338) |

`tsl_ann` (smb_simple), `ta_ann`/`ref%ta_ann` (glacial smb) and `pr_ann` are
not exported in phase 1 (§10). Capability flags live in `caps_t` (§3.1);
`has_ocn` is new (yelmox implied it by `.not. has_ocn_shelf`).

### 2.4 surface (regional; chion)

| Imports | units | rank | tsem | req | source |
|---|---|---|---|---|---|
| `z_srf`, `H_ice`, `z_bed`, `z_sl` | m | 2 | I | R | hub |
| `tas`, `pr` | K, kg m-2 yr-1 | 3 (month) | I | R | climate, lapse transforms (§4.1) |

| Exports | units | rank | tsem | cons | req | yelmox (`surface_chion%ann`) |
|---|---|---|---|---|---|---|
| `smb` | kg m-2 yr-1 | 2 | M | ext | R | `smb` [mm w.e./yr], factor 1 |
| `T_srf` | K | 2 | M | int | R | `tsrf` |
| `melt`, `runoff`, `refrz` | kg m-2 yr-1 | 2 | M | ext | O | same names, factor 1 |
| `albedo` | 1 | 2 | M | int | O | `alb_s` |

### 2.5 shelf (regional; marine_shelf)

| Imports | units | rank | tsem | req | source |
|---|---|---|---|---|---|
| `H_ice`, `z_bed`, `f_grnd`, `z_sl`, `z_srf` | m, m, 1, m, m | 2 | I | R | hub |
| `to_ann`, `so_ann`, `dto_ann` | K, psu, K | 3 (depth) | I | set A | climate (`has_ocn`) |
| `to_shlf`, `so_shlf`, `dto_shlf`, `dso_shlf` | K, psu, K, psu | 2 | I | set B | climate/ocean (`has_ocn_shelf`) |

Contract (checked at init): exactly one of set A, set B is complete and
active, else stop "shelf `shelf_grl` needs ocean forcing: climate `clim_grl`
declares neither profiles (has_ocn) nor shelf-base values (has_ocn_shelf)".
The shelf owns the profile → shelf-base interpolation (decision 6).

| Exports | units | rank | tsem | cons | req | yelmox (`mshlf%now`) |
|---|---|---|---|---|---|---|
| `bmb_shlf` | kg m-2 yr-1 | 2 | M | ext | R | `bmb_shlf` [m ie/yr] `*rho_ice` |
| `T_shlf` | K | 2 | I | int | R | `T_shlf` |
| `tf_shlf` | K | 2 | I | int | R | `tf_shlf` |

`mask_ocn` stays internal (from `f_grnd` and its regions reference, §8).

### 2.6 geothermal, sediments (regional; data)

Exports `Q_geo` [W m-2] (yelmox `gthrm%now%ghf` [mW m-2], datasrc `units_conv`
`*1e-3`) and `H_sed` [m] (`sed%now%H`): rank 2, INSTANT, intensive, required.
Cadence off, put once in `init_state`, `con%restart = .false.`.

### 2.7 forcing (global, 0-D; tsgen) and hub

forcing: one export per coefficient of `[forcing] coefs`, name = coefficient
name, units = its `units`, rank 0, INSTANT, intensive; imports `M_ice` (Gt,
from ice, lag 1, summed over domains) and one 0-D import per derived
coefficient, named like it, from its producer (§7).

hub (domain, not a component): `hub_exports` = `z_bed`, `H_ice`, `z_srf`,
`f_ice`, `f_grnd`, `z_sl` (INSTANT, every refresh) and `z_srf_ref` (once);
`hub_imports` = the mirror set (`z_bed`, `H_ice`, `z_srf`, `f_ice`, `f_grnd`,
`z_sl`) when `grid_hub == grid(ice)`, else the anomaly set (`dz_bed`,
`dH_ice`, `z_sl`) (§5).

## 3. Abstract types and wrappers

### 3.1 `kryos_component`

```fortran
integer, parameter, public :: PASS_INIT = 0, PASS_GEOMETRY = 1, PASS_BOUNDARY = 2, PASS_STATE = 3, &
                              PASS_SPINUP = 4, PASS_FINAL = 5, PASS_RESTART = 6   ! SPINUP: phase 2 (WP1)

type, public :: caps_t                          ! set by the backend in init
    logical :: has_ocn = .false., has_ocn_shelf = .false., has_smb = .false., has_Qd = .false.
end type

type, public :: constraints_t                   ! decision 12; one function per backend module
    logical           :: linked = .true.        ! .false. in the stub
    character(len=16) :: build_flag = ""        ! "yelmo=1"
    logical           :: regional = .true., global = .false.   ! allowed scopes
    integer           :: max_instances = 0      ! 0: unlimited
    character(len=64) :: fixed_grid = "", domain = ""          ! "": any
    logical           :: restart = .true.       ! .false.: no backend restart file
end type

! axis_t: name ("month" | "depth"), units, v(:). census_t: instances per "kind:backend" (run).
! grid_ref_t, init_ctx_t, step_time_t: core-design §2, plus init_ctx_t%axes(:) (producers'
! axes) and step_time_t%dt (the instance's cadence step [yr]; dtt for every step).

type, abstract, public :: component_t
    character(len=64)   :: name = "", kind = "", backend = "", grid = "", domain = ""
    character(len=16)   :: scope = "regional"   ! "regional" | "global"
    integer             :: prec = sp            ! kind of the backend's arrays (§9)
    integer             :: nx = 0, ny = 0
    type(cadence_t)     :: cadence
    type(caps_t)        :: caps
    type(constraints_t) :: con
    type(field_state_t) :: imports, exports
    type(axis_t), allocatable :: axes(:)
contains
    procedure(c_init),   deferred :: init         ! (ctx): read group; set prec, nx, ny, caps, opt flags
    procedure :: init_state => component_init_state      ! (st, pass): default nothing
    procedure :: prepare => component_prepare            ! (st): imports -> backend boundary; default nothing (phase 2)
    procedure(c_update), deferred :: update       ! (st): integrate st%t0 -> st%t1; put every export
    procedure(c_rst),    deferred :: restart_write, restart_read   ! (dir, t)
    procedure :: restart_files => component_restart_files  ! () -> base names restart_write creates in `dir`
    procedure :: write_output => component_write_output  ! (file, t, stream, init): exports, core §6
    procedure :: end => component_end                    ! (): default nothing
end type
```

Free procedures with `class(component_t)` dummies, used only inside component
modules: `backend_admit(c, census)` enforces `con` with named stops (decision
12): "instance `ice_grl` needs backend `yelmo`: rebuild kryos with `yelmo=1`",
"backend `rembo` allows one instance; a second was requested by domain
`Antarctica`", "backend `rembo` runs only on grid `GRL-16KM`", wrong scope or
domain; then counts the instance. `component_close(c, st)` =
`fields_finalize_interval(c%exports, st%t1)`. `kryos_mask` (§8).

`prepare(st)` (phase 2, WP1) puts the instance's imports into the backend's
boundary in the backend's units and precision, so that an experiment hook
(`HOOK_ICE_BND`) can change the boundary before `update`. Only the run
sequence calls it, through `domain_prepare(dom, k, clk, init)`, once per step
when the instance is due, before the hook and `update` (and at the cold start
before `PASS_STATE`). The ice is the only backend that implements it (yelmo:
`bnd_from_imports`); `update` of such a backend checks that `prepare` ran for
the step. `PASS_SPINUP` runs after `PASS_STATE` and the hook between them:
the cold-start equilibration of the instance (yelmo: `init_equil_time`,
`init_time_thrm`) and the exports of the initial state.

`restart_files(self) result(files)` (`character(len=PATH_LEN)`, allocatable)
returns the base names (no folder) of the files `restart_write` creates in the
folder it is given; kryos_restart lists them in the manifest and checks they
exist. Default `[<name>_restart.nc]` when `con%restart`, else none; a backend
whose files differ overrides it (tsgen forcing, Yelmo sidecars). Called only
by the `<comp>_restart_files(slot)` wrappers.

### 3.2 ice, in full

```fortran
module kryos_ice
    use kryos_component
    use fesm_fields, only: fields_init, fields_add, fields_lock, FIELD_INSTANT, FIELD_MEAN, &
                           FIELD_INTENSIVE, FIELD_EXTENSIVE
    implicit none
    private                                          ! public: the types and the wrappers below

    type :: ice_opt_t                                ! optional fields the backend supports (set in init)
        logical :: T_shlf = .false., tf_shlf = .false., Q_geo = .false., H_sed = .false., Qd = .false.
        logical :: uxy_s = .false., dz_bed_dt_corr = .false., z_bed_ref = .false.
    end type

    type, abstract, extends(component_t) :: ice_t
        type(ice_opt_t) :: opt
        integer :: hi_z_bed = 0, hi_z_sl = 0, hi_smb = 0, hi_T_srf = 0, hi_bmb_shlf = 0   ! imports
        integer :: hi_T_shlf = 0, hi_tf_shlf = 0, hi_Q_geo = 0, hi_H_sed = 0, hi_Qd = 0
        integer :: h_H_ice = 0, h_z_srf = 0, h_z_bed = 0, h_z_sl = 0, h_f_ice = 0, h_f_grnd = 0
        integer :: h_dHidt = 0, h_smb = 0, h_bmb = 0, h_cmb = 0, h_mb_net = 0, h_M_ice = 0
        integer :: h_ux_s = 0, h_uy_s = 0, h_dz_bed_dt_corr = 0, h_z_bed_ref = 0
    end type

    type :: ice_slot_t                               ! what domain_t holds
        class(ice_t), allocatable :: c
    end type

    interface                                        ! factory, submodule kryos_ice_factory
        module subroutine ice_create(ice, ctx, census)
            type(ice_slot_t), intent(out)   :: ice
            type(init_ctx_t), intent(in)    :: ctx  ! instance, backend, grid, domain, scope
            type(census_t),   intent(inout) :: census
        end subroutine
    end interface

contains                                             ! wrappers: the only call sites of the bindings
    subroutine ice_init(ice, ctx)
        type(ice_slot_t), intent(inout) :: ice
        type(init_ctx_t), intent(in)    :: ctx
        call ice%c%init(ctx)
        call ice_declare(ice%c)
    end subroutine

    subroutine ice_init_state(ice, st, pass)
        type(ice_slot_t),  intent(inout) :: ice
        type(step_time_t), intent(in)    :: st          ! t0 = t1 = time_init (core §3.2)
        integer,           intent(in)    :: pass
        type(step_time_t) :: s
        s = st;  s%t0 = st%t1 - st%dt                   ! virtual previous step for MEAN puts (§1)
        call ice%c%init_state(s, pass)
        call component_close(ice%c, s)
    end subroutine

    subroutine ice_update(ice, st)
        type(ice_slot_t),  intent(inout) :: ice
        type(step_time_t), intent(in)    :: st
        if (st%t1 == st%t0) return                      ! idempotent at dt = 0
        call ice%c%update(st)
        call component_close(ice%c, st)
    end subroutine
    ! ice_present(ice) = allocated(ice%c); ice_restart_write/read(ice, dir, t),
    ! ice_restart_files(ice) result(files) (character(len=PATH_LEN), allocatable),
    ! ice_write_output(ice, file, t, stream, init), ice_end(ice): one line each.

    subroutine ice_declare(c)                        ! the vocabulary of §2.1, once for all backends
        class(ice_t), intent(inout) :: c
        call fields_init(c%imports, trim(c%name)//":imports")
        call fields_add(c%imports, "smb", "kg m-2 yr-1", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_smb, &
                        grid=c%grid, dims=[c%nx, c%ny])
        if (c%opt%Q_geo) call fields_add(c%imports, "Q_geo", "W m-2", FIELD_INSTANT, FIELD_INTENSIVE, &
                                         c%hi_Q_geo, grid=c%grid, dims=[c%nx, c%ny])
        ! ... every import of §2.1
        call fields_lock(c%imports)
        call fields_init(c%exports, trim(c%name)//":exports")
        call fields_add(c%exports, "dHidt", "m yr-1", FIELD_MEAN, FIELD_EXTENSIVE, c%h_dHidt, &
                        grid=c%grid, dims=[c%nx, c%ny])
        call fields_add(c%exports, "M_ice", "Gt", FIELD_INSTANT, FIELD_EXTENSIVE, c%h_M_ice)   ! 0-D
        ! ... every export of §2.1
        call fields_lock(c%exports)
    end subroutine
end module kryos_ice

submodule (kryos_ice) kryos_ice_factory              ! the only place naming ice backends
    use kryos_ice_yelmo, only: yelmo_ice_t, yelmo_ice_constraints
contains
    module subroutine ice_create(ice, ctx, census)
        type(ice_slot_t), intent(out)   :: ice
        type(init_ctx_t), intent(in)    :: ctx
        type(census_t),   intent(inout) :: census
        select case (trim(ctx%backend))
        case ("yelmo");  allocate(yelmo_ice_t :: ice%c);  ice%c%con = yelmo_ice_constraints()
        ! phase 2: case ("data"); allocate(data_ice_t :: ice%c); ice%c%con = data_ice_constraints()
        case default
            call fesm_error("ice_create", "unknown ice backend '"//trim(ctx%backend)//"'", &
                            "instance "//trim(ctx%instance)//"; known: yelmo")
        end select
        ice%c%name = ctx%instance;  ice%c%kind = "ice";  ice%c%backend = ctx%backend
        ice%c%grid = ctx%grid_name; ice%c%domain = ctx%domain; ice%c%scope = ctx%scope
        call backend_admit(ice%c, census)            ! a stub stops here, before init
    end subroutine
end submodule
```

The backend file (`src/backends/ice/kryos_ice_yelmo.f90`, adapters-design)
holds one extension type with one `contains` block (`init`, `init_state`,
`update`, `restart_write`, `restart_read`, `write_output` for Yelmo's own
streams) and `yelmo_ice_constraints()` (`constraints_t(build_flag =
"yelmo=1")`); its `update` works by handle (`fields_get(self%imports,
self%hi_z_bed, self%ylmo%bnd%z_bed)`, `fields_put(self%exports, self%h_dHidt,
dhdt, t0=st%t0, t1=st%t1)`). The stub module has the same name and type, the
constraints with `linked = .false.`, and bindings that stop (safety net). The
factory stays a submodule of the component module (it `use`s the backends,
which `use` the component: no cycle). A host adds its own extension with
`ice_adopt(ice, proto, con, ctx, census)` (`allocate(ice%c, mold=proto)`),
without editing kryos. Experiments reach `ice%c` only by `select type`.

### 3.3 The other components, by difference

| Module | Abstract type, slot | Extra bindings and wrappers | Handles (exports / imports, §2) | Phase 1 backend; constraints |
|---|---|---|---|---|
| `kryos_geo` | `geo_t`, `geo_slot_t`; one type for both scopes, declaration by `scope` | `reduce(st)` (default no-op; global part), wrappers `geo_reduce`, `geo_scatter` (§6) | regional `h_dz_bed`, `h_dz_ss`, `h_bsl`, `h_dV_bsl`, `h_z_bed_ref`, `h_dz_bed_dt` / `hi_H_ice`, `hi_dz_bed_dt_corr`, `hi_bsl`, `hi_A_ocn`; global `h_bsl`, `h_A_ocn` / `hi_dV_bsl` | `fastisostasy`: two extension types in one file, `fastiso_geo_t` (regional, `isos_class`) and `fastiso_bsl_t` (global, `bsl_class`); the factory picks by (backend, scope); regional and global allowed, one global instance |
| `kryos_climate` | `climate_t`, `climate_slot_t` | — | `h_tas`, `h_pr`, `h_z_srf_clim`, `h_lapse`, `h_beta_p`, ocean and capability handles, `h_derived(:)` / geometry `hi_*`, `hi_coef(:)` with `coef_names(:)` | `snapshots`; regional |
| `kryos_surface` | `surface_t`, `surface_slot_t` | `init_state` used (chion spin-up in `PASS_BOUNDARY`, one cycle in `PASS_FINAL`) | `h_smb`, `h_T_srf`, optionals / `hi_z_srf`, `hi_H_ice`, `hi_z_bed`, `hi_z_sl`, `hi_tas`, `hi_pr` | `chion` (`chion=1`) |
| `kryos_shelf` | `shelf_t`, `shelf_slot_t` | — | `h_bmb_shlf`, `h_T_shlf`, `h_tf_shlf` / geometry, set A, set B | `marine_shelf` |
| `kryos_geothermal`, `kryos_sediments` | `geothermal_t`, `sediments_t` | `update` no-op (cadence off) | `h_Q_geo` / —; `h_H_sed` / — | `data` (`component_t` + `datasrc_t`, in the component module, design §13); `restart = .false.` |
| forcing (module name: Q4) | `forcing_t`, `forcing_slot_t` | `kill()`, `restart_due()` → `forcing_kill`, `forcing_restart_due` | `h_coef(:)` / `hi_M_ice`, `hi_derived(:)` | `tsgen`; global only, one instance, no stub |

The rembo1 case (phase 2) is only a constraints function:
`constraints_t(build_flag="rembo=1", max_instances=1, fixed_grid="GRL-16KM")`.

## 4. Wiring

### 4.1 Rows and transforms

```fortran
integer, parameter :: T_NONE = 0          ! dst = remap(src)
integer, parameter :: T_ADD_REF = 1       ! dst = ref_dst%<aux> + remap(src)                  (decision 4)
integer, parameter :: T_SUB_REF = 2       ! dst = remap(src - ref_src%<aux>)                 (hub anomalies)
integer, parameter :: T_ADD_SCALAR = 3    ! dst = src_state%<aux> (0-D) + remap(src)          (z_sl = bsl + dz_ss)
integer, parameter :: T_LAPSE_T = 4       ! tas re-elevated to the consumer's z_srf
integer, parameter :: T_LAPSE_P = 5       ! pr  re-elevated to the consumer's z_srf
integer, parameter :: LAG_INIT = -1       ! exchanged once after init (static fields)

type :: wire_t                            ! kryos_wiring; kinds: ice geo climate surface shelf
    character(len=16) :: src_kind = "", dst_kind = ""   ! geothermal sediments hub; globe: forcing geo
    character(len=FIELD_NAME_LEN) :: src_field = "", dst_field = "", aux = ""   ! aux: ref field or 0-D src field
    character(len=8)  :: method = ""      ! "" (from cons) | con | bilin | auto
    integer           :: lag = 0, transform = T_NONE    ! lag: 0 | 1 | LAG_INIT
    logical           :: gather = .false. ! N domains -> one global 0-D field, summed in domain-name order
    logical  :: active = .false.          ! resolved by wire_resolve:
    integer  :: hs = 0, hd = 0, ha = 0, prec_src = 0, prec_dst = 0
    real(dp), allocatable :: work(:,:,:), pre(:,:,:), aux3(:,:,:,:)   ! persistent workspaces
end type
pure function wire(src_kind, src_field, dst_kind, dst_field, method, lag, transform, aux, gather) result(w)
```

A closed enum evaluated by `kryos_wiring`: no user code, no scale factors (the
exchange never converts units). `T_LAPSE_*` are the identity when the two grids
are equal (yelmox returns before the lapse step, coupling:529); otherwise, per
month m, with `zc = bilin(z_srf_clim)`, `bp = bilin(beta_p)`, `lp =
bilin(lapse)` in `aux3` and `z` the consumer's current `z_srf` import (the core
runs hub → surface before climate → surface):

```
dz = z - zc;   tas = bilin(tas) - lp(:,:,m)*dz;   pr = bilin(pr)*exp(-bp*lp(:,:,m)*dz)
```

`lp` of a uniform field rounds back to the scalar `lapse(m)` in `sp`.

### 4.2 The phase 1 table (`kryos_wiring`, data)

Method `auto` = coupler `"auto"` without `norm` = yelmox
`remap_method_smooth` (bilin to refine, con otherwise; same choice for the
square projected cells of every phase 1 grid; core R3). Lag 0: the producer has
run this step; 1: from the previous step; I: once at init.

| src | field | dst | field | method | lag | transform | yelmox |
|---|---|---|---|---|---|---|---|
| ice | `H_ice` | geo | `H_ice` | auto | 1 | | coupling:146 |
| ice | `z_bed_ref` | geo | `z_bed_ref` | auto | I | | kryos_startup.f90:236 (`isos_init_ref`) |
| ice | `H_ice_ref` | geo | `H_ice_ref` | auto | I | | kryos_startup.f90:238 |
| globe geo | `bsl`, `A_ocn` | geo | `bsl`, `A_ocn` | — | 0 | | `bsl_update` before `isos_update` |
| geo | `dV_bsl` | globe geo | `dV_bsl` | — | 0 | gather Σ | fastisostasy.f90:839 |
| geo | `dz_bed` | ice | `z_bed` | auto | 0 | ADD_REF `z_bed` | coupling:222,225 |
| geo | `dz_ss` | ice | `z_sl` | auto | 0 | ADD_SCALAR `bsl` | coupling:223,226 |
| ice | `H_ice`, `z_bed`, `z_srf`, `f_ice`, `f_grnd`, `z_sl` | hub | same (mirror set) | "" (same grid) | 0 | | coupling:557-562 |
| ice | `z_bed` | hub | `dz_bed` | bilin | 0 | SUB_REF `z_bed` | coupling:565 |
| ice | `H_ice` | hub | `dH_ice` | bilin | 0 | SUB_REF `H_ice` | coupling:564,566 |
| ice | `z_sl` | hub | `z_sl` | bilin | 0 | | coupling:567 |
| hub | `z_srf`, `H_ice`, `z_bed`, `f_grnd`, `z_sl` | climate | same | auto | 0 | | coupling:430-434 |
| hub | `z_srf_ref` | climate | `z_srf_ref` | auto | I | | coupling:435 (static) |
| hub | `z_srf`, `H_ice`, `z_bed`, `z_sl` | surface | same | auto | 0 | | coupling:487-490 |
| hub | `H_ice`, `z_bed`, `f_grnd`, `z_sl`, `z_srf` | shelf | same | bilin | 0 | | coupling:592-596 |
| climate | `tas` | surface | `tas` | bilin | 0 | LAPSE_T | coupling:527-535 |
| climate | `pr` | surface | `pr` | bilin | 0 | LAPSE_P | coupling:528-536 |
| climate | `to_ann`, `so_ann`, `dto_ann` | shelf | same | bilin | 0 | | coupling:610-613 |
| climate | `to_shlf`, `so_shlf`, `dto_shlf`, `dso_shlf` | shelf | same | bilin | 0 | | coupling:600-603 |
| surface | `smb`, `T_srf` | ice | `smb`, `T_srf` | con | 1 | | coupling:314-315 |
| shelf | `bmb_shlf`, `T_shlf` | ice | same | con | 1 | | coupling:353-354 |
| climate | `Qd` | ice | `Qd` | con | 1 | | coupling:182 |
| geothermal | `Q_geo` | ice | `Q_geo` | con | I | | kryos.f90:356 |
| sediments | `H_sed` | ice | `H_sed` | con | I | | kryos.f90:353 |
| ice | `M_ice` | globe forcing | `M_ice` | — | 1 | gather Σ | yelmox.f90:185 |
| climate | `<derived>` | globe forcing | `<derived>` | — | 1 | | yelmox_climate.f90:338 |
| globe forcing | `<name>` | climate | `coef_<name>` | — | 0 | | `coefs_value` in `terms_eval` |

The mirror set and the anomaly set are alternatives, chosen by
`grid_hub == grid(ice)` at resolution. The two forcing rows with `<...>` are
rules: the destination's declared names select the rows (a climate import
`coef_alpha` reads forcing `alpha`; a forcing import `dTo` reads the climate
export `dTo`). Not in phase 1: `climate.smb, T_srf → surface` (backend
`surface = "climate"`), `shelf.tf_shlf → ice.tf_shlf`, all of §10.

### 4.3 `wire_exchange`

```fortran
subroutine wire_exchange(w, src, dst, cpl, t, t_last, ref_src, ref_dst)
    ! w inout; src in; dst inout; cpl inout; t = clock time, t_last = producer's
    ! cadence_t_last (hub: hub_t); ref_src, ref_dst: dom%refs of the two grids (optional)
    associate (fs => src%f(w%hs), fd => dst%f(w%hd))
    call lag_check(w, fs, t, t_last)                    ! §4.4
    if (fs%rank == 0) then
        call fields_put(dst, w%hd, fs%v(1,1,1), t0=fs%t0, t1=fs%t1);  return
    end if
    if (w%transform == T_SUB_REF) then                  ! anomaly in the source's precision
        call sub_ref(w%pre, fs%v, ref_src, w%aux, w%prec_src)
        call remap(cpl, w%pre, fs%grid, w%work, fd%grid, method=w%method)
    else                                                ! identity copy when grids are equal
        call remap(cpl, fs%v, fs%grid, w%work, fd%grid, method=w%method)
    end if
    ! no rounding: the consumer narrows once, after its own conversion (§9)
    select case (w%transform)                           ! in the destination's precision (§9)
    case (T_ADD_REF);    call add_ref(w%work, ref_dst, w%aux, w%prec_dst)
    case (T_ADD_SCALAR); call add_scalar(w%work, src%f(w%ha)%v(1,1,1), w%prec_dst)
    case (T_LAPSE_T, T_LAPSE_P); call lapse(w, src, dst, cpl, w%prec_dst)
    end select
    call fields_put(dst, w%hd, w%work, t0=fs%t0, t1=fs%t1)
    end associate
end subroutine
```

The core's `exchange(dom, src, dst, t)` picks states, references and `t_last`
and calls this per active row in table order. `gather_into(globe, doms, kind,
t)` sums, in domain-name order, the 0-D values published at `t` (stale ones
count 0) and puts the sum.

### 4.4 Resolution at init (`wire_resolve`), log, lag check

Each row is active iff both fields are declared, else inactive with a reason.
Then, stopping with a named message on failure:

1. **Gaps**: every required import has exactly one active row (forcing rule
   included): "ice `ice_grl` import `bmb_shlf` has no producer: row
   shelf.bmb_shlf inactive (domain `Greenland` has no shelf instance)"; two
   active rows into one field also stop.
2. `fields_check_link` per active row; rank 3: same `zdim`, identical axis.
3. No masked field across two grids; MEAN → INSTANT rows need equal cadences.
4. Component contracts: shelf ocean set (§2.5); a regional `fastisostasy` geo
   needs a global geo part; phase 2 capability rules (`surface = "climate"`
   needs `has_smb`).
5. `coupler_prime` every remapped pair (core); with the frozen cache a new pair
   stops (regression §5 row 1).

Logged once (`LOG_INFO`, after `fields_describe` of every state): one line per
row, `src.field [grid] -> dst.field [grid]  method (resolved)  lag  transform
active|inactive: reason`. The resolved `auto` needs a public
`coupler_method(cpl, src, dst, method)` (wrapper of the private
`resolve_method`; one-line fesm-utils addition). Nothing is logged per
exchange.

Lag check (as core §4): the published `t1` equals `t_last` (the producer put
what it produced), and `t1 ≤ t` with the producer not due (lag 0) or `t1 < t`
(lag 1); during the cold start (`t = t0`) every row counts as lag 0;
`LAG_INIT` rows run once. Rows run where the core sequence calls `exchange`;
the lag column documents and enforces that order, it does not schedule.

## 5. The hub

### 5.1 What it holds

Per domain: `htopo_class` (`sp` in phase 1, `dp` after the gate as a
rebaseline, regression Q2) with the reference `ref` (z_bed, H_ice, z_srf,
z_bed_sd from FesmData on `grid_hub`), regions, basins and the current
geometry; `hub_imports`, `hub_exports` (§2.7), `hub_t`; `refs(:)`, the
reference on each instance grid, made once by `con` through the coupler's `sp`
path and stored widened (core R4), only for grids in use (no map pair yelmox
lacks).

### 5.2 Refresh after the ice step (`domain_hub_update`)

Rows ice → hub (§4.2) fill `hub_imports`, then:

```
if (grid_hub == grid(ice)) then          ! mirror: the hub is the ice geometry, fractional f_grnd kept
    z_bed, H_ice, z_srf, f_grnd, z_sl, f_ice = the mirror imports
else                                     ! hi-res reference + the ice's anomalies (bilin)
    dz_bed = bilin(z_bed_ice - refs(ice)%z_bed)     ! T_SUB_REF in sp, as yelmox
    dH_ice = bilin(H_ice_ice - refs(ice)%H_ice)     ! refs(ice)%H_ice = con(ref%H_ice) = yelmox H_ice_ref_y
    z_sl   = bilin(z_sl_ice)
    call htopo_update(hub, dz_bed, dH_ice, z_sl)    ! z_bed = ref + dz; H = max(ref + dH, 0);
                                                    ! f_grnd in {0,1} by flotation; z_srf
    f_ice  = merge(1, 0, H_ice > 0)
end if
put hub_exports at t; hub_t = t
```

This is `couple_yelmo_to_htopo` (coupling:541-570) with Yelmo's
`bnd%z_bed_ref` replaced by `refs(ice)%z_bed`, the same values (Yelmo's
`z_bed_ref` is set from it at init, decision 4). The hub's bed comes through
the ice (the bed the ice actually used; CISM's un-folded `z_bed`), not from geo
directly (Q3). At init `htopo_init` sets the current geometry to the
reference and `hub_exports` (with `z_srf_ref`) is put at `t0`.

### 5.3 Reference consistency (decision 4)

`domain_check_reference(dom)` after geo's `PASS_GEOMETRY` and after every
restore: geo's `z_bed_ref` export against `refs(geo)%z_bed`, and the ice's
`z_bed_ref` (if declared) against `refs(ice)%z_bed`; log mean and max |diff|;
stop above 1 m (yelmox `check_isostasy_reference`, same message: on a restart
the geo or ice restart came from another topography). On a cold start both are
exact by construction.

## 6. geo: two levels in phase 1

Instances: one regional `fastiso_geo_t` per domain (`[geo_grl]`, on the geo
grid, `isos_class`), one global `fastiso_bsl_t` (`[geo_glob]`, `bsl_class`,
modes `const | file | fastiso` = decision 5's `const | file | own`). Per step
(core §3.4, with the corrections of the header):

```
geo_update(globe%geo, st)                    ! bsl_update(bsl, st%t_eval): const/file value; A_ocean from bsl
do d, if geo due:
    exchange(dom, "ice", "geo", t)           ! H_ice, dz_bed_dt_corr: H(t_{n-1}), lag 1
    exchange_globe(globe, "geo", dom, "geo", t)   ! bsl, A_ocn of this step
    geo_update(dom%geo, st)                  ! adapter: local bsl_class (bsl_now = bsl, A_ocean_now = A_ocn,
                                             ! method as the global's); isos_update(isos, H, t, bsl_loc, dwdt_corr);
                                             ! puts dz_bed = w+we, dz_ss = z_ss - now%bsl, dV_bsl = now%deltaV_bsl
gather_into(globe, doms, "geo", t)           ! dV_bsl: Σ over domains whose geo ran at t, name order
geo_reduce(globe%geo, st)                    ! fastiso: bsl_now = bsl_now - real(Σ dV, wp)/A_ocean_now; put bsl
do d, if geo ran at t:
    exchange_globe(globe, "geo", dom, "geo", t)   ! the reduced bsl into the regional import
    geo_scatter(dom%geo, st)                 ! wrapper: regional export bsl := import bsl
```

Consumers form `z_sl = bsl + remap(dz_ss)` from the regional exports, so
domain rows never read the globe. A domain whose geo did not run keeps its
exports, including the `bsl` it last solved with (yelmox uses
`isos%now%bsl`, which changes only inside `isos_update`).

Bit identity, one domain (regression §5 row 10, core R5): `isos_update`
computes in `sp` `bsl_loc = bsl - deltaV/A` and `z_ss` with it; `geo_reduce`
computes `bsl - real(dV_1, sp)/A` (`real` exact), the same `sp` expression, so
the reduced `bsl` equals `isos%now%bsl` and `z_sl` at the ice equals
`isos%now%bsl + dz_ss_y`. For N > 1 the sum is order-independent; each
domain's `z_ss` saw only its own contribution within the step (the others'
arrive one step later; logged when N > 1). VILMA2 (phase 4) uses the same
calls: `geo_update` no-op, gather merges `H_ice` (`coupler_merge_*`),
`geo_reduce` solves, scatter rows split `dz_bed`, `dz_ss` (`coupler_split`).

## 7. forcing

- Global, 0-D, every step, `clock = "forcing"`: `st%t_eval` = `time_rel` at
  t1, as `coefs_update(coefs, ts%time_rel, var)`.
- Backend `tsgen` ports `coefs_class` under a new module name (Q4): `[forcing]
  coefs`, `[coef_<name>]` through the par library (core §5); tsgen
  coefficients advance with `var = M_ice` (lag 1, Σ over domains); derived
  ones copy their import; all are put at t1.
- Derived coefficients: the producer reads `from`, `factor`, `region` from the
  same `[coef_<name>]` group and exports a 0-D field named like the
  coefficient. The climate uses its own value inside its update (yelmox
  `coefs_set`, then `terms_eval(derived=.true.)`), so the forcing's lag-1 copy
  serves other consumers (OBM, phase 3). Two producing domains: init stop
  until phase 3.
- Snapshots: its terms name coefficients in `init` (`self%coef_names`); the
  climate layer declares 0-D imports `coef_<name>` (`hi_coef(:)`); the forcing
  rule wires them, a missing one is a named gap ("climate `clim_grl` term
  `lgm` needs coefficient `alpha`: not in `[forcing] coefs`"); `update` reads
  `fields_get(self%imports, self%hi_coef(k), f)`, exact (`sp` widened).
- Restart: `forcing_restart.nc` (every tsgen state; core §7) and
  `forcing_fields.nc` (derived values). Wrappers `forcing_kill` (tsgen kill
  switch) and `forcing_restart_due` (`restart_every_df`) for the run. Output:
  exports to `forcing_ts.nc` (`coef_<name>`); `df_dt`, `dv_dt` as a backend
  stream.

## 8. Masks and fractions

Fractions cross, integer masks never. `f_ice`: ice-covered fraction of a
cell; `f_grnd`: grounded fraction, ice-free land above sea level counting as
grounded (Yelmo, `htopo`). The canonical mask, in `kryos_component`:

```fortran
elemental integer function kryos_mask(f_ice, f_grnd)     ! majority class
    ! MASK_GRND f_ice >= 0.5 and f_grnd >= 0.5 | MASK_FLOAT f_ice >= 0.5 |
    ! MASK_LAND f_grnd >= 0.5 | MASK_OCEAN otherwise
```

On another grid the fractions are remapped and the function applied there.
Phase 1 consumers:

| Consumer | Needs | Served by |
|---|---|---|
| marine_shelf | `f_grnd` (its `mask_ocn` = `f_grnd` + its regions' `mask_ocn_ref`) | hub → shelf `f_grnd` (bilin) |
| FastIsostasy | its own masks from `H_ice` and `z_ss` | nothing extra |
| chion | land/ice from `z_bed`, `H_ice`, `z_sl` | hub geometry |
| snapshots (esm in phase 2) | `f_grnd`, `z_sl` | hub → climate |
| hub output, future consumers (pico, OBM freshwater, CMIP) | four-value mask | `kryos_mask` on `hub_exports`; written as `mask` in the hub file |

## 9. Precision of transforms and bit identity

Rule (decision 2026-10-25, replacing "a remapped value is rounded to the
destination's `prec`"): a row never rounds; the consumer narrows once, after
its own conversion. A `T_SUB_REF` difference is formed in the source's `prec`;
every other transform is evaluated in the destination's `prec` with the yelmox
operand order (`sp` and `dp` kernels), an `sp` kernel narrowing its operands:
they are `sp` values widened, or `dp` remaps of them, which narrowed are the
coupler's `sp` remap (R1). `T_LAPSE_P` computes its factor in `sp` and applies
it to the unnarrowed `pr` (which carries a unit conversion the consumer
reverses; `dz = 0` is exact). Rounding every row to the consumer's `sp` broke
the unit round trips (`pr`, `bmb_shlf`, `Q_geo`: the `dp` product was rounded
before the consumer divided the factor back, C9 test). With `sp` operands,
`sp` arithmetic is what yelmox did, so:

| Where | yelmox | kryos | Equal because |
|---|---|---|---|
| geo → ice `z_bed` | `z_bed_ref + dz_bed_y` (sp) | `real(ref,sp) + real(remap(dz),sp)` | same `sp` operands and op (core R4: refs via the sp path) |
| geo → ice `z_sl` | `isos%now%bsl + dz_ss_y` | geo `bsl` (reduced, §6) + remapped `dz_ss` | §6 |
| ice → hub anomalies (Antarctica) | `bnd%z_bed - bnd%z_bed_ref`, sp, then sp bilin | SUB_REF in sp, dp bilin, narrowed by `htopo_update`'s sp | R1 (sp remap = dp remap + one rounding) |
| hub refresh | `htopo_update` in sp | same routine, sp hub | unchanged code |
| climate → surface | copy (identical grids in both gates) | identity | — |
| unit round trips | — | ×ρ then ÷ρ, ×360 then ÷360, ×1e-3 then ×1e3 | exact for `sp` values (regression §3.3) |

Remaining non-identity: `bmb_shlf` ×ρ and ÷ρ around the shelf → ice `con` in
Antarctica (regression §5 row 5, core Q2). This rule answers the hub part of
core Q2 without moving any sum into an adapter (Q1). After the gate, `prec =
dp` for the hub is a declaration change.

## 10. In yelmox, not in phase 1

| Field / edge (yelmox) | Phase 1 | Goes to (phase 2 unless noted) |
|---|---|---|
| `relax` (`ytopo topo_rel` ramp), `opt` cb_ref (`optimize_cb_ref`) | dropped (gate: `relax = opt = False`, regression Q1) | `experiments/relax`, `opt_cb_ref` (yelmo plugins; `HOOK_STEP_BEGIN`, cb_ref guess at `HOOK_COLD_PRE_STATE`) |
| `opt` tf_corr (marine_shelf `tf_corr`; reads Yelmo `H_ice`, `H_grnd`, `dHidt`, pd `H_ice`) | dropped | marine_shelf plugin (decision 6); needs ice `dHidt`, `f_grnd` exports and pd data |
| NEGIS `cb_ref` | dropped | `experiments/negis` |
| glacial smb scaling (`ta_ann` now/ref → ice smb) | dropped | `experiments/glacial_smb` on the surface → ice row; climate exports `ta_ann`, `ta_ann_ref` |
| `lim_pd_ice` (−4 m ie/yr outside pd ice) | dropped | yelmo plugin at `HOOK_ICE_BND`, after `prepare` |
| cold-start recipes (recon, marine ice, kill shelves) | only `init_equil_time` (core §3.2) | plugins at `HOOK_COLD_POST_STATE` (core §3.9) |
| thermal equilibration (`init_time_thrm`) | dropped | yelmo key `init_time_thrm` in `PASS_SPINUP` (WP1) |
| `ice_bnd` (`H_ice_ref`, `mask_ice`, `tau_relax` from a history, coef-scaled) | dropped | yelmo plugin reading a `datasrc_t` |
| smb_simple (`tsl_ann`; Yelmo `H_ice_ref` reach-in) | dropped | backend; climate export `tsl_ann`; mask from `refs(surface)%H_ice` |
| `surface = "climate"` (`has_smb`), esm (`has_ocn_shelf`, `has_Qd`, `Qd`), rembo, smbpal | dropped | backends; `Qd` per-cell units revisited with esm |
| `shelf.tf_shlf → ice` | no row | with CISM/PISM (decision 6) |
| bsl `mixed` (file + accumulated fastiso) | not offered (unused in every yelmox config) | only if a run needs it |
| OBM: Yelmo → fwf, `obm_north/south` derived coefficients | — | phase 3, `ocean` component from ice exports `mb_net`, `smb`, `bmb`, `cmb`, `dHidt` |
| CMIP output, `isos_ts` framework diagnostics, climate `ref` anomaly output | — | post-processing; backend streams |
| `with_* = False` | — | `data` backends |

## 11. Open questions for Alex

1. **Transform precision** (§9): evaluate transforms in the declared precision
   of the arrays yelmox held (`sp` kernels next to `dp`, remap results rounded
   to the consumer's precision), which makes the Antarctica hub anomalies
   bit-identical and leaves only `bmb_shlf` (proposed); or keep transforms in
   `dp` and accept ulp differences there (core Q2 as written)?
2. **Lapse rate on identical grids**: yelmox re-elevates `tas`/`pr` only when
   climate and surface grids differ; on one grid the surface-elevation change
   since the last climate update (up to `dt_clim`) is ignored. Reproduce in
   phase 1 (proposed, needed for the gate) and apply it on every grid
   afterwards as a rebaseline, or keep yelmox's behaviour?
3. **geo reads the ice, not the hub**: design §5 says geo now reads the hub;
   decision 4 and the gate need `ice.H_ice → geo` (Antarctica: geo on the
   32 km ice grid, hub 16 km; reading the hub's reconstructed `H_ice` would
   change the load beyond round-off, and `dz_bed_dt_corr` is not on the hub).
   Correct §5 to "geo reads the ice exports" (proposed)?
4. **Module name**: `kryos_forcing` is taken by yelmox `libs/kryos_forcing.f90`
   (design §1: no reuse while both build together). Name the component
   module `kryos_forcing_comp` (proposed) or another name? (Core-design's
   `kryos_output` collides with yelmox's `kryos_output` in the same way.)
5. **Experiment hooks** (decided, phase 2 WP1, `phase2-plan.md`): `prepare` on
   `component_t` and a core `PASS_SPINUP` (§3.1, core §3.2, §3.9), so that
   `lim_pd_ice` can act on the converted boundary and the cold-start recipes
   sit between the state and the equilibration; both bit-neutral on the
   phase 1 gates (same calls, same order).
