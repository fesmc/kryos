# kryos: design proposal

Status: design settled 2026-10-24; no code yet. Evidence: the four studies in
`study/` (yelmox, CLIMBER-X, fesm-utils, the backends). Section 15 records the
thirteen decisions, each with its reasoning and the rejected alternatives.

## 1. What kryos is

A Fortran 2008 framework for cryosphere modelling on one or more regional
domains, with optional global components. Every component (ice sheet included)
is an interchangeable backend behind a small, model-agnostic interface, and
every component has a `data` backend that reads its outputs from files. One
static library (`libkryos.a` plus `.mod` files, the archive-target convention
of the stack) that other programs can link, one driver program for N domains
that is itself only a client of the library, analysis tooling in Julia.

kryos replaces yelmox. The in-tree prototype `yelmox/libs/kryos*.f90` is the
reference for what must be generalised; its module names must not be reused
(the two may be built together during migration).

## 2. Vocabulary

| Noun | Meaning |
|---|---|
| **component** | A role with a fixed, canonical set of import and export fields: `ice`, `geo` (solid earth + sea level, after CLIMBER-X), `climate`, `surface` (SMB + surface temperature), `shelf` (sub-shelf melt), `ocean`, `geothermal`, `sediments`, `forcing` (run-level scalars). |
| **backend** | One implementation of a component, chosen at runtime by name: a prognostic model (`yelmo`, `cism`, `pism`, `fastisostasy`, `vilma2`, `chion`, `smbpal`, `rembo`, `marine_shelf`, `obm`), a parametrisation (`smb_simple`), or `data`. Optional libraries sit behind same-interface stubs. |
| **instance** | A component with a backend, a grid and a scope. Named in the par file; the name is its namelist group, output file prefix and restart file name. |
| **domain** | A regional scope: a named hub grid, the hi-res geometry hub, per-component grids, its regional instances, its own map cache, and its footprint on every global grid. |
| **globe** | The one global scope: global grids, the global geometry (reference topography + merged ice from the domains + background ice outside them), and the global instances (`vilma2`, `obm`, the barystatic sea level, the forcing scalars). |
| **run** | Clock and scheduler, the globe, the domains, output schedules, the restart bundle, and the coupling sequence in the driver. |

Scope is a property of the instance, not of the component type (section 6).

## 3. Architecture

```
run
 ├─ clock + scheduler        one integer-tick clock; a cadence per instance
 ├─ globe                    global grids; global geometry hub; footprints
 │    ├─ forcing             scalars at the forcing time: indices, CO2, orbit, bsl (if prescribed)
 │    ├─ geo       (global)  vilma2 | data | none        ← gather loads, solve, scatter dz_bed/dz_ss/bsl
 │    └─ ocean     (global)  obm | data | none           ← gather FWF, solve, scatter T/S
 ├─ domain(1..N)             e.g. "Greenland"
 │    ├─ hub (htopo)         hi-res geometry truth: reference + anomalies from the models
 │    ├─ coupler             grids by name, lazy cached maps (fesm-utils)
 │    ├─ ice                 yelmo | cism | pism | data
 │    ├─ geo      (regional) fastisostasy | data        (its global part: bsl reduction)
 │    ├─ climate             snapshots | esm | rembo | data
 │    ├─ surface             chion | smbpal | smb_simple | climate | data
 │    ├─ shelf               marine_shelf | pico | data   (owns ocean → shelf-base interpolation)
 │    ├─ geothermal, sediments   data
 │    └─ experiments         opt, relax, negis, ... (plugins; backend-specific)
 ├─ exchange                 field registry (name, units, grid, time semantics); wiring tables per component pair
 ├─ output                   per instance, native grid, from its export list
 └─ restart                  bundle folder per time, file per stateful instance, manifest
```

Dependency flow is one-directional: fesm-utils → kryos core → components →
backends → driver. Backends never `use` the core's run or domain types; the
core never `use`s a backend module.

## 4. The component contract

```fortran
type, abstract :: component_t
    character(len=64) :: name, kind, backend, grid     ! instance, role, backend, native grid tag
    character(len=16) :: scope                          ! "regional" | "global"
    type(field_state_t) :: imports, exports             ! declared at init (section 5)
    type(cadence_t)     :: cadence                      ! dt, phase, t_last
    type(caps_t)        :: caps                         ! capability flags (has_smb, has_ocn_shelf, ...)
contains
    procedure(c_init),    deferred :: init          ! (par_file, group, grid, time, cnst)
    procedure(c_update),  deferred :: update        ! (t_target): integrate from t_last to t_target; idempotent at dt=0
    procedure(c_rst_w),   deferred :: restart_write ! (dir, time)
    procedure(c_rst_r),   deferred :: restart_read  ! (dir, time)
    procedure                      :: write_output  ! default: write exports on the native grid from a variable table
    procedure                      :: end
end type
```

- One abstract type per component extends `component_t` with the canonical
  field handles (`h_H_ice`, `h_z_srf`, ...) so backends fill them by handle, not
  by name lookup, in `update`.
- A **factory** per component allocates the concrete type from the backend
  string. An unlinked optional backend is a stub extension whose `init` stops
  with "rebuild with `cism=1`".
- **Global instances** implement two extra steps, `gather(domain)` and
  `scatter(domain)`, around `update`; regional instances do not (section 6).
- **Two-step init** is allowed (Yelmo needs topography before its state): the
  domain initialises every instance's geometry from the hub first, lets the
  surface and climate produce initial boundary fields, then completes the
  state (`init_state`).
- **Capabilities** are declared at init and checked against the wiring: an
  `ocean` or `climate` with `has_ocn_shelf` lets `shelf` skip its own interpolation; a `climate`
  with `has_smb` allows `surface = "climate"`. Invalid combinations fail at
  init with a named reason (yelmox `climate_caps_class` pattern).
- **Backend internals** are reachable only through `select type` in an
  experiment plugin (`experiments/opt_cb_ref.f90` requires `yelmo_ice_t`
  and stops otherwise). The core never does this.
- **Precision**: each backend keeps its own `wp`; the contract is `dp`. Each
  backend takes the run's `phys_const_class` if it can; the adapter calls
  `phys_const_compare` and logs the differences otherwise (VILMA2, CISM, rembo).

## 5. Exchange

The field registry (`fesm_fields`, new in fesm-utils) holds `dp` buffers with
name, canonical units, grid tag, validity interval, time semantics
(instantaneous | mean over the interval | accumulated), and a conservation flag
(extensive fields are remapped `con`, intensive ones `bilin`). Names are resolved
to integer handles at init.

Canonical vocabulary (full table in `study/backends.md` section A):

| Component | Imports | Exports |
|---|---|---|
| ice | `z_bed`, `z_sl` (2D), `smb`, `T_srf`, `bmb_shlf`; opt `T_shlf`/`tf_shlf`, `Q_geo`, `H_sed`, `Qd` | `H_ice`, `z_srf`, `f_ice`, `f_grnd`; interval-mean `dHidt`, `smb`, `bmb`, `cmb`, `mb_net`; opt `ux_s`, `uy_s`, `Q_b` |
| geo | reference `z_bed_ref`, `H_ice_ref` (init); `H_ice` | `dz_bed`, `dz_ss`, scalar `bsl`; opt `dz_bed_dt` |
| climate | geometry; forcing scalars | monthly `tas`, `pr`, the `z_srf` they refer to, lapse rate; opt ocean `T`,`S` profiles or shelf-base values, `smb`, `T_srf`, `Qd` |
| surface | geometry; `tas`, `pr`; opt radiation, wind, humidity (chion) | `smb`, `T_srf`; opt `melt`, `runoff`, `refrz`, `albedo` |
| shelf | geometry; ocean `T`,`S` | `bmb_shlf`, `T_shlf`, `tf_shlf` |
| ocean | freshwater flux per domain | `T`,`S` profiles or box values |
| geothermal, sediments | — | `Q_geo`, `H_sed` |
| forcing | ice volume (feedback) | named scalars: indices, CO2, orbit, `bsl` if prescribed |

Conventions: geometry in m, time in years (`dp`), temperatures in K, heat flux
in W m-2 positive into the ice, every mass flux in `kg m-2 yr-1` and
"+ = ice gains" (adapters flip CISM/PISM melt and CISM `bheatflx`), fractions
in [0,1] instead of integer masks (the core derives one canonical four-value
mask), a `valid` mask instead of in-band missing values.

`kg m-2 yr-1` (decided 2026-10-24) is density-free and per year, and it is
numerically identical to the mm w.e./yr that yelmox already exchanges, so the
phase 1 bit-identity gate holds. `rho_ice` enters only inside each ice-sheet
adapter (Yelmo, CISM: m ice/yr or m ice/s), with that backend's own density;
`sec_year` enters only at the chion, CISM and PISM boundaries. Rates per
second were rejected: not for precision (relative precision is
magnitude-independent; only time in seconds would suffer, and time stays in
years), but because every Fortran backend in the stack is per year, the round
trip is not bit-exact, and the exchange layer must not carry a year-length
assumption.

Wiring between components is a **table per component pair, fixed in code**
(source field → target field, remap method, lag): `geo → ice:
z_bed = hub%z_bed_ref + con(dz_bed)`, `surface → ice: smb con, lag 1`,
and so on. The driver sequence calls `exchange(dom, "surface", "ice")`.
This keeps coupling readable and non-configurable by users, while adding a
backend touches no wiring (yelmox's `couple_*_to_yelmo` code is replaced by
data, but not by user configuration).

The **hub** (`htopo`, moved to `dp`) stays the geometry source of truth per
domain: reference topography from FesmData, refreshed after each ice step with
the ice model's anomalies, read by climate, surface and shelf (and now also by
geo, which today reads Yelmo directly). Slow, smooth fields (`dz_bed`,
`dz_ss`) always cross grids as anomalies; the hub holds the reference.

## 6. Regional versus global

**Rule:** an instance is regional (lives in a domain, on a domain grid) or
global (lives in the globe, on a global grid). Which it is depends on the
backend, not on the component:

| Component | Regional backends | Global backends | Notes |
|---|---|---|---|
| ice | yelmo, cism, pism, data | data (e.g. ICE-6G as the background) | a global prognostic ice model is not planned but fits |
| geo | fastisostasy, data | vilma2, data | fastisostasy per domain + a global reduction for `bsl` |
| climate | rembo, snapshots, esm, data | data (global products read once, split to domains) | snapshots/esm today read global files per domain; either works |
| surface | chion, smbpal, smb_simple, data | — | column physics; grid is the domain's |
| shelf | marine_shelf, pico, data | — | |
| ocean | data | obm, data | |
| geothermal, sediments | data | data | |
| forcing | — | tsgen, data | 0-D |

**Mechanism** (the VILMA2 plan's load → solve → land, made general):

1. `gather`: each domain exports its contribution to the global instance's
   grid. Extensive fields (ice load, freshwater flux) are remapped `con` with
   **full-target-area normalisation** and a **coverage fraction** `f_d`; the
   global field is `(1 − Σ f_d)·background + Σ f_d·contribution_d`. The
   background is the globe's own data instance (reference ice, or an ice
   history). Overlapping domains are resolved by owner weights, not by call
   order. Volumes are checked (`coupler_integral`) and logged.
2. `update`: the global instance integrates once per run step. For
   FastIsostasy the global part is only the barystatic sum
   `bsl += Σ_d ΔV_d / A_ocn` (order-independent, unlike today's sequential
   mutation).
3. `scatter`: each domain imports the global exports on its hub, as anomalies
   (`dz_bed`, `dz_ss`) added to the hub reference, or as remapped fields
   (ocean `T`, climate anomalies).

The globe keeps a global geometry hub (reference topography with the domain
references spliced in, merged ice) so that global components see a complete
Earth; a run with no global instance has an empty globe and pays nothing.

What is necessarily regional: anything whose physics is defined on the ice
domain's projected grid (ice dynamics, shelf melt, SMB columns). What is
necessarily global: the sea-level equation, orbital/CO2/index scalars, a box
ocean. Everything else is a backend choice, and the framework needs both scopes
for most components.

## 7. Time

- One clock (`fesm_clock`): integer ticks in `int64` internally, `dp` years at
  every interface, a declared calendar. Removes the `mod(nint(t*100))` cadence
  arithmetic, the 1e-6 output-time tolerance and the chunk-boundary rules.
- One **forcing time** `time_rel` (yelmox's), with the two named exceptions
  (fixed orbit; calendar-driven ESM periods).
- Each instance has a **cadence** (`dt`, phase, `t_last`) and is called with
  `update(t_target)`; adapters loop fixed-step backends (CISM, chion, smbpal)
  and compute `dt` for advance-by-dt backends (VILMA2, obm). `update` is
  idempotent at `dt = 0`.
- **Lag is declared** in the wiring table, not implied by call order. Default:
  bed and sea level end-of-step, surface and shelf fluxes from the previous
  step (today's scheme), written down.
- Restarts only at times that are multiples of every cadence; the scheduler
  state (`t_last` per instance) and the exchange accumulators are part of the
  bundle, so a restart is a restore and never a step.

## 8. The data backend

`datasrc_t` (new in fesm-utils, on `varslice`): one varslice per export field,
modes static | climatology | transient, time axis and units read from the file
and checked against the field's declared units, served into the export state
with the right validity interval. Every component's `data` backend is
`component_t + datasrc_t` and a few lines. This replaces the piecemeal cases of
today (`ice_bnd`, `with_ice_sheet = False`, `bsl method = file`, geothermal,
sediments, CLIMBER-X `fake_*`) with one mechanism, and gives for free the
isostasy-only, climate-only and SMB-only run types.

## 9. Configuration, output, restart, build

- **Par files**: `nml` with a defaults file per backend as schema
  (`<backend>_defaults.nml`, synchronised copies for external libraries as
  chion/VILMA2 do) and mandatory `nml_validate`. One group per instance, named
  by the instance: `[Greenland]` lists its instances; `[ice_grl]` has
  `backend = "yelmo"`, `grid = "GRL-16KM"`, `dt = 1.0`; backend physics groups
  are named in it. No hemisphere suffixes threaded through every group; a
  second domain is a second set of instance groups. Shared library groups
  (snapshots, indices) via `par_library` as now.
- **Output**: each instance writes its export fields on its native grid from a
  variable table (`variable_io`), on `timeout` schedules, with fixed per-field
  precision; backend-internal diagnostics stay the backend's own file. CMIP
  output becomes a post-processing table over the ice exports.
- **Restart**: `restart-<time>/` with one file per stateful instance and a
  manifest (instance, backend, grid, file, git hash); global instances at the
  top, domains in subfolders. Derived state (hubs) is rebuilt.
- **Build**: configme manifest with pinned dependencies (keep); every backend
  library optional with a stub (`yelmo=1 cism=0 pism=0 vilma=0 rembo=0 chion=1`
  defaults); generated dependency rules (no hand-written per-object lists);
  archive targets as dependencies; `make check` that runs tests. PISM via a
  C shim and `bind(C)`; a string-keyed C API for kryos itself (Yelmo C API
  precedent) so a Julia host (Yelmo.jl, Kryonomos.jl) can drive the same
  components later.
- **Tooling** (grid generation, forcing builder, regression diffs against
  yelmox) in Julia.

## 10. What we keep, and why

| From | Pattern | Why it earned its place |
|---|---|---|
| yelmox | Hub as geometry truth, updated with anomalies | keeps hi-res detail; decouples component grids from the ice grid; makes the ice model one producer among others |
| yelmox | Only anomalies cross the isostasy boundary, with a reference check | smooth fields survive remapping; catches restart mismatches |
| yelmox | Grids as strings, lazy method-typed maps cached to disk | moving a component to another grid is a par-file change |
| yelmox | Neutral output type + capability flags + init-time contracts (`climate_out_class`) | consumers are backend-agnostic; bad combinations fail at init, named |
| yelmox | Run-level vs domain-level ownership (`bsl`, coefficients, OBM) | first explicit regional/global split; becomes `globe` vs `domain` |
| yelmox | Coefficient registry: named scalars, derived values, restart, kill switch | the 0-D case of the exchange; becomes the `forcing` component |
| yelmox | Forcing as reference + coefficient × pattern; one forcing time | one vocabulary replaced six mechanisms; killed the two-clock bug class |
| yelmox | Restart = restore, never step; bundle folder per time | bit-exact restarts; chunked runs |
| yelmox | Per-module output on the native grid; same-interface stubs; configme manifest; `defaults_file` + `nml_validate` + par library | |
| yelmox | Bit-identity gates on Levante for every refactor phase | the working method; kryos keeps a regression harness against yelmox |
| CLIMBER-X | Neutral ice exchange type with canonical units, converted once in the adapter | consumers never change when a backend is added |
| CLIMBER-X | Domain = named grid from a registry + bundle of regional components on it | domains are data, not code |
| CLIMBER-X | Background + regional overlay, coverage masks deciding regional vs global fallback | the cleanest existing answer to the regional/global split |
| CLIMBER-X | Global solid Earth with a narrow contract (`h_ice` in; `rsl`, `z_bed` out), model owns its exotic grid and remap (VILMA2) | host never learns spectral details; backends interchangeable |
| CLIMBER-X | Regional SMB/BMB can run on geometry from data when there is no ice model | offline experiments for free |
| CLIMBER-X | Physical-constant comparison at coupling time | silent constant mismatches are a classic drift source |

## 11. What we shed, and why

| Legacy | Why it goes |
|---|---|
| Ice model at the centre: `kryos_domain%yelmo`, `couple_*_to_yelmo`, `wp` from Yelmo, dx in "Yelmo units", every kryos object compiled against Yelmo | the ice sheet becomes one component; a neutral interface is the precondition for CISM/PISM and for data ice sheets |
| All backends as fields of one type; string dispatch in 4–6 places per component | one polymorphic slot per component; a backend touches one file |
| Cross-component reach-ins (climate takes `marshelf_class`; steps write `mshlf%now`; smb_simple reads Yelmo's `H_ice_ref`; isostasy reads Yelmo not the hub) | hidden edges in the dependency graph; all go through the hub and the exchange |
| Framework-written component output (`isos_write_step`, surface writers) | the instance writes its exports; no framework edit per backend |
| Two drivers (`yelmox`, `yelmox_bipolar`) with diverging loops, restart layouts and `restore_bsl` | one driver for N domains; the OBM bridge becomes a global `ocean` instance |
| Experiment logic in the framework (NEGIS, glacial SMB scaling, `lim_pd_ice`, marine-ice start, recon, relax/opt with hard-coded numbers) | plugins with explicit backend requirements, outside the core |
| `domain_ctl` grab bag; required-key groups without defaults; integer method codes; suffix-based multi-domain par files | per-instance groups with schemas |
| Float cadence arithmetic; cadences inside components (REMBO, chion annual per call); lag implied by call order | integer clock, per-instance cadence, declared lag |
| Sequential mutation of the shared `bsl` by each domain | order-independent global reduction |
| Hard-coded paths (`maps/`, `input/yelmo_phys_const.nml`, output names); per-library `parse_path` copies; per-library precision kinds with ad-hoc conversion | one path/keyword service; `dp` exchange with conversion in adapters |
| CLIMBER-X: flat 200-field hub + 4600-line pairwise coupler with physics inside; link-time dummy modules with drifting type copies; global mutable control/timer flags used inside backends; hard-coded geography | the reasons kryos uses typed per-instance states, a runtime `data` backend and a passed-in clock |

## 12. Additions to fesm-utils

fesm-utils stays a library of mechanisms with no knowledge of specific
components. Candidates, in the order to do them:

| # | Addition | Replaces / why |
|---|---|---|
| 1 | `fesm_log` / `fesm_error(proc, msg, detail [, stat])` + `to_str` | ten private copies; the team's own "fifth module" rule is exceeded |
| 2 | `fesm_fields`: `field_t`, `field_state_t` (`dp` buffers, units, grid tag, validity, time semantics, conservation flag, handles) | nothing exists; the exchange of section 5 |
| 3 | `coupler` extensions: coverage fraction and full-target-area normalisation in `map_field`, `coupler_merge`/`coupler_split`/`coupler_integral`, Gauss cell bounds from Gauss weights, grid hash in the cache key, identity shortcut, mask/fill pass-through, no fixed `GRID_MAX`/`MAP_MAX` | the global ↔ regional gaps of `study/fesm-utils.md` section 2.7; removes VILMA2's rescale |
| 4 | `fesm_clock`: integer-tick `clock_t` + `cadence_t` with a convenient, well-documented interface (decision 8); `tstep_class` and `timeout` become views on it | section 7 |
| 5 | `datasrc_t` on `varslice`; varslice reads the time axis and units from the file | section 8 |
| 6 | restart bundle helper (paths + manifest) | yelmox convention made reusable |
| 7 | `outfile_t`: `variable_io` + `timeout` + `ncio` driven by a field state | per-instance output |
| 8 | `htopo` in `dp`, mass-conserving `H_ice` | hub precision |
| 9 | Prefixed module names (`fesm_*`) with aliases for the old ones | `precision`, `constants`, `timer`, `coupler` collide with CLIMBER-X and FastIsostasy today, hidden by `-zmuldefs` |
| 10 | `make check` + CI | before anything downstream depends on new pieces |

Kept in kryos until a second host wants them: `component_t`, the component
vocabularies, the wiring tables, the domain/globe topology and the sequence.
Adoption work, not new code: `phys_constants` in VILMA2 and CLIMBER-X.

## 13. Repository layout

```
kryos/
  src/core/        kryos_component, kryos_fields_vocab, kryos_domain, kryos_globe, kryos_run,
                   kryos_wiring, kryos_scheduler, kryos_restart, kryos_output
  src/components/  one module per component: abstract type + handles + data backend + factory
                   ice/  geo/  climate/  surface/  shelf/  ocean/  geothermal/  sediments/  forcing/
  src/backends/    <component>/<backend>.f90 and <backend>_stub.f90
                   ice/yelmo, ice/cism, ice/pism (C shim), geo/fastisostasy,
                   geo/vilma2, surface/chion, surface/smbpal, surface/smb_simple, climate/snapshots,
                   climate/esm, climate/rembo, shelf/marine_shelf, ocean/obm
  src/experiments/ opt_cb_ref, relax, negis, glacial_smb, marine_ice_start, recon_start
  drivers/kryos.f90
  par/  input/<backend>_defaults.nml  maps/  tests/  tools/ (Julia)  docs/
```

## 14. Migration

Each phase is gated by bit-identity (or documented differences) against yelmox
on the Levante test cases.

0. fesm-utils additions 1–5, 9, 10 (independent of kryos).
1. kryos core + `ice` (yelmo), `geo` (fastisostasy), `climate`
   (snapshots), `surface` (chion), `shelf` (marine_shelf), geothermal and
   sediments as data; one domain; the hub; output and restart. Target:
   reproduce `yelmox_Greenland` and `yelmox_Antarctica`.
2. `data` backends for every component; esm and rembo climate; smbpal and
   smb_simple; experiments plugins. Target: all yelmox base configs.
3. Globe: `forcing`, N domains, `ocean = obm`. Target: `yelmox_bipolar`.
4. `geo = vilma2` (VILMA2's own gaps: export `u` and `N` separately,
   barystatic definition, `libvilma.a`, FesmData viscosity format).
5. `ice = cism` (serial), then `pism` (C shim).
6. yelmox frozen; CLIMBER-X may adopt `ice` + `geo` later.

## 15. Decisions

All settled 2026-10-24 with Alex. Each entry: the decision, why, what was
rejected.

1. **Backend mechanism**: decided (2026-10-24), abstract types + factory. One
   abstract type per component extends `component_t` with the deferred
   contract and the field handles; each backend is one file with one extension
   type that `use`s the contract and its own library; a factory allocates the
   dynamic type from the backend string once, at init (`allocate(yelmo_ice_t
   :: ice)`), so no procedure branches on a backend name. An unlinked optional
   library is a stub extension whose `init` stops with "rebuild with
   `cism=1`". A host may register its own extension (CLIMBER-X: SICOPOLIS)
   without editing kryos. Experiment plugins reach backend internals only via
   `select type` and stop on a mismatch. Rules to contain compiler risk
   (gfortran ≥ 9, ifort ≥ 19, ifx): no polymorphic arrays (wrapper type for
   collections), no `class(*)` in hot paths, no `final` on types with
   allocatable polymorphic components. Rejected: `select case` + stubs (the
   yelmox climate pattern): the branch repeats in every procedure, every
   backend is always compiled in, and a host cannot add a backend. A
   string-keyed C API (Yelmo C API style) for PISM and a Julia host is
   orthogonal and sits on top.
2. **Exchange precision**: decided (2026-10-24), `dp` buffers; each backend
   converts at its own boundary and keeps its internal `wp` (half the stack is
   already `dp`; volume and area sums over 1e6 cells need it). Rejected: a
   kryos `wp = sp` matching Yelmo.
3. **Mass-flux units**: decided, `kg m-2 yr-1` (section 5). Alternatives
   considered: `kg m-2 s-1` (CF/ISMIP, chion and PISM native; rejected, see
   section 5) and `m ice yr-1` (Yelmo's; needs one agreed `rho_ice` in the
   exchange).
4. **Bed and sea level**: decided (2026-10-24). `geo` exports anomalies
   `dz_bed`, `dz_ss` relative to the reference state it was initialised with,
   plus `bsl`. The domain owns the reference topography: on the hub, and on
   every component grid as a copy made once at init by the same conservative
   map. Each consumer forms full fields on its own grid,
   `z_bed = z_bed_ref(grid) + remap(dz_bed)` and `z_sl = bsl + remap(dz_ss)`,
   all in one absolute datum (the geoid of the reference topography). `geo` is
   initialised from the domain reference, so the yelmox reference check becomes
   a check that `geo` and hub agree. Yelmo's internal `z_bed_ref` is set from
   the same reference at init, so phase 1 is bit-identical; CISM, PISM and
   `data` need no reference concept. `ice → geo` is total `H_ice` (both
   backends apply flotation); Yelmo's `dzbdt_corr` becomes an optional export
   that only FastIsostasy declares it accepts. Under `vilma2` the anomalies
   arrive on the global grid and are split onto each hub; under `geo = data` a
   bed or rsl history is differenced against the reference so it enters as
   anomalies too. Why anomalies: displacements are long-wavelength and survive
   a coarse geo grid; absolute beds would smooth troughs and sills. Rejected:
   anomalies into each backend's own reference (yelmox today; needs a
   reference in every backend and a runtime equality check), and the
   CLIMBER-X datum (bed relative to local sea surface, `z_sl = 0`: loses
   absolute elevations by metres, ties the meaning of `z_bed` to geo
   internals, no place for a sea-level curve).
5. **Solid earth and sea level**: decided (2026-10-24), one component `geo`
   (named after CLIMBER-X's) with two levels: a per-domain part (`dz_bed`,
   `dz_ss` on a domain grid, or on the global grid for VILMA2) and a global
   part (global load, the solve, `bsl`), stepped gather → update → scatter.
   FastIsostasy's global part is only the order-independent barystatic
   reduction; VILMA2's domain part is a copy. The `bsl` source is a component
   option (`own | file | const`), so bed from FastIsostasy with sea level from
   a curve stays expressible. Why one component: the sea-level equation
   produces `u`, `N` and `bsl` from one solve, so two components would have
   to share one backend state. Rejected: one flat component (a regional
   instance cannot own a global scalar; today's driver-owned `bsl_class`), and
   separate `isostasy` and `sealevel` components. `geo` is the only two-level
   component in the design.
6. **Shelf melt**: decided (2026-10-24), always a separate `shelf`
   component; every ice backend takes `bmb_shlf` as an external field (Yelmo
   `bnd%bmb_shlf`, CISM `whichbmlt_float = 4`, PISM external ocean model), and
   the internal schemes of CISM and PISM are not used. The `shelf` component
   is responsible for obtaining `bmb`: kryos owns the ocean → shelf-base
   interpolation (profiles on `depth`, or 2D shelf-base values from a
   `climate`/`ocean` with `has_ocn_shelf`, which also absorbs the ESM
   backend's current interpolation and removes that reach-in); the melt-law
   backend (`marine_shelf` laws, `pico`, `data`) receives shelf-base `T`, `S`,
   `tf` and geometry and returns `bmb_shlf`. `tf_corr` becomes an experiment
   plugin on the `marine_shelf` backend. Why: one path, and the same melt law
   under every ice backend, which an intercomparison needs. An `own_bmb`
   capability on `ice` (route ocean fields to the ice backend, skip `shelf`)
   was considered and deferred; it can be added later as one more capability
   and wiring row without touching the above.
7. **Wiring**: decided (2026-10-24), fixed tables in code per component pair
   (`wire(src, field, dst, field, lag)`), executed by a generic `exchange`;
   the table is logged at init and init stops with a named gap when a
   required import has no producer. Why: every legitimate rewiring is already
   an instance or capability choice (`surface = "climate"`, `has_ocn_shelf`,
   `geo = "data"`), so a par-file connection schema would add only the ability
   to make wrong runs silently. Non-standard edges (NEGIS, glacial SMB
   scaling) are backend-specific plugins, not connections. The table is
   already data; it can be exposed later if a real need appears.
8. **Time base**: decided (2026-10-24), integer ticks. An `int64` count of a
   fixed run-wide resolution (`tick`, a fraction 1/N of a year chosen so that
   1, 0.1, 0.01, 1/12 and 1/360 yr are whole numbers) from an absolute origin.
   Every `dt`, output time and restart time is converted at init and must be
   an exact multiple, or init stops with the offending value. "Due" is an
   integer modulo against absolute time, so chunks with different `dtt` align
   and a restart is valid exactly when every cadence divides the restart tick.
   Years are `dp` at every interface, computed as `n*tick`, never
   accumulated; backends never see ticks. Calendar day-of-year for sub-annual
   components comes from the run's declared calendar, not from the tick.
   Why: removes the tolerance and modulo bug class of the v2.4 sweep
   (`cadence_due` on `nint(time*100)`, `timeout` 1e-6 matching, chunk
   alignment rules). Rejected: float years with Kahan and tolerances (status
   quo); a step counter `time_init + n*dtt` (exact but breaks when `dtt`
   changes between chunks). Built as a clean fesm-utils module (`fesm_clock`,
   row 4 of section 12) with a convenient interface: construct from years and
   a resolution, `due(cadence)`, `interval(cadence)` in years, `is_output_time`,
   `is_restart_time`, `to_years`/`from_years`, restart of its state, and a
   readable log of the schedule at init. `tstep_class` becomes a view on it.
9. **Where `component_t` lives**: decided (2026-10-24), kryos first
   (`src/core`). Its mechanisms (`fesm_fields`, `fesm_clock`, `fesm_log`) go
   to fesm-utils now; the per-component abstract types carry field
   vocabularies and stay in kryos in any case. Why: one consumer until
   CLIMBER-X adopts `ice` + `geo`, and it will then depend on `libkryos`
   for those contracts anyway; the lifecycle contract will change during
   phases 1–3 and must not require a fesm-utils release and pin bumps across
   the stack each time; fesm-utils keeps its mechanisms-only principle.
   Revisit at phase 3 and promote (one module move) if a second host wants
   it. Rejected for now: `fesm_component` from the start, which would freeze
   the contract before any backend implements it.
10. **Driver**: decided (2026-10-24), one program. The program is a few
    lines (`run_init`, loop `run_step`, `run_end`); the coupling sequence is
    written out once in the library module `kryos_run`, looping over
    domains at each step, so tests and a future Julia or C host call the same
    sequence. Global instances own their per-domain ports (the obm instance
    reads `obm_box = "north"` from each domain's group), which removes the
    reason `yelmox_bipolar` exists as a separate driver. A single domain is
    N = 1 with an empty globe; offline SMB, geo-only and climate-only runs
    are `data` backends, not drivers. Rule: no experiment conditionals in
    `run_step`; experiments are plugins at fixed hook points. Rejected: one
    program per configuration, which in yelmox produced two drivers that
    drifted (calendar, CMIP output, restart layout, `restore_bsl`, hand-
    interleaved order).
11. **CISM's scalar sea level**: decided (2026-10-24), fold now. The CISM
    adapter sets `eus = bsl` and `topg = z_bed − dz_ss`, so `eus − topg =
    z_sl − z_bed` and flotation, the marine-connection mask and the calving
    tests are exact; grounded elevations inside CISM are shifted by `dz_ss`
    (metres to tens of metres), whose gradient is negligible for the driving
    stress. The adapter un-folds on export (`z_bed = topg + dz_ss`, `z_srf =
    usrf + dz_ss`), so kryos never sees the shift; only CISM's own output
    files do. Init logs that CISM runs with folded sea level and each geo
    update logs the maximum `|dz_ss|`. A 2D sea level is proposed upstream;
    when accepted, the fold is removed from one adapter. Also required of the
    adapter: trigger the `topg` halo update and `marine_connection_mask`
    recompute after each external change (they only run inside CISM's
    isostasy branch). Rejected: waiting for the upstream change, which would
    hold phase 5 on a CISM release cycle for a first use (comparison with
    Yelmo) where the fold is adequate.
12. **rembo1**: decided (2026-10-24), accepted as a special case and not
    refactored now: one instance per executable, on the Greenland grid it was
    compiled for, behind the stub pattern (`rembo=1`). The constraint is
    handled by a **general mechanism**, not by rembo-specific code: every
    backend declares its constraints at registration (`max_instances`, a
    fixed grid name or none, required domain, no restart, ...), and the
    factory enforces them at init with a named stop ("backend `rembo` allows
    one instance; a second was requested by domain `Antarctica`", "backend
    `rembo` runs only on grid `GRL-16KM`"). Any future model that can only
    work one way (fixed grid, single instance, one hemisphere) uses the same
    declaration. A refactor of rembo1 into an instantiable type with a runtime
    grid stays possible later and would only relax its declared constraints.
13. **Names**: decided (2026-10-24). `globe` for the global scope (`earth`
    collides with solid earth, `global` is an adjective, `run` is the whole
    simulation); `instance` for a configured component (component + backend +
    grid + scope, one par-file group). Component names follow CLIMBER-X
    where possible: `ice` (an ice sheet or a glacier; not `ice_sheet`),
    `geo`, `climate`, `surface` (not `smb`: it exports `T_srf` and may be a
    snowpack), `shelf`, `ocean`, `geothermal`, `sediments`, `forcing`.
    Backends carry their library's name (`yelmo`, `cism`, `pism`,
    `fastisostasy`, `vilma2`, `chion`, `smbpal`, `smb_simple`, `rembo`,
    `snapshots`, `esm`, `marine_shelf`, `pico`, `obm`, `data`). Module
    prefixes `kryos_` and `fesm_`, with the old bare fesm-utils names kept as
    deprecated aliases. New types end in `_t`; existing library `_class`
    types stay. Lifecycle: `init`, `init_state`, `update`, `restart_write`,
    `restart_read`, `write_output`, `end`. Files: `<instance>.nc`,
    `<instance>_ts.nc`, `restart-<time>/<instance>_restart.nc`,
    `input/<backend>_defaults.nml`.
