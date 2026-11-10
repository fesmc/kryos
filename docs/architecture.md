---
title: "Architecture"
---

How a kryos run is put together: the components and their backends, the instances named in the par file, the
domains and the globe, the fields they exchange, the clock and the coupling sequence, and the output and restart
files. How several regional domains share one global solid Earth is the subject of
[Global hub and regional domains](global-hub.md).

The design notes behind this page are in the repository: [`docs/design.md`](https://github.com/fesmc/kryos/blob/main/docs/design.md), [`docs/dev/core-design.md`](https://github.com/fesmc/kryos/blob/main/docs/dev/core-design.md) and
[`docs/dev/contracts-design.md`](https://github.com/fesmc/kryos/blob/main/docs/dev/contracts-design.md) (see [Developer notes](developer-notes.md)). Where they differ from the code, the code
wins; this page follows the code.

## Overview

```{mermaid}
flowchart TB
    run["run (kryos_run)<br/>clock, schedules, output, restart"]
    globe["globe (kryos_globe)<br/>forcing, geo (global part),<br/>global ice field (kryos_globe_ice)"]
    d1["domain North<br/>hub NH-32KM, coupler,<br/>ice_nor, geo_nor"]
    d2["domain Antarctica<br/>hub ANT-16KM, coupler,<br/>ice_ant, geo_ant"]
    d3["domain ...<br/>hub, coupler,<br/>ice, geo, climate, surface, shelf, ..."]
    run --> globe
    run --> d1
    run --> d2
    run --> d3
    globe <-->|"gather / scatter<br/>(wiring rows, globe_ice)"| d1
    globe <--> d2
    globe <--> d3
```

- A **run** owns one clock, the output schedules, the restart bundles and the coupling sequence
  ([`kryos_run.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_run.f90),
  [`kryos_run_sequence.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_run_sequence.f90)).
- The **globe** holds the instances that belong to the run, not to a domain: the `forcing` (run-level scalars) and
  the global part of `geo` ([`kryos_globe.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_globe.f90)).
  A global geo resolved on the host's loads (VILMA2) adds the global grids, their coupler and the global ice field
  ([`kryos_globe_ice.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_globe_ice.f90)).
- A **domain** is a hub grid with its own component instances, its coupler (grids and map cache) and the reference
  topography on every grid in use
  ([`kryos_domain.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_domain.f90)).
- Instances exchange **fields** only through a fixed **wiring table**
  ([`kryos_wiring.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_wiring.f90)): source field, target
  field, remap method, lag and transform, in code, not in the par file.

The dependency flow is one-directional: fesm-utils → core configuration and `component_t` → components and their
backends (factories as submodules) → wiring → domain, globe → output, restart → run. The core never `use`s a backend
module; a backend never sees the run, the domain or the globe.

## Components and backends

A **component** is a role with a fixed vocabulary of imports and exports (`src/components/<comp>/`). A **backend** is
one implementation of it, chosen at run time by the `backend` key of the instance group; the factory of the component
(`kryos_<comp>_factory.f90`) allocates it. An optional library that is not built is a stub with the same names whose
`init` stops with "rebuild with `<flag>=1`".

| Component | Role | Backends (factory) | Scope |
|---|---|---|---|
| `ice` | ice sheet: thickness, surface, grounding; mass-balance diagnostics | `yelmo`, `plastic` ([PlasticIceSheet](backend-plastic.md)) | regional |
| `geo` | solid Earth and sea level: bed and sea-surface change | `fastisostasy` (regional, plus its global part for `bsl`), `vilma2` ([VILMA2](backend-vilma2.md): global Earth, regional ports) | regional and global |
| `climate` | atmosphere and ocean fields | `snapshots` | regional |
| `surface` | surface mass balance and surface temperature | `chion`, `smb_simple` ([smb_simple](backend-smb-simple.md)) | regional |
| `shelf` | sub-shelf melt | `marine_shelf` | regional |
| `geothermal` | geothermal heat flux | `data` | regional |
| `sediments` | sediment thickness | `data` | regional |
| `forcing` | named run-level scalars (coefficients) | `tsgen` | global |

Scope is a property of the instance, not of the component: the geo factory picks the global VILMA2 Earth
(`vilma_earth_t`) for scope `global` and the regional port (`vilma_geo_t`) for scope `regional`
([`kryos_geo_vilma2.f90`](https://github.com/fesmc/kryos/blob/main/src/backends/geo/kryos_geo_vilma2.f90)).

Every backend extends `component_t`
([`kryos_component.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_component.f90)) through its component's
abstract type. The bindings the core calls:

| Binding | When |
|---|---|
| `init(ctx)` | once; the init context carries the instance's grid, the domain's regions on it, the reference topography on it, the run's physical constants and calendar year, and the axes of the instances made before |
| `init_state(st, pass)` | the cold-start passes `PASS_INIT`, `PASS_GEOMETRY`, `PASS_BOUNDARY`, `PASS_STATE`, `PASS_SPINUP`, `PASS_FINAL`, and `PASS_RESTART` after a restore |
| `prepare(st)` | before `update` when due: imports into the backend's own boundary (Yelmo), so that an experiment hook can modify it |
| `update(st)` | integrate `t0 → t1` and put every export |
| `restart_write`, `restart_read`, `restart_files` | the backend's own restart files (none for a stateless backend) |
| `write_output` | the default writer (exports on the native grid), plus backend streams |
| `flotation` | ice only, optional (capability `has_flotation`): the ice model's grounding rule, used by the global gather |

Experiments (`src/experiments/`: `relax`, `opt_cb_ref`, `opt_tf`, `recon`) are plugins between the domain and the
run, called at four fixed hooks (`HOOK_STEP_BEGIN`, `HOOK_ICE_BND`, `HOOK_COLD_PRE_STATE`, `HOOK_COLD_POST_STATE`).
A domain lists them in its `experiments` key.

## Instances and par groups

An **instance** is a component with a backend, a grid and a scope. It is named in the par file, and its name is its
namelist group, the prefix of its output files and the name of its restart files. Instance names are unique in a run
(`run_init` stops on a duplicate).

The groups of a run (example: [`par/kryos_NH-ANT-PAT_plastic_vilma2.nml`](https://github.com/fesmc/kryos/blob/main/par/kryos_NH-ANT-PAT_plastic_vilma2.nml)):

| Group | Template | Names |
|---|---|---|
| `[kryos]` | `[kryos]` | the domains, the globe group, the clock group, the restart bundle and its schedule, `out_dir`, `map_fldr`, the physical constants, library files |
| `[clock]` | `[clock]` | `time_init`, `time_end`, `dtt`, `method`, calendar |
| `[out_1D]`, `[out_2D]`, `[out_rst]`, ... | `[schedule]` | output and restart schedules |
| `[globe]` | `[globe]` | the `forcing` instance and the global `geo` instance |
| `[North]`, `[Antarctica]`, ... | `[domain]` | `grid_hub`, the `topo` and `regions` groups, one instance per component (`""` = absent), `experiments` |
| `[ice_nor]`, `[geo_glob]`, ... | `[instance]` + `[<component>_<backend>]` | `backend`, `grid`, `dt`, `phase`, `clock`, `out_*`, plus the backend's keys |
| `[vilma]`, `[globe_ice]`, `[topo_nor]`, `[X_z_bed]`, ... | library or datasrc schemas | backend physics groups and data field groups, named by a key of an instance or of the globe |

A component runs in a domain if and only if the domain names an instance for it.

```fortran
&globe
    geo = "geo_glob"
/
&North
    grid_hub = "NH-32KM"
    topo     = "topo_nor"
    regions  = "regions_nor"
    ice      = "ice_nor"
    geo      = "geo_nor"
/
&ice_nor
    backend = "plastic"
    grid    = ""          ! the hub grid, NH-32KM
    dt      = 0.0         ! every step
    clock   = "model"
    out_2D  = "out_2D"
    ...                   ! the keys of [ice_plastic]
/
```

## The par-file template

[`input/kryos_defaults.nml`](https://github.com/fesmc/kryos/blob/main/input/kryos_defaults.nml) holds the schema of
every group: one template per group (`[kryos]`, `[clock]`, `[schedule]`, `[globe]`, `[domain]`, `[instance]`,
`[experiment]`) and one per backend (`[ice_plastic]`, `[geo_vilma2]`, `[surface_smb_simple]`, ...), plus the groups of
code kryos hosts itself (`[orbit]`, `[globe_ice]`, `[smb_simple]`). Every group of the par file is validated against
its template (`config_validate`); an instance group against `[instance]` and `[<component>_<backend>]` together (a key
is valid if either declares it). The template also gives the defaults of the keys a group leaves out. Keys marked
"required" must be set.

Backend physics groups of external libraries (`[yelmo]`, `[isos]`, `[chion]`, `[vilma]`, ...) keep their own schemas
in `input/<library>_defaults.nml`; data field groups follow `input/datasrc_defaults.nml`. Instance, domain and run
groups must be in the par file; named library groups may come from the files of `[kryos] library`. All keys are listed
in [Par files](reference-par-files.md).

## Domains

A domain is a regional scope: a named hub grid, the hub, per-component grids, its regional instances and its own
coupler.

- **The hub** (fesm-utils `htopo`) holds the prepared reference topography (`[topo]`: `z_bed`, `z_srf`, `H_ice` from
  the file of `path`, with `z_bed_f_sd`, `remove_lakes`, `H_ice_min`, `grad_lim_zb`), the regions, and the current
  geometry. After every ice step it follows the ice (`domain_hub_update`): on the ice grid it mirrors the ice; on a
  finer hub it adds the ice's anomalies of `z_bed` and `H_ice` (bilinear) to its own reference and recomputes the
  grounding and surface. Climate, surface and shelf read the geometry from the hub. The hub is never stored in a
  restart: it is derived from the ice.
- **Grids.** Every instance runs on its own grid (`grid` key); a blank `grid` follows the default chain
  (`GRID_DEFAULT`): ice ← hub, geo ← ice, climate ← ice, surface ← climate, shelf ← hub, geothermal and sediments ←
  ice.
- **The coupler** (fesm-utils `coupler`) resolves grids by name from `grid_<name>.txt` in the map folder
  (`[kryos] map_fldr`, default `maps`) and builds each map once, caching it in that folder. A remap between two
  identical grids is a copy.
- **References.** The hub prepares the reference topography once; the reference on every instance grid in use is
  made from it by a conservative remap and kept in `dom%refs`. Every component takes its reference from there: the
  ice's `z_bed_ref`, geo's `z_bed_ref`. `domain_check_reference` compares geo's and the ice's `z_bed_ref` exports with
  the domain's reference and stops above 1 m (after geo's `PASS_GEOMETRY` and after every restore).

Instances are made in `INIT_ORDER` (ice, geo, geothermal, sediments, climate, surface, shelf: a producer of an axis
before its consumers) and restored in `RESTORE_ORDER` (ice, geo, climate, surface, shelf, geothermal, sediments).

## Fields: imports and exports

Each instance declares an import and an export field state (`fesm_fields`) at `init`. A field has a name, canonical
units, a grid, a rank (0-D, 2-D or 3-D with an axis such as `month` or `depth`), a validity interval `[t0, t1]`, a
time semantics and a conservation flag. Values are `dp` in the field states; each instance declares the precision of
its own arrays (`prec`), and a consumer narrows once.

| Time semantics | Meaning |
|---|---|
| `INSTANT` | a state at `t1` (geometry, bed, sea level, temperatures) |
| `MEAN` | the mean over the producer's interval `[t0, t1]` (rates: `dHidt`, `smb`, `bmb`, `cmb`, `mb_net` of the ice) |

A `MEAN` export is finalised at `t1` after the backend's call (`component_close`). A row from a `MEAN` export into an
`INSTANT` import needs equal cadences on both sides (checked by `wire_resolve`). At the cold start a `MEAN` export
put in `init_state` covers a virtual previous step `[t0 − dt, t0]`, so that lag-1 consumers find `t1 = t0` at the
first step.

- **Optional fields.** A backend declares an optional export or import only if it supports it (`self%opt%<name>`,
  capabilities `self%caps`); a row that touches an undeclared field is inactive.
- **Optional imports.** A backend may mark an import optional (`component_optional`): it may stay unwired. The ice
  marks `T_shlf`, `tf_shlf`, `Q_geo`, `H_sed`, `Qd` optional when declared; the plastic backend also `z_bed`, `z_sl`,
  `T_srf`, `bmb_shlf` and (except with `drive = "smb"` or the inverted mass balance) `smb`, so that a domain with
  only ice and geo is valid. Until `z_bed` and `z_sl` are published it uses the reference bed and `z_sl = 0`.
- **Core-filled imports.** A global geo resolved on the host's loads marks the imports the core fills outside the
  rows (`component_core_filled`: the external ice and the background reference, filled by `gia_gather` and
  `globe_ice`).
- **Diagnostics.** An `INSTANT` export declared with `restart = .false.` is not stored in a restart bundle; its
  backend puts it again at `PASS_RESTART`. No active row may read one.

## Wiring

The wiring table (`wiring_table` in
[`kryos_wiring.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_wiring.f90)) is fixed in code. A row is
`src.field → dst.field`, a remap method (`con`, `bilin`, `auto`, or `""` from the field's conservation flag), a lag
and a transform:

| Transform | Result |
|---|---|
| none | `dst = remap(src)` |
| `T_ADD_REF` | `dst = ref_dst + remap(src)`: `z_bed` of the ice = its reference bed + geo's `dz_bed` |
| `T_SUB_REF` | `dst = remap(src − ref_src)`: the ice's anomalies onto a finer hub |
| `T_ADD_SCALAR` | `dst = src_0D + remap(src)`: `z_sl` of the ice = `dphi` + geo's `dz_ss` |
| `T_LAPSE_T`, `T_LAPSE_P` | `tas`, `pr` re-elevated to the consumer's surface |
| `T_GIA` | the core's gather into a host-resolved global geo; resolved and lag-checked, never run by the exchange |

| Lag | The consumer in step `t_{n−1} → t_n` sees |
|---|---|
| `0` | the producer's state of this step (geo → ice `z_bed`, `z_sl`; hub → climate, surface, shelf; global geo → regional geo) |
| `1` | the producer's state of the previous step (ice → geo `H_ice`; surface → ice `smb`, `T_srf`; shelf → ice; ice → forcing `M_ice`; ice → global geo, the gather) |
| `LAG_INIT` | exchanged once after `PASS_INIT` (static fields: the ice's `z_bed_ref`, `H_ice_ref` → geo; `Q_geo`, `H_sed` → ice; the hub's `z_srf_ref` → climate) |

At init `wire_resolve` binds the rows of each domain (with the globe's nodes): a row is active if both fields are
declared. It stops, naming the field and the row, on a required import without an active row (a gap), on two active
rows into one import, on a link mismatch (units, rank, axis), on a masked field across two grids, on `MEAN → INSTANT`
with unequal cadences, and on a diagnostic that an active row reads. Every active pair of grids is primed in the
coupler. One line per row is logged.

At every exchange the lag is checked: the published `t1` must equal the producer's last update, and be at most `t`
(lag 0) or before `t` (lag 1). A sequence edit that swapped two calls therefore stops the run instead of changing the
result. The rows run where the sequence calls `exchange(dom, src, dst, t)` (one source and destination, in table
order) or `exchange_into(dom, dst, t)` (every row into `dst`, sources in the order geo, surface, shelf, climate, then
the static ones).

## Clock and cadences

One clock (`fesm_clock`) counts integer ticks (`[clock] ticks_per_year`, default 36000) and reports years in `dp`.
`[clock]` gives `time_init`, `time_end`, the run step `dtt`, and the forcing time (`method`: `rel`, `cal` or `const`;
`time_ref`, `calendar`, `year`).

Each instance has a **cadence** from its keys:

| Key | Meaning |
|---|---|
| `dt` | `0`: every step; `> 0`: on that cadence grid [yr]; `< 0`: never after init (static data) |
| `phase` | offset of the cadence grid [yr] |
| `clock` | `model` (the model time) or `forcing` (`t_eval` is the forcing time `time_rel` at `t1`) |

An instance is called with `step_time_t` `(t0, t1, t_eval, dtt, tick)`: its own interval since its last update.
`clock_start` stops at init if a restart time or `time_end` is off any cadence grid. A regional VILMA2 port takes the
cadence of the global Earth (see [Global hub and regional domains](global-hub.md)).

## The run sequence

```fortran
call run_init(run, par_file)            ! config, clock, globe, domains, wiring, schedules; cold start or restore
do while (.not. run_done(run))
    call run_step(run)                  ! advance the clock, the sequence, output and restart at the new time
end do
call run_end(run)
```

### Init (`run_init`)

1. Configuration, log level, lossless switch, physical constants (`[kryos] phys_const` of `phys_const_file`).
2. The clock.
3. The globe: forcing, then the global geo (and, for a host-resolved geo, the global ice field).
4. Per domain: groups, grids, hub, coupler, references; the instances (factory, `init`, cadence); for a
   host-resolved global geo the ports get its cadence and every domain's coupler its grid. Then `gia_check`, the
   wiring of every domain, the experiments.
5. Schedules and output streams.
6. `PASS_INIT` of every instance (globe first, then each domain in `INIT_ORDER`) and the static rows; then the restore
   from `[kryos] restart`, or the cold start.
7. Describe, `clock_start`, output at `time_init`.

### Cold start

Per domain, in the order of yelmox's `domain_init_state`
([`run_cold_start`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_run_sequence.f90)):

- a) bed and sea level for every domain: rows ice → geo, global `bsl` → geo, geo `PASS_GEOMETRY`,
  `domain_check_reference`; then the global reduction of the initial volume change (`reduce_bsl`);
- rows geo → ice, ice `PASS_GEOMETRY`;
- b) boundary fields from the hub: climate, surface and shelf `PASS_BOUNDARY` (each after its rows from the hub and
  the climate), then every row into the ice and the ice's `prepare`;
- c) the ice state: `HOOK_COLD_PRE_STATE`, ice `PASS_STATE`, `HOOK_COLD_POST_STATE`, ice `PASS_SPINUP`;
- with VILMA2: the seed of the Earth under the final t0 load, then c') geo → ice and ice `PASS_GEOMETRY` again;
- d) the hub on the final geometry and one more surface cycle (`PASS_FINAL`), its rows into the ice.

### One step (`run_sequence`)

```{mermaid}
flowchart TB
    s1["1. globe: gather (VILMA2) or bsl_update;<br/>forcing: M_ice gathered, coefficients updated"]
    s2["2. HOOK_STEP_BEGIN, every domain"]
    s3["3. geo, every domain when due:<br/>ice → geo, globe → geo, geo update;<br/>then reduce_bsl"]
    s4["4. per domain: rows into ice, prepare, HOOK_ICE_BND,<br/>ice update, hub update,<br/>climate, surface, shelf when due"]
    s5["5. kill switch, output, restart bundle"]
    s1 --> s2 --> s3 --> s4 --> s5
```

1. **Globe, start of the step.** A host-resolved global geo, when due: `gia_gather` (the ice of the previous step),
   then the Earth's update; FastIsostasy: `bsl_update`. The forcing: `M_ice` of every domain gathered (lag 1), the
   climate's derived coefficients (lag 1), then the forcing update.
2. **Hooks** `HOOK_STEP_BEGIN` of every domain's experiments.
3. **Geo**, every domain on its cadence: rows ice → geo (`H_ice` of the previous step), the globe → geo rows, the geo
   update. Then `reduce_bsl`: the domains' `dV_bsl` gathered in name order, `bsl −= Σ dV / A_ocn`, the reduced `bsl`
   back into each geo (`geo_scatter`). With VILMA2 there is no `dV_bsl` source and it returns.
4. **Per domain**: `exchange_into(ice)` (geo with lag 0, surface and shelf with lag 1), and when the ice is due its
   `prepare`, `HOOK_ICE_BND` and `update`; the hub update; climate (rows from the hub, coefficients from the forcing),
   surface (rows from the hub, the climate and the ice) and shelf (rows from the hub and the climate), each when due.
5. **End of step** (`run_step`): the forcing's kill switch, the due outputs, the restart bundle when its schedule is
   due (always at `time_end`), a forcing-increment bundle when due.

Presence is tested (`domain_present`), never a backend name: a domain without climate, surface or shelf skips those
calls.

### Restore

`run_restore` restores, it never steps: the clock and every cadence (a continuation: the bundle time is `time_init`;
a branch with instances stops), the manifest against the run's instances, the globe (forcing, then geo), then per
domain every instance in `RESTORE_ORDER` (its backend files, its imports and exports as the bundle holds them,
`PASS_RESTART`), `domain_check_reference` and the hub. A regional FastIsostasy geo first gets its rows from the
restored ice and the global `bsl`; no other row runs.

## Output

A stream is an instance's `out_<stream>` key whose value names a schedule group (not `"none"`):

| Key | File | Content |
|---|---|---|
| `out_2D` | `<instance>.nc` | rank-2 and rank-3 exports on the native grid (`xc`, `yc` in km; a global grid on `lon`, `lat`) |
| `out_1D` | `<instance>_ts.nc` | the 0-D exports |
| `out_<stream>` of the backend | `<instance>_<stream>.nc` | backend streams, e.g. `out_plastic`, `out_vilma`, `out_gia`, `out_yelmo` |
| domain `out_2D` | `<domain>_hub.nc` | the hub |
| `[globe_ice] out_2D` | `globe_ice.nc` | the global geometry hub |
| (first 1D schedule) | `timing.txt` | wall time per stage |

`out_vars` restricts the exports of `2D` and `1D` (`""` = all). Files go to `[kryos] out_dir`.

The **variable tables** `input/kryos-variables-<component>.md` (exports) and
`input/kryos-variables-<component>-<backend>.md` (backend streams; `kryos-variables-hub.md`,
`kryos-variables-globe_ice.md`) give per variable its dimensions, units, long name and `prec`. The `prec` column is
the only source of rounding (absolute, `r` relative, `-` lossless). A rank-2 or rank-3 export without a row stops the
writer, and a row whose units differ from the export's stops; 0-D series are lossless. `[kryos] write_lossless =
True` writes every file lossless (the gates set it). Restart files are always lossless.

## Restart bundles

A bundle is a folder per restart time, written when the `[kryos] restart_out` schedule is due and always at
`time_end` ([`kryos_restart.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_restart.f90)):

```
restart-<time>/
  clock_restart.nc                       the clock and the cadence of every instance that steps
  <inst>_imports.nc  <inst>_exports.nc   every global instance (forcing, geo_glob)
  <inst>_restart.nc ...                  a global instance's backend files (e.g. geo_glob_restart.nc)
  <domain>/
    <inst>_imports.nc  <inst>_exports.nc every instance of the domain
    <inst>_*restart*.nc                  its backend files (stateful backends)
  manifest.nml                           [bundle] + one group per instance; written last
```

- The field files hold the state: every import, every export an active row reads, and every `MEAN` export.
  Diagnostics are not stored and are put again at `PASS_RESTART`.
- A backend declares the files its `restart_write` creates (`restart_files`); the writer checks they exist, the
  restore checks that every file the manifest lists exists.
- The manifest lists each instance's scope, domain, component, backend, grid, files and git version. A restore stops
  on an instance missing on either side or with another component, backend or grid; a changed `dt` is allowed.
- The hub and the global ice field are not stored: both are derived.
- A run continues from a bundle with `[kryos] restart = "<run dir>/restart-<time>"` and `[clock] time_init` set to the
  bundle time; a chunk of a long run is the same par file with these keys and `[clock] time_end` overridden (runme
  `-p clock.time_init=... clock.time_end=... kryos.restart=...`).

## Source map

| File | Content |
|---|---|
| [`src/core/kryos_config.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_config.f90) | par file, templates, validation, library lookup, cadence keys |
| [`src/core/kryos_component.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_component.f90) | `component_t`, `init_ctx_t`, `step_time_t`, passes, default writer |
| [`src/core/kryos_wiring.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_wiring.f90) | the wiring table, resolution, exchange, lag check |
| [`src/core/kryos_domain.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_domain.f90) | domains: grids, hub, coupler, references, instance slots, exchanges |
| [`src/core/kryos_globe.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_globe.f90) | the globe: forcing, global geo, global coupler, globe rows |
| [`src/core/kryos_globe_ice.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_globe_ice.f90) | the global ice field and reference |
| [`src/core/kryos_run.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_run.f90) | `run_t`, `run_init`, `run_step`, `run_end` |
| [`src/core/kryos_run_sequence.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_run_sequence.f90) | cold start, restore, one step |
| [`src/core/kryos_run_gia.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_run_gia.f90) | reference, gather and seed of a host-resolved global geo |
| [`src/core/kryos_output.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_output.f90), [`kryos_outfile.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_outfile.f90) | output streams and files |
| [`src/core/kryos_restart.f90`](https://github.com/fesmc/kryos/blob/main/src/core/kryos_restart.f90) | restart bundles and the manifest |
