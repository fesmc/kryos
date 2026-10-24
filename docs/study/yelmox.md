# Study: yelmox (v2.4) as input to kryos

Scope: yelmox `main` at 62bc7a6 (v2.4 + VILMA2 plan), read 2026-10-24. All
paths are relative to `/Users/alrobi001/models/yelmox` unless noted. Line
numbers refer to that state.

The existing in-tree "kryos" (`libs/kryos*.f90`) is a first generalization of
yelmox: one `kryos_domain` bundling Yelmo and its components, a hi-res geometry
hub, a per-domain coupler and one module per concern. It is still
Yelmo-centred: Yelmo is a field of the domain, the coupling routines are named
`couple_*_to_yelmo`, and most framework modules compile against Yelmo.

---

## 1. Architecture as it stands

### 1.1 Layering

```
fesm-utils   coupler (grids by name, map cache), htopo (hub), regions, nml, ncio,
             varslice, tsgen, timestepping, timeout, timer, phys_constants
   |
libs/        component wrappers: marine_shelf(+pico), smbpal/, smb_simple, surface_chion,
             geothermal, sediments, snapshots, esm_forcing, climate_rembo(_stub),
             climate_out, climate_terms, climate_mods, yelmox_climate, ice_bnd, orbit_forcing
   |
libs/kryos*  kryos (domain type, ctl, init, remap), kryos_coupling (step_*/couple_*),
             kryos_startup (cold start, restart bundles), kryos_output,
             kryos_regions (experiment physics), kryos_forcing (coefficient registry)
   |
drivers      yelmox/yelmox.f90 (1 domain), yelmox_bipolar/yelmox_bipolar.f90 (+ obm_coupling.f90)
```

External models: Yelmo (`yelmo/`), FastIsostasy (`FastIsostasy/`), chion,
insol, rembo1 (optional), all as symlinked sibling repos.

### 1.2 `kryos_domain` (libs/kryos.f90:163-184)

| Field | Type | Role |
|---|---|---|
| `yelmo` | `yelmo_class` | ice sheet (concrete type, not abstract) |
| `mshlf` | `marshelf_class` | sub-shelf melt |
| `isos` | `isos_class` (FastIsostasy) | bedrock + sea surface |
| `cl` | `yelmox_climate_class` | climate backend state (snapshots/esm/rembo all embedded) |
| `clim` | `climate_out_class` | neutral climate output (now/ref + capability flags) |
| `smb`, `smbs`, `schn` | smbpal, smb_simple, chion wrappers | three surface backends, **all** present in the type, one used |
| `sed`, `gthrm` | sediments, geothermal | static data, read once |
| `cnst` | `phys_const_class` | physical constants, loaded once, handed to every component |
| `topo` | `htopo_class` | hi-res geometry hub |
| `cpl` | `coupler_class` | this domain's grid table + map cache |
| `rlx`, `opt`, `ngs`, `gsmb` | relax/optimization/NEGIS/glacial-smb params | experiment and spin-up machinery |
| `ibnd` | `ice_bnd_class` | prescribed ice in Yelmo's border cells |
| `ctl` | `domain_ctl` | `[domain]`, `[comps]`, `[sim]`, `[output]` |
| `mask_recon` | logical(:,:) | recon cold start |

`domain_ctl` (kryos.f90:55-128) mixes four concerns: component selection
(`with_*`, `climate`, `surface_method`, `dt_clim`, `dt_isos`), cold-start
recipe (`init_method`, `recon_*`, `init_marine_H`, ...), experiment switches
(`relax`, `opt`, `use_negis`, `scale_glacial_smb`, `lim_pd_ice`), grids
(`grid_hub/ice/isos/clim/surface/shelf` + derived `dx_*`, `south`) and output
switches.

Run-level objects are deliberately **not** in the domain: `bsl_class`
(FastIsostasy barystatic sea level) and `coefs_class` (forcing coefficients)
are owned by the driver (yelmox.f90:38,62; bipolar:64-65), and the OBM by the
bipolar driver.

### 1.3 The hub (`htopo`, fesm-utils/src/htopo.f90)

- Lives on `grid_hub`, the finest grid of the setup. Holds the static
  reference (`ref`: z_bed, H_ice, z_srf, z_bed_sd from FesmData), the regions
  (`reg`), and the *current* geometry `z_bed, H_ice, z_srf, f_grnd, z_sl`.
- It is the source from which every component is populated, Yelmo included:
  `domain_init` remaps `topo%ref` conservatively to `grid_ice` and passes it as
  both Yelmo's initial and present-day topography (kryos.f90:303-313).
- Refreshed after the ice step by `couple_yelmo_to_htopo`
  (kryos_coupling.f90:541-570): identical grid -> mirror Yelmo; finer hub ->
  keep the hi-res reference and add Yelmo's *anomalies* (bed displacement
  `z_bed - z_bed_ref`, ice change from the hub reference as Yelmo received it),
  refined bilinearly; `htopo_update` recomputes flotation and surface on the
  hub.
- Climate, surface and shelf read geometry from the hub only
  (step_climate 430-435, step_surface 480-490, step_shelf 592-596).
  **Isostasy does not**: it reads `yelmo%tpo%now%H_ice` and
  `yelmo%bnd%dzbdt_corr` directly (step_isostasy 146-147).
- The hub is never written to restart; it is rebuilt from the restored models
  (kryos_startup.f90:100, 404-406).

### 1.4 Coupler and grids

- `coupler_class` per domain (`dom%cpl`). Grids are named strings resolved from
  `maps/grid_<name>.txt` (`MAP_FLDR = "maps"`, kryos.f90:50). Maps are
  directed, method-typed edges (`con`, `bilin`, `nn`, ...), built lazily,
  cached in memory and on disk (docs/multigrid.md "Conceptual model").
- `remap(dom, var_src, src, var_dst, dst, method)` (kryos.f90:714-750)
  short-circuits to a copy when `src == dst` and allocates the destination.
- `remap_method_smooth(dx_src, dx_dst)` (687-703): bilinear to refine,
  conservative to coarsen. Spacings are expressed in "Yelmo dx units"
  (`dx_hub = yelmo%grd%G%dx * ratio`, kryos.f90:301, 344-345, 362, 377, 420).
- Grid defaults chain (kryos.f90:574-578): `grid_ice <- grid_hub`,
  `grid_isos <- grid_ice`, `grid_clim <- grid_ice`, `grid_surface <-
  grid_clim`, `grid_shelf <- grid_hub`.
- Regions reach each component on its grid by nearest neighbour
  (`domain_regions_load`, kryos.f90:201-213); each component builds its own
  masks from selection expressions in its own group.

### 1.5 Stepping order (yelmox.f90:176-208)

```
tstep_update(ts, dtt)
bsl_update(bsl, time_rel)                      ! run level
coefs_update(coefs, time_rel, var = V_ice[Gt]) ! run level; V_ice read from dom%yelmo%reg
step_relax          ! Yelmo ytopo params ramp
step_optimize       ! cb_ref (Yelmo grid), tf_corr (shelf grid, in mshlf%now)
step_isostasy       ! dt_isos cadence; H_ice(t_{n-1}) -> isos_update to t_n
couple_to_yelmo     ! isostasy (this step), surface + shelf + Qd (from last step)
step_ice_bnd        ! H_ice_ref in border strips at time_rel
step_icesheet       ! yelmo_update(t_n)
couple_yelmo_to_htopo
step_climate        ! dt_clim cadence, on grid_clim, geometry from hub
step_surface        ! every step, on grid_surface
step_shelf          ! every step, on grid_shelf
[forcing-increment restart] [kill switch]
```

Output/restart are written at the top of the loop for the current time, the
final state always written (yelmox.f90:147-174).

### 1.6 Data flow

```
             run level (driver)
   bsl_class ──bsl_now──────────────┐            coefs_class (scalars by name)
      ▲  (fastiso: isos_update       │             │  ▲ derived: air_anomaly (climate),
      │   writes bsl_now back)       │             │  │          obm_north/south (OBM)
 ────────────────────────────────────┼─────────────┼──┼──────── domain ──────────────
                                     ▼             ▼  │
  Yelmo (grid_ice) ──H_ice, dzbdt_corr──► FastIsostasy (grid_isos)
     ▲    │                     [con/bilin]        │ w+we, z_ss-bsl  (anomalies only)
     │    │                                        ▼
     │    │                          couple_isostasy_to_yelmo: z_bed = z_bed_ref + dz,
     │    │                                                    z_sl = bsl + dz_ss
     │    └─z_bed-z_bed_ref, H_ice-H_ref, z_sl [bilin]──► htopo hub (grid_hub)
     │                                                          │ z_srf,H_ice,z_bed,f_grnd,z_sl,
     │                                                          │ z_srf_ref  [con/bilin]
     │                                                          ▼
     │                       climate backend (grid_clim) ◄── coefs (terms)
     │                         snapshots | esm | rembo     (esm also reads mshlf params)
     │                         -> climate_out: tas,pr (monthly), ta_ann, tsl_ann, pr_ann,
     │                            z_srf, beta_p, lapse(12); to_ann/so_ann(depth) or
     │                            T_shlf/S_shlf/dT/dS; smb,tsrf (has_smb); Qd (has_Qd)
     │            ┌───────────────────────┼──────────────────────────┐
     │            ▼ tas,pr [bilin + lapse]▼ to/so(z) or T/S_shlf     ▼ Qd [con]
     │   surface (grid_surface)     marine_shelf (grid_shelf)        │
     │   chion|smbpal|smb_simple    + hub geometry [bilin]           │
     │   + hub geometry             (smb_simple also reads           │
     │   (smb_simple: H_ice_ref     yelmo H_ice_ref)                 │
     │    from Yelmo grid)                                           │
     │            │ smb, tsrf [con]         │ bmb_shlf, T_shlf [con] │
     └────────────┴─────── couple_to_yelmo (next step) ──────────────┘
       + Q_geo (geothermal), H_sed (sediments): written once at init
       + H_ice_ref in fixed/relax cells (ice_bnd, time_rel)
       + optional: glacial smb scaling (clim ta_ann vs ref), lim_pd_ice (-4 m/yr)
```

Fields crossing boundaries:

| From -> To | Fields | Grid hop | Method | Where |
|---|---|---|---|---|
| Yelmo -> isostasy | `H_ice`, `dzbdt_corr` | ice -> isos | smooth | coupling:146-147 |
| isostasy -> Yelmo | `w+we`, `z_ss-bsl` (anomalies) | isos -> ice | smooth | coupling:222-226 |
| Yelmo -> hub | `z_bed-z_bed_ref`, `H_ice-H_ref`, `z_sl` | ice -> hub | bilin | coupling:564-568 |
| hub -> climate | z_srf, H_ice, z_bed, f_grnd, z_sl, z_srf_ref | hub -> clim | smooth | coupling:430-435 |
| climate -> surface | tas, pr (monthly), z_srf, beta_p, lapse; tsl_ann (smb_simple) | clim -> surface | bilin + lapse | coupling:508-539, 481 |
| hub -> surface | z_srf, H_ice, z_bed, z_sl | hub -> surface | smooth | coupling:487-490 |
| climate -> shelf | to_ann, so_ann, dto (3D) or T/S/dT/dS_shlf (2D) | clim -> shelf | bilin | coupling:598-613 |
| hub -> shelf | H_ice, z_bed, f_grnd, z_sl, z_srf | hub -> shelf | bilin | coupling:592-596 |
| surface -> Yelmo | smb (mm w.e. -> m ie), tsrf | surface -> ice | con | coupling:306-322 |
| shelf -> Yelmo | bmb_shlf, T_shlf | shelf -> ice | con | coupling:353-356 |
| climate -> Yelmo | Qd | clim -> ice | con | coupling:182-183 |
| climate (smb_direct) -> Yelmo | smb, tsrf | clim -> surface -> ice | bilin, con | coupling:306-309 |
| geothermal/sediments -> Yelmo | Q_geo, H_sed | ice (read with remap) | once | kryos.f90:353,356 |
| Yelmo -> coefs | V_ice (Gt) | scalar | | yelmox.f90:185 |
| climate -> coefs | air_anomaly | scalar | | yelmox_climate.f90:337-340 |
| Yelmo -> OBM (bipolar) | fwf from mb_net, smb, bmb, cmb, H_ice, dHidt, f_grnd | scalar | | obm_coupling.f90:246-253 |
| OBM -> climate (bipolar) | tn, ts -> derived coefs -> ocean terms | scalar | | obm_coupling.f90:261-276 |

### 1.7 Coefficient registry (`kryos_forcing`, libs/kryos_forcing.f90)

- One `coefs_class` per run, driver-owned, shared by all domains (header 1-25).
- A coefficient is a named scalar of the forcing time. Two kinds: `tsgen`
  (const, ramps, sin, series from file normalised to [0,1], PI/PID feedback on
  the ice volume) advanced once per step (`coefs_update`, 162-176); `derived`,
  set by a producer with `coefs_set` (178-195): `air_anomaly` (climate),
  `obm_north/south` (OBM).
- Consumers read by name (`coefs_value`, 254-259): the climate terms
  (`climate_terms`), `ice_bnd` (scaled), the OBM forcing.
- It also owns run control that is really forcing-driven: kill switch
  (`coefs_kill`, 283-293), forcing-increment restarts (`restart_every_df`,
  295-326), restarts of tsgen state (`coef_<name>_restart.nc`, 328-388), 1D
  output appended to `yelmo_ts.nc` (390-427; kryos_output.f90:182).
- Group lookup: par file first, then the `[forcing] var_defs` library files
  (`var_defs_find`, libs/yelmox_nml.f90); unset keys from
  `input/yelmox_defaults.nml` templates via `defaults_group`.

On top of it, the climate forcing is "reference + sum of terms", each term =
coefficient x pattern (snapshots difference, uniform, replace) on an explicit
field list (libs/climate_terms.f90 header; docs/_forcing-terms.md), followed by
static transforms (`taper`, `bias`, libs/climate_mods.f90), all contracts
checked at init against `climate_caps_class`.

### 1.8 Regions

`kryos_regions.f90` is misnamed: only `domain_regions_init` (20-46) is about
regions (it sets Yelmo's internal `outfldr` and regional file names). The rest
is experiment physics: `domain_init_marine_ice` (48-59, hard-coded 600/800 m,
-500 m), `negis_update_cb_ref` (61-105, hard-coded -11 ka/-4 ka/0 ka ramp
times), `calc_glacial_smb` (107-137). The regions machinery itself is
fesm-utils `regions` (selection expressions like
`"region:Greenland & ~zone:open_ocean"`).

### 1.9 Startup and restart (`kryos_startup.f90`)

Cold start (`domain_init_state`, 131-195):
1. isostasy reference + state from Yelmo's z_bed_ref/H_ice_ref, check, land on
   Yelmo, `calc_ytopo_diagnostic` (218-250);
2. `step_ice_bnd`; hub refresh; climate (`init=.true.`, rembo equilibrates);
   surface (`init=.true.`, snowpack spin-up); shelf; `couple_to_yelmo`;
3. cb_ref initial guess / NEGIS; `yelmo_init_state(robin-cold)`;
4. `domain_init_ice`: kill shelves, marine ice, `none|equil|recon|recon_ref`,
   thermal equilibration (252-350, Yelmo-specific `yelmo_update_equil`,
   `beta_min` hack at 343-347);
5. hub refresh, surface update again on final geometry.

Restart (`domain_startup`, 57-104): restore bundle; rebuild hub; climate
re-derived with `force=.true.`; surface re-derived only if stateless.

Bundles: folder `restart-<kyr>-kyr/` (389-400), one file per stateful
component (`yelmo_restart.nc`, `isos_restart.nc`, `marine_shelf.nc`,
`chion_restart.nc`/`smbpal_restart.nc`, `rembo_restart.nc`) plus run-level
`bsl_restart.nc` and `coef_<name>_restart.nc` (106-129). Bipolar: run-level
files at top, one subfolder per domain, plus `obm_restart.nc`
(yelmox_bipolar.f90:368-382). Restarts are bit-exact continuations: nothing is
stepped at restart.

### 1.10 Output (`kryos_output.f90`)

- One file per module, on its own grid, grid recorded in the file
  (`yelmo.nc`, `yelmo_sm.nc`, `htopo.nc`, `isos.nc`, `mshlf.nc`,
  `<surface_method>.nc`, `snap.nc|esm.nc|rembo.nc`; 1D `*_ts.nc`).
- Schedules from fesm-utils `timeout` (`[tm_2D]`, `[tm_2Dsm]`, `[tm_1D]`,
  `[tm_rst]`).
- Fixed per-field precision (bit rounding) with a global lossless override.
- The writers for isostasy, shelf and the three surface models live **in the
  framework** (260-361), reading component internals (`isos%out%w`,
  `mshlf%now%tf_shlf`, `sc%ann%PDDs`, ...).
- CMIP output (`cmip_output.f90`) is Yelmo + marine_shelf specific and
  silently requires `grid_shelf == grid_ice` (header 1-5).

---

## 2. Component interfaces today

### 2.1 Summary

| Component | init | update | restart | output | Backend switch |
|---|---|---|---|---|---|
| Ice sheet (Yelmo) | `yelmo_init(yelmo, filename, grid_def, time, domain, grid_name, group, reg, topo_init, topo_pd, cnst)` (kryos.f90:310) + `yelmo_init_state`, `domain_init_ice` | `yelmo_update(yelmo, time)` (coupling:386) | `yelmo_restart_write`; read = `yelmo_restart_read_topo_bnd` + `yelmo_restart_init` (startup:493-497) | `yelmo_write_init/step`, `yelmo_regions_write`, CMIP | none (hard-wired); `with_ice_sheet` only skips the update |
| Isostasy (FastIsostasy) | `isos_init(isos, path_par, group, nx, ny, dx, dy, cnst, domain, grid_name)` (kryos.f90:348) + `isos_init_ref`, `isos_init_state` | `isos_update(isos, H_ice, time, bsl, dwdt_corr)` (coupling:149) | `isos_restart_write`; read via `isos_init_state` with `use_restart` | framework writes `out%z_bed,z_ss,w,we` | none (VILMA2 planned as runtime string) |
| Sea level | `bsl_init(bsl, path_par, time_rel)` | `bsl_update(bsl, time_rel)`; under fastiso written by each `isos_update` | `bsl_restart_write/read` | via `isos_ts.nc` | `[barysealevel] method` string: const, file, fastiso, mixed |
| Climate | `climate_init(cl, method, filename, domain, grid, time, reg, south, sfx, timeline_group, smb_direct, coefs)` (yelmox_climate:108) | `climate_update(cl, out, ts, z_srf, H_ice, z_bed, f_grnd, z_sl, z_srf_ref, dtt, mshlf, coefs, init)` (262) | `climate_restart_write` (rembo only; reads its own restart from its own par file) | `climate_write_init/2D/1D` | runtime `[comps] climate` = snapshots, esm, rembo; rembo also build flag `rembo=1` with stub |
| Surface | chion: `surface_chion_init(sc, filename, lats, group, chion_group, orbit_group, cnst)`; smbpal: `smbpal_init(smb, filename, x, y, lats, group, itm_group, orbit_group)`; smb_simple: `smb_simple_init(smbs, filename, x, y, lat, group, units, init, cnst)` + `smb_simple_set_mask(H_ice_ref)` | chion: `surface_chion_update(sc, t2m, pr, z_srf, z_bed, H_ice, z_sl, time_rel)` (+`_spinup`); smbpal: `smbpal_update_monthly(smb, t2m, pr, z_srf, H_ice, time_bp)` (+`_equil`); smb_simple: `smb_simple_update(smbs, z_srf, t_sl)` | chion, smbpal restart incl. annual outputs; smb_simple stateless | framework writers (kryos_output:296-361) | runtime `[comps] surface_method` = chion, smbpal, smb_simple, climate |
| Shelf | `marshelf_init(mshlf, filename, group, nx, ny, domain, grid_name, reg, xc, yc, cnst)` (marine_shelf:628) | `marshelf_update_shelf(mshlf, H_ice, z_bed, f_grnd, z_sl, dx, depth, to, so, dto)` + `marshelf_update(mshlf, H_ice, z_bed, f_grnd, z_sl, dx, z_srf)` | `marshelf_restart_write/read` | framework writer (279-294) | `[marine_shelf] bmb_method` string: pico, lin, quad, quad-nl, *-slope, anom |
| Geothermal | `geothermal_init(gthrm, filename, grid, domain, group)` | none | none | none (Yelmo writes Q_geo) | `use_obs` bool + `ghf_const` |
| Sediments | `sediments_init(sed, filename, grid, domain, group)` | none | none | none | integer `method` 0/1 |
| ice_bnd | `ice_bnd_init(ibnd, filename, group, domain, grid, coefs)` | `ice_bnd_update(ibnd, H_ice_ref, mask_ice, tau_relax, time, coefs)` | stateless | via `H_ice_ref` in yelmo.nc | `method` = file, scaled |
| Orbit | `orbit_ctl_load(ctl, filename, group)` (per surface model) | `orbit_pars_now(ctl, is, time_rel)` | stateless | none | `method` = forcing_time, fixed_time, fixed_orbit |
| Coefficients | `coefs_init(coefs, path_par, time)` | `coefs_update(coefs, time, var)`, `coefs_set` | `coefs_restart_write/read` | `coefs_write_step` | per coefficient `method` |
| OBM (bipolar) | `obm_init`, `obm_ctl_load`, `obm_coupling_init` | `obm_update(obox, dt, name)`, `obm_exchange` | `write_obm_restart` | `write_obm_update` | `obm_name` (nautilus only) |

### 2.2 Observations per component

**Ice sheet.** No interface layer: the domain holds `type(yelmo_class)` and the
framework reads and writes `yelmo%tpo%now`, `yelmo%bnd`, `yelmo%dyn`,
`yelmo%dta%pd`, `yelmo%reg`, `yelmo%grd`, `yelmo%par` directly. Inputs it
receives are exactly the `bnd%` fields set in `couple_*_to_yelmo` plus
`H_sed`, `Q_geo`, `H_ice_ref`; the full Yelmo boundary type has more
(`fmb_shlf`, `tf_shlf`, `enh_srf`, `calv_mask`, ...; yelmo/src/yelmo_defs.f90
~889-947). `with_ice_sheet = False` still requires Yelmo to be initialized (its
geometry is the static state).

**Isostasy.** Clean contract at the coupling seam: only anomalies cross grids
(`couple_isostasy_to_yelmo`, coupling:186-227, with a long rationale), plus a
consistency check of the reference bed (229-281). But the FastIsostasy types
(`isos_class`, `bsl_class`) appear in the framework and drivers, the output
writer reads `isos%out%w/we`, `dwdt_corr` is a Yelmo-specific correction
field, and under `fastiso` each domain's `isos_update` mutates the shared
`bsl%bsl_now` (FastIsostasy src/fastisostasy.f90:839-849).

**Climate.** The best-designed seam: one `climate_out_class` read by
everyone, capability flags (`has_ocn_shelf`, `has_smb`, `has_Qd`) and init
contracts. Weak points:
- `yelmox_climate_class` embeds all backends' state (snapshots, esm, esm_ctl,
  rembo) in one type (yelmox_climate.f90:78-91); dispatch by `select case` in
  every routine.
- `climate_update` takes a `marshelf_class` (283): esm interpolates the ocean
  to the shelf base using the shelf's parameters (`extrap_shlf`,
  `marshelf_interp_shelf`, 799-847) and `esm_forcing.f90` uses `marine_shelf`
  (line 11). The climate computes a shelf-base field on grid_clim geometry,
  then the shelf remaps it to grid_shelf.
- `climate_init` needs to know the surface choice (`smb_direct`) and the esm
  timeline group (`timeline_group`, read in `esm_init` 757-762).
- REMBO keeps module-global state, runs on its compiled grid (must equal
  grid_clim), reads its own par file and restart, has its own internal
  cadences (climate_rembo.f90:6-9).
- The atmosphere carries `z_srf`, `lapse(12)`, `beta_p` so consumers can
  re-elevate it: good, but it makes the "climate interface" include a
  downscaling model.

**Surface.** Three wrappers with three different signatures; the framework
dispatches in four places (init kryos.f90:384-413, step coupling:477-505,
couple 302-319, output kryos_output:110-118, restart startup:429-432,
516-519). Units are smbpal's (mm w.e./yr), converted to Yelmo m ie/yr in the
coupler using Yelmo's constant `bnd%c%conv_we_ie` (coupling:321). smb_simple
reads Yelmo's `H_ice_ref` for its target mask (kryos.f90:402-404,
coupling:482-483) and needs `tsl_ann` from the climate. Insolation is computed
inside each surface model (each owns an `insol_class` and an `[orbit<sfx>]`).

**Shelf.** Two-call update (`marshelf_update_shelf` to get the shelf-base
ocean, then `marshelf_update`). With `has_ocn_shelf` the framework writes
directly into `mshlf%now%T_shlf/S_shlf/dT_shlf/dS_shlf` (coupling:600-603).
`tf_corr` (an optimization control) lives in `mshlf%now` and is written by
`step_optimize` (104-114). Own `wp = sp` definition (marine_shelf.f90:24).

**Geothermal, sediments.** Data-only; read on the Yelmo grid at init and
assigned to `yelmo%bnd` inside `domain_init` (kryos.f90:352-356), not through
a coupling step. Each has its own `parse_path` for `{domain}/{grid_name}`.

**ice_bnd.** The only existing "prescribed data" backend: a varslice ice
history at `time_rel` (or a coefficient times one slice) into Yelmo's
`H_ice_ref` in fixed/relaxation cells. It uses `yelmo_defs:MASK_ICE_FIXED`
(ice_bnd.f90:29), so it compiles against Yelmo.

**Orbit.** Not a component: a config type loaded separately by each surface
model per domain (orbit_forcing.f90:40-65).

### 2.3 Yelmo-specific leaks (inventory)

| Leak | Where |
|---|---|
| `wp` taken from `yelmo` throughout | kryos.f90:28, kryos_coupling.f90:9, kryos_startup.f90:9, kryos_output.f90:9, kryos_regions.f90:7, obm_coupling.f90:34, drivers |
| `yelmo_class` field in the domain | kryos.f90:164 |
| `yelmo_defs` mask constants | kryos.f90:29, kryos_startup.f90:12, kryos_regions.f90:8, ice_bnd.f90:29 |
| `bnd%` writes | coupling:183, 225-226, 321-336, 355-356, 369-370; kryos.f90:353,356; startup:308-328 |
| internal params mutated | `dyn%par%till_method = -1` (kryos.f90:462, 539), `tpo%par%topo_rel*` (relax), `dyn%par%beta_min` (startup:343-347) |
| Yelmo routines called by the framework | `calc_cb_ref`, `calc_ytopo_diagnostic`, `yelmo_update_equil`, `yelmo_init_topo`, `yelmo_restart_read_topo_bnd`, `yelmo_regions_update`, `optimize_cb_ref`, `relax_update` (startup:9-16, coupling:17-18) |
| Driver reads Yelmo internals | `dom%yelmo%reg%V_ice * dom%yelmo%bnd%c%rho_ice` (yelmox.f90:185; bipolar:209-210), `yelmo%grd%G%nx` in prints |
| Output | `YELMO_VARS_2D` names (kryos_output:31-40), `yelmo_regions_write`, CMIP |
| Grid spacing units | all `dx_*` in "Yelmo dx units" |
| OBM freshwater flux | reads `tpo%now%mb_net, smb, bmb, cmb, dHidt` and `bnd%c` (obm_coupling:246-253) |
| Region output | `domain_regions_init` sets `yelmo%outfldr`, `yelmo%reg%fnm` |
| Physical constants file | `input/yelmo_phys_const.nml` (kryos.f90:53); `phys_const` set in both `[domain]` and `[yelmo]` |
| Build | every kryos object compiled with `$(INC_YELMO) $(INC_ISOSTASY) $(INC_CHION) $(INC_INSOL)` (config/Makefile_yelmox.mk:63-90) |

---

## 3. What is good and worth keeping

| Pattern | Where | Why it works |
|---|---|---|
| **Hub as geometry source of truth, updated with anomalies** | htopo, coupling:541-570 | Components never see a coarse round trip of the hi-res bed/ice; the finest information is kept and only the model's *change* is refined. Decouples component grids from the ice-model grid. Makes the ice model one producer of geometry among others, which is the kryos intent. |
| **Only anomalies cross the isostasy boundary** | coupling:186-227 | Displacements are long-wavelength and grid-insensitive; absolute fields would smooth away troughs and sills. Paired with an explicit reference check (229-281) that catches restart mismatches. Generalizes to any slow, smooth field (VILMA2 plan keeps it). |
| **Grids as strings + lazy, cached, method-typed maps** | fesm-utils coupler, kryos.f90:714-750 | Moving a component to another grid is a par-file change. Identity is free. Disk cache amortizes expensive hi-res maps. Grid table on disk (`maps/grid_<name>.txt`) is the single definition. |
| **Grid default chain** | kryos.f90:574-578 | A minimal setup names one grid; refinements are opt-in. |
| **Coarsen conservatively, refine bilinearly** | `remap_method_smooth` 687-703 | Mass-like fields keep their integral; smooth fields don't become blocky. Should be a coupler policy, keyed on the field's nature rather than Yelmo dx. |
| **Re-elevation of the atmosphere across grids** | `remap_atm_to_surface` 508-539, `climate_out` lapse/beta_p | Lets climate run coarse while the surface keeps elevation detail. |
| **Neutral output type + capability flags + init-time contracts** | `climate_out_class`, `climate_caps_class`, `terms_check`, `climate_mods_init` | Downstream code is backend-agnostic; invalid combinations fail at init with a named reason instead of silently. This is the template for every kryos component interface. |
| **Run-level vs domain-level ownership** | bsl, coefs, OBM in drivers | The first explicit separation of global and regional state; a second domain shares them without code changes. |
| **Coefficient registry (named scalars, shared, restartable, derived)** | kryos_forcing.f90 | A generic scalar bus: forcing indices, feedback controllers, and model-produced scalars (OBM, air anomaly) through one mechanism with restart, output and kill-switch. Natural place for global scalars (GMT, CO2, bsl, orbit). |
| **Forcing as an equation (reference + coefficient x pattern)** | climate_terms.f90, docs/_forcing-terms.md | The par file transcribes physics; one vocabulary replaced six ad-hoc mechanisms (docs/_forcing-design.md section 2). |
| **One forcing time** | `ts%time_rel` for climate, orbit, bsl, ice_bnd | Eliminated the LIS insolation bug class (two clocks). Exceptions are named (fixed orbit, esm calendar). |
| **Restart bundle: folder per time, file per stateful component, run-level files on top, domain subfolders** | kryos_startup.f90:106-129, 402-530; bipolar:368-382 | Exact continuation, chunked runs (docs/run-chunks.md), forcing-increment branching. Derived state (hub, stateless climate) is rebuilt, not stored. |
| **Restart = restore, never step** | startup:57-104 | Bit-exact restarts; spin-up logic cannot accidentally rerun. |
| **Per-module output on native grid, grid in the file, fixed precision per field** | kryos_output.f90 | No interpolation artifacts in output; compressible; each component's file stands alone. |
| **Step primitives that are no-ops when inactive; sequence written out in the driver** | kryos_coupling.f90 | The coupling order is readable in one place; switching components off needs no branching in the driver. |
| **Optional dependency behind a same-interface stub** | climate_rembo.f90 / climate_rembo_stub.f90, config/common.mk:63-73 | Default build has no rembo dependency; selecting it at runtime without the build fails at init with an actionable message. VILMA2 and CISM plans reuse it. |
| **Regions on every component grid; components make their own masks from expressions** | kryos.f90:201-213, fesm-utils regions | No precomputed per-grid mask files; one regions dataset serves all grids. |
| **Physical constants loaded once, handed to all components** | kryos.f90:278 | Consistency of rho_ice, rho_sw, g across components. Should be run-level. |
| **Defaults files + `nml_validate` + library lookup** | yelmox_nml.f90, input/yelmox_defaults.nml, input/forcing_db.nml | Par files list only differences; typos fail; shared groups live once. runme `par_defaults`/`par_library` let `-p` add groups. |
| **Domain suffix on groups** | `domain_init(..., group_suffix)` | Several domains in one file, `-p` can target one. (Mechanism good; see 4 for limits.) |
| **Per-component wall-time accounting** | fesm-utils timer, yelmox.f90:142-144 | Cheap, made the `dt_isos` cadence decision obvious (CHANGELOG v2.4). |
| **Prescribed-data component prototype** | ice_bnd.f90 | varslice + forcing time + optional coefficient scaling is exactly the "data backend" pattern kryos wants for every component. |

---

## 4. Legacy to shed

### 4.1 Structural

1. **Ice model at the centre.** `kryos_domain%yelmo`, `couple_*_to_yelmo`,
   `remap` signature on `kryos_domain`, Yelmo dx units, `wp` from Yelmo. Every
   component's output is landed on the ice grid by framework code that knows
   Yelmo's `bnd`. A neutral ice interface (inputs/outputs as in the CISM plan,
   docs/_cism-backend-plan.md "Plan" step 1) is the precondition for
   everything else.
2. **All backends as fields of one type.** `smb`, `smbs`, `schn` in the domain;
   snapshots/esm/rembo in `yelmox_climate_class`. Dispatch by string in 4-6
   places per component. Replace with one polymorphic slot per component (or
   one dispatch point) so adding a backend touches one file.
3. **Cross-component reach-ins.** climate takes `marshelf_class`; step_shelf
   writes `mshlf%now%T_shlf`; step_optimize writes `mshlf%now%tf_corr`;
   smb_simple reads Yelmo's `H_ice_ref`; glacial smb scaling reads climate
   internals and Yelmo latitude; isostasy reads Yelmo directly instead of the
   hub. Each is a hidden edge in the dependency graph.
4. **Component output written by the framework.** `isos_write_step`,
   `mshlf_write_step`, three surface writers in kryos_output. Output knowledge
   belongs to the component (or to a declared field list), otherwise every new
   backend needs framework edits.
5. **Two drivers duplicating the loop** (yelmox.f90 vs yelmox_bipolar.f90).
   Bipolar differs in: no `run_step`/`calendar` (bipolar:99), no CMIP output,
   explicit `dom_north`/`dom_south` instead of an array, interleaved order (all
   isostasy before OBM), different restart layout, `restore_bsl` flag to paper
   over the bsl restore asymmetry (startup:57-104), different timer table. The
   asymmetry argument (north <-> fn, south <-> fs) is an OBM detail that a
   generic "global component with per-domain ports" would absorb.
6. **Experiment logic in the framework.** NEGIS (kryos.f90:135-145, 436-485;
   kryos_regions.f90:61-105, hard-coded ka), glacial smb scaling, `lim_pd_ice`
   (hard-coded -4 m/yr, coupling:334-337), `init_marine_H` (600/800/-500 m),
   recon cold start (`beta_min = 100`, 200 yr, 50 m thresholds), relax/opt
   tied to Yelmo internals and the shelf. These belong to the ice backend
   (spin-up tooling) or to an experiment plugin layer, not to the domain type
   and its ctl.
7. **`domain_ctl` as a grab bag** (component selection + grids + cold start +
   experiment switches + output). Kryos should split component config (per
   component), run control (driver), and experiment hooks.
8. **Cadences inside components.** REMBO's own `dtime_emb/dtime_smb`,
   FastIsostasy's internal `dt_*`, chion's annual cycle per call regardless of
   `dtt`. The framework can't reason about them.

### 4.2 Conventions

- **Required-key namelist groups without defaults** for framework groups
  (`[domain]`, `[comps]`, `[sim]`, `[output]`: "every key is required",
  kryos.f90:56-59) next to defaulted forcing groups: two regimes.
- **Integer method codes** (`isos.method = 3`, `sed.method = 0`), mixed
  `.true.`/`True`, duplicated `phys_const` in `[domain]` and `[yelmo]`.
- **Suffix-based multi-domain configs**: one 1403-line par file for bipolar;
  shared vs suffixed groups decided by convention (Yelmo physics shared,
  `[barysealevel]` unsuffixed); per-domain physics needs pointer-field
  re-wiring (bipolar nml header). A per-domain file (or table) would scale
  better.
- **Hard-coded paths/names**: `maps/` (kryos.f90:50),
  `input/yelmo_phys_const.nml` (53), `input/yelmox_defaults.nml`
  (yelmox_nml.f90:17), output names, restart folder pattern, `./` output dir,
  hydro masks under `ice_data/` (bipolar nml) while everything else is
  `fesmdata/`, run-dir links set by runme (`input`, `ice_data`, `fesmdata`,
  `isostasy_data`, `maps`).
- **Each component re-implements path keyword substitution** (`parse_path` in
  geothermal, sediments, marine_shelf).
- **Per-library precision kinds** (`smbpal_precision`, local `wp = sp` in
  marine_shelf, geothermal, sediments, esm_forcing; `dp` in rembo/VILMA2) with
  ad-hoc `real(..., wp)` conversions.
- **Compatibility shims**: `series_<name>_restart.nc` fallback
  (kryos_forcing.f90:361-370), `[coupling]` -> `[comps]` renames, docs that
  still mention `tsf` and an old `kryos_domain` (docs/multigrid.md).
- **Ad-hoc error handling**: `stop 1`, `error stop 1`, bare `stop` mixed.

### 4.3 Known bugs from the v2.4 sweep that are structural

From docs/_v24-release-sweep.md (some since fixed; the classes remain):
effective output switches resolved in several places (item 10), derived
coefficients with no single producer when two domains share a registry
("bipolar both domains set dTo", item 48), run-level objects restored
asymmetrically between drivers (item 51), `counter_restart` not restored,
`fwf_definition` ignored in `ice2ocean` (item 7), coordinate units passed
wrongly to smb_simple (item 5, km labelled m). All are symptoms of implicit
contracts between framework and component.

---

## 5. Time stepping and the coupling lag

### 5.1 Who owns time

- The driver owns `tstep_class` (fesm-utils timestepping) and `dtt`
  (yelmox.f90:86). `ts` carries several clocks: `time` (model), `time_rel`
  (forcing time; `const` holds it fixed), `time_elapsed` (since run start),
  `time_cal`.
- Clock use is component-specific:

| Consumer | Clock |
|---|---|
| Yelmo, isostasy (`isos_update`) | `ts%time` (coupling:149, 386) |
| bsl, coefs, snapshots, rembo insolation, orbit, ice_bnd, NEGIS, surface (orbit) | `ts%time_rel` |
| esm | `ts%time` (calendar), periods read from the timeline group |
| `dt_clim`, `dt_isos`, `dt_cmip` cadences; relax/opt windows | `ts%time_elapsed` |

  Isostasy integrates on model time while bsl follows the forcing time: in a
  `const` run these diverge by design, but the mismatch is implicit.

### 5.2 Cadences

- `cadence_due(time, dt)` (kryos.f90:215-222): `mod(nint(time*100),
  nint(dt*100)) == 0` on `time_elapsed`. Silent if `dt` is not a multiple of
  `dtt` (it just fires at the common multiples); 0.01-yr resolution; `dt <= 0`
  disables.
- Climate is held between updates; isostasy holds the bed and integrates the
  whole interval in one call; surface and shelf run every step; chion runs one
  annual cycle per call regardless of `dtt`.
- Cadences counted from run start: chunk boundaries must be multiples of
  `dt_clim`, `dt_isos` and output intervals, and relax/opt windows restart per
  chunk (docs/run-chunks.md "What changes after the first chunk").
- On restart, the climate is re-derived (`force=.true.`) and equals the
  original only if the restart time falls on the `dt_clim` cadence
  (startup:76-81). Derived coefficients (`air_anomaly`) only refresh on climate
  updates.

### 5.3 The coupling lag

Within step n (t_{n-1} -> t_n):

| Input to Yelmo | Computed from | At |
|---|---|---|
| z_bed, z_sl | isostasy integrated to t_n with load H(t_{n-1}) | t_n |
| smb, T_srf | climate/surface at end of step n-1, geometry H(t_{n-1}), forcing time_rel(t_{n-1}) | t_{n-1} |
| bmb_shlf, T_shlf | shelf at end of step n-1 | t_{n-1} |
| Qd | climate at end of step n-1 | t_{n-1} |
| H_ice_ref (ice_bnd) | time_rel(t_n) | t_n |

So the bed is "end-of-interval" while the surface forcing is
"start-of-interval": a mixed explicit scheme. In bipolar the OBM steps with the
previous step's freshwater/atmosphere, then `obm_exchange` updates the
derived coefficients and re-applies the ocean terms before the shelf step
(yelmox_bipolar.f90:233, 267-273), so OBM->shelf has no lag but
ice->OBM->ice does.

### 5.4 Fragile points

- Lag is implicit in call order, not declared; the comment on
  `couple_to_yelmo` (coupling:160-170) is the only documentation.
- Mixed clocks per component (5.1) with no declared clock in the interface.
- Cadence arithmetic on `time_elapsed` breaks with non-multiples and across
  chunks.
- Restart correctness depends on cadence alignment (climate) and on the order
  of run-level restores (`restore_bsl`).
- `fastiso` bsl: each domain's `isos_update` adds its contribution in sequence,
  so within one step the second domain sees the first domain's update
  (order-dependent).
- Components with internal sub-stepping (REMBO, FastIsostasy adaptive,
  chion daily) are opaque to the driver.

A from-scratch design should own one clock, have each component declare its
clock (model or forcing), cadence and which inputs it accepts lagged, and let
the scheduler derive the order and the restart-safe points.

---

## 6. Regional vs global

### 6.1 What is implicitly global today

| Object | Scope today | Notes |
|---|---|---|
| `bsl_class` | run (driver) | FastIsostasy type; `[barysealevel]` unsuffixed; under `fastiso` integrated from every domain's load in sequence |
| `coefs_class` | run (driver) | shared registry; derived coefficients can be set by any domain (last writer wins) |
| OBM | run (bipolar driver) | global box model with per-hemisphere ports (fn/fs, tn/ts) |
| Timeline, output schedules | run | shared `[ctrl]`, `[tm_*]` |
| Orbit / insolation | **per domain, per surface model** | physically global; `[orbit<sfx>]` duplicated; each surface model builds its own insol table |
| Physical constants | per domain | physically per planet |
| ESM / paleo global products | per domain | each domain reads global files (`{grid_global}`) and remaps online |
| Yelmo physics groups | shared by convention in bipolar | not a model property, a config shortcut |

Regional: Yelmo, FastIsostasy (regional by construction: padded FFT domain),
marine shelf, surface, climate backends (snapshots/esm/rembo all produce
regional fields), hub, coupler, regions.

### 6.2 How bipolar couples two domains

- Two full `kryos_domain`s, each with its own coupler, hub and grids
  (yelmox_bipolar.f90:71, 309-344); no grid shared between them.
- Shared: `bsl` (each domain's `isos_update` adds its load change), `coefs`
  (e.g. `gmt`), OBM.
- Exchanges through scalars only (obm_coupling.f90): atm->OBM from a
  coefficient (`gain * coef`), ice->OBM freshwater flux computed from Yelmo
  mass-balance fields on the Yelmo grid within a hydrographic mask, OBM->ice
  via derived coefficients `obm_north/obm_south` consumed by each domain's
  ocean terms (`climate_apply_ocean` re-applies them after the climate step).
- Ordering is hand-interleaved in the driver (isostasy both -> OBM -> ice,
  climate, surface both -> exchange -> shelf both).

### 6.3 What a truly global isostasy/sea level needs (from docs/_vilma2-backend-plan.md)

- **Two-level component**: `isos_run_class` (driver-owned, one per run;
  absorbs `bsl_class`) and `isos_domain_class` (per domain: only the
  `dz_bed`, `dz_ss` on `grid_isos` that the ice model needs).
- **Three-phase step** so N domains feed one solve:
  `isostasy_load` (per domain: remap H_ice to the global grid, splice into the
  load) -> `isostasy_solve` (once) -> `isostasy_land` (per domain: global
  anomalies back to `grid_isos`). FastIsostasy maps to load=update, solve=no-op.
- **A global grid** (Gauss `GLOBAL-T64`) known to the coupler, maps from each
  domain grid to it (per-domain couplers would each hold a domain->global map;
  a run-level coupler is cleaner). Lon-lat global products remapped by the
  coupler. `remap_method_smooth` in Yelmo dx units is meaningless here: load
  `con`, fields back `bilin`.
- **Global load splice**: `H_global = f*H_domain + (1-f)*H_ext`, with `f` the
  domain footprint fraction; ice outside the domains from a fixed reference or a
  global ice history (ICE-6G_C), read on its own cadence.
- **Global reference topography** with domain references spliced in;
  reference check becomes footprint-masked on the global grid.
- **One cadence per run** (`[isostasy] dt`) replacing per-domain `dt_isos`.
- **Run-level restart** (`vilma_restart.nc`), run-level output (`vilma.nc`),
  per-domain neutral output (`z_bed, z_ss, dz_bed, dz_ss, bsl`).
- **bsl definition from the global model** (barystatic rise as ocean mean of
  C*rsl), `z_sl` with geoid and rotation: gravitationally consistent.
- Backend-side changes (export `dz_bed`/`dz_ss` separately, `libvilma.a`,
  FesmData viscosity format) and a stub (`vilma=1`).

Generalization for kryos: this load/solve/land split is the generic pattern
for any global component (climate model, ocean box model, sea level) serving
regional ones. The OBM in bipolar is a hand-coded instance of it with scalar
"grids" (boxes). The coefficient registry is the degenerate 0-D case.

---

## 7. Build system and dependencies

| Element | Where | Assessment |
|---|---|---|
| configme manifest with pinned deps | `.configme/manifest.toml` (`fesm-utils`, `yelmo`, `FastIsostasy`, `rembo1`, `chion`, `insol` at `stable`/tags) | Keep. Reproducible installs; `.install.sh` is a generated replay script (clones, symlinks `fesm-utils` into every component, per-package configure and builds). Superseded the submodule plan (docs/_submodules-plan.md). |
| Generated compiler/machine header in `Makefile` | Makefile top block ("Generated by configme") | Keep: one place for FC/flags/netCDF per machine. |
| `config/common.mk` | dependency roots, serial/omp variants, `LFLAGS` | Keep the idea; but all libraries are linked always (common.mk:97): Yelmo, FastIsostasy, chion, insol are mandatory. |
| Stub for optional deps | `rembo ?= 0` (common.mk:63-73), `climate_rembo_stub.f90` | Keep and generalize: every non-core backend (Yelmo, PISM, CISM, VILMA2, chion, REMBO...) should be optional with a stub, so kryos can be built with any subset. |
| Archive dependency rule | end of `config/Makefile_yelmox.mk`: objects depend on `libyelmo.a`, `libisostasy.a`, `libfesmutils.a` | Keep: prevents stale-type `.mod` layouts (fesmc/FastHydrology#10). |
| Sanity check for heavy fesm-utils deps | Makefile:110-130 (`fesmutils-deps-check`, LIS/FFTW) | Keep: fails early with instructions. |
| Recursive builds of each dep | `yelmo-static`, `isostasy-static`, ... | OK, but transitive deps (FastHydrology, elsa, tracer) must be named on the link line because static archives don't carry them (common.mk:33-45). |
| Hand-written per-object dependency rules | config/Makefile_yelmox.mk | Shed: easy to break, and the kryos objects list every backend as a dependency (63-90), so the core can't be built without them. Use generated dependencies (makedepf90/fortdepend) or CMake/fpm. |
| Single include dir for all `.mod`s | `libyelmox/include` | Risk of name collisions (`htopo.mod`, `insolation.mod` stale, CHANGELOG v2.4 "Changed"). |
| runme `info.json` | `.runme/info.json` | Keep: exe aliases, `par_defaults`, `par_library`, run-dir links. |
| `obm_coupling.f90` in the driver dir | yelmox_bipolar/ | A bridge compiled only for one program; generalizes to "coupling modules for global components". |

---

## 8. Other things a from-scratch design should know

1. **The CISM plan already scoped the neutral ice interface** (inputs: z_bed,
   z_sl, smb, T_srf, bmb_shlf, T_shlf/t_ocn/s_ocn, Q_geo, H_sed, Qd, mask_ice,
   regions/basins; outputs: H_ice, z_srf, z_base, f_grnd, masks, dHidt,
   smb/bmb/cmb/mb_net, calving, ux_s/uy_s, z_bed_ref/H_ice_ref, dzbdt_corr) and
   the open decision "abstract type vs select case + stubs"
   (docs/_cism-backend-plan.md, "Open decisions" 2). It recommends an abstract
   type so host programs (CLIMBER-X with SICOPOLIS) can add backends.
   CLIMBER-X `src/ice/ice_def.f90` is a second precedent.
2. **CISM practicalities**: INI config files, its own netCDF I/O and restarts,
   `stop`/`mpi_abort` on errors, module-global registry (max 8 instances),
   topg/eus halo + marine mask recompute needed after external changes, SMB in
   m ice/s, temperatures in degC, projection-area integrals. A kryos ice
   interface must allow backends that own their I/O and restarts.
3. **PISM** (not studied in yelmox) will have similar traits (PETSc/MPI, own
   config/IO). Plan for components that are MPI-parallel while kryos is
   OpenMP/serial.
4. **Unit conventions are converted in the coupler** (smb mm w.e./yr ->
   m ie/yr with Yelmo's constant). A kryos interface should fix SI (or
   declared) units per field.
5. **Precision**: Yelmo sp, VILMA2/REMBO dp, others local kinds. Pick one
   interface precision and convert at backend boundaries only.
6. **"Data" backends exist piecemeal**: `ice_bnd` (prescribed ice in strips),
   `with_ice_sheet = False` (static geometry), `bsl method = file|const`,
   geothermal/sediments (pure data), snapshots and esm (data + anomaly logic),
   smb_simple (parametric). varslice (time-sliced, keyword paths, online remap)
   is the reader to build a uniform data backend on. A generic
   "prescribed component: read field X at forcing time from group G" would
   cover ice (H_ice, z_srf), isostasy (z_bed anomaly), smb, bmb, climate.
7. **Forcing builder and check scripts** (`scripts/forcing-builder/`,
   `scripts/check_par_nml.py`) consume the par-file schema; a kryos schema
   should be machine-readable from the start (defaults files are a good
   basis).
8. **Climate grid limitation**: one `grid_clim` for both the high-res
   reference and the coarse transient forcing (kryos.f90:107-111,
   docs/multigrid.md open issue 3). A kryos climate interface should allow
   reference and anomaly on different grids.
9. **Hemisphere is derived from the climate grid latitude** (kryos.f90:365)
   and used for seasons and lapse rates; a global/climate component would need
   it per cell, not per domain.
10. **The OBM shows a global component with heterogeneous ports**: scalar
    inputs per hemisphere, scalar outputs consumed via coefficients. A generic
    port/registry concept (field or scalar, with a grid or a box id) would
    remove the bespoke bridge.
11. **Testing**: `tests/test_coefs.f90`, `tests/test_timestepping.f90`,
    `tests/diff_nc.py`, bipolar `run_tests.sh`, and Levante bit-identity checks
    for refactors (docs/_par-files-plan.md "Verification"). Bit-identical
    refactor phases have been the working method; kryos should keep a
    regression harness against yelmox outputs.
12. **Docs drift**: docs/multigrid.md lists an outdated `kryos_domain` and
    "transient forcing (`tsf`)"; CISM plan line numbers are pre-split. Plans in
    `docs/_*.md` are the real design record; kryos should keep ADR-style
    decision docs close to the code.
