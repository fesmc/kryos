# kryos phase 4 — plan

Status: revised 2026-10-29 with Alex's decisions (§17, Q1-Q14 all decided); no code. Phase 4
(`design.md` §14) is `geo = vilma2`: the VILMA2 global solid Earth and sea level
coupled to Yelmo, including its restart. Target (the science case): a running
VILMA2-Yelmo setup with two active Yelmo domains, NH (`yelmox_North`, NH-32KM) and
ANT (32 km ice grid, hub ANT-16KM as in the shipped base); climate and surface from
snapshots + tsgen + chion, or from `smb_simple`. No `obm`, no bipolar setup (much
later). The machinery stays general (N domains, any backends). Phase 4 starts once
phase 3's two-domain work is in (`K3-sci` registered).

Refs: `design`, `core`, `contracts`, `adapters`, `regression`, `phase2` =
`docs/design.md`, `docs/dev/*-design.md`, `docs/dev/phase2-plan.md`; `yx-plan` =
`~/models/yelmox/docs/_vilma2-backend-plan.md` (superseded as a yelmox plan, its
physics carries over); `ht` = the approved VILMA2 change 4 design
(`vilma2-change4-design-approved.md`, scratchpad of the planning session);
`phase3` = `docs/dev/phase3-plan.md` (§10 here and §3.5 there agree). Code references are to
kryos `main` f57e826 (code of 93aab84; line numbers updated from cbabc7b: phase 2 WP1 moved the
run files) unless noted; `vilma2 main` = d8041e8, which
has changes 1-4 (host terms merged and pushed). The VILMA2 host API is the README
section "Host-resolved terms" and `doc/design.md` §15.

## 0. Summary

1. Two instances on two levels (design D5): a **global** `geo_glob` (`vilma_earth_t`,
   owns the `solid_earth`, lives on VILMA2's Gauss grid) and one **regional port**
   per domain (`geo_nh`, `geo_ant`, `vilma_geo_t`, stateless; it holds the scattered
   anomalies between global steps and is what the unchanged `geo → ice` rows read).
2. The loads reach the Earth by a **gather** run by the core, not by wiring rows:
   per domain, on the ice grid, the host-resolved terms `I_g`, `C`, `F` (Yelmo's
   flotation, external ice outside the footprint), conservatively merged to the
   Gauss grid with the coverage `w` (`coupler_merge_*`, which exists), summed over
   NH + ANT; the external ice for the `1 − w` part from a history.
3. The **reference** is built at run time from the same code: background X
   (`IceHistory-PaleoMIST` at 0 ka, `z_bed` and `H_ice` on `GLOBAL-0.5DEG`, a datasrc
   of the globe) plus each domain's present day (`dta%pd`) gathered as in a step. No
   FesmData product is needed; the only FesmData gap is a grid (below).
4. The **scatter** splits `dz_bed`, `dz_ss` (anomalies) and `dphi` to each domain's
   geo grid by the wiring rows `geo_g → geo` (bilinear, in the domain's coupler); `z_sl = dphi + dz_ss` at the ice. This
   needs one small contract change: a scalar `dphi` beside `bsl` (Q1).
5. **Time**: VILMA2 steps on the geo cadence (10 yr, revisited after G7), loads at the
   previous step (lag 1, as FastIsostasy), the ramp is VILMA2's. Cold start: VILMA2 is
   initialised at the reference before the domains, the domains start on the
   reference bed, then the t0 load is gathered, VILMA2 is spun up (its own spin-up,
   holding the t0 load) and seeded (`dt = 0`) between blocks c) and d), the anomalies
   are landed and the domains reconcile (§8.1).
6. **Restart**: VILMA2's `vilma_restart_write/read` file is the backend file of
   `geo_glob` in the bundle; the reference is rebuilt before the read from static
   ice exports that `PASS_INIT` already publishes, so it is identical by
   construction and VILMA2 checks it (`check_reference`). Gate: chunk == cont.
7. **Phase 3 and phase 4**: phase 3 delivers two domains in one run, the `bsl` reduction for N domains
   and the kryos self-baselines `K3-sci` (§10). The
   whole global side is phase 4: grids and a coupler in the globe, rank-2 fields on
   global instances, rank-2 globe ↔ domain exchange, footprints, the global geometry hub
   (assembly, checks, output), a generic globe step, global grid names in the manifest.
8. **New findings** that change the inputs (§1.4): FesmData's `GLOBAL-T170` is 512 × 256,
   VILMA2's `lmax = 170` grid is 680 × 342; Yelmo's `bnd%H_ice_ref` is not a fixed
   reference (the ice adapter must export `dta%pd`); the contract's scalar `bsl`
   cannot carry VILMA2's `Δφ`.
9. **Output precision**: every phase 4 output (global sea-level file, VILMA2
   diagnostics on the Gauss grid, per-domain geo port output) uses ncio rounding
   (`prec_abs`/`prec_rel` per variable); restart files stay lossless (§9.1).
10. Work: 9 packages (§12), 2 L, 5 M, 2 S. Gates: no yelmox baseline, so: standalone vs
    kryos with prescribed ice, mass and `bsl` budgets, chunk == cont, an ANT
    deglaciation, the cost of T170 + 3D viscosity (§13, §14).
11. Q1-Q14 are decided (§17), none open. FastIsostasy + VILMA2
    is a later package (§18).

## 1. Scope, inputs, findings

### 1.1 In and out

In: `geo = vilma2` on the globe and as the domains' port; the global grids, coupler,
footprints and geometry hub; the gather, reference, scatter, time stepping, spin-up,
restart, output (rounded, §9.1); the packages, tests and gates. The science case has two
active Yelmo domains (NH-32KM, ANT-32KM); the code is written for N domains and any
backends.
Out: `obm` and the bipolar setup (much later); FastIsostasy together with VILMA2
(VILMA2 change 5, §18); `bsl` from a file (`geo = data`, phase 2 WP6 deferred); other
ice backends (phase 5).

### 1.2 Settled (not reopened)

Global GIA is built in kryos only. The yelmox plan's physics and decisions carry
over: global ice field, footprints (`mask_ice /= MASK_ICE_NONE`, fixed), host-resolved
terms on the host grid with the host's flotation, reference from `dta%pd`, `ice_ext`
(`ref | file`, anomaly against the history's own present day), densities from the host,
init and restart sequence, VILMA2 as in CLIMBER-X (`lmax = 170`, PREM, Bagge structure
6, 3D viscosity). The global reference topography is assembled at run time
(`design.md` §6 "Mechanism", §14). VILMA2 main d8041e8 has changes 1-4: `dz_bed`, `dz_ss`,
`dphi`, `bsl` (ocean mean of `C·rsl`) and `bsl_ive`, densities as parameters, the Bagge
reader, `make vilma-static`, and the host terms (`ht`; README "Host-resolved terms",
`doc/design.md` §15): Gauss passthrough of the weighted `w`, `I_g`, `C`, `F`, one global
`w` summed by the host, `h_ice`/`h_ice_eq` = the `1 − w` part, the restart carries the
reference and the last host terms.

Translation of the yelmox plan into kryos terms:

| yx-plan | kryos |
|---|---|
| `kryos_sealevel` (run), `[sealevel] method = "vilma2"` | `geo_glob` (global geo instance, `[globe] geo`), `kryos_globe_ice` (the global ice field), `gia_*` run procedures |
| `kryos_isostasy` (domain), `[isostasy] method = "global"` | the domain's geo instance with backend `vilma2` (regional port) |
| `dta%pd` | ice exports `z_bed_ref`, `H_ice_pd` (static; §4.4) |
| `ice_ext*`, `z_bed_ref_file` | `[globe_ice]` group (§3.4) |
| `vilma_restart.nc` | `geo_glob_restart.nc` in the bundle (§8.3) |
| `sealevel.nc`, `sealevel_ts.nc` | `geo_glob.nc`, `geo_glob_ts.nc` (§9) |
| `barystatic`, `const`, `file` methods | not needed in phase 4 (fastiso `bsl` stays; `geo = data`, phase 2 WP6) |
| yelmox_bipolar, `obm` | out (much later) |
| yelmox phase 4 (FastIsostasy + VILMA2, channels) | later package (§18) |

### 1.3 States checked

kryos `main` f57e826 (rebased onto it); vilma2 main d8041e8 (host
terms merged); fesm-utils `phase0-kryos` 5662899 (the checkout kryos links; `grid_init_gaussian`
defaults to "Spherical Earth"); climber-x
478b00a; FesmData local registry (IceHistory, SolidEarth, Topo, Grids); yelmo (local);
`phase3-plan.md` as revised 2026-10-29 (two domains, no obm, §10).

### 1.4 Findings of this review

| # | Finding | Evidence | Consequence |
|---|---|---|---|
| F1 | FesmData `GLOBAL-T170` is 512 × 256 ("as used by VILMA"); VILMA2's default de-aliased grid for `lmax = 170` is `nlat = 2 lmax + 2 = 342`, `nphi = 4 lmax = 680` | `FesmData/docs/domains.qmd:87`; `vilma2 src/vilma_coupling.f90:690-708`, `grid_is_gauss` `:710-730` | kryos names and writes its own Gauss grid `GAUSS-L<lmax>` (Q9); never `GLOBAL-T170` |
| F2 | Global products that exist: `IceHistory-PaleoMIST` (`z_topo`, `z_bed`, `H_ice`, `dz_sl`, `f_ice`, 80-0 ka every 2.5 kyr) on `GLOBAL-0.25/0.5/1/2/5DEG`; also ICE5G, ICE6G_C, GLAC1D (+nn9927, nn9894), PRISM4, PaleoMIST_max. `Topo` has **no** Global product (domains only). `SolidEarth-Bagge2021` exists on `GLOBAL-T170` (`log10_visc` on `depth` × `structure`) | `FesmData/registry/Global/IceHistory.toml`, `SolidEarth.toml`; `Topo/README.md` | X default = PaleoMIST at 0 ka (§5.2); the viscosity file exists |
| F3 | Yelmo's `bnd%H_ice_ref`/`bnd%z_bed_ref` equal `dta%pd` only at init (`yelmo_ice.f90:1267-1268`); `ice_bnd` and recon overwrite `H_ice_ref` later (phase 2 WP4, WP7: "re-published every step") | yelmo `yelmo_ice.f90`; `phase2-plan.md` WP4, WP7; kryos `kryos_ice_yelmo.f90:411-412` exports `bnd%` values | the ice contract needs a static `H_ice_pd` (phase 2 WP3 already adds it, copy of `dta%pd%H_ice`) and a static `footprint` |
| F4 | fesm-utils already has `coupler_merge_begin/add/end`, `coupler_split`, `coupler_integral`, `coupler_set_area`, `coupler_owner_reduce`, and `grid_init_gaussian` (Gauss grid from `nlat`, `nlon`; planet default "Spherical Earth" since fesm-utils 5662899, so the cell areas are the quadrature areas `R²·dλ·w_j`) | `coupler.f90:100-104, 229, 295-440, 466`; `coords/coordinates.f90:271` | no fesm-utils change is required; design §12 row 3 can be marked done |
| F5 | `geo_declare` declares only 0-D fields for a global geo; `globe_ctx` sets no grid; `exchange_globe` stops on a non-0-D destination; `GLOBE_ORDER` is `forcing, geo`; the geo regional contract always declares `H_ice` | `kryos_geo.f90:203-207`; `kryos_globe.f90:135, 523-531`; `kryos_restart.f90:78` | the rank-2 globe is phase 4 (WP1, WP4; §10) |
| F6 | The contract's `z_sl = bsl + remap(dz_ss)` (`ADD_SCALAR`, aux `bsl`) needs a scalar that, added to the anomaly, gives the sea-surface change. VILMA2's `dz_ss = N + Δφ`; its `bsl` is the ocean mean of `C·S`, which differs from `Δφ` by the ocean mean of `C(N − u)` (cm to m) | `kryos_wiring.f90:116`; `vilma_coupling.f90:609-639`; yx-plan "One budget" | Q1: add `dphi` |
| F7 | VILMA2's `dt = 0` branch solves with `ice0`, not the new ice (documented known issue, not fixed by change 4, `ht` §0) | `ht` §0, §8 Q4 | the seed must follow a spin-up under the same load (§7.4) |
| F8 | VILMA2's restart time label is `se%time` (VILMA2's own clock from init, advanced by spin-up passes), not the model time; a restart file holds one slice; the folder is made with `mkdir -p` inside | `vilma_io.f90:67-81, 419-473`; climber-x `src/geo/vilma2.F90:66-67, 158` | write with `se%time`, read the last slice (as CLIMBER-X); no time argument |
| F9 | Measured cost (Levante, `ifx -Ofast`): lmax 128, 3D Bagge, 261 coupling steps of 100 yr: solver 1585 s (6.1 s per call), setup 172 s, 10 GB; fitted exponent in `lmax` 2.40 (3D) | `vilma2 doc/performance.md` | T170 estimate ≈ 12 s per call, ≈ 20 GB (§14) |
| F10 | `config/common.mk:179-181` stops on `vilma=1`; the `vilma` switch and the stamp already exist (`:103, 107`) | `config/common.mk` | WP2 replaces the stop by the real wiring |
| F11 | The cold start runs the geo state before the ice state (`run_cold_start` block a), `kryos_run_sequence.f90:56-72`), but the t0 load is final only after `PASS_STATE` (recon, equil) | `kryos_run_sequence.f90:74-112` | two-pass cold start (§8.1, Q8) |
| F12 | Yelmo `H_ice` is the cell mean (`V_ice = Σ H_ice·area`, `yelmo_regions.f90:365`), so `I_g = H_ice·f_grnd` needs no `f_ice` factor | yelmo | as in yx-plan |
| F13 | VILMA2's writers (`vilma_io.f90`: `put2d`, `put_hor`, `put3d`, `put_scalar`, `put_rotmem`) call `nc_write` with units and long name only: no `prec_abs`/`prec_rel`; `input/vilma-variables.md` has 5 columns (fesm-utils `variable_io` accepts a 6th `prec`, lossless without); restart (`vilma_restart_write`) and diagnostics (`vilma_write_step`, `vilma_write_horizontal`) share the same `put*` | `vilma2 src/vilma_io.f90:195-236, 367-438`; `fesm-utils src/variable_io.f90:20-21, 170-237`; `ncio.f90:529-560` | VILMA2 cannot round today; Q13, §9.1 |

## 2. Architecture

### 2.1 Instances (science case)

| Instance | Scope | Component | Backend | Grid | Role |
|---|---|---|---|---|---|
| `geo_glob` | global | geo | `vilma2` (`vilma_earth_t`) | `GAUSS-L170` | the Earth: solve, exports `dz_bed`, `dz_ss`, `dphi`, `bsl`, `bsl_ive`, diagnostics |
| `geo_nh`, `geo_ant` | regional | geo | `vilma2` (`vilma_geo_t`, the port) | default chain: the ice grid | holds the split anomalies; feeds `geo → ice` |
| `ice_nh`, `ice_ant` | regional | ice | yelmo | NH-32KM; ANT-32KM (hub ANT-16KM, as in the shipped base) | loads (gather source) and consumers |
| `forcing` | global | forcing | tsgen | — | shared by both domains (phase 3); independent of phase 4 |
| `[globe_ice]` | globe | — (a par group read by the globe, Q5) | datasrc | `GLOBAL-0.5DEG` | background X and the external history |

### 2.2 Modules and files

New:

| File | Content | Size |
|---|---|---|
| `src/backends/geo/kryos_geo_vilma2.f90`, `…_stub.f90` | `vilma_earth_t`, `vilma_geo_t`, `vilma2_*_constraints()`; stub with the same names, `vilma=0` | M + S |
| `src/core/kryos_globe_ice.f90` | `globe_ice_t`, the global geometry hub (all phase 4, §10): global grids (background, Gauss) and their coupler, per-domain footprints and coverage, background and external ice (datasrc), reference assembly (spliced and VILMA2's), the per-domain host terms, merges, checks, diagnostics, `globe_ice.nc`. Knows no domain type: callers pass arrays and the coupler | L |
| `src/core/kryos_run_gia.f90` | submodule of `kryos_run`: `gia_reference`, `gia_gather`, `gia_seed` (the scatter is the rows `geo_g → geo`, `bsl_to_geo`; loops over `run%dom` in `run%dom_order`, as `reduce_bsl`, `gather_forcing` in `kryos_run_sequence.f90:338-388`) | M |
| `tests/test_globe_ice.f90`, `tests/test_geo_vilma2.f90`, `tests/test_gia_sequence.f90` | §13 | M |

Changed: `kryos_geo.f90` (opts, declaration), `kryos_geo_factory.f90` (backend `vilma2`
by scope), `kryos_ice.f90` + `kryos_ice_yelmo.f90` (exports, flotation binding),
`kryos_globe.f90` (`globe%cpl`, grid, rank 2), `kryos_wiring.f90` (rows, `T_GIA`, `dphi`),
`kryos_run.f90` / `kryos_run_sequence.f90` (§7, §8), `kryos_restart.f90` (manifest `grid`
of a global instance; `GLOBE_ORDER` stays `forcing, geo`), `kryos_output.f90` /
`kryos_outfile.f90` (lon-lat writer with rounding, §9.1), `kryos_config.f90` templates,
`input/kryos_defaults.nml`, `input/kryos-variables-geo.md` (+ the `prec` column for the
new variables), `input/sync.txt` (+ `vilma_defaults.nml` from `VILMAROOT`),
`config/common.mk`, `config/Makefile`, `.configme/manifest.toml`.

### 2.3 Vocabulary

**geo, global part** (`vilma_earth_t`; rank 2 on `GAUSS-L<lmax>`, the weighted
fields per unit Gauss-cell area as in `ht` §1):

| Import | Units | Meaning | Filled by |
|---|---|---|---|
| `H_ice_ext` | m | external ice, mean over the `1 − w` part | `globe_ice` |
| `host_w` | 1 | Σ over domains of the coverage of the host region | gather |
| `host_Ig` | m | `w · H · f_g` | gather |
| `host_C` | 1 | `w · (1 − f_g)` | gather |
| `host_F` | m | `w · ζ⁽⁰⁾ (C − C_ref)` | gather |
| `z_bed_ref`, `H_ice_ref` (static) | m | background X on the Gauss grid (the `1 − w` part) | `globe_ice` |
| `host_Ig_ref`, `host_C_ref` (static) | m, 1 | reference host terms | gather, once |

| Export | Units | Meaning |
|---|---|---|
| `dz_bed` (rank 2) | m | `u`, the bed change |
| `dz_ss` (rank 2) | m | `N`, the sea-surface change **relative to `dphi`** (VILMA2's `dz_ss − Δφ`) |
| `C_ocn` (rank 2) | 1 | total ocean function (diagnostic) |
| `dphi` | m | the uniform sea-surface shift `Δφ` |
| `bsl` | m | ocean mean of `C · rsl` |
| `bsl_ive` | m | ice-volume equivalent (diagnostic) |
| `bsl_ive_host` | m | the host part of `bsl_ive` (diagnostic; 0 without host terms) |
| `A_ocn` | m2 | `∫ C dA` (diagnostic) |
| `mass_resid` | 1 | worst SLE mass residual of the step |

`rsl = dz_ss + dphi − dz_bed` is a diagnostic of the output, not an export.

**geo, regional port** (`vilma_geo_t`): imports `dz_bed`, `dz_ss`, `dphi`, `bsl`
(from the globe) and the static `z_bed_ref`, `H_ice_ref` (from the ice, existing
rows); exports `dz_bed`, `dz_ss`, `dphi`, `bsl`, `z_bed_ref` on the geo grid. It has
no `H_ice` import (the existing row `ice.H_ice → geo.H_ice` is then inactive), no
`dV_bsl`, so `reduce_bsl` finds no gather source and returns.
`z_bed_ref` is the ice's, passed through, so `domain_check_reference` still compares
it with `refs(geo grid)`.

**ice additions** (contract, optional; the yelmo adapter sets them; also needed by
phase 2 WP3 except `footprint`): static exports `footprint` (1 where Yelmo holds
ice, `mask_ice /= MASK_ICE_NONE` at the end of init) and `H_ice_pd` (`dta%pd%H_ice`);
and a binding `flotation(H, z_bed, z_sl) → f_grnd` (capability `has_flotation`, wrapper
`ice_flotation(ice, …)` like the other wrappers; decided, Q2), implemented by the yelmo
adapter with Yelmo's own routines: `calc_H_grnd` (no `f_ice`) and
`calc_f_grnd_subgrid_area` or `_linear` per the run's `gl_sep`, on the instance grid with
the run's `rho_ice`, `rho_sw`.

### 2.4 Wiring

Added or changed rows (`kryos_wiring.f90:100-141`, one table, fixed in code):

| src | field | dst | field | method | lag | transform |
|---|---|---|---|---|---|---|
| geo_g | `dz_bed`, `dz_ss` | geo | same | bilin (anomalies) | 0 | — (rank 2, remapped in the domain's coupler) |
| geo_g | `dphi` | geo | `dphi` | — | 0 | — |
| geo_g | `bsl` | geo | `bsl` | — | 0 | existing row |
| geo | `dz_ss` | ice | `z_sl` | auto | 0 | `ADD_SCALAR`, aux **`dphi`** (was `bsl`, Q1) |
| ice | `H_ice`, `f_grnd`, `z_bed`, `z_sl` | geo_g | `host_w`, `host_Ig`, `host_C`, `host_F` | — | 1 | `T_GIA`, gather |
| ice | `z_bed_ref`, `H_ice_pd`, `footprint` | geo_g | `host_Ig_ref`, `host_C_ref` | — | `LAG_INIT` | `T_GIA`, gather |

`T_GIA` rows exist so that `wire_resolve` finds a producer for every required import
(gaps stop with a name), logs them, and checks the lags; `wire_exchange` never runs them
(it stops if asked): the arithmetic is `gia_gather` (§4), because it needs flotation,
the external ice and the merge, which are not a remap plus a transform (D3). The
imports `H_ice_ext`, `z_bed_ref`, `H_ice_ref` of the global geo have no component
producer; `wire_resolve` skips the gap check for the imports of a global geo that
declares `opt%host` and logs "filled by gia_gather / globe_ice".
The rank-2 globe rows use the domain's coupler (`dom%cpl`): the Gauss grid is
registered there (`coupler_add_grid`) and the `GAUSS → geo grid` bilinear maps are
primed by `wire_prime` like any other pair.

## 3. The `vilma2` geo backend

### 3.1 Adapter

`vilma_earth_t extends geo_t` (scope global; `prec = dp`: VILMA2 is `real64`, so no
narrowing anywhere) holds `type(solid_earth) :: se`. Lifecycle (`core` §3, adapters §9):

| Binding | Action |
|---|---|
| `init(ctx)` | read `[geo_glob]` (keys §3.4); `vilma_par_load(se%par, par_file, defaults_file = "input/vilma_defaults.nml", group = <vilma>)` (as climber-x `vilma2.F90:54-56`); **set `se%par%rho_ice`, `rho_water` from `ctx%cnst` (`rho_ice`, `rho_sw`)** and log both pairs; `phys_const_compare` for `r_earth`, `g`; require `solver = "v2"` and `sle_fixed_ocean = .false.` (`scheme = "fe"` is needed only by change 5); build the Gauss grid (§3.5); declare fields |
| `init_state(PASS_INIT)` | nothing put |
| `init_state(PASS_GEOMETRY)` | after `gia_reference`: `solid_earth_init(se, z_bed_eq = z_bed_ref, h_ice_eq = H_ice_ref, grid = gauss grid, w = host_w, I_g = host_Ig_ref, C = host_C_ref)` (no `h_ice_init`: `ht` Q5); stop if VILMA2 built a remap (`se%remap`), i.e. the grid was not recognised as its Gauss grid; put the exports at the reference (`dz_bed = dz_ss = dphi = bsl = 0`, `C_ocn` = reference) so the domains find them (§8.1) |
| `init_state(PASS_STATE)` (cold only) | `solid_earth_spinup(se, H_ice_ext, w, I_g, C, F)` (a no-op unless `equil_time_max > 0` or `pre_spinup_1d`), then `solid_earth_update(se, H_ice_ext, 0, w, I_g, C, F)` (the seed), put the exports |
| `update(st)` | `dt = st%t1 − st%t0` [yr, `dp`]; `solid_earth_update(se, H_ice_ext, dt, host_w, host_Ig, host_C, host_F)`; put `dz_bed = se%dz_bed`, `dz_ss = se%dz_ss − se%dphi`, `C_ocn = se%C` (VILMA2's host-grid outputs, never `se%gg`: the Gauss state is north first, the grid ascends), `dphi`, `bsl`, `bsl_ive`, `A_ocn` (∫C over VILMA2's sphere), `mass_resid` at `st%t1` |
| `restart_write(dir, t)` | `vilma_restart_write(se, se%time, filename = "<name>_restart.nc", folder = dir)` |
| `restart_read(dir, t)` | `vilma_restart_read(se, file)` into the model built by `PASS_GEOMETRY` (last slice) |
| `write_output(stream)` | default writer for the exports; backend streams `gia` and `vilma` (§9) |
| `end` | `solid_earth_finalize(se)` |

`PASS_RESTART` puts nothing (adapters A5): the exports are the bundle's, and
`vilma_restart_read` ends with `solid_earth_sync_host`, so the Gauss state and the
exports agree.

`vilma_geo_t extends geo_t` (scope regional; `con%restart = .false.`): `init` takes its
cadence from the global (a `dt` key that differs stops); `init_state(PASS_GEOMETRY)` puts
the reference exports; `init_state(PASS_STATE)` and `update(st)` do the same thing: copy
the imports `dz_bed`, `dz_ss`, `dphi`, `bsl` to the exports with their validity, and
`z_bed_ref` from its import. The core fills the imports with the scatter rows just
before (§6).

### 3.2 Stub, factory, constraints

`kryos_geo_vilma2_stub.f90` keeps the module and type names, `linked = .false.`,
bindings stop "rebuild with `vilma=1`" (as `kryos_geo_fastisostasy_stub.f90`). The geo
factory (`kryos_geo_factory.f90:20-35`) gets `case ("vilma2")` choosing `vilma_earth_t`
for scope `global`, `vilma_geo_t` for `regional`. Constraints: global
`constraints_t(build_flag="vilma=1", regional=.false., global=.true., max_instances=1)`;
regional `build_flag="vilma=1"`, `restart=.false.`. Init-time contract checks (named
stops): a domain geo with `vilma2` needs the global geo with `vilma2` (and the reverse,
until §18); the global geo needs every domain's `ice` to declare `has_flotation` and the
`footprint`, `H_ice_pd` exports; footprints of the domains must not overlap (merge
check); `fastisostasy` + `vilma2` together stops with "not implemented (phase 4 later
package)".

### 3.3 Build

`config/common.mk`: remove the stop at `:179-181`; add `VILMAROOT = vilma2`,
`INC_VILMA = -I${VILMAROOT}/libvilma/include`, `LIB_VILMA = -L… -lvilma`,
`VILMA_ARCHIVE = ${VILMAROOT}/libvilma/include/libvilma.a`, and SHTns from fesm-utils:
`LIB_SHTNS = -L${FESMUTILSROOT}/SHTns/shtns-{serial,omp}/lib -lshtns[_omp]` chosen with
`openmp` like the LIS/FFTW lines and appended to `LIB_EXTERNAL` before FFTW (`LFLAGS` is
`LIB_BACKENDS LIB_FESMUTILS LIB_NC LIB_EXTERNAL`, `config/common.mk:203`: backends first,
as VILMA2 needs). The `vilma=1` branch mirrors the `fastisostasy=1`
branch (`:109-122`): `kryos_geo_vilma2_obj`, `INC_BACKENDS`, `LIB_BACKENDS` with `-lvilma`
(already before `-lfesmutils` in `LFLAGS`), `dep_archives
+= $(VILMA_ARCHIVE)` (the file target, so kryos objects rebuild when libvilma changes),
`backend_roots`, `build_info_libs += vilma:$(VILMAROOT)`, `backend_static += vilma-static`,
`need_lis_fftw = 1`, `test_skip` for the vilma tests when off. `config/Makefile` gets
`vilma-static: fesmutils-static` → `$(MAKE) -C $(VILMAROOT) vilma-static openmp=$(openmp)`
(SHTns is the one-time `make shtns` of fesm-utils, like its `libs`; `deps-check` names it when missing). `input/sync.txt` gets
`vilma_defaults.nml VILMAROOT input/vilma_defaults.nml`. `.configme/manifest.toml` gets
`vilma2:main` (decided, Q10: main d8041e8 has changes 1-4, the API is the README "Host-resolved
terms" and `doc/design.md` §15; no merge to wait for); `vilma2`'s own manifest pins fesm-utils
`phase0-kryos`, as kryos. `openmp` must match across kryos, fesm-utils, vilma2.

### 3.4 Configuration

```fortran
&globe        geo = "geo_glob"   forcing = "forcing" /
&geo_glob                              ! the global part
    backend  = "vilma2"
    dt       = 10.0                    ! [yr] geo cadence (Q6, revisit after G7); clock = "model"
    allow_no_spinup = .false.          ! a run without VILMA2's spin-up (pre_spinup_1d, equil_time_max > 0) stops
                                       ! (gia_check) unless .true.: the seed solves the previous load, so the t0 load
                                       ! enters with the first call (F7); tests and present-day starts only
    vilma    = "vilma"                 ! the &vilma group; defaults input/vilma_defaults.nml
    ice_bg   = "globe_ice"             ! the global ice field group
/
&vilma                                 ! overrides only, as climber-x geo_par.nml:124-145
    lmax = 170   earth = "PREM"   n_layer = 11   r_bot = …   r_top = …   eta = …   rheology = …
    l_visc_3d = .true.   visc_3d_file = "fesmdata/Global/GLOBAL-T170/GLOBAL-T170_SolidEarth-Bagge2021.nc"
    visc_3d_structure = 6
    pre_spinup_1d = .true.   equil_time_max = …          ! Q7: VILMA2's own spin-up; value set after G2c
/
&globe_ice
    grid      = "GLOBAL-0.5DEG"        ! background grid (Q3, decided)
    z_bed_ref = "X_z_bed"   H_ice_ref = "X_H_ice"       ! datasrc field groups (below)
    ext       = "file"                 ! "ref" | "file": ice outside the footprints
    ext_field = "X_dH_ice"             ! datasrc group of the history (anomaly)
    ext_tref  = 0.0                    ! [yr] model time of the history's present day
    out_2D    = "none"                 ! schedule of globe_ice.nc (hub output, §9)
/
&geo_nh   backend = "vilma2"   /       ! the ports; grid follows the chain (ice grid)
&geo_ant  backend = "vilma2"   /
```

As built (WP4/WP5): `[geo_vilma2]` also has `out_gia` and `out_vilma` (schedules of `<inst>_gia.nc`, the host imports as VILMA2 got them, and
`<inst>_vilma.nc`, VILMA2's own diagnostics) and `allow_no_spinup`; `[globe_ice]` has `target` (the global grid of the host terms for the
observer without VILMA2; with `geo = vilma2` it is the Gauss grid). The Earth's 2D files are written on `lon`, `lat`.
The densities are not keys: they are the run's (`phys_const`). `[globe_ice]` follows a new
template `[globe_ice]`; the three field groups follow `datasrc_defaults.nml`:

- `X_z_bed`, `X_H_ice`: `{grid_name}_IceHistory-PaleoMIST.nc`, variables `z_bed`, `H_ice`,
  `mode = "transient"` evaluated once at `ext_tref` (a slice of the history; if
  `fesm_datasrc` has no static-slice read, this is it);
- `X_dH_ice`: the same file, variable `H_ice`, `mode = "transient"`, `anomaly_ref = "slice"`
  with the reference time `ext_tref`, `out_of_range = "hold"` after the last time and an
  error before the first (the yx-plan rule; `fesm_datasrc` serves exactly this,
  `fesm_datasrc.f90:1-35`).

Templates for the instance groups: `[geo_vilma2]` in `input/kryos_defaults.nml` with `vilma`,
`ice_bg`; the regional port has none.

### 3.5 Gauss grid

At `geo_glob` init: `grid_init_gaussian(g, "GAUSS-L<lmax>", nlon = nphi, nlat = nlat,
lon180 = .false., x0 = 0)` with `nlat = 2 lmax + 2`,
`nphi = 4 lmax` (or the `nlat`, `nphi` of `[vilma]` when set), the same numbers as
`build_sht` (`vilma_coupling.f90:690-708`); `planet` is left at its default, "Spherical Earth"
(fesm-utils 5662899), so the cell areas of the grid are the quadrature areas `R²·dλ·w_j`;
`grid_cdo_write` of its description to
`maps/grid_GAUSS-L<lmax>.txt` (explicit `yvals`, so any coupler resolves it by name);
`coupler_add_grid` in the globe's coupler: its default model area (`grid%area·xy_conv²`) is
already the quadrature area, so no `coupler_set_area` for this grid (the adapter checks the
planet radius of the grid against `se%par%r_earth` and stops on a mismatch). `gauss_axes_rows` (VILMA2, 1e-6° on
every axis value) is the check that kryos's grid is VILMA2's: the adapter stops if
`se%remap` is true after init. Longitudes start at 0 (SHTns). Latitudes ascend; VILMA2's
state is north first, so the passthrough reverses the rows both ways (`se%host_rows`).

## 4. The global ice field (gather)

Run by `gia_gather(run)` when the global geo is due, at the start of the step (§7), on the
exports the ice published at the previous step (lag 1): `wire_resolve`'s lag check
(`t1 == cadence_t_last`) is applied to every source field.

### 4.1 Per step

1. **External ice on the background grid** (`globe_ice`): `H_ext(t_load) = max(0, H_X +
   ΔH_hist(t_load))` for `ext = "file"` (`ΔH_hist` from the `anomaly_ref = slice` field; the
   clipped mass is logged), `H_X` for `"ref"`. Evaluated at every geo step: no separate
   cadence, so the restart is exact (D7). With the default (X = the history's own slice at
   `ext_tref`) `H_ext(t) = H_hist(t)` exactly.
2. `H_ice_ext` on the Gauss grid: `H_ext` remapped `con` (extensive) from the background
   grid. In partial cells (`w < 1`) this is the whole-cell mean applied to the `1 − w` part
   (D17; masked variant later).
3. **Per domain d**, name order, on its ice grid `i` (one pass over the grid):
   - `ext_i` = `H_ext` remapped (`con`) from the background grid to grid `i`;
   - `m` = ice `footprint`; `H`, `f_Y` = ice `H_ice`, `f_grnd`; `z_bed`, `z_sl` = the ice's
     exports (the bed and sea level Yelmo used);
   - `H_tot = m·H + (1 − m)·ext_i`;
   - `f_g = m·f_Y + (1 − m)·flot(H_tot, z_bed, z_sl)`, `flot` = `ice_flotation` (Yelmo's rule,
     §4.3), evaluated once on the whole grid;
   - `I = H_tot·f_g`, `Cc = 1 − f_g`, `F = ζ⁽⁰⁾·(f_g,ref − f_g)` with `ζ⁽⁰⁾ = z_bed_ref − z_sl_ref`
     (`z_sl_ref = 0`) and `f_g,ref` from the reference (§5.1);
   - `coupler_merge_add` into three accumulators: `I` (`extensive`), `Cc` (`intensive`: the
     amount is `coverage × mean`, which is `ht`'s `w·(1 − f_g)`), `F` (`extensive`).
4. `coupler_merge_end` with background 0 for each: `host_Ig`, `host_C`, `host_F`, and
   `host_w = frac_tot` (the summed coverage; overlap above `overlap_tol` stops naming the
   domains). The result is `ht`'s form (weighted fields per unit target area) with no
   further step.
5. Checks and logs (§4.5), then the globe geo's imports are put (validity `[t_load, t_load]`).

One ice cadence across the domains (decided at the WP4 checkpoint): the load time is the ice's last update
(`cadence_t_last`, lag 1; the cold start: t0), and `gia_gather` stops when two domains' ice instances were last
updated at different times: the Earth has one load time.

### 4.2 Footprint

`footprint` is static (published at `PASS_INIT`). `gia_gather` compares it to the first
gathered one (and `host_w` to the reference `w`) every step: a change stops naming the
domain (settled: the footprint is fixed). If the recon
recipe (phase 2 WP4) can change `mask_ice`, the first check after the cold start catches it.

### 4.3 Flotation

Outside the footprint (and in the reference) the rule must be Yelmo's, so it is the ice
instance's binding `flotation` (Q2): kryos calls the wrapper `ice_flotation(dom%ice, H, z_bed,
z_sl, f_g)`; the core never sees Yelmo. Test (§13): on Yelmo's own geometry the binding
reproduces `tpo%f_grnd` bit for bit where Yelmo uses no pinning points (otherwise the
differences are listed); this makes a footprint cell and a non-footprint cell consistent.
Densities are Yelmo's `bnd%c` (the run's), which are the ones VILMA2 got (§3.1).

### 4.4 Source areas

The extensive merges use the source grid's model area. The coupler default for a projected
grid is the plane `dx·dy` (what Yelmo integrates), which differs from the true area by up to
about 5% on the polar stereographic grids (yx-plan). D6: `coupler_set_area` of each ice grid
in the globe's coupler is set to the true (ellipsoidal) cell area `grid%area · xy_conv²`, so
the load per true area is right and the global integrals are spherical (the Gauss grid needs no
call: its default area is the quadrature area, §3.5); Yelmo's `M_ice`
differs from the Gauss-grid grounded ice by that factor and the mass check (§4.5) reports
both. Decided (Q11).

Metrics at the poles (decided, WP3 checkpoint): the amounts (`host_Ig`, `host_F`, the merged ice)
are per true area, so the merges conserve Yelmo's mass; the weights `w` and `C` are the coupler's
geometric coverage on the sphere (map metric). Near the poles an ellipsoidal cell is up to 0.9 %
larger than the sphere's, so in rim cells (0 < w < 1) the amount per area and `w` times the cell
mean differ by up to that much (measured on ANT-32KM and GRL-16KM against GAUSS-L170: 0.9 %;
spherical grids: 8e-5); in the interior they agree to round-off, and the merge closes to 4e-14.

### 4.5 Checks and logs (every gather)

- finite; `0 ≤ host_C ≤ host_w + tol`; `host_Ig ≥ 0`; `host_w ≤ 1 + tol` (VILMA2's
  `solid_earth_host_check` repeats it, with its own 1e-9; it stays as it is). **Producer side**: the merge makes
  `host_w` and `host_C` from conservative maps, whose weights are exact to the round-off of their polygon areas
  (`host_w` = 1 + 2.8e-9 at the polar row of GAUSS-L170), not of the machine; the gather therefore clips `host_w` to
  [0, 1] and `host_C` to [0, `host_w`] after the merge, logs the largest excess (`kryos_globe_ice`, every gather and
  the reference), and stops if an excess exceeds `HOST_CLIP_TOL` = 1e-6 (named parameter in `kryos_globe_ice.f90`: that
  is not round-off). The conservation of the hub's merged reference and ice is held to `HUB_TOL` = 1e-9 (observed
  2.6e-10 reference, 1.5e-10 ice), the same bound as the closure of a gather;
- conservation: `∫ host_Ig dA` on the Gauss grid (`coupler_integral`) against `Σ_d ∫ m·H·f_g dA`
  and the external part on the source grids, relative difference ≤ 1e-9 (merge exactness);
- per domain and external: grounded ice volume and its sea-level equivalent
  `−(ρ_i/ρ_w) ΔV_g / A_ocn` against the reference (Gt and m ESL), for the log and
  `geo_glob_ts.nc`; the clipped external mass;
- at the reference: `A_ocn_ref = ∫ C_tot,ref dA` against `cnst%area_seasurf` (1% logged).

## 5. The reference

### 5.1 Assembly (once, `gia_reference`)

After every instance's `PASS_INIT` and the static rows (`run_pass_init`,
`kryos_run.f90:438-454`), before the cold start and before any restore:

1. Background X on the background grid: `z_bed_X`, `H_X` (datasrc at `ext_tref`), remapped to
   the Gauss grid (`con`): `z_bed_ref`, `H_ice_ref` of the global geo (the `1 − w` part,
   what VILMA2's `z_bed_eq`, `h_ice_eq` mean with host terms, `ht` §1, §4).
2. Per domain on its ice grid: `H_ref,tot = m·H_ice_pd + (1 − m)·ext_ref,i` (`ext_ref,i` = X ice
   remapped `con` to the grid), `z_bed = z_bed_ref` (the ice's static export = `dta%pd%z_bed`) in
   every host cell, `f_g,ref = flot(H_ref,tot, z_bed_ref, 0)`, `I_ref = H_ref,tot f_g,ref`,
   `C_ref = 1 − f_g,ref`; merged exactly as in §4.1 (the same procedure with `H_ice_pd`,
   `z_bed_ref`, `z_sl = 0` as inputs, so reference and step cannot drift apart): `host_w`,
   `host_Ig_ref`, `host_C_ref`; `F_ref ≡ 0` (`ht` Q2).
3. The **spliced global reference topography** (`z_bed`: intensive merge; `H_ice`: extensive
   merge, background X) for the global geometry hub (§5.4) and diagnostics:
   `(1 − Σf_d)·X + Σ f_d·ref_d`, with `ref_d` the domain's reference of item 2 (`H_ref,tot` and
   `z_bed_ref` on its whole grid: `H_ice_pd` inside the footprint, X elsewhere; decided, WP3
   checkpoint: one reference definition for VILMA2 and the hub, D18). `design.md` §6 item 1 says
   "present day" for `ref_d`; outside the footprint that is X here. VILMA2 itself gets the
   background part only (item 1): giving it the spliced ice would count the domain's present
   day twice in partial cells.

Consistency self-test (§13): at the reference state (`H = H_ice_pd`, `z_bed = z_bed_ref`,
`z_sl = 0`) one gather returns `ΔI_g = 0`, `ΔC = 0`, `F = 0` to round-off.

### 5.2 Which X

| Option | Fields | Exists | Comment |
|---|---|---|---|
| **IceHistory-PaleoMIST at 0 ka**, `GLOBAL-0.5DEG` (**decided**, Q3) | `z_bed`, `H_ice` | yes (0.25, 0.5, 1, 2, 5°) | same file as the history: `H_ext(t_ref) = X` by construction; history back to 80 ka (Q4) |
| same on `GLOBAL-0.25DEG` | | yes | native resolution, 4× data; not chosen |
| ICE6G_C or GLAC1D at 0 ka | `z_topo` (+ `H_ice`) | yes (0.5°, 1°…) | not chosen; `z_topo` includes ice base under shelves, not `z_bed`: a conversion would be needed |
| a global Topo product (BedMachine + GEBCO) | `z_bed`, `H_ice` | **no** (Topo is domain-only) | settled: not a FesmData product |

### 5.3 Consistency with `ice_ext`

Outside the footprints the load is `X + [H_hist(t) − H_hist(t_ref)]` (clipped): a mismatch
between the history's present day and X cannot create a load. With the default, X is the
history's slice, so the anomaly is the history itself. Inside the host region but outside
the footprint the same external field is interpolated to the ice grid (§4.1 step 3), and in
the reference the external part is X. The ice of the footprint is Yelmo's at every time,
whatever the history says there.

### 5.4 The global geometry hub

All of it is phase 4 (phase 3 delivers two domains only, §10); it is `globe_ice_t`
(`kryos_globe_ice.f90`), so one code path serves the hub and VILMA2: the global grids
(background `GLOBAL-0.5DEG`, `GAUSS-L<lmax>`) and `globe%cpl`, the footprint and coverage
`f_d` of every domain on each global grid (owner weights `o_d` for overlapping domains,
`coupler_owner_reduce`; the science case has none, overlap above `overlap_tol` stops), the
model areas (§4.4), the assembly of §5.1 and the merge of the domains' ice of §4.1, and the
checks of §4.5. With `geo = vilma2` the hub is a by-product of `gia_reference`/`gia_gather`;
with any other geo (e.g. `fastisostasy`) it runs **observer-only** when `[globe_ice]` is
present: it assembles and writes, nothing feeds back (the gate G1b). Output: `globe_ice.nc`
on the background grid at `[globe_ice] out_2D`: the spliced reference (`z_bed`, `H_ice`), the
merged ice of the step (lag 1, as `ice → geo`), the coverage `Σ f_d`; rounded (§9.1). No
restart (derived; rebuilt like the domain hubs). No `[globe_ice]` group: no cost.

As built (WP3): `globe_ice_add_domain(..., priority)` gives overlapping domains owner weights
(descending priority, ties by ice grid name; pairs that share no global cell are skipped, so no
ice-to-ice map is needed for disjoint domains); domains without priority must not overlap (the merge
stops naming them). The file variables are `z_bed_ref`, `H_ice_ref` (the spliced reference, static),
`cov` (owned coverage, at most 1), `H_ice` and `t_load` (one record per output); the table is
`input/kryos-variables-globe_ice.md`. The hub's merged reference and ice are checked against the
background outside the domains plus the owned ice of the domains (1e-9, logged at the reference and at
every write); `tools/bin/kryos-globe-ice.jl` recomputes the merge offline from the coupler's map files.

## 6. Scatter

The scatter after the global update, for every domain (as built: `bsl_to_geo` runs the resolved rows, then
`domain_update`; no separate procedure):

1. rows `geo_g → geo` into the port's imports: `dz_bed`, `dz_ss` (`bilin` remap through the domain's coupler,
   in which the Gauss grid is registered; no references: anomalies stay anomalies), `dphi`, `bsl`;
2. `domain_update(dom, C_GEO)` (the port copies imports to exports at `t1`);
3. the existing rows (`kryos_wiring.f90:114-116`) `geo.dz_bed → ice.z_bed` (`ADD_REF`),
   `geo.dz_ss → ice.z_sl` (`ADD_SCALAR` with `dphi`), run in `exchange_into(dom, "ice")` at
   step 4 as today.

The geo grid of the port defaults to the ice grid (`GRID_DEFAULT`: geo ← ice), so each field
is remapped once (Gauss → ice grid); the hub follows through the ice as today
(`domain_hub_update`): the anomalies reach the hub resolution through the ice's `z_bed`, which
is what design D4 means by "split onto each hub" in kryos terms (design.md change, §11). A
port on the hub grid (`[geo_ant] grid = "ANT-16KM"`) is a par-file choice, at one more remap
(geo → ice `auto`). `bsl` goes to every domain (diagnostics, consumers of the barystatic level).
The domain's `bsl` export is VILMA2's ocean mean; the sea level the ice sees is
`dphi + dz_ss`.

## 7. Time stepping

### 7.1 Step order (`run_sequence`, `kryos_run_sequence.f90:219-299`)

```
1. globe: [new] if the global geo needs a gather (opt%host) and is due: gia_gather(run)   ! ice at t_{n-1}
          globe_update(geo)                       ! fastiso: bsl_update; vilma2: solid_earth_update
          forcing as today
3. geo:   rows geo_g -> geo (bsl_to_geo, rank 2 for a port), port update, when due
          reduce_bsl: no source for vilma2, returns
4. per domain: unchanged (exchange_into(ice) takes z_bed, z_sl of this step; lag 0)
```

For the fastiso runs nothing changes (the gather and scatter are skipped when the geo has no
`opt%host`), which the phase 1-3 gates verify bit for bit (G1).

### 7.2 Cadence and lag

VILMA2 has no cadence of its own; the adapter passes `dt = st%t1 − st%t0` of the instance's
cadence interval. `geo_glob` `dt = 10` (as FastIsostasy `dt_isos` and CLIMBER-X
`n_year_geo`, `nml/control.nml:52`), a multiple of `dtt`, and restarts only on cadence times
(`clock_start`). Lag 1: the load of the call at step `n` is the ice of step `n − 1`, so
VILMA2's load history is one geo interval behind the ice, like FastIsostasy's. VILMA2 ramps
linearly from the previous call's load to this one (`vilma_coupling.f90:369-467`, host terms
likewise, `ht` §3). The cost of the cadence is §14 (Q6).

### 7.3 `bsl` and `z_sl`

`z_sl` at the ice = `dphi + remap(dz_ss)`, lag 0 (this step's solve), as in phase 1. `bsl`
is a diagnostic and a scattered scalar.

### 7.4 Spin-up and the seed

`solid_earth_spinup` relaxes VILMA2's memory under the t0 load while holding the reference as
datum (`vilma_coupling.f90:492-554`; `pre_spinup_1d`, `equil_time_max`, `equil_rate_tol` of
`[vilma]`); host terms are held over all passes (`ht` §1). It consumes VILMA2's clock
(`se%time`) but not the model's. The seed `solid_earth_update(…, dt = 0, …)` follows it with
the same load: F7 (the `dt = 0` branch solves with `ice0`, i.e. the spin-up load). Without a
spin-up the seed solves the reference and the first call ramps the whole t0 load over one
interval: only for tests and present-day starts.

Decided (Q7): VILMA2's own spin-up, holding the t0 load (`pre_spinup_1d = .true.`); the
value of `equil_time_max` (and `equil_rate_tol`) stays a parameter, set after a test (gate
G2c, §13.2): the test varies `equil_time_max` on the 80 ka start and reports the bed rate
at the seed and the drift of the first intervals.

## 8. Init and restart

### 8.1 Cold start (two-pass)

`run_init`: the globe is initialised **before** the domains (`kryos_run.f90:134`), so
`geo_glob` cannot build VILMA2 at `init`; it does at `PASS_GEOMETRY`.

```
run_pass_init:       PASS_INIT of globe, then of each domain; static rows      (as today)
[new] gia_reference(run)                               §5.1: statics -> globe geo imports (w, X, refs)
[new] globe_pass(geo, PASS_GEOMETRY)                   solid_earth_init at the reference; exports 0
run_cold_start, per domain, block a):
      exchange(ice -> geo, init)                       no active rows for the port
      [changed] globe_to_geo rows (init)               dz = dphi = bsl = 0 into the ports   (was bsl_to_geo)
      domain_pass(geo, PASS_GEOMETRY); domain_check_reference
   reduce_bsl                                          returns
per domain, blocks b), c):                             as today, on the reference bed
[new] gia_seed(run), after c) of every domain          (c) includes the cold-start hooks and the ice PASS_SPINUP):
                                                       gather at t0 (ice final); globe_pass(geo, PASS_STATE)
                                                       = spin-up + seed; rows geo_g -> geo; ports PASS_STATE
per domain, [new] block c'):                           exchange(geo -> ice); domain_pass(ice, PASS_GEOMETRY)
                                                       (bnd bed/sea level, calc_ytopo_diagnostic); hub update
per domain, block d):                                  surface PASS_FINAL on the final geometry; as today
```

The domain loop of `run_cold_start` (`kryos_run_sequence.f90:74-112`) is split after block c)
(a second loop with c') and d)); per domain the order is unchanged, so N = 1 fastiso runs
are bit-identical. Why this order (decided, Q8): the t0 load is final only after `PASS_STATE`
(recon, equil) and the ice `PASS_SPINUP`; yx-plan step 3 ("after all domains") and its
reconcile; moving block d) after the reconcile makes the surface's final cycle see the final
geometry (the yx-plan did not).
Ice thickness is kept fixed through the reconcile; `z_srf` changes by the bed shift.

### 8.2 Restore

`run_restore` (`kryos_run_sequence.f90:119-189`): after the clock, cadences and manifest,
the globe step 3 becomes: `gia_reference`; `globe_pass(geo, PASS_GEOMETRY)` (VILMA2 built at
the reference, memory zero); `restart_read_global(geo)` (`vilma_restart_read` of the bundle
file, which `check_reference`s the static fields against the model just built, tolerance 1e-3,
and restores the memory, `σ_n`, rotation, `dt_try`, the ramp's last host terms, `se%time`);
the exports and imports are the bundle's (`fields_read`). Domains in `RESTORE_ORDER`: the
ports are stateless (their exports restored). Nothing is gathered or stepped.
The reference is identical because it is a function of files and of static exports that
`PASS_INIT` publishes from the freshly initialised Yelmo, independent of the state;
`geo_glob_restart.nc` records the reference fields too. A cross-resolution restore is not
possible: the manifest `grid` (`GAUSS-L<lmax>`) differs and stops.

### 8.3 Bundle and manifest

```
restart-<t>/
  geo_glob_restart.nc          VILMA2: memory, σ_n, rotation, reference, last host terms, se%time    (new)
  geo_glob_imports.nc  geo_glob_exports.nc     fields on the Gauss grid (the last gathered load, outputs)
  <domain>/geo_<dom>_imports.nc  geo_<dom>_exports.nc        (no backend file: restart = .false.)
```

`restart_files` of `vilma_earth_t` is the default (`<name>_restart.nc`); manifest entry
scope `global`, backend `vilma2`, grid `GAUSS-L170` (needs the global instance grid of F5,
WP4), `git` from `build_info_libs`. `GLOBE_ORDER` stays `forcing`, `geo` (geo last: it needs
the domains' `PASS_INIT`, which `run_pass_init` already did).
Size of the file: the 3D memory is the large term (not measured; G7; the restart schedule is
decided after G7, Q12). Restart files stay lossless (§9.1).

### 8.4 Exactness (the gate)

chunk == cont for every file of the run and every bundle, `kryos-compare --chunks`. What
carries it: stateless external ice and reference; `se%time`, `dt_try` and the ramp state
restored; the gather reads only restored ice exports; the first step after a restore does
not recompute the global layer's past. What can break it: OpenMP reduction order (first
the `cont2 == cont` determinism run, same binary, same threads, as regression §1.2), a
Yelmo restart that is not exact (phase 2 Q9, fixed in S4), VILMA2's fast-math local builds
(gate flags `-O2 -fp-model precise`; bit identity there is by dispatch, `ht` §0).

## 9. Output

| File | Content | Rounded |
|---|---|---|
| `geo_glob.nc` | exports on the Gauss grid at `out_2D`: `dz_bed`, `dz_ss`, `C_ocn`; stream `gia`: the global ice field as VILMA2 got it (`host_w`, `host_Ig`, `host_C`, `H_ice_ext`, total effective ice) | yes |
| `geo_glob_ts.nc` | `bsl`, `bsl_ive`, `bsl_ive_host`, `dphi`, `A_ocn`, `mass_resid`; per domain and external: grounded ice change [Gt, m ESL]; clipped external mass; VILMA2 timers (`t_solver`, `t_remap`, `stepper%n_accept`) | no (time series stay lossless) |
| `geo_glob_vilma.nc` (stream `vilma`, optional) | VILMA2's own `vilma_write_step` file (`rsl`, `z_bed`, …) for comparisons | yes (Q13, decided) |
| `geo_<dom>.nc`, `geo_<dom>_ts.nc` | the port's exports on its grid (the default writer) | 2D yes, ts no |
| `globe_ice.nc` | the hub output of §5.4 (reference, merged ice, coverage) | yes |
| `geo_glob_restart.nc`, bundle `*_imports.nc`, `*_exports.nc` | restart (§8.3) | no (lossless) |

The default 2D writer writes `xc`, `yc` in km (`kryos_outfile.f90:175-181`); a global lon-lat
grid needs `lon`, `lat` in degrees (writer extension, WP5; the Gauss grid has explicit
latitudes, written from the grid object).

### 9.1 Precision (decided, Q12 addition)

Every phase 4 output except restarts and time series is written with ncio bit rounding,
`prec_abs` / `prec_rel` per variable, as in the yelmox ncio precision rollout (`ncio` at
fesm-utils 4020c7f and later; the 6th `prec` column of the `variable_io` tables; `nc_lossless`
set from `[kryos] write_lossless` stays the global debug switch, which `ctrl.write_lossless`
maps to). Restart files, bundle field files, time series, grids and debug output stay lossless.

- **kryos side** (no new mechanism): the default writer of `geo_glob`, the ports and
  `globe_ice.nc` take the `prec` column of `input/kryos-variables-geo.md` (and the new rows
  for the gia stream and `globe_ice`); `outfile_var2d_dp/sp` already pass `prec_abs`,
  `prec_rel` (`kryos_outfile.f90:212-230`). Proposed values: `dz_bed`, `dz_ss`, `z_bed_ref`,
  `H_ice_ext`, `host_Ig`, `H_ice`, `z_bed` 1e-3 m; `dz_bed_dt` 1e-4; `C_ocn`, `host_w`,
  `host_C` 1e-4 (fractions); `host_F` 1e-3 m; `rsl` 1e-3 m.
- **VILMA2 side** (F13, **decided Q13 a**): `vilma_io.f90` does not round. The own diagnostics file
  (`geo_glob_vilma.nc`, stream `vilma`, run through `vilma_write_step` and
  `vilma_write_horizontal`: `rsl`, `z_bed`, `h_ice`, `C_ocean`, `u_east`, …) needs the
  feature in VILMA2: a `prec` column in `input/vilma-variables.md` (read by the existing
  `variable_io`), `prec_abs=v%prec_abs, prec_rel=v%prec_rel` in `put2d`, `put_hor`, `put3d`,
  and a `lossless` argument that the restart path passes `.true.` (the restart shares the
  `put*`; the yelmo rollout does the same with its `lossless` optional). About 40 lines and a
  test (restart round trip unchanged and lossless; rounded variable within tolerance). VILMA2
  gets a `prec` column and `prec_abs`/`prec_rel` in its `put*`, with a lossless restart path
  (Q13); the standalone driver rounds too, one mechanism.
- **Tests** (G4 extension): rounded and lossless runs give identical restart bundles; max
  error / tolerance ≤ 1 for every rounded variable of the files above; file sizes reported
  (G7).

## 10. Phase 3 and phase 4

Phase 4 starts once phase 3's two-domain work is in (Alex, 2026-10-29). Phase 3 (`phase3-plan.md`)
delivers two domains in one run (NH-32KM + ANT), the `bsl` reduction for N domains, shipped forcing,
`smb_simple`, forcing-increment bundles, output precision and the kryos self-baselines `K3-sci`;
the global geometry hub and everything global is phase 4. Phase 2 WP2 + WP3 (gate `ant-full`) and
WP4 (gate `nor`) are done before phase 3 `sci`. This section and `phase3` §3.5 agree.

### 10.1 What phase 4 needs from phase 3

| # | Need | Used by | Phase 3 |
|---|---|---|---|
| P1 | Two domains in one run: NH (`yelmox_North`, NH-32KM, the `nor` recon cold start) and ANT (32 km ice grid, hub ANT-16KM as in the shipped base), disjoint, name order; `run_check_names`, loops, one par file (assembled from the two converted files) | everything; G1, G4-G6 | WP1; `sci-const`, `sci` (N10) |
| P2 | The `bsl` reduction for N domains, order-independent (`reduce_bsl`), fastiso | G1 (fastiso identity of the two-domain runs); for `vilma2` `reduce_bsl` returns | WP1 item 3 (order swap, budget) |
| P3 | Shared tsgen forcing: one `forcing` instance read by both domains (coefficients, `gather_forcing` in name order); the shipped forcing of each base (North `alpha` + derived `dTo`, ANT none) | science case climate | Q16 (`scaled` is Later) |
| P4 | N-domain output and restart: the bundle with the globe at the top and a folder per domain, the manifest | §8.3 adds the global instance's grid | unchanged `main` mechanism, gated by `sci` (chunk) |
| P5 | Forcing-increment bundles (`restart_every_df`) | transients (G5, production) | WP4, `sci-dfm` |
| P6 | Output precision in the default writer: the `prec` column, ncio rounding, restart lossless, `write_lossless` switch; backend streams with tables | §9.1 reuses the mechanism for the new variables | WP3, `sci-nl` |
| P7 | `smb_simple` and the snapshots + tsgen + chion path | the two surface/climate options of the science case | WP2 (Q17), `sci-smb` |
| P8 | The NH recon cold start and the `vm-m16` restart fix | G4 chunk == cont with NH | phase 2 WP4 (gate `nor`) |
| P9 | Kryos self-baselines of the two-domain fastiso run: `K3-sci` (files and bundles of `sci` and `sci-const` at 1000, kryos commit of the phase 3 sign-off); no yelmox baseline for it | G1, G1b, G6 drift reference | §3.1, registered at sign-off |

P1-P9 block the Levante gates, not the harness tests: WP1-WP3 and WP6 run on the Mac with a
mock ice (§12).

### 10.2 What phase 3 does not provide (phase 4 builds it)

Phase 3 no longer has a globe geometry hub or rank-2 globe items; all are phase 4:

| Item | Phase 4 package |
|---|---|
| globe coupler with several global grids (background, Gauss), `coupler_set_area` per ice grid, domains' grids resolvable by name (N1) | WP4 (`globe%cpl`) |
| a global instance has a grid: `ctx%grid_name`, `nx`, `ny`, manifest `grid` (N2) | WP4; `vilma_earth_t` sets its grid at `init` and the globe copies it |
| rank-2 fields on a global geo (N3) | WP1 (declaration by `opt%host`) |
| rank-2 globe ↔ domain rows, gather by merge, scatter by the wiring rows (bilinear, domain coupler), lag check (N4) | WP4 (`gia_gather`, `bsl_to_geo`; `exchange_globe` rank 2) |
| footprint and coverage of each domain on each global grid, owner weights (N5) | WP3 (`globe_ice`) |
| global geometry hub: reference assembly, merged ice, checks, `globe_ice.nc` (N6) | WP3, §5.4; gate G1b |
| generic global step in `run_sequence` (N7) | WP4 (§7.1; the hooks are `gia_*`) |
| lon-lat output (N8) | WP5 |
| `GLOBE_ORDER` (N9) | unchanged (`forcing`, `geo`); no `ocean` |

### 10.3 Phase 3 scope that phase 4 does not use

`ocean`/`obm`, freshwater coupling, `ocean_glob`, `scaled` coefficients, globe ↔ globe rows,
blocks 3b and 5b, the bipolar mode of the converter and any yelmox two-domain reference are Later
in `phase3` (§7). Nothing here assumes them; the N-domain `bsl` gather (`dV_bsl`) is unused by
`vilma2`, the `M_ice` gather (`gather_forcing`, name order) is shared.

## 11. Changes outside the phase

**design.md** (applied in the commit "docs: phase 3/4 edits"):
1. §5 vocabulary: geo exports add `dphi`; `dz_ss` relative to `dphi`; ice exports add `footprint`,
   `H_ice_pd` and the `flotation` capability; geo global imports/exports rank 2 for `vilma2`.
2. §6 "Mechanism" item 1: the gather may carry host-resolved terms (VILMA2) and not only
   `con`-remapped fields; item 3: `z_sl = dphi + dz_ss`; last paragraph: the reference for a
   global component comes from the ice's present day (`dta%pd`), not the hub reference.
3. §14 phase 4: dependencies (`vilma2:main`, `libvilma.a`), the global geometry hub (grids, footprints,
   assembly, checks, output) is phase 4 (§14 item 3 already says phase 3 is two domains), the
   two-pass cold start, FastIsostasy + VILMA2 as a later package, the Gauss grid note (F1), output
   rounding (§9.1).
4. §15 D4: "split onto each hub" → onto the domain's geo grid (ice grid by default), the hub through
   the ice; D5: add the regional port and `dphi`; D12-style note that the background of the globe is
   a datasrc (Q5).
5. §12 row 3: mark done (merge, split, integral, set_area, Gaussian grid).
6. §3 diagram: `globe_ice` under the globe.

**contracts-design.md**: §2.1 ice (`footprint`, `H_ice_pd`, binding), §2.2 geo (`dphi`; the global
vilma2 vocabulary of §2.3 here), §3.3 geo row (two `vilma2` types), §4.2 table (§2.4 here), §6 last
paragraph ("VILMA2 (phase 4) uses the same calls" → gather by merge, `geo_reduce` is the solve,
scatter by split). **core-design.md**: §3.2 (two-pass cold start), §3.3 (reference before the geo
restore), §3.4 (§7.1). **adapters-design.md**: new §12 `geo = vilma2`; yelmo adapter additions.

**VILMA2**: change 4 is merged and pushed (main d8041e8); `libvilma.a` rebuilt against
`phase0-kryos`. No further API change is required for the coupling (change 5, later). One small
change for output rounding (F13, §9.1, Q13): `prec` column, `prec_abs`/`prec_rel` in the `put*`,
restart lossless.
**fesm-utils**: none required. To verify in WP3: `fesm_datasrc` reads a slice of a time-dependent
file as a static field (else use `transient` at `ext_tref`), the time axis "years ago" of the
IceHistory files against the model clock, `coupler_merge` handles a projected source and a
Gaussian target with the set areas. **FesmData**: none required; optional, a Gaussian grid of
VILMA2's size if products on it are wanted (Q9). **yelmo adapter**: §2.3.

## 12. Work packages

Sizes: S ≲ 300 lines, M ≲ 1000, L more. Gate names §13.

| WP | What | Size | Depends | Gate |
|---|---|---|---|---|
| WP0 | Stack S5 = S4 + `vilma=1`: `vilma2:main` (d8041e8, Q10) in `.configme/manifest.toml`, `make shtns` in fesm-utils (serial and omp) on Levante and the Mac, `tools/gate/stack.toml`; frozen maps for the new pairs (§13.1); data check of the files this plan reads | M (Levante time) | — | G0 |
| WP1 | Contracts: geo opts, rank-2 global declaration, `dphi` (fastiso exports `dphi = bsl`; the row aux changes), `footprint`, `H_ice_pd`, `flotation` binding (yelmo adapter + wrapper), rows and `T_GIA`, `wire_resolve` rule for host imports | M | — | G1 (fastiso identity), G2a |
| WP2 | `vilma2` adapter, stub, factory, build (§3), `[geo_vilma2]`, `[globe_ice]` templates, Gauss grid | M | WP0, WP1 | G2a, G2b |
| WP3 | `kryos_globe_ice`: global grids, footprints and coverage, background and external datasrc, reference, host terms, merges, checks, diagnostics, the hub and `globe_ice.nc` (§5.4) | L | WP1 | G1b, G2a |
| WP4 | Run sequence: `gia_*`, two-pass cold start, restore, step; `globe%cpl`, the global instance grid, rank-2 globe rows (all of §10.2) | L | WP2, WP3, phase 3 (two domains) | G2b, G3 |
| WP5 | Restart bundle and manifest for the global grid; output streams, lon-lat writer, output rounding (kryos side, and the VILMA2 change, Q13) | S-M | WP4 | G4 |
| WP6 | Tests: mock ice instance (`ice_adopt`), unit tests, the VILMA2 standalone comparison scripts (Julia) | M | from WP2 | G2 |
| WP7 | Levante gates G1-G7, cost table, baselines registered | M (compute) | WP4, WP5 | G1-G7 |
| WP8 | Science-case par files (`par/`): NH + ANT, `vilma2`, snapshots + tsgen + chion and the `smb_simple` variant; plots | S | WP7 | — |

Critical path: WP0 ∥ WP1 → WP2 ∥ WP3 → WP4 → WP5 ∥ WP6 → WP7 → WP8. Outside kryos: phase 3's
two-domain work (§10.1; WP4 integration and the Levante gates wait for it), phase 2 `nor` and
`ant-full` (NH and ANT domains that restart exactly). VILMA2 is ready (d8041e8; the output
rounding change, Q13, is small and part of WP5). WP1-WP3 and WP6 (harness) need neither phase 3 nor Levante:
they run on the Mac with `vilma=1 yelmo=0` using a mock ice.

## 13. Tests and gates

No yelmox baseline exists. The references are VILMA2 itself, analytic cases, budgets, the phase 3
kryos self-baselines `K3-sci` (fastiso side) and, once G4 passes, `K4-*` (registered as in regression §3.5).

### 13.1 Unit tests (Mac, `make check`)

- `test_globe_ice`: external history (anomaly of the slice is zero at `ext_tref`; clip, hold-last,
  stop-before-first); merge conservation (`coupler_integral` vs source integrals, 1e-9); overlap
  stop; footprint change stop; **reference self-test** (§5.1: gather at the reference state gives
  zero anomalies); reference assembly twice is bit identical; the two accumulators of `I` and `C`
  reproduce `ht`'s weighted fields on a case with partial cells.
- `test_geo_vilma2` (`vilma=1`, mock ice): init at the reference, `update`, `dt = 0` seed;
  restart round trip (write, rebuild, read, continue: bitwise); the port copies exactly;
  `se%remap` stop on a wrong grid; densities set from `cnst`; the stub stops with the rebuild message.
- `test_gia_sequence`: two mock domains (disjoint footprints), cold start order, restore order,
  `reduce_bsl` returns, fastiso instances unaffected.
- Yelmo-linked (`yelmo=1`): `flotation` binding vs `tpo%f_grnd` on the initial geometry of grl, ant,
  nor (bitwise, or the listed differences); `H_ice_pd == dta%pd%H_ice`; `footprint` == `mask_ice /=
  NONE`.
- `test_globe_ice` also covers the hub (§5.4): coverage and owner weights with two fake overlapping
  domains, the merged reference and ice conserve volume (1e-9), `globe_ice.nc` read back equals
  the offline merge, domain-order swap bit-identical.
- Precision (§9.1): `test_outfile` rounding on a Gauss-grid variable (within tolerance, lossless
  switch restores exact), VILMA2's `test_io_prec` (rounded diagnostics, restart exact, Q13).
- New map pairs to freeze (Levante `mapgen`, as WP0 of phase 2): `GAUSS-L170 ↔ {NH-32KM, ANT-32KM}`
  (con for the merge, bilin for the scatter), `GAUSS-L170 ↔ ANT-16KM` and `ANT-16KM ↔ ANT-32KM`
  (the shipped ANT hub; a port on the hub grid, and the hub update), `GLOBAL-0.5DEG →
  {GAUSS-L170, NH-32KM, ANT-32KM}` (con).

### 13.2 Gates (Levante unless stated; flags `-O2 -fp-model precise`, production `-Ofast` once as information)

| Gate | Test | Pass |
|---|---|---|
| G0 | S5 builds with `vilma=1` and `vilma=0`; `make check` passes both | builds, tests (the vilma tests skipped when off) |
| G1 | no regression: K(S5, `vilma=1`, fastiso configs) == K(S4) for grl, ant, ant-full, nor **and `K3-sci`** (the phase 3 `sci` and `sci-const` files and bundles at 1000, P9; no yelmox baseline); the `dphi` change and the split cold-start loop are bit-neutral | bit identical |
| G1b | **hub, observer-only**: the phase 3 `sci` run (two domains, fastiso, S5 binary) with `[globe_ice]` set (§5.4): (1) bit-identical to the run without it; (2) merged reference and ice conserve volume (`coupler_integral` against the domains' ice and the background outside the footprints, 1e-9 relative, the round-off of the conservative maps: observed 2.6e-10 and 1.5e-10; logged at init and every output); (3) `cont2` determinism and domain-order swap bit-identical; (4) `globe_ice.nc` equals the merge computed offline (Julia, `KryosTools`) | (1), (3) bit identical; (2) 1e-9; (4) equal |
| G2a | **VILMA2 standalone vs kryos, same load**, analytic: a grounded land-ice cap (no coast in the host region), mock ice; the standalone driver `vilma.x` on the Gauss-grid load vs `geo_glob` | `dz_bed`, `dz_ss`, `bsl` equal to the remap accuracy (proposed 1e-6 relative; set after the first run) |
| G2b | same, realistic: mock ice replaying a PaleoMIST deglaciation on the NH and ANT grids vs VILMA2 standalone on the same history. Differences are the host-resolved flotation: reported against `lmax` 85, 128, 170; not an equality gate | difference falls with `lmax` (convergence report), `bsl` within a proposed 2%; thresholds (set 2026-10-31, `gate-reports/p4-g2b.md`): `bsl` 1 % of the range, `dz_bed` host rms 2 % of max, rim and outside 1 %, not larger at lmax 170 than at 85, lost-cell uplift median ≥ 0.9 and within 0.05 of the standalone; S11 PASS |
| G2c | **spin-up** (Q7): the 80 ka cold start with `equil_time_max` = 0 and a few values (and `equil_rate_tol`); reports per value the bed rate and `dz_bed` at the seed, the drift of `dz_bed` over the first intervals, the spin-up cost | `equil_time_max` chosen from the report (Alex) |
| G3 | **mass and `bsl` budgets** per step: (i) merge conservation (§4.5), (ii) VILMA2 `mass_resid` ≤ `sle_tol`-level, (iii) `∫ C·rsl dA = −(ρ_i/ρ_w) ∫ ΔI_g dA + ∫ F dA` recomputed by kryos from the exports and the gathered terms (`F` the host part and the own part's fill, which VILMA2 exports as its fill integral), (iv) `bsl_ive_host`, the host part of `bsl_ive` (`−(ρ_i/ρ_w) ∫(I_g − I_g,ref) dA / ∫C dA` over the host columns, which Yelmo owns with its flotation), against the sea-level equivalent of the Yelmo grounded-ice change (the sum over the domains), compared **as volumes** (`bsl_ive_host·A_ocn` of VILMA2's solve against the equivalent times the sea surface area of the constants; as sea levels they differ by the change of the ocean area, 3.4 % at an LGM load). `bsl_ive` as a whole is not compared: its `1 − w` part counts only the grounded external ice, by VILMA2's flotation, while the external equivalent of the log counts all PaleoMIST ice of that part, floating included. The difference `bsl_ive − bsl_ive_host` and the floating external ice (external equivalent minus own part, as volumes) are reported for information, (v) `A_ocn_ref` against 3.62e14 m² | (i) 1e-9, (iii) 1e-6 relative; (iv) 1e-6 as volumes (on `bsl_ive_host`; the former comparison of `bsl_ive` had 1e-5, observed 4.3e-7); (ii), (v) reported with the proposed bounds 1e-6, 1%. Until VILMA2's fill integral lands, (iii) is reported as the match of the recomputed term with VILMA2's `bsl_ive·A_ocn` (observed 8e-8) |
| G4 | **chunk == cont**: harness (mock ice, 2 + 2 steps vs 4, bundles and files) and the science case on Levante (`cont`, `cont2`, `chunk`, `kryos-compare --chunks`, 200-500 yr), including `geo_glob_restart.nc`; precondition `cont2 == cont`; and the rounding check of §9.1: a rounded and a lossless run give identical restarts, max error / tolerance ≤ 1 on every rounded output variable | bit identical, as phase 1 and 2 |
| G5 | **ANT deglaciation**: `yelmox_Antarctica` base with VILMA2 global, NH as external PaleoMIST only (one active Yelmo domain), over a transient window with a large ESL change (21 → 0 ka if the base forcing allows; to be settled in WP0's data check); checks: signs and magnitudes of `dz_bed` (uplift under loss), G3 budgets, `bsl` against the ESL of the ANT ice change; **documented-difference report** against the same run with `geo = fastisostasy` (different Earth, far field, rotation) and a VILMA2 standalone replay of the recorded ANT load | budgets pass; the report is read by Alex |
| G6 | science case smoke: both domains active (NH-32KM + ANT), climate and surface from snapshots + tsgen + chion, and again with `smb_simple`; ≥ 2 kyr of the `sci` configuration with `geo = vilma2`; no stop, G3 each step, restart at the middle; drift report against `K3-sci` (the fastiso `sci` files over their 1 kyr) | stable, budgets, report |
| G7 | cost (§14), restart size and write time, output sizes rounded vs lossless | table; decisions on `dt` (Q6) and the restart schedule (Q12) |

Registered baselines after G4 and G6: `K4-ant`, `K4-sci` bundles and files (`K4-sci-smb` for the `smb_simple` variant).

## 14. Cost of T170 + 3D viscosity

From `vilma2 doc/performance.md` (Levante, `ifx -Ofast`, 64-128 threads, F9): lmax 128, 3D:
6.1 s per call at 100-yr steps, setup 172 s, 10 GB. With the fitted `lmax^2.40` the estimate
for T170 is ≈ 12 s per call (≈ 2.0×), setup ≈ 6 min, memory ≈ 20-25 GB (cubic scaling assumed).
A call at `dt = 10` has fewer sub-steps than one at 100 yr (the explicit Maxwell ceiling, `cfl`, is
of order 100 yr for the mantle layers), so it may cost less, but the SLE solve per call stays. Budget
if the cost per call holds: 120 kyr at 10 yr = 12 000 calls ≈ 40 h of VILMA2; at 50 yr ≈ 8 h; at
100 yr ≈ 4 h. G7 measures, on the science case, per call at `dt` = 10, 20, 50, 100: wall of
`gia_gather` (should be ms), VILMA2 `update`, spin-up passes, setup, peak RSS, restart file size,
the VILMA2 share of the step wall, with `OMP_PROC_BIND=spread` and ≤ 64 threads (VILMA2
recommendation; Yelmo and VILMA2 share the thread count). Output: the table and a
recommendation for the production `dt` and `lmax` (T170 is settled; the table says what it costs).

## 15. Risks

| Risk | Mitigation |
|---|---|
| Flotation outside the footprint differs from Yelmo's (pinning points, `f_ice`) and gives spurious `ΔI_g`, `F` at t0 | the reference self-test and the binding test; one code path for reference and step |
| One-interval lag of the load, and a seed solved with `ice0` | documented; spin-up under the same load; G2b measures against a standalone run |
| Cost: 12 s per call at T170 3D, ≥ 20 GB | G7; `dt` up to 50 yr, ports hold anomalies (linear-in-time interpolation of the anomalies between calls is an option, not planned) |
| Restart exactness with OpenMP and 3D memory | `cont2` first; same thread count; bundle compare |
| Restart file size and write time of the 3D memory | G7 measures; restart schedule |
| Host region cut through an ice sheet or a basin at the edge | external mean over the whole cell for the `1 − w` part (D17) |
| NH recon changes `mask_ice` after `PASS_INIT` | the first footprint check stops |
| PaleoMIST covers 80-0 ka only (decided for now, Q4); no Eurasia in GLAC-1D nn9927 | runs of 80 ka or shorter; a longer history is a later datasrc change |
| Rounded output mistaken for restart data, or the VILMA2 `put*` shared by restart and diagnostics rounding a restart | restart paths pass `lossless`; G4 compares rounded vs lossless restarts (§9.1) |
| Closed depressions below sea level count as ocean (Yelmo `f_grnd`, VILMA2, FastIsostasy alike) | limitation as in the yx-plan; connectivity test later |

## 16. Decisions taken in this plan (reviewable)

D1 Global `geo_glob` + regional ports; the core gathers from the ice instances, the port is a
receiver. D2 `dz_ss` = `N`, scalar `dphi` (Q1). D3 The gather is a dedicated core procedure, with
`T_GIA` rows only for resolution, logging and lag checks. D4 Gather sources are contract exports at
the previous step plus static ones; flotation through the ice binding (Q2). D5 VILMA2's reference is
the background X on the Gauss grid; the spliced reference is the hub's and diagnostic. D6 Source areas
are true spherical (§4.4). D7 The external ice is evaluated at every geo step; no separate cadence.
D8 The background is a datasrc group of the globe, not an instance (Q5). D9 `geo_glob` cadence 10 yr as default (Q6). D10 Two-pass cold start (Q8). D11 The Gauss grid is
`GAUSS-L<lmax>`, built and written by kryos (Q9). D12 The footprint is fixed; a change stops the run.
D13 Densities and `r_earth` come from `phys_const`; the `[vilma]` keys are overridden and logged.
D14 `geo_glob_restart.nc` is written with the label `se%time` and read as the last slice. D15 The
first kryos results are the baselines (no yelmox run exists). D16 The port is on the ice grid (chain
default). D17 The `1 − w` part uses the whole-cell mean of the external field. D18 In host cells outside the footprint the bed, in the reference and in a step, is the domain's (`z_bed_ref`, the ice's `z_bed`) and the ice is external (X in the reference, `H_ext(t)` in a step); the yx-plan does not say. D19 The global geometry hub is `globe_ice_t` (§5.4): one code path for the hub and VILMA2, observer-only without `geo = vilma2`; phase 3 delivers two domains only. D20 Output rounding by ncio for every phase 4 output except restarts and time series (§9.1). D21 `vilma2:main` is pinned; kryos stays general in N domains and backends, the science case is NH + ANT.

## 17. Decisions

Decided by Alex, 2026-10-29 (numbering as in the draft; `Was:` is the draft's question).

1. **Decided: add `dphi` to the geo contract; FastIsostasy exports `dphi = bsl`** (bit-neutral,
   G1; the row `geo.dz_ss → ice.z_sl` adds `dphi`; `dz_ss` is the anomaly relative to it; `bsl`
   stays the barystatic diagnostic) (WP1). Was: `dphi`, or the VILMA2 port exporting `Δφ` as `bsl`.
2. **Decided: an ice-component flotation binding that calls Yelmo's own routines**
   (`has_flotation`, `ice_flotation`) (§4.3, WP1). Was: that, or a kryos function in the core.
3. **Decided: background X = `IceHistory-PaleoMIST` at 0 ka on `GLOBAL-0.5DEG`** (§5.2). Was:
   that, the 0.25° grid, or ICE6G_C/GLAC1D.
4. **Decided: PaleoMIST back to 80 ka for now** (runs of 80 ka or shorter; a longer history
   later). Was: PaleoMIST, GLAC-1D nn9927 for 120 ka, or an LGM start.
5. **Decided: the background ice is a datasrc of the globe** (`[globe_ice]`, no new component
   backend). Was: datasrc, or an `ice` instance of scope global with backend `data`.
6. **Decided: geo cadence 10 yr; revisit after G7** (§7.2, §14). Was: 10 yr, 20-50 yr, or
   10 yr with interpolated anomalies.
7. **Decided: VILMA2's own spin-up, holding the t0 load.** `equil_time_max` stays a parameter,
   set after a test: gate G2c (§13.2), not an open question. Was: VILMA2's spin-up, none, or a
   pre-history run.
8. **Decided: the two-pass order, the seed between cold-start blocks c) and d)** (§8.1). Was: that,
   the yx-plan order, or a seed before the ice `PASS_STATE`.
9. **Decided: kryos builds and names `GAUSS-L170` itself** (§3.5). Was: that, a FesmData grid, or
   `nlat = 256`, `nphi = 512` in VILMA2.
10. **Decided: pin `vilma2:main`** (d8041e8: changes 1-4; host terms merged and pushed; API in the
    README "Host-resolved terms" and `doc/design.md` §15) (§3.3). Was: merge `host-terms` then pin
    main, or pin the branch.
11. **Decided: true spherical cell area for the gather** (§4.4). Was: that, or the plane `dx·dy`.
12. **Decided: restart schedule and 3D-memory size after G7** (G7 reports size and write time;
    the options of the draft stay: the common schedule, or a thinner one for `geo_glob` with
    compression). **Addition: output precision** (§9.1): every phase 4 output (global sea-level
    file, VILMA2 diagnostics on the Gauss grid, per-domain geo port output) uses ncio rounding,
    restarts lossless.

Decided later the same day:

13. **Decided: VILMA2 gets output rounding (a): a `prec` column in its variable table and
    `prec_abs`/`prec_rel` in its `put*`, with a lossless restart path** (F13, §9.1, WP5; ~40
    lines + a test; the `vilma` stream and the standalone driver round too, one mechanism).
    Was: (a) that (recommended), or (b) kryos writes all global output itself (the default
    writer already rounds), the `vilma` stream is dropped and VILMA2's diagnostics file stays a
    standalone-driver feature.
14. **Decided: `globe_ice.nc` is on the background grid** (`GLOBAL-0.5DEG`: spliced reference,
    merged ice, coverage; the grid the data are native on, small) (§5.4). Was: (a) that
    (recommended), (b) the Gauss grid (the grid VILMA2 got, no extra remap), or (c) both.

## 18. Later: FastIsostasy + VILMA2

Not in phase 4 (VILMA2 change 5: response channels, `scheme = "fe"`; FastIsostasy change for an
external sea surface and far-field bed). In the structure above it is: a regional `fastiso_geo_t`
in the domain next to the global `vilma_earth_t`; the global geo exports one extra rank-2 field per
domain (`dz_bed_far`, the response outside the FastIsostasy region) through the same scatter rows;
the regional instance reads `dz_ss` from the global instead of its own geoid and adds the far field
to its bed; the load of FastIsostasy is the ice of the domain it already reads. The init stop of
§3.2 names the missing package.
