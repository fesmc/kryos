# kryos phase 2 — plan

Status: reviewed 2026-10-26; all open questions (§5) decided. No code yet. Phase 1 is done: K == R2 bit for
bit, grl and ant, 1 kyr + chunks + 15 kyr (`gate-log.md`).

Phase 2 (`design.md` §14): "`data` backends for every component; esm and rembo
climate; smbpal and smb_simple; experiments plugins. Target: all yelmox base
configs." smbpal is dropped (Q14: deprecated, chion is its successor).

Refs: `yx` = yelmox `topo-prep` (`~/models/yelmox/.claude/worktrees/topo-prep`,
74f3a2b = v2.4 + R2 commit). `contracts`, `core`, `adapters`, `regression` =
`docs/dev/*-design.md`.

## 1. Scope

### 1.1 Bases × features

Columns: GRL `yelmox_Greenland`, ANT `yelmox_Antarctica` (shipped:
relax = opt = True), GRP `GreenlandPaleo`, LIS, NOR `North`, PYR `Pyrenees`,
SRG, eANT/eGRL `esm_*`, rGRL `rembo_Greenland` (all `yx/yelmox/*.nml`).
`1` = done in phase 1, `x` = phase 2 work, `·` = not used.

| Feature | GRL | ANT | GRP | LIS | NOR | PYR | SRG | eANT | eGRL | rGRL |
|---|---|---|---|---|---|---|---|---|---|---|
| hub / ice grid | GRL-16KM | ANT-16/32KM | GRL-PAL-16KM | LIS-32KM | NH-32KM | PYR-250M | SRG-250M | ANT-8KM | GRL-8KM | GRL-16KM |
| converter refuses today (`tools/test/test_convert.jl:261-274`) | — | `sim.relax`, `sim.opt` | `bias.active` | `surface_method` | `init_method` | `with_shelf` | `with_isostasy` | `climate` | `climate` | `climate` |
| ice = yelmo | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| geo = fastisostasy | 1 | 1 | 1 | 1 | 1 | x (geo 2KM ≠ ice 250M) | · | 1 | 1 | 1 |
| geo = data (static bed) | · | · | · | · | · | · | x | · | · | · |
| bsl method | fastiso 1 | 1 | file x¹ | 1 | 1 | 1 | const x¹ | 1 | 1 | 1 |
| climate = snapshots | 1 | 1 | 1 + bias x | 1 | 1 | 1 (`uniform`, dt_clim 50) | 1 (MAR groups) | · | · | ocean only x |
| climate = esm | · | · | · | · | · | · | · | x | x | · |
| climate = rembo | · | · | · | · | · | · | · | · | · | x |
| climate → surface lapse, grids differ | · | · | · | · | · | x² | · | · | · | · |
| surface = chion | 1 | 1 | 1 | · | 1 | 1 | 1 | 1 | · | · |
| surface = smb_simple | · | · | · | x | · | · | · | · | · | · |
| surface = climate (`has_smb`) | · | · | · | · | · | · | · | · | x | x |
| shelf = marine_shelf (profiles) | 1 | 1 | 1 | 1 | 1 | · | · | · | · | 1 |
| shelf, shelf-base input (`has_ocn_shelf`) | · | · | · | · | · | · | · | x | x | · |
| shelf = data (none) | · | · | · | · | · | x | x | · | · | · |
| sediments method 1 (file) | · | · | · | x¹ | x¹ | · | · | · | x¹ | x¹ |
| ice `Qd` from climate | · | · | · | · | · | · | · | x | x | · |
| cold start | none 1 | equil 1 | none 1 | recon x | recon x | equil 1 | equil 1 | equil 1 | equil 1 | equil 1 |
| plugin relax | · | x | · | · | · | · | · | x | · | · |
| plugin opt cb_ref (L21) | · | x | · | · | · | · | · | x | x | · |
| plugin opt tf_corr (L21-points) | · | x | · | · | · | · | · | x | x | · |
| plugin lim_pd_ice | · | · | · | · | · | · | · | · | · | x |
| clock from a timeline group (`run_step = spinup`), calendar | · | · | · | · | · | · | · | x | x | · |
| tsgen `rel`, transient forcing time | · | · | · | · | · | x (−120 kyr) | · | · | · | · |
| Yelmo calving `vm-m16` (restart not exact, `gate-log.md:306`) | · | · | x | x | x | · | · | · | · | x |

¹ code exists from phase 1 (adapters §3, §7), never gated. ² first use of
`T_LAPSE_*` with `dz ≠ 0` (contracts §4.1).

### 1.2 Ported in phase 2 but used by no base (feature gates, §3.3)

NEGIS cb_ref, glacial SMB scaling, `ice_bnd`, `init_marine_H`,
`init_kill_shelves`, `init_time_thrm > 0`, `recon_ref`
(all off in every base, §1.1 extraction; `contracts` §10, `adapters` §11).
pico stays refused until the E30 fix (`regression` §2.5); no base uses it.

### 1.3 Out of phase 2

| Phase | What |
|---|---|
| 3 | globe with global grids, N domains, `ocean = obm`, `yelmox_bipolar` (`yx/yelmox_bipolar/`), bsl `mixed`, forcing-increment bundles (`core` §7), two domains producing one derived coefficient (`contracts` §7), suffixed multi-domain par files |
| 4 | `geo = vilma2` |
| 5 | `ice = cism`, then `pism`; `shelf.tf_shlf → ice`, `own_bmb` (D6) |
| — | `ocean = data` (no `ocean` component before phase 3); geothermal `f_stdev` (0 in every base) |

## 2. Work packages

Each: what, decisions (recommendation in **bold**), depends on, gate.
Sizes: S ≲ 300 lines, M ≲ 1000, L more (ported yelmox lines in brackets).

### WP0 Baselines and infrastructure — M (Levante time)

- yelmox `topo-prep` commit 2: the four `[topo]` keys of R2 in the other 8
  bases. Today `prepared=.true.` is unconditional (`yx/libs/kryos.f90:315`)
  but only grl/ant set the keys, so on `topo-prep` the other bases skip all
  reference cleaning. **Decided (Q1): per-base values from each base's own
  Yelmo keys** (`grad_lim_zb` from ytopo: PYR 1.0, others 0.5; `z_bed_f_sd`
  from yelmo_init_topo: 0 everywhere; `remove_lakes = True`, `H_ice_min = 1`
  as grl/ant). No base smooths (`smooth_* < 1`), so the hub needs no
  smoothing.
- Stack S4 = S3 + `rembo=1` (rembo1 v1.3 is in the stack but not linked,
  `gate-log.md:37`).
- Maps: one `mapgen` R0 run per new base on `v24`, then freeze (as phase 1,
  `gate-log.md:103-118`). New pairs: PYR-2KM↔PYR-250M, SRG-16KM↔SRG-250M,
  the GLOBAL→grid maps of each climate.
- Data check before WP gates (§4). **Decided (Q17):** an agent resolves every
  file each base opens and checks it against the gate's links (`fesmdata`,
  `ice_data`, `isostasy_data`, `maps`, `input`); per-base table. Grid/topo
  products are expected in fesmdata, esm data in `ice_data`; missing files
  found in `/work/ba1442/fesmdata` are linked, never copied, and logged in
  `gate-log.md`; the shared store is never modified; anything else is
  reported to Alex.
- Converter (`tools/data/yelmox2kryos.toml`): lift each refusal with its WP;
  `dt_clim ≤ dtt` → `dt = 0`; `ctrl.run_step` → `[kryos] clock`;
  `--check-par` (`regression` §2.7, still missing).
- Gate: R1(S4) == R0 for every base; R2 vs R0 documented per base.

### WP1 Experiment plugin mechanism (core) — M

What: `src/experiments/`, hooks in `kryos_run_sequence` (today "hooks: none",
`src/core/kryos_run_sequence.f90:233`).

Decisions:

1. Form: **Decided (Q2): abstract `experiment_t` (deferred `init`,
   `hook(dom, id, st)`, optional restart), factory by `kind`, an adopt for
   host-provided experiments** — the D1 pattern. **Held in `run_t`, one
   `experiment_slot_t(:)` list per domain** (not in `domain_t`: that would be
   a circular module dependency with `hook(dom, ...)`); the run sequence
   calls `hook(run%domains(d), id, st)`.
2. Layer: **experiments sit between domain and run, take `domain_t`, reach
   backends only by `select type` and stop on a mismatch** (D1, `design.md`
   §4). Recon needs the domain sequence (hub, climate, surface, below), so
   ice-only plugins are not enough.
3. Hooks (yelmox order, `yx/yelmox/yelmox.f90:180-199`,
   `yx/libs/kryos_startup.f90:131-195`):

   | Hook | Where | Users |
   |---|---|---|
   | `HOOK_STEP_BEGIN` | after forcing, before geo | relax, opt cb_ref, tf_corr |
   | `HOOK_ICE_BND` | after `exchange_into(ice)` and a new ice binding `prepare` (imports → backend boundary, already converted), before `ice_update` | lim_pd_ice, glacial_smb, negis, ice_bnd |
   | `HOOK_COLD_PRE_STATE` | before ice `PASS_STATE` | opt `cf_init`, negis |
   | `HOOK_COLD_POST_STATE` | after `yelmo_init_state`, before the equilibration | kill_shelves, marine_H, recon |

   **Decided (Q3): a `prepare(st)` binding on `component_t` (empty default;
   imports → backend boundary, converted), called by the wrapper before the
   hook and `update`; a general core pass `PASS_SPINUP` after `PASS_STATE`
   (cold-start equilibration; yelmo: equil, thrm, then the exports),
   `HOOK_COLD_POST_STATE` between them. Bit-neutral on grl/ant (same calls,
   same order).** Needed for bit
   identity: lim_pd_ice subtracts 4 m ie/yr after the `sp` conversion
   (`yx/libs/kryos_coupling.f90:334-337`); doing it on the kg m-2 yr-1 import
   is not the same `sp` arithmetic.
4. Config: **`[<domain>] experiments = "relax" "opt"`, one group per
   experiment named by it, with `kind`**; converter maps `[sim]` switches and
   `[relax]`, `[opt]`, `[negis]`, ... groups.
5. Time windows: yelmox uses `time_elapsed`, restarted per chunk
   (`kryos_coupling.f90:45,75,95`). **Decided (Q4): both modes.** A
   continuation restart keeps the run's origin, so windows run on elapsed
   time continuous across chunks (`core` §7: K chunk == K cont, Y chunk ≠ Y
   cont after the first chunk, documented §3). A branch restart (WP13) only
   gives the starting state of a new experiment, possibly at another time:
   the origin, and so the windows, start anew.
6. State: **plugins stateless by default** (cb_ref in Yelmo's restart,
   tf_corr in marine_shelf's, `adapters` §6.4; relax is a function of time).

Depends: —. Gate: unit tests (hook order, `select type` stop, restart round
trip); bit-neutral on grl/ant (no plugin active).

### WP2 relax, opt cb_ref — S (glue ~150; `yelmo/libs/ice_optimization.f90` as library)

- relax: `relax_update` on `ytopo` par at `HOOK_STEP_BEGIN`.
- opt cb_ref: `optimize_cb_ref` on Yelmo `dyn%now%cb_ref` with `dta%pd`
  (`kryos_coupling.f90:75-84`); `cf_init` at `HOOK_COLD_PRE_STATE`
  (`kryos_startup.f90:197-216`).
- Requires `yelmo_ice_t`. Depends WP1. Gate: ANT shipped (§3.1).

### WP3 opt tf_corr — S

- **Decided (Q5): a general calibration, contract fields only, no `select
  type`** (`kryos_coupling.f90:95-115`):
  1. The shelf component owns the thermal-forcing corrections: `tf_corr`
     (optimised; optional start from file, `tf_correction`/`tf_path`) and
     `tf_corr_basin` (per-region constants) move from marine_shelf to
     `kryos_shelf_ocean`, written to the instance's restart, handed to the
     backend with the shelf-base ocean inputs. A shelf backend adds them to
     its thermal forcing; marine_shelf keeps `(T_shlf - T_fp) + tf_corr +
     tf_corr_basin` (bit-identical).
  2. Optional ice exports `H_grnd`, `H_ice_pd`, `H_grnd_pd` (yelmo: copies of
     `tpo%now%H_grnd`, `dta%pd%H_ice`, `dta%pd%H_grnd`; not `H_ice_ref`).
  3. Plugin `opt_tf`: reads `H_ice`, `H_grnd`, `dHidt`, `H_ice_pd`,
     `H_grnd_pd`, remaps ice → shelf through the domain coupler with yelmox's
     method (narrowed to `sp`), updates the component's `tf_corr`.
  4. Check: the `dHidt` export (MEAN over the ice call) equals yelmox's
     instantaneous `tpo%now%dHidt` when the ice step equals the coupling
     step; else an instantaneous export. Gate ANT confirms.
- Depends WP1. Gate: ANT shipped.

### WP4 Cold-start recipes — M (~150: `kryos_startup.f90:252-350`, `kryos_regions.f90:48-59`)

- recon (LIS, NOR), recon_ref, marine_H, kill_shelves at
  `HOOK_COLD_POST_STATE`; `none`/`equil`/`init_time_thrm` stay yelmo adapter
  keys, run in `PASS_SPINUP`.
- recon re-runs hub → climate (forced) → surface → ice mid-recipe
  (`kryos_startup.f90:332-337`): **the core exposes one
  `domain_refresh_boundary(dom, t0)` the plugin calls**.
- recon changes Yelmo `H_ice_ref` (`:327-330`): the ice `H_ice_ref` export
  is re-published after the recipe (consumers: smb_simple, WP7).
- Depends WP1. Gate: NOR (chion), LIS (with WP7).

### WP5 Ice-boundary plugins — S each (~90 + `ice_bnd.f90` 158)

| Plugin | Base | Needs |
|---|---|---|
| lim_pd_ice | rGRL | Yelmo `dta%pd%H_ice` (`kryos_coupling.f90:334-337`) |
| glacial_smb | none | climate exports `ta_ann`, `ta_ann_ref` (bilin clim → ice), `kryos_regions.f90:107-137` |
| negis | none | Yelmo `till_method = -1`, regions `NEGIS:*` (`kryos_coupling.f90:381`) |
| ice_bnd | none | datasrc of an ice history, optional forcing coefficient (`scaled`) |

- **Decided (Q13): `ice_bnd` is a yelmo plugin with its own datasrc**
  (`fesm_datasrc`, `contracts` §10), setting Yelmo `H_ice_ref` in the strip
  cells at `HOOK_ICE_BND`. A data ice-history instance + yelmo `H_ice_ref`
  import (`adapters` §11) only if another consumer appears.
- Depends WP1. Gate: rGRL (lim_pd_ice); feature gates (§3.3).

### WP6 `data` backends — M

What a `data` backend is (`design.md` §8): `component_t + datasrc_t` in the
component module, one field group per export, modes `const | static |
climatology | transient` (as geothermal/sediments now).

| Component | Exports served | Base use |
|---|---|---|
| shelf | `bmb_shlf`, `T_shlf`, `tf_shlf` | PYR, SRG (`const`, = `with_shelf = False`) |
| geo | `dz_bed`, `dz_ss`, `bsl`, `z_bed_ref` | SRG (`const` 0, = `with_isostasy = False`); history mode (D4: difference against the dataset's own reference slice) |
| surface | `smb`, `T_srf` (+ optionals) | none (SMB-forced runs) |
| climate | `tas`, `pr`, `z_srf_clim`, `lapse`, `beta_p`; optional ocean, `smb`/`T_srf` (`has_smb`) | none |
| ice | `H_ice` from a file; geometry from geo imports by flotation | none (geo-only, offline climate/SMB) |

Decisions:

1. Absent component: **Decided (Q6): a `data` instance with `const`** (wiring
   stays complete, D7 gaps still stop); the converter writes it from
   `with_shelf = False` / `with_isostasy = False`.
2. MEAN exports from a file: **the datasrc value at `t1`, put over `[t0, t1]`**
   (documented as sampled, not averaged).
3. ice = data: **Decided (Q7): thickness only; `z_bed`, `z_sl` from geo like a model;
   `z_srf`, `f_grnd` by flotation with the run's `cnst`; `smb`, `bmb`, `cmb`,
   `mb_net` made optional exports in the ice contract (today R, `contracts`
   §2.1)**; `M_ice` from `H_ice`.
4. SRG identity: check Yelmo's initial `bnd%z_sl` (expected 0) equals the
   data `bsl` before the gate.

Depends: —. Gate: PYR, SRG (shelf, geo); replay gate for ice, climate,
surface (§3.3).

### WP7 smb_simple + climate `tsl_ann` — M (`smb_simple.f90` 985)

- Backend `surface = smb_simple` (private module, stateless, no restart).
- Climate vocabulary adds `tsl_ann` (snapshots `has_tsl_ann`,
  `src/backends/climate/kryos_climate_snapshots.f90:132`); with WP5, `ta_ann`,
  `ta_ann_ref`.
- Mask: yelmox reads Yelmo `bnd%H_ice_ref` (bilin ice → surface) every call
  (`kryos_coupling.f90:481-483`). **Decided (Q8): import `H_ice_ref` from the
  ice export** (bilin ice → surface; follows recon and ice_bnd,
  bit-identical); the yelmo adapter puts `H_ice_ref` every step (today once,
  as static). LIS sets `mask_file`, so either
  passes LIS.
- Depends WP4 for LIS. Gate: LIS.

### WP8 surface = climate — S

- Pass-through backend on `grid_surface`: imports climate `smb`, `T_srf`
  (bilin, `kryos_coupling.f90:302-305`), exports them; surface → ice `con` as
  today. Capability check `has_smb` (`contracts` §4.4 item 4).
- Depends WP9 or WP10. Gate: eGRL, rGRL.

### WP9 climate = esm — L (`esm_forcing.f90` 892 + `yelmox_climate.f90:721-857`)

Decisions:

1. Shelf-base ocean (D6 removes the `marshelf_class` reach-in,
   `yelmox_climate.f90:842-849`, `esm_forcing.f90:572,785`). **Decided
   (Q10): the shelf component owns it.** A first, bit-neutral refactor moves
   the "ocean at the shelf base" step out of the marine_shelf backend into
   `src/components/shelf/kryos_shelf_ocean.f90`: shelf-base depth
   (`interp_depth`), vertical weights (`calc_shelf_variable_mean/layer/depth`,
   `interp_method`), cavity extrapolation (`ocn_variable_extrapolation`,
   `extrap_shlf`), code unchanged; keys read by the component; run by the
   shelf wrapper before the backend's update, so every backend gets
   shelf-base values (skipped when the climate exports them,
   `has_ocn_shelf`). `find_ocean` (melt mask) stays in marine_shelf. Gate:
   grl/ant bit-identical. Then esm exports reference, anomaly and variability
   profiles (`depth`); the component interpolates each and sums in yelmox's
   order (`T_shlf = interp(to_ref) + interp(dto) + interp(dto_var)`), a third
   shelf input set. Bit-identical when climate grid = shelf grid (both esm
   bases). A later move of the generic routines to fesm-utils stays open.
2. `Qd` (m3 s-1 per cell, `contracts` §2.1): a per-cell amount under `con`
   remap is wrong across grids. **Decided (Q11): a per-cell amount class in
   fesm_fields** (beside intensive, extensive). Identity map: passed through
   (bit-identical, both esm bases). Other maps: the wiring divides by the
   source cell area, remaps `con`, multiplies by the target cell area. `Qd`
   keeps m3 s-1 per cell in both components. Needs a fesm-utils change and
   cell areas of both grids in the wiring (same new stack as the Yelmo
   vm-m16 fix). (A per-area flux with adapter conversions is not
   bit-identical: x/a*a can differ from x by an ulp.)
3. Timeline: **`[kryos] clock = "spinup"`; `time_ref/hist/proj/esm_ref`,
   `clim_var`, `clim_seed` move to the esm instance group; `clock = "model"`**
   (calendar years, `core` §4 exception; `yelmox_climate.f90:757-767`).
4. Branch restore (spinup → transient at another time) stops today
   (`src/core/kryos_run_sequence.f90:142-147`): **implement in phase 2**
   (field times moved to `time_init`), not needed by the base gates.
5. `random_seed` (global state, `:767`): per-instance seed; `clim_var` is off
   in both bases.

Depends: WP8 (eGRL), WP1–3 (opt/relax in both). Gate: eANT, eGRL.

### WP10 climate = rembo — M (wrapper 95 + branch ~80 + snapshots ocean reuse)

- Constraints only (decision 12, `contracts` §3.3): `build_flag = "rembo=1"`,
  `max_instances = 1`, `fixed_grid = "GRL-16KM"`. Ocean from the snapshots
  private modules (ocean-only reference, `yelmox_climate.f90:183-200`);
  `dT_sum` from the uniform terms. rembo's internal cadences (`dtime_emb`,
  `dtime_smb`) stay inside (special case).
- Restart: rembo1 restores only its precipitation correction (`anf_dat = 3`,
  `rembo1/src/climate.f90:703-757`); yelmox re-steps it on restart
  (`kryos_startup.f90:75,101`), which kryos forbids (`adapters` §1.4).
  **Decided (Q12): make rembo1 restart-exact**: its prognostic state
  (day-level snowpack, albedo and previous-day surface state, `day(...)` in
  `rembo_main.f90`, with `ppcorr`) in a rembo1 restart write/read; a rembo1
  change and a new stack, like the Yelmo vm-m16 fix.
- Depends WP0 (S4), WP8, WP5 (lim_pd_ice). Gate: rGRL (cont and chunk).

### WP11 smbpal — dropped

- **Decided (Q14): not ported.** smbpal is deprecated; chion is its
  successor (ITM, PDD). Removed from `design.md`; the converter refuses
  `surface_method = "smbpal"` as deprecated.

### WP12 Climate transforms bias, taper — S

- Already ported (`src/backends/climate/snapshots/kryos_snap_mods.f90`);
  converter maps `[bias]`, `[taper]` to the climate `bias`, `taper` keys and
  lifts the refusal. Gate: GRP.

### WP13 Core items left open by phase 1 — S each

- Branch restore (WP9 item 4).
- Clock year (`tropical`) vs `cnst%sec_year`: unify (`core` §9, "a phase 2
  item"). Answer-changing for chion/snapshots `pr` (they use `clk%sec_year`);
  **do it with the WP14 rebaseline**.
- Lapse transform with `dz ≠ 0` (PYR) and geo on a coarser grid than the ice
  (PYR: `auto` both ways): first uses, no new code expected.
- MEAN imports for slower consumers (`contracts` §1): no base needs them
  (every MEAN row has equal cadences); **leave for phase 3**.

### WP14 `dp` hub rebaseline — S code, M gate work

- `prec = dp` for the hub (`contracts` §9). Expected bit-neutral where the hub
  is the ice grid (mirror: every base except ANT, hub 16 km / ice 32 km; to
  verify, the reference preparation also runs in the hub's precision);
  changes ANT (anomaly path) at round-off. Registered rebaseline
  (`regression` §3.5).
- **Decided (Q16): last package of phase 2, after every yelmox gate passes;
  then all gates switch to kryos self-baselines at once**, with the
  clock-year unification (WP13).

### WP15 CMIP output — deferred (after phase 2)

- `cmip_output.f90` needs Yelmo internals (`taub`, `ux_b`, `uz_*`, `T_ice`,
  `fmb`, `z_base`) and the shelf's `tf_shlf`, not only ice exports.
- **Decided (Q15): Julia post-processing (`KryosTools`) over the Yelmo
  backend's file + `tf_shlf`, written at `dt_cmip`** (design §9). Not built
  in phase 2 (proof of concept first); no gate.

## 3. Gates

### 3.1 Configurations and baselines

Per base: the phase 1 `-p` (`regression` §1.1: lossless, `tm_2D.dt=100
tm_2Dsm.dt=10 tm_1D.dt=10 tm_rst.dt=500`, 1 kyr from `time_init`), on stack
S4, gate flags; runs `cont`, `cont2`, `chunk`.

| Gate | Base, `-p` extras | Baseline | Packages |
|---|---|---|---|
| `ant-full` (the named phase 2 gate, `regression` §1.1) | shipped ANT | R2 (`yx` 74f3a2b) | WP1–3 |
| `grp` | GRP | R2' (WP0 commit) | WP12 |
| `pyr` | PYR, `time_end=-119e3` | R2' | WP6 |
| `srg` | SRG, `time_end=200`, `tm_rst.dt=100` (dtt 5) | R2' | WP6 |
| `nor` | NOR | R2' | WP4 |
| `lis` | LIS | R2' | WP4, WP7 |
| `rgrl` | rGRL | R2' (S4, `rembo=1`) | WP5, WP8, WP10 |
| `eant` | eANT (`[spinup]` clock) | R2' | WP1–3, WP9 |
| `egrl` | eGRL | R2' | WP2, WP3, WP8, WP9 |
| `grl`, `ant` | phase 1 gates, kept | R2 | regression of everything |

Baseline rule: R0 (v2.4, `v24`) == R1 (v2.4, S4) per base first; R2' =
`topo-prep` + WP0 commit; R2' vs R0 documented as for grl/ant
(`gate-log.md:175-220`). Gate 1: K == R2'.

### 3.2 Expected identity and accepted differences

- Expected: E/C classes as in phase 1, every output file and bundle.
- Documented, not failures:

| Where | Why | Handling |
|---|---|---|
| Y chunk vs Y cont (ANT shipped, eANT, eGRL) | yelmox restarts `time_elapsed` per chunk (relax/opt windows) | compare K chunk vs K cont only |
| chunk vs cont, GRP, LIS, NOR, rGRL | Yelmo `vm-m16` calving restart not exact (`gate-log.md:306-308`): `eps_eig_1/2` not in Yelmo's restart | Decided (Q9): fix in Yelmo before these gates (new stack); then expected identical |
| `bmb_shlf` ties on 16 ↔ 32 km | `adapters` §1.2 (did not occur in phase 1) | registered if they occur |
| esm shelf-base values | interpolation moves to the shelf (WP9) | exact only on identical grids (both bases); else registered |

### 3.3 Gates without a yelmox base

- Feature gates: a base + `-p` switch that yelmox runs: grl
  `sim.use_negis=True`, grl `sim.scale_glacial_smb=True`, grp
  `ice_bnd.active=True`, grl `sim.init_marine_H=True`, `init_kill_shelves`,
  `init_time_thrm`.
- Replay gate (data ice, climate, surface; no yelmox counterpart): a K run's
  lossless exports at every step, read back by a `data` instance; the
  downstream components must equal the original run bit for bit.

### 3.4 After WP14

`dp` hub registered as rebaseline; yelmox comparisons become drift reports;
all gates re-baselined to kryos at the WP14 commit (`regression` §3.5).

## 4. Order (critical path first)

| # | Package | Size | Unblocks |
|---|---|---|---|
| 1 | WP0 baselines, maps, data check (parallel with 2) | M | every gate |
| 2 | WP1 plugin mechanism | M | 3, 4, 7, 8 |
| 3 | WP2 + WP3 relax, opt | S + S | **`ant-full`** |
| 4 | WP12 bias | S | `grp` |
| 5 | WP6 shelf/geo `const` | S | `pyr`, `srg` |
| 6 | WP4 recon + WP7 smb_simple | M + M | `nor`, `lis` |
| 7 | WP8 + WP10 + WP5 lim_pd_ice | S + M + S | `rgrl` |
| 8 | WP9 esm (+ WP13 branch restore) | L | `eant`, `egrl` |
| 9 | WP5 rest, WP6 rest (ice/climate/surface data) | M each | feature and replay gates |
| 10 | WP14 `dp` hub + clock-year unification, rebaseline | S/M | phase 2 sign-off |

Sign-off: one `docs/gate.md` row per gate of §3.1 green, plus the §3.3 gates.

## 5. Open questions for Alex

1. **Decided: per-base values from each base's Yelmo keys** (WP0). Was: R2 for the other 8 bases: add the grl/ant `[topo]` keys (`z_bed_f_sd = 0`,
   `remove_lakes`, `H_ice_min = 1`, `grad_lim_zb = 0.5`) to every base in a
   second `topo-prep` commit (**proposed**), or per-base values?
2. **Decided: abstract `experiment_t`, slots in `run_t` per domain** (WP1). Was: Plugin form: abstract `experiment_t` + factory + slots, taking `domain_t`
   (**proposed**), or plain per-experiment fields?
3. **Decided: `prepare` on `component_t`, general `PASS_SPINUP`** (WP1). Was: Split ice `update` into `prepare` + `update`, and yelmo `PASS_STATE` into
   `PASS_STATE` + `PASS_SPINUP`, for `HOOK_ICE_BND` and the recipes
   (**proposed**)?
4. **Decided: both** (continuous for continuation restarts, new origin for
   branch restarts; WP1 item 5). Was: plugin time windows on elapsed time continuous across chunks (**proposed**,
   differs from yelmox chunks)?
5. **Decided: general calibration; corrections owned by the shelf component, optional ice exports `H_grnd`, `H_ice_pd`, `H_grnd_pd`** (WP3). Was: tf_corr reads Yelmo `H_grnd` and `pd` data by `select type` (**proposed**,
   bit-identical), or only ice exports (not identical)?
6. **Decided: `data const` instances** (WP6). Was: Absent shelf/geo as `data const` (**proposed**) or optional components?
7. **Decided: thickness only + flotation; mass fluxes optional** (WP6). Was: `ice = data`: thickness only + flotation, mass fluxes optional in the ice
   contract (**proposed**)?
8. **Decided: the ice `H_ice_ref` export, put every step** (WP7). Was: smb_simple mask from the ice `H_ice_ref` export (**proposed**) or
   `refs(surface)%H_ice`?
9. **Decided: fix Yelmo first** (cause: `eps_eig_1/2` missing from the
   restart). Was: `vm-m16` bases (GRP, LIS, NOR, rGRL): skip chunk vs cont there until Yelmo
   is restart-exact (**proposed**), or fix Yelmo first (a Yelmo change and a
   new stack)?
10. **Decided: the shelf component (`kryos_shelf_ocean`)** interpolates
    (WP9 item 1). Was: esm ocean: shelf interpolates esm's reference and anomaly profiles
    separately (**proposed**), or esm keeps the interpolation with the
    shelf's parameters?
11. **Decided: per-cell amount class in fesm_fields** (WP9 item 2). Was: `Qd` as a per-area flux in the exchange (**proposed**) or per cell?
12. **Decided: make rembo1 restart-exact** (WP10). Was: rembo: no restart, restore stops (**proposed**), or non-exact restarts
    registered?
13. **Decided: yelmo plugin with its own datasrc** (WP5). Was: `ice_bnd`: yelmo plugin with its own datasrc (**proposed**) or a data
    ice-history instance + yelmo `H_ice_ref` import?
14. **Decided: dropped** (deprecated; chion is its successor). Was: smbpal: port in phase 2 with a feature gate (**proposed**) or defer?
15. **Decided: Julia post-processing, deferred after phase 2.** Was: CMIP: Julia post-processing (**proposed**) or a yelmo adapter stream?
16. **Decided: last package.** Was: `dp` hub (+ clock-year unification): last package of phase 2
    (**proposed**), or now?
17. **Decided: agent check, links only** (WP0). Was: Data on Levante for the 8 new bases (GRL-PAL, LIS, NH, PYR, SRG, ANT-8KM,
    GRL-8KM, ISMIP7, ESM, rembo inputs): who checks, and may missing products
    be synced from `/work/ba1442/fesmdata`?
