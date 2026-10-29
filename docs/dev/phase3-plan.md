# kryos phase 3 — plan

Status: reviewed 2026-10-29, scope narrowed the same day (Alex): two domains, no
yelmox reference for them. All questions (§5) are decided, Q16-Q18 on 2026-10-29
(Q19: order with phase 2). No code.
Phase 3 (`design.md` §14) was "Globe: `forcing`, N domains, `ocean = obm`.
Target: `yelmox_bipolar`". It is now **two domains in one run** for the science
case: NH (`yelmox_North`, NH-32KM) and ANT (as the shipped Antarctica: hub
ANT-16KM, ice ANT-32KM), climate and surface from snapshots + tsgen + chion or
smb_simple, geo = FastIsostasy per domain. No `obm`, no bipolar setup, no globe
geometry hub: those are "Later" (§7, findings kept) and phase 4 (`geo = vilma2`,
`docs/dev/phase4-plan.md`). There is no yelmox run, par file or commit for the
two-domain configuration: the gates are kryos self-consistency, anchored on the
single-domain gates already proven against yelmox.

Summary:

- The core already loops over N domains and has never run N = 2 (row tests only,
  `tests/test_globe.f90`). Phase 3 runs it: a par-file assembly from the two
  converted single-domain files (WP1), the `bsl` reduction over two domains, a first
  two-domain run.
- Gates (§3): `sci-const` (bsl constant, domains independent) must equal two
  single-domain kryos runs bit for bit per domain, outputs and bundles included; the
  single-domain runs are gated against yelmox (`ant`, `nor`), so N = 2 is tied to
  yelmox transitively. `sci` (bsl coupled): `cont2 == cont`, `chunk == cont`, domain
  order swap identical, `bsl` budget closed. Single runs: `sci-smb`, `sci-dfm`,
  `sci-nl`, `sci-plug` (after phase 2 `ant-full`). About 17 runs.
- Forcing: shipped forcing of each base (North `alpha` const + derived `dTo`, ANT
  none), so `scaled` coefficients are not needed and go to Later (Q16, decided).
- smb_simple is back (surface alternative): a work package here (WP2), carried over from
  phase 2 WP7 (Q17).
- New: per-variable output precision (WP3). The default writer already rounds from
  the variable tables and `[kryos] write_lossless` overrides it; the package closes
  the gaps (backend writers with hard-coded precisions) and adds a non-lossless gate.
- Also: forcing-increment bundles (WP4), gate tooling for N domains (WP5).
- What phase 4 gets: two domains in one run, the `bsl` reduction, and registered
  kryos self-baselines of the two-domain run (`K3-sci`). Everything global (grids,
  coupler, rank-2 globe rows, hub, footprints) is phase 4's (§3.5, checked against
  `phase4-plan.md` §10).

Refs: `design`, `core`, `contracts`, `adapters`, `regression`, `p2` = `docs/design.md`
and `docs/dev/{core,contracts,adapters,regression}-design.md`,
`phase2-plan.md`; `p4` = `phase4-plan.md` (§10 there and §3.5 here agree). Code
references are kryos `main` f57e826 (code of 93aab84; includes phase 2 WP1, the experiment
mechanism). fesm-utils = `phase0-kryos` 5662899.

## 1. Scope

### 1.1 Science configuration, feature by feature

`sci` = NH + ANT in one kryos run, obm off, relax/opt off (Q5), fastiso `bsl`.

| # | Feature | yelmox (single domain) | kryos today | Gap, package |
|---|---|---|---|---|
| 1 | N domains in one run | one domain per run; the bipolar driver (two domains, suffixed groups) is not used as a reference | `run_t%dom(:)`, `[kryos] domains` (`kryos_run.f90:76,141`), unique instance names (`run_check_names`, `:303`), loops over domains in `run_cold_start` and `run_sequence`; converted par files are single-domain with unsuffixed physics and snapshot groups, which collide in one file | never run N = 2; par-file assembly: WP1 |
| 2 | Per-domain physics | North and Antarctica differ in `ytopo`, `ycalv`, `ydyn`, `ytill`, `yhyd` (`yelmox_North.nml` vs `yelmox_Antarctica.nml`) | each ice instance names its own physics groups | assembly renames and repoints: WP1 |
| 3 | Forcing | North `alpha` (const 1), `dTo` derived (`from = air_anomaly`), terms `lgm`, `ocean`; ANT no coefficients | `forcing` = tsgen: tsgen and derived; one producer per derived coefficient (a second stops, `contracts` §7) | none (test of the stop with N = 2: WP1); `scaled`: Later |
| 4 | `bsl` from two domains | (single domain only) | gather `dV_bsl` (Σ in name order, `gather_into`), `bsl -= ΣdV/A` (`kryos_geo_fastisostasy.f90:456-467`, `reduce_bsl`, `kryos_run_sequence.f90:338`), every domain gets the full reduced `bsl`; design decision 5: order-independent | gated by order swap and budget (§3.2) |
| 5 | Surface alternative | `surface_method = "smb_simple"`; mask from Yelmo `bnd%H_ice_ref` | not ported (chion only) | WP2 |
| 6 | Cold start | North `init_method = "recon"`, ANT `equil` | phase 2 WP4 (recon, gate `nor`) | phase 2; `sci` waits for it |
| 7 | Output and restart | `restart_every_df` bundles | flat `<instance>*.nc`, bundle with a folder per domain, manifest (`core` §6, §7); increments deferred to phase 3 (`core` §7) | WP4 |
| 8 | Output precision | ncio `prec_abs`/`prec_rel`, `nc_set_lossless` from `ctrl.write_lossless` | default writer rounds from `input/kryos-variables-*.md` (`prec` column); `[kryos] write_lossless`; backend writers partly hard-coded | WP3 |

### 1.2 The step

Unchanged from `main`: without an ocean the domains are independent. `run_sequence`: forcing,
hooks, geo per domain, `reduce_bsl`, then ice, hub, climate, surface, shelf per domain. The only
coupling between the domains is the shared `bsl` and the `M_ice` gathered for the forcing (a sum of two
terms: commutative). With `bsl` constant nothing couples them at all, which is what `sci-const` uses.
No new block, no sequence change (blocks 3b and 5b of the first draft are Later).

### 1.3 The science-configuration par file

Built from the kryos conversions of the two shipped bases (`par/kryos_*.nml`, converter
`tools/bin/yelmox2kryos.jl`; NH from the phase 2 `nor` conversion), not from a yelmox
bipolar file. Findings from reading the bases:

- north = `yelmox_North.nml`: NH-32KM for hub and every component grid, `sim` with
  `init_method = "recon"` (ICE-6G_C at -26 kyr, `recon_regions = region:North_America`),
  `[forcing]` `alpha`, `dTo`, terms `lgm`, `ocean`;
- south = `yelmox_Antarctica.nml`: hub and shelf ANT-16KM, ice ANT-32KM, no forcing coefficients,
  `relax`/`opt` off (Q5);
- physics groups differ per base (§1.1 row 2): per-domain groups, not shared;
- the clock: `[clock]` is run-wide: `time_init = 0`; North's `tstep_const = -21e3` and the
  `rel` method give the same result with `alpha` constant, so one clock serves both;
- instance names `ice_nor`/`ice_ant`... are the single-domain names, so single and two-domain files compare
  by name.

The `bsl` method is the `[barysealevel]` group: `fastiso` (the shipped value) for `sci`,
`const` for `sci-const` and for its two single-domain references.

### 1.4 Not in phase 3

Moved to Later (§7): `ocean` component and obm, freshwater coupling, box → ice and the
climate `refresh`, globe ↔ globe rows, atm2obm, `scaled` coefficients (Q16), the bipolar mode
of the converter and the `yelmox_bipolar` gate ladder, any yelmox two-domain reference. To phase 4: the
global geometry hub, globe grids and coupler, footprints, rank-2 globe rows, `geo = vilma2`.
Out of scope: `ocean = data`, `ice = data` (the frozen-ice forcing tests of `bip_forcing`), bsl
`mixed` (stays refused), a second `ocean` instance.

## 2. Work packages

Each: what, decisions (recommendation in **bold**), depends on, gate. Sizes as
phase 2: S ≲ 300 lines, M ≲ 1000, L more.

### WP0 Infrastructure for the kryos-only gates — S

No yelmox commits, no yelmox baselines, no new stack.

- Stack S4 as it is (`gate-log.md` WP0); the kryos binary of the gate is built as for `nor`.
- Maps: none new expected (NH-32KM: hub = ice grid, `bilin_GLOBAL-1DEG_NH-32KM` made in phase 2;
  ANT-16KM ↔ ANT-32KM and the ANT climate maps from phase 1; `$G/maps/S4`); verify against the
  par file, `mapgen` R0 for anything missing, freeze.
- Data check (an agent, as `p2` Q17): every file the par file opens against the gate links; the NH
  and ANT bases are covered by `gate-reports/wp0-data-check.md`, so only the composition is new.
- The par files `par/kryos_NH-ANT.nml` (assembled, WP1) and, for the references of `sci-const`, the
  single-domain files with `[barysealevel] method = "const"`.
- Depends: phase 2 WP4 (`nor` conversion and gate). Gate: the data check; the single-domain gates
  `ant`, `nor` still pass on the gate binary.

### WP1 N domains in one run — M

1. **Par-file assembly** (`tools/src/Convert.jl`, a `merge` step over two converted files; the single-domain
   converter and its refusal of suffixed groups stay as they are): instance groups already carry the domain
   suffix; the other domain-specific groups (`topo`, `regions`, `snapshots`, `yelmo`, `ytopo`, `ycalv`, ...,
   `snapshot_*`, `isos`) that exist in both files are renamed with the suffix and the instance keys that name
   them are repointed; run-wide groups (`[kryos]`, `[clock]`, `[out_*]`, `[globe]`, `[geo_glob]`,
   `[barysealevel]`, `[forcing]` with the union of coefficients, `[coef_*]`) are taken once and must agree
   (a differing value stops, naming both); `[kryos] domains` lists both. **Decided (Q18): a merge step**;
   alternative (not taken): a hand-written file checked by `--check-par`. `--check-par`
   (`regression` §2.7, still missing) covers the merged groups.
2. First N = 2 run, local (the Mac test harness, grids of the local tests, 80 yr): the three checks of §3.2 at
   small scale (const equality with two single-domain runs; order swap; chunk 40 + 40 == 80).
3. `bsl` with N = 2: no code expected (gather, reduce and scatter exist, `reduce_bsl`); the **budget check**
   is new: a test procedure and `KryosTools` function recompute `bsl(t) - bsl(0)` from the geo exports
   (`dV_bsl` of every domain, `A_ocn`) and compare with the exported `bsl` (WP5).
4. Test of the stop for two producers of one derived coefficient (two climates with `dTo`): the message names
   the coefficient and suggests a coefficient of its own per domain (**keep the stop**, Q10). Two active
   rows into one import already stop (`contracts` §4.4).
5. Plugins (phase 2 WP1) per domain: the experiment lists are per domain (`[<domain>] experiments`); south
   `relax` + `opt` windows on elapsed time, continuous across chunks (phase 2 WP1 item 5).
6. Output tree: no kryos change (flat, unique instance names, `core` §6).

Depends: WP0; phase 2 WP1-3 for item 5 (`sci-plug`). Gate: `sci-const`, `sci`.

### WP2 smb_simple and climate `tsl_ann` — M (`smb_simple.f90` 985)

Phase 2 WP7, carried over (**decided, Q17: a package of phase 3, gated by unit tests
against yelmox's expressions plus `sci-smb`; phase 2 WP7 moves here**; the text below replaces it):

- Backend `surface = smb_simple` (private module, stateless, no restart).
- Climate vocabulary adds `tsl_ann` (snapshots `has_tsl_ann`,
  `src/backends/climate/kryos_climate_snapshots.f90:132`); `ta_ann`, `ta_ann_ref` only if the
  NH parameters use them (glacial scaling is off in every base).
- Mask: yelmox reads Yelmo `bnd%H_ice_ref` (bilin ice → surface) every call
  (`kryos_coupling.f90:481-483`). **Decided (`p2` Q8): import `H_ice_ref` from the ice
  export** (bilin ice → surface; follows recon, bit-identical); the yelmo adapter puts
  `H_ice_ref` every step (today once, static). With no `mask_file` (the NH case) this is the path
  that matters: after recon `H_ice_ref` is the recon ice (phase 2 WP4 re-publishes it).
- NH parameters: **the `[smb_simple]` group of `yelmox_North.nml`** (the NH base itself: `mask_file = ""`, mask from
  `H_ice_ref`; identical to the Greenland group; LIS differs only in `co2_const` 190, inert with `a3 = 0`, and its
  Batchelor `mask_file`/`mask_idx`). Both shipped values of `a1`, `beta0`, `c0`, `L_acc` are the LIS tuning ("tuned to NH smbpal (LGM)").
  The converter carries the group when `comps.surface_method = "smb_simple"` (`sci-smb` is `nor`/`sci` with that override).
- **Done (branch `p3-smb`)**: backend `kryos_surface_smb_simple.f90` (kernels `kryos_smb_simple`), `tsl_ann` export, `H_ice_ref` put
  every step, rows `climate.tsl_ann` and `ice.H_ice_ref` into the surface, converter, `tests/test_surface_smb_simple.f90`
  (bit-identical to the verbatim yelmox module) and the `smb_round_trip` of `test_run_sequence` (cont = chunk, GRL-16KM). `sci-smb` not run.
- Anchor (Q17): no yelmox run of smb_simple on NH-32KM exists and the LIS gate is deferred, so the gate is
  unit tests against the literal yelmox expressions plus `sci-smb`; a later LIS gate gives the run-level anchor.
- Depends: phase 2 WP4 (recon and `H_ice_ref`). Gate: unit tests, `sci-smb` (§3.2).

### WP3 Output precision — S

What exists (checked in `main`):

- Default writer (`component_write_output`, `kryos_component.f90`): ranks 0, 2, 3 written with
  the `prec` of the variable-table row (`input/kryos-variables-<component>.md`, column `prec`:
  `1e-3` absolute, `1e-4r` relative, `1e-3,1e-4r` both, `-` lossless) through ncio
  `prec_abs`/`prec_rel` (`variable_io` `var_io_write`, fesm-utils v1.4+). Every row of
  the eight tables has a value. A 0-D export without a row (a coefficient the par file names)
  is written lossless (`scalar_row`).
- Global override: `[kryos] write_lossless` → `nc_set_lossless` (`kryos_output.f90:115`,
  `kryos_run.f90:123`): ncio ignores all precisions. The converter maps `ctrl.write_lossless`
  to it (`yelmox2kryos.toml:478`); every gate `-p` list sets it, so the gates are unaffected
  and the bit-identity checks never see a rounded value.
- Restarts are lossless by construction: `fesm_fields_restart_write` writes `dp` without
  precision; the Yelmo restart forces lossless on its table; a rounded run restarts like a
  lossless one (yelmox verification 2026-10-21: rounded vs lossless restarts identical,
  max error / tolerance 0.977, size -31% ANT-32KM, -33% GRL-16KM;
  `~/.claude/projects/-Users-alrobi001-models-yelmox/memory/ncio-precision-rollout.md`).
- Backend writers: Yelmo's own tables (yelmo `topo-prepared` contains the ncio-precision
  merge a8d0f65b); the hub (`htopo_write_*` library, `f_ice` 1e-4 at `kryos_domain.f90:1197`).
  **Hard-coded**: FastIsostasy file `prec_abs=1e-3` (`kryos_geo_fastisostasy.f90:372`),
  snapshots (`kryos_climate_snapshots.f90:513-525`, `snapshots/kryos_snap_snapshots.f90:503`).

Package:

1. Audit of every file kryos writes against the yelmox values (the default list of the
   ncio-precision design: m 1e-3, m yr-1 1e-3 or 1e-4r, K 1e-3, fractions 1e-4, Pa 1e-4r,
   smb 1e-4... as in the tables): a list of the writers and their precisions in
   `docs/running.md`.
2. **Move the hard-coded precisions of backend streams into tables** (FastIsostasy `w`,
   `we`, snapshots, chion/hub where kryos owns the call), read by the backend's writer
   from a `kryos-variables-<component>-<backend>.md` like the instance tables; unlisted
   backend variable: lossless (as the default writer).
3. Time series (0-D, `*_ts.nc`, `forcing_ts.nc`, `timing.txt`), grids and debug output stay
   lossless (user decision 2026-10-21, yelmox); only listed rows round.
4. Tests (`test_output.f90`): each table row round-trips within its tolerance; the same
   run with `write_lossless` gives exact values; a restart bundle from a rounded run equals
   the lossless one; a table row without a unit or with a bad `prec` string stops.
5. **Gate `sci-nl`** (§3.2): the same `sci` run with `write_lossless` off and on, compared against
   each other: restart bundles identical, every rounded file within its row's tolerance, sizes reported
   (target about -30%). The gate `-p` lists keep `ctrl.write_lossless=True`; the lossless override
   is the key of the gate, no new mechanism.

Depends: —. Gate: unit tests; `sci-nl`.

### WP4 Forcing-increment bundles — S

`restart_every_df` (`coefs_restart_due`, `coefs_restart_fldr`): the forcing exports `restart_due` (binding
exists, `forcing_restart_due`), `run_step` writes a bundle `restart-<n>` (counter as yelmox,
`yelmox_bipolar.f90:286-290`). `core` §7 deferred it to this phase; 0 in every shipped config. Gate
`sci-dfm`: a ramp on `alpha` with `restart_every_df` set; a restart from one `restart-<n>` bundle continues
bit-identically (§3.2); also the first N = 2 run with a moving coefficient.

### WP5 Gate tooling for N domains — M

- `tools/gate/sci.toml` and variants (instances of both domains, `geo_glob`); a `Gate.jl` mode for
  kryos-against-kryos comparisons (no yelmox side, no variable map): per-domain comparison of a
  two-domain run against single-domain runs (instance files by name, the domain's folders of the bundles,
  globe files except the aggregates over domains), chunk and restart checks as phase 1.
- Domain order swap (`[kryos] domains` reversed) run and comparison of all files and bundles.
- `bsl` budget check (WP1 item 3) as a `KryosTools` function and a gate step.
- `--nl` comparison for `sci-nl`: lossless vs rounded (tolerance from the variable tables).
- `-p` with the merged group names (`barysealevel.method=const`, `comps_nor.surface_method=smb_simple`...).
- Registration of the reference files `K3-sci` (`regression` §3.5).
- Depends: WP0. Gate: the tool's own tests.

### WP6 Documentation — S

`docs/running.md` (two-domain runs, precision), `docs/build.md`, and the changes listed in §6.

## 3. Gates

### 3.1 References

No yelmox reference for the two-domain run. The references are the single-domain kryos runs of the same binary
(for `sci-const`), the run itself (for `sci`: its own `cont`), and the unit tests. Their yelmox anchor:
`ant` (and `ant-full`) and `nor`, gated K == R2/R2' in phase 1 and 2, which keep passing (WP1 does not touch
single-domain paths). After sign-off the `sci` and `sci-const` files and bundles at 1000 are registered as
`K3-sci` (kryos self-baseline, `regression` §3.5), the baseline for later kryos phases.

### 3.2 Gates and checks

Common `-p`: `ctrl.time_end=1000 ctrl.write_lossless=True tm_2D.dt=100 tm_2Dsm.dt=10 tm_1D.dt=10
tm_rst.dt=500`, 1 kyr from `time_init`, stack S4, gate flags (`-O2 -fp-model precise`). A two-domain run is
~15 min (NOR 11, ANT 3.5, `gate-log.md`).

| Gate | Setup | Checks |
|---|---|---|
| `sci-const` | `sci` with `[barysealevel] method = "const"`; references: the single-domain runs of NH (`nor` configuration) and ANT (`ant` configuration, relax/opt off), same binary, same `bsl` method | the domains are independent, so per domain **bit-identical**: every instance file, the domain's folders of the bundles at 500 and 1000, and the globe files except the aggregates over domains (`M_ice` of the forcing); 3 runs |
| `sci` | NH + ANT, `bsl` fastiso | `cont2 == cont`; `chunk == cont`, restarts included (`--restarts`, bundles at 500 and 1000); the same run with the domain order swapped in `[kryos] domains` identical in every file and bundle (order-independent reduction); **`bsl` budget closes**: `bsl(t) - bsl(0)` = `-Σ_d Σ dV_bsl_d / A_ocn` recomputed from the exports, to round-off of the `sp` sums (bound set from the first run); 4 runs (`cont`, `cont2`, `chunk`, swap) |
| `sci-smb` | `sci` + `comps_nor.surface_method=smb_simple` | tier B; WP2 unit tests |
| `sci-dfm` | `sci` + `alpha` ramp + `restart_every_df` | tier B, with a restart from one `restart-<n>` bundle |
| `sci-nl` | `sci` with `write_lossless` off | restart bundles identical to those of `sci`; every rounded file within the tolerance of its table row; sizes reported |
| `sci-plug` | the shipped south, `relax` and `opt` on; **after phase 2 `ant-full`** | tier B |

Tier B: a `cont` run (bundles at 500 and 1000) and the second half restarted from the 500 bundle (K chunk ==
K cont, bit for bit), no stop, `bsl` budget closed. About 17 runs in all. The expected-error and unit checks
of §3.4 run on the Mac.

### 3.3 Expected identity and accepted differences

| Where | Why | Handling |
|---|---|---|
| `sci-const` globe files | `M_ice` of the forcing aggregates two domains | excluded from the per-domain comparison (listed in `sci.toml`) |
| `sci-const` vs `ant`/`nor` gates | the single-domain gates run `bsl` fastiso (own volume only), `sci-const` runs `const`: the chain is code identity plus the `bsl` source; no yelmox `const` run (accepted, Q18) | Q18 |
| `sci-nl` | rounding by design | within the table tolerance; restarts exact |
| `bmb_shlf` ties, ANT shelf 16 km ↔ ice 32 km | `adapters` §1.2 (did not occur in phase 1) | registered if they occur |
| `timing.txt` | not comparable | not compared |

A difference not in the table is a failure with a named cause (`regression` §5).

### 3.4 Tests without a yelmox counterpart

Unit tests (Mac): two-domain run, domain order swap, const equality at small scale (WP1); the `bsl` budget; the stop for two
`dTo` producers; smb_simple expressions against the literal yelmox code on random `sp` fields (as
`test_conversions`); precision round trip (WP3); increment bundle round trip (WP4).

### 3.5 What phase 4 needs from phase 3 (checked against `p4` §10)

| `p4` §10 | Need | Phase 3 |
|---|---|---|
| N10 | two domains in one run (NH + ANT, disjoint, name order), `domain_check_reference` per domain, NH recon and the vm-m16 fix | **provided** (WP1, `sci-const`, `sci`; recon and fix from phase 2) |
| `p4` G1 | no regression of fastiso configs; the `dphi` change and the split cold start bit-neutral | **provided as kryos self-baselines**: `K3-sci` (files and bundles of `sci` and `sci-const` at 1000, registered at phase 3 sign-off). `p4` G1 reads: K(S5, `vilma=1`, fastiso configs) == K(S4) for grl, ant, ant-full, nor **and `K3-sci`** (no yelmox baseline, no phase 3 bipolar gate) |
| `p4` G6 | drift report against the fastiso baseline of the science case | against `K3-sci` (`sci`); no obm in G6 |
| `bsl` | order-independent reduction for N domains (fastiso runs stay; `reduce_bsl` returns for vilma2) | **provided** (WP1 item 3: budget and order swap) |
| N1, N2, N3, N4, N5, N7, N8 | globe coupler and grids, global instance grid, rank-2 fields on a global geo, rank-2 globe ↔ domain rows, footprints, generic global step, lon-lat writer | **not provided**: `p4` WP4 builds them (its own fallback column) |
| N6 | global geometry hub (spliced reference and merged ice as one function) | **not provided**: `p4` §5.1 item 3 builds it for VILMA2; a hub as a product of the globe has no consumer before obm |
| N9 | `GLOBE_ORDER` generalised | **not provided**: today `forcing, geo` (`kryos_restart.f90:78`); `p4` adds the global geo's names, `ocean` joins Later |
| (`p4` §2.1) | `ocean_glob` instance, `[globe] ocean`, obm in the science case | **dropped** (also in `p4`, §6) |

Phase 4's WP4 grows accordingly (it already lists "N1-N7 items phase 3 does not cover").

## 4. Order

| # | Package | Size | Unblocks |
|---|---|---|---|
| 1 | WP0 infrastructure (parallel with all) | S | every gate |
| 2 | WP3 output precision (independent) | S | `sci-nl` |
| 3 | WP1 N = 2 run, par-file assembly, budget | M | `sci-const`, `sci` |
| 4 | WP5 gate tooling (parallel with 3) | M | all gates |
| 5 | WP4 increments | S | `sci-dfm` |
| 6 | WP2 smb_simple (needs phase 2 WP4) | M | `sci-smb` |
| 7 | WP6 docs | S | sign-off |

Critical path: phase 2 sign-off (`nor`, recon) → WP0 → WP1 → `sci-const`, `sci`. Phase 2 WP2 + WP3
(relax/opt, gate `ant-full`) and WP4 (recon, gate `nor`) are done before phase 3 `sci` (Q19), so
`sci-plug` needs no wait; the plugins stay off in `sci-const` and `sci`. Sign-off: `sci-const`, `sci`, `sci-smb`,
`sci-dfm`, `sci-nl`, `sci-plug` green and `K3-sci` registered; `grl`, `ant`, `ant-full`, `nor` still green.
Then phase 4.

## 5. Decisions

**Decided** (Alex, 2026-10-29, in the style of `phase2-plan.md` §5; the numbers keep the
first draft's; 16-19 decided on review):

1. **Decided: no yelmox baseline for the two-domain run: no yelmox pin, commits, par file or stack.** The
   gates are kryos self-consistency anchored on the single-domain gates (§3). Was: yelmox dev pinned at a
   commit on the v24 pins or a new stack, a chain of yelmox baselines with preparatory yelmox commits.
2. **Out of scope: where obm lives** (obm is Later, §7).
3. **Decided: `bsl` with two domains is the order-independent reduction of design decision 5, gated by
   the domain order swap and the budget** (WP1 item 3, §3.2). Was: a preparatory yelmox commit
   making yelmox's reduction order-independent (superseded: no yelmox reference).
4. **Decided: a reduced gate set**: `sci-const` (equals two single-domain runs), `sci` (determinism, chunk,
   order swap, budget), single runs for variants (§3.2). Was: the ladder of eight setups; then full checks on
   one configuration plus yelmox-referenced variants.
5. **Decided: shipped south with relax/opt off first; the plugins-on rung after phase 2
   `ant-full`** (`sci-plug`). Was: the same, or only the shipped base.
6. **Out of scope: the freshwater port** (Later, §7).
7. **Out of scope: the box → ice path** (Later, §7).
8. **Out of scope: the global background reference** (phase 4).
9. **Decided: the whole global geometry hub moves to phase 4** (§3.5). Was: grids,
   footprints and assembly in phase 3.
10. **Decided: two producers of one derived coefficient: keep the stop** (WP1 item 4). Was:
    keep the stop, or a per-domain implicit suffix.
11. **Decided: `ice = data` and `ocean = data` are out of scope** (§1.4). Was: non-goals or
    `ice = data` brought forward.
12. **Decided: the science configuration is NH (`yelmox_North`, NH-32KM) + ANT (32 km ice grid), obm off,
    freshwater coupling off; no yelmox run of it.** The first answer (a yelmox bipolar dev run with a new par
    file, two preparatory commits) was superseded the same day; its findings are kept in Later (§7). Was: a bipolar par file
    with NH + ANT as the first rung of phase 4.
13. **Decided: forcing-increment bundles in phase 3** (WP4). Was: implement or defer.
14. **Decided: `component_t` stays in kryos** (unless CLIMBER-X asks). Was: stay or promote.
15. **Decided (conditional): `scaled` coefficients only if the science-case baseline par file
    uses them (`dTa_*` from gmt); otherwise deferred.** The par file of §1.3 does not
    (shipped forcing of both bases), so `scaled` is Later; see Q16.

16. **Decided: forcing of the science configuration is each base's shipped forcing; `scaled` deferred.**
    North `alpha` const 1 + derived `dTo` (terms `lgm`, `ocean`), ANT none; no `scaled`, fewest new
    keys; `scaled` stays Later (§7, an S package: `libs/kryos_forcing.f90:139-147,217-245`; needs a
    present-day reference snapshot for an NH domain that starts from an LGM reconstruction). Was: (a)
    that (recommended), or (b) the bipolar gmt-driven forcing (`dTa_*` = factor × gmt), which
    needs `scaled`.
17. **Decided: smb_simple is a package of phase 3 (WP2), gated by unit tests against yelmox's
    expressions plus `sci-smb`** (no run-level yelmox anchor, §WP2). Was: (a) that
    (recommended), or (b) phase 2 WP7 un-deferred and gated on LIS as written (LIS maps, data and
    baseline, i.e. phase 2 WP0 for LIS; gives the yelmox anchor).
18. **Decided: the two-domain par file is assembled by a merge step (WP1 item 1), and the
    const-`bsl` anchor is accepted (no yelmox `const` runs).** `const` only removes the gather and the
    reduction, and the unit tests exercise it, so the transitive anchor holds for everything except
    the `const` source (§3.3). Was: (i) a merge step (recommended) or a hand-written file checked by
    `--check-par`; (ii) accept the anchor (recommended), or one single-domain yelmox run per base with
    `barysealevel.method=const` against the existing yelmox code (`topo-prep`, a par switch, no new
    commit) and a kryos run of the same.
19. **Decided: phase 2 WP2 + WP3 (relax/opt, gate `ant-full`) and WP4 (recon, gate `nor`) are done
    before phase 3 `sci`** (§4). Was: `sci-plug` waiting for `ant-full` while the rest ran with the
    plugins off (decision 5 keeps the plugins-off first rung).

## 6. Changes needed in other docs (applied in the commit "docs: phase 3/4 edits")

- `design.md` §14 phase 3: "N domains; `smb_simple`; output precision; forcing-increment
  bundles. Target: NH-32KM + ANT in one kryos run, gated by kryos self-consistency against the
  single-domain gates (no yelmox reference)"; the `ocean = obm`, freshwater and bipolar target
  (`yelmox_bipolar`) move to a later phase (after phase 4); the globe grid/hub parts of §6 (last
  paragraph) and §14 phase 4 stay with phase 4; §6 mechanism: the bsl reduction with two domains is
  order-independent and gated by the order swap and the budget; §5 table and decision 10 (ports
  are domain-side scalars) stay as the Later design; §15 decision 9 settled (`component_t` stays).
- `core-design.md` §7: forcing-increment bundles done; §6 "As implemented": precision
  tables for backend streams, `write_lossless` the lossless override; §3.4 sequence and §2
  `globe_t` unchanged (blocks 3b, 5b, globe rows are Later).
- `contracts-design.md` §7: the two-producers stop kept (test added); §2 vocabulary
  `tsl_ann` (WP2); §10 row OBM "Later".
- `regression-design.md`: §2.5 the suffixed multi-domain refusal stays (the bipolar mode of the
  converter is Later; two-domain par files are assembled from single-domain ones); §1.3 and §3: a
  kryos-against-kryos gate type (per-domain equality with single-domain runs, order swap, budget) beside
  the yelmox-referenced ones; §3.5 the `K3-sci` registration; §9 row 9: the non-lossless pair `sci-nl`.
- `phase2-plan.md` top "Deferred" list: WP7 (smb_simple) moves to phase 3 WP2 (Q17);
  §1.3 row 3 (phase 3 list): "N domains, smb_simple, increments, output precision; obm, freshwater,
  `yelmox_bipolar`, bsl `mixed`, two producers, suffixed multi-domain par files → Later (bipolar/obm)";
  §4 sign-off line "then phase 3 (N domains: NH + ANT) and phase 4 (`geo = vilma2`)"; WP13 "MEAN
  imports for slower consumers: leave for phase 3" → no consumer (the freshwater port is Later).
- `phase4-plan.md` (applied there): §0 item 7 and §10 header ("assumed
  from phase 3: globe grids, gather/scatter, global hub") → only N10 and the `bsl`
  reduction come from phase 3 (§3.5); N1-N9 are phase 4's (WP4 grows); §2.1 `ocean_glob`,
  §3.4 `ocean = "ocean_glob"` and `[globe] ocean`, G6 ("+ obm + bipolar forcing") and the
  title line ("with `obm`") drop obm; **G1 "the phase 3 bipolar gate" → the registered phase 3
  self-baselines `K3-sci` (files and bundles of `sci` and `sci-const` at 1000, kryos commit of the
  phase 3 sign-off): K(S5, `vilma=1`, fastiso configs) == K(S4) for grl, ant, ant-full, nor and `K3-sci`**;
  G6 drift report against `K3-sci`; §8.3 `GLOBE_ORDER` "the phase 3 names (forcing, ocean, geo)" →
  forcing, geo.

## 7. Later (bipolar/obm)

Not phase 3, not scheduled; the findings of the first draft, kept so the work is not lost.
Refs as in the first draft (`yb` = dev e31c408).

| Item | Findings and references |
|---|---|
| `ocean` component + `obm` backend (L) | `src/components/ocean/` (the `kryos_forcing` pattern), backend = port of `libs/obm/` (693 lines: `obm.f90` 253, `nautilus.f90` 257, `obm_defs.f90` 76, `ice2ocean.f90` 74, `ocean2ice.f90` 33), all `sp`. Contract (0-D, global): imports `dT_<role>` from forcing coefficients named by `[obm_forcing]`, `fwf_<box>` from the domains; exports `T_north`, `T_south` (+ `_init`) to the forcing, the rest output only; boxes are the vocabulary. Bit-identity hazards ported literally with a unit test per expression: `tau_star_south` mixed `sp`/`dp` (`(a/tau)*dt`, `obm_coupling.f90:361-368`), the hysteresis `factor` (`rate*dt`, `r8_normal_ab` noise from `random_number`: global state, per-instance seed as `p2` WP9 item 5, `hyst_on = False` in every setup). Restart `ocean_restart.nc` = `write_obm_restart`; `GLOBE_ORDER` (`kryos_restart.f90:78`) becomes forcing, ocean, geo. Output `<instance>_obm.nc` keeps `write_obm_update`'s variables. Constraints: `max_instances = 1`, global, no build flag. Where obm lives: copy at the pin, yelmox frozen once the gate passes (recommended); a shared repository in the configme stack; or diverging copies (not recommended). Sergio changes `libs/obm` (`639516e`, `fc4f247`). |
| Freshwater coupling | `couple_fwf_north/south` (`obm_coupling.f90:98-99,330-333`), `calc_fwf` (`ice2ocean.f90`; `dVdt = -sum(dHidt)*dx*dy*conv`, whole domain; `mb`, `ablation` inside the hydrographic mask, `{grid_name}_HYDRO-BASINS_zickfeld2004.nc`, read whenever `couple_fwf_*` is on). Port options: (a) a core **domain port**, a 0-D `fwf` [Sv] from the ice exports (`dHidt`, `mb_net`, `smb`, `bmb`, `cmb`, `f_grnd`, `dx*dy`), row `domain.fwf → ocean.fwf_<box>` with `<box>` = key `obm_box` of the domain group (recommended); (b) a 0-D export in every ice adapter; (c) the backend reads 2D exports. Precision: narrow to `sp`, sum in `sp` in column-major order with yelmox's operand order; divide the `rho_ice`-scaled exports by `rho_ice`. `dHidt`: yelmox reads the instantaneous value, kryos exports a MEAN: equal when the ice step equals the coupling step (`p2` WP3 item 4; stop at init otherwise). |
| atm2obm, globe ↔ globe rows | `coupling_atm2obm` (`:347-372`): `tstar_i = tstar_i_init + dT_i`, vapour fluxes from `gmt`. Involves no domain: `globe_t` gets its own resolved rows (rows live in each domain's table today, `kryos_globe.f90`), resolved per run so a zero-domain run (`bip_obm_only`) works. |
| Box → ice and the climate `refresh` | `coupling_obm2ism` (`:402-417`): box `tn`, `ts` + init values → derived coefficients `obm_north/south` (`coefs_set_obm`, `libs/kryos_forcing.f90:278-293`) → ocean terms re-applied (`climate_apply_ocean`, `libs/yelmox_climate.f90:347-358`) after both climate updates, before the shelf. Port: a second producer kind `ocean` (exports `T_north`, `T_north_init`, the forcing applies `anomaly` in `sp`), a new `refresh(st)` binding on `component_t` (climate snapshots adapter: `terms_eval(derived=.true.)`, `terms_apply_ocn`, re-put `to_ann`, `so_ann`, `dto_ann`), a second `fields_put` of an INSTANT export at the same `t1` (to check in `fesm_fields`). Alternatives: ocean exports a shelf-base field (not identical); ocean before climate (changes lags). Sequence: block 3b (`ocean_update`, inputs as received at the end of the previous step), block 5b (ports, rows, ocean → forcing, `climate_refresh`), per-domain loop split into A (ice, hub, climate, surface) and B (shelf); cold start: ocean `PASS_INIT` and the `ocean → forcing` rows at `init` before any domain; restore: ocean after forcing. |
| `scaled` coefficients (S-M) | factor × a tsgen coefficient, restart-free (`libs/kryos_forcing.f90:139-147,217-245`): `dTa_*` = factor × `gmt` (`input/forcing_db.nml`); ports to `kryos_forcing_tsgen.f90` and the `snap_coef_table_t`; `coefs_varies`, `coefs_check_used` (init stop for unused coefficients). Needed when the science-case par file is gmt-driven (Q16 b). |
| Init checks of the bipolar setup | `obm_coupling_init` (`:139-233`): box coefficient of the other hemisphere, box and air ocean terms in one domain (`check_domain_ocean`, `:235-277`), box on but atmosphere varying without `atm2obm` (`check_box_forced`, `:279-296`), unused coefficients; named messages so the `chk_*` cases of `run_tests.sh` carry over. |
| Bipolar mode of the converter, `yelmox_bipolar` gate ladder | Converter: lifts `refuse_group` for suffixed groups (`tools/data/yelmox2kryos.toml:87-98`), domains from `domain_north`/`domain_south`, `[ctrl] active_north/south` → `[kryos] domains`, obm keys accepted only when off. Ladder: base `yelmox_bipolar_Bipolar.nml`, setups of `params.sh`: `bip_uncoupled`, `bip_obm_only` (zero domains), `bip_obm_to_ice`, `bip_gis_to_obm`, `bip_ais_to_obm`, `bip_gis_obm`, `bip_ais_obm`, `bip_full`, modifiers `bip_ocean_replace`, `bip_lag_south`, `bip_gmt_ramp` (`coef_gmt.method=ramp-time f_min=-6 f_max=0 df_sign=-1 dt_ramp=60`: gmt is `const` 0 in the base, so a ramp moves the box); Greenland GRL-16KM + ANT-32KM (hub ANT-32KM there); rungs a-h/full/var/dfm of the first draft (commit 9b942c8 of this file), tier A/B, ~100 runs, yelmox baselines on a prepared yelmox branch. `bip_forcing` (frozen-ice tests) needs `ice = data`: out of scope. |
| yelmox two-domain reference (dropped 2026-10-29) | Findings from reading yelmox dev e31c408: `yelmox_bipolar` runs NH + ANT with obm and freshwater off as a par-file change (the base par file is the uncoupled setup `bip_uncoupled`; `obm_coupling_init` checks only what is on); the par file needs north groups from `yelmox_North.nml`, south groups from `yelmox_Antarctica.nml` (the bipolar `domain_south` has hub ANT-32KM), per-domain physics groups (`nml_ytopo`... repointed in `[yelmo_<sfx>]`), `[forcing] coefs = alpha dTo`. Code: reference preparation (`prepared = .true.` at `libs/kryos.f90`, cherry-pick of `topo-prep` 74f3a2b, 5ac89c8, dev changed 16 lines there) and an order-independent `bsl` reduction (each domain's `isos_update` on a copy of `bsl`, then `bsl_now -= (Σ ΔV in sp, name order)/A_ocean_now`, the expression of `kryos_geo_fastisostasy.f90:465`; yelmox mutates the shared `bsl` sequentially, `libs/kryos_coupling.f90:119-150`). v24 pins already contain the ncio precision code (fesm-utils v1.4, yelmo v2.4); dev output `forcing_ts.nc` replaces the coefficients of `yelmo_ts.nc` (`df48108`). |
