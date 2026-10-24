# Regression harness and configuration migration — design note

Status: proposal (phase 1), no code. Companion notes: `core-design.md` (par-file
schema, instance groups, defaults files) and `adapters-design.md` (backend
adapters, unit conversions), written in parallel; where this note names a group,
key or file that those notes own, it is marked *(core)* or *(adapters)* and
their choice wins.

Scope: how kryos phase 1 is shown to reproduce `yelmox_Greenland` and
`yelmox_Antarctica` bit for bit (design.md §14), and how that stays true: the
reference runs, the conversion of yelmox par files, the comparison tool, the
automation on Levante, and the classification of differences.

**Hard constraint.** "Bit-identical" means every compared value is equal bit
for bit, in every output file and restart bundle, at every record, between
yelmox and kryos built from the same dependency commits with the same compiler
and flags. Anything else is a difference with a named cause (§5), registered
(§3.4) or fixed.

Facts this note relies on (checked 2026-10-24):

- yelmox `v2.4` pins fesm-utils `v1.4`, yelmo `v2.4`, FastIsostasy
  `v0.3-fesmc`, rembo1 `v1.3`, chion `v0.2.0`, insol `v1.2.0`. yelmox dev
  after v2.4 changes only docs and the manifest.
- fesm-utils `phase0-kryos` branches from `v1.4` (27 commits ahead).
- Levante: host alias `levante` (user a271116), work tree `/work/ba1442/robinson/`,
  FesmData products at `/work/ba1442/fesmdata`. Production flags for ifx on
  Levante are `-Ofast -march=znver3` (configme `dkrz_levante.mk`).
- This Mac: yelmox is built (`yelmox.x`, gfortran, 2026-10-19), with the yelmo,
  fesm-utils, FastIsostasy and chion archives. The FesmData v2 products the
  bases read are **not** here (`~/data/fesmdata`: ANT-16KM topography, ANT-32KM
  PMIP, the test domains, NH-16KM, Global; no Greenland; no `fesmdata` link in
  yelmox); `ice_data` is the v1 layout; `isostasy_data` is complete. So the
  full gate runs on Levante only, until the products are synced (§4.5).
- Yelmo, fesm-utils (`htopo`, coupler sp path) and marine_shelf work in `sp`;
  the sp `map_field` accumulates in `dp` and rounds once
  (`map_field_grid_grid_sp`).


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

## 1. Reference runs

### 1.1 Configurations

Two gate configurations, both from yelmox tag `v2.4`, cold start at forcing time
0 (`tstep_method = "const"`, as in the bases), `dtt = 10`, `dt_clim = dt_isos =
10`, so every component updates every step:

| Gate | Base | `-p` |
|---|---|---|
| `grl` | `yelmox/yelmox_Greenland.nml` | `ctrl.time_end=1000 ctrl.write_lossless=True tm_2D.dt=100 tm_2Dsm.dt=10 tm_1D.dt=10 tm_rst.dt=500` |
| `ant` | `yelmox/yelmox_Antarctica.nml` | the same + `sim.relax=False sim.opt=False` (Q1) |

Why these numbers. 1000 yr is 100 coupling steps with chion's daily cycle,
two restart writes and the strong transient of a cold start; a 1-ulp
difference anywhere is amplified by Yelmo's adaptive time stepping within tens
of steps. Cost: grl-paleo chunks do ~5 kyr per 30 min on 16 threads, so one
gate run is ~10 min. `tm_2Dsm.dt=10` writes
`yelmo_sm.nc` every step, which locates the first divergence to one step.
`write_lossless` is required: per-field bit rounding (`nc_bitround`) hides
small differences and does not commute with a unit conversion (§3.2).
Timeouts are counted from the start of the run in yelmox and on absolute
time in kryos (decision 8); with `time_init = 0` and chunk boundaries that are
multiples of every `dt` they coincide.

`ant` drops relaxation and optimisation because they are experiment plugins
(`src/experiments/relax`, `opt_cb_ref`), which design §14 puts in phase 2. The
shipped base (`relax = opt = True`) becomes the phase 2 gate. `init_method =
"equil"` stays: it is a cold-start option of the yelmo backend
(`yelmo_update_equil`), not an experiment.

The base length (15 kyr) is run once at the end of phase 1, not per commit.

### 1.2 Run set per gate

| Run | What | Purpose |
|---|---|---|
| `cont` | 0 → 1000 in one job | the comparison target |
| `chunk` | 0 → 500, then 500 → 1000 from `restart-0.500-kyr` (yelmox: `scripts/run_chunks.sh`) | restart round trip: records after 500 must equal `cont` |
| `cont2` | `cont` again, once per new stack | determinism (same binary, same threads): must equal `cont`, else no gate is meaningful |

kryos runs `cont` and `chunk`. The comparisons are `K.cont` vs `Y.cont` (the
gate), `K.chunk` vs `K.cont` and `Y.chunk` vs `Y.cont` (restart exactness of
each model), and the restart bundles at 500 and 1000 of `K.cont` vs `Y.cont`.
A kryos run never starts from a yelmox bundle (different layouts); the round
trip is tested inside each model.

`chunk` vs `cont` (`Compare.compare_chunks`, the gate's `chunks`): chunk `k`
covers `(t_k, t_k+1]`, the first chunk also `t_0`. Every `cont` record in that
interval must be in chunk `k` and equal (E), and the `cont` bundle at `t_k+1`
must equal the bundle chunk `k` writes at `t_k+1`. The record at `t_k` and the
bundle a chunk writes at its own start are not compared: they are written
before the first step, and the diagnostics in them that a restart does not
restore are still being recomputed (yelmox v2.4: `visc_bar`, `mb_clip_tot`, `pc_tau*`,
`z_bed_filt`, … at t = 500; prognostic state exact, gate log 2026-10-24).

### 1.3 Dependency versions: the same fesm-utils for both

yelmox pins fesm-utils `v1.4` (stable); kryos needs `phase0-kryos`. The gate
builds both models against **one** set of dependency checkouts (one commit of
each library, built once, linked into both installs), so the reference chain
is:

1. **R0**: yelmox `v2.4` with its own tag pins (stack `v24`).
2. **R1**: yelmox `v2.4` with fesm-utils at commit `F` of `phase0-kryos`, every
   other library at its v2.4 tag (stack `S`). Gate 0: `R1 == R0`, every file,
   bundle and variable. This proves that the phase 0 changes do not alter
   yelmox's answers for these two configurations.
3. **K**: kryos built on stack `S` (the same checkouts and archives). Gate 1:
   `K == R1`.

4. **R2** (decision 2026-10-24, Alex): Yelmo is fed its topography as given,
   as in CLIMBER-X. The cleaning Yelmo did itself (bed roughness scaling,
   englacial lakes, thin ice, bed gradient limit) moves into the hub
   (fesm-utils `htopo`, `[topo]` keys `z_bed_f_sd`, `remove_lakes`,
   `H_ice_min`, `grad_lim_zb`, all off by default), done once on the hub
   grid, and Yelmo gets `yelmo_init(..., prepared=.true.)`. Stack `S2` = `S`
   with those fesm-utils and yelmo commits (both off by default: `R1(S2) ==
   R0` first); R2 = yelmox `v2.4` + a preparatory commit that turns them on
   (`[topo]` of the gate configs, `prepared=.true.`). Gate 1 becomes
   `K == R2`; `R2` vs `R0` is documented (one thin-ice rule instead of Yelmo's
   two, the gradient limit on the hub grid: Antarctica 16 km instead of 32 km).

Whenever a library moves (a new fesm-utils commit, or a yelmo change kryos
needs), the stack gets a new id, R1 is rerun on it and must equal R0 (or the
change is a registered deviation, §3.4). If yelmox v2.4 does not build against
`F`, the fix goes into fesm-utils (aliases, restored signatures), never into
yelmox.

What phase 0 has already proven, and what it has not:

| Change in `phase0-kryos` | Evidence | Status |
|---|---|---|
| coupler/`map_field` apply surface (containers, pass-through, `norm`/`frac2` absent) | `test_coupler_bitident` + `bitident.sh` against the pre-change head; cached maps 409/409 identical | proven |
| cell edges from one source (`x_bnds`) | `grid_area` bit-identical (CHANGELOG) | proven |
| identity shortcut `remap(src == dst)` | not identical for a caller that self-remaps; yelmox short-circuits before the coupler (`c5df773`) | proven by grep only; R1 confirms |
| polar stereographic on a sphere | GRL/ANT grids are WGS84 (`semi_major_axis = 6378137`), unchanged | proven |
| legacy `varslice` path | `test_varslice_legacy`, byte-identical | proven for groups with a `time` key; a group without one now takes the file's time axis |
| `fesmdata_flux_conv`, `temp_ref` via `units_conv` | checked bit for bit | proven |
| **`nml` string vectors split at blanks outside quotes; comment after a group name** | none on yelmox files | **unproven**: every vector key of the staged par files, `input/forcing_db.nml`, `input/clim_db.nml` and both defaults files must parse the same (§4.3, `test_nml_compat`) |
| **newly generated maps** (prune `cutoff·max(w)`, spherical clip fixes) | changes generated maps by design | **avoided**: no map is generated in R1 or K (frozen maps, §5) |
| **ncio restart read** (`|x| < 1e-30` → 0, NaN → -9999) | ncio unchanged between v1.4 and `F` | same convention in both models; new kryos restart state is covered by the `chunk` runs |

So the recommendation is: rebuild yelmox v2.4 against `F`, show `R1 == R0`
before kryos is compared to anything, and use R1 as the reference.

### 1.4 Where the references live

Under `$G = /work/ba1442/robinson/kryos-gate` (a variable, `GATE_ROOT`):

```
$G/stacks/<stack>/            yelmox and kryos installs, shared dependency checkouts
$G/maps/<stack>/              frozen map cache (read-only, §5)
$G/ref/<gate>/yx-<sha7>/<stack>/{cont,cont2,chunk}/    yelmox runs (R0: stack v24)
$G/runs/<gate>/kx-<sha7>/<stack>/{cont,chunk}/         kryos runs
$G/reports/<date>_kx-<sha7>_<stack>.{md,json}
```

`<stack>` is the first 8 characters of the hash of the resolved
`tools/gate/stack.toml` (§4.1); `v24` for R0. A reference is made read-only
(`chmod -R a-w`) once its gate passes. References are kept until the end of
phase 3.

### 1.5 Optional ladder

If the core allows a domain without `geo` or `shelf` (the hub reference as the
static bed, *(core)*), the gate can also be climbed in steps against yelmox
runs with the same components off (`comps.with_isostasy=False`,
`comps.with_shelf=False`): ice + climate + surface, then + shelf, then + geo.
Useful while components land; not required for the phase 1 sign-off.

## 2. Configuration migration

### 2.1 Tool

A Julia package in `kryos/tools/KryosTools/` (Julia ≥ 1.10), modules
`Namelist` (reader), `Convert` (yelmox → kryos), `Compare` (§3), `Gate`
(staging, manifests, reports). Dependencies: NCDatasets, TOML and Printf
(stdlib), JSON; CairoMakie as a weak dependency (package extension, as in
Pagos.jl), so the comparison runs on Levante without loading Makie. Scripts in
`tools/bin/`: `yelmox2kryos.jl`, `kryos-compare.jl`.

The namelist reader exists: `yelmox/scripts/forcing-builder/Namelist.jl` reads
files as fesm-utils `nml.f90` does (comments kept) and `Resolve.jl` resolves a
base + library + defaults + runme `-p` (`g.k=v` scalar, `g.k==a,b` vector), with
the source of every value. Proposed home: FesmUtils.jl (`FesmUtils.Namelist`),
used by the forcing builder and by KryosTools (Q5). It must follow the new
quoted-vector rule of fesm-utils `nml`.

```bash
julia --project=tools tools/bin/yelmox2kryos.jl --yelmox $YELMOXROOT \
    --par $G/ref/grl/yx-<sha>/<stack>/cont/yelmox_Greenland.nml \
    --out par/gate/grl.nml --report par/gate/grl_conversion.md
```

### 2.2 Input: the staged par file

The converter reads the par file **as staged by runme in the reference run
directory** (what yelmox actually read, `-p` already applied), plus `input/` of
the yelmox commit (defaults, `forcing_db.nml`, `clim_db.nml`). So runme's
override semantics are not reimplemented on the gate path. A `base + -p` mode
(through `Resolve.parse_overrides`) exists for interactive use.

### 2.3 Output: fully expanded

The converted par file sets **every** key of every group it writes, from the
resolved yelmox value (par file, else library, else yelmox's defaults file).
Library groups referenced by name (`snapshot_*`, `coef_*`, `term_*`, the
`clim_db` varslice groups) are copied in. The gate par file is therefore
self-contained and independent of kryos' defaults files, which then only serve
`nml_validate`: a kryos default that differs from yelmox's cannot change a gate
answer (the v2.4 sweep found the yelmox copy of `yelmo_defaults.nml` out of
sync with yelmo, A1; this removes that class). A `--keep-library` mode writes
references instead, for user-facing conversions.

A header records the yelmox commit, the source file, the `-p` list, the
converter version and the sha256 of every input file the par file names
(series, timeout files, `yelmo_phys_const.nml`); the gate compares these hashes
between the two run directories (§4.2).

### 2.4 Mapping rules

Group and key names on the kryos side are *(core)*; the rules are:

| yelmox | kryos |
|---|---|
| `[ctrl]` timeline (`tstep_method`, `tstep_const`, `time_init`, `time_end`, `dtt`, `calendar*`), `restart`, `write_lossless` | run/clock group *(core)*; `dtt` → run step and the `ice` instance `dt` |
| `[comps]` `with_*`, `climate`, `surface_method`, `dt_clim`, `dt_isos` | the domain's instance list; backend of each instance; `dt` of `climate_<d>`, `geo_<d>` |
| `[domain]` `name`, `phys_const`, `grid_*`, `basins` | domain group `[Greenland]`; instance `grid` (`grid_ice = ""` → `grid_hub`) |
| `[topo]`, `[regions]` | domain hub reference and regions *(core)* |
| `[tm_2D]`, `[tm_2Dsm]`, `[tm_1D]`, `[tm_rst]` | output schedules (2D, Yelmo small set, `_ts`), restart schedule (`schedule_t` keys = `timeout` keys) |
| `[output]` `write_*` | per-instance output switches |
| `[yelmo]`, `[ytopo]`, `[ycalv]`, `[ydyn]`, `[ytill]`, `[yhyd]`, `[ymat]`, `[ytherm]`, `[yelmo_masks]`, `[yelmo_init_topo]`, `[yelmo_data]` | physics groups of `ice_<d>` (yelmo), verbatim |
| `[sim]` `init_method` (`none`, `equil`), `init_equil_time`, `init_marine_H`, `init_kill_shelves`, `init_time_thrm` | yelmo backend cold-start keys *(adapters)* |
| `[isos]` | physics group of `geo_<d>` (fastisostasy), verbatim |
| `[barysealevel]` | `bsl` source of `geo` (`own | file | const`) *(core)* |
| `[forcing]`, `[coef_*]`, `[term_*]`, `[snapshot_*]`, `[snapshots]`, `[orbit]`, clim_db groups | physics groups of `climate_<d>` (snapshots), schema unchanged; owner of the coefficients in phase 1 *(core)* |
| `[chion]`, `[surface_chion]`, `[itm_chion]`, `[pdd_chion]`, `[chion_const]` | physics groups of `surface_<d>` (chion), verbatim |
| `[marine_shelf]` | physics group of `shelf_<d>` (marine_shelf), verbatim |
| `[ghf]`, `[sed]` | `geothermal_<d>`, `sediments_<d>` with backend `data` (`datasrc` groups: `obs_path` → static field; `use_obs = False` → const) |
| `{domain}`, `{grid_name}` in paths | kept; resolved by the kryos path service *(core)* |

Instance suffixes: `grl`, `ant` (from `[domain] name`).

Groups of inactive options are dropped with a line in the report
(`[smbpal]`, `[smb_simple]`, `[itm]`, `[pico]` when unused; `[relax]`, `[opt]`,
`[negis]` when their switch is off).

### 2.5 Refusals (phase 1)

The converter stops, naming the key and the phase that brings it:
`climate` ≠ `snapshots` (esm, rembo: phase 2); `surface_method` ≠ `chion`
(phase 2); any `with_* = False` (data backends, phase 2); `sim.relax`,
`sim.opt`, `use_negis`, `scale_glacial_smb`, `lim_pd_ice` true (plugins, phase
2); `init_method = recon | recon_ref`; `[ice_bnd] active`; `[bias]`/`[taper]`
active; `bmb_method = "pico"` or `*-slope` (broken, sweep E30); `write_cmip =
True` (CMIP post-processing table, later); `ctrl.restart` ≠ `"None"`; any
suffixed multi-domain file (phase 3).

### 2.6 Checks

1. **Accounting.** Every key of the resolved yelmox configuration appears in the
   report as mapped (destination group.key), dropped (reason) or refused. An
   unaccounted key is an error.
2. **Value round trip.** The converted file, read back with `Namelist`, gives
   the resolved yelmox value for every mapped key (type and value).
3. **Validation.** `kryos.x --check-par <file>`: reads the par file, runs
   `nml_validate` on every instance and backend group against its defaults
   file, runs the factory's backend constraints (decision 12) and the clock
   checks (every `dt`, output and restart time a whole number of ticks), and
   exits 0 or 1 without reading data. A flag on the one driver (decision 10),
   not a second program.
4. All eleven yelmox bases go through the converter in `make check-tools`:
   `yelmox_Greenland` and `yelmox_Antarctica` (with `sim.relax=False
   sim.opt=False`) convert and validate; the others are refused with the right
   phase. This needs no data and runs locally.

## 3. Comparison tool

### 3.1 Usage

```bash
julia --project=tools tools/bin/kryos-compare.jl REF RUN \
    --map tools/gate/map_phase1.toml --deviations tools/gate/deviations.toml \
    --gate grl [--restarts] [--json out.json] [--plot VAR ...]
```

`REF` and `RUN` are run directories (yelmox or kryos; the map says which file
pattern each side uses; `--map identity` compares two runs of the same model).
Exit 0 pass, 1 fail, 2 error (missing file or record, unreadable map).
`tests/diff_nc.py` in yelmox is not suitable: it drops values outside (0, 1e4)
and works with tolerances.

### 3.2 Mapping table

`tools/gate/map_phase1.toml`, one row per file or variable, `{ice}` etc. the
instance names of the gate. File and variable names on the kryos side are
*(adapters)*; the table is the contract and changes with them. The gate needs
the Yelmo backend to keep its own output (Yelmo's writer, the same
`YELMO_VARS_2D` and `YELMO_VARS_2D_SM` sets, the same schedules): this file
carries the full Yelmo state at `==` and is the most sensitive single check.

| yelmox | kryos | Class |
|---|---|---|
| `yelmo.nc` (all) | `{ice}_yelmo.nc` (backend file) | E |
| `yelmo_sm.nc` (all) | `{ice}_yelmo_sm.nc` | E |
| `yelmo_ts.nc` (Yelmo regions) | `{ice}_yelmo_ts.nc` | E |
| `yelmo_ts.nc:coef_<name>` | coefficient owner's `_ts` *(core)* | E |
| `yelmo.nc:H_ice`, `z_srf`, `dHidt` | `{ice}.nc:H_ice`, `z_srf`, `dHidt` (exports, m, m yr-1) | C: `Float64(x)` |
| `yelmo.nc:smb`, `bmb`, `cmb`, `mb_net` (m ice yr-1) | `{ice}.nc:smb`, ... (kg m-2 yr-1) | C: `Float64(x) * 910.0` (`rho_ice`, `[Earth]`) |
| `htopo.nc:z_bed`, `H_ice`, `z_srf`, `f_grnd`, `z_sl` | hub file *(core)* | E (sp hub, Q2) |
| `isos.nc:z_bed`, `z_ss` | geo backend file `{geo}_fastisostasy.nc` | E |
| `isos.nc:w`, `we` | `{geo}.nc:dz_bed` | C: `Float64(w + we)` (sum in FastIsostasy's kind) |
| `isos_ts.nc:bsl` | `{geo}_ts.nc:bsl` | E |
| `isos_ts.nc:z_bed_mean`, `w_mean`, ... | — | A (framework diagnostics, shed, §11) |
| `chion.nc:smb`, `melt`, `runoff`, `refrz` (mm w.e. yr-1) | `{surface}.nc:smb`, ... (kg m-2 yr-1) | C: `Float64(x)` (factor 1) |
| `chion.nc:tsrf`, `alb_s` | `{surface}.nc:T_srf`, `albedo` | C: `Float64(x)` |
| `chion.nc:pr`, `sf`, `t2m`, `S`, `PDDs`, `H_snow`, `mask` | chion backend file, if kept *(adapters)*, else A | E / A |
| `mshlf.nc:bmb_shlf` (m ice yr-1) | `{shelf}.nc:bmb_shlf` (kg m-2 yr-1) | C: `Float64(x) * 910.0` |
| `mshlf.nc:T_shlf`, `tf_shlf` | `{shelf}.nc:T_shlf`, `tf_shlf` | C: `Float64(x)` |
| `snap.nc` (`tas`, `pr` mm d-1, ocean `to`, `so`, `z_srf`, terms) | `{climate}.nc` | E for K, psu, m; C: `Float64(x) * 360.0` for `pr` (calendar days, *(adapters)*) |
| restart `yelmo_restart.nc` | `restart-<t>/{ice}_restart.nc` | E |
| restart `isos_restart.nc` | `restart-<t>/{geo}_restart.nc` | E |
| restart `marine_shelf.nc` | `restart-<t>/{shelf}_restart.nc` | E |
| restart `chion_restart.nc` | `restart-<t>/{surface}_restart.nc` | E |
| restart `bsl_restart.nc`, `coef_<name>_restart.nc` | global part of geo, coefficient owner *(core)* | E |
| — | manifest, scheduler state, exchange buffers | N |

Every variable of every yelmox file must match a row (A included): a new
yelmox variable without a row is a failure, so the table stays complete.
Variable names inside backend files are unchanged (same writer).

### 3.3 Classes

| Class | Pass when |
|---|---|
| E exact | same type and shape; equal bit for bit (all NaN payloads equal; `-0.0` vs `+0.0` reported as "signed zero", a failure) |
| C conversion-exact | `run` equals, bit for bit, the declared expression applied to `ref` in the declared precision (`Float64(x_sp) * 910.0`), and the `units` attribute equals the declared canonical unit |
| U ulp | `|run − f(ref)| ≤ n·ulp` in the run's type, `n` declared (default 1); for a derived value the two models cannot compute in the same order |
| T tolerance | `|d| ≤ atol + rtol·|ref|`; only through a registered deviation |
| A absent | yelmox variable with no counterpart; listed in the report, never compared |
| N new | kryos-only; listed; covered by `K.chunk` vs `K.cont` |

C is exact because every conversion in phase 1 multiplies an `sp` value by a
constant with a short mantissa (910, 360, 86400, 1000): the product is exact in
`dp`, and dividing by the same constant gives back the `sp` value, so the
adapter's inverse conversion restores the bits Yelmo got in yelmox. The
exception is a conversion and its inverse separated by a remap (§5, row 5).

The time axes are compared first: every record time of a `ref` file must exist
in the `run` file, equal bit for bit. Static variables (lat, lon, area, masks)
are compared like the others.

### 3.4 Report

Per variable: class, status, cells differing, max abs, max rel, max ulp, first
differing record (time and index), location (i, j) of the maximum. Per gate:
the **first divergence**, the earliest time at which any compared variable
differs, with all variables that differ at that time, ordered by the coupling
sequence (geo, ice, hub, climate, surface, shelf): the root-cause hint, since
after the first difference everything diverges. Markdown and JSON. `--plot
VAR` draws, with CairoMakie, yelmox, kryos and the difference at the first
differing record and at the last, and the time series of max |difference|.

### 3.5 Deliberate non-identity

`tools/gate/deviations.toml`, one entry per accepted difference:

```toml
[[deviation]]
id        = "chion-diag-mean"
rows      = ["chion.nc:PDDs"]
class     = "U"            # or "T" with atol, rtol; or "skip"
scope     = "local"        # local | rebaseline
since     = "<kryos sha>"
changelog = "CHANGELOG.md#..."   # must exist
approved  = "Alex, <date>"
```

- **local**: output-only differences that do not feed back (a diagnostic, a
  label). Only the listed rows are relaxed; every other row must still pass,
  otherwise the deviation is not local and the gate fails.
- **rebaseline**: an answer-changing difference (a bug fix, the `dp` hub).
  The yelmox gate closes at the last identical kryos commit (recorded in
  `docs/gate.md`); from then on kryos is gated against its own runs at the
  `since` commit (`--map identity`, bit for bit), and the yelmox comparison
  becomes a drift report of integrated quantities (`V_ice`, `A_ice`, `V_af`,
  `bsl`) without pass/fail. Phase 1 aims for no rebaseline before sign-off.

The tool lists deviations that matched nothing; they are removed.

## 4. Automation

### 4.1 Stack

`tools/gate/stack.toml`: commit of yelmox, kryos and each library (fesm-utils,
yelmo with FastHydrology/elsa/tracer, FastIsostasy, chion, insol, rembo1), the
configme and runme versions, compiler and netCDF modules, flags,
`OMP_NUM_THREADS = 16`, partition. One stack, two installs sharing the
dependency checkouts and archives:

```bash
configme install yelmox:v2.4+fesm-utils:phase0-kryos -m dkrz_levante -c ifx \
    --dir $G/stacks/$S/yelmox
configme install kryos -m dkrz_levante -c ifx --dir $G/stacks/$S/kryos \
    --link fesm-utils=$G/stacks/$S/yelmox/fesm-utils --link yelmo=... \
    --link FastIsostasy=... --link chion=... --link insol=...
```

then each checkout is set to the exact commit of `stack.toml` (detached) and
reconfigured; `gate.sh` refuses a dirty checkout or one off its pin. configme's
kryos orchestrator needs `data_link = ["ice_data", "fesmdata",
"isostasy_data"]` as yelmox has, and kryos a `.runme/info.json` (exe alias
`kryos`, `par_defaults`, `par_library`, `links`).

**Flags.** Both models and every library are built for the gate with
`DFLAGS_NODEBUG="-O2 -fp-model precise -march=znver3 -traceback"` (on the make
command line, so it reaches the recursive dependency builds): `-Ofast` allows
FMA contraction and reassociation that depend on the shape of the source, and
kryos rewrites all the glue code (conversions, anomaly sums, wrappers) that
yelmox had in its own files. Library code compiled once is the same in both
models either way; the glue is not. The production `-Ofast` pair is compared
once per phase, as information (Q3).

### 4.2 `tools/gate/gate.sh`

`make gate` in kryos calls it. Subcommands (all is the default):

| Step | Does |
|---|---|
| `stack` | installs or verifies the stack, builds `yelmox.x` and `kryos.x`, writes `stack.lock.json` |
| `maps` | seeds `$G/maps/<stack>` from yelmox's `maps/` (grid descriptions and cached files); after R0 has run, makes it read-only |
| `ref <gate>` | if the R1 folder is missing: stages `cont`, `cont2` (new stack), `chunk` with `runme -s` and submits with `sbatch --parsable` (the `run_chunks.sh` pattern); links the frozen maps |
| `run <gate>` | converts the staged R1 par file (§2), validates it (`--check-par`), stages and submits the kryos `cont` and `chunk` runs; checks that the input-file hashes of both run directories agree |
| `compare` | one job on the shared partition with `--dependency=afterok:` on all run jobs: `kryos-compare.jl` for R1 vs R0 (new stack), `cont2` vs `cont`, `chunk` vs `cont` (both models), K vs R1 (outputs and `--restarts`) |
| `record` | prints the `docs/gate.md` row |

Model runs use runme as yelmox does (`-q compute -w 00:30:00 --omp 16`); the
compare job uses a template `tools/gate/compare.sbatch` (account `ba1442`,
partition `shared`, 1 task, 4 cpus, 16 GB, 30 min, Julia from juliaup in
`$HOME`). The whole gate is submitted at once and needs no waiting.

### 4.3 In the repository

- `make check`: the Fortran unit tests (`tests/test_*.f90`), including
  `test_conversions.f90`: each adapter conversion against its literal
  expression from the map table, on random `sp` values with zeros, signed
  zeros, subnormals, extremes and missing values; and the remap rows of §5
  (the `sp` remap path against `dp` remap + one rounding).
- `make check-tools`: `Pkg.test()` of KryosTools (comparison on synthetic file
  pairs, one per class and failure mode) and the converter on all yelmox bases
  (needs `YELMOXROOT`).
- `test_nml_compat` (in fesm-utils, `make check`): reads every key of a list of
  yelmox files as a string vector with `nml` and prints it; `bitident.sh`-style
  comparison against `v1.4`. Run on the staged gate par files and the yelmox
  library and defaults files. Closes the unproven row of §1.3 for all yelmox
  configurations, not only the two gates.

### 4.4 Records

`docs/gate.md`: one row per official gate (Levante only): date, gate, kryos,
yelmox, stack id (and the commits it resolves to), flags, length, outcome,
first divergence (if any), deviations used, report path on Levante. The JSON
reports stay on Levante under `$G/reports`. A phase is signed off by a row with
both gates green.

### 4.5 Local (this Mac)

What runs now, without data: `make check`, `make check-tools`, the converter
and `--check-par` on every yelmox base, `test_nml_compat`, and the comparison
tool on R0 vs R1 copies once they exist. NorthTest (NHT-16KM/32KM: topography,
regions, basins) is enough for hub and regions tests, not for a coupled run:
phase 1 has no `data` climate or surface.

A local short gate (100 yr, gfortran `-O2 -ffp-contract=off`, 4 threads,
`GATE_MACHINE=macbook`) needs the GRL-16KM, ANT-16KM and ANT-32KM products
synced from `/work/ba1442/fesmdata` to `~/data/fesmdata` (plus a `fesmdata`
link in yelmox and the Antarctica `ice_data` files); the grids are small
(106×181 to 381×381). Local gates speed up debugging; only Levante rows count.

## 5. Expected sources of non-identity

| # | Source | Phase 1 handling |
|---|---|---|
| 1 | Newly generated maps (prune criterion, spherical clip, new map pairs kryos might request, e.g. for `z_bed_ref` on component grids) | Both models run on the frozen yelmox map cache, read-only: an attempt to generate a map stops the run and fails the gate. kryos uses exactly yelmox's (src, dst, method) pairs in phase 1. |
| 2 | Identity remaps | yelmox short-circuits; kryos goes through `coupler remap`, which copies when `src == dst` (phase 0). The cached `scrip-*_GRL-16KM_GRL-16KM.nc` files must never be applied. All grids are identical in `grl`; in `ant` hub/shelf are ANT-16KM, ice/geo/climate/surface ANT-32KM. |
| 3 | `sp` ↔ `dp` exchange | Upcast is exact. The adapter rounds to `sp` at the point yelmox did: a remapped field is rounded once (`real(v, sp)`, as `map_field_grid_grid_sp`), then combined in `sp` in yelmox's order (`z_bed = z_bed_ref + dz_bed`, `smb * conv_we_ie * 1e-3`). Rounding the sum instead of the remapped term is a double rounding and differs. `test_conversions`. |
| 4 | `kg m-2 yr-1` vs m ice yr-1, mm w.e. yr-1, mm d-1 | Exact when the value is `sp` and the factor short (§3.3): smb ×1, bmb ×910 and ÷910, pr ×360 and ÷360, then ÷86400 as yelmox. `rho_ice`, `rho_w`, `sec_year`, `sec_day` from `input/yelmo_phys_const.nml [Earth]`, passed to each backend as yelmox passes them. `test_conversions`. |
| 5 | A conversion and its inverse separated by a remap | Only `ant` `bmb_shlf` (shelf ANT-16KM → ice ANT-32KM, `con`): `con(x·ρ)/ρ` differs from `con(x)` by ~1 ulp in `dp`, so the `sp` result can flip 1 ulp when the exact value is within that of a rounding midpoint (order 1e-9 per cell and update). Quantified before the first Antarctica gate by a test on the recorded lossless `mshlf.nc` of R1: count flips over all records. If any occur, it is a decision for Alex (accept as a rebaseline, or remap that field in m ice yr-1 inside the shelf→ice wiring). |
| 6 | Hub precision | `htopo` is `sp` today; design §5 moves it to `dp` (§12 row 8), which changes every consumer of hub geometry. Phase 1 keeps the `sp` hub until the gate is green; the `dp` hub is then a rebaseline (Q2). |
| 7 | ncio read flush (`|x| < 1e-30` → 0) on restart | The same backend restart files are read by the same ncio in both models. New kryos state (exchange buffers, scheduler) is read through ncio too and only feeds values that yelmox also restored through ncio (Yelmo `bnd`). Covered by `chunk` vs `cont`. |
| 8 | Lossless restarts | Backend restart writers are unchanged (Yelmo forces lossless on its table); `fesm_fields_restart_write` writes `dp` without rounding. E on the bundles at 500 and 1000. |
| 9 | Output precision rounding | Gate runs are lossless. Backend files keep their writers and per-field precision; instance export files declare theirs in the variable table. Once per phase a non-lossless `grl` pair is compared: E rows must still pass. |
| 10 | `bsl` reduction | One domain: `Σ_d` over one term is that term; the expression must be yelmox's (`bsl_update` order). E on `bsl`. |
| 11 | Clock | Times are `n·tick` in kryos and accumulated `dtt` in yelmox; with integer `dtt` both are exact integers in `dp`. Gate `dt` values stay integer years. Time axes compared first. |
| 12 | Defaults drift | Removed by the fully expanded par file (§2.3). |
| 13 | `nml` parsing changes | `test_nml_compat` and R1 == R0. |
| 14 | Compiler | Gate flags without fast-math; same compiler, modules and archives for both (§4.1). |
| 15 | OpenMP | Same thread count, partition and node type for both; `cont2` checks run-to-run determinism once per stack. |

## 6. Order of work

1. KryosTools `Namelist` + `Compare` with tests (local, now).
2. Levante stacks `v24` and `S`: R0, R1, `cont2`, `chunk`; R1 == R0 is the
   phase 0 sign-off and needs no kryos code. `test_nml_compat`; the
   `bmb_shlf` flip count (§5 row 5).
3. Converter once `core-design.md` fixes the schema; `--check-par`.
4. kryos gates as components land (ladder §1.5), then both full gates, the
   `docs/gate.md` row, sign-off; then the `dp` hub rebaseline, 15-kyr run.

## 7. Open questions

1. **Antarctica gate**: phase 1 compares `ant` with `sim.relax=False
   sim.opt=False` (proposed; the plugins are phase 2), or the `relax` and
   `opt_cb_ref` plugins move into phase 1 so the shipped base is the target?
2. **Hub precision**: keep `htopo` in `sp` through the phase 1 gate and move it
   to `dp` afterwards as a registered rebaseline (proposed), or `dp` now and
   give up bit identity for every hub consumer?
3. **Gate flags**: gate builds with `-O2 -fp-model precise` (proposed) and the
   production `-Ofast` pair as information only, or must the production flags
   be bit-identical too?
4. **Access and data**: can the gate be run from this environment (`ssh
   levante` non-interactively, `sbatch` as a271116), or does Alex start
   `gate.sh`? And may the GRL-16KM/ANT-16KM/ANT-32KM FesmData products be
   synced to `~/data/fesmdata` for local short gates?
5. **Namelist reader**: move `Namelist.jl`/`Resolve.jl` from yelmox's forcing
   builder to FesmUtils.jl, used by both (proposed), or keep a copy in kryos
   tools?
