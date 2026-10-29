# Running kryos

## Build

See [build.md](build.md): link or clone the dependencies, `configme config`,
then `make deps-static` and `make kryos` (the driver `libkryos/bin/kryos.x`).
`make check` runs the tests.

## Configure

A run reads one par file (namelist). The shipped ones are converted from the
yelmox configurations (yelmox `topo-prep`):

| Par file | Source | Notes |
|---|---|---|
| `par/kryos_Greenland.nml` | `yelmox/yelmox_Greenland.nml` | 15 kyr, yelmox output schedules |
| `par/kryos_Antarctica.nml` | `yelmox/yelmox_Antarctica.nml` | `-p sim.relax=False sim.opt=False` (relaxation and optimisation are phase 2) |

Each has its conversion report next to it (`par/kryos_<domain>_conversion.md`:
every yelmox key mapped, default-equal or dropped, with the reason). Do not
edit a converted file: change the yelmox source and convert again
(`tools/bin/yelmox2kryos.jl`, Julia project `tools/`).

Groups and defaults:

- kryos groups (`[kryos]`, the clock, schedules, `[globe]`, domains,
  instances) are validated against the templates in `input/kryos_defaults.nml`
  (`[instance]` plus `[<component>_<backend>]` for an instance); a key left out
  takes the template's value.
- Backend library groups (`[yelmo]`, `[ytopo]`, ..., `[chion]`) take their
  defaults from the library's defaults file in `input/` (`yelmo_defaults.nml`,
  `chion_defaults.nml`); data field groups from `datasrc_defaults.nml`;
  snapshot, term and coefficient groups from `snapshots_defaults.nml` and
  `forcing_defaults.nml`. Groups whose library has no defaults file (`[topo]`,
  `[regions]`, `[isos]`, `[barysealevel]`, `[marine_shelf]`) have no
  template: set every key the library reads. `[vilma]` (`geo = vilma2`) takes
  its defaults from `input/vilma_defaults.nml`: the group sets overrides only,
  and its name is fixed (VILMA2 looks the defaults up under it).
- Named library groups (`snapshot_*`, `term_*`, `coef_*`, data field groups)
  are looked up in the par file, then in the files of `[kryos] library`
  (namelist files of such groups; kryos ships none, the converted par files
  hold their groups).

So a par file may list only what differs from these defaults and libraries.
The converted files are fully expanded instead: every key of every group is
set, `library = ""`, and the file needs no defaults (it stays valid when a
default changes; the header records the sha256 of the input files it was made
against).

## input/

The run directory's `input/` is the repository's `input/` (runme links it).
It holds kryos' own files and synchronised copies of the files the libraries
read.

kryos' own: `kryos_defaults.nml`, `forcing_defaults.nml`,
`snapshots_defaults.nml`, `kryos-variables-<component>.md` and the tables of the
backend streams (`kryos-variables-<component>-<backend>.md`,
`kryos-variables-hub.md`; precisions of the output, see Output), `phys_const.nml`
(the `[kryos] phys_const_file` default, read by `phys_const_load`; the run's
one set of constants, passed to every backend), and the series of the
configurations: `sl_deglaciation.dat` (Greenland `[barysealevel] sl_path`),
`alpha_orbital_r2020.dat` (Greenland `[coef_alpha]`),
`alpha_combined_125kyr_interp.dat` (the yelmox default glacial index).

Synchronised copies, listed in `input/sync.txt` with their source:

| Copies | Source |
|---|---|
| `yelmo_defaults.nml`, `yelmo_phys_const.nml`, `yelmo-variables-*.md` | `$(YELMOROOT)/input/` |
| `elsa_defaults.nml`, `tracer_defaults.nml` | `$(YELMOROOT)/elsa/input/`, `$(YELMOROOT)/tracer/input/` (read when Yelmo's tracers are on) |
| `chion_defaults.nml` | `$(CHIONROOT)/input/` |
| `datasrc_defaults.nml` | `$(FESMUTILSROOT)/par/` |
| `vilma_defaults.nml`, `vilma-variables.md` | `$(VILMAROOT)/input/` (the `[vilma]` defaults and the variable table `vilma_io` reads, `geo = vilma2`) |
| `LA2004.INSOLN.5Ma.txt`, `LA2004.INSOLP.1Ma.txt` | `$(INSOLROOT)/input/` (`[orbit] fldr`) |

Yelmo reads `sec_year` from `input/yelmo_phys_const.nml`
(`ybound_define_physical_constants`) even though kryos passes it the run's
constants, so the file must be present (its constants are currently identical to
`phys_const.nml`).

`make sync-input` copies every source over its copy; `make check` first runs
`make input-check`, which fails naming each copy that differs from its source
(or a missing source). Update a library, then `make sync-input` and commit the
copies.

## Run

With runme (`.runme/info.json`: executable `kryos` = `libkryos/bin/kryos.x`,
par file `par/kryos_Greenland.nml`, defaults `input/kryos_defaults.nml`; links
`input`, `ice_data`, `fesmdata`, `isostasy_data`, `maps` from the repository
root into the run directory). Once per checkout: `runme config init`, then set
`hpc` and `account` in `.runme/config.toml`.

```bash
# local, 4 threads, 1000 years
runme -r -e kryos -n par/kryos_Greenland.nml -o output/grl-test --omp 4 \
      -p clock.time_end=1000

# submitted (Levante; queue aliases: runme queues)
runme -rs -q compute -w 24:00:00 -e kryos -n par/kryos_Antarctica.nml \
      -o output/ant --omp 16

# ensemble over a Yelmo parameter
runme -rs -q compute -w 24:00:00 -e kryos -n par/kryos_Greenland.nml \
      -o output/grl-ens --omp 16 -p ydyn.beta_q=0.2,0.5,1.0
```

`-p group.key=value` overrides a key of the par file (`bsl` is an alias of
`barysealevel`); the clock group is the one `[kryos] clock` names (`clock`).
Without runme: from a run directory with the links above,
`./kryos.x <par file>`.

## Output

Files go to `[kryos] out_dir` (`./`), one stream per instance and schedule
group (`out_2D`, `out_1D`, ... name schedule groups such as `[out_2D]`):

| File | Content |
|---|---|
| `<instance>.nc` | 2D/3D exports of the instance (`out_2D`) |
| `<instance>_ts.nc` | 0-D exports (`out_1D`) |
| `<instance>_<stream>.nc` | a backend's own streams, e.g. `ice_grl_yelmo.nc`, `ice_grl_yelmo_sm.nc`, `ice_grl_yelmo_ts.nc` |
| `<domain>_hub.nc` | the domain's hub topography (`[<domain>] out_2D`) |
| `timing.txt` | timing table (on the first 1D schedule) |

Instance names are unique in a run, so the files need no domain prefix.

### Precision

Output is rounded per variable (ncio bit rounding, `float`), as yelmox does; the
restart bundle is never rounded. The `prec` column of the variable tables in
`input/` is the only source of rounding: `-` lossless, `1e-3` absolute (in the
units of the variable), `1e-4r` relative, `1e-3,1e-4r` both. A variable a table
does not list is written lossless; a row without a unit, a bad `prec` entry, or a
unit that differs from the data's stops the run. The tables:

| Table | Read by | Content |
|---|---|---|
| `kryos-variables-<component>.md` | the default writer (`component_write_output`) | exports of an instance of that component kind |
| `kryos-variables-<component>-<backend>.md`, `kryos-variables-hub.md` | the backend's own writer | `geo-fastisostasy`, `climate-snapshots`, `surface-chion`, `shelf-marine_shelf`: the variables of that backend's 2D stream; `hub`: `z_bed`, `H_ice`, `z_srf`, `f_grnd`, `z_sl` (passed to `htopo_write_step`) and `f_ice` of the hub file |
| `yelmo-variables-*.md` | Yelmo (synchronised copies) | the Yelmo streams |

`[kryos] write_lossless = True` (ncio `nc_set_lossless`) ignores every precision:
the gates set it, so their files are never rounded. The default tolerances are
yelmox's: m 1e-3, m yr-1 1e-3 and 1e-4 relative or 1e-4 absolute (rates), K 1e-3,
fractions 1e-4, mass balances 1e-1 kg m-2 yr-1, precipitation, insolation and
degree days 1e-4 relative.

Every file kryos writes, and its precision (audited against the yelmox writers
of `libs/kryos_output.f90`, `yelmox_climate.f90`, `snapshots.f90`,
`climate_mods.f90`: the values of the backend streams are identical):

| File | Written by | Precision |
|---|---|---|
| `<instance>.nc` (ice, geothermal, sediments) | default writer | table of the component kind |
| `<instance>_ts.nc` (0-D exports), `geo_glob_ts.nc` | default writer | lossless (`-` in every 0-D row) |
| `forcing_ts.nc` | tsgen adapter (`nc_write`) | lossless |
| `<instance>.nc` of a regional geo, fastisostasy (`z_bed`, `z_ss`, `w`, `we`) | FastIsostasy adapter | `geo-fastisostasy`: 1e-3 |
| `<instance>_ts.nc` of a regional geo, fastisostasy | FastIsostasy adapter | lossless |
| `<instance>.nc` of a climate, snapshots | snapshots adapter, snapshots and transform modules | `climate-snapshots`: temperature 1e-3, precipitation and salinity 1e-4 relative; static transform fields (`taper_w`, `bias_*_ann`) and `term_<name>` lossless |
| `<instance>.nc` of a surface, chion | chion adapter | `surface-chion`: smb, melt, runoff, refreezing, snowpack 1e-1, T 1e-3, `pr`, `sf`, `S`, `PDDs` 1e-4 relative, albedo 1e-4, `mask` lossless |
| `<instance>.nc` of a shelf, marine_shelf | marine_shelf adapter | `shelf-marine_shelf`: `bmb_shlf` 1e-4, `T_shlf`, `tf_shlf` 1e-3 |
| `<domain>_hub.nc` | fesm-utils `htopo_write_step` (`z_bed`, `H_ice`, `z_srf`, `f_grnd`, `z_sl`), kryos (`f_ice`; `mask` lossless) | `hub`: `z_bed`, `H_ice`, `z_srf`, `z_sl` 1e-3, `f_grnd` and `f_ice` 1e-4 |
| `<instance>_yelmo.nc`, `<instance>_yelmo_sm.nc` | Yelmo | `yelmo-variables-*.md` |
| `<instance>_yelmo_ts*.nc` (global and region series) | Yelmo | lossless |
| `timing.txt` | timing table | text |
| `restart-<time>/` (every file) | `fesm_fields` restart, backends | lossless (`dp`; Yelmo restart forced lossless; chion, marine_shelf, FastIsostasy restarts written without precision) |

Rounded and lossless runs have identical restart bundles; a rounded file is
within the tolerance of its row of the lossless one (the `sci-nl` gate checks it
on the science configuration). The output files of the 80 yr GRL-16KM run of
`test_run_sequence` (no Yelmo streams) are 37% smaller rounded (the instance
files -24 to -57%; 22.3 MB lossless, 14.1 MB rounded; with the Yelmo streams
`yelmo` and `yelmo_sm` -35%: 32.4 MB, 21.0 MB); the bundles are the same size.

## Restart

When the schedule named by `[kryos] restart_out` is due, and always at
`time_end`, the run writes a bundle `<out_dir>/restart-<time>/` (`<time>` in
years, shortest exact decimal, e.g. `restart-15000`):

```
restart-15000/
  clock_restart.nc                    clock and the cadence of every instance that steps
  forcing_imports.nc  forcing_exports.nc  forcing_<coef>_restart.nc
  geo_glob_imports.nc  geo_glob_exports.nc  geo_glob_restart.nc
  Greenland/
    <inst>_imports.nc  <inst>_exports.nc  every instance of the domain
    <inst>_restart.nc                     stateful backends (Yelmo: + elsa/tracer sidecars when active)
  manifest.nml                        [bundle] + one group per instance; written last
```

`manifest.nml` holds `[bundle]` (`time`, `ticks`, `ticks_per_year`,
`kryos_git`, `fesm_utils_git`, `domains`) and one group per instance (`scope`,
`domain`, `component`, `backend`, `grid`, `files`, `git`). A bundle without
it is incomplete.

Restart from a bundle with `[kryos] restart` (a folder relative to the run
directory, or absolute) and the clock continued, e.g. the second chunk of a
run:

```bash
runme -rs -q compute -w 24:00:00 -e kryos -n par/kryos_Greenland.nml \
      -o output/grl-2 --omp 16 \
      -p clock.time_init=15000 clock.time_end=30000 \
         kryos.restart=../grl/restart-15000
```

The restore stops, naming the instance, if the run's instances differ from
the manifest's (scope, domain, component, backend, grid) or a listed file is
missing; a changed `dt` is allowed. Restart times and `time_end` must be on
every instance's cadence grid (checked at init).

### Forcing-increment bundles

A tsgen coefficient with `[coef_<name>] restart_every_df = <df>` (> 0; 0 = off, the
default) also writes a bundle `<out_dir>/restart-<n>/` (`n` = 1, 2, ... in the order
written; the same content as `restart-<time>`, its `manifest.nml` has the time) each
time it has moved by more than `df` since its last bundle (either direction; the
first mark is the initial value). It is written after the step's output, on a step
whose time is on every cadence grid (a trigger between two such steps fires on the
next one). The count is stored in the bundles: a run restored from `restart-3` numbers on from 4. A bundle
is never overwritten: a run that would write into an existing `restart-<n>` (restored
from an older bundle than the last in its `out_dir`) stops; give it another `out_dir`.
A derived coefficient cannot carry `restart_every_df`.
