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
  template: set every key the library reads.
- Named library groups (`snapshot_*`, `term_*`, `coef_*`, data field groups)
  are looked up in the par file, then in the files of `[kryos] library`
  (e.g. `"input/clim_db.nml" "input/forcing_db.nml"`).

So a par file may list only what differs from these defaults and libraries.
The converted files are fully expanded instead: every key of every group is
set, `library = ""`, and the file needs no defaults (it stays valid when a
default changes; the header records the sha256 of the input files it was made
against).

## input/

The run directory's `input/` is the repository's `input/` (runme links it).
It holds kryos' own files (`kryos_defaults.nml`, `forcing_defaults.nml`,
`snapshots_defaults.nml`, `phys_const.nml`, `kryos-variables-<component>.md`)
and synchronised copies of the files the libraries and the configurations
read, listed in `input/sync.txt` with their source:

| Copies | Source |
|---|---|
| `yelmo_defaults.nml`, `yelmo_phys_const.nml`, `yelmo-variables-*.md` | `$(YELMOROOT)/input/` |
| `elsa_defaults.nml`, `tracer_defaults.nml` | `$(YELMOROOT)/elsa/input/`, `$(YELMOROOT)/tracer/input/` |
| `phys_const.nml` | `$(YELMOROOT)/input/yelmo_phys_const.nml` |
| `chion_defaults.nml` | `$(CHIONROOT)/input/` |
| `datasrc_defaults.nml` | `$(FESMUTILSROOT)/par/` |
| `LA2004.INSOLN.5Ma.txt`, `LA2004.INSOLP.1Ma.txt` | `$(INSOLROOT)/input/` |
| `forcing_db.nml`, `clim_db.nml`, `sl_deglaciation.dat`, `alpha_orbital_r2020.dat`, `alpha_combined_125kyr_interp.dat` | `$(YELMOXROOT)/input/` (no library source) |

`YELMOXROOT` defaults to a `yelmox` link at the repository root, like the
dependencies. `make sync-input` copies every source over its copy;
`make check` first runs `make input-check`, which fails naming each copy that
differs from its source (or a missing source). Update a library, then
`make sync-input` and commit the copies.

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

Instance names are unique in a run, so the files need no domain prefix. The
variable tables are `input/kryos-variables-<component>.md`;
`[kryos] write_lossless = True` writes full precision.

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
