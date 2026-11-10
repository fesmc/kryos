---
title: "Getting started"
---

From an empty folder on DKRZ Levante to a run of the plastic ice sheet (ice backend `plastic`, the Fortran
library [PlasticIceSheet](https://github.com/fesmc/PlasticIceSheet)) coupled to the global VILMA2 Earth (geo
backend `vilma2`). Everything here runs on Levante: builds, tests, runs and analysis. The data paths are those of
project `ba1442`.

- **Program:** `kryos.x`, built with `make kryos plastic=1 vilma=1 yelmo=0`.
- **Config:** a par file of [par/](https://github.com/fesmc/kryos/tree/main/par), e.g.
  `par/kryos_North_plastic_vilma2_b1.nml` (see [Par files](reference-par-files.md)).
- **Data:** FesmData products under `/work/ba1442/fesmdata` (read-only), kryos-made inputs under
  `/work/ba1442/kryos_data`.

## Environment

Load the modules in every shell that builds (e.g. in `~/.bashrc`; they are not Levante's defaults):

```bash
export MODULEPATH=/sw/spack-workplace/spack-0.20.1/modules:$MODULEPATH
module load intel-oneapi-compilers/2023.2.1-gcc-11.2.0 \
    netcdf-c/4.9.3-openmpi-5.0.8-intel-2021.5.0 netcdf-fortran/4.6.2-openmpi-5.0.8-intel-2021.5.0 \
    python3/2023.01-gcc-11.2.0 git/2.43.7-gcc-11.2.0
```

configme's machine file `dkrz_levante` hardcodes the matching netCDF prefixes (checked 2026-11-10). The
production compiler is Intel `ifx` (configme machine `dkrz_levante`, compiler `ifx`).

configme writes the Makefiles; install it once, from its `fortdeps` branch (its `main` has no kryos package):

```bash
pip install git+https://github.com/fesmc/configme@fortdeps
export PATH="${PATH}:${HOME}/.local/bin"      # if `configme` is not found
```

## The code

The dependencies are checkouts, or symlinks to checkouts, at the kryos root (the slots are gitignored). The pins are in
[.configme/manifest.toml](https://github.com/fesmc/kryos/blob/main/.configme/manifest.toml):

| Checkout | Repository | Ref | Needed for |
|---|---|---|---|
| `kryos` | [fesmc/kryos](https://github.com/fesmc/kryos) | `main` | everything |
| `fesm-utils` | [fesmc/fesm-utils](https://github.com/fesmc/fesm-utils) | `phase0-kryos` | always (ncio, grids, coupler, htopo; LIS, FFTW, SHTns) |
| `vilma2` | [fesmc/vilma2](https://github.com/fesmc/vilma2) | `main` | `vilma=1`; and `make check` in any build (two `input/` copies come from it) |
| `PlasticIceSheet` | [fesmc/PlasticIceSheet](https://github.com/fesmc/PlasticIceSheet) | `main` | `plastic=1`; not in the manifest: clone it by hand |
| `yelmo` | [fesmc/yelmo](https://github.com/fesmc/yelmo) | `stable` | `yelmo=1`; and `make check` (the Yelmo copies in `input/`) |
| `FastIsostasy` | [palma-ice/FastIsostasy](https://github.com/palma-ice/FastIsostasy) | `stable` | `fastisostasy=1` (default on) |
| `chion` | [fesmc/chion](https://github.com/fesmc/chion) | `stable` | `chion=1` (default on) |
| `insol` | [fesmc/insol](https://github.com/fesmc/insol) | `stable` | `insol=1` (default on; chion needs it) |

FastHydrology, elsa and tracer live inside the yelmo checkout.

Clone kryos and its dependencies at the refs of the table (`https://github.com/...` without a GitHub SSH key):

```bash
cd ~/models
git clone git@github.com:fesmc/kryos.git
cd kryos
git clone -b phase0-kryos git@github.com:fesmc/fesm-utils.git
git clone git@github.com:fesmc/vilma2.git
git clone git@github.com:fesmc/PlasticIceSheet.git
git clone -b stable git@github.com:fesmc/yelmo.git
git clone -b stable git@github.com:palma-ice/FastIsostasy.git
git clone -b stable git@github.com:fesmc/chion.git
git clone -b stable git@github.com:fesmc/insol.git
```

For checkouts that already exist, or clones next to kryos, symlink them instead (`ln -s ../fesm-utils fesm-utils`
at the kryos root).

`configme install kryos` covers only part of this: it clones the seed set insol, fesm-utils, yelmo, FastIsostasy,
FastHydrology and chion and builds fesm-utils, insol and chion, but does not clone vilma2 or PlasticIceSheet (the
manifest only pins refs).

## Configure

```bash
cd ~/models/kryos
configme config kryos --only -m dkrz_levante -c ifx
```

This writes the root `Makefile` (the template
[config/Makefile](https://github.com/fesmc/kryos/blob/main/config/Makefile), the compiler and machine fragments and
netCDF) and `config/deps.mk`, the per-object dependency rules generated from the `module`, `use` and `include`
statements under `src/` (configme >= 0.16). `--only` leaves the dependencies' own Makefiles alone. Rerun it, or
`make depend`, after adding a source file or a `use`. The dependency wiring is in
[config/common.mk](https://github.com/fesmc/kryos/blob/main/config/common.mk).

## Build

`make` alone lists the targets and the current switches. `openmp=1` is the default; fesm-utils, VILMA2,
PlasticIceSheet and kryos must all use the same `openmp=`.

### Dependencies (once)

fesm-utils' external libraries, one variant, serially (no `-j`; SHTns links the FFTW built here, so FFTW comes
first):

```bash
cd ~/models/kryos/fesm-utils
make fftw  openmp=1
make lis   openmp=1
make shtns openmp=1
```

Then fesm-utils and the backend libraries the build links, from the kryos root:

```bash
cd ~/models/kryos
make deps-static vilma=1 yelmo=0
```

`deps-static` runs `fesmutils-static` and the `<backend>-static` target of every backend switched on:
here FastIsostasy, chion, insol (on by default) and VILMA2 (`make -C vilma2 vilma-static openmp=1`). With
`plastic=1` it would also run `plastic-static`, but that passes only `config=` to PlasticIceSheet (not `openmp=`
or `lis=`), so build PlasticIceSheet by hand:

```bash
cd ~/models/PlasticIceSheet
L=$HOME/models/kryos/fesm-utils/lis/lis-omp     # the LIS kryos links (openmp=1)
make static config=dkrz_levante_ifx openmp=1 lis=1 LISROOT=$L
make check  config=dkrz_levante_ifx openmp=1 lis=1 LISROOT=$L
```

- `make static` writes `libplastic/include/libplastic.a` and the `.mod` files; `make check` builds and runs the
  library's tests (`test_core`, `test_closures`, `test_surface`, `test_steady`, `test_inversion`).
- Pass the same switches to `check` as to `static`: with `lis=1` the tests link `-llis`.
- `lis=1` (ifx only): Lis for the linear solves of the steady SMB-driven solver. Without `LISROOT` the
  `dkrz_levante_ifx` fragment points at a personal fesm-utils build; pass the LIS of your kryos checkout, the one
  kryos puts on its link line (`LIB_PLASTIC` includes `LIB_LIS`).
- `make clean` between compilers: the `.mod` files are compiler-specific.
- API, options and tests: PlasticIceSheet
  [docs/fortran-port.md](https://github.com/fesmc/PlasticIceSheet/blob/main/docs/fortran-port.md).

`make deps-check plastic=1 vilma=1 yelmo=0` fails, naming the command that builds it, for every archive that is
missing (e.g. SHTns: `make -C fesm-utils shtns openmp=1`).

### kryos

```bash
cd ~/models/kryos
make kryos plastic=1 vilma=1 yelmo=0      # -> libkryos/bin/kryos.x
make check plastic=1 vilma=1 yelmo=0      # the tests
```

- The driver target is `kryos`. It builds `libkryos/include/libkryos.a` and relinks
  `libkryos/bin/kryos.x` on every call. `make libkryos/bin/kryos.x` does **not** relink the driver: the only rule for
  `libkryos/bin/%.x` builds the test programs from `tests/%.f90`.
- `yelmo=0` is enough for the plastic runs: the Yelmo adapter becomes its stub, and the factory stops with "rebuild
  kryos with yelmo=1" if a par file asks for it. The other switches stay at their defaults.
- One binary built with `plastic=1 vilma=1` runs every plastic par file, the fixed-bed one too.
- Objects are stamped with the switches (`libkryos/include/.build-openmp1-...`): changing any switch rebuilds every
  object. Pass the same switches to `make kryos` and `make check`.
- `make check` first runs `make input-check`, which compares the synchronised copies in `input/` with their sources
  in the yelmo, chion, insol, fesm-utils and vilma2 checkouts: those must be present whatever the switches.
- With `plastic=1` and `vilma=1` it runs, among the others, `test_plastic_margin`, `test_ice_plastic` and
  `test_geo_vilma2`. The long tests (`make check-long`) need `yelmo=1`: with `yelmo=0` all three are skipped.

Targets:

| Target | Builds |
|---|---|
| `make kryos` | `libkryos/bin/kryos.x`, the driver ([drivers/kryos.f90](https://github.com/fesmc/kryos/blob/main/drivers/kryos.f90)) |
| `make kryos-static` | `libkryos/include/libkryos.a` and the `.mod`/`.smod` files in `libkryos/include` |
| `make check` | `make input-check`, then every `tests/test_*.f90` but the long ones, linked against libkryos, and runs them; fails if one exits non-zero |
| `make check-long` | the long tests: `test_run_sequence` (GRL-16KM, every backend, restart round trips, about 4.5 min), `test_run_domains` (GRL-16KM + ANT-32KM, about 8 min; needs Julia with the project of `tools/`, else it reports OPEN and passes), `test_gia_sequence` (VILMA2 with GRL-16KM + ANT-32KM at lmax 16, about 6 min; needs `vilma=1`) |
| `make sync-input` | copies the library sources over the synchronised copies in `input/` (`input/sync.txt`) |
| `make input-check` | fails naming every copy in `input/` that differs from its source |
| `make deps-check` | fails with instructions if a dependency archive is missing |
| `make deps-static` | builds fesm-utils and the enabled backend libraries in their own trees |
| `make depend` | regenerates `config/deps.mk` (`configme deps .`) |
| `make clean` / `make cleanall` | kryos objects and executables / also the dependencies |

Switches:

| Switch | Default | Effect |
|---|---|---|
| `openmp=0\|1` | 1 | OpenMP flag; links fesm-utils `include-serial` or `include-omp` and the matching LIS, FFTW, SHTns |
| `debug=1\|2` | 0 | debug or profiling flags of the compiler fragment |
| `pic=1` | 0 | `-fPIC` |
| `plastic` | 0 | ice backend `plastic`: `PlasticIceSheet/libplastic/include/libplastic.a` and LIS; off: stub, `test_ice_plastic` not built |
| `vilma` | 0 | geo backend `vilma2` (the global Earth and its regional ports): `vilma2/libvilma/include/libvilma.a`, SHTns of fesm-utils; off: stub, `test_geo_vilma2` and `test_gia_sequence` not built |
| `yelmo` | 1 | ice backend `yelmo` (+ FastHydrology, elsa, tracer, LIS, FFTW); off: stub, the Yelmo experiments are stubs too |
| `fastisostasy` | 1 | geo backend `fastisostasy`; off: stub |
| `chion`, `insol` | 1 | surface backend `chion` (needs `insol=1`); off: stub |
| `cism`, `pism`, `rembo` | 0 | stop: the adapters do not exist yet |

Build notes:

- All backend include dirs are on every compile line, after `libkryos/include` and fesm-utils, so a module name
  that exists in two libraries (`precision` in fesm-utils and FastIsostasy) resolves to fesm-utils. kryos' own modules
  are prefixed `kryos_`.
- `kryos-static` does not descend into the dependencies: it builds against what is built, and a missing archive stops
  the build with the command that builds it.
- A host links kryos as yelmox links yelmo: `-Ikryos/libkryos/include -Lkryos/libkryos/include -lkryos` plus the
  backend libraries. The archive is rewritten from scratch, so its mtime changes only when an object does.
- The marine_shelf and pico code of yelmox is compiled in kryos (`src/backends/shelf/marine_shelf/`); it needs no
  switch.

How the generated dependencies work and how to add a backend:
[docs/dev/build.md](https://github.com/fesmc/kryos/blob/main/docs/dev/build.md).

## The Levante tools

The scripts of [tools/levante/](https://github.com/fesmc/kryos/tree/main/tools/levante) share one environment,
[env.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/env.sh), which every script sources:

| Variable | Default | |
|---|---|---|
| `KRYOS` | the checkout the script is in | |
| `KRYOS_WORK` | `/work/ba1442/$USER/kryos-work` | your writable work folder: runs, scores, calibrations |
| `FESMDATA`, `FESMDATA_SOURCES` | `/work/ba1442/fesmdata`, `/work/ba1442/fesmdata-sources` | FesmData and its original data, read-only; `FESMDATA` is linked as `$KRYOS/fesmdata` |
| `KRYOS_DATA` | `/work/ba1442/kryos_data` | kryos reference data, read-only: refcorr topographies, the released tau_b products (`plastic/`), the shared remapping weights (`maps/`); linked as `$KRYOS/kryos_data` |
| `VILMA2` | `$KRYOS/vilma2` | the vilma2 checkout of the standalone runs (`bin/vilma.x`, `input/`) |
| `PM_DIR`, `PM_REF` | `$KRYOS_WORK/paleomist`, the reference runs | your standalone VILMA2 runs; the read-only reference runs the scores compare with |
| `ACCOUNT` | the `account` of `.runme/config.toml`, else `ba1442` | the Slurm account of the small Julia and analysis jobs of the tools (model runs take runme's) |
| `KRYOS_RUNME` | `runme` | the runme command: the model runs are staged and submitted with runme (below) |
| `KRYOS_JULIA` | juliaup's `julia +1.13` | the julia command (by full path: Levante's python module puts an old julia first on `PATH`) |

Override any of them in the environment or in `~/.kryos.env`. The pages write the paths with these variables
(`$KRYOS_WORK/...`); to use them in a shell, source the file (`. tools/levante/env.sh`). Then `kryos_jl <script.jl> ...`
runs julia with the tools' project.

Once, in a login shell (it needs your own juliaup with Julia 1.13:
`curl -fsSL https://install.julialang.org | sh -s -- --yes && juliaup add 1.13`, and runme:
`pip install --user git+https://github.com/fesmc/runme`):

```bash
tools/levante/setup.sh
```

It checks the paths and the builds, and sets up the checkout:

- the links `fesmdata -> $FESMDATA` and `kryos_data -> $KRYOS_DATA` at the kryos root (runme links them into every
  run folder);
- the grid files and remapping weights of `$KRYOS_DATA/maps` copied into the checkout's `maps/`, the maps folder of
  every run (so the runs reuse the weights);
- runme: installed and >= 0.5.27 (else `runme update`), and `.runme/config.toml` created from
  `.runme/config.default.toml` (`hpc = "dkrz_levante"`, `account = "ba1442"`; edit it for another account);
- `ice_data` and `isostasy_data`: read only by the Yelmo par files. Without them runme warns that the link is
  missing, which is harmless for the plastic runs;
- the Julia project of the tools, `tools/levante` (`Project.toml`, `Manifest.toml`: CairoMakie, NCDatasets,
  DelimitedFiles on Julia 1.13), instantiated in your own depot (a fresh instantiate took 3.5 min). Every Julia tool
  runs with `julia --project=tools/levante` (or `kryos_jl`).

It is safe to rerun.

## Data

A run reads its inputs through links in the run directory; the par files name paths relative to it
(`fesmdata/{domain}/{grid_name}/...`, with `{domain}` and `{grid_name}` replaced by kryos). runme makes the links
from those of the checkout ([.runme/info.json](https://github.com/fesmc/kryos/blob/main/.runme/info.json) `links`):

| Link | Target on Levante | What |
|---|---|---|
| `fesmdata` | `/work/ba1442/fesmdata` | FesmData products, `<Domain>/<GRID>/<GRID>_<PRODUCT>.nc`: topography, regions, PaleoMIST ice history and tau_b domains, the Bagge et al. (2021) viscosity |
| `kryos_data` | `/work/ba1442/kryos_data` | kryos-made inputs, e.g. `refcorr/PaleoMIST-sci/` (the present-day reference correction r, [tools/levante/gia_refcorr.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/gia_refcorr.sh)), the released tau_b products `plastic/` |
| `input` | the checkout's `input/` | defaults, variable tables, constants, synchronised copies (below) |
| `maps` | the checkout's `maps/` | grid descriptions `grid_<GRID>.txt` and the cached remapping weights (`*.nc`) |
| `ice_data`, `isostasy_data` | `/work/ba1442/ice_data`, `/work/ba1442/isostasy_data` | the yelmox data of the Yelmo configurations, linked when the checkout has them; the plastic par files do not read them |

**FesmData is read-only.** Never write into `/work/ba1442/fesmdata`, nor through a link into it (a copy of a link
is a link: copy with symlinks followed). It is shared by every run of the project; fetch or update products only with
FesmData's own tool (`julia fesmdata.jl`, see the FesmData
[download page](https://fesmc.github.io/FesmData/download.html)). FesmData products are the source of the inputs.

What the plastic par files read:

| File | Par files |
|---|---|
| `fesmdata/{domain}/{grid_name}/{grid_name}_IceHistory-PaleoMIST.nc` | every plastic run: the margin chronology (`f_ice`, 33 slices, 80 ka to present) |
| `fesmdata/{domain}/{grid_name}/{grid_name}_IceHistory-PaleoMIST_taub.nc` | every plastic run: Gowan's tau_b domains and their series |
| `fesmdata/{domain}/{grid_name}/{grid_name}_TOPO-GEBCO2026.nc`, `_TOPO-BedMachine-v4.nc` | the reference topography without correction (`kryos_North_plastic.nml`, the `_b1` runs) |
| `fesmdata/{domain}/{grid_name}/{grid_name}_REGIONS.nc` | every plastic run |
| `fesmdata/Global/{grid_name}/{grid_name}_IceHistory-PaleoMIST.nc` | the VILMA2 runs: the ice outside the plastic domains (`[globe_ice]`, GLOBAL-0.5DEG) |
| `fesmdata/Global/GLOBAL-T170/GLOBAL-T170_SolidEarth-Bagge2021.nc` | the science-Earth runs (3-D viscosity) |
| `kryos_data/refcorr/PaleoMIST-sci/...` | the science-Earth runs: the reference bed today's + r |

The maps folder must hold `grid_<GRID>.txt` for every hub and instance grid of the run; the first run writes the
remapping weights into it. kryos writes only VILMA2's Gauss grid file, `grid_GAUSS-L<lmax>.txt`. The repository's
[maps/](https://github.com/fesmc/kryos/tree/main/maps) has the grid files of ANT-16KM, ANT-32KM, GRL-16KM, NH-32KM,
PAT-PAL-8KM and GLOBAL-0.5DEG; for any other grid copy FesmData's description as it is
(`fesmdata/<Domain>/<GRID>/grid_<GRID>.txt`). The maps folder of every run is the checkout's `maps/` (runme's link):
`tools/levante/setup.sh` copies the grid files and the remapping weights of our runs from `$KRYOS_DATA/maps` into it,
so the plastic/VILMA2 runs need no new weights, and the first run of a new grid writes its weights there.

### input/

The run directory's `input/` is the repository's `input/`. kryos' own files: `kryos_defaults.nml` (the templates
of the par file groups), `forcing_defaults.nml`, `snapshots_defaults.nml`, `phys_const.nml` (the run's one set of
constants, `[kryos] phys_const_file`), the variable tables `kryos-variables-*.md` (with the output precision, see
Output) and the series of the yelmox configurations (`sl_deglaciation.dat`, `alpha_orbital_r2020.dat`,
`alpha_combined_125kyr_interp.dat`).

The synchronised copies of the files the libraries read, listed in
[input/sync.txt](https://github.com/fesmc/kryos/blob/main/input/sync.txt) with their source:

| Copies | Source |
|---|---|
| `yelmo_defaults.nml`, `yelmo_phys_const.nml`, `yelmo-variables-*.md` | `yelmo/input/` |
| `elsa_defaults.nml`, `tracer_defaults.nml` | `yelmo/elsa/input/`, `yelmo/tracer/input/` |
| `chion_defaults.nml` | `chion/input/` |
| `datasrc_defaults.nml` | `fesm-utils/par/` |
| `vilma_defaults.nml`, `vilma-variables.md` | `vilma2/input/` (the `[vilma]` defaults and the variable table of VILMA2's file) |
| `LA2004.INSOLN.5Ma.txt`, `LA2004.INSOLP.1Ma.txt` | `insol/input/` |

After updating a library: `make sync-input`, then commit the copies. Yelmo reads `sec_year` from
`input/yelmo_phys_const.nml` even though kryos passes it the run's constants, so the file must be present.

## A first run

### Running with runme

Every model run (`kryos.x`, the standalone `vilma.x`, the offline fits) is staged and submitted with
[runme](https://github.com/fesmc/runme), the fesmc tool that yelmox and climber-x use too. From a project folder
(here the kryos checkout, with its [.runme/](https://github.com/fesmc/kryos/tree/main/.runme)) runme creates the run
folder, copies the executable and the par file into it, writes the `-p` values into the copied par file, links the
input folders, writes the job script and submits it. The checkout's project:

- [.runme/info.json](https://github.com/fesmc/kryos/blob/main/.runme/info.json): the executables `kryos`
  (`libkryos/bin/kryos.x`, the default), `taub_fit` (`libkryos/bin/plastic_taub_fit.x`) and `closure_fit`
  (`libkryos/bin/plastic_closure_fit.x`); the defaults `input/kryos_defaults.nml`; the links `input`, `fesmdata`,
  `kryos_data`, `maps`, `ice_data`, `isostasy_data` (Data, above).
- `.runme/submit_slurm`, `submit_slurm_omp`: the job template (one node, one task, `--omp` cpus;
  `OMP_PROC_BIND=true`, `OMP_STACKSIZE=512M`, `OMP_MAX_ACTIVE_LEVELS=1`, `ulimit -s unlimited`). The job runs
  `./kryos.x <par file>` directly (one OpenMP task, no `srun`) in the run folder and then writes `status`, `<exit code> <seconds>`, which the
  calibration and ensemble tools read.
- `.runme/config.toml` (from `config.default.toml`, by `tools/levante/setup.sh`): `hpc = "dkrz_levante"`,
  `account = "ba1442"`, `omp = 16`.

Run runme from the checkout. A run of the science Earth (NH-ANT-PAT, 3-D) and a short 1-D (`_b1`) run:

```bash
. $HOME/models/kryos/tools/levante/env.sh       # KRYOS, KRYOS_WORK
cd $KRYOS
runme -r -s -q shared -w 08:00:00 --omp 64 -m 60G -o $KRYOS_WORK/runs/sci-3d -n par/kryos_NH-ANT-PAT_plastic_vilma2.nml
runme -r -s -q shared -w 01:00:00 --omp 16 -m 30G -o $KRYOS_WORK/runs/b1-short -n par/kryos_NH-ANT_plastic_vilma2_b1.nml \
    -p clock.time_end=-60e3
```

| Option | |
|---|---|
| `-o <folder>` | the run folder (an ensemble: the folder of the member folders `0`, `1`, ...) |
| `-n <par file>` | the par file, copied into the run folder under its own name and passed to the executable |
| `-e <exe>` | the executable alias of `info.json` (default `kryos`) |
| `-p group.key=value ...` | keys written into the copied par file; the key must exist in the par file or in `input/kryos_defaults.nml` (a key only in the defaults is added to its group); a comma list makes an ensemble, `key==a,b,c` or `'key=[a,b,c]'` a vector |
| `-s` | write the job script `job.submit` (with `-r`: submit it) |
| `-r` | run: with `-s` submit the job, without it run in the background on the login node |
| `-q`, `-w`, `--omp`, `-m` | the queue alias of `dkrz_levante` (`shared`, `compute`), the wall time, the OpenMP threads, the memory |
| `-i <table>`, `-j <ids>`, `--array [N]` | an ensemble from a parameter table, a subset of its members, and (with `-s`) its submission as one Slurm job array with at most N tasks at once (runme >= 0.5.27; writes `array.submit`, `array.dirs` and `array.json`, the job id, into `-o`) |

Without `-r` runme only stages the folder (and, with `-s`, writes `job.submit`: `sbatch job.submit` in the folder
submits it later). A folder that exists is reused: runme removes its `*.nc` and par files first.

Keep at most 30 jobs queued or running at a time on Levante; throttle job arrays (runme `--array N`).

### The run directory

```
<run>/
  kryos.x                a copy of libkryos/bin/kryos.x
  <par file>.nml         the copy of the par file, with the -p values
  input      -> $KRYOS/input
  fesmdata   -> /work/ba1442/fesmdata
  kryos_data -> /work/ba1442/kryos_data
  maps       -> $KRYOS/maps
  job.submit             the Slurm job (-s)
  runme.json             runme's record: the command, the -p values, the git revision of the checkout, the state
  params.txt, info.txt   the -p values as a table
  out.out, out.err       the log of the job
  status                 "<exit code> <seconds>", written when kryos.x ends
```

Copying `kryos.x` keeps the run independent of later builds. `kryos.x` takes one argument, the par file. The members
of a calibration or ensemble have the same layout, staged by runme from the calibration folder
([Calibration and ensembles](experiment-calibration.md)).

### Which par file first

| Par file | What | Cost (Levante) |
|---|---|---|
| `par/kryos_North_plastic.nml` | NH-32KM plastic on a fixed bed (today's, `z_sl = 0`), no globe: the smoke run | 25 s on 4 cores (800 solves) |
| `par/kryos_North_plastic_vilma2_b1.nml` | NH-32KM plastic with VILMA2 on Gowan's 1-D Earth and B1's reference | 13.6 min on 16 threads |
| `par/kryos_NH-ANT_plastic_vilma2_b1.nml` | as `_b1` plus a plastic instance on ANT-16KM | 13.4 min on 16 threads; 11-25 min on 8 threads (calibration members, reduced output) |
| `par/kryos_NH-ANT_plastic_vilma2.nml`, `par/kryos_NH-ANT-PAT_plastic_vilma2.nml` | the science Earth (Bagge et al. 2021 3-D viscosity, refcorr reference); `-PAT` adds PAT-PAL-8KM | |

All run 80 ka to present (`[clock] time_init = -80e3`, `time_end = 0`, `dtt = 100`). Costs from
[docs/dev/plastic-backend.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-backend.md). VILMA2
dominates and barely scales: 8 threads per run is the cheapest per run on the shared partition.

### Changing keys

Set keys with `-p group.key=value` when the run is staged (above). runme types the value by the key's current value
(a string is quoted for you) and stops on a key that neither the par file nor `input/kryos_defaults.nml` has. `bsl`
is an alias of `barysealevel`. The Yelmo configurations work the same way:

```bash
runme -rs -q compute -w 24:00:00 -e kryos -n par/kryos_Antarctica.nml -o output/ant --omp 16
runme -r -e kryos -n par/kryos_Greenland.nml -o output/grl-test --omp 4 -p clock.time_end=1000
```

To change the par file of a staged folder in place, edit it, or use
[tools/levante/set_nml.py](https://github.com/fesmc/kryos/blob/main/tools/levante/set_nml.py), which sets existing
keys only (a key the group does not have stops it; quote strings yourself):

```bash
python3 tools/levante/set_nml.py <run>/<par file>.nml clock.time_end=-70e3 'kryos.restart_out="none"'
```

## Output

Files go to `[kryos] out_dir` (`./`), one stream per instance and schedule group (`out_1D`, `out_2D`, ... name
schedule groups of the par file). Instance names are unique in a run, so the files have no domain prefix. For
`par/kryos_NH-ANT_plastic_vilma2_b1.nml`:

| File | Content |
|---|---|
| `ice_nor.nc`, `ice_ant.nc` | 2D exports of the ice instances (`out_2D`, every 2.5 kyr: `H_ice`, `z_srf`, `z_bed`, `z_sl`, `f_ice`, `smb`, ...) |
| `ice_nor_ts.nc`, `ice_ant_ts.nc` | their 0-D exports (`out_1D`, every 100 yr: `M_ice`, ...) |
| `ice_nor_plastic.nc`, `ice_ant_plastic.nc` | the backend's own stream `plastic` (`out_plastic`): `H_ice`, `z_srf`, `z_bed`, `z_sl`, `mask`, `sdist`, `tau_b`, `taub_domain`, and `converged`, `n_outer_used`, `stress_res_median`, `stress_res_max` of the last solve |
| `geo_glob.nc`, `geo_glob_ts.nc` | the global Earth's exports: `dz_bed`, `dz_ss`, `C_ocn` on the Gauss grid; `bsl`, `bsl_ive`, `mass_resid`, ... |
| `geo_glob_vilma.nc` | VILMA2's own diagnostics (`out_vilma`): `rsl`, `z_bed`, `dz_bed`, `bsl`, `bsl_ive`, ...; the file the GAPSLIP scoring reads |
| `geo_glob_gia.nc` | the global ice field as VILMA2 got it (`out_gia`; off in these files) |
| `geo_nor.nc`, `geo_nor_ts.nc`, ... | the regional geo ports on the ice grids |
| `globe_ice.nc` | the global ice field of the hub (`[globe_ice] out_2D`; off in these files) |
| `<domain>_hub.nc` | the domain's hub topography (`[<domain>] out_2D`; off in these files) |
| `timing.txt` | the timing table |
| `restart-<time>/` | restart bundles (below) |

The variables of each stream: the tables `input/kryos-variables-<component>.md` (exports),
`input/kryos-variables-ice-plastic.md` (the `plastic` stream), `input/kryos-variables-geo-vilma2.md` (the `gia`
stream) and `input/vilma-variables.md` (VILMA2's file).

### Precision

Output is rounded per variable (ncio bit rounding, `float`), as in yelmox; restart bundles are never rounded. The
`prec` column of the variable tables in `input/` is the only source of rounding: `-` lossless, `1e-3` absolute (in the
units of the variable), `1e-4r` relative, `1e-3,1e-4r` both. A variable a table does not list is written lossless;
a row without a unit, a bad `prec` entry, or a unit that differs from the data's stops the run. 0-D rows are lossless,
so the time series are.

| Table | Read by | Examples |
|---|---|---|
| `kryos-variables-<component>.md` | the default writer: `<instance>.nc`, `<instance>_ts.nc` | `ice`: `H_ice`, `z_bed` 1e-1 m, `f_ice` 1e-4, `smb` 1e-1 kg m-2 yr-1; `geo`: `dz_bed` 1e-2 m |
| `kryos-variables-ice-plastic.md` | the plastic backend: `<instance>_plastic.nc` | `H_ice`, `sdist`, `tau_b` 1e-1; `mask` 1e-4; `div_q` 1e-4 m yr-1 |
| `kryos-variables-geo-vilma2.md` | the vilma2 backend: `<instance>_gia.nc` | `H_ice_ext` 1e-3 m, `host_w` 1e-4 |
| `vilma-variables.md` | VILMA2: `<instance>_vilma.nc` | `rsl`, `z_bed`, `dz_bed` 1e-2 m; `bsl`, `bsl_ive` lossless |
| `kryos-variables-<component>-<backend>.md`, `kryos-variables-hub.md` | the other backends' streams, the hub file | `geo-fastisostasy`, `climate-snapshots`, `surface-chion`, `shelf-marine_shelf`; `hub` |
| `yelmo-variables-*.md` | Yelmo | the Yelmo streams |

`[kryos] write_lossless = True` ignores every precision (the gates set it). Rounded and lossless runs have identical
restart bundles, and a rounded file is within the tolerance of its row of the lossless one. On the 80-yr GRL-16KM run
of `test_run_sequence` the output files are 37 % smaller rounded.

## Restarts

When the schedule named by `[kryos] restart_out` is due (`out_rst`, every 10 kyr in the plastic par files), and
always at `time_end`, the run writes a bundle `<out_dir>/restart-<time>/` (`<time>` in years, shortest exact
decimal, e.g. `restart--40000` at 40 ka):

```
restart-<time>/
  clock_restart.nc                     the clock and the cadence of every instance that steps
  geo_glob_imports.nc  geo_glob_exports.nc  geo_glob_restart.nc    the global Earth (VILMA2's state)
  North/
    <inst>_imports.nc  <inst>_exports.nc   every instance of the domain
    <inst>_restart.nc                      stateful backends
  manifest.nml                         [bundle] + one group per instance; written last
```

A run with a forcing instance also writes `forcing_imports.nc`, `forcing_exports.nc` and
`forcing_<coef>_restart.nc` at the top level; a run without one (the plastic par files) writes no `forcing_*.nc`, and
the top level holds `clock_restart.nc`, the `geo_glob_*.nc` files and `manifest.nml`.

- The plastic backend's `<instance>_restart.nc` holds `H`, `w` and the mask `w` belongs to (the library's warm
  start), `H_prev` and `z_srf`, lossless; the read checks the grid and the time.
- `manifest.nml` holds `[bundle]` (`time`, `ticks`, `ticks_per_year`, `kryos_git`, `fesm_utils_git`, `domains`)
  and one group per instance (`scope`, `domain`, `component`, `backend`, `grid`, `files`, `git`). A bundle without it
  is incomplete.

Restart with `[kryos] restart` (a folder relative to the run directory, or absolute) and the clock continued, e.g.
a second chunk:

```bash
runme -r -s -q shared -w 08:00:00 --omp 64 -m 60G -o $KRYOS_WORK/runs/chunk2 -n par/kryos_NH-ANT-PAT_plastic_vilma2.nml \
    -p clock.time_init=-40e3 kryos.restart=../chunk1/restart--40000
```

(the Yelmo Greenland run likewise: `-p clock.time_init=15000 clock.time_end=30000 kryos.restart=../grl/restart-15000`).

- The restore stops, naming the instance, if the run's instances differ from the manifest's (scope, domain,
  component, backend, grid) or a listed file is missing; a changed `dt` is allowed.
- Restart times and `time_end` must be on every instance's cadence grid (checked at init).
- The fixed-bed plastic run (`mode = "flat"`, constant tau_b) restarted from its 40 ka bundle equals the continuous
  run bit for bit.

Forcing-increment bundles: a tsgen coefficient with `[coef_<name>] restart_every_df = <df>` (> 0) also writes
`<out_dir>/restart-<n>/` (`n` = 1, 2, ...) each time it has moved by more than `df` since its last bundle. A bundle
is never overwritten: a run that would write into an existing `restart-<n>` stops. The plastic par files have no
forcing instance.

## Next

- How a run is put together: [Architecture](architecture.qmd), [Global hub and regional domains](global-hub.qmd).
- The backends: [PlasticIceSheet](backend-plastic.md), [VILMA2](backend-vilma2.md),
  [smb_simple](backend-smb-simple.md).
- The experiments: [Reproducing Gowan et al. (2021)](experiment-gowan.md),
  [Simulating our own ice](experiment-plastic.md), [Calibration and ensembles](experiment-calibration.md).
- The par files and the Levante tools: [Par files](reference-par-files.md), [Tools](reference-tools.md).
