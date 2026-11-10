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

The gate builds used the modules of the Levante login profile: `intel-oneapi-compilers/2023.2.1`,
`netcdf-c/4.9.3`, `netcdf-fortran/4.6.2` (openmpi-5.0.8, intel-2021.5.0) and `python3/2023.01`
([gate-log](https://github.com/fesmc/kryos/blob/main/docs/dev/gate-log.md)). The production compiler is Intel
`ifx` (configme machine `dkrz_levante`, compiler `ifx`).

configme configures and clones the stack; install it once:

```bash
pip install git+https://github.com/fesmc/configme
export PATH="${PATH}:${HOME}/.local/bin"      # if `configme` is not found
```

<!-- TODO(check): the module names and versions above are from the phase 0 gate log (2026-10); confirm they are still the login-profile defaults. -->

## The code

kryos and its dependencies are sibling checkouts, or symlinks, at the kryos root (gitignored). The pins are in
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

Clone and configure the stack with configme, then add PlasticIceSheet:

```bash
cd ~/models
configme install kryos -m dkrz_levante -c ifx --dir ~/models/kryos   # -d https without a GitHub SSH key
git clone git@github.com:fesmc/PlasticIceSheet.git
ln -s ../PlasticIceSheet kryos/PlasticIceSheet
```

<!-- TODO(check): that `configme install kryos` clones vilma2 (it is in the repo's manifest, but not in configme's seed default_packages for kryos) and which dependencies it builds. -->

For checkouts that already exist, symlink them instead (`ln -s ../fesm-utils fesm-utils` at the kryos root, or
`configme install kryos --link fesm-utils=/abs/path`).

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

## Data

A run reads its inputs through links in the run directory; the par files name paths relative to it
(`fesmdata/{domain}/{grid_name}/...`, with `{domain}` and `{grid_name}` replaced by kryos).

| Link | Target on Levante | What |
|---|---|---|
| `fesmdata` | `/work/ba1442/fesmdata` | FesmData products, `<Domain>/<GRID>/<GRID>_<PRODUCT>.nc`: topography, regions, PaleoMIST ice history and tau_b domains, the Bagge et al. (2021) viscosity |
| `kryos_data` | `/work/ba1442/kryos_data` | kryos-made inputs, e.g. `refcorr/PaleoMIST-sci/` (the present-day reference correction r, [tools/levante/gia_refcorr.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/gia_refcorr.sh)) |
| `input` | the repository's `input/` | defaults, variable tables, constants, synchronised copies (below) |
| `maps` | a writable folder of your own | grid descriptions `grid_<GRID>.txt` and the cached remapping weights (`*.nc`) |
| `ice_data`, `isostasy_data` | `/work/ba1442/ice_data`, `/work/ba1442/isostasy_data` | the yelmox data of the Yelmo configurations; the plastic par files do not read them |

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
remapping weights into it (VILMA2 also writes its `grid_GAUSS-L<lmax>.txt`). The repository's
[maps/](https://github.com/fesmc/kryos/tree/main/maps) has `grid_ANT-16KM.txt`, `grid_ANT-32KM.txt` and
`grid_GRL-16KM.txt`; FesmData has the description of each grid in its grid folder
(`fesmdata/<Domain>/<GRID>/grid_<GRID>.txt`). `par/kryos_NH-ANT-PAT_plastic_vilma2.nml` needs `grid_PAT-PAL-8KM.txt`. For the plastic/VILMA2 runs link the maps folder
of our runs, `/work/ba1442/robinson/kryos-ice-plastic-maps`: it has the grid files of NH-32KM, ANT-16KM, PAT-PAL-8KM,
GLOBAL-0.5DEG and GAUSS-L170 and all their remapping weights.

<!-- TODO(check): that FesmData's grid_<GRID>.txt can be copied into a new maps folder as it is. -->

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

### The run directory

```
<run>/
  kryos.x                a copy of libkryos/bin/kryos.x
  par.nml                a copy of the par file (with any overrides)
  input/                 the repository's input/ (a link, or a copy to freeze it)
  fesmdata   -> /work/ba1442/fesmdata
  kryos_data -> /work/ba1442/kryos_data
  maps       -> <your maps folder>
  job.submit
```

This is the layout the Levante tools stage, e.g. a member of the tau_b calibration
([tools/levante/plastic_esmda.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_esmda.sh)) and
[tools/levante/calib.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/calib.sh). Copying `kryos.x`
keeps the run independent of later builds. `kryos.x` takes one argument, the par file.

```bash
K=$HOME/models/kryos
R=<your work folder>/runs/nor-plastic-b1
mkdir -p $R && cd $R
cp $K/libkryos/bin/kryos.x .
cp $K/par/kryos_North_plastic_vilma2_b1.nml par.nml
ln -s $K/input input
ln -s /work/ba1442/fesmdata fesmdata
ln -s /work/ba1442/kryos_data kryos_data
ln -s <your maps folder> maps
git -C $K rev-parse HEAD > kryos_commit.txt
```

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

### The job

A job script after those of the calibration (8 threads, 10 GB, shared partition):

```bash
#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=ba1442
#SBATCH --time=01:00:00
#SBATCH --job-name=kryos-plastic
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --output=./out.out
#SBATCH --error=./out.err
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=true
export OMP_STACKSIZE=512M
export OMP_MAX_ACTIVE_LEVELS=1
ulimit -s unlimited
srun ./kryos.x par.nml
```

```bash
cd $R && sbatch job.submit
```

Keep at most 30 jobs queued or running at a time on Levante; throttle job arrays (`--array=...%N`).

### Changing keys

Edit `par.nml` in the run directory, or set keys with `set_nml.py` as the Levante tools do (a string value keeps
its quotes):

```bash
python3 tools/levante/set_nml.py par.nml clock.time_end=-70e3 'kryos.restart_out="none"' 'ice_nor.out_plastic="none"'
```

[tools/levante/set_nml.py](https://github.com/fesmc/kryos/blob/main/tools/levante/set_nml.py) sets existing keys only (a key
the group does not have stops it). For the plastic runs
[tools/levante/plastic_run.sh](https://github.com/fesmc/kryos/blob/main/tools/levante/plastic_run.sh) stages and submits a
run folder in one step ([Simulating our own ice](experiment-plastic.md)).

### With runme

runme stages and submits from the repository root
([.runme/info.json](https://github.com/fesmc/kryos/blob/main/.runme/info.json): executable `kryos` =
`libkryos/bin/kryos.x`, defaults `input/kryos_defaults.nml`). Once per checkout: `runme config init`, then set `hpc`
and `account` in `.runme/config.toml`.

```bash
runme -rs -q compute -w 24:00:00 -e kryos -n par/kryos_Antarctica.nml -o output/ant --omp 16
runme -r -e kryos -n par/kryos_Greenland.nml -o output/grl-test --omp 4 -p clock.time_end=1000
```

`-p group.key=value` overrides a key (`bsl` is an alias of `barysealevel`); a comma list makes an ensemble
(`-p ydyn.beta_q=0.2,0.5,1.0`). runme links `input`, `ice_data`, `fesmdata`, `isostasy_data` and `maps` from the
repository root into the run directory, but not `kryos_data`: the par files that read `kryos_data/` (the
science-Earth runs) need that link added by hand.

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
`forcing_<coef>_restart.nc` at the top level.

<!-- TODO(check): that a run without a forcing instance (the plastic par files) writes no forcing_*.nc into the bundle. -->

- The plastic backend's `<instance>_restart.nc` holds `H`, `w` and the mask `w` belongs to (the library's warm
  start), `H_prev` and `z_srf`, lossless; the read checks the grid and the time.
- `manifest.nml` holds `[bundle]` (`time`, `ticks`, `ticks_per_year`, `kryos_git`, `fesm_utils_git`, `domains`)
  and one group per instance (`scope`, `domain`, `component`, `backend`, `grid`, `files`, `git`). A bundle without it
  is incomplete.

Restart with `[kryos] restart` (a folder relative to the run directory, or absolute) and the clock continued, e.g.
a second chunk:

```bash
python3 set_nml.py par.nml clock.time_init=-40e3 'kryos.restart="../chunk1/restart--40000"'
```

With runme: `-p clock.time_init=15000 clock.time_end=30000 kryos.restart=../grl/restart-15000`.

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

- How a run is put together: [Architecture](architecture.md), [Global hub and regional domains](global-hub.md).
- The backends: [PlasticIceSheet](backend-plastic.md), [VILMA2](backend-vilma2.md),
  [smb_simple](backend-smb-simple.md).
- The experiments: [Reproducing Gowan et al. (2021)](experiment-gowan.md),
  [Simulating our own ice](experiment-plastic.md), [Calibration and ensembles](experiment-calibration.md).
- The par files and the Levante tools: [Par files](reference-par-files.md), [Tools](reference-tools.md).
