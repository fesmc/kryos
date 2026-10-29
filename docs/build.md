# Building kryos

kryos follows the build convention of the stack (yelmo, chion, yelmox): a
`config/Makefile` template that configme fills in for a machine and compiler,
dependency wiring in `config/common.mk`, a static library as the product. What
is new: the per-object dependency rules are generated (`config/deps.mk`), and
every backend library is optional behind a stub.

## Configure

Dependencies are sibling checkouts or symlinks at the repo root:
`fesm-utils`, `yelmo`, `FastIsostasy`, `chion`, `insol`, `vilma2` (pinned in
`.configme/manifest.toml`; gitignored). `configme install kryos` clones them;
for existing checkouts, symlink them (`ln -s ../fesm-utils fesm-utils`).

```bash
configme config kryos --only -m macbook -c gfortran
```

This writes the root `Makefile` (template + compiler/machine fragments +
netCDF) and `config/deps.mk` (needs configme >= 0.16). `--only` leaves the
dependencies' own Makefiles alone.

`vilma=1` (geo backend `vilma2`) needs `make vilma-static` in `vilma2/` (or
`make vilma-static` here, which builds fesm-utils first) and the SHTns that
fesm-utils builds once (`make -C fesm-utils shtns openmp=0|1`; `make deps-check`
names it when missing). VILMA2, fesm-utils and kryos must use the same `openmp=`.
`make check` needs the `vilma2` checkout whatever the switch: two of the
synchronised copies in `input/` (`vilma_defaults.nml`, `vilma-variables.md`) come from it.

## Targets

| Target | Builds |
|---|---|
| `make kryos-static` | `libkryos/include/libkryos.a` and the `.mod`/`.smod` files in `libkryos/include` |
| `make kryos` | `libkryos/bin/kryos.x`, the driver (`drivers/kryos.f90`) |
| `make check` | `make input-check`, then every `tests/test_*.f90`, linked against libkryos, then runs them; fails if one exits non-zero |
| `make check-long` | the long tests, left out of `make check`: `tests/test_run_domains.f90` (two domains, GRL-16KM + ANT-32KM, about 8 min; needs Julia with the project of `tools/`, else it reports OPEN and passes; its header lists the checks) |
| `make sync-input` | copies the library sources over the synchronised copies in `input/` (`input/sync.txt`; [running.md](running.md)) |
| `make input-check` | fails naming every copy in `input/` that differs from its source |
| `make deps-check` | fails with instructions if a dependency archive is missing |
| `make deps-static` | builds fesm-utils and the enabled backend libraries in their own trees |
| `make depend` | regenerates `config/deps.mk` (`configme deps .`) |
| `make clean` / `make cleanall` | kryos objects and executables / also the dependencies |
| `make` / `make usage` | lists the targets and the current switches |

A host links kryos as yelmox links yelmo: `-Ikryos/libkryos/include
-Lkryos/libkryos/include -lkryos` and the same backend libraries, and makes its
objects depend on `libkryos.a` so a type change rebuilds them. The archive is
rewritten from scratch, so its mtime changes only when an object does.

`kryos-static` does not descend into the dependencies: it builds against
whatever is built, and a missing archive stops the build with the command that
builds it. Use `make deps-static` (or each dependency's own make) first.

## Switches

| Switch | Default | Effect |
|---|---|---|
| `openmp=0\|1` | 1 | OpenMP flag; links fesm-utils `include-serial` or `include-omp` (and the matching LIS/FFTW). The host must use the same variant. |
| `pic=1` | 0 | `-fPIC`, for a later shared library or C API |
| `debug=1\|2` | 0 | debug or profiling flags from the compiler fragment |
| `yelmo` | 1 | ice backend `yelmo`: adapter, `INC_YELMO`/`LIB_YELMO` (+ FastHydrology, elsa, tracer, LIS, FFTW) |
| `fastisostasy` | 1 | geo backend `fastisostasy` (regional and global part): adapter, `INC_ISOSTASY`/`LIB_ISOSTASY` (+ LIS, FFTW); off: stub, and `tests/test_geo_fastisostasy.f90` is not built (`test_skip`) |
| `chion`, `insol` | 1 | surface backend `chion`: `chion=1` compiles the chion host (with its orbit) into the adapter and needs `insol=1`; off: stub |
| `vilma` | 0 | geo backend `vilma2` (the global Earth and the regional port): adapter, `INC_VILMA`/`LIB_VILMA` (`libvilma.a` of `vilma2/`), SHTns of fesm-utils (`LIB_SHTNS`; LIS, FFTW); off: stub, and `tests/test_geo_vilma2.f90` is not built (`test_skip`), `tests/test_geo_vilma2_stub.f90` is |
| `cism`, `pism`, `rembo` | 0 | stop with the phase in which the adapter arrives |

Objects are stamped with the switch settings
(`libkryos/include/.build-openmp1-...-rembo0`); changing any switch rebuilds
every object, since a stub replaces a module others were compiled against.

## Generated dependencies

`config/deps.mk` is written by configme from the `module`, `submodule`, `use`
and `include` statements under `src/` (kryos is registered with
`[orchestrator.fortran_deps] src = ["src"]`, and `config/Makefile` includes
the file). It contains:

- `vpath` for every source directory, so one pattern rule
  `$(objdir)/%.o: %.f90` compiles every file;
- `deps_objs`, the library objects (the Makefile uses it as `kryos_objs`);
- one rule per object, `$(objdir)/kryos_ice.o: $(objdir)/kryos_component.o`;
- for a module defined in two files (an adapter and its stub), a check that
  `<module>_obj` is set, and `$(<module>_obj)` wherever the module is used.

Modules not defined under `src/` (fesm-utils, yelmo, intrinsics) are external
and get no rule. Object names must be unique across `src/` (one flat object
dir). Rerun `configme config` or `make depend` after adding a file or a `use`.

## Adding a backend

For backend `foo` of component `ice`:

1. `src/backends/ice/kryos_ice_foo.f90`: `module kryos_ice_foo`, type
   `foo_ice_t` extending `ice_t`, `use`s the library.
2. `src/backends/ice/kryos_ice_foo_stub.f90`: the same module and type, no
   library; every procedure stops with "rebuild kryos with foo=1".
3. The factory (`src/components/ice/kryos_ice_factory.f90`): `use
   kryos_ice_foo` and a `case ("foo")` that allocates `foo_ice_t`.
4. `config/common.mk`: `FOOROOT`, `INC_FOO`, `LIB_FOO`, `FOO_ARCHIVE`, a
   `foo ?= 0` switch whose `ifeq` block sets `kryos_ice_foo_obj` to the
   adapter or the stub and, when on, appends to `INC_BACKENDS`,
   `LIB_BACKENDS`, `dep_archives`, `backend_roots`, `backend_static`; add
   `foo` to `backend_config`.
5. `config/Makefile`: a missing-archive rule for `$(FOO_ARCHIVE)` and a
   `foo-static` target; add `foo` to `usage`.
6. `configme config` (or `make depend`), then `make check foo=0` and
   `make check foo=1`.

## Notes

- All backend include dirs are on every compile line, after
  `libkryos/include` and fesm-utils, so a bare module name that exists in
  two libraries (`precision` in fesm-utils and FastIsostasy) resolves to
  fesm-utils. kryos' own modules are all prefixed `kryos_`.
- The marine_shelf and pico code of yelmox is compiled in kryos (`src/backends/shelf/marine_shelf/`, modules `kryos_mshlf*`); it needs no switch. The chion host and orbit modules are inside `kryos_surface_chion.f90`.
- The yelmo adapter compiles against `yelmo.mod` and holds a `yelmo_class`,
  but the placeholder calls nothing in libyelmo, so no libyelmo member is
  linked yet.
- macOS ships GNU make 3.81, whose timestamps have one-second resolution: a
  switch changed within the same second as the last compile may not rebuild
  every object. `make clean` if in doubt.
