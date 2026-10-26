# kryos documentation

Design settled 2026-10-24 (thirteen decisions in design.md section 15). Build
system and source skeleton in place; components are placeholders.

- [build.md](build.md): configure, targets, switches, generated dependencies,
  adding a backend with a stub.
- [running.md](running.md): par files and defaults, the synchronised `input/`,
  running with runme, output files, restart bundles.
- [design.md](design.md): the proposal — vocabulary, architecture, component
  contract, exchange, regional vs global, time, data backends, what is kept
  from yelmox and CLIMBER-X and why, what is shed and why, additions to
  fesm-utils, repository layout, migration, the thirteen decisions.
- `study/`: the evidence, one document per codebase read, with file:line
  references as of 2026-10-24:
  - [yelmox.md](study/yelmox.md): the current framework and its in-tree
    `kryos*` prototype; interfaces per component; leaks; what worked; legacy.
  - [climber-x.md](study/climber-x.md): the `ice_model` interface, the `cmn`
    hub, regional domains inside a global model, geo/sea level/VILMA, dummies.
  - [fesm-utils.md](study/fesm-utils.md): package inventory, the grid/map/coupler
    stack and its global↔regional gaps, time, I/O, nml, candidates to add.
  - [backends.md](study/backends.md): Yelmo, CISM, PISM, FastIsostasy, VILMA2,
    chion, smbpal, smb_simple, rembo1, insol, marine_shelf/PICO, obm,
    geothermal, sediments, FastHydrology: needs, provides, clocks, units,
    signs, instantiability; the minimum interface per component.
