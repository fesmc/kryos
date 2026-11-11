---
title: "kryos"
---

**kryos** is a component-based framework for ice sheets and the solid Earth. A run is a set of *domains*, each a hub
grid with its own component instances (ice, surface, climate, geo, ...), coupled through one global solid Earth. Each
component is a contract with interchangeable *backends*: Yelmo or the perfectly-plastic ice sheet of
[PlasticIceSheet](https://github.com/fesmc/PlasticIceSheet) for the ice; VILMA2 or FastIsostasy for the Earth;
snapshots, chion or `smb_simple` for climate and surface. Outside the domains, the global Earth sees an ice history
read from file. kryos is the successor of [YelmoX](https://fesmc.github.io/yelmox/).

These pages concentrate on the set-up we use for the last glacial cycle: Gowan et al.'s (2021) PaleoMIST ice
history, imposed on VILMA2 alone (the reproduction of Gowan's run) or simulated with the plastic ice sheet on several
regional domains (Northern Hemisphere, Antarctica, Patagonia) coupled to the global VILMA2 Earth with a 1-D or a 3-D
viscosity structure, and scored against the GAPSLIP 2.1 sea-level index points.

## Where to start

- [Getting started](getting-started.md): the code and its dependencies on Levante, building, the data, a first run.
- [Architecture](architecture.qmd) and [Global hub and regional domains](global-hub.qmd): how a run is put
  together and how the regional ice reaches the global Earth and back.
- Backends: [PlasticIceSheet](backend-plastic.md), [VILMA2](backend-vilma2.md),
  [PaleoEmulator](backend-paleoemu.md), [smb_simple](backend-smb-simple.md).

## Experiments

| Page | What it reproduces |
|------|--------------------|
| [Reproducing Gowan et al. (2021)](experiment-gowan.md) | The PaleoMIST ice imposed on standalone VILMA2 (1-D Earth), the GAPSLIP comparison |
| [Simulating our own ice](experiment-plastic.md) | The plastic ice sheet online with VILMA2, 1-D or 3-D Earth, two or three domains |
| [Calibration and ensembles](experiment-calibration.md) | Present-day tau_b factors, ES-MDA, the weighted ensemble |
| [Analysis and figures](experiment-analysis.md) | Maps, ensemble spread, site comparisons, residual correlation |
| [Implied mass balance](experiment-implied-mb.md) | The mass balance a reconstructed ice sheet needs (in development) |

## Rules of the road

- Builds, tests, runs and analysis run on DKRZ Levante; the FesmData store (`/work/ba1442/fesmdata`) is read-only:
  never write into it except through FesmData's own tools.
- At most 30 jobs queued or running at a time on Levante (throttle job arrays).
- Plots and analysis are in Julia with CairoMakie.

Source: [github.com/fesmc/kryos](https://github.com/fesmc/kryos). Developer notes (design, gate reports, plans) are
in the repository's `docs/dev/` and `docs/study/`: see [Developer notes](developer-notes.md).
