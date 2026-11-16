---
title: "kryos"
---

**kryos** is a component-based framework for ice sheets and the solid Earth. A run is a set of *domains*, each a grid
with its own component instances (ice, surface, climate, geo), coupled through one global solid Earth. Each component is
a contract with interchangeable *backends*: Yelmo or the perfectly-plastic ice sheet of
[PlasticIceSheet](https://github.com/fesmc/PlasticIceSheet) for the ice; VILMA2 or FastIsostasy for the Earth;
snapshots, chion or `smb_simple` for climate and surface. Outside the domains the global Earth sees an ice history read
from file. kryos is the successor of [YelmoX](https://fesmc.github.io/yelmox/).

The site has two parts:

- **The science**: the reconstruction of ice sheets and sea level of the last glacial cycle made with kryos, the
  material of a paper in preparation. [Start here](science/index.md).
- **The model and how to use it**: getting started, the how-to pages, how kryos works, the backends and the reference.

## Where to start

| you want to | read |
|---|---|
| understand the science and where it stands | [Reconstructing the last glacial cycle](science/index.md), then [Results](science/results.md) |
| install kryos and make a first run | [Getting started](getting-started.md) |
| run imposed ice, our ice online, an ensemble or the Earth search | the How-to pages: [imposed ice](howto/imposed-ice.md), [online ice](howto/online-ice.md), [ensembles](howto/ensembles.md), [the Earth search](howto/earth-search.md) |
| understand how a run is put together | [Architecture](architecture.qmd), [Global hub and regional domains](global-hub.qmd) |
| configure a backend | [PlasticIceSheet](backend-plastic.md), [VILMA2](backend-vilma2.md), [smb_simple](backend-smb-simple.md) |
| look up a par file or a tool | [Par files](reference-par-files.md), [Tools](reference-tools.md) |

## Rules of the road

- Builds, runs and analysis run on the HPC machines, DKRZ Levante and AWI's albedo; a laptop is for editing.
- Runs, scores and data go in the project's shared work folder (`KRYOS_WORK`: `/work/ba1442/kryos-work` on Levante,
  `/albedo/work/projects/p_forclima/kryos-work` on albedo), never in a personal folder.
- The FesmData store is read-only: never write into it except through FesmData's own tools.
- Jobs queued or running at a time, your own: at most 50 on albedo, 30 on Levante.
- Analysis and plots in Julia with CairoMakie.

Source: [github.com/fesmc/kryos](https://github.com/fesmc/kryos). Design notes, gate reports and plans are in the
repository's `docs/dev/` and `docs/study/`: see [Developer notes](developer-notes.md).
