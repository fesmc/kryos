---
title: "Developer notes"
---

The design notes, plans, studies and gate reports behind kryos live in the repository, not on this site: they record decisions
and checks with file and line references as of their date, and are kept as written. The pages of this site describe the code
as it is; where they disagree with a note, the code and these pages are current. Status of each note as of 2026-11-16.

## The framework

| Note | Subject | Status |
|------|---------|--------|
| [design.md](https://github.com/fesmc/kryos/blob/main/docs/design.md) | The design proposal: vocabulary, architecture, contracts, exchange, regional vs global, time, the decisions | settled 2026-10-24; implemented |
| [dev/core-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/core-design.md) | The core (phase 1) | proposal of 2026-10-24; implemented |
| [dev/contracts-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/contracts-design.md) | Component contracts, wiring and the hub | proposal of 2026-10-24; implemented |
| [dev/adapters-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/adapters-design.md) | Backend adapters | design of 2026-10-24; implemented |
| [dev/regression-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/regression-design.md) | Regression harness and configuration migration | proposal (phase 1); used by the gates |
| [dev/build.md](https://github.com/fesmc/kryos/blob/main/docs/dev/build.md) | The build: generated dependencies, adding a backend | current |

## Plastic ice and VILMA2

| Note | Subject | Status |
|------|---------|--------|
| [dev/plastic-backend.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-backend.md) | The ice backend `plastic`: keys, tests, Stage 1, the ES-MDA calibrations, the 3-D Earth, the implied mass balance | record; ES-MDA retired |
| [dev/vilma2-paleomist.md](https://github.com/fesmc/kryos/blob/main/docs/dev/vilma2-paleomist.md) | VILMA2 standalone under the PaleoMIST ice history (B1, B3), the present-day correction r | record; B1 superseded by I1 |
| [dev/plastic-ensemble-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-ensemble-design.md) | The weighted ensemble: prior, likelihood, weights, GNSS, the constraints on the ice | implemented (ens-v1, ens-v2) |
| [dev/solid-earth-literature.md](https://github.com/fesmc/kryos/blob/main/docs/dev/solid-earth-literature.md) | Solid-Earth structure for GIA: Bagge 2021, Lloyd 2024, Pan 2022, seismic viscosity, regional constraints, GNSS | literature review, 2026-11-12 |
| [dev/earth-search-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/earth-search-design.md) | The broad Earth search with emulation: parameters, waves A, B and B2, the emulator, the scores | decided; wave A done, wave B emulated, wave B2 running |
| [dev/final-ensemble-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/final-ensemble-design.md) | Our ice on sampled Earths: a wide τ_b prior, an ice-only screen, RSL from step responses | proposed; decisions open |
| [dev/plastic-smb-inversion-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-smb-inversion-design.md) | The implied mass balance by joint inversion | implemented at present day, paused |
| [dev/plastic-margin-sources.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-margin-sources.md) | Margin targets from original sources | not pursued (2026-11-09) |

## Plans and gates

| Note | Subject | Status |
|------|---------|--------|
| [dev/phase2-plan.md](https://github.com/fesmc/kryos/blob/main/docs/dev/phase2-plan.md), [phase3](https://github.com/fesmc/kryos/blob/main/docs/dev/phase3-plan.md), [phase4](https://github.com/fesmc/kryos/blob/main/docs/dev/phase4-plan.md) | Phase plans | reviewed 2026-10-26 to 10-29 |
| [dev/gate-log.md](https://github.com/fesmc/kryos/blob/main/docs/dev/gate-log.md), [dev/gate-reports/](https://github.com/fesmc/kryos/tree/main/docs/dev/gate-reports) | Phase 0 sign-off on Levante and the regression reports (yelmox against kryos) | record |
| [dev/calib-cf.md](https://github.com/fesmc/kryos/blob/main/docs/dev/calib-cf.md), [ant-friction](https://github.com/fesmc/kryos/blob/main/docs/dev/ant-friction.md), [setup-diagnosis](https://github.com/fesmc/kryos/blob/main/docs/dev/setup-diagnosis.md), [setup-diagnosis-tables](https://github.com/fesmc/kryos/blob/main/docs/dev/setup-diagnosis-tables.md) | Yelmo basal friction and the setup of the Yelmo science case | record |

## Studies

The codebases read for the design (2026-10-24): [yelmox](https://github.com/fesmc/kryos/blob/main/docs/study/yelmox.md),
[CLIMBER-X](https://github.com/fesmc/kryos/blob/main/docs/study/climber-x.md),
[fesm-utils](https://github.com/fesmc/kryos/blob/main/docs/study/fesm-utils.md),
[the candidate backends](https://github.com/fesmc/kryos/blob/main/docs/study/backends.md).

## This site

The pages are Markdown in `docs/`, rendered with [Quarto](https://quarto.org) (`docs/_quarto.yml`); the workflow
`.github/workflows/docs.yml` renders them on every pull request and publishes them to GitHub Pages on every push to `main`
that touches `docs/`. The framework pages (Getting started, Architecture, Global hub, the backends, Reference) and the how-to
pages (`docs/howto/`) belong to kryos. The science pages (`docs/science/`) are the reconstruction of the last glacial cycle that
kryos serves, kept apart so that they can move to the repository of the paper. To preview locally: `quarto preview docs`.
