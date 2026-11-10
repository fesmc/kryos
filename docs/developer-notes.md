---
title: "Developer notes"
---

The design notes, plans and gate reports behind kryos live in the repository, not on this site: they record decisions
and checks with file and line references as of their date, and are kept as written. The pages of this site describe
the code as it is; where they disagree with a note, the code and these pages are current.

## Design

| Note | Subject |
|------|---------|
| [design.md](https://github.com/fesmc/kryos/blob/main/docs/design.md) | The design proposal: vocabulary, architecture, contracts, exchange, regional vs global, time, the decisions |
| [dev/core-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/core-design.md) | The core (phase 1) |
| [dev/contracts-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/contracts-design.md) | Component contracts, wiring and the hub |
| [dev/adapters-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/adapters-design.md) | Backend adapters |
| [dev/regression-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/regression-design.md) | Regression harness and configuration migration |

## Plastic ice and VILMA2

| Note | Subject |
|------|---------|
| [dev/plastic-backend.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-backend.md) | The ice backend `plastic`: keys, tests, calibration results, the implied mass balance |
| [dev/vilma2-paleomist.md](https://github.com/fesmc/kryos/blob/main/docs/dev/vilma2-paleomist.md) | VILMA2 under the PaleoMIST ice history (Gowan et al. 2021) |
| [dev/plastic-margin-sources.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-margin-sources.md) | Margin targets from original sources (future option) |
| [dev/plastic-ensemble-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-ensemble-design.md) | The weighted ensemble (design draft) |
| [dev/plastic-smb-inversion-design.md](https://github.com/fesmc/kryos/blob/main/docs/dev/plastic-smb-inversion-design.md) | The implied mass balance by joint inversion (design draft) |

## Plans and gates

| Note | Subject |
|------|---------|
| [dev/phase2-plan.md](https://github.com/fesmc/kryos/blob/main/docs/dev/phase2-plan.md), [phase3](https://github.com/fesmc/kryos/blob/main/docs/dev/phase3-plan.md), [phase4](https://github.com/fesmc/kryos/blob/main/docs/dev/phase4-plan.md) | Phase plans |
| [dev/gate-log.md](https://github.com/fesmc/kryos/blob/main/docs/dev/gate-log.md) | Phase 0 sign-off on Levante |
| [dev/calib-cf.md](https://github.com/fesmc/kryos/blob/main/docs/dev/calib-cf.md), [ant-friction](https://github.com/fesmc/kryos/blob/main/docs/dev/ant-friction.md), [setup-diagnosis](https://github.com/fesmc/kryos/blob/main/docs/dev/setup-diagnosis.md) | Yelmo basal friction and setup of the science case |

## Studies

The codebases read for the design (2026-10-24): [yelmox](https://github.com/fesmc/kryos/blob/main/docs/study/yelmox.md),
[CLIMBER-X](https://github.com/fesmc/kryos/blob/main/docs/study/climber-x.md),
[fesm-utils](https://github.com/fesmc/kryos/blob/main/docs/study/fesm-utils.md),
[the candidate backends](https://github.com/fesmc/kryos/blob/main/docs/study/backends.md).

## This site

The pages are Markdown in `docs/`, rendered with [Quarto](https://quarto.org) (`docs/_quarto.yml`); the workflow
`.github/workflows/docs.yml` renders them on every pull request and publishes them to GitHub Pages on every push to
`main` that touches `docs/`. To preview locally: `quarto preview docs`.
