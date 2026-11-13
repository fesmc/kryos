# experiments

The launch scripts of the glacial-cycle simulations of
[docs/simulations.md](../docs/simulations.md): every row of its ladder is made by one call of one script here, and
scored by `scores.sh`. The workflow, the concepts and the open questions:
[docs/plastic-workflow.md](../docs/plastic-workflow.md).

All scripts run on Levante or albedo from a kryos checkout installed with `configme install kryos` and built with
`make kryos plastic=1 vilma=1 yelmo=0` (docs/getting-started.md). They source `tools/levante/env.sh`, which sets
`KRYOS_WORK` (your writable folder, where everything below goes), `FESMDATA`, `KRYOS_DATA`, `VILMA2` and the queues.

| ids | script | makes |
|---|---|---|
| P1, P2 | none | the published reconstructions are FesmData products (`IceHistory-PaleoMIST`, `IceHistory-GLAC1D`), scored as they are |
| I1-I8 | `refs-ice.sh <id>` | an ice history imposed on VILMA2 with one Earth: `$KRYOS_WORK/kryos_data/refcorr/<ice>_<earth>` |
| T1, D1, D3 | `refs-tau.sh <ens-v2> <id>` | kryos online from ens-v2's m000: `$KRYOS_WORK/runs/pm-online/{gowan,pdfac}`, `$KRYOS_WORK/runs/res16/m000-NH16` |
| D2, E2 | `ens.sh v2` | ens-v2, 50 members (m000 = D2): `$KRYOS_WORK/ens-v2` |
| E1 | `ens.sh v1` | ens-v1, 25 members (on Levante) |
| Earth-only | `earth-pca.sh basis v1`, `sample`, `run` | the PCA basis of the expert Earths, 30 sampled Earths, one pass each: `.../refcorr/ice-PMIST_E3D-PCA-E<nnn>` |
| S0-S4 | `earth-search.sh sample-a`, `run-a`, `score`, `sample-b`, `run-b 1d`, `run-b 3d`, `score b` | the Earth search: wave A, 200 1-D Earths under PaleoMIST and GLAC-1D (402 cases); wave B, 50 3-D Earths and the 1-D runs of their profiles (200 cases); scores and emulators: `$W/{earths,earthsB,kryos_data/refcorr,score,emulate,emulateB}` |
| scores | `scores.sh <ensemble folder>` | `<folder>/iter_01/scores.csv`: every member and every reference of the folder's config |

Order for a full rebuild: `refs-ice.sh all` and `ens.sh v2` (independent), then `refs-tau.sh <ens-v2> all` and
`earth-pca.sh` (they need ens-v2's m000 and the basis), then `scores.sh <ens-v2>` once the runs are done. The references
the scores read are the `[[references]]` of `tools/levante/plastic_ensemble_v2.toml` (copied into the ensemble folder at
its start); their `run` paths are the folders these scripts write.

Each `gia_refcorr.sh` case records its inputs and the code versions in `<case>/provenance.txt`, and each kryos run its
build in its output (`kryos_build_info`).
