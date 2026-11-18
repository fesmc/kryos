# HPC environment of kryos (Levante, albedo)

| file | what |
|---|---|
| `env.sh` | the environment of a kryos checkout on Levante or AWI's albedo: `KRYOS`, `KRYOS_WORK` (the project's shared work folder), `FESMDATA`, `KRYOS_DATA`, `VILMA2`, `PM_DIR`, `PM_REF`, `ACCOUNT`, the queues, `KRYOS_RUNME`, `KRYOS_JULIA`; override any of them in the environment or in `~/.kryos.env` |
| `setup.sh` | the one-time setup of the checkout: checks the paths and the builds, links `fesmdata` and `kryos_data`, seeds `maps/` with the shared weights, runme and `.runme/config.toml` |
| `set_nml.py` | sets existing keys of a fully expanded par file in place (`set_nml.py <par file> group.key=value ...`) |
| `archive/` | the gate, regression and Yelmo-calibration scripts of the development phases, kept as they were run (their paths are those of their time; not maintained) |

```bash
tools/levante/setup.sh
```

Model runs (`kryos.x`, `vilma.x`) are staged and submitted with [runme](https://github.com/fesmc/runme) from the
checkout (its `.runme/`), e.g.

```bash
runme -r -s -q shared -w 08:00:00 --omp 64 -m 60G -o $KRYOS_WORK/runs/sci-3d -n par/kryos_NH-ANT-PAT_plastic_vilma2.nml
```

The experiments, scores, ensembles and analysis tools of the glacial-cycle reconstruction are in
[awi-esc/ice-earth-recon](https://github.com/awi-esc/ice-earth-recon) since kryos v0.1.0; its `tools/env.sh` sources
this `env.sh`.
