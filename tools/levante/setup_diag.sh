#!/bin/bash
#SBATCH --account=ba1442 --partition=shared --time=00:30:00 --mem=16G --cpus-per-task=4 --job-name=jl-setup-diag -o /work/ba1442/robinson/kryos-gate/setup-diag/diag.%j.log
# setup-diag (docs/dev/setup-diagnosis.md): sbatch setup_diag.sh -> tools/levante/setup_diag.jl (copied to $G/setup-diag), tables and maps to $G/setup-diag/out.
G=/work/ba1442/robinson/kryos-gate; export JULIA_DEPOT_PATH=$G/julia-depot
$HOME/.juliaup/bin/julia +1.13 --startup-file=no --project=$G/julia-plot $G/setup-diag/setup_diag.jl $G/setup-diag/out
