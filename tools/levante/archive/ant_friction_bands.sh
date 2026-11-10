#!/bin/bash
#SBATCH --account=ba1442 --partition=shared --time=00:30:00 --mem=8G --cpus-per-task=2 --job-name=jl-af-bands -o /work/ba1442/robinson/kryos-gate/calib/out/af-bands.%j.log
# ant-friction bands (docs/dev/ant-friction.md): sbatch ant_friction_bands.sh <out.md> <run name> ...  (run names in the calib run root; script copied to $G/calib)
G=/work/ba1442/robinson/kryos-gate; R=$G/runs/calib/kx-624f686/S9; export JULIA_DEPOT_PATH=$G/julia-depot
out=$1; shift
$HOME/.juliaup/bin/julia +1.13 --startup-file=no --project=$G/julia-plot $G/calib/ant_friction_bands.jl $out "${@/#/$R/}"
