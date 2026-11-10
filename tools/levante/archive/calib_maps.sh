#!/bin/bash
#SBATCH --account=ba1442 --partition=shared --time=00:30:00 --mem=16G --cpus-per-task=4 --job-name=jl-calib-maps -o /work/ba1442/robinson/kryos-gate/calib/maps/maps.%j.log
# calib-cf maps (docs/dev/calib-cf.md): sbatch calib_maps.sh -> tools/calib_maps.jl on every metrics output of $G/calib/out, PNGs to $G/calib/maps.
# Plot environment $G/julia-plot (CairoMakie, NCDatasets, DelimitedFiles), created by $G/julia-plot/install.sh (an sbatch job):
#   julia +1.13 --project=$G/julia-plot -e 'using Pkg; Pkg.add(["CairoMakie", "NCDatasets", "DelimitedFiles"]); Pkg.precompile()'
G=/work/ba1442/robinson/kryos-gate; export JULIA_DEPOT_PATH=$G/julia-depot
$HOME/.juliaup/bin/julia +1.13 --startup-file=no --project=$G/julia-plot $G/calib/calib_maps.jl $G/calib/out $G/calib/maps
