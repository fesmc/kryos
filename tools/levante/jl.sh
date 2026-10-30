#!/bin/bash
# Run a Julia analysis script of this tools project as a small job on the shared partition (WP7c reports).
#   jl.sh <out file> <script.jl> [args ...]       (the project: the tools of the p4w7 checkout (FesmUtils.jl resolves there, as cmp_g4.sh); the depot of the gate)
G=/work/ba1442/robinson/kryos-gate
OUT=$1; shift
sbatch --parsable -A ba1442 -p shared -t 01:00:00 -c 4 --mem=16G -J jl-$(basename $1 .jl) -o $OUT --wrap \
  "export PATH=\$HOME/.juliaup/bin:\$HOME/.local/bin:\$PATH JULIA_DEPOT_PATH=$G/julia-depot; julia +1.13 --project=$G/p4w7/kryos/tools $*"
