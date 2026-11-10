#!/bin/bash
# jl.sh <out file> <script.jl> [args ...]: a Julia script of tools/levante (its project, env.sh) as a small job on the
# shared partition (1 h, 4 cpus, 16 GB; JL_TIME, JL_MEM override). Prints the job id.
set -euo pipefail
_h=$(dirname "$(readlink -f "$0")"); . "$_h/env.sh"
[ $# -ge 2 ] || { sed -n '2,3p' "$0"; exit 1; }
OUT=$(readlink -f "$1"); shift
sbatch --parsable -A "$ACCOUNT" -p shared -t "${JL_TIME:-01:00:00}" -c 4 --mem="${JL_MEM:-16G}" -J "jl-$(basename "$1" .jl)" \
    -o "$OUT" --wrap "$KRYOS_JULIA --project=$_h $*"
