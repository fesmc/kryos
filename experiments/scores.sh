#!/bin/bash
# The scores of docs/simulations.md: every member of an ensemble and every reference of its config, against GAPSLIP 2.1
# RSL (age <= 25 ka), the ESL stack and GNSS uplift, under the ensemble's error model (tools/levante/plastic_ensemble.jl
# scores).
#
#   experiments/scores.sh <ensemble folder> [<wave>]      default wave 1; a job of about 10 min
#
# The references and the groups are those of <folder>/config.toml, copied from tools/levante/plastic_ensemble_v<k>.toml
# when the ensemble started; a reference whose run is not finished is left out and named. Output:
# <folder>/iter_<wave>/scores.csv, one row per series, group and window (columns of docs/simulations.md: Scores).
set -euo pipefail
X=$(cd "$(dirname "$(readlink -f "$0")")" && pwd); T=$X/../tools/levante; . "$T/env.sh"
[ $# -ge 1 ] || { sed -n '2,11p' "$0"; exit 1; }
C=$(realpath "$1"); k=${2:-1}
sbatch --parsable -A $ACCOUNT $KRYOS_PART -t 00:30:00 -c 8 --mem=32G -J ens-scores -o "$C/scores-$k.out" --wrap \
    "$KRYOS_JULIA --project=$T $T/plastic_ensemble.jl scores $C $k"
