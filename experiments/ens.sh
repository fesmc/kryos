#!/bin/bash
# The weighted ensembles of tau_b (docs/simulations.md: E1 ens-v1, E2 ens-v2; docs/dev/plastic-ensemble-design.md): the
# members, their runs and the weights, from the config in tools/levante (plastic_ensemble_v<k>.toml).
#
#   experiments/ens.sh <v1|v2> [<folder>]      default folder $KRYOS_WORK/ens-<v>; it must not exist
#
# One init job draws the members (seeded) and submits them as a runme job array, then the weights after them
# (tools/levante/plastic_ensemble.sh). ens-v2: 50 members, about 3 h each on albedo; m000 is the prior mean (D2). Needs
# libkryos/bin/kryos.x built with plastic=1 vilma=1. The scores: experiments/scores.sh <folder>.
set -euo pipefail
X=$(cd "$(dirname "$(readlink -f "$0")")" && pwd); T=$X/../tools/levante; . "$T/env.sh"
[ $# -ge 1 ] || { sed -n '2,10p' "$0"; exit 1; }
v=$1; C=${2:-$KRYOS_WORK/ens-$v}
cfg=$T/plastic_ensemble_$v.toml; [ -f "$cfg" ] || { echo "no $cfg"; exit 1; }
bash "$T/plastic_ensemble.sh" start "$cfg" "$C"
