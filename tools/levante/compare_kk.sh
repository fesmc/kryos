#!/bin/bash -l
# The checks of a kryos-against-kryos gate (tools/gate/<gate>.toml, kind = "kryos") over its runs; report in $G/reports/KK.
#   [TOOLS=<checkout with the Julia project>] compare_kk.sh <gate> <kryos sha> [stack=S5b]
# The runs are $G/runs/<gate>/kx-<sha>/<stack>/<run>; one line per gate in summary.txt.
set -uo pipefail
module load python3 >/dev/null 2>&1 || true
export PATH=$HOME/.juliaup/bin:$HOME/.local/bin:$PATH
G=/work/ba1442/robinson/kryos-gate; TOOLS=${TOOLS:-$G/kryos}
export JULIA_DEPOT_PATH=$G/julia-depot
g=$1; SHA=$2; STACK=${3:-S5b}
OUT=$G/reports/KK; mkdir -p $OUT
RUNS=$G/runs/$g/kx-$SHA/$STACK
julia +1.13 --project=$TOOLS/tools $TOOLS/tools/bin/kryos-gate-kk.jl check $g $RUNS --md $OUT/$g-$SHA-$STACK.md > $OUT/$g-$SHA-$STACK.log 2>&1
rc=$?
echo "$g | $SHA | $STACK | $(sed -n 1p $OUT/$g-$SHA-$STACK.md) | exit $rc" >> $OUT/summary.txt
