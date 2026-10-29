#!/bin/bash -l
# Phase 2 WP2-WP4 gates (phase2-plan): kryos K against R2 (yelmox topo-prep, stack S5), reports in $G/reports/P2.
#   compare_p2.sh <gate> <kryos sha> [stack=S5] [R2 root=$G/ref-S5]     gate: ant-cf | ant-full | nor
# K runs in $G/runs/<gate>/kx-<sha>/<stack>, R2 in <R2 root>/<gate>/yx-<sha>/<stack>/cont. Per gate: K cont and
# cont2 against R2 (mapped, outputs and restart bundles), cont2 against cont (identity) and K chunk against K cont
# (--chunks; yelmox restarts the window time per chunk, so R2 chunk is not compared). One line per comparison
# in summary.txt. The tools are those of the kryos checkout, which must be at <sha>.
set -uo pipefail
module load python3 >/dev/null 2>&1 || true
export PATH=$HOME/.juliaup/bin:$HOME/.local/bin:$PATH
G=/work/ba1442/robinson/kryos-gate; K=$G/kryos
export JULIA_DEPOT_PATH=$G/julia-depot
T=$K/tools
g=$1; SHA=$2; STACK=${3:-S5}; R2ROOT=${4:-$G/ref-S5}
OUT=$G/reports/P2; mkdir -p $OUT
kc () {   # kc <name> <ref> <run> <opts...>
    local name=$1 ref=$2 run=$3; shift 3
    if [ ! -d "$ref" ] || [ ! -d "$run" ]; then echo "$name | MISSING | $ref $run" >> $OUT/summary.txt; return; fi
    julia +1.13 --project=$T $T/bin/kryos-compare.jl "$ref" "$run" --md $OUT/$name.md "$@" > $OUT/$name.log 2>&1
    echo "$name | $(sed -n 1p $OUT/$name.md | sed 's/# kryos-compare: //') | $(grep -m1 '^- compared:' $OUT/$name.md | sed 's/- compared: //')" >> $OUT/summary.txt
}
R2=$(ls -d $R2ROOT/$g/yx-*/$STACK | head -1); KR=$G/runs/$g/kx-$SHA/$STACK
M="--map $T/gate/map_phase1.toml --deviations $T/gate/deviations.toml --gate $T/gate/$g.toml"
kc $g-K-cont-vs-R2    $R2/cont $KR/cont  $M --restarts
kc $g-K-cont2-vs-R2   $R2/cont $KR/cont2 $M --restarts
kc $g-K-cont2-vs-cont $KR/cont $KR/cont2 --restarts
kc $g-K-chunk-vs-cont $KR/cont $KR/chunk --chunks --gate $T/gate/$g.toml
echo "== $(date -Is) done $g $SHA $STACK" >> $OUT/summary.txt
