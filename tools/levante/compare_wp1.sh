#!/bin/bash -l
# WP1 gate (phase2-plan WP1): kryos K against R2 on grl and ant, reports in $G/reports/WP1.
#   compare_wp1.sh <kryos sha> [stack]     (K runs in $G/runs/<gate>/kx-<sha>/<stack>, stack default S4)
# Per gate: K cont and cont2 against R2 (yelmox 74f3a2b, stack S2; mapped, outputs and restart bundles),
# cont2 against cont and chunk against cont (identity). One line per comparison in summary.txt.
set -uo pipefail
module load python3 >/dev/null 2>&1 || true
export PATH=$HOME/.juliaup/bin:$HOME/.local/bin:$PATH
G=/work/ba1442/robinson/kryos-gate; K=$G/kryos
export JULIA_DEPOT_PATH=$G/julia-depot
T=$K/tools
SHA=$1; STACK=${2:-S4}
OUT=$G/reports/WP1; mkdir -p $OUT
kc () {   # kc <name> <ref> <run> <opts...>
    local name=$1 ref=$2 run=$3; shift 3
    if [ ! -d "$ref" ] || [ ! -d "$run" ]; then echo "$name | MISSING | $ref $run" >> $OUT/summary.txt; return; fi
    julia +1.13 --project=$T $T/bin/kryos-compare.jl "$ref" "$run" --md $OUT/$name.md "$@" > $OUT/$name.log 2>&1
    echo "$name | $(sed -n 1p $OUT/$name.md | sed 's/# kryos-compare: //') | $(grep -m1 '^- compared:' $OUT/$name.md | sed 's/- compared: //')" >> $OUT/summary.txt
}
for g in grl ant; do
    R2=$G/ref/$g/yx-74f3a2b/S2; KR=$G/runs/$g/kx-$SHA/$STACK
    M="--map $T/gate/map_phase1.toml --deviations $T/gate/deviations.toml --gate $T/gate/$g.toml"
    kc $g-K-cont-vs-R2    $R2/cont $KR/cont  $M --restarts
    kc $g-K-cont2-vs-R2   $R2/cont $KR/cont2 $M --restarts
    kc $g-K-cont2-vs-cont $KR/cont $KR/cont2 --restarts
    kc $g-K-chunk-vs-cont $KR/cont $KR/chunk --chunks --gate $T/gate/$g.toml
done
echo "== $(date -Is) done $SHA $STACK" >> $OUT/summary.txt
