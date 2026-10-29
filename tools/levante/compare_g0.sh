#!/bin/bash -l
# Phase 4 WP0 gate G0 (phase4-plan §13.2): S5 identity against S4, grl and ant, reports in $G/reports/P4WP0.
#   compare_g0.sh <kryos sha S5> [<kryos sha S4>=f81e004] [<R2 root S5>=$G/ref-S5]
# K(S5, vilma=0) cont vs K(S4) cont (code of main == f81e004 apart from docs), and R2(S5) cont vs R2(S4) cont
# (yelmox topo-prep 5ac89c8, rebuilt in S5); outputs and restart bundles (--restarts). One line per comparison in summary.txt.
set -uo pipefail
module load python3 >/dev/null 2>&1 || true
export PATH=$HOME/.juliaup/bin:$HOME/.local/bin:$PATH
G=/work/ba1442/robinson/kryos-gate; K=$G/kryos
export JULIA_DEPOT_PATH=$G/julia-depot
T=$K/tools
SHA5=$1; SHA4=${2:-f81e004}; R2ROOT=${3:-$G/ref-S5}
OUT=$G/reports/P4WP0; mkdir -p $OUT
kc () {   # kc <name> <ref> <run> <opts...>
    local name=$1 ref=$2 run=$3; shift 3
    if [ ! -d "$ref" ] || [ ! -d "$run" ]; then echo "$name | MISSING | $ref $run" >> $OUT/summary.txt; return; fi
    julia +1.13 --project=$T $T/bin/kryos-compare.jl "$ref" "$run" --md $OUT/$name.md "$@" > $OUT/$name.log 2>&1
    echo "$name | $(sed -n 1p $OUT/$name.md | sed 's/# kryos-compare: //') | $(grep -m1 '^- compared:' $OUT/$name.md | sed 's/- compared: //')" >> $OUT/summary.txt
}
for g in grl ant; do
    kc $g-KS5-vs-KS4  $G/runs/$g/kx-$SHA4/S4/cont  $G/runs/$g/kx-$SHA5/S5/cont  --restarts
    kc $g-R2S5-vs-R2S4 $G/ref/$g/yx-5ac89c8/S4/cont $R2ROOT/$g/yx-5ac89c8/S5/cont --restarts
done
echo "== $(date -Is) done K(S5) $SHA5 vs K(S4) $SHA4" >> $OUT/summary.txt
