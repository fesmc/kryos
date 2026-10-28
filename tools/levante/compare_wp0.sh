#!/bin/bash -l
# WP0 comparisons (phase2-plan §3.1) with kryos-compare (identity mode), reports in $G/reports/WP0.
#   compare_wp0.sh ident            S4 identity: R1(S4) vs R0, R2(S4) vs R2, K(S4) vs K(S3) (grl, ant, cont, --restarts)
#   compare_wp0.sh base <gate>...   per gate: R1 vs R0, cont2 vs cont, chunk vs cont (R0, R1, R2'), R1 chunk vs R0 chunk,
#                                   R2' vs R0 (cont; differences expected)
# One line per comparison in $G/reports/WP0/summary.txt: name, outcome, compared/failing counts.
set -uo pipefail
module load python3 >/dev/null 2>&1 || true
export PATH=$HOME/.juliaup/bin:$HOME/.local/bin:$PATH
G=/work/ba1442/robinson/kryos-gate; K=$G/kryos
export JULIA_DEPOT_PATH=$G/julia-depot
T=${TOOLS:-$K/tools}                      # KryosTools checkout (project + bin/kryos-compare.jl)
GATES=${GATEDIR:-$T/gate/wp0}             # gate files with chunks / restart_times
OUT=$G/reports/WP0; mkdir -p $OUT
kc () {   # kc <name> <ref> <run> <opts...>
    local name=$1 ref=$2 run=$3; shift 3
    if [ ! -d "$ref" ] || [ ! -d "$run" ]; then echo "$name | MISSING | $ref $run" >> $OUT/summary.txt; return; fi
    julia +1.13 --project=$T $T/bin/kryos-compare.jl "$ref" "$run" --md $OUT/$name.md "$@" > $OUT/$name.log 2>&1
    echo "$name | $(sed -n 1p $OUT/$name.md | sed 's/# kryos-compare: //') | $(grep -m1 '^- compared:' $OUT/$name.md | sed 's/- compared: //')" >> $OUT/summary.txt
}
case $1 in
ident)
    for g in grl ant; do
        R0=$G/ref/$g/yx-918d996/v24; R1=$G/ref/$g/yx-918d996/S4
        R2o=$G/ref/$g/yx-74f3a2b/S3; R2n=$G/ref/$g/yx-5ac89c8/S4
        kc $g-R1S4-vs-R0 $R0/cont $R1/cont --restarts
        kc $g-R2S4-vs-R2 $R2o/cont $R2n/cont --restarts
        kc $g-KS4-vs-KS3 $G/runs/$g/kx-af8f1b9/S3/cont $G/runs/$g/kx-af8f1b9/S4/cont --restarts
    done ;;
base)
    shift
    for g in "$@"; do
        GT=$GATES/$g.toml
        R0=$G/ref/$g/yx-918d996/v24; R1=$G/ref/$g/yx-918d996/S4; R2=$G/ref/$g/yx-5ac89c8/S4
        kc $g-R1-vs-R0-cont       $R0/cont  $R1/cont  --restarts
        kc $g-R0-cont2-vs-cont    $R0/cont  $R0/cont2 --restarts
        kc $g-R1-cont2-vs-cont    $R1/cont  $R1/cont2 --restarts
        kc $g-R2-cont2-vs-cont    $R2/cont  $R2/cont2 --restarts
        kc $g-R0-chunk-vs-cont    $R0/cont  $R0/chunk --chunks --gate $GT
        kc $g-R1-chunk-vs-cont    $R1/cont  $R1/chunk --chunks --gate $GT
        kc $g-R2-chunk-vs-cont    $R2/cont  $R2/chunk --chunks --gate $GT
        kc $g-R1-vs-R0-c00        $R0/chunk/c00 $R1/chunk/c00 --restarts
        kc $g-R1-vs-R0-c01        $R0/chunk/c01 $R1/chunk/c01 --restarts
        kc $g-R2-vs-R0-cont       $R0/cont  $R2/cont
    done ;;
esac
echo "== $(date -Is) done $*" >> $OUT/summary.txt
