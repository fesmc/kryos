#!/bin/bash -l
# Phase 4 WP7 runs on stack S7 (setup_S7.sh), kryos at the S7 checkouts. Prints the job ids. KSHA: the commit the checkouts are at.
#   gate_S7.sh g1    G1: K(S7, vilma=1) cont for grl, ant, ant-full, nor (maps S5 files) and the two-domain sci, sci-const cont
# Compare with compare_S7.sh.
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate; S=$G/stacks/S7; P=$G/p4w7
KSHA=$(git -C $S/kryos rev-parse --short=7 HEAD)
RK=${RK:-$P/levante/run_kryos.sh}; RKK=${RKK:-$P/levante/run_kryos_kk.sh}
M=${MAPS:-$G/maps/S5}
case ${1:-} in
  g1)
    for g in grl ant ant-full nor; do
        K=$S/kryos KSHA=$KSHA VILMA=1 STACK=S7 MAPS=$M $RK $g cont
    done
    for g in sci sci-const; do
        K=$S/kryos KSHA=$KSHA VILMA=1 STACK=S7 MAPS=$M TOOLS=$P/kryos $RKK $g cont
    done ;;
  *) echo "usage: gate_S7.sh g1"; exit 1 ;;
esac
