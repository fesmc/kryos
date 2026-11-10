#!/bin/bash -l
# Phase 4 S8 runs on stack S8 (setup_S8.sh), kryos at the S8 checkouts. Prints the job ids. KSHA: the commit the checkouts are at.
#   gate_S7.sh g1    G1: K(S7, vilma=1) cont for grl, ant, ant-full, nor (maps S5 files) and the two-domain sci, sci-const cont
#   gate_S8.sh g4    G4: the science run with geo = vilma2 (cont, cont2, chunk at -79800, lossless), 400 yr
# Compare with compare_S8.sh (g1) and cmp_g4.sh (SHA=$KSHA STK=S8).
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate; S=$G/stacks/S8; P=$G/p4w7
KSHA=$(git -C $S/kryos rev-parse --short=7 HEAD)
RK=${RK:-$P/levante/run_kryos.sh}; RKK=${RKK:-$P/levante/run_kryos_kk.sh}
M=${MAPS:-$G/maps/S5}
case ${1:-} in
  g1)
    for g in grl ant ant-full nor; do
        K=$S/kryos KSHA=$KSHA VILMA=1 STACK=S8 MAPS=$M $RK $g cont
    done
    for g in sci sci-const; do
        K=$S/kryos KSHA=$KSHA VILMA=1 STACK=S8 MAPS=$M TOOLS=$P/kryos $RKK $g cont
    done ;;
  g4)
    STK=S8 K=$S/kryos $P/levante/g34_vilma2.sh ;;
  *) echo "usage: gate_S8.sh g1|g4"; exit 1 ;;
esac
