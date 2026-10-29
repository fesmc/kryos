#!/bin/bash -l
# Phase 4 WP7 runs on stack S6 (setup_S6.sh), kryos cd85c14. Prints the job ids.
#   gate_S6.sh id    identity: K(S6, vilma=0) cont for grl, ant (maps S5) and R2'(S6) cont for grl, ant
#   gate_S6.sh g1    G1: K(S6, vilma=1) cont for grl, ant, ant-full, nor (maps S6) and the two-domain sci, sci-const cont
# Compare with compare_S6.sh. RK, RKK: the repo's run_kryos.sh, run_kryos_kk.sh.
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate; S=$G/stacks/S6; KSHA=cd85c14; P=$G/p4w7
RK=${RK:-$P/levante/run_kryos.sh}; RKK=${RKK:-$P/levante/run_kryos_kk.sh}
M=${MAPS:-$G/maps/S5}      # the maps of the fastiso configurations: the S5 files (maps/S6 = maps/S5 + the new pairs)
case ${1:-} in
  id)
    for g in grl ant; do
        K=$S/kryos-v0 KSHA=$KSHA VILMA=0 STACK=S6 RLABEL=S6-v0 MAPS=$M $RK $g cont
        OUTROOT=$G/ref-S6 $G/tools/submit_ref.sh S6 $g cont $G/maps/S5 | sed "s/^/R2 $g /"
    done ;;
  g1)
    for g in grl ant ant-full nor; do
        K=$S/kryos KSHA=$KSHA VILMA=1 STACK=S6 MAPS=$M $RK $g cont
    done
    for g in sci sci-const; do
        K=$S/kryos KSHA=$KSHA VILMA=1 STACK=S6 MAPS=$M TOOLS=$P/kryos $RKK $g cont
    done ;;
  *) echo "usage: gate_S6.sh id|g1"; exit 1 ;;
esac
