#!/bin/bash -l
# Register the kryos self-baseline K3-sci (regression-design 3.5, phase3-plan 3.1): the output files and restart bundles of the
# two-domain runs `cont` of gates sci and sci-const at 1000 yr, read-only, at the kryos commit of the sign-off.
#   [G=$GATE_ROOT] register_k3.sh <kryos sha> [stack=S5b]
# Refuses unless $G/reports/KK/summary.txt has a PASS line for both gates at <sha> and <stack>, and when the baseline exists.
# Result: $G/ref/K3-sci/kx-<sha>/<stack>/{sci,sci-const}/cont/ (nc files, bundles, par.nml, kryos_commit.txt; not the links to input
# and data, not the executable, not the logs), and registered.txt beside them. Later phases compare against it with
# kryos-compare.jl --map identity --restarts (phase 4 G1).
set -euo pipefail
G=${G:-/work/ba1442/robinson/kryos-gate}
SHA=$1; STACK=${2:-S5b}
SUM=$G/reports/KK/summary.txt
for g in sci sci-const; do
    grep -qE "^$g \| $SHA \| $STACK \| .*: PASS" $SUM 2>/dev/null || { echo "no PASS of gate $g at $SHA on $STACK in $SUM: not registering"; exit 1; }
    [ -d $G/runs/$g/kx-$SHA/$STACK/cont ] || { echo "no run $G/runs/$g/kx-$SHA/$STACK/cont"; exit 1; }
done
DEST=$G/ref/K3-sci/kx-$SHA/$STACK
[ -e $DEST ] && { echo "exists: $DEST"; exit 1; }
for g in sci sci-const; do
    SRC=$G/runs/$g/kx-$SHA/$STACK/cont
    mkdir -p $DEST/$g/cont
    tar -C $SRC -cf - --exclude=input --exclude=maps --exclude=fesmdata --exclude=ice_data --exclude=isostasy_data \
        --exclude=fesm-utils --exclude=kryos.x --exclude='*.log' --exclude=job.submit . | tar -C $DEST/$g/cont -xf -
    n=$(find $DEST/$g/cont -name '*.nc' | wc -l)
    [ "$n" -gt 0 ] || { echo "no nc files copied from $SRC"; exit 1; }
    echo "$g: $n nc files"
done
{ echo "K3-sci: kryos $SHA, stack $STACK, registered $(date -Is) by $(id -un)"; echo "gates sci, sci-const passed (reports/KK/summary.txt)"; } > $DEST/registered.txt
chmod -R a-w $DEST
echo "registered $DEST"
