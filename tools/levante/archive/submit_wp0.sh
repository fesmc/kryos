#!/bin/bash -l
# Submit the phase 2 WP0 baselines (phase2-plan §3.1), per gate: cont, cont2, chunk.
#   submit_wp0.sh R0|R1|R2 "<gates>" ["<runs>"]
#     R0: yelmox v2.4 on stack v24 (default exe of the stack)
#     R1: yelmox v2.4 (918d996) on S4, exe $G/stacks/S4/bin/R1
#     R2: yelmox topo-prep (5ac89c8, R2') on S4, exe $G/stacks/S4/bin/R2
# The yelmox work tree of the stack is put at the commit of the baseline (par files); logs the job ids
# in $G/reports/jobs-wp0.txt as "<baseline> <gate> <run> <jobids>".
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate
B=$1; GATES=$2; RUNS=${3:-"cont cont2 chunk"}
case $B in
    R0) STACK=v24; COMMIT=918d996; unset EXE ;;
    R1) STACK=S4;  COMMIT=918d996; export EXE=$G/stacks/S4/bin/R1/yelmox.x ;;
    R2) STACK=S4;  COMMIT=5ac89c8; export EXE=$G/stacks/S4/bin/R2/yelmox.x ;;
    *) echo "unknown baseline $B" >&2; exit 1 ;;
esac
git -C $G/stacks/$STACK/yelmox checkout -q $COMMIT
for g in $GATES; do
    # job size per gate (shared partition: 940 MB per CPU)
    WALL=01:00:00; NCPU=
    case $g in
        grl|ant) WALL=00:30:00 ;;
        pyr|srg|egrl) WALL=04:00:00; NCPU=32 ;;
        eant) WALL=04:00:00; NCPU=64 ;;
    esac
    export WALL NCPU
    for r in $RUNS; do
        out=$($G/tools/submit_ref_wp0.sh $STACK $g $r) || { echo "FAILED $B $g $r" | tee -a $G/reports/jobs-wp0.txt; continue; }
        echo "$B $g $out" | tee -a $G/reports/jobs-wp0.txt
    done
done
