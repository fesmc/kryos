#!/bin/bash -l
# S5b identity gate (setup_S5b.sh): the cont runs of grl and ant with K (kryos 5a663a3, vilma=0, built on S5b) and with R2'
# (yelmox topo-prep, S5b build), staged and submitted as the S5 runs were (maps S5). Prints the job ids.
#   gate_S5b.sh            then: sbatch compare_S5b.sh (after the jobs)
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate; S=$G/stacks/S5b; KSHA=5a663a3
RK=${RK:-$G/tools/run_kryos.sh}   # the repo's tools/levante/run_kryos.sh (K, KSHA overrides)
for g in grl ant; do
    K=$S/kryos KSHA=$KSHA STACK=S5b MAPS=$G/maps/S5 $RK $g cont
    OUTROOT=$G/ref-S5b $G/tools/submit_ref.sh S5b $g cont $G/maps/S5 | sed "s/^/R2 $g /"
done
