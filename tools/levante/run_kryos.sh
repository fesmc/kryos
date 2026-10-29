#!/bin/bash -l
# Update kryos on Levante, build kryos.x with the gate flags, stage and submit a kryos gate run.
# The job template is that of the base gate (ant-cf, ant-full: ant).
#   [STACK=S5] [SRC=<git url>] [BUNDLE=kryos.bundle] [REF=main] [MAPS=$G/maps/S] run_kryos.sh <gate> cont|cont2|chunk|base   (STACK default S3)
# cont2: a second cont run (determinism). The code is the REF branch of SRC (default $G/bundles/$BUNDLE;
# SRC=git@github.com:fesmc/kryos.git pulls from origin). MAPS is the folder linked as <rundir>/maps.
# chunk: c00 0 -> 500, c01 500 -> 1000 from c00/restart-500 (afterok), as submit_ref.sh.
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate; K=$G/kryos
GATE=$1; RUN=$2; STACK=${STACK:-S3}   # stack the kryos dependency links point at
BUNDLE=${BUNDLE:-kryos.bundle}; REF=${REF:-main}; SRC=${SRC:-$G/bundles/$BUNDLE}; MAPS=${MAPS:-$G/maps/S}
GF="-O2 -fp-model precise -march=znver3 -traceback"
CM="env PYTHONPATH=$G/tools/configme-src/src python3 -m configme"
cd $K
git fetch -q $SRC $REF:refs/remotes/bundle/main
git checkout -q --detach bundle/main
SHA=$(git rev-parse --short=7 HEAD)
$CM config kryos --only -m dkrz_levante -c ifx > /dev/null; $CM deps . > /dev/null
make kryos openmp=1 DFLAGS_NODEBUG="$GF" > $G/stacks/kryos-build-$SHA.log 2>&1 || { tail -30 $G/stacks/kryos-build-$SHA.log; exit 1; }
OUT=$G/runs/$GATE/kx-$SHA/$STACK/$RUN
[ -e $OUT ] && { echo "exists: $OUT"; exit 1; }

stage () {   # stage <dir> <jobname> <set_nml args...>; echo the folder
    local dir=$1 name=$2; shift 2
    mkdir -p $dir/input; cd $dir
    for f in $G/stacks/$STACK/yelmox/input/*; do ln -sf $f input/; done
    for f in $K/input/*; do ln -sf $f input/; done
    ln -sf $K/fesm-utils/par/datasrc_defaults.nml input/datasrc_defaults.nml
    for d in fesmdata ice_data isostasy_data; do ln -s /work/ba1442/$d $d; done
    ln -s $K/fesm-utils fesm-utils; ln -s $MAPS maps
    cp $K/par/gate/$GATE.nml par.nml; cp $K/libkryos/bin/kryos.x .
    [ $# -gt 0 ] && python3 $G/tools/set_nml.py par.nml "$@"
    git -C $K rev-parse HEAD > kryos_commit.txt
    sed -e "s/--job-name=.*/--job-name=$name/" -e 's|^srun .*|srun ./kryos.x par.nml|' \
        $G/ref/${GATE%%-*}/yx-74f3a2b/S2/cont/job.submit > job.submit
}
case $RUN in
  cont|cont2)
    stage $OUT $GATE-kx-$RUN
    echo "$GATE $RUN $SHA $(cd $OUT && sbatch --parsable job.submit) $OUT" ;;
  base)                                            # 15 kyr (regression-design 1.1), 2 h limit
    stage $OUT $GATE-kx-base clock.time_end=15000.0
    sed -i 's/--time=.*/--time=02:00:00/' $OUT/job.submit
    echo "$GATE base $SHA $(cd $OUT && sbatch --parsable job.submit) $OUT" ;;
  chunk)
    stage $OUT/c00 $GATE-kx-c00 clock.time_init=0.0 clock.time_end=500.0
    j0=$(cd $OUT/c00 && sbatch --parsable job.submit)
    stage $OUT/c01 $GATE-kx-c01 clock.time_init=500.0 clock.time_end=1000.0 "kryos.restart=\"$OUT/c00/restart-500\""
    j1=$(cd $OUT/c01 && sbatch --parsable --dependency=afterok:$j0 job.submit)
    echo "$GATE chunk $SHA $j0 $j1 $OUT" ;;
esac
