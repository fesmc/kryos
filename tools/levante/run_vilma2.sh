#!/bin/bash -l
# Stage and submit a run of the science case with geo = vilma2 on stack S6 (phase 4 WP7), kryos built with vilma=1.
#   [K=$S/kryos] [MAPS=$G/maps/S6] [PAR=par/kryos_NH-ANT_vilma2.nml (relative to K, or absolute)] [THREADS=64] [JOBTIME=02:00:00] [PART=compute] \
#   [CT0= CTM= CTE=] run_vilma2.sh <gate> cont|cont2|chunk|mapgen [set_nml args: group.key=value ...]
# cont, cont2: one run, mapgen: a run in a writable map folder $G/maps/gen/S6 (copy of $MAPS) that creates the missing maps;
# chunk: c00 CT0 -> CTM and c01 CTM -> CTE from c00/restart-<CTM> (afterok). The run folder is $G/runs/<gate>/kx-<sha>/S6/<run>.
# The node is exclusive (VILMA2 at lmax 170 needs 20-25 GB); OMP_PROC_BIND=spread. The extra args are set in the par file
# (existing keys only: clock.time_init, clock.time_end, geo_glob.out_gia, ...).
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate; S=$G/stacks/${STK:-S6}; K=${K:-$S/kryos}; STACK=${STK:-S6}
GATE=$1; RUN=$2; shift 2
MAPS=${MAPS:-$G/maps/S6}; PAR=${PAR:-par/kryos_NH-ANT_vilma2.nml}; THREADS=${THREADS:-64}; JOBTIME=${JOBTIME:-02:00:00}
PART=${PART:-compute}
[ -x $K/libkryos/bin/kryos.x ] || { echo "no $K/libkryos/bin/kryos.x"; exit 1; }
SHA=$(git -C $K rev-parse --short=7 HEAD)
ROOT=$G/runs/$GATE/kx-$SHA/$STACK
OUT=$ROOT/$RUN
[ -e $OUT ] && { echo "exists: $OUT"; exit 1; }
if [ "$RUN" = mapgen ]; then
    MAPS=$G/maps/gen/S6-$GATE-$SHA
    [ -e $MAPS ] && { echo "exists: $MAPS"; exit 1; }
    mkdir -p $MAPS; cp -p $G/maps/S5/* $MAPS/; chmod -R u+w $MAPS
    cp -p /work/ba1442/fesmdata/Global/GLOBAL-0.5DEG/grid_GLOBAL-0.5DEG.txt $MAPS/    # the background grid of [globe_ice] (FesmData)
    [ -d $G/maps/S6 ] && { cp -pn $G/maps/S6/* $MAPS/ 2>/dev/null || true; }
fi

stage () {   # stage <dir> <jobname> <set_nml args...>
    local dir=$1 name=$2; shift 2
    mkdir -p $dir/input; cd $dir
    for f in $G/stacks/$STACK/yelmox/input/*; do ln -sf $f input/; done
    for f in $K/input/*; do ln -sf $f input/; done
    ln -sf $K/fesm-utils/par/datasrc_defaults.nml input/datasrc_defaults.nml
    for d in fesmdata ice_data isostasy_data; do ln -s /work/ba1442/$d $d; done
    ln -s $K/fesm-utils fesm-utils; ln -s $MAPS maps
    case $PAR in /*) cp $PAR par.nml ;; *) cp $K/$PAR par.nml ;; esac; cp $K/libkryos/bin/kryos.x .
    [ $# -gt 0 ] && python3 $G/tools/set_nml.py par.nml "$@"
    git -C $K rev-parse HEAD > kryos_commit.txt
    cat > job.submit <<JOB
#!/bin/bash
#SBATCH --qos=normal
#SBATCH --partition=$PART
#SBATCH --time=$JOBTIME
#SBATCH --job-name=$name
#SBATCH --account=ba1442
#SBATCH --output=./out.out
#SBATCH --error=./out.err
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=$THREADS
#SBATCH --exclusive
#SBATCH --no-requeue
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OMP_NUM_THREADS=$THREADS
export OMP_MAX_ACTIVE_LEVELS=1
export OMP_STACKSIZE=512M
export MKL_NUM_THREADS=$THREADS
ulimit -s unlimited
ulimit -v unlimited
echo "node \$(hostname) \$(date -Is)"
/usr/bin/time -v -o time.txt srun ./kryos.x par.nml
rc=\$?
echo "end \$(date -Is) exit \$rc"
exit \$rc
JOB
}
case $RUN in
  cont|cont2|mapgen|lossless)
    stage $OUT $GATE-$RUN "$@"
    echo "$GATE $RUN $SHA $(cd $OUT && sbatch --parsable job.submit) $OUT" ;;
  chunk)
    : ${CT0:?} ${CTM:?} ${CTE:?}
    stage $OUT/c00 $GATE-c00 clock.time_init=$CT0 clock.time_end=$CTM "$@"
    j0=$(cd $OUT/c00 && sbatch --parsable job.submit)
    RD=$(python3 -c "t=float('$CTM'); print(int(t) if t==int(t) else t)")
    stage $OUT/c01 $GATE-c01 clock.time_init=$CTM clock.time_end=$CTE "kryos.restart=\"$OUT/c00/restart-$RD\"" "$@"
    j1=$(cd $OUT/c01 && sbatch --parsable --dependency=afterok:$j0 job.submit)
    echo "$GATE chunk $SHA $j0 $j1 $OUT" ;;
  *) echo "unknown run $RUN"; exit 1 ;;
esac
