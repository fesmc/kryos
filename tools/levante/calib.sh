#!/bin/bash -l
# Stage and submit a calib-cf run (docs/dev/calib-cf.md) on stack S9: kryos at $S/kryos-v0 (vilma=0, FastIsostasy), the par of par/calib/ copied to $G/calib/par.
#   [AFTER=<jobid>] [THREADS=64] [JOBTIME=03:00:00] [PART=compute] [TAG=] calib.sh <Antarctica|North> <pd|lgm> [set_nml args ...]
# lgm restarts from the end of the pd run of the same TAG (restart-$RSTT, default 30000; kryos.restart), afterok AFTER when given.
# Run folder $G/runs/calib/kx-<sha>/S9/<domain>-<case><TAG>. Exclusive node of the compute partition (PART), THREADS threads (the shared queue had a wait of hours).
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate; S=$G/stacks/S9; K=$S/kryos-v0; STACK=S9
DOM=$1; CASE=$2; shift 2
TAG=${TAG:-}; THREADS=${THREADS:-64}; JOBTIME=${JOBTIME:-03:00:00}; PART=${PART:-compute}; MAPS=${MAPS:-$G/maps/gen/calib}
SHA=$(git -C $K rev-parse --short=7 HEAD)
ROOT=$G/runs/calib/kx-$SHA/$STACK; OUT=$ROOT/$DOM-$CASE$TAG
[ -x $K/libkryos/bin/kryos.x ] || { echo "no kryos.x"; exit 1; }
[ -e $OUT ] && { echo "exists: $OUT"; exit 1; }
OV=("$@")
[ $CASE = lgm ] && OV+=("kryos.restart=\"$ROOT/$DOM-pd$TAG/restart-${RSTT:-30000}\"")
mkdir -p $OUT/input; cd $OUT
for f in $G/stacks/$STACK/yelmox/input/*; do ln -sf $f input/; done
for f in $K/input/*; do ln -sf $f input/; done
ln -sf $K/fesm-utils/par/datasrc_defaults.nml input/datasrc_defaults.nml
printf '# Time\tRSL\n-1.0e6\t-120.0\n1.0e6\t-120.0\n' > input/sl_const_lgm120.dat
for d in fesmdata ice_data isostasy_data; do ln -s /work/ba1442/$d $d; done
ln -s $K/fesm-utils fesm-utils; ln -s $MAPS maps
cp $G/calib/par/$DOM-$CASE.nml par.nml; cp $K/libkryos/bin/kryos.x .
[ ${#OV[@]} -gt 0 ] && python3 $G/tools/set_nml.py par.nml "${OV[@]}"
git -C $K rev-parse HEAD > kryos_commit.txt
cat > job.submit <<JOB
#!/bin/bash
#SBATCH --qos=normal
#SBATCH --partition=$PART
#SBATCH --time=$JOBTIME
#SBATCH --job-name=cal-$DOM-$CASE$TAG
#SBATCH --account=ba1442
#SBATCH --output=./out.out
#SBATCH --error=./out.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$THREADS
#SBATCH --exclusive
#SBATCH --no-requeue
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OMP_NUM_THREADS=$THREADS
export OMP_MAX_ACTIVE_LEVELS=1
export OMP_STACKSIZE=512M
ulimit -s unlimited
echo "node \$(hostname) \$(date -Is)"
/usr/bin/time -v -o time.txt srun ./kryos.x par.nml
rc=\$?
echo "end \$(date -Is) exit \$rc"
exit \$rc
JOB
echo "$DOM-$CASE$TAG $SHA $(sbatch --parsable ${AFTER:+--dependency=afterok:$AFTER} job.submit) $OUT"
