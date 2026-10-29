#!/bin/bash
#SBATCH --partition=compute
#SBATCH --account=ba1442
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --exclusive
#SBATCH --job-name=g2a
# G2a: tests/g2a_cap.f90 (env G2A_SPINUP=0: no spin-up in both runs) (the VILMA2 adapter on a land-ice cap, lmax 170, the science-case &vilma) and vilma.x on the same load,
# then g2a_compare.jl. sbatch --export=ALL,NSTEP=6,TAG=cap g2a.sh. Run folder $G/runs/g2a/kx-<sha>/S6/<TAG>.
set -uo pipefail
module load python3 >/dev/null 2>&1 || true
export PATH=$HOME/.juliaup/bin:$HOME/.local/bin:$PATH JULIA_DEPOT_PATH=/work/ba1442/robinson/kryos-gate/julia-depot
G=/work/ba1442/robinson/kryos-gate; S=$G/stacks/S6; K=$S/kryos
SHA=$(git -C $K rev-parse --short=7 HEAD); OUT=$G/runs/g2a/kx-$SHA/S6/${TAG:-cap}
[ -e $OUT ] && { echo "exists: $OUT"; exit 1; }
mkdir -p $OUT/maps $OUT/input; cd $OUT
for f in $S/yelmox/input/*; do ln -sf $f input/; done
for f in $K/input/*; do ln -sf $f input/; done
for d in fesmdata ice_data isostasy_data; do ln -s /work/ba1442/$d $d; done
ln -s $K/fesm-utils fesm-utils
cp $K/libkryos/bin/g2a_cap.x .; cp $S/vilma2/bin/vilma.x .
export OMP_NUM_THREADS=64 OMP_PROC_BIND=spread OMP_PLACES=cores OMP_STACKSIZE=512M
ulimit -s unlimited
echo "== adapter $(date -Is)"
/usr/bin/time -v ./g2a_cap.x ${NSTEP:-6} > adapter.log 2>&1; echo "adapter exit $?"
grep "^g2a:" adapter.log
grep -E "Maximum resident|Elapsed" adapter.log
echo "== standalone $(date -Is)"
/usr/bin/time -v ./vilma.x vilma_g2a.nml > vilma.log 2>&1; echo "vilma.x exit $?"
grep -E "PROFILE setup|solid_earth_update  =|Maximum resident|Elapsed" vilma.log
echo "== compare $(date -Is)"
julia +1.13 --project=$G/p4w7/kryos/tools $G/p4w7/levante/g2a_compare.jl $OUT --md $OUT/g2a.md
