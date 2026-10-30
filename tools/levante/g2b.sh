#!/bin/bash
#SBATCH --partition=compute
#SBATCH --account=ba1442
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --exclusive
#SBATCH --job-name=g2b
# G2b: the mock-ice replay of PaleoMIST through the kryos globe and VILMA2 (tests/g2b_replay.f90, MODE=kryos) or the standalone vilma.x
# on the same history (MODE=vilma). sbatch --export=ALL,LMAX=85,NSTEP=2100,MODE=kryos,TAG=l085 g2b.sh
# Env: LMAX (170), NSTEP (2100), MODE (kryos|vilma), TAG, SPIN (1; 0: no spin-up, smoke), SHA (label of the harness), RATE_TOL (1e-5), TIME_MAX (2e4), NOUT (10),
# FESM (the fesmdata folder; g2b_visc_mirror.sh makes one with the 3D viscosity mirrored north-south). The diagnosis runs used
# -p shared, no --exclusive, --mem=24G, -t 01:10:00 (L85).
# Run folder $G/runs/g2b/kx-$SHA/S9/$TAG/{kryos,vilma}. Workspace $G/p4g2b: kryos (a copy of the S9 checkout, library as built there, plus the harness),
# vilma2 (a copy of the S9 clone, bin/vilma.x), data (g2b_prep.jl), maps (writable; the S6 maps linked, the Gauss maps of each lmax made once).
set -uo pipefail
G=/work/ba1442/robinson/kryos-gate; W=$G/p4g2b; S=$G/stacks/S9; K=$W/kryos
LMAX=${LMAX:-170}; NSTEP=${NSTEP:-2100}; MODE=${MODE:?MODE=kryos|vilma}; TAG=${TAG:-l$LMAX}; SPIN=${SPIN:-1}; SHA=${SHA:-wip}
OUT=$G/runs/g2b/kx-$SHA/S9/$TAG/$MODE
[ -e $OUT ] && { echo "exists: $OUT"; exit 1; }
mkdir -p $OUT/input; cd $OUT
for f in $S/yelmox/input/*; do ln -sf $f input/; done
for f in $K/input/*; do ln -sf $f input/; done
for d in ice_data isostasy_data; do ln -s /work/ba1442/$d $d; done
ln -s ${FESM:-/work/ba1442/fesmdata} fesmdata
ln -s $K/fesm-utils fesm-utils
export OMP_NUM_THREADS=64 OMP_PROC_BIND=spread OMP_PLACES=cores OMP_STACKSIZE=512M
ulimit -s unlimited
export G2B_DATA=$W/data G2B_SPINUP=$SPIN G2B_RATE_TOL=${RATE_TOL:-1.0e-5} G2B_TIME_MAX=${TIME_MAX:-2.0e4} G2B_NOUT=${NOUT:-10}
if [ $MODE = kryos ]; then
  ln -s $W/maps maps
  echo "== kryos g2b_replay lmax $LMAX nstep $NSTEP $(date -Is)"
  /usr/bin/time -v $K/libkryos/bin/g2b_replay.x $LMAX $NSTEP > run.log 2>&1; echo "g2b_replay exit $?"
  grep -E "^g2b:" run.log | grep -v " step "
  grep -E "Maximum resident|Elapsed" run.log
else
  mkdir maps
  G2B_NMLONLY=1 $K/libkryos/bin/g2b_replay.x $LMAX $NSTEP > nml.log 2>&1 || { cat nml.log; exit 1; }
  echo "== vilma.x lmax $LMAX nstep $NSTEP $(date -Is)"
  /usr/bin/time -v $W/vilma2/bin/vilma.x vilma_g2b.nml > run.log 2>&1; echo "vilma.x exit $?"
  grep -E "PROFILE setup|solid_earth_update  =|vilma_write_step|Maximum resident|Elapsed" run.log
fi
echo "== done $(date -Is)"
