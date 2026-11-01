#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=ba1442
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --job-name=pm_iterref
# pm-ref: iterate the present-day reference so that the modelled 0 ka bed equals today's (SELEN's present topography as final
# condition). Iteration n runs B1 (pm_run.sh, i_eq = 3) with the reference correction r_n (r_0 = 0: i_eq = 1), then r_{n+1} = rsl(0 ka)
# of that run (pm_prep.jl rslref), until max |r_n - rsl_n(0 ka)| < TOL m or NMAX runs.
#   sbatch --export=ALL,BASE=I910,TOL=1.0,NMAX=8 pm_iterref.sh      (other pm_run.sh settings pass through, default RHOI=910)
set -uo pipefail
G=/work/ba1442/robinson/kryos-gate; P=$G/runs/paleomist; T=$(dirname $(readlink -f $0)); [ -f $T/pm_run.sh ] || T=$P/tools_ref
export PATH=$HOME/.juliaup/bin:$PATH JULIA_DEPOT_PATH=$G/julia-depot
BASE=${BASE:?}; TOL=${TOL:-1.0}; NMAX=${NMAX:-8}; ROOT=${ROOT:-$P/ref}; mkdir -p $ROOT
export ROOT LMAX=${LMAX:-170} EARTH=${EARTH:-1d} FORCE=${FORCE:-B} TMAX=${TMAX:-2.0e6} RHOI=${RHOI:-910.0} NTHR=${NTHR:-32}
prev=""
for n in $(seq 0 $((NMAX - 1))); do
  RUN=${BASE}_it$n RSLREF=$prev bash $T/pm_run.sh
  r=$ROOT/${BASE}_rslref_$((n + 1)).nc
  julia +1.13 --project=$G/p4w7/kryos/tools $T/pm_prep.jl rslref $ROOT/${BASE}_it$n/vilma_out.nc $r $prev | tee -a $ROOT/${BASE}_iter.log
  m=$(tail -1 $ROOT/${BASE}_iter.log | sed 's/.*max |d| = \([0-9.]*\).*/\1/')
  [ -n "$prev" ] && awk -v m=$m -v t=$TOL 'BEGIN{exit !(m < t)}' && { echo "converged after $((n + 1)) runs"; break; }
  prev=$r
done
