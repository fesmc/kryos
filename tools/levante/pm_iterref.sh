#!/bin/bash
# pm-ref: iterate the present-day reference of a standalone VILMA2 PaleoMIST run so that the modelled 0 ka bed equals
# today's (SELEN's present topography as final condition). Iteration n runs pm_run.sh (i_eq = 3) with the reference
# correction r_n (r_0 = 0: i_eq = 1), then r_{n+1} = rsl(0 ka) of that run (pm_prep.jl rslref), until
# max |r_n - rsl_n(0 ka)| < TOL m or NMAX runs. Submits itself (shared partition, 3 h per run x NMAX):
#   BASE=I910 TOL=1.0 NMAX=8 tools/levante/pm_iterref.sh      (other pm_run.sh settings pass through, default RHOI=910)
# Runs and the log in ROOT (default $PM_DIR/ref).
set -uo pipefail
_h=$(dirname "$(readlink -f "$0")"); [ -f "$_h/env.sh" ] || _h=$KRYOS/tools/levante; . "$_h/env.sh"
BASE=${BASE:?}; TOL=${TOL:-1.0}; NMAX=${NMAX:-8}; ROOT=${ROOT:-$PM_DIR/ref}; mkdir -p $ROOT
export ROOT LMAX=${LMAX:-170} EARTH=${EARTH:-1d} FORCE=${FORCE:-B} TMAX=${TMAX:-2.0e6} RHOI=${RHOI:-910.0} NTHR=${NTHR:-32}
if [ "${KRYOS_JOB:-}" != 1 ]; then
    sbatch --parsable -A $ACCOUNT -p shared -t $(( 3 * NMAX )):00:00 -c $NTHR --mem=16G -J pm-iterref-$BASE \
        -o $ROOT/${BASE}_iterref-%j.out --export=ALL,KRYOS_JOB=1 "$_h/pm_iterref.sh"
    exit
fi
prev=""
for n in $(seq 0 $((NMAX - 1))); do
  RUN=${BASE}_it$n RSLREF=$prev bash $_h/pm_run.sh
  r=$ROOT/${BASE}_rslref_$((n + 1)).nc
  kryos_jl $_h/pm_prep.jl rslref $ROOT/${BASE}_it$n/vilma_out.nc $r $prev | tee -a $ROOT/${BASE}_iter.log
  m=$(tail -1 $ROOT/${BASE}_iter.log | sed 's/.*max |d| = \([0-9.]*\).*/\1/')
  [ -n "$prev" ] && awk -v m=$m -v t=$TOL 'BEGIN{exit !(m < t)}' && { echo "converged after $((n + 1)) runs"; break; }
  prev=$r
done
