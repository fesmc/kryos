#!/bin/bash
# pm-ref: iterate the present-day reference of a standalone VILMA2 PaleoMIST run so that the modelled 0 ka bed equals
# today's (SELEN's present topography as final condition). Iteration n is a run of a case (tools/levante/paleomist/pm_<EARTH>.nml),
# staged by runme in the PaleoMIST project ROOT (pm_project.sh), with the reference correction r_n (r_0 = 0: i_eq = 1, then
# i_eq = 3 and rsl_restart_file = r_n), then r_{n+1} = rsl(0 ka) of that run (pm_prep.jl rslref), until
# max |r_n - rsl_n(0 ka)| < TOL m or NMAX runs. Submits itself (shared partition, 3 h per run x NMAX):
#   BASE=I910 TOL=1.0 NMAX=8 tools/levante/pm_iterref.sh
# Env: BASE (the run names <BASE>_it<n>), EARTH (1d), FORCE (B), LMAX (170), TMAX (spin-up cap, 2.0e6), RHOI (910.0),
#      PM_P (more runme -p settings, e.g. "vilma.rotation_k_s=0.94349"), NTHR (32), TOL (1.0), NMAX (8), ROOT ($PM_DIR/ref).
# Runs and the log <BASE>_iter.log in ROOT.
set -uo pipefail
_h=$(dirname "$(readlink -f "$0")"); [ -f "$_h/env.sh" ] || _h=$KRYOS/tools/levante; . "$_h/env.sh"
BASE=${BASE:?}; TOL=${TOL:-1.0}; NMAX=${NMAX:-8}; ROOT=${ROOT:-$PM_DIR/ref}
EARTH=${EARTH:-1d}; FORCE=${FORCE:-B}; LMAX=${LMAX:-170}; TMAX=${TMAX:-2.0e6}; RHOI=${RHOI:-910.0}; NTHR=${NTHR:-32}; PM_P=${PM_P:-}
CASE=$_h/paleomist/pm_$EARTH.nml
[ -f "$CASE" ] || { echo "no case $CASE (EARTH: 1d | l80 | bm | 3d | t3)"; exit 1; }
if [ "${KRYOS_JOB:-}" != 1 ]; then
    bash "$_h/pm_project.sh" "$ROOT" || exit 1
    sbatch --parsable -A $ACCOUNT $KRYOS_PART -t $(( 3 * NMAX )):00:00 -c $NTHR --mem=16G -J pm-iterref-$BASE \
        -o $ROOT/${BASE}_iterref-%j.out --export=ALL,KRYOS_JOB=1 "$_h/pm_iterref.sh"
    exit
fi
export OMP_NUM_THREADS=$NTHR OMP_PROC_BIND=true OMP_STACKSIZE=512M OMP_MAX_ACTIVE_LEVELS=1
ulimit -s unlimited
prev=""
for n in $(seq 0 $((NMAX - 1))); do
  R=${BASE}_it$n
  ref=""; [ -n "$prev" ] && ref="ctl.i_eq=3 ctl.rsl_restart_file=$prev"
  (cd $ROOT && $KRYOS_RUNME -e vilma -n $CASE -o $R -p ctl.file_forcing=pm_data/pm_force_$FORCE.nc vilma.lmax=$LMAX \
       vilma.equil_time_max=$TMAX vilma.rho_ice=$RHOI $ref $PM_P > /dev/null) || { echo "runme failed for $R"; exit 1; }
  echo "== $R $(date -Is) on $(hostname)" | tee -a $ROOT/${BASE}_iter.log
  (cd $ROOT/$R && t0=$(date +%s) && ./vilma.x $(basename $CASE) > out.out 2> out.err; echo "$? $(( $(date +%s) - t0 ))" > status)
  read rc sec < $ROOT/$R/status; [ "$rc" = 0 ] || { echo "$R failed (exit $rc): $ROOT/$R/out.err" | tee -a $ROOT/${BASE}_iter.log; exit 1; }
  r=$ROOT/${BASE}_rslref_$((n + 1)).nc
  kryos_jl $_h/pm_prep.jl rslref $ROOT/$R/vilma_out.nc $r $prev | tee -a $ROOT/${BASE}_iter.log
  m=$(tail -1 $ROOT/${BASE}_iter.log | sed 's/.*max |d| = \([0-9.]*\).*/\1/')
  [ -n "$prev" ] && awk -v m=$m -v t=$TOL 'BEGIN{exit !(m < t)}' && { echo "converged after $((n + 1)) runs"; break; }
  prev=$r
done
