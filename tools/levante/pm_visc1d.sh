#!/bin/bash
# The 1-D viscosity profile of VILMA2's pre_spinup_1d for EARTH=bm of pm_run.sh (the science Earth collapsed to 1-D):
#   tools/levante/pm_visc1d.sh build           pm_visc1d.f90 against libvilma.a of VILMA2, with the flags of its vilma.x
#                                              (make -n vilma; run in a login shell with the Intel modules)
#   tools/levante/pm_visc1d.sh <3-D run>       submits: pm_visc1d.x with the namelist of a 3-D run of pm_run.sh (e.g.
#                                              $PM_DIR/B3_l170) -> $PM_DIR/data/pm_visc1d_bagge6.csv, then pm_prep.jl visc1d
#                                              -> $PM_DIR/data/pm_visc_bm.nc (a laterally uniform FesmData-form file)
set -euo pipefail
_h=$(dirname "$(readlink -f "$0")"); [ -f "$_h/env.sh" ] || _h=$KRYOS/tools/levante; . "$_h/env.sh"
X=$PM_DIR/data/pm_visc1d.x
if [ "${1:-}" = build ]; then
  mkdir -p $PM_DIR/data
  cmd=$(make -s -C $VILMA2 -n vilma | sed -n '/src\/vilma_main.f90/{N;s/\\\n//;p}')
  cmd=${cmd/src\/vilma_main.f90/$_h/pm_visc1d.f90}; cmd=${cmd/-o bin\/vilma.x/-o $X}
  echo "$cmd"
  (cd $VILMA2 && eval "$cmd")
  echo "built $X"; exit 0
fi
if [ "${KRYOS_JOB:-}" != 1 ]; then
  [ $# -ge 1 ] || { sed -n '2,7p' "$0"; exit 1; }
  [ -x $X ] || { echo "no $X: tools/levante/pm_visc1d.sh build"; exit 1; }
  R=$(readlink -f "$1"); [ -f $R/pm.nml ] || { echo "no $R/pm.nml (a 3-D run of pm_run.sh)"; exit 1; }
  sbatch --parsable -A $ACCOUNT -p shared -t 00:30:00 -c 16 --mem=40G -J pm-visc1d -o $PM_DIR/data/pm_visc1d-%j.out \
      --export=ALL,KRYOS_JOB=1 "$_h/pm_visc1d.sh" "$R"
  exit
fi
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd "$1"
$X pm.nml $PM_DIR/data/pm_visc1d_bagge6.csv
kryos_jl $_h/pm_prep.jl visc1d $PM_DIR/data/pm_visc1d_bagge6.csv $PM_DIR/data/pm_visc_bm.nc
