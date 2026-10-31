#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=ba1442
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --job-name=pm-visc1d
# PaleoMIST test of the standalone VILMA2: the 1-D profile of pre_spinup_1d for EARTH=bm of pm_run.sh.
#   bash -l pm_visc1d.sh build ; sbatch pm_visc1d.sh        (from $P/tools, with the Intel modules of the login shell)
# build: pm_visc1d.f90 against libvilma.a of the S11 stack, with the flags of its vilma.x (make -n vilma);
# run:   pm_visc1d.x with the namelist of the 3-D run (B3_l170/pm.nml) -> $P/data/pm_visc1d_bagge6.csv, then
#        pm_prep.jl visc1d -> $P/data/pm_visc_bm.nc (the profile as a laterally uniform FesmData-form file).
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate; P=$G/runs/paleomist; V=$G/stacks/S11/vilma2; T=$P/tools
X=$P/data/pm_visc1d.x
if [ "${1:-}" = build ]; then
  cmd=$(make -s -C $V -n vilma | sed -n '/src\/vilma_main.f90/{N;s/\\\n//;p}')
  cmd=${cmd/src\/vilma_main.f90/$T/pm_visc1d.f90}; cmd=${cmd/-o bin\/vilma.x/-o $X}
  echo "$cmd"
  (cd $V && eval "$cmd")
  echo "built $X"; exit 0
fi
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd $P/B3_l170
$X pm.nml $P/data/pm_visc1d_bagge6.csv
export PATH=$HOME/.juliaup/bin:$PATH JULIA_DEPOT_PATH=$G/julia-depot
julia +1.13 --project=$G/p4w7/kryos/tools $T/pm_prep.jl visc1d $P/data/pm_visc1d_bagge6.csv $P/data/pm_visc_bm.nc
