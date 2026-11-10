#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=ba1442
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --job-name=taub-fit
#SBATCH --output=fit-%x-%j.out
# Stage 1 of the tau_b calibration (tests/plastic_taub_fit.f90): the per-domain present-day factors of one grid.
#   cd <fit folder>; K=<kryos checkout> sbatch <K>/tools/levante/taub_fit/fit.sh <K>/tools/levante/taub_fit/fit_NH.nml
# The namelists here are those of the released fits: fit_NH.nml and fit_ANT.nml (runs/taub-fit, v2), fit_PAT.nml
# (added in runs/taub-fit-v3, whose factor file holds the three grids). Then the report and the factor file:
#   julia --project=<julia-plot> <K>/tools/levante/plastic_taub_fit_report.jl <out folder> fit_NH-32KM.nc fit_ANT-16KM.nc ...
# Needs: make libkryos/bin/plastic_taub_fit.x plastic=1
echo cpus $SLURM_CPUS_PER_TASK; export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16} OMP_STACKSIZE=512M OMP_MAX_ACTIVE_LEVELS=1
ulimit -s unlimited
K=${K:-/work/ba1442/robinson/kryos-ice-plastic}
/usr/bin/time -f "%e s %M kB" "$K/libkryos/bin/plastic_taub_fit.x" "$1"
