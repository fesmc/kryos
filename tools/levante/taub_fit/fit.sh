#!/bin/bash
# Stage 1 of the tau_b calibration (tests/plastic_taub_fit.f90): the per-domain present-day factors of one grid.
#   tools/levante/taub_fit/fit.sh <fit folder> <fit.nml> ...      one job per namelist, run in <fit folder>
#   tools/levante/taub_fit/fit.sh $KRYOS_WORK/taub-fit tools/levante/taub_fit/fit_NH.nml tools/levante/taub_fit/fit_ANT.nml
# The namelists here are those of the released factor files: fit_NH.nml and fit_ANT.nml (v2), fit_PAT.nml (added in v3,
# whose factor file holds the three grids; $KRYOS_DATA/plastic). {fesmdata} in a namelist is replaced by FESMDATA. Then
# the report and the factor file of all grids (lambda by the discrepancy principle or the L-curve corner):
#   kryos_jl tools/levante/plastic_taub_fit_report.jl <fit folder>/report <fit folder>/fit_NH-32KM.nc ...
# Needs: make libkryos/bin/plastic_taub_fit.x plastic=1. Each fit takes up to 2 h on 16 threads (ANT-16KM).
set -euo pipefail
_h=$(dirname "$(readlink -f "$0")")/..; [ -f "$_h/env.sh" ] || _h=$KRYOS/tools/levante; . "$_h/env.sh"
if [ "${KRYOS_JOB:-}" = 1 ]; then
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16} OMP_STACKSIZE=512M OMP_MAX_ACTIVE_LEVELS=1
    ulimit -s unlimited
    /usr/bin/time -f "%e s %M kB" "$KRYOS/libkryos/bin/plastic_taub_fit.x" "$1"
    exit
fi
[ $# -ge 2 ] || { sed -n '2,3p' "$0"; exit 1; }
X=$KRYOS/libkryos/bin/plastic_taub_fit.x
[ -x "$X" ] || { echo "no $X: make libkryos/bin/plastic_taub_fit.x plastic=1"; exit 1; }
D=$1; shift; mkdir -p "$D"; D=$(cd "$D" && pwd)
for nml in "$@"; do
    n=$(basename "$nml")
    sed -e "s|{fesmdata}|$FESMDATA|g" "$nml" > "$D/$n"
    (cd "$D" && sbatch --parsable -A "$ACCOUNT" -p shared -t 02:00:00 -c 16 --mem=24G -J "taub-fit" -o "fit-${n%.nml}-%j.out" \
        --export=ALL,KRYOS_JOB=1 "$_h/taub_fit/fit.sh" "$D/$n")
done
