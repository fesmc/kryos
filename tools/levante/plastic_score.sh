#!/bin/bash
# plastic_score.sh <kryos run folder> <label> [<score folder>]: the run's VILMA2 stream (geo_glob_vilma.nc) against the
# GAPSLIP 2.1 index points with pm_compare.jl, beside the standalone references B1 (B1_l170, rho_ice 931) and I910_it0 (B1 at
# 910 kg m-3) of PM_REF, as in docs/dev/vilma2-paleomist.md. Submits itself as a small job; the report goes to
# <score folder>/out (default $KRYOS_WORK/score/<label>/out: pm_indicators.md and the tables of pm_compare.jl).
set -euo pipefail
_h=$(dirname "$(readlink -f "$0")"); [ -f "$_h/env.sh" ] || _h=$KRYOS/tools/levante; . "$_h/env.sh"
if [ "${KRYOS_JOB:-}" = 1 ]; then
    RUN=$1; S=$3
    kryos_jl "$_h/pm_compare.jl" "$FESMDATA" "$S" "$S/out" B1="$PM_REF/B1_l170" B1_170="$PM_REF/ref/I910_it0" B3="$RUN"
    exit
fi
[ $# -ge 2 ] || { sed -n '2,5p' "$0"; exit 1; }
RUN=$(readlink -f "$1"); LBL=$2; S=${3:-$KRYOS_WORK/score/$LBL}
[ -f "$RUN/geo_glob_vilma.nc" ] || { echo "no $RUN/geo_glob_vilma.nc (the run's VILMA2 stream: out_vilma of geo_glob)"; exit 1; }
for r in B1_l170 ref/I910_it0; do [ -d "$PM_REF/$r" ] || { echo "no reference $PM_REF/$r (PM_REF)"; exit 1; }; done
mkdir -p "$S"; S=$(readlink -f "$S")
sbatch --parsable -A "$ACCOUNT" $KRYOS_PART -t 01:00:00 -c 4 --mem=32G -J "score-$LBL" -o "$S/score-%j.out" \
    --export=ALL,KRYOS_JOB=1 "$_h/plastic_score.sh" "$RUN" "$LBL" "$S"
