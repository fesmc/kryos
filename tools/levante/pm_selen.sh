#!/bin/bash
# PaleoMIST's own SELEN sea level at points, from the Stokes coefficients of the PANGAEA bundle (Gowan et al. 2021,
# doi:10.1594/PANGAEA.905800), with Gowan's program calc_sl/rslc.f90 (SELEN 2.9, GPL; run-standard.sh), built here unchanged. It is the
# check of the grid synthesis of pm_prep.jl selen (the program itself needs ~40 ms per point: 3 h for the 0.5 degree grid).
#   CASE=72_73_74_75/ehgr tools/levante/pm_selen.sh      (submits itself; the standard version: the final ice runs NA 72, EU 73, ANT 74, PAT 75)
#   CASE=72_73_74_75_h/ehgr tools/levante/pm_selen.sh    (the 500-yr version, ice load linearly interpolated; its own data.inc)
# Points: the sites of pm_compare.jl and a few more (test_sl.txt: lon lat). Output: $W/rsl_spreadsheet.dat (age lon lat rsl, rsl = S(t) - S(0)).
set -euo pipefail
_h=$(dirname "$(readlink -f "$0")"); [ -f "$_h/env.sh" ] || _h=$KRYOS/tools/levante; . "$_h/env.sh"
B=$FESMDATA_SOURCES/gowan2021_paleomist/ice_reconstruction
CASE=${CASE:?72_73_74_75/ehgr | 72_73_74_75_h/ehgr}
W=$PM_DIR/selen/${CASE//\//_}
if [ "${KRYOS_JOB:-}" != 1 ]; then
    mkdir -p $W
    sbatch --parsable -A $ACCOUNT -p shared -t 00:30:00 -n 1 -c 1 --mem=4G -J pm_selen -o $W/job-%j.out \
        --export=ALL,KRYOS_JOB=1 "$_h/pm_selen.sh"
    exit
fi
mkdir -p $W; cd $W
cp $B/calc_sl/HARMONICS.f90 $B/calc_sl/SHTOOLS.f90 $B/calc_sl/rslc.f90 .
# data.inc: the case folder's own if it has one (the 500-yr version: NN = 160, DELTA = 0.5), else calc_sl's (NN = 32, DELTA = 2.5)
d=$B/Stokes_Coefficients/${CASE%%/*}/data.inc; [ -f $d ] || d=$B/calc_sl/data.inc
cp $d data.inc
ln -sf $B/Stokes_Coefficients/$CASE/shs-ice5glgm-gsc-3it-rot-128-v3.bin .
cat > test_sl.txt <<EOF
275.0 60.0
20.0 63.0
40.0 75.0
185.0 -80.0
310.0 -77.0
300.4 13.1
0.25 -89.75
123.75 -45.25
200.25 30.75
EOF
gfortran -o shtools.o -c -O3 SHTOOLS.f90
gfortran -o harmonics.o -c -O3 HARMONICS.f90
gfortran -o calc_sl -O3 rslc.f90 harmonics.o shtools.o
echo "== calc_sl $CASE ($d) $(date -Is)"
./calc_sl > calc_sl.log
echo "== done $(date -Is): $(wc -l < rsl_spreadsheet.dat) lines"
