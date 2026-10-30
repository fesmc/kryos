#!/bin/bash -l
# The S11 gates (VILMA2 row order) on stack S11 (setup_S11.sh); scripts from the stack's kryos checkout. Prints the job ids.
#   gate_S11.sh g0     make check and make check-long of $S/kryos (vilma=1), make check of $S/kryos-v0
#   gate_S11.sh g2a    G2a: the cap, adapter vs vilma.x on the ascending Gauss inputs and on the lon-lat inputs (g2a.sh)
#   gate_S11.sh g2b    G2b: the PaleoMIST replay, kryos and vilma.x, lmax 85, 128, 170 (g2b.sh)
#   gate_S11.sh g4     G3 and G4: sci-vilma2 cont, cont2, chunk, lossless (400 yr from -80 ka; g34_vilma2.sh, as S10)
#   gate_S11.sh g6     G6: the science par (spin-up 1e-5 / 2e4), 2 kyr cont and chunk at -79000; smb_simple 1 kyr cont and chunk at -79500
#   gate_S11.sh g2c    G2c: the chosen spin-up (1e-5, 2e4 yr) with the science par, 100 yr (run e2e4_t1em5)
# Compare with compare_S11.sh.
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate; export STK=S11; S=$G/stacks/S11; K=$S/kryos; L=$K/tools/levante
RV=$L/run_vilma2.sh
mkdir -p $G/reports/S11
case ${1:-} in
  g0)
    sbatch --parsable --export=ALL,STK=S11,V=1,TGT="check check-long",LOGTAG=vilma1 $L/check_kryos.sh
    sbatch --parsable --export=ALL,STK=S11,V=0,LOGTAG=vilma0 $L/check_kryos.sh ;;
  g2a)
    sbatch --parsable --export=ALL,STK=S11,NSTEP=6,TAG=cap -o $G/reports/S11/g2a-%j.out $L/g2a.sh ;;
  g2b)
    for l in 85 128 170; do
        t=$(printf "l%03d" $l); tm=01:00:00; [ $l = 128 ] && tm=01:30:00; [ $l = 170 ] && tm=02:30:00
        for m in kryos vilma; do
            sbatch --parsable -t $tm --export=ALL,STK=S11,LMAX=$l,NSTEP=2100,MODE=$m,TAG=$t -o $G/reports/S11/g2b-$t-$m-%j.out $L/g2b.sh
        done
    done ;;
  g4)
    for r in cont cont2 lossless; do STK=S11 K=$K PAR=$G/p4w7/kryos_NH-ANT_vilma2.nml SINGLE=1 JOBTIME=01:30:00 \
        $RV sci-vilma2 $r clock.time_end=-79600 out_2D.dt=50 out_2Dsm.dt=50 out_rst.dt=200 'geo_glob.out_gia="out_2D"' 'geo_glob.out_vilma="out_2D"' \
        $([ $r = lossless ] && echo kryos.write_lossless=True); done
    STK=S11 K=$K PAR=$G/p4w7/kryos_NH-ANT_vilma2.nml JOBTIME=01:30:00 CT0=-80000 CTM=-79800 CTE=-79600 \
        $RV sci-vilma2 chunk clock.time_end=-79600 out_2D.dt=50 out_2Dsm.dt=50 out_rst.dt=200 'geo_glob.out_gia="out_2D"' 'geo_glob.out_vilma="out_2D"' ;;
  g6)
    COM=(out_1D.dt=10 out_2D.dt=100 out_2Dsm.dt=50 'geo_glob.out_gia="out_2Dsm"' 'geo_glob.out_vilma="out_2Dsm"' kryos.write_lossless=True)
    export JOBTIME=04:00:00
    PAR=$K/par/kryos_NH-ANT_vilma2.nml SINGLE=1 $RV g6 cont "${COM[@]}" clock.time_end=-78000 out_rst.dt=1000
    PAR=$K/par/kryos_NH-ANT_vilma2.nml CT0=-80000 CTM=-79000 CTE=-78000 $RV g6 chunk "${COM[@]}" out_rst.dt=1000
    PAR=$K/par/gate/g6-smb.nml SINGLE=1 $RV g6-smb cont "${COM[@]}" clock.time_end=-79000 out_rst.dt=500
    PAR=$K/par/gate/g6-smb.nml CT0=-80000 CTM=-79500 CTE=-79000 $RV g6-smb chunk "${COM[@]}" out_rst.dt=500 ;;
  g2c)
    PAR=$K/par/kryos_NH-ANT_vilma2.nml SINGLE=1 JOBTIME=01:30:00 $RV g2c e2e4_t1em5 clock.time_end=-79900 out_2D.dt=10 out_1D.dt=10 \
        'geo_glob.out_gia="out_2D"' 'geo_glob.out_vilma="out_2D"' kryos.write_lossless=True ;;
  *) echo "usage: gate_S11.sh g0|g2a|g2b|g4|g6|g2c"; exit 1 ;;
esac
