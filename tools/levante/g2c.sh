#!/bin/bash -l
# Gate G2c (phase 4 WP7c): the spin-up grid of the science case with geo = vilma2, 100 yr from the cold start at -80 ka, stack S8.
#   g2c.sh <equil_time_max [yr]> <equil_rate_tol [m/yr]> [more set_nml args]     prints the job id
# The par file is the science case with `equil_rate_tol` added to [vilma] (a copy in $G/p4w7c/par; the key is VILMA2's own default 1e-3).
# 2D output every 10 yr (the first 10 coupling intervals), streams gia/vilma on, lossless (the differences are small), no restart written.
# Run folder: $G/runs/g2c/kx-<sha>/S8/e<T>_t<tol>.
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate; P=$G/p4w7c; export STK=S8
EQ=$1; TOL=$2; shift 2
PARDIR=$P/par; mkdir -p $PARDIR
NAME=e$(python3 -c "print(f'{float(\"$EQ\"):.0e}'.replace('+0','').replace('+',''))")_t$(python3 -c "print(f'{float(\"$TOL\"):.0e}'.replace('-0','m').replace('-','m'))")
PAR=$PARDIR/kryos_NH-ANT_vilma2_$NAME.nml
awk -v tol="$TOL" '{print} /^&vilma/{inv=1} inv && /equil_time_max/{print "    equil_rate_tol    = " tol "     ! [m/yr] G2c"; inv=0}' $G/stacks/S8/kryos/par/kryos_NH-ANT_vilma2.nml > $PAR
export PAR JOBTIME=${JOBTIME:-01:30:00}
OV=(clock.time_end=-79900 out_2D.dt=10 out_1D.dt=10 'geo_glob.out_gia="out_2D"' 'geo_glob.out_vilma="out_2D"' kryos.write_lossless=True vilma.equil_time_max=$EQ)
SINGLE=1 exec $P/levante/run_vilma2.sh g2c $NAME "${OV[@]}" "$@"
