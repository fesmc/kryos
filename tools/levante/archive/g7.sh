#!/bin/bash -l
# Gate G7 (phase 4 WP7c): the cost of the science case with geo = vilma2 at geo cadence dt = 10, 20, 50, 100 yr, 400 yr from -80 ka, stack S8.
#   g7.sh <dt>        prints the job id. Output as a production run would write it (rounded): 2D every 100 yr (the streams gia and vilma off, as the science par file), bundle at -79800.
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate; P=$G/p4w7c; export STK=S8
DT=$1
export PAR=$G/stacks/S8/kryos/par/kryos_NH-ANT_vilma2.nml JOBTIME=${JOBTIME:-01:30:00}
SINGLE=1 $P/levante/run_vilma2.sh g7 dt$DT clock.time_end=-79600 out_2D.dt=100 out_2Dsm.dt=100 out_rst.dt=200 \
    geo_glob.dt=$DT geo_nor.dt=$DT geo_ant.dt=$DT
