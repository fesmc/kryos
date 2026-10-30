#!/bin/bash -l
# Gate G5 (phase 4 WP7c): ANT deglaciation from PaleoMIST at -20 ka to 0, Antarctica the only active domain, stack S8.
#   g5.sh vilma2|fastiso <run name> [set_nml args ...]       prints the job id
# Par: par/gate/g5-vilma2.nml, par/gate/g5-fastiso.nml (copies in $G/p4w7c/par; tools/bin/kryos-ant-par.jl). Defaults of the
# full run: 2D and restart schedules below, streams gia and vilma (lossless) every 100 yr on the schedule out_2Dsm.
# Run folder $G/runs/g5/kx-<sha>/S8/<run name>. THREADS, JOBTIME as run_vilma2.sh.
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate; P=$G/p4w7c; export STK=S8
KIND=$1; NAME=$2; shift 2
export PAR=$P/par/g5-$KIND.nml JOBTIME=${JOBTIME:-04:00:00}
OV=(out_2D.dt=500 out_2Dsm.dt=100 out_1D.dt=10 out_rst.dt=5000 kryos.write_lossless=True)
[ $KIND = vilma2 ] && OV+=('geo_glob.out_gia="out_2Dsm"' 'geo_glob.out_vilma="out_2Dsm"')
SINGLE=1 exec $P/levante/run_vilma2.sh g5 $KIND-$NAME "${OV[@]}" "$@"
