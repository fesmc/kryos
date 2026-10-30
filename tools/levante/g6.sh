#!/bin/bash -l
# Gate G6 (phase 4 WP7c): the science case smoke, both domains, geo = vilma2, from -80 ka, stack S9 (STK=S8 for the first runs). Prints the job ids.
#   g6.sh cont|chunk            2 kyr (restart bundle at the middle, -79000)
#   g6.sh smb-cont|smb-chunk    the smb_simple variant (par/gate/g6-smb.nml, tools/bin/kryos-smb-par.jl; the -p list of the conversion note is not enough), 1 kyr (restart at -79500)
# Output: 1-D every 10 yr, 2D every 100 yr, the streams gia and vilma (lossless, G3 (iii)) every 50 yr, bundles at the middle and the end.
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate; P=$G/p4w7c; export STK=${STK:-S9}
export PAR=${PAR:-$P/par/kryos_NH-ANT_vilma2.nml} JOBTIME=${JOBTIME:-04:00:00}   # the science par of p4-wp7c (spin-up 1e-5 / 2e4 yr), copied to $P/par; the spin-up adds about 1 h to the cold start
COM=(out_1D.dt=10 out_2D.dt=100 out_2Dsm.dt=50 'geo_glob.out_gia="out_2Dsm"' 'geo_glob.out_vilma="out_2Dsm"' kryos.write_lossless=True)
case ${1:-} in
  cont)       SINGLE=1 $P/levante/run_vilma2.sh g6 cont "${COM[@]}" clock.time_end=-78000 out_rst.dt=1000 ;;
  chunk)      CT0=-80000 CTM=-79000 CTE=-78000 $P/levante/run_vilma2.sh g6 chunk "${COM[@]}" out_rst.dt=1000 ;;
  smb-cont)   PAR=$P/par/g6-smb.nml SINGLE=1 $P/levante/run_vilma2.sh g6-smb cont "${COM[@]}" clock.time_end=-79000 out_rst.dt=500 ;;
  smb-chunk)  PAR=$P/par/g6-smb.nml CT0=-80000 CTM=-79500 CTE=-79000 $P/levante/run_vilma2.sh g6-smb chunk "${COM[@]}" out_rst.dt=500 ;;
  *) echo "usage: g6.sh cont|chunk|smb-cont|smb-chunk"; exit 1 ;;
esac
