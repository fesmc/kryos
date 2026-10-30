#!/bin/bash -l
# Gate G6 (phase 4 WP7c): the science case smoke, both domains, geo = vilma2, from -80 ka, stack S8. Prints the job ids.
#   g6.sh cont|chunk            2 kyr (restart bundle at the middle, -79000)
#   g6.sh smb-cont|smb-chunk    the smb_simple variant (conversion note), 1 kyr (restart at -79500)
# Output: 1-D every 10 yr, 2D every 100 yr, the streams gia and vilma (lossless, G3 (iii)) every 50 yr, bundles at the middle and the end.
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate; P=$G/p4w7c; export STK=S8
export PAR=$G/stacks/S8/kryos/par/kryos_NH-ANT_vilma2.nml JOBTIME=${JOBTIME:-02:00:00}
COM=(out_1D.dt=10 out_2D.dt=100 out_2Dsm.dt=50 'geo_glob.out_gia="out_2Dsm"' 'geo_glob.out_vilma="out_2Dsm"' kryos.write_lossless=True)
SMB=('surface_nor.backend="smb_simple"' 'surface_ant.backend="smb_simple"')
case ${1:-} in
  cont)       SINGLE=1 $P/levante/run_vilma2.sh g6 cont "${COM[@]}" clock.time_end=-78000 out_rst.dt=1000 ;;
  chunk)      CT0=-80000 CTM=-79000 CTE=-78000 $P/levante/run_vilma2.sh g6 chunk "${COM[@]}" out_rst.dt=1000 ;;
  smb-cont)   SINGLE=1 $P/levante/run_vilma2.sh g6-smb cont "${COM[@]}" "${SMB[@]}" clock.time_end=-79000 out_rst.dt=500 ;;
  smb-chunk)  CT0=-80000 CTM=-79500 CTE=-79000 $P/levante/run_vilma2.sh g6-smb chunk "${COM[@]}" "${SMB[@]}" out_rst.dt=500 ;;
  *) echo "usage: g6.sh cont|chunk|smb-cont|smb-chunk"; exit 1 ;;
esac
