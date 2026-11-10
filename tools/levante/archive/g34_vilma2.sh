#!/bin/bash -l
# Gates G3 and G4 of phase 4 WP7: the runs of the science case with geo = vilma2 (gate sci-vilma2), 400 yr from -80 ka.
#   g34_vilma2.sh [cont cont2 chunk lossless]       (default: all four); prints the job ids
# Overrides (set_nml): 400 yr, 2D and VILMA2 streams every 50 yr, bundles every 200 yr (chunk 0 -> -79800 -> -79600).
# lossless: the same as cont with [kryos] write_lossless (G4 rounding check).
set -euo pipefail
P=/work/ba1442/robinson/kryos-gate/p4w7
export PAR=$P/kryos_NH-ANT_vilma2.nml JOBTIME=${JOBTIME:-01:30:00}
OV=(clock.time_end=-79600 out_2D.dt=50 out_2Dsm.dt=50 out_rst.dt=200 'geo_glob.out_gia="out_2D"' 'geo_glob.out_vilma="out_2D"')
for r in ${@:-cont cont2 chunk lossless}; do
  case $r in
    cont|cont2) $P/levante/run_vilma2.sh sci-vilma2 $r "${OV[@]}" ;;
    chunk)      CT0=-80000 CTM=-79800 CTE=-79600 $P/levante/run_vilma2.sh sci-vilma2 chunk "${OV[@]}" ;;
    lossless)   $P/levante/run_vilma2.sh sci-vilma2 lossless "${OV[@]}" kryos.write_lossless=True ;;
  esac
done
