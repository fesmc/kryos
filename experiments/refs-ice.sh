#!/bin/bash
# The imposed-ice rungs of the ladder (docs/simulations.md, I1-I8): an ice history imposed on VILMA2 with one Earth, its
# present-day correction r iterated until the 0 ka bed is today's (tools/levante/gia_refcorr.sh: one job per case).
#
#   experiments/refs-ice.sh <id> ...        ids: I1 ... I8, or all
#
# A case is <ice>_<earth> (the simulation `<ice> | <earth>` of docs/simulations.md) under $KRYOS_WORK/kryos_data/refcorr; the
# scores read its converged run (the reference `run` of tools/levante/plastic_ensemble_v2.toml). Restartable: rerunning a
# finished case only redoes its reference files. Run times on albedo: 1-D 3-5 min per iteration, VILMA2's 1Db about 17 min,
# 3-D about 2 h; 3-5 iterations.
# The GLAC-1D histories (FesmData IceHistory-GLAC1D_80k, v2.5.0) have no bed: today's comes from PaleoMIST (ZREF, XZREF).
# Before that FesmData release, point GLAC1D_GLOBAL at a local build of the product (its Global folder).
set -euo pipefail
X=$(cd "$(dirname "$(readlink -f "$0")")" && pwd); T=$X/../tools/levante; . "$T/env.sh"
[ $# -ge 1 ] || { sed -n '2,13p' "$0"; exit 1; }
F=$FESMDATA; G=${GLAC1D_GLOBAL:-$F/Global}
PM=$F/Global/GLOBAL-0.25DEG/GLOBAL-0.25DEG_IceHistory-PaleoMIST.nc
P1=$KRYOS/par/kryos_NH-ANT_plastic_vilma2_b1.nml                     # E1D-PMIST: Gowan et al. (2021)'s 1-D Earth
PB=$KRYOS/par/kryos_NH-ANT-PAT_plastic_vilma2_bagge1d.nml            # E1D-B21-6: Bagge 6's radial mean (VILMA2's 1Db)
PL=$KRYOS/par/kryos_NH-ANT-PAT_plastic_vilma2_bagge1d_layers.nml     # E1D-B21-6-layers: the same as plain 1-D layer means
P3=$KRYOS/par/kryos_NH-ANT-PAT_plastic_vilma2.nml                    # E3D-B21-6: Bagge et al. (2021) structure 6
PARS=$KRYOS_WORK/experiments/par; mkdir -p "$PARS"

# case <name> <par> <ice: PMIST | GLAC1D_80k | GLAC1D_nn9927_80k> [<cpus>]
case_() {
    local c=$1 par=$2 ice=$3 ncpu=${4:-128}
    if [ "$ice" = PMIST ]; then
        env -u ZREF -u XZREF HIST=$PM CASE=$c PAR=$par NCPU=$ncpu bash "$T/gia_refcorr.sh"
    else
        HIST=$G/GLOBAL-0.25DEG/GLOBAL-0.25DEG_IceHistory-$ice.nc XHIST=$G/GLOBAL-0.5DEG/GLOBAL-0.5DEG_IceHistory-$ice.nc \
            ZREF=$PM CASE=$c PAR=$par NCPU=$ncpu bash "$T/gia_refcorr.sh"
    fi | sed "s/^/$c: job /"
}
# the par file of Bagge's structure k: the 3-D science par with visc_3d_structure = k
b21() { local p=$PARS/E3D-B21-$1.nml; cp "$P3" "$p"; python3 "$T/set_nml.py" "$p" "vilma.visc_3d_structure=$1"; echo "$p"; }

ids=("$@"); [ "$1" = all ] && ids=(I1 I2 I3 I4 I5 I6 I7 I8)
for id in "${ids[@]}"; do
    case $id in
    I1) case_ ice-PMIST_E1D-PMIST "$P1" PMIST 16 ;;
    I2) case_ ice-PMIST_E1D-B21-6 "$PB" PMIST
        case_ ice-PMIST_E1D-B21-6-layers "$PL" PMIST 16 ;;
    I3) case_ ice-PMIST_E3D-B21-6 "$P3" PMIST ;;
    I4) for k in 1 2 3 4 5 7 8; do case_ ice-PMIST_E3D-B21-$k "$(b21 $k)" PMIST; done ;;
    I5) case_ ice-GLAC1D_E1D-PMIST "$P1" GLAC1D_80k 16 ;;
    I6) case_ ice-GLAC1D_E1D-B21-6-layers "$PL" GLAC1D_80k 16 ;;
    I7) case_ ice-GLAC1D_E3D-B21-6 "$P3" GLAC1D_80k ;;
    I8) case_ ice-GLAC1D-9927_E1D-PMIST "$P1" GLAC1D_nn9927_80k 16
        case_ ice-GLAC1D-9927_E3D-B21-6 "$P3" GLAC1D_nn9927_80k ;;
    *)  echo "unknown id $id (I1 ... I8, all)"; exit 1 ;;
    esac
done
