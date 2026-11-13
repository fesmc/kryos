#!/bin/bash
# The online rungs of the ladder (docs/simulations.md, T1, D1, D3): kryos with PlasticIceSheet every step on the bed that
# VILMA2 deforms, from the par file of ens-v2's member m000 (D2: Gowan et al. (2021)'s tau_b times the present-day factors,
# Bagge 6, the refcorr reference of PaleoMIST-sci).
#
#   experiments/refs-tau.sh <ens-v2 folder> <id> ...      ids: T1, D1, D3, or all
#
#   T1  tau-PMIST | E1D-PMIST              m000 on Gowan's 1-D Earth without the present-day factors and scales
#   D1  tau-PMIST-PD | E1D-PMIST           m000 on Gowan's 1-D Earth (the factors kept)
#   D3  tau-PMIST-PD | E3D-B21-6 @NH16     m000 with the Northern Hemisphere on NH-16KM instead of NH-32KM
#
# The runs go to $KRYOS_WORK/runs/pm-online/{gowan,pdfac} and $KRYOS_WORK/runs/res16/m000-NH16 (runme, queue $KRYOS_QUEUE);
# about 3 h each on albedo (64 threads), D3 about 3.3 h. The ensemble folder (experiments/ens.sh v2) gives m000's par file,
# the runme setup and kryos.x.
set -euo pipefail
X=$(cd "$(dirname "$(readlink -f "$0")")" && pwd); T=$X/../tools/levante; . "$T/env.sh"
[ $# -ge 2 ] || { sed -n '2,15p' "$0"; exit 1; }
E=$(realpath "$1"); shift
M0=$E/iter_01/0/par.nml; [ -f "$M0" ] || { echo "no $M0: run experiments/ens.sh v2 first"; exit 1; }

# a run folder set up as the ensemble's: runme config, input/, kryos.x and the data links
setup() { mkdir -p "$1"; cd "$1"; cp -r "$E/.runme" "$E/input" .; cp "$E/kryos.x" .
          for l in fesmdata kryos_data maps; do ln -sfn "$(readlink -f "$E/$l")" $l; done; }
submit() { runme -r -s -e kryos -n "$1" -o "$2" -q "$KRYOS_QUEUE" -w 11:00:00 --omp 64 -m 30G --jobname "$3" 2>&1 | tail -1; }

# Gowan's 1-D Earth (par/kryos_NH-ANT_plastic_vilma2_b1.nml) in the &vilma group of a 3-D par file
EARTH=(
 "vilma.n_layer=12"
 "vilma.r_top=6371.0e3, 6291.0e3, 6251.0e3, 6151.0e3, 5971.0e3, 5771.0e3, 5701.0e3, 5600.0e3, 4943.0e3, 4287.0e3, 3630.0e3, 3480.0e3"
 "vilma.r_bot=6291.0e3, 6251.0e3, 6151.0e3, 5971.0e3, 5771.0e3, 5701.0e3, 5600.0e3, 4943.0e3, 4287.0e3, 3630.0e3, 3480.0e3, 0.0"
 "vilma.eta=1.0e40, 1.0e40, 4.0e20, 4.0e20, 4.0e20, 4.0e20, 4.0e22, 4.0e22, 4.0e22, 4.0e22, 4.0e22, 0.0"
 "vilma.rheology=0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2"
 "vilma.l_visc_3d=.false." "vilma.visc_shift_file=''" "vilma.pre_spinup_1d=.false." "vilma.equil_time_max=2.0e6" "vilma.equil_rate_tol=1.0e-5")
NOFAC=(); for i in nor ant pat; do NOFAC+=("ice_$i.taub_domain_factor_file=''" "ice_$i.taub_scale_file=''"); done

ids=("$@"); [ "$1" = all ] && ids=(T1 D1 D3)
for id in "${ids[@]}"; do
    case $id in
    T1|D1)
        R=$KRYOS_WORK/runs/pm-online; ( setup "$R"
        c=pdfac; [ $id = T1 ] && c=gowan
        cp "$M0" par_$c.nml; python3 "$T/set_nml.py" par_$c.nml "${EARTH[@]}"
        [ $id = T1 ] && python3 "$T/set_nml.py" par_$c.nml "${NOFAC[@]}"
        submit par_$c.nml $c pm-$c ) ;;
    D3)
        R=$KRYOS_WORK/runs/res16; ( setup "$R"
        sed "s/'NH-32KM'/'NH-16KM'/g" "$M0" > par.nml
        [ "$(grep -c "'NH-16KM'" par.nml)" = 2 ] || { echo "D3: expected grid and grid_hub NH-16KM in par.nml"; exit 1; }
        submit par.nml m000-NH16 res16 ) ;;
    *)  echo "unknown id $id (T1, D1, D3, all)"; exit 1 ;;
    esac
done
