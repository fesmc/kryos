#!/bin/bash
# The Earth-only ensemble earth-pca-v1 (docs/simulations.md: ice-PMIST | E3D-PCA-E<nnn>; docs/dev/earth-search-design.md):
# PaleoMIST's ice imposed on 30 3-D Earths sampled from a PCA basis of the expert viscosity fields (the 8 structures of
# Bagge et al. 2021, Lloyd et al. 2024, Pan et al. 2022), each run once with the present-day correction r of PaleoMIST-sci.
#
#   experiments/earth-pca.sh basis <v1|v2>     the basis (a job, about 30 min): $KRYOS_WORK/kryos_data/plastic/earth_basis_<v>.nc
#   experiments/earth-pca.sh sample            the 30 members of v1 (seed 1): $KRYOS_WORK/earth-pca-v1/E<nnn>.{nc,nml}, members.csv
#   experiments/earth-pca.sh run [<m1> <m2>]   the runs of members m1..m2 (default 1..30), one 3-D pass each (about 2 h)
#
# v1 weighs the fields equally; v2 by study (Bagge's 8 count as one) with a common lid (22.5). The ensemble used v1.
# PAN_LITHO is the lithosphere thickness of Pan et al. (2022) (default: albedo's isostasy_data copy); the field itself is
# VILMA2's input/pan2022.nc. The runs go to $KRYOS_WORK/kryos_data/refcorr/ice-PMIST_E3D-PCA-E<nnn>/it0 (R0: no iteration,
# no reference files).
set -euo pipefail
X=$(cd "$(dirname "$(readlink -f "$0")")" && pwd); T=$X/../tools/levante; . "$T/env.sh"
[ $# -ge 1 ] || { sed -n '2,13p' "$0"; exit 1; }
J="$KRYOS_JULIA --project=$T"
B=$KRYOS_WORK/kryos_data/plastic; S=$KRYOS_WORK/earth-pca-v1
PAN_LITHO=${PAN_LITHO:-/albedo/work/projects/p_forclima/isostasy_data/earth_structure/lithothickness/pan2022.nc}

case $1 in
basis)
    v=${2:?v1 or v2}; w=""; [ $v = v2 ] && w="weights=study lid=22.5"
    mkdir -p "$B"
    sbatch --parsable -A $ACCOUNT $KRYOS_PART -t 02:00:00 -c 32 --mem=120G -J earth-basis-$v -o "$B/earth_basis_$v.out" --wrap \
        "export JULIA_NUM_THREADS=32; $J $T/plastic_ensemble.jl earth-basis $B/earth_basis_$v.nc $VILMA2/input/pan2022.nc $PAN_LITHO kmax=6 $w && \
         $J $T/plastic_earth_basis_plot.jl $B/earth_basis_$v.nc $B/earth_basis_$v.png" ;;
sample)
    $J "$T/plastic_ensemble.jl" earth-sample "$B/earth_basis_v1.nc" "$KRYOS/par/kryos_NH-ANT-PAT_plastic_vilma2.nml" "$S" n=30 K=4 seed=1 ;;
run)
    m1=${2:-1}; m2=${3:-30}
    export R0=$KRYOS_DATA/refcorr/PaleoMIST-sci/r_3.nc DOMAINS="North:NH-32KM:North/NH-32KM/NH-32KM_TOPO-GEBCO2026.nc"
    for m in $(seq -f "%03g" $m1 $m2); do
        echo "E$m: job $(CASE=ice-PMIST_E3D-PCA-E$m PAR=$S/E$m.nml bash "$T/gia_refcorr.sh")"
    done ;;
*)  sed -n '2,13p' "$0"; exit 1 ;;
esac
