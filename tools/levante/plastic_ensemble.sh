#!/bin/bash
# The jobs of a weighted ensemble (plastic_ensemble.jl; docs/dev/plastic-ensemble-design.md):
#   plastic_ensemble.sh start <config.toml> <ensemble folder>   init (wave 1) as a job, then its members and the weights
#   plastic_ensemble.sh submit <ensemble folder> <k>            the members of wave k (job array) and the weights after them
# The members run as those of a calibration (one runme job array, plastic_esmda.sh). Internal: plastic_ensemble.sh init|weights ...
set -u
T=${ESMDA_TOOLS:-$(cd "$(dirname "$(readlink -f "$0")")" && pwd)}    # the jobs (copies of this script) find the tools here
export ESMDA_TOOLS=$T
. "$T/env.sh"
J="$KRYOS_JULIA --project=$T"
small="-A $ACCOUNT $KRYOS_PART -t 02:00:00 -c 4 --mem=32G"
it() { printf "%s/iter_%02d" "$1" "$2"; }
. "$T/members.sh"

case $1 in
start)
    CFG=$(realpath "$2"); C=$(realpath -m "$3")
    [ -e "$C" ] && { echo "exists: $C"; exit 1; }
    [ -x "$KRYOS/libkryos/bin/kryos.x" ] || { echo "no $KRYOS/libkryos/bin/kryos.x: make kryos plastic=1 vilma=1"; exit 1; }
    mkdir -p "$C"
    sbatch --parsable $small -J ens-init -o "$C/init.out" "$T/plastic_ensemble.sh" init "$CFG" "$C" ;;
init)
    $J "$T/plastic_ensemble.jl" init "$2" "$3" || exit 1
    bash "$T/plastic_ensemble.sh" submit "$3" 1 ;;
submit)
    C=$2; k=$3; I=$(it "$C" "$k")
    jid=$(submit_members "$C" "$I" "ens-$k") || exit 1
    wid=$(sbatch --parsable $small -J "ens-w$k" -d "afterany:$jid" -o "$I/weights.out" "$T/plastic_ensemble.sh" weights "$C" "$k") || exit 1
    echo "$(date -Is) wave $k: members $jid (array $(paste -sd, "$I/members.txt")), weights $wid" | tee -a "$C/jobs.txt" ;;
weights)
    $J "$T/plastic_ensemble.jl" weights "$2" "$3" ;;
*)
    sed -n 2,6p "$0"; exit 1 ;;
esac
