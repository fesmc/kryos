#!/bin/bash
# The jobs of a weighted ensemble (plastic_ensemble.jl; docs/dev/plastic-ensemble-design.md):
#   plastic_ensemble.sh start <config.toml> <ensemble folder>   init (wave 1) as a job, then its members and the weights
#   plastic_ensemble.sh submit <ensemble folder> <k>            the members of wave k (job array) and the weights after them
# The members run as those of a calibration (plastic_esmda.sh member). Internal: plastic_ensemble.sh init|weights ...
set -u
G=/work/ba1442/robinson/kryos-gate
T=${ESMDA_TOOLS:-$(cd "$(dirname "$0")" && pwd)}
export PATH=$HOME/.juliaup/bin:$PATH JULIA_DEPOT_PATH=$G/julia-depot ESMDA_TOOLS=$T
J="julia +1.13 --project=$G/p4w7/kryos/tools"
small="-A ba1442 -p shared -t 02:00:00 -c 4 --mem=32G"
it() { printf "%s/iter_%02d" "$1" "$2"; }

case $1 in
start)
    CFG=$(realpath "$2"); C=$(realpath -m "$3")
    [ -e "$C" ] && { echo "exists: $C"; exit 1; }
    mkdir -p "$C"
    sbatch --parsable $small -J ens-init -o "$C/init.out" "$T/plastic_ensemble.sh" init "$CFG" "$C" ;;
init)
    $J "$T/plastic_ensemble.jl" init "$2" "$3" || exit 1
    bash "$T/plastic_ensemble.sh" submit "$3" 1 ;;
submit)
    C=$2; k=$3; I=$(it "$C" "$k")
    source "$C/slurm.env"
    ids=$(paste -sd, "$I/members.txt")
    thr=""; [ "$THROTTLE" -gt 0 ] && thr="%$THROTTLE"
    jid=$(sbatch --parsable -A "$ACCT" -p "$PART" -t "$TIME" -N 1 -n 1 -c "$THREADS" --mem="$MEM" -J "ens-$k" \
          --array="$ids$thr" -o "$I/slurm-%a.out" "$T/plastic_esmda.sh" member "$C" "$k") || exit 1
    wid=$(sbatch --parsable $small -J "ens-w$k" -d "afterany:$jid" -o "$I/weights.out" "$T/plastic_ensemble.sh" weights "$C" "$k") || exit 1
    echo "$(date -Is) wave $k: members $jid (array $ids), weights $wid" | tee -a "$C/jobs.txt" ;;
weights)
    $J "$T/plastic_ensemble.jl" weights "$2" "$3" ;;
*)
    sed -n 2,6p "$0"; exit 1 ;;
esac
