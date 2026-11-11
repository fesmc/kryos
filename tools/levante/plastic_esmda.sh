#!/bin/bash
# The jobs of the tau_b calibration by ES-MDA (plastic_esmda.jl; docs/dev/plastic-backend.md, "Calibration"):
#   plastic_esmda.sh start <config.toml> <calibration folder>   init (prior ensemble) as a job, then iteration 1
#   plastic_esmda.sh submit <calibration folder> <k>            the members of iteration k (job array) and its analysis after them
# The analysis scores the members, updates and submits iteration k + 1 until the folder's `next` says done. Internal:
#   plastic_esmda.sh init|analyse ...                            the bodies of the jobs
# The members of an iteration are one job array staged and submitted by runme from the folder (a frozen runme project:
# plastic_esmda.jl init): iter_kk/<m> with par.nml, the multipliers of scale/mNNN.nc, the links, the run and `status` (exit
# code, seconds; the folder's job template). A failed member is dropped by the analysis.
set -u
T=${ESMDA_TOOLS:-$(cd "$(dirname "$(readlink -f "$0")")" && pwd)}    # the jobs (copies of this script) find the tools here
export ESMDA_TOOLS=$T
. "$T/env.sh"
J="$KRYOS_JULIA --project=$T"
small="-A $ACCOUNT $KRYOS_PART -t 01:00:00 -c 4 --mem=32G"
it() { printf "%s/iter_%02d" "$1" "$2"; }
. "$T/members.sh"

case $1 in
start)
    CFG=$(realpath "$2"); C=$(realpath -m "$3")
    [ -e "$C" ] && { echo "exists: $C"; exit 1; }
    [ -x "$KRYOS/libkryos/bin/kryos.x" ] || { echo "no $KRYOS/libkryos/bin/kryos.x: make kryos plastic=1 vilma=1"; exit 1; }
    mkdir -p "$C"
    sbatch --parsable $small -J esmda-init -o "$C/init.out" "$T/plastic_esmda.sh" init "$CFG" "$C" ;;
init)
    $J "$T/plastic_esmda.jl" init "$2" "$3" || exit 1
    bash "$T/plastic_esmda.sh" submit "$3" 1 ;;
submit)
    C=$2; k=$3; I=$(it "$C" "$k")
    jid=$(submit_members "$C" "$I" "esmda-$k") || exit 1
    aid=$(sbatch --parsable $small -J "esmda-up$k" -d "afterany:$jid" -o "$I/update.out" "$T/plastic_esmda.sh" analyse "$C" "$k") || exit 1
    echo "$(date -Is) iteration $k: members $jid (array $(paste -sd, "$I/members.txt")), analysis $aid" | tee -a "$C/jobs.txt" ;;
analyse)
    C=$2; k=$3
    $J "$T/plastic_esmda.jl" update "$C" "$k" || exit 1
    nxt=$(cat "$C/next")
    [ "$nxt" = done ] || bash "$T/plastic_esmda.sh" submit "$C" "$nxt" ;;
*)
    sed -n 2,9p "$0"; exit 1 ;;
esac
