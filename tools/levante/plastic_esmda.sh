#!/bin/bash
# The jobs of the tau_b calibration by ES-MDA (plastic_esmda.jl; docs/dev/plastic-backend.md, "Calibration"):
#   plastic_esmda.sh start <config.toml> <calibration folder>   init (prior ensemble) as a job, then iteration 1
#   plastic_esmda.sh submit <calibration folder> <k>            the members of iteration k (job array) and its analysis after them
# The analysis scores the members, updates and submits iteration k + 1 until the folder's `next` says done. Internal:
#   plastic_esmda.sh init|member|analyse ...                     the bodies of the jobs
# A member: its folder iter_kk/mNNN (taub_scale.nc, par.nml), links to the frozen input/ and kryos.x of the folder, the run, and
# `status` (exit code, seconds). A failed member is dropped by the analysis.
set -u
G=/work/ba1442/robinson/kryos-gate
T=${ESMDA_TOOLS:-$(cd "$(dirname "$0")" && pwd)}
export PATH=$HOME/.juliaup/bin:$PATH JULIA_DEPOT_PATH=$G/julia-depot ESMDA_TOOLS=$T
J="julia +1.13 --project=$G/p4w7/kryos/tools"
small="-A ba1442 -p shared -t 01:00:00 -c 4 --mem=32G"
it() { printf "%s/iter_%02d" "$1" "$2"; }

case $1 in
start)
    CFG=$(realpath "$2"); C=$(realpath -m "$3")
    [ -e "$C" ] && { echo "exists: $C"; exit 1; }
    mkdir -p "$C"
    sbatch --parsable $small -J esmda-init -o "$C/init.out" "$T/plastic_esmda.sh" init "$CFG" "$C" ;;
init)
    $J "$T/plastic_esmda.jl" init "$2" "$3" || exit 1
    bash "$T/plastic_esmda.sh" submit "$3" 1 ;;
submit)
    C=$2; k=$3; I=$(it "$C" "$k")
    source "$C/slurm.env"
    ids=$(paste -sd, "$I/members.txt")
    thr=""; [ "$THROTTLE" -gt 0 ] && thr="%$THROTTLE"
    jid=$(sbatch --parsable -A "$ACCT" -p "$PART" -t "$TIME" -N 1 -n 1 -c "$THREADS" --mem="$MEM" -J "esmda-$k" \
          --array="$ids$thr" -o "$I/slurm-%a.out" "$T/plastic_esmda.sh" member "$C" "$k") || exit 1
    aid=$(sbatch --parsable $small -J "esmda-up$k" -d "afterany:$jid" -o "$I/update.out" "$T/plastic_esmda.sh" analyse "$C" "$k") || exit 1
    echo "$(date -Is) iteration $k: members $jid (array $ids), analysis $aid" | tee -a "$C/jobs.txt" ;;
member)
    C=$2; k=$3; source "$C/slurm.env"
    D=$(it "$C" "$k")/$(printf "m%03d" "$SLURM_ARRAY_TASK_ID")
    cd "$D" || exit 1
    ln -sfn "$C/input" input; ln -sfn /work/ba1442/fesmdata fesmdata; ln -sfn /work/ba1442/kryos_data kryos_data; ln -sfn "$MAPS" maps
    export OMP_PROC_BIND=true OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OMP_STACKSIZE=512M OMP_MAX_ACTIVE_LEVELS=1
    ulimit -s unlimited
    t0=$(date +%s)
    "$C/kryos.x" par.nml > out.out 2> out.err
    rc=$?
    echo "$rc $(( $(date +%s) - t0 ))" > status
    exit $rc ;;
analyse)
    C=$2; k=$3
    $J "$T/plastic_esmda.jl" update "$C" "$k" || exit 1
    nxt=$(cat "$C/next")
    [ "$nxt" = done ] || bash "$T/plastic_esmda.sh" submit "$C" "$nxt" ;;
*)
    sed -n 2,9p "$0"; exit 1 ;;
esac
