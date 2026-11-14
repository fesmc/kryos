#!/bin/bash
# Many small gia_refcorr.sh cases in few jobs (the 1-D Earths of the Earth search, docs/dev/earth-search-design.md): each job is a
# pool of NWORK workers on its own cores, each worker runs the next case of the list that is neither done nor claimed, until the
# list is done or the job's time is nearly up. For when the number of jobs is limited and a case needs a fraction of a node.
#   LIST=cases.txt NJOB=4 tools/levante/gia_refcorr_pool.sh        (submits NJOB jobs; resubmit the same list to go on)
# The list: one case per line, "CASE PAR [VAR=value ...]" (the env of gia_refcorr.sh for that case: HIST, XHIST, R0, ...); # comments.
# Env: LIST; NJOB (1); NWORK (8) workers per job, each on NCPU/NWORK cores (taskset) with as many OpenMP threads; NCPU (128),
#   MEM (64G), TIME (12:00:00), PART (partition and QOS: KRYOS_PART; KRYOS_PART_LONG for longer jobs, TIME in hours, e.g. 36:00:00); MARGIN [min] (60): a worker claims no case after TIME - MARGIN, nor when
#   less than 1.2 times the longest finished case of the list is left (3-D cases: 6-11 h on 32 cores); KD and the
#   rest of gia_refcorr.sh's env for every case (KEEP=last for many cases). Case output: $KD/refcorr/<CASE>/pool.out.
# A case is done when its iter.log ends in "== done" (or "not converged"/"failed": then it is not retried; delete its folder to rerun).
# A claim (the link <case>/claim -> the job id) of a job that no longer runs is taken over, so a resubmission resumes the cases of a job
# that was cut off (gia_refcorr.sh restarts after the last finished iteration).
set -uo pipefail
_h=$(dirname "$(readlink -f "$0")"); [ -f "$_h/env.sh" ] || _h=$KRYOS/tools/levante; . "$_h/env.sh"
LIST=$(readlink -f "${LIST:?}"); KD=${KD:-$KRYOS_WORK/kryos_data}
NWORK=${NWORK:-8}; NCPU=${NCPU:-128}; TIME=${TIME:-12:00:00}; MARGIN=${MARGIN:-60}
if [ "${KRYOS_JOB:-}" != 1 ]; then
    mkdir -p $KD/refcorr
    for j in $(seq 1 ${NJOB:-1}); do
        sbatch --parsable -A $ACCOUNT ${PART:-$KRYOS_PART} -t $TIME -c $NCPU --mem=${MEM:-64G} -J gia_refcorr_pool \
            -o $KD/refcorr/pool-%j.out --export=ALL,KRYOS_JOB=1,LIST=$LIST,KD=$KD "$_h/gia_refcorr_pool.sh"
    done
    exit
fi
# every case would fail at once without these (a job submitted from a shell without the login environment)
command -v ${KRYOS_RUNME%% *} > /dev/null || { echo "no ${KRYOS_RUNME%% *} on PATH: submit from a login shell (bash -l)"; exit 1; }
[ -x $VILMA2/bin/vilma.x ] || { echo "no $VILMA2/bin/vilma.x"; exit 1; }
IFS=: read -r hh mm ss <<< "$TIME"; jend=$(( $(date +%s) + 10#$hh * 3600 + 10#$mm * 60 )); deadline=$(( jend - MARGIN * 60 ))
# the longest case of the list so far [s] (start to "== done" in its iter.log): a worker claims a case only with 1.2 times that left
longest() {
    local c rest t0 t1 m=0
    while read -r c rest; do
        [ -z "$c" ] || [[ $c == \#* ]] && continue
        f=$KD/refcorr/$c/iter.log; grep -q "^== done" $f 2> /dev/null || continue
        t0=$(sed -n '1s/.* \([0-9-]*T[0-9:+-]*\)$/\1/p' $f); t1=$(grep -m1 "^== done" $f | grep -o "[0-9-]*T[0-9:+-]*$")
        [ -n "$t0" ] && [ -n "$t1" ] && m=$(( $(date -d $t1 +%s) - $(date -d $t0 +%s) > m ? $(date -d $t1 +%s) - $(date -d $t0 +%s) : m ))
    done < $LIST
    echo $m
}
# the cores of the job, NCPU/NWORK per worker
cpus=(); for r in $(taskset -cp $$ | sed 's/.*: //' | tr ',' ' '); do
    case $r in *-*) cpus+=($(seq ${r%-*} ${r#*-})) ;; *) cpus+=($r) ;; esac
done
per=$(( ${#cpus[@]} / NWORK )); [ $per -ge 1 ] || { echo "$NWORK workers on ${#cpus[@]} cores"; exit 1; }
echo "== pool $SLURM_JOB_ID on $(hostname): $NWORK workers x $per cores, list $LIST, claims until $(date -d @$deadline -Is)"
finished() { tail -1 $1/iter.log 2> /dev/null | grep -q "^== done\|^not converged\|failed"; }
claim() {   # claim the case folder $1: the link <case>/claim -> <job id> (atomic); a claim of a job that no longer runs is taken over
    local c=$1/claim j
    mkdir -p $1; ln -sT $SLURM_JOB_ID $c 2> /dev/null && return 0
    [ -L $c ] || return 1
    j=$(readlink $c) && squeue -h -j $j 2> /dev/null | grep -q . && return 1
    rm -f $c; ln -sT $SLURM_JOB_ID $c 2> /dev/null
}
worker() {
    local w=$1 set case par vars
    set=$(IFS=,; echo "${cpus[*]:$(( w * per )):$per}")
    while read -r case par vars; do
        [ -z "$case" ] || [[ $case == \#* ]] && continue
        [ $(date +%s) -lt $deadline ] && [ $(( $(date +%s) + 12 * $(longest) / 10 )) -lt $jend ] || break
        O=$KD/refcorr/$case; finished $O && continue; claim $O || continue
        echo "worker $w: $case $(date -Is)"
        env KRYOS_JOB=1 CASE=$case PAR=$par NCPU=$per NTHR=$per $vars taskset -c $set bash $_h/gia_refcorr.sh > $O/pool.out 2>&1
        echo "worker $w: $case exit $? $(date -Is) ($(tail -1 $O/iter.log 2> /dev/null))"
        rm -f $O/claim
    done < $LIST
}
for w in $(seq 0 $(( NWORK - 1 ))); do worker $w & sleep 2; done
wait
echo "== pool $SLURM_JOB_ID done $(date -Is)"
