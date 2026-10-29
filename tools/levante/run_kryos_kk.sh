#!/bin/bash -l
# Stage and submit a run of a kryos-against-kryos gate (tools/gate/<gate>.toml, kind = "kryos"; phase 3 WP5).
#   [STACK=S5b] [K=<kryos checkout>] [KSHA=<commit>] [MAPS=$G/maps/S5] [TOOLS=<checkout with the Julia project>] \
#       run_kryos_kk.sh <gate> <run>
# <run> is a [[run]] of the gate file: its par file (par/gate/<par>.nml) with its overrides `kp` (merged group names,
# `kryos.domains` reversed for `swap`); `chunk = true`: c00 0 -> time(1) and c01 from c00/restart-<t> (afterok), as
# run_kryos.sh. The code is KSHA in K (built here with the gate flags). Runs go to $G/runs/<gate>/kx-<sha>/<stack>/<run>.
# The job template is that of the `nor` reference runs (time limit raised by JOBTIME, default 01:00:00).
set -euo pipefail
module load python3 >/dev/null 2>&1 || true
export PATH=$HOME/.juliaup/bin:$HOME/.local/bin:$PATH
G=/work/ba1442/robinson/kryos-gate; K=${K:-$G/kryos}
GATE=$1; RUN=$2; STACK=${STACK:-S5b}; MAPS=${MAPS:-$G/maps/S5}; TOOLS=${TOOLS:-$K}
GF="-O2 -fp-model precise -march=znver3 -traceback"
CM="env PYTHONPATH=$G/tools/configme-src/src python3 -m configme"
JL="env JULIA_DEPOT_PATH=$G/julia-depot julia +1.13 --project=$TOOLS/tools $TOOLS/tools/bin/kryos-gate-kk.jl"
cd $K
[ -n "${KSHA:-}" ] && git checkout -q --detach $KSHA
SHA=$(git rev-parse --short=7 HEAD)
if [ ! -x $K/libkryos/bin/kryos.x ] || [ "$(cat $K/libkryos/bin/.kryos_sha 2>/dev/null)" != "$SHA" ]; then
    $CM config kryos --only -m dkrz_levante -c ifx > /dev/null; $CM deps . > /dev/null
    BLOG=$G/stacks/kryos-build-$SHA.log; [ "$K" = "$G/kryos" ] || BLOG=$K/build.log
    make kryos openmp=1 DFLAGS_NODEBUG="$GF" > $BLOG 2>&1 || { tail -30 $BLOG; exit 1; }
    echo $SHA > $K/libkryos/bin/.kryos_sha
fi
OUT=$G/runs/$GATE/kx-$SHA/$STACK/$RUN
[ -e $OUT ] && { echo "exists: $OUT"; exit 1; }
STAGE=$($JL stage $GATE $RUN)
PAR=$(echo "$STAGE" | sed -n 's/^PAR=//p')
mapfile -t KP < <(echo "$STAGE" | sed -n 's/^KP=//p')
CHUNK=$(echo "$STAGE" | sed -n 's/^CHUNK=//p')

stage () {   # stage <dir> <jobname> <set_nml args...>
    local dir=$1 name=$2; shift 2
    mkdir -p $dir/input; cd $dir
    for f in $G/stacks/$STACK/yelmox/input/*; do ln -sf $f input/; done
    for f in $K/input/*; do ln -sf $f input/; done
    ln -sf $K/fesm-utils/par/datasrc_defaults.nml input/datasrc_defaults.nml
    for d in fesmdata ice_data isostasy_data; do ln -s /work/ba1442/$d $d; done
    ln -s $K/fesm-utils fesm-utils; ln -s $MAPS maps
    cp $K/par/gate/$PAR.nml par.nml; cp $K/libkryos/bin/kryos.x .
    local args=("${KP[@]}" "$@")
    [ ${#args[@]} -gt 0 ] && python3 $G/tools/set_nml.py par.nml "${args[@]}"
    git -C $K rev-parse HEAD > kryos_commit.txt
    sed -e "s/--job-name=.*/--job-name=$name/" -e "s/--time=.*/--time=${JOBTIME:-01:00:00}/" -e 's|^srun .*|srun ./kryos.x par.nml|' \
        $(ls $G/ref/nor/yx-*/S*/cont/job.submit | head -1) > job.submit
}
if [ -n "$CHUNK" ]; then
    read -r T0 TM TE <<< "$CHUNK"
    stage $OUT/c00 $GATE-kx-c00 clock.time_init=$T0 clock.time_end=$TM
    j0=$(cd $OUT/c00 && sbatch --parsable job.submit)
    RD=$(python3 -c "t=float('$TM'); print(int(t) if t==int(t) else t)")
    stage $OUT/c01 $GATE-kx-c01 clock.time_init=$TM clock.time_end=$TE "kryos.restart=\"$OUT/c00/restart-$RD\""
    j1=$(cd $OUT/c01 && sbatch --parsable --dependency=afterok:$j0 job.submit)
    echo "$GATE $RUN $SHA $j0 $j1 $OUT"
else
    stage $OUT $GATE-kx-$RUN
    echo "$GATE $RUN $SHA $(cd $OUT && sbatch --parsable job.submit) $OUT"
fi
