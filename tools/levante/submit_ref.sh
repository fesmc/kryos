#!/bin/bash -l
# Stage and submit one yelmox reference run of the kryos gate (regression-design §1).
#   submit_ref.sh <stack> <gate> <run> [MAPS]
#     stack: v24 | S | S2 | S3 | S4        (the stack under $G/stacks)
#     gate:  grl | ant (phase 1) | grp | pyr | srg | nor | lis | rgrl | eant | egrl (phase 2, phase2-plan §3.1)
#     run:   cont | cont2 | chunk | mapgen | base (15 kyr, grl/ant only, regression-design 1.1)
#     MAPS:  map folder linked as <rundir>/maps (default $G/maps/S for grl/ant, $G/maps/S4 for the phase 2 gates)
# Environment: EXE (runme -e executable, default yelmox = libyelmox/bin/yelmox.x of the stack), PART (partition queue alias, default shared), WALL (hh:mm:ss, default 01:00:00),
#              NCPU (CPUs allocated, default 16; `shared` gives 940 MB per CPU, so more memory = more CPUs;
#              `runme --cpus` does not change the template, so --cpus-per-task is edited and OMP threads stay 16),
#              OUTROOT (default $G/ref).
# Prints "<run> <jobid> [<jobid>]" on stdout.
# chunk mirrors yelmox scripts/run_chunks.sh (runme -s per chunk in <dir>/cNN,
# time_init/time_end per chunk, ctrl.restart = previous chunk's bundle, sbatch
# afterok), staged by hand so the maps link can be repointed before submission.
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate
STACK=$1; GATE=$2; RUN=$3
ROOT=$G/stacks/$STACK/yelmox
SHA=$(git -C "$ROOT" rev-parse --short=7 HEAD)
PART=${PART:-shared}; WALL=${WALL:-01:00:00}; NCPU=${NCPU:-}
OUTROOT=${OUTROOT:-$G/ref}; EXE=${EXE:-yelmox}

# Gate table: par file, group holding the timeline, time_init, time_end, chunk boundary, restart dt, extra -p.
RST=500; EXTRA=()
case $GATE in
    grl)  PAR=yelmox/yelmox_Greenland.nml;        TG=ctrl;   TI=0;      TE=1000;   TM=500 ;;
    ant)  PAR=yelmox/yelmox_Antarctica.nml;       TG=ctrl;   TI=0;      TE=1000;   TM=500;   EXTRA=(sim.relax=False sim.opt=False) ;;
    grp)  PAR=yelmox/yelmox_GreenlandPaleo.nml;   TG=ctrl;   TI=0;      TE=1000;   TM=500 ;;
    nor)  PAR=yelmox/yelmox_North.nml;            TG=ctrl;   TI=0;      TE=1000;   TM=500 ;;
    lis)  PAR=yelmox/yelmox_LIS.nml;              TG=ctrl;   TI=0;      TE=1000;   TM=500 ;;
    rgrl) PAR=yelmox/yelmox_rembo_Greenland.nml;  TG=ctrl;   TI=0;      TE=1000;   TM=500 ;;
    pyr)  PAR=yelmox/yelmox_Pyrenees.nml;         TG=ctrl;   TI=-120e3; TE=-119e3; TM=-119.5e3 ;;
    srg)  PAR=yelmox/yelmox_SRG.nml;              TG=ctrl;   TI=0;      TE=200;    TM=100;   RST=100 ;;
    eant) PAR=yelmox/yelmox_esm_Antarctica.nml;   TG=spinup; TI=0;      TE=1000;   TM=500 ;;
    egrl) PAR=yelmox/yelmox_esm_Greenland.nml;    TG=spinup; TI=0;      TE=1000;   TM=500 ;;
    *) echo "unknown gate $GATE" >&2; exit 1 ;;
esac
case $GATE in grl|ant) MAPS=${4:-$G/maps/S} ;; *) MAPS=${4:-$G/maps/S4} ;; esac
RUNOPTS=(-q "$PART" -w "$WALL" --omp 16)
[ "$RUN" = base ] && { RUNOPTS=(-q "$PART" -w 02:00:00 --omp 16); TE=15000; }
P=(ctrl.write_lossless=True tm_2D.dt=100 tm_2Dsm.dt=10 tm_1D.dt=10 tm_rst.dt=$RST "${EXTRA[@]}")

# bundle folder name of a time, as yelmox (kryos_startup.f90 restart_bundle_dir): restart-<kyr, f20.3>-kyr
bundle () { printf "restart-%.3f-kyr" "$(awk -v t="$1" 'BEGIN{printf "%.6f", t/1000}')"; }

if [ "$RUN" = mapgen ]; then
    OUT=$G/maps/gen/$GATE
else
    OUT=$OUTROOT/$GATE/yx-$SHA/$STACK/$RUN
fi
[ -e "$OUT" ] && { echo "exists: $OUT" >&2; exit 1; }
cd "$ROOT"

# stage_submit <dir> <jobname> <dep jobid or ""> <params...>; repoint maps; submit; echo jobid
stage_submit () {
    local dir=$1 name=$2 dep=$3; shift 3
    runme -s "${RUNOPTS[@]}" --jobname "$name" -e "$EXE" -n "$PAR" -o "$dir" -p "$@" > "$dir.runme.log" 2>&1 \
        || { cat "$dir.runme.log" >&2; exit 1; }
    if [ -n "$NCPU" ]; then
        sed -i -e "s/^#SBATCH --cpus-per-task=.*/#SBATCH --cpus-per-task=${NCPU}/" \
               -e 's/^export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK/export OMP_NUM_THREADS=16/' "$dir/job.submit"
    fi
    rm "$dir/maps"
    ln -s "$MAPS" "$dir/maps"
    if [ -n "$dep" ]; then
        (cd "$dir" && sbatch --parsable --dependency=afterok:"$dep" job.submit)
    else
        (cd "$dir" && sbatch --parsable job.submit)
    fi
}

mkdir -p "$(dirname "$OUT")"
case $RUN in
    cont|cont2)
        j=$(stage_submit "$OUT" "$GATE-$STACK-$RUN" "" $TG.time_end=$TE "${P[@]}")
        echo "$RUN $j" ;;
    base)
        j=$(stage_submit "$OUT" "$GATE-$STACK-base" "" $TG.time_end=$TE "${P[@]}")
        echo "$RUN $j" ;;
    mapgen)
        TEN=$(awk -v t="$TI" 'BEGIN{printf "%g", t+10}')
        j=$(stage_submit "$OUT" "$GATE-mapgen" "" $TG.time_end=$TEN "${P[@]}")
        echo "$RUN $j" ;;
    chunk)
        mkdir -p "$OUT"
        j0=$(stage_submit "$OUT/c00" "$GATE-$STACK-chunk0" "" $TG.time_init=$TI $TG.time_end=$TM "${P[@]}")
        j1=$(stage_submit "$OUT/c01" "$GATE-$STACK-chunk1" "$j0" $TG.time_init=$TM $TG.time_end=$TE "${P[@]}" \
             ctrl.restart="$OUT/c00/$(bundle $TM)")
        echo "$RUN $j0 $j1" ;;
    *) echo "unknown run $RUN" >&2; exit 1 ;;
esac
