#!/bin/bash
# smb_ensemble.sh: an ensemble of the SMB tuning on imposed ice (ice backend data), on Levante or albedo (env.sh):
#
#   tools/levante/smb_ensemble.sh submit <ens dir> <params.txt> [group.key=value ...]
#       copies the parameter file to <ens dir>/params.txt (runme's format: a header of group.key, a row per member) and
#       submits, as two runme job arrays from the checkout, every member's
#         pm/<i>: par/kryos_NH-ANT_data_paleoemu.nml, an LGM snapshot from -22 to -20 ka (the emulator's filters start in
#                 equilibrium; 20 ka within a few % of the transient from -50 ka, about 3 min on 32 threads), or with
#                 TRANSIENT=1 from -50 ka to 0 (about 20 min);
#         pd/<i>: par/kryos_GRL-ANT_data_paleoemu_pd.nml, the present day on the reference topography (about 3 min);
#       the key=value arguments override keys of both par files for every member.
#   tools/levante/smb_ensemble.sh post <ens dir> [smb_fit.jl options]
#       submits itself ($KRYOS_PART, 2 h): smb_basin_budget.jl on every pm member (20 ka; the LGM window 19-23 ka for a
#       transient), then smb_fit.jl into <ens dir>/fit (default options: lgm=-20:-20, or -23:-19 for a transient;
#       exclude=Cordillera_and_Alaska,British_Isles). Prints the job id.
#
# Queues: albedo 30min (12h for a transient), Levante shared; OMP threads 32 (SMB_OMP). Each member's run folder holds
# its par file with the member's values; runme writes params_ensemble.txt.
set -euo pipefail
_h=${SMB_TOOLS:-$(dirname "$(readlink -f "$0")")}; . "$_h/env.sh"      # SMB_TOOLS: in the job ($0 is Slurm's copy)
[ $# -ge 2 ] || { sed -n '2,18p' "$0"; exit 1; }
mode=$1; E=$(readlink -f "$2"); shift 2
omp=${SMB_OMP:-32}
case $mode in
submit)
    [ $# -ge 1 ] || { echo "submit needs <params.txt>"; exit 1; }
    P=$(readlink -f "$1"); shift
    ov=("$@"); pdp=(); [ ${#ov[@]} -gt 0 ] && pdp=(-p "${ov[@]}")
    mkdir -p "$E"; cp "$P" "$E/params.txt"
    n=$(( $(grep -cv '^\s*$' "$E/params.txt") - 1 ))
    if [ "${TRANSIENT:-0}" = 1 ]; then
        clock="clock.time_init=-50e3"; echo transient > "$E/mode"
        case $KRYOS_MACHINE in albedo) qpm="-q 12h -w 01:30:00" ;; *) qpm="-q shared -w 01:30:00" ;; esac
    else
        clock="clock.time_init=-22e3 clock.time_end=-20e3 out_2D.dt=1000.0"; echo snapshot > "$E/mode"
        case $KRYOS_MACHINE in albedo) qpm="-q 30min" ;; *) qpm="-q shared -w 00:30:00" ;; esac
    fi
    case $KRYOS_MACHINE in albedo) qpd="-q 30min" ;; *) qpd="-q shared -w 00:30:00" ;; esac
    cd "$KRYOS"
    $KRYOS_RUNME -r -s -e kryos -n par/kryos_NH-ANT_data_paleoemu.nml -o "$E/pm" -i "$E/params.txt" --array "$n" \
        -p $clock ${ov[@]+"${ov[@]}"} $qpm --omp $omp -m 60G --jobname smb-pm | tail -1
    $KRYOS_RUNME -r -s -e kryos -n par/kryos_GRL-ANT_data_paleoemu_pd.nml -o "$E/pd" -i "$E/params.txt" --array "$n" \
        ${pdp[@]+"${pdp[@]}"} $qpd --omp $omp -m 60G --jobname smb-pd | tail -1
    echo "$n members in $E (pm: $(cat "$E/mode")); then: tools/levante/smb_ensemble.sh post $E" ;;
post)
    if [ -z "${SLURM_JOB_ID:-}" ]; then
        sbatch --parsable -A "$ACCOUNT" $KRYOS_PART -t 02:00:00 -c 4 --mem=16G -J smb-post -o "$E/post.out" \
            --export=ALL,SMB_TOOLS="$_h" "$_h/smb_ensemble.sh" post "$E" "$@"
        exit 0
    fi
    if [ "$(cat "$E/mode" 2> /dev/null)" = transient ]; then win="tmin=-23 tmax=-19"; lgm="lgm=-23:-19"; else win="tmin=-20 tmax=-20"; lgm="lgm=-20:-20"; fi
    for d in "$E"/pm/[0-9]*; do
        grep -q "kryos_run: end" "$d/out.out" 2> /dev/null || { echo "$d has not finished"; exit 1; }
        kryos_jl "$_h/smb_basin_budget.jl" "$d" "$d/budget" $win > "$d/budget.log" 2>&1
    done
    for d in "$E"/pd/[0-9]*; do
        grep -q "kryos_run: end" "$d/out.out" 2> /dev/null || { echo "$d has not finished"; exit 1; }
    done
    opts=(${@+"$@"}); printf '%s\n' "${opts[@]}" | grep -q '^lgm=' || opts+=("$lgm")
    printf '%s\n' "${opts[@]}" | grep -q '^exclude=' || opts+=("exclude=Cordillera_and_Alaska,British_Isles")
    kryos_jl "$_h/smb_fit.jl" "$E" "$E/fit" "${opts[@]}" ;;
*) echo "mode: submit or post"; exit 1 ;;
esac
