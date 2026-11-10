# Sourced by plastic_esmda.sh and plastic_ensemble.sh (after env.sh).
# the members of an iteration (members.txt) as one job array, staged and submitted by runme from the folder; prints the job id
submit_members() {
    local C=$1 I=$2 name=$3
    source "$C/slurm.env"
    (cd "$C" && $KRYOS_RUNME -r -s -e kryos -n "$I/par.nml" -o "$I" -i "$I/scale_files.txt" -j "$(paste -sd, "$I/members.txt")" \
        --array "$THROTTLE" -q "$PART" -w "$TIME" --omp "$THREADS" -m "$MEM" --jobname "$name") > "$I/runme.out" 2>&1 \
        || { tail -5 "$I/runme.out" >&2; return 1; }
    sed -n 's/.*"job_id": "\([0-9]*\)".*/\1/p' "$I/array.json"
}
