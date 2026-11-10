# The environment of the Levante tools (tools/levante): sourced by every script here; the Julia tools read the same
# variables. Any of them can be set beforehand, in the environment or in ~/.kryos.env (sourced first; another file with
# KRYOS_ENV_FILE). Defaults:
#
#   KRYOS          the kryos checkout (the one this file is in)
#   KRYOS_WORK     your writable work folder: run folders, scores, calibrations (/work/ba1442/$USER/kryos-work)
#   KRYOS_MAPS     the maps folder of the runs: grid descriptions and remapping weights; runs write new weights into it
#                  ($KRYOS_WORK/maps; created by kryos_maps_init from the checkout's maps/ and KRYOS_MAPS_SEED)
#   KRYOS_MAPS_SEED  a folder of existing grid files and weights to copy from, read-only ($KRYOS_DATA/maps)
#   FESMDATA       the FesmData store, read-only: never write into it except with FesmData's own tools
#                  (/work/ba1442/fesmdata)
#   FESMDATA_SOURCES  the original data FesmData is built from, read-only (/work/ba1442/fesmdata-sources)
#   KRYOS_DATA     kryos reference data, read-only for runs: refcorr topographies, the released tau_b products
#                  ($KRYOS_DATA/plastic) (/work/ba1442/kryos_data)
#   VILMA2         the vilma2 checkout of the standalone runs (vilma.x, input/) ($KRYOS/vilma2)
#   PM_DIR         your folder of standalone VILMA2 PaleoMIST runs ($KRYOS_WORK/paleomist)
#   PM_REF         the reference standalone runs (B1_l170, ref/I910_it0, B3_l170, the inputs data/): read-only, the
#                  scoring references ($KRYOS_DATA/paleomist)
#   ACCOUNT        the Slurm account (ba1442)
#   KRYOS_JULIA    the julia command (julia; with juliaup "julia +1.13" if that channel is installed)
#
# Functions: kryos_jl <script.jl> [args]   julia with the project of tools/levante (run tools/levante/setup.sh once);
#            kryos_maps_init               create KRYOS_MAPS if missing (grid files of the checkout, seeded weights);
#            kryos_link_data <folder>      the links fesmdata, kryos_data, input, maps of a run folder.

_kryos_env_file=${KRYOS_ENV_FILE:-$HOME/.kryos.env}
[ -f "$_kryos_env_file" ] && . "$_kryos_env_file"
unset _kryos_env_file

if [ -z "${KRYOS:-}" ]; then
    KRYOS=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../.." && pwd)
fi
export KRYOS
export KRYOS_WORK=${KRYOS_WORK:-/work/ba1442/$USER/kryos-work}
export KRYOS_MAPS=${KRYOS_MAPS:-$KRYOS_WORK/maps}
export FESMDATA=${FESMDATA:-/work/ba1442/fesmdata}
export FESMDATA_SOURCES=${FESMDATA_SOURCES:-/work/ba1442/fesmdata-sources}
export KRYOS_DATA=${KRYOS_DATA:-/work/ba1442/kryos_data}
export KRYOS_MAPS_SEED=${KRYOS_MAPS_SEED:-$KRYOS_DATA/maps}
export VILMA2=${VILMA2:-$KRYOS/vilma2}
export PM_DIR=${PM_DIR:-$KRYOS_WORK/paleomist}
export PM_REF=${PM_REF:-$KRYOS_DATA/paleomist}
export ACCOUNT=${ACCOUNT:-ba1442}
if [ -z "${KRYOS_JULIA:-}" ]; then
    # juliaup's launcher by its full path: Levante's python module (mambaforge) puts an old julia first on PATH
    KRYOS_JULIA=julia
    _ju=$HOME/.juliaup/bin
    if [ -x "$_ju/julia" ] && "$_ju/juliaup" status 2> /dev/null | grep -qE "^ *\*? *1\.13 "; then KRYOS_JULIA="$_ju/julia +1.13"
    elif [ -x "$_ju/julia" ]; then KRYOS_JULIA=$_ju/julia; fi
    unset _ju
fi
export KRYOS_JULIA

kryos_jl() {
    $KRYOS_JULIA --project="$KRYOS/tools/levante" "$@"
}

kryos_maps_init() {
    [ -d "$KRYOS_MAPS" ] && return 0
    mkdir -p "$KRYOS_MAPS"
    cp -n "$KRYOS"/maps/grid_*.txt "$KRYOS_MAPS"/ 2> /dev/null || true
    if [ -d "$KRYOS_MAPS_SEED" ]; then
        cp -n "$KRYOS_MAPS_SEED"/grid_*.txt "$KRYOS_MAPS_SEED"/scrip-*.nc "$KRYOS_MAPS"/ 2> /dev/null || true
    fi
    echo "kryos_maps_init: $KRYOS_MAPS ($(ls "$KRYOS_MAPS" | wc -l) files)"
}

kryos_link_data() {
    local d=$1
    kryos_maps_init
    ln -sfn "$FESMDATA" "$d/fesmdata"
    ln -sfn "$KRYOS_DATA" "$d/kryos_data"
    ln -sfn "$KRYOS/input" "$d/input"
    ln -sfn "$KRYOS_MAPS" "$d/maps"
}
