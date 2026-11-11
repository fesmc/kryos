# The environment of the Levante tools (tools/levante): sourced by every script here; the Julia tools read the same
# variables. Any of them can be set beforehand, in the environment or in ~/.kryos.env (sourced first; another file with
# KRYOS_ENV_FILE). Defaults:
#
#   KRYOS          the kryos checkout (the one this file is in)
#   KRYOS_WORK     your writable work folder: run folders, scores, calibrations (/work/ba1442/$USER/kryos-work)
#   FESMDATA       the FesmData store, read-only: never write into it except with FesmData's own tools
#                  (/work/ba1442/fesmdata); linked as $KRYOS/fesmdata by setup.sh
#   FESMDATA_SOURCES  the original data FesmData is built from, read-only (/work/ba1442/fesmdata-sources)
#   KRYOS_DATA     kryos reference data, read-only for runs: refcorr topographies, the released tau_b products
#                  ($KRYOS_DATA/plastic), the seed of the maps folder ($KRYOS_DATA/maps) (/work/ba1442/kryos_data); linked
#                  as $KRYOS/kryos_data by setup.sh
#   VILMA2         the vilma2 checkout of the standalone runs (vilma.x, input/, its .runme) ($KRYOS/vilma2)
#   PM_DIR         your folder of standalone VILMA2 PaleoMIST runs ($KRYOS_WORK/paleomist)
#   PM_REF         the reference standalone runs (B1_l170, ref/I910_it0, B3_l170, the inputs data/): read-only, the
#                  scoring references ($KRYOS_DATA/paleomist)
#   ACCOUNT        the Slurm account of the small jobs of the tools (the account of runme's .runme/config.toml, else ba1442)
#   KRYOS_RUNME    the runme command (runme): model runs are staged and submitted with runme (fesmc/runme)
#   KRYOS_JULIA    the julia command (julia; with juliaup "julia +1.13" if that channel is installed)
#
# Model runs (kryos.x, vilma.x) go through runme from the checkout (its .runme: the links input, fesmdata, kryos_data,
# maps into every run folder). Functions: kryos_jl <script.jl> [args]   julia with the project of tools/levante (run
# tools/levante/setup.sh once).

_kryos_env_file=${KRYOS_ENV_FILE:-$HOME/.kryos.env}
[ -f "$_kryos_env_file" ] && . "$_kryos_env_file"
unset _kryos_env_file

if [ -z "${KRYOS:-}" ]; then
    KRYOS=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../.." && pwd)
fi
export KRYOS
export KRYOS_WORK=${KRYOS_WORK:-/work/ba1442/$USER/kryos-work}
export FESMDATA=${FESMDATA:-/work/ba1442/fesmdata}
export FESMDATA_SOURCES=${FESMDATA_SOURCES:-/work/ba1442/fesmdata-sources}
export KRYOS_DATA=${KRYOS_DATA:-/work/ba1442/kryos_data}
export VILMA2=${VILMA2:-$KRYOS/vilma2}
export PM_DIR=${PM_DIR:-$KRYOS_WORK/paleomist}
export PM_REF=${PM_REF:-$KRYOS_DATA/paleomist}
if [ -z "${ACCOUNT:-}" ]; then
    ACCOUNT=$(sed -n 's/^ *account *= *"\([^"]*\)".*/\1/p' "$KRYOS/.runme/config.toml" 2> /dev/null | head -1)
    [ -n "$ACCOUNT" ] && [ "$ACCOUNT" != CHANGEME ] || ACCOUNT=ba1442
fi
export ACCOUNT
export KRYOS_RUNME=${KRYOS_RUNME:-runme}
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
