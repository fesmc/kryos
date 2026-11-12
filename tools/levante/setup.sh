#!/bin/bash
# setup.sh: the one-time setup of the HPC tools, run in a login shell on Levante or albedo from anywhere:
#   tools/levante/setup.sh
# Checks the environment (env.sh: the paths, their access); links the data stores into the checkout (fesmdata, kryos_data:
# runme links them into every run folder) and seeds the checkout's maps/ with the shared remapping weights
# ($KRYOS_DATA/maps); runme (fesmc/runme, >= 0.5.27) and its .runme/config.toml (from .runme/config.default.toml); Julia
# (juliaup with the 1.13 channel) and the Julia project of the tools (tools/levante/Project.toml, Manifest.toml: CairoMakie,
# NCDatasets, DelimitedFiles; in your own depot, JULIA_DEPOT_PATH or ~/.julia); creates KRYOS_WORK; reports what is built.
# Safe to rerun.
set -uo pipefail
_h=$(dirname "$(readlink -f "$0")"); . "$_h/env.sh"
ok() { printf "  ok    %s\n" "$*"; }
no() { printf "  MISSING  %s\n" "$*"; bad=1; }
bad=0
echo "kryos checkout: $KRYOS (machine $KRYOS_MACHINE)"
echo "paths (env.sh; override in the environment or ~/.kryos.env):"
for v in FESMDATA KRYOS_DATA PM_REF; do      # KRYOS_DATA/paleoemu: the inputs of the paleoemu par files
    d=${!v}; [ -r "$d" ] && ok "$v = $d (readable)" || no "$v = $d (not readable)"
done
mkdir -p "$KRYOS_WORK" 2> /dev/null && [ -w "$KRYOS_WORK" ] && ok "KRYOS_WORK = $KRYOS_WORK (writable)" || no "KRYOS_WORK = $KRYOS_WORK (cannot write)"
echo "checkout links (runme links them into every run folder):"
for l in fesmdata:$FESMDATA kryos_data:$KRYOS_DATA; do
    n=${l%%:*}; t=${l#*:}
    if [ ! -e "$KRYOS/$n" ] && [ ! -L "$KRYOS/$n" ]; then ln -s "$t" "$KRYOS/$n" && ok "$n -> $t (created)"
    elif [ "$(readlink -f "$KRYOS/$n")" = "$(readlink -f "$t")" ]; then ok "$n -> $t"
    else no "$n -> $(readlink "$KRYOS/$n"), not $t: remove it, or set the variable to it"; fi
done
for n in ice_data isostasy_data; do      # linked by runme when present; read by the Yelmo par files only
    [ -e "$KRYOS/$n" ] && ok "$n -> $(readlink "$KRYOS/$n")" || echo "  -     $n (only the Yelmo par files read it; runme warns that it is missing, harmless for the plastic runs)"
done
nw=$(ls "$KRYOS/maps"/scrip-*.nc 2> /dev/null | wc -l)
cp -n "$KRYOS_DATA"/maps/grid_*.txt "$KRYOS_DATA"/maps/scrip-*.nc "$KRYOS/maps"/ 2> /dev/null
ok "maps/: $(ls "$KRYOS/maps"/grid_*.txt | wc -l) grids, $(ls "$KRYOS/maps"/scrip-*.nc 2> /dev/null | wc -l) weights ($(( $(ls "$KRYOS/maps"/scrip-*.nc 2> /dev/null | wc -l) - nw )) copied from $KRYOS_DATA/maps)"
echo "runme:"
if ! command -v ${KRYOS_RUNME%% *} > /dev/null 2>&1; then
    no "runme: pip install --user git+https://github.com/fesmc/runme  (then rerun this script)"
else
    rv=$($KRYOS_RUNME --version 2> /dev/null | awk '{print $2}')
    if printf '%s\n0.5.27\n' "$rv" | sort -V -C 2> /dev/null && [ "$rv" != 0.5.27 ]; then
        no "runme $rv: the ensemble tools need >= 0.5.27 (--array): runme update"
    else ok "runme $rv"; fi
    if [ ! -f "$KRYOS/.runme/config.toml" ]; then
        case $KRYOS_MACHINE in
            albedo) sed -e 's/^hpc = "dkrz_levante"/hpc = "awi_albedo"/' -e "s/^account = \"ba1442\"/account = \"$ACCOUNT\"/" \
                        "$KRYOS/.runme/config.default.toml" > "$KRYOS/.runme/config.toml" ;;
            *) cp "$KRYOS/.runme/config.default.toml" "$KRYOS/.runme/config.toml" ;;
        esac && ok ".runme/config.toml (created from config.default.toml for $KRYOS_MACHINE, account $ACCOUNT; edit it for another account)"
    else ok ".runme/config.toml (account $ACCOUNT)"; fi
fi
echo "builds:"
[ -x "$KRYOS/libkryos/bin/kryos.x" ] && ok "kryos.x" || no "kryos.x: make kryos plastic=1 vilma=1 yelmo=0 (getting-started)"
[ -x "$KRYOS/libkryos/bin/plastic_taub_fit.x" ] && ok "plastic_taub_fit.x" || echo "  -     plastic_taub_fit.x (for Stage 1 only: make libkryos/bin/plastic_taub_fit.x plastic=1)"
[ -x "$VILMA2/bin/vilma.x" ] && ok "vilma.x ($VILMA2)" || echo "  -     vilma.x (for the standalone runs only: make vilma in $VILMA2)"
echo "julia:"
if ! command -v ${KRYOS_JULIA%% *} > /dev/null 2>&1; then
    no "julia: install juliaup and Julia 1.13 (your own, in your home):"
    echo "          curl -fsSL https://install.julialang.org | sh -s -- --yes && juliaup add 1.13"
    echo "        then rerun this script (KRYOS_JULIA=\"julia +1.13\" is picked up when the channel exists)"
else
    v=$($KRYOS_JULIA -e 'print(VERSION)' 2> /dev/null)
    case $v in 1.13*) ok "$KRYOS_JULIA ($v)" ;; *) no "$KRYOS_JULIA is $v: the tools' Manifest is for 1.13 (juliaup add 1.13)" ;; esac
    echo "  instantiating tools/levante (first time: some minutes) in ${JULIA_DEPOT_PATH:-$HOME/.julia} ..."
    if $KRYOS_JULIA --project="$_h" -e 'using Pkg; Pkg.instantiate(); using CairoMakie, NCDatasets' > /dev/null 2>&1; then
        ok "Julia project tools/levante"
    else
        no "Julia project tools/levante: run  $KRYOS_JULIA --project=$_h -e 'using Pkg; Pkg.instantiate()'  to see why"
    fi
fi
[ $bad = 0 ] && echo "setup: all set" || { echo "setup: fix the MISSING items above"; exit 1; }
