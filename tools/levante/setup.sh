#!/bin/bash
# setup.sh: the one-time setup of the Levante tools, run in a login shell on Levante from anywhere:
#   tools/levante/setup.sh
# Checks the environment (env.sh: the paths, their access), Julia (juliaup with the 1.13 channel), instantiates the
# Julia project of the tools (tools/levante/Project.toml, Manifest.toml: CairoMakie, NCDatasets, DelimitedFiles; in your
# own depot, JULIA_DEPOT_PATH or ~/.julia), creates KRYOS_WORK and KRYOS_MAPS, and reports what is built. Safe to rerun.
set -uo pipefail
_h=$(dirname "$(readlink -f "$0")"); . "$_h/env.sh"
ok() { printf "  ok    %s\n" "$*"; }
no() { printf "  MISSING  %s\n" "$*"; bad=1; }
bad=0
echo "kryos checkout: $KRYOS"
echo "paths (env.sh; override in the environment or ~/.kryos.env):"
for v in FESMDATA KRYOS_DATA PM_REF KRYOS_MAPS_SEED; do
    d=${!v}; [ -r "$d" ] && ok "$v = $d (readable)" || no "$v = $d (not readable)"
done
mkdir -p "$KRYOS_WORK" 2> /dev/null && [ -w "$KRYOS_WORK" ] && ok "KRYOS_WORK = $KRYOS_WORK (writable)" || no "KRYOS_WORK = $KRYOS_WORK (cannot write)"
kryos_maps_init > /dev/null && ok "KRYOS_MAPS = $KRYOS_MAPS ($(ls "$KRYOS_MAPS" | wc -l) files)"
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
