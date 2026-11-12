#!/bin/bash -l
# stack.sh: a kryos stack on AWI's albedo, from scratch: the kryos checkout and the dependencies of the backends that
# BUILD turns on, at the refs of its .configme/manifest.toml, configured for albedo by configme (awi_albedo, ifx: the
# machine's flags come from configme, the repositories keep no machine files), linked and built; then
# tools/levante/setup.sh (the data links, maps, runme, the Julia project of the tools).
#
#   tools/albedo/stack.sh <dir> [<kryos ref>]      (from a kryos checkout; or curl it and run it anywhere)
#
#   <dir>        a new folder, e.g. /albedo/work/user/$USER/kryos-plastic; it gets kryos/ (the checkout, <kryos ref>,
#                default main) and the packages: fesm-utils, chion, insol always; yelmo (yelmo=1), FastIsostasy
#                (fastisostasy=1), vilma2 (vilma=1), PlasticIceSheet (plastic=1), PaleoEmulator (paleoemu=1). An existing
#                checkout in <dir> is kept as it is.
#   BUILD        the make flags of the backends (default "yelmo=0 fastisostasy=0 vilma=1 plastic=1 openmp=1": the plastic
#                ice backend with VILMA2); a backend is on when BUILD has <flag>=1
#   CONFIGME     the configme command, one that knows the kryos orchestrator and its packages (configme >= 0.16; default
#                configme). Before its release: CONFIGME="env PYTHONPATH=<configme checkout>/src python3 -m configme"
#
# The data come from the shared stores of p_forclima (tools/levante/env.sh): FesmData ($FESMDATA, kept up to date
# with FesmData's `julia fesmdata.jl mirror`, shared/machines/albedo.env) and kryos_data ($KRYOS_DATA, a copy of
# Levante's /work/ba1442/kryos_data: maps, paleoemu, plastic, refcorr, paleomist).
# fesm-utils is built serially (its Makefile has no module dependencies: -j races); its external libraries (FFTW, Lis,
# SHTns: make libs, slow, one-time) only when a backend needs them (yelmo, fastisostasy, vilma, plastic).
set -euo pipefail
[ $# -ge 1 ] || { sed -n '2,17p' "$0"; exit 1; }
W=$1; KREF=${2:-main}
BUILD=${BUILD:-"yelmo=0 fastisostasy=0 vilma=1 plastic=1 openmp=1"}
CONFIGME=${CONFIGME:-configme}
on() { [[ " $BUILD " == *" $1=1 "* ]]; }
omp=0; on openmp && omp=1
command -v git > /dev/null 2>&1 || module load git > /dev/null 2>&1
mkdir -p "$W"; cd "$W"
exec > >(tee -a "$W/stack.log") 2>&1
echo "== $(date -Is) kryos $KREF in $W ($BUILD)"
[ -d kryos ] || git clone -q -b "$KREF" git@github.com:fesmc/kryos.git kryos

# the packages: name, the backend flag that turns it on (- always), its make target in kryos (- none), linked to fesm-utils
pkgs=()
while read -r p flag target link; do
    [ "$flag" = - ] || on "$flag" || continue
    pkgs+=("$p:$target:$link")
done << 'LIST'
fesm-utils      -             -                 0
chion           -             chion-static      1
insol           -             insol-static      1
yelmo           yelmo         yelmo-static      1
FastIsostasy    fastisostasy  isostasy-static   1
vilma2          vilma         vilma-static      1
PlasticIceSheet plastic       plastic-static    1
PaleoEmulator   paleoemu      paleoemu-static   0
LIST
ref() { sed -n "s/^ *\"$1:\([^\"]*\)\".*/\1/p" kryos/.configme/manifest.toml; }
for e in "${pkgs[@]}"; do
    p=${e%%:*}; r=$(ref $p); [ -n "$r" ] || { echo "no $p in kryos/.configme/manifest.toml"; exit 1; }
    [ -d $p ] || git clone -q -b "$r" git@github.com:fesmc/$p.git $p
    printf "%-16s %-14s %s\n" $p "$r" "$(git -C $p log --oneline -1)"
done

echo "== $(date -Is) configure (configme: awi_albedo, ifx)"
for e in "${pkgs[@]}"; do
    IFS=: read -r p target link <<< "$e"
    [ "$link" = 1 ] && ln -sfn "$W/fesm-utils" "$W/$p/fesm-utils"
    [ -f "$W/$p/config/Makefile" ] && (cd $p && $CONFIGME config $p --only -m awi_albedo -c ifx > "$W/config-$p.log")
    ln -sfn "$W/$p" "$W/kryos/$p"
done
cd kryos
$CONFIGME config kryos --only -m awi_albedo -c ifx > "$W/config-kryos.log"
$CONFIGME deps . > /dev/null

echo "== $(date -Is) build"
make fesmutils-static openmp=$omp > build.log 2>&1
if on yelmo || on fastisostasy || on vilma || on plastic; then
    make -C fesm-utils libs openmp=$omp >> build.log 2>&1 || { tail -30 build.log; exit 1; }
fi
for e in "${pkgs[@]}"; do
    IFS=: read -r p target link <<< "$e"
    [ "$target" = - ] || make $target $BUILD >> build.log 2>&1 || { echo "$target failed:"; tail -30 build.log; exit 1; }
done
make kryos $BUILD -j8 >> build.log 2>&1 || { tail -30 build.log; exit 1; }
# VILMA2's standalone driver (vilma2/bin/vilma.x): the runs of tools/levante/gia_refcorr.sh (present-day correction, references)
if on vilma; then make -C vilma2 vilma openmp=$omp >> build.log 2>&1 || { tail -30 build.log; exit 1; }; fi
ls -l libkryos/bin/kryos.x
echo "== $(date -Is) tools/levante/setup.sh"
tools/levante/setup.sh || true
echo "== $(date -Is) done"
