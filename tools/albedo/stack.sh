#!/bin/bash -l
# stack.sh: a kryos stack on AWI's albedo, from scratch: the kryos checkout and the dependencies at the refs of its
# .configme/manifest.toml, configured for albedo (configme, awi_albedo, ifx), linked and built; then
# tools/levante/setup.sh (the data links, maps, runme, the Julia project of the tools).
#
#   tools/albedo/stack.sh <dir> [<kryos ref>]      (from a kryos checkout; or curl it and run it anywhere)
#
#   <dir>        a new folder, e.g. /albedo/work/user/$USER/kryos-work; it gets kryos/ (the checkout, <kryos ref>, default
#                main) and fesm-utils/ chion/ insol/ PaleoEmulator/ (the manifest's refs; yelmo, FastIsostasy and vilma2
#                are left out: the build is yelmo=0 fastisostasy=0 vilma=0 plastic=0)
#   BUILD        the make flags of the backends (default "yelmo=0 fastisostasy=0 vilma=0 plastic=0 paleoemu=1 openmp=1")
#   CONFIGME     the configme command, one that knows the kryos orchestrator (configme >= 0.16; default configme).
#                Before its release: CONFIGME="env PYTHONPATH=<configme checkout>/src python3 -m configme"
#
# The data come from the shared stores of p_forclima (tools/levante/env.sh): FesmData ($FESMDATA, kept up to date
# with FesmData's `julia fesmdata.jl mirror`, shared/machines/albedo.env) and kryos_data ($KRYOS_DATA, a copy of
# Levante's /work/ba1442/kryos_data: maps, paleoemu, plastic, refcorr, paleomist).
# fesm-utils is built serially (its Makefile has no module dependencies: -j races).
set -euo pipefail
[ $# -ge 1 ] || { sed -n '2,8p' "$0"; exit 1; }
W=$1; KREF=${2:-main}
BUILD=${BUILD:-"yelmo=0 fastisostasy=0 vilma=0 plastic=0 paleoemu=1 openmp=1"}
CONFIGME=${CONFIGME:-configme}
command -v git > /dev/null 2>&1 || module load git > /dev/null 2>&1
mkdir -p "$W"; cd "$W"
exec > >(tee -a "$W/stack.log") 2>&1
echo "== $(date -Is) kryos $KREF in $W"
[ -d kryos ] || git clone -q -b "$KREF" git@github.com:fesmc/kryos.git kryos
ref() { sed -n "s/^ *\"$1:\([^\"]*\)\".*/\1/p" kryos/.configme/manifest.toml; }
for p in fesm-utils chion insol PaleoEmulator; do
    r=$(ref $p); [ -n "$r" ] || { echo "no $p in kryos/.configme/manifest.toml"; exit 1; }
    [ -d $p ] || git clone -q -b "$r" git@github.com:fesmc/$p.git $p
    printf "%-14s %-14s %s\n" $p "$r" "$(git -C $p log --oneline -1)"
done
for p in fesm-utils chion insol; do (cd $p && $CONFIGME config $p --only -m awi_albedo -c ifx > "$W/config-$p.log"); done
for p in chion insol; do ln -sfn "$W/fesm-utils" "$W/$p/fesm-utils"; done
for p in fesm-utils chion insol PaleoEmulator; do ln -sfn "$W/$p" "$W/kryos/$p"; done
cd kryos
$CONFIGME config kryos --only -m awi_albedo -c ifx > "$W/config-kryos.log"
$CONFIGME deps . > /dev/null
echo "== $(date -Is) build ($BUILD)"
make fesmutils-static $BUILD > build.log 2>&1
make chion-static $BUILD -j8 >> build.log 2>&1
make insol-static $BUILD >> build.log 2>&1
make paleoemu-static $BUILD >> build.log 2>&1
make kryos $BUILD -j8 >> build.log 2>&1 || { tail -30 build.log; exit 1; }
ls -l libkryos/bin/kryos.x
echo "== $(date -Is) tools/levante/setup.sh"
tools/levante/setup.sh || true
echo "== $(date -Is) done"
