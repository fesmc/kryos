#!/bin/bash -l
# stack.sh: a kryos stack on AWI's albedo, from scratch: configme install kryos (kryos and every dependency at the refs
# of its .configme/manifest.toml, configured for albedo: awi_albedo, ifx; fesm-utils and its external libraries built),
# make kryos with the backends of BUILD (it builds their libraries first), VILMA2's standalone driver; then
# tools/levante/setup.sh (the data links, maps, runme, the Julia project of the tools).
#
#   tools/albedo/stack.sh <dir> [<kryos ref>]      (from a kryos checkout; or curl it and run it anywhere)
#
#   <dir>        a new folder, e.g. /albedo/work/user/$USER/kryos-plastic; it gets kryos/ (<kryos ref>, default main) with
#                the dependency checkouts inside. An existing install is updated by configme, not replaced.
#   BUILD        the make flags of the backends (default "yelmo=0 fastisostasy=0 vilma=1 plastic=1 openmp=1": the plastic
#                ice backend with VILMA2)
#   CONFIGME     the configme command (configme >= 0.16.1: vilma2 and PlasticIceSheet in the kryos set; default configme)
#
# The data come from the shared stores of p_forclima (tools/levante/env.sh): FesmData ($FESMDATA, kept up to date
# with FesmData's `julia fesmdata.jl mirror`, shared/machines/albedo.env) and kryos_data ($KRYOS_DATA, a copy of
# Levante's /work/ba1442/kryos_data: maps, paleoemu, plastic, refcorr, paleomist).
set -euo pipefail
[ $# -ge 1 ] || { sed -n '2,13p' "$0"; exit 1; }
W=$1; KREF=${2:-main}
BUILD=${BUILD:-"yelmo=0 fastisostasy=0 vilma=1 plastic=1 openmp=1"}
CONFIGME=${CONFIGME:-configme}
omp=0; [[ " $BUILD " == *" openmp=1 "* ]] && omp=1
command -v git > /dev/null 2>&1 || module load git > /dev/null 2>&1
mkdir -p "$W"; cd "$W"
exec > >(tee -a "$W/stack.log") 2>&1
echo "== $(date -Is) configme install kryos:$KREF in $W"
$CONFIGME install kryos:$KREF -m awi_albedo -c ifx --build-deps > "$W/install.log" 2>&1 || { tail -30 "$W/install.log"; exit 1; }
cd kryos
for p in kryos fesm-utils vilma2 PlasticIceSheet yelmo FastIsostasy chion insol; do
    d=.; [ $p = kryos ] || d=$p
    printf "%-16s %s\n" $p "$(git -C $d log --oneline -1)"
done

echo "== $(date -Is) build ($BUILD)"
make kryos $BUILD -j8 > build.log 2>&1 || { tail -30 build.log; exit 1; }
# VILMA2's standalone driver (vilma2/bin/vilma.x): the runs of tools/levante/gia_refcorr.sh (present-day correction, references)
if [[ " $BUILD " == *" vilma=1 "* ]]; then make -C vilma2 vilma openmp=$omp >> build.log 2>&1 || { tail -30 build.log; exit 1; }; fi
ls -l libkryos/bin/kryos.x
echo "== $(date -Is) tools/levante/setup.sh"
tools/levante/setup.sh || true
echo "== $(date -Is) done"
