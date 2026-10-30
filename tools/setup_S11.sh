#!/bin/bash -l
# Stack S11 = S10 + vilma2 host-orient (the host interface owns the row order: the passthrough of a Gauss host with ascending
# latitude reverses the rows both ways, host-grid ocean function se%C; the legacy Gauss-grid reader takes the row order from the
# file's lat axis) + kryos vilma-orient (the adapter exports VILMA2's host-grid outputs, never se%gg), both from bundles (not
# pushed): $G/bundles/vilma2-host-orient.bundle (VBRANCH), $G/bundles/kryos-vilma-orient.bundle (KBRANCH).
# fesm-utils, yelmo and the other packages as S10 (afae7a0; not rebuilt). Every vilma2 restart written before S11 holds mirrored
# state and is invalid (vilma_restart_read stops on its reference unless the reference is north-south symmetric).
# kryos checkouts, linked to S11 and built with the gate flags:
#   $S/kryos     vilma=1 (libvilma, SHTns); also libkryos/bin/hub_areas.x and g2a_cap.x
#   $S/kryos-v0  vilma=0
# As setup_S10.sh without the fesm-utils, SHTns and yelmox steps and the backward-compatibility variants; the changes are the
# vilma2 step (bundle), KBRANCH, the stack name and S10.
# Not repeatable: S11 must not exist. Run as a batch job: sbatch --account=ba1442 --partition=shared --time=02:00:00 --cpus-per-task=8 --wrap="bash -l <abs path of this file>"
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate
S=$G/stacks/S11; R=$S/yelmox; GF="-O2 -fp-model precise -march=znver3 -traceback"
VIL=c17d513; VBRANCH=host-orient; VB=$G/bundles/vilma2-host-orient.bundle
KBRANCH=vilma-orient; KB=$G/bundles/kryos-vilma-orient.bundle
CM="env PYTHONPATH=$G/tools/configme-src/src python3 -m configme"
exec > >(tee "$G/stacks/s11-setup.log") 2>&1
echo "== $(date -Is) copy S10 -> S11"
[ -e $S ] && { echo "exists: $S"; exit 1; }
cp -a $G/stacks/S10 $S
chmod -R u+w $S/stack.lock 2>/dev/null || true
rm -f $S/stack.lock $S/vilma-link-test.* $S/build-vilma2*.log $S/check-*.log
rm -rf $S/vilma2 $S/kryos $S/kryos-v0 $S/kryos-x $S/kryos-d
echo "== $(date -Is) vilma2 $VIL ($VBRANCH)"
git clone -q -b $VBRANCH $VB $S/vilma2
cd $S/vilma2
git cat-file -e $VIL^{commit}
git checkout -q $VIL
ln -s $R/fesm-utils fesm-utils
$CM config vilma2 --only -m dkrz_levante -c ifx > $S/vilma2-config.log 2>&1
make vilma-static openmp=1 DFLAGS_NODEBUG="$GF" > $S/build-vilma2.log 2>&1
ls -l libvilma/include/libvilma.a
for p in . fesm-utils yelmo yelmo/FastHydrology yelmo/elsa yelmo/tracer FastIsostasy chion insol rembo1 ../vilma2; do
  printf "%-20s %s\n" "$p" $(git -C $R/$p rev-parse HEAD); done > $S/stack.lock
for C in kryos kryos-v0; do
  D=$S/$C; V=1; [ $C = kryos-v0 ] && V=0
  echo "== $(date -Is) kryos $KBRANCH on S11, vilma=$V ($D)"
  git clone -q -b $KBRANCH $KB $D
  KSHA=$(git -C $D rev-parse --short=7 HEAD)
  for p in fesm-utils yelmo FastIsostasy chion insol; do ln -sfn $R/$p $D/$p; done
  ln -sfn $S/vilma2 $D/vilma2
  for d in ice_data fesmdata isostasy_data; do ln -sfn /work/ba1442/$d $D/$d; done
  (cd $D && $CM config kryos --only -m dkrz_levante -c ifx > /dev/null && $CM deps . > /dev/null)
  (cd $D && make kryos openmp=1 vilma=$V DFLAGS_NODEBUG="$GF" > $D/build.log 2>&1) || { tail -30 $D/build.log; exit 1; }
  ls -l $D/libkryos/bin/kryos.x
  sha256sum $D/libkryos/bin/kryos.x
  echo $KSHA-vilma$V > $D/libkryos/bin/.kryos_sha
  [ $C = kryos ] && (cd $D && make libkryos/bin/hub_areas.x libkryos/bin/g2a_cap.x openmp=1 vilma=1 DFLAGS_NODEBUG="$GF" > $D/build-tools.log 2>&1)
done
echo "== $(date -Is) done"
