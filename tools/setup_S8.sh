#!/bin/bash -l
# Stack S8 = S7 + fesm-utils phase0-kryos a13cd38 (ncio returns the values as stored: the underflow flush on read is gone; every dependent
# package is rebuilt) + vilma2 main 4419636 (the exact-restart tests with values below 1e-30 in the memory; the library is unchanged against
# d2006c8, rebuilt) + kryos main e4453dc (from GitHub, no bundle).
# yelmox R2' (topo-prep 5ac89c8, rembo=1) is rebuilt into $G/stacks/S8/bin/R2.
# kryos has two checkouts in the stack, both linked to S8 and built with the gate flags:
#   $S/kryos     vilma=1 (libvilma, SHTns)
#   $S/kryos-v0  vilma=0
# Also built in $S/kryos: libkryos/bin/hub_areas.x.
# As setup_S7.sh; the changes are the pins FU, VIL, KSHA, the stack name, the source stack S7, and the kryos clone from origin.
# Not repeatable: S8 must not exist. Run with: setsid nohup bash -l <this file> < /dev/null > /dev/null 2>&1 &
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate
S=$G/stacks/S8; R=$S/yelmox; GF="-O2 -fp-model precise -march=znver3 -traceback"
FU=a13cd38; VIL=4419636; YX=5ac89c8; KSHA=e4453dc
CM="env PYTHONPATH=$G/tools/configme-src/src python3 -m configme"
exec > >(tee "$G/stacks/s8-setup.log") 2>&1
echo "== $(date -Is) copy S7 -> S8"
[ -e $S ] && { echo "exists: $S"; exit 1; }
cp -a $G/stacks/S7 $S
chmod -R u+w $S/stack.lock 2>/dev/null || true
rm -f $S/stack.lock $S/build-gate*.log $S/vilma-link-test.* $S/build-vilma2*.log $S/check-*.log
rm -rf $S/vilma2 $S/bin/R2 $S/kryos $S/kryos-v0
echo "== $(date -Is) fesm-utils $FU"
cd $R/fesm-utils
git fetch -q origin phase0-kryos
git cat-file -e $FU^{commit}
git checkout -q $FU
git merge-base --is-ancestor $FU origin/phase0-kryos || { echo "$FU is not on origin/phase0-kryos"; exit 1; }
$CM config fesm-utils -m dkrz_levante -c ifx > $S/fesmutils-config.log 2>&1
echo "== $(date -Is) SHTns serial, omp (own flag scheme; links the stack's fftw)"
make shtns openmp=0 > $S/build-shtns-serial.log 2>&1
make shtns openmp=1 > $S/build-shtns-omp.log 2>&1
ls -l SHTns/shtns-serial/lib SHTns/shtns-omp/lib
echo "== $(date -Is) build R2'(S8): yelmox topo-prep $YX, rembo=1 (every package cleaned and rebuilt)"
cd $R
git checkout -q $YX
mkdir -p $S/bin/R2
MAKE_EXTRA="rembo=1" SKIP_LIBS=1 $G/tools/build_stack_wp0.sh $S
mv $S/build-gate-fortran.log $S/build-gate-R2.log
rm -rf $S/bin/R1 $S/bin/R2/yelmox.x
cp -p libyelmox/bin/yelmox.x $S/bin/R2/yelmox.x
sha256sum $S/bin/R2/yelmox.x
echo "== $(date -Is) fesm-utils serial library (rebuilt: the clean above removed it)"
make -C $R/fesm-utils fesmutils-static openmp=0 DFLAGS_NODEBUG="$GF" > $S/build-fesmutils-serial.log 2>&1
echo "== $(date -Is) vilma2 $VIL"
git clone -q git@github.com:fesmc/vilma2.git $S/vilma2
cd $S/vilma2
git checkout -q $VIL
ln -s $R/fesm-utils fesm-utils
$CM config vilma2 --only -m dkrz_levante -c ifx > $S/vilma2-config.log 2>&1
make vilma-static openmp=1 DFLAGS_NODEBUG="$GF" > $S/build-vilma2.log 2>&1
ls -l libvilma/include/libvilma.a
for p in . fesm-utils yelmo yelmo/FastHydrology yelmo/elsa yelmo/tracer FastIsostasy chion insol rembo1 ../vilma2; do
  printf "%-20s %s\n" "$p" $(git -C $R/$p rev-parse HEAD); done > $S/stack.lock
for V in 1 0; do
  D=$S/kryos; [ $V = 0 ] && D=$S/kryos-v0
  echo "== $(date -Is) kryos $KSHA on S8, vilma=$V ($D)"
  rm -rf $D
  git clone -q git@github.com:fesmc/kryos.git $D
  git -C $D checkout -q --detach $KSHA
  for p in fesm-utils yelmo FastIsostasy chion insol; do ln -sfn $R/$p $D/$p; done
  ln -sfn $S/vilma2 $D/vilma2
  for d in ice_data fesmdata isostasy_data; do ln -sfn /work/ba1442/$d $D/$d; done
  (cd $D && $CM config kryos --only -m dkrz_levante -c ifx > /dev/null && $CM deps . > /dev/null)
  (cd $D && make kryos openmp=1 vilma=$V DFLAGS_NODEBUG="$GF" > $D/build.log 2>&1) || { tail -30 $D/build.log; exit 1; }
  ls -l $D/libkryos/bin/kryos.x
  sha256sum $D/libkryos/bin/kryos.x
  echo $KSHA-vilma$V > $D/libkryos/bin/.kryos_sha
  [ $V = 1 ] && (cd $D && make libkryos/bin/hub_areas.x openmp=1 vilma=1 DFLAGS_NODEBUG="$GF" > $D/build-tools.log 2>&1)
done
echo "== $(date -Is) done"
