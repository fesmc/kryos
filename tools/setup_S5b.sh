#!/bin/bash -l
# Stack S5b = S5 + fesm-utils phase0-kryos 3ba56a9 (sphere radius override, "years ago" time units, coupler_merge_end log threshold;
# a coords/varslice/coupler change, so every dependent package is rebuilt) + vilma2 main 2676879 (libvilma.a, openmp variant),
# rebuilt on it. yelmox R2' (topo-prep 5ac89c8, rembo=1) is rebuilt into $G/stacks/S5b/bin/R2 for the R2(S5b) == R2(S5) check.
# kryos has its own checkout $S/kryos (at KSHA, linked to S5b) so that $G/kryos and the S5 links are not touched.
# As setup_S5.sh; the changes to S5 are the pins FU, VIL, KSHA, the stack name and the kryos checkout.
# Not repeatable: S5b must not exist. Run with: setsid nohup bash -l <this file> < /dev/null > /dev/null 2>&1 &
set -euo pipefail
G=/work/ba1442/kryos-gate; G=/work/ba1442/robinson/kryos-gate
S=$G/stacks/S5b; R=$S/yelmox; GF="-O2 -fp-model precise -march=znver3 -traceback"
FU=3ba56a9; VIL=2676879; YX=5ac89c8; KSHA=5a663a3
CM="env PYTHONPATH=$G/tools/configme-src/src python3 -m configme"
exec > >(tee "$G/stacks/s5b-setup.log") 2>&1
echo "== $(date -Is) copy S5 -> S5b"
[ -e $S ] && { echo "exists: $S"; exit 1; }
cp -a $G/stacks/S5 $S
chmod -R u+w $S/stack.lock 2>/dev/null || true; rm -f $S/stack.lock $S/build-gate*.log $S/vilma-link-test.* $S/build-vilma2*.log
rm -rf $S/vilma2 $S/bin/R2
echo "== $(date -Is) fesm-utils $FU"
cd $R/fesm-utils
git fetch -q origin phase0-kryos
git cat-file -e $FU^{commit}
git checkout -q $FU
git merge-base --is-ancestor $FU origin/phase0-kryos || { echo "$FU is not on origin/phase0-kryos"; exit 1; }
configme config fesm-utils -m dkrz_levante -c ifx > $S/fesmutils-config.log 2>&1
echo "== $(date -Is) SHTns serial, omp (own flag scheme; links the stack's fftw)"
make shtns openmp=0 > $S/build-shtns-serial.log 2>&1
make shtns openmp=1 > $S/build-shtns-omp.log 2>&1
ls -l SHTns/shtns-serial/lib SHTns/shtns-omp/lib
echo "== $(date -Is) build R2'(S5b): yelmox topo-prep $YX, rembo=1 (every package cleaned and rebuilt)"
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
echo "== $(date -Is) kryos $KSHA on S5b ($S/kryos; built by tools/levante/run_kryos.sh with K=$S/kryos KSHA=$KSHA)"
git clone -q $G/kryos $S/kryos
git -C $S/kryos checkout -q --detach $KSHA
for p in fesm-utils yelmo FastIsostasy chion insol; do ln -sfn $R/$p $S/kryos/$p; done
ln -sfn $S/vilma2 $S/kryos/vilma2
for d in ice_data fesmdata isostasy_data; do ln -sfn /work/ba1442/$d $S/kryos/$d; done
(cd $S/kryos && $CM config kryos --only -m dkrz_levante -c ifx > /dev/null && $CM deps . > /dev/null)
echo "== $(date -Is) done"
