#!/bin/bash -l
# Stack S10 = S9 + fesm-utils fields-diag afae7a0 (fields_add restart = .false.: diagnostics left out of the restart files, ignored
# and logged when an older file holds one; field_t changes, so every dependent package is rebuilt, yelmo included) + kryos rst-state
# (state-only restore: the diagnostics declared and put again on PASS_RESTART), both from bundles (not pushed):
# $G/bundles/fu-fields-diag.bundle (45e09c2..fields-diag), $G/bundles/kryos-rst-state.bundle (KBRANCH). vilma2 main e449d58 (S9: 4419636;
# bsl_ive_host, which kryos main exports since f4ec006).
# yelmox R2' (topo-prep 5ac89c8, rembo=1) is rebuilt into $G/stacks/S10/bin/R2 (the step that rebuilds every package).
# kryos checkouts, all linked to S10 and built with the gate flags:
#   $S/kryos     vilma=1 (libvilma, SHTns); also libkryos/bin/hub_areas.x
#   $S/kryos-v0  vilma=0
#   $S/kryos-x   vilma=1, + one extra diagnostic export (bsl_ive_x), local commit of p4-s10-diag-extra.diff (backward-compatibility test)
#   $S/kryos-d   vilma=1, - the diagnostic bsl_ive_host, local commit of p4-s10-diag-drop.diff (backward-compatibility test)
# As setup_S9.sh; the changes are the pins FU (bundle), no yelmo step (30361f96 as S9), KBRANCH (bundle), the variants, the stack name and S9.
# Not repeatable: S10 must not exist. Run as a batch job: sbatch --account=ba1442 --partition=shared --time=02:00:00 --cpus-per-task=8 --wrap="bash -l <abs path of this file>"
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate
S=$G/stacks/S10; R=$S/yelmox; GF="-O2 -fp-model precise -march=znver3 -traceback"
FU=afae7a0; FB=$G/bundles/fu-fields-diag.bundle; VIL=e449d58; YX=5ac89c8; KBRANCH=rst-state; KB=$G/bundles/kryos-rst-state.bundle
PD=$G/rst10            # the variant patches (docs/dev/gate-reports/p4-s10-diag-*.diff, copied there)
CM="env PYTHONPATH=$G/tools/configme-src/src python3 -m configme"
exec > >(tee "$G/stacks/s10-setup.log") 2>&1
echo "== $(date -Is) copy S9 -> S10"
[ -e $S ] && { echo "exists: $S"; exit 1; }
cp -a $G/stacks/S9 $S
chmod -R u+w $S/stack.lock 2>/dev/null || true
rm -f $S/stack.lock $S/build-gate*.log $S/vilma-link-test.* $S/build-vilma2*.log $S/check-*.log
rm -rf $S/vilma2 $S/bin/R2 $S/kryos $S/kryos-v0
echo "== $(date -Is) fesm-utils $FU"
cd $R/fesm-utils
git fetch -q $FB fields-diag:refs/remotes/b10/fields-diag
git cat-file -e $FU^{commit}
git checkout -q $FU
git merge-base --is-ancestor $FU b10/fields-diag || { echo "$FU is not on fields-diag"; exit 1; }
$CM config fesm-utils -m dkrz_levante -c ifx > $S/fesmutils-config.log 2>&1
echo "== $(date -Is) SHTns serial, omp (own flag scheme; links the stack's fftw)"
make shtns openmp=0 > $S/build-shtns-serial.log 2>&1
make shtns openmp=1 > $S/build-shtns-omp.log 2>&1
ls -l SHTns/shtns-serial/lib SHTns/shtns-omp/lib
echo "== $(date -Is) build R2'(S10): yelmox topo-prep $YX, rembo=1 (every package cleaned and rebuilt)"
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
for C in kryos kryos-v0 kryos-x kryos-d; do
  D=$S/$C; V=1; [ $C = kryos-v0 ] && V=0
  echo "== $(date -Is) kryos $KBRANCH on S10, vilma=$V ($D)"
  rm -rf $D
  git clone -q -b $KBRANCH $KB $D
  case $C in
    kryos-x) git -C $D apply $PD/p4-s10-diag-extra.diff; git -C $D commit -q -am "S10 test variant: extra diagnostic bsl_ive_x" ;;
    kryos-d) git -C $D apply $PD/p4-s10-diag-drop.diff;  git -C $D commit -q -am "S10 test variant: diagnostic bsl_ive_host dropped" ;;
  esac
  KSHA=$(git -C $D rev-parse --short=7 HEAD)
  for p in fesm-utils yelmo FastIsostasy chion insol; do ln -sfn $R/$p $D/$p; done
  ln -sfn $S/vilma2 $D/vilma2
  for d in ice_data fesmdata isostasy_data; do ln -sfn /work/ba1442/$d $D/$d; done
  (cd $D && $CM config kryos --only -m dkrz_levante -c ifx > /dev/null && $CM deps . > /dev/null)
  (cd $D && make kryos openmp=1 vilma=$V DFLAGS_NODEBUG="$GF" > $D/build.log 2>&1) || { tail -30 $D/build.log; exit 1; }
  ls -l $D/libkryos/bin/kryos.x
  sha256sum $D/libkryos/bin/kryos.x
  echo $KSHA-vilma$V > $D/libkryos/bin/.kryos_sha
  [ $C = kryos ] && (cd $D && make libkryos/bin/hub_areas.x openmp=1 vilma=1 DFLAGS_NODEBUG="$GF" > $D/build-tools.log 2>&1)
done
echo "== $(date -Is) done"
