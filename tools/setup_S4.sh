#!/bin/bash -l
# Stack S4 = S3 + rembo=1 (rembo1 v1.3 linked into yelmox) + Yelmo vm-m16 restart fix
# (yelmo topo-prepared ea90fee7: eps_eig_1/2 recomputed on restart).
# Builds both yelmox binaries of the phase 2 baselines (gate flags, rembo=1) and keeps them:
#   $G/stacks/S4/bin/R1/yelmox.x   yelmox v2.4 (918d996): R1
#   $G/stacks/S4/bin/R2/yelmox.x   yelmox topo-prep 5ac89c8 (R2' = 74f3a2b + [topo] keys of the other 8 bases)
# runme copies the exe at staging, so runs are staged from these copies (tools/levante/submit_wp0.sh).
# Then kryos is relinked to S4 (code unchanged) and rebuilt (clean: yelmo objects changed).
# Bundles (in $G/bundles): ym4.bundle (yelmo topo-prepared, b5ab5d2c..ea90fee7), yx2.bundle (yelmox 918d996..topo-prep).
set -euo pipefail
G=/work/ba1442/robinson/kryos-gate; S=$G/stacks/S4; R=$S/yelmox; GF="-O2 -fp-model precise -march=znver3 -traceback"
exec > >(tee "$G/stacks/s4-setup.log") 2>&1
echo "== $(date -Is) copy S3 -> S4"
[ -e $S ] && { echo "exists: $S"; exit 1; }
cp -a $G/stacks/S3 $S
chmod -R u+w $S/stack.lock 2>/dev/null || true; rm -f $S/stack.lock $S/build-gate*.log
cd $R/yelmo
git fetch -q $G/bundles/ym4.bundle topo-prepared:refs/remotes/b4/topo-prepared
git checkout -q ea90fee7
cd $R
git fetch -q $G/bundles/yx2.bundle topo-prep:refs/remotes/b4/topo-prep
git rev-parse --short=7 refs/remotes/b4/topo-prep
mkdir -p $S/bin/R1 $S/bin/R2
echo "== $(date -Is) build R1(S4): yelmox v2.4, rembo=1"
git checkout -q 918d996
MAKE_EXTRA="rembo=1" SKIP_LIBS=1 $G/tools/build_stack_wp0.sh $S
mv $S/build-gate-fortran.log $S/build-gate-R1.log
cp -p libyelmox/bin/yelmox.x $S/bin/R1/yelmox.x
echo "== $(date -Is) build R2'(S4): yelmox topo-prep, rembo=1"
git checkout -q 5ac89c8
make clean > /dev/null
make yelmox openmp=1 rembo=1 DFLAGS_NODEBUG="$GF" > $S/build-gate-R2.log 2>&1
cp -p libyelmox/bin/yelmox.x $S/bin/R2/yelmox.x
sha256sum $S/bin/R*/yelmox.x
for p in . fesm-utils yelmo yelmo/FastHydrology yelmo/elsa yelmo/tracer FastIsostasy chion insol rembo1; do
  printf "%-20s %s\n" "$p" $(git -C $R/$p rev-parse HEAD); done > $S/stack.lock
echo "linked rembo: $(grep -c 'rembo' $S/build-gate-R2.log) lines; climate_rembo.o: $(ls $R/libyelmox/include/climate_rembo.o 2>&1 | tail -1)"
echo "== $(date -Is) kryos on S4"
K=$G/kryos
for p in fesm-utils yelmo FastIsostasy chion insol; do ln -sfn $S/yelmox/$p $K/$p; done
make -C $K clean > /dev/null 2>&1 || true
echo "== $(date -Is) done"
