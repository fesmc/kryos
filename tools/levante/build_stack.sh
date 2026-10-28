#!/bin/bash -l
# Gate build of one yelmox stack: every package cleaned, then rebuilt with the
# gate flags given on the make command line (they reach the recursive builds).
#   build_stack.sh <stack-dir>     e.g. $G/stacks/v24
# Environment: SKIP_LIBS=1 skips the fesm-utils externals; MAKE_EXTRA="rembo=1" adds make arguments.
set -euo pipefail
S=$1
R=$S/yelmox
GF="-O2 -fp-model precise -march=znver3 -traceback"
LOG=$S/build-gate${SKIP_LIBS:+-fortran}.log
exec > >(tee "$LOG") 2>&1
module list 2>&1
echo "== $(date -Is) gate build of $R"
echo "DFLAGS_NODEBUG=$GF  MAKE_EXTRA=${MAKE_EXTRA:-}"

# fesm-utils external libraries (omp variant, the one yelmox openmp=1 links):
# their flags are CFLAGS_LIB/FFLAGS_LIB of config/libs.mk, not DFLAGS_NODEBUG.
cd "$R/fesm-utils"
if [ -z "${SKIP_LIBS:-}" ]; then
make -j16 fftw openmp=1 CFLAGS_LIB="$GF -fPIC" FFLAGS_LIB="$GF -fPIC"
make -j16 lis openmp=1 CFLAGS_LIB="$GF -fPIC" FFLAGS_LIB="$GF -fPIC"
fi
echo "== $(hostname) $(date -Is) Fortran builds"

# Fortran objects of every package
for d in fesm-utils yelmo yelmo/FastHydrology yelmo/elsa yelmo/tracer FastIsostasy chion insol rembo1 .; do
    echo "== make clean in $d"
    make -C "$R/$d" clean || echo "   (no clean target / failed in $d)"
done

cd "$R"
make yelmox openmp=1 ${MAKE_EXTRA:-} DFLAGS_NODEBUG="$GF"
ls -l libyelmox/bin/yelmox.x
echo "== $(date -Is) done"
