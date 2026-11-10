#!/bin/bash
#SBATCH --partition=compute
#SBATCH --account=ba1442
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --job-name=kryos-check
# G0: `make check` of a kryos checkout of stack S6 with the gate flags. sbatch --export=ALL,V=1|0 check_kryos.sh
# TGT: the make targets (default check; S10: "check check-long").
# V=1: $S/kryos (vilma=1), V=0: $S/kryos-v0 (vilma=0). Log: $S/check-<LOGTAG, default vilma<V>>.log
set -uo pipefail
module load python3 >/dev/null 2>&1 || true
export PATH=$HOME/.juliaup/bin:$HOME/.local/bin:$PATH
G=/work/ba1442/robinson/kryos-gate; S=$G/stacks/${STK:-S6}; GF="-O2 -fp-model precise -march=znver3 -traceback"
export JULIA_DEPOT_PATH=$G/julia-depot KRYOS_TOOLS=$G/p4w7/kryos/tools OMP_NUM_THREADS=16 OMP_STACKSIZE=256M
ulimit -s unlimited
D=$S/kryos; [ "${V:-1}" = 0 ] && D=$S/kryos-v0
D=${CK:-$D}; LT=${LOGTAG:-vilma${V:-1}}      # CK: another checkout (the re-check of the test fixes), LOGTAG: its log name
cd $D
git rev-parse HEAD
make -k ${TGT:-check} openmp=1 vilma=${V:-1} DFLAGS_NODEBUG="$GF" > $S/check-$LT.log 2>&1
rc=$?
echo "make ${TGT:-check} vilma=${V:-1} exit $rc"; tail -4 $S/check-$LT.log
if [ $rc -ne 0 ]; then    # make -k built what it could: run every test executable that exists, one line each
    echo "== individual runs of the built tests"
    for t in libkryos/bin/test_*.x; do
        ./$t > $S/check-$LT-$(basename $t .x).log 2>&1; echo "$(basename $t .x) exit $?"
    done
fi
exit $rc
