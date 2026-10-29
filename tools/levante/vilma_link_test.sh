#!/bin/bash -l
# Link test of libvilma.a on Levante (phase 4 WP0): vilma2's own host-API tests compile and link against the archive
# with the host link line (-lvilma, -lfesmutils, SHTns, FFTW, netCDF), gate flags, openmp variant of the stack.
#   vilma_link_test.sh <stack-dir>        e.g. $G/stacks/S5      (log: <stack-dir>/vilma-link-test.log, job id on stdout)
# Builds on the login node, runs in a short job (16 threads, 30 min, account ba1442). Tests: params (defaults file),
# host_terms, host_bitident, restart_host (host API), io_prec (output precision), remap (Gauss grid from fesm-utils grid_init_gaussian),
# coupling (lmax 16, 3D viscosity of input/bagge2021.nc) and test_visc_fesmdata (the FesmData Bagge2021_visc.nc, read only).
set -euo pipefail
S=$1; V=$S/vilma2; GF="-O2 -fp-model precise -march=znver3 -traceback"
VISC=${VISC:-/work/ba1442/fesmdata-prepared/Bagge2021_visc/Bagge2021_visc.nc}
REFV=${REFV:-$V/input/bagge2021.nc}      # vendored in vilma2 (isostasy_data has no bagge2021.nc)
T="test_params test_host_terms test_host_bitident test_restart_host test_remap test_coupling test_io_prec"
cd $V
make $T test_visc_fesmdata openmp=1 DFLAGS_NODEBUG="$GF" > $S/build-vilma2-tests.log 2>&1 || { tail -30 $S/build-vilma2-tests.log; exit 1; }
cat > $S/vilma-link-test.job <<JOB
#!/bin/bash -l
#SBATCH --job-name=vilma-link
#SBATCH --account=ba1442
#SBATCH --partition=compute
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=$S/vilma-link-test.log
#SBATCH --error=$S/vilma-link-test.log
export OMP_NUM_THREADS=16 OMP_STACKSIZE=256M
ulimit -s unlimited
cd $V
echo "libvilma.a: \$(ls -l libvilma/include/libvilma.a)"
for t in $T; do echo "--- \$t ---"; ./bin/\$t.x || { echo "FAILED \$t"; exit 1; }; done
echo "--- test_visc_fesmdata (read only: $VISC, $REFV) ---"
./bin/test_visc_fesmdata.x $VISC $REFV 6 || echo "FAILED test_visc_fesmdata"
echo "== all link tests done"
JOB
sbatch --parsable $S/vilma-link-test.job
