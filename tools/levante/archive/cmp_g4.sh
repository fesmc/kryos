#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=ba1442
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --job-name=cmp-g4
#SBATCH --output=/work/ba1442/robinson/kryos-gate/reports/cmp-g4-%j.log
# G4: sci-vilma2 (cont, cont2, chunk): cont2 == cont and chunk == cont, bundles included (kryos-gate-kk.jl check);
# rounded cont against the lossless run: bundles identical, every rounded file within its table tolerance (round_check.jl).
set -uo pipefail
module load python3 >/dev/null 2>&1 || true
export PATH=$HOME/.juliaup/bin:$HOME/.local/bin:$PATH JULIA_DEPOT_PATH=/work/ba1442/robinson/kryos-gate/julia-depot
G=/work/ba1442/robinson/kryos-gate; T=$G/p4w7/kryos/tools; SHA=${SHA:-cd85c14}; STK=${STK:-S6}; OUT=$G/reports/$STK; mkdir -p $OUT
R=$G/runs/sci-vilma2/kx-$SHA/$STK
echo "== cont2, chunk"
julia +1.13 --project=$T $T/bin/kryos-gate-kk.jl check sci-vilma2 $R --md $OUT/sci-vilma2-kk.md > $OUT/sci-vilma2-kk.log 2>&1
echo "exit $?"; sed -n 1,30p $OUT/sci-vilma2-kk.md
echo "== rounding"
julia +1.13 --project=$T $G/p4w7/levante/round_check.jl $R/cont $R/lossless $G/stacks/$STK/kryos/input --md $OUT/sci-vilma2-round.md 2>&1 | tail -60
