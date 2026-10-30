#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=ba1442
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --job-name=cmp-S11
#SBATCH --output=/work/ba1442/robinson/kryos-gate/reports/S11/cmp-S11-%j.log
# sbatch --export=ALL compare_S11.sh g2b|g4|g3|g6|g2c; reports and summary.txt in $G/reports/S11 (scripts of the stack checkout).
#   g2b: g2b_compare.jl over l085 l128 l170 (and the S9 diagnosis runs for the before/after table), g2b_diag.jl for l085.
#   g4:  sci-vilma2 cont2 == cont, chunk == cont, rounding (cmp_g4.sh).
#   g3:  g3_budget.jl on sci-vilma2 (cont with its lossless twin) and on the g6 cont (lossless itself).
#   g6:  chunk == cont for g6 and g6-smb (kryos-gate-kk.jl, bundles included), drift against K3-sci (g6_report.jl).
#   g2c: g2c_report.jl on the S11 spin-up run.
set -uo pipefail
module load python3 >/dev/null 2>&1 || true
export PATH=$HOME/.juliaup/bin:$HOME/.local/bin:$PATH
G=/work/ba1442/robinson/kryos-gate; T=$G/p4w7/kryos/tools; S=$G/stacks/S11; K=$S/kryos; L=$K/tools/levante
KSHA=${KSHA:-$(cut -d- -f1 $K/libkryos/bin/.kryos_sha)}   # the commit the stack's kryos was built at (git is not on the compute nodes' PATH)
S11T=$G/s11/tools                                       # scripts fixed after the stack was set up (g3_budget.jl, g2a_compare.jl)
export JULIA_DEPOT_PATH=$G/julia-depot
OUT=$G/reports/S11; mkdir -p $OUT
J="julia +1.13 --project=$T"
KK="$J $T/bin/kryos-gate-kk.jl check"
K3=$G/ref/K3-sci/kx-2d80b43/S5b/sci/cont
case ${1:-} in
  g2b)
    B=$G/runs/g2b/kx-$KSHA/S11
    $J $L/g2b_compare.jl $B l085 l128 l170 --md $OUT/g2b.md --csv $OUT/g2b-csv > $OUT/g2b.log 2>&1
    echo "g2b | exit $? | g2b.md" >> $OUT/summary.txt
    $J $L/g2b_diag.jl $B l085 $G/p4g2b/maps $OUT/g2b-diag > $OUT/g2b-diag.log 2>&1
    echo "g2b diag l085 | exit $? | g2b-diag/" >> $OUT/summary.txt ;;
  g4)
    SHA=$KSHA STK=S11 bash $L/cmp_g4.sh > $OUT/g4-vilma2.log 2>&1
    echo "sci-vilma2 | $(sed -n 1p $OUT/sci-vilma2-kk.md) | g4-vilma2.log" >> $OUT/summary.txt ;;
  g3)
    R=$G/runs/sci-vilma2/kx-$KSHA/S11
    $J $S11T/g3_budget.jl $R/cont $R/lossless $K/libkryos/bin/hub_areas.x --md $OUT/g3-sci-vilma2.md > $OUT/g3-sci-vilma2.log 2>&1
    echo "g3 sci-vilma2 | exit $? | g3-sci-vilma2.md" >> $OUT/summary.txt
    R=$G/runs/g6/kx-$KSHA/S11
    $J $S11T/g3_budget.jl $R/cont $R/cont $K/libkryos/bin/hub_areas.x --md $OUT/g3-g6.md > $OUT/g3-g6.log 2>&1
    echo "g3 g6 | exit $? | g3-g6.md" >> $OUT/summary.txt ;;
  g6)
    for g in g6 g6-smb; do
        $KK $K/tools/gate/$g.toml $G/runs/$g/kx-$KSHA/S11 --md $OUT/$g-kk.md > $OUT/$g-kk.log 2>&1
        echo "$g chunk == cont | exit $? | $(sed -n 1p $OUT/$g-kk.md)" >> $OUT/summary.txt
        $J $L/g6_report.jl $G/runs/$g/kx-$KSHA/S11/cont $K3 --md $OUT/$g-drift.md > $OUT/$g-drift.log 2>&1
        echo "$g drift | exit $? | $g-drift.md" >> $OUT/summary.txt
    done ;;
  g2c)
    $J $L/g2c_report.jl $G/runs/g2c/kx-$KSHA/S11 --md $OUT/g2c.md > $OUT/g2c.log 2>&1
    echo "g2c | exit $? | g2c.md" >> $OUT/summary.txt ;;
  *) echo "usage: compare_S11.sh g2b|g4|g3|g6|g2c"; exit 1 ;;
esac
echo "== $(date -Is) done $1" >> $OUT/summary.txt
