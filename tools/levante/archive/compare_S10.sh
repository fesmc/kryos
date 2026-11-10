#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=ba1442
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --job-name=cmp-S10
#SBATCH --output=/work/ba1442/robinson/kryos-gate/reports/cmp-S10-%j.log
# sbatch --export=ALL,KSHA=<sha> compare_S10.sh g1|g4|bc; reports and summary.txt in $G/reports/S10.
#   g1: K(S10, vilma=1) against K(S9) 624f686 (grl, ant, ant-full, nor) and sci, sci-const against K3-sci, kryos-compare --restarts.
#       The bundles of S10 lack the diagnostics (regional geo bsl, shelf tf_shlf): absent, by design.
#   g4: sci-vilma2 cont2 == cont, chunk == cont, rounding (cmp_g4.sh) and the sci chunk == cont (kryos-gate-kk.jl, that check only).
#   bc: tools/gate/sci-vilma2-bc.toml over $G/runs/sci-vilma2-bc/kx-$KSHA/S10 (gate_S10.sh bc).
set -uo pipefail
module load python3 >/dev/null 2>&1 || true
export PATH=$HOME/.juliaup/bin:$HOME/.local/bin:$PATH
G=/work/ba1442/robinson/kryos-gate; T=$G/p4w7/kryos/tools; KSHA=${KSHA:?}; KS9=624f686
export JULIA_DEPOT_PATH=$G/julia-depot
OUT=$G/reports/S10; mkdir -p $OUT
KK="julia +1.13 --project=$T $T/bin/kryos-gate-kk.jl check"
kc () {   # kc <name> <ref> <run>
    if [ ! -d "$2" ] || [ ! -d "$3" ]; then echo "$1 | MISSING | $2 $3" >> $OUT/summary.txt; return; fi
    julia +1.13 --project=$T $T/bin/kryos-compare.jl "$2" "$3" --md $OUT/$1.md --restarts > $OUT/$1.log 2>&1
    echo "$1 | $(sed -n 1p $OUT/$1.md | sed 's/# kryos-compare: //') | $(grep -m1 '^- compared:' $OUT/$1.md | sed 's/- compared: //')" >> $OUT/summary.txt
}
case ${1:-} in
  g1)
    for g in grl ant ant-full nor; do
        kc $g-K-S10-vs-S9     $G/runs/$g/kx-$KS9/S9/cont  $G/runs/$g/kx-$KSHA/S10/cont
    done
    for g in sci sci-const; do
        kc $g-K-S10-vs-K3sci  $G/ref/K3-sci/kx-2d80b43/S5b/$g/cont $G/runs/$g/kx-$KSHA/S10/cont
    done ;;
  g4)
    SHA=$KSHA STK=S10 bash $G/p4w7/levante/cmp_g4.sh > $OUT/g4-vilma2.log 2>&1
    echo "sci-vilma2 | $(sed -n 1p $OUT/sci-vilma2-kk.md) | g4-vilma2.log" >> $OUT/summary.txt
    $KK sci $G/runs/sci/kx-$KSHA/S10 --only "chunk vs cont (chunks)" --md $OUT/sci-chunk.md > $OUT/sci-chunk.log 2>&1
    echo "sci chunk | exit $? | $(sed -n 1p $OUT/sci-chunk.md)" >> $OUT/summary.txt ;;
  bc)
    $KK $G/rst10/sci-vilma2-bc.toml $G/runs/sci-vilma2-bc/kx-$KSHA/S10 --md $OUT/sci-vilma2-bc.md > $OUT/sci-vilma2-bc.log 2>&1
    echo "sci-vilma2-bc | exit $? | $(sed -n 1p $OUT/sci-vilma2-bc.md)" >> $OUT/summary.txt ;;
  *) echo "usage: KSHA=<sha> compare_S10.sh g1|g4|bc"; exit 1 ;;
esac
echo "== $(date -Is) done $1" >> $OUT/summary.txt
