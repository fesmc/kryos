#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=ba1442
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --job-name=cmp-S9
#SBATCH --output=/work/ba1442/robinson/kryos-gate/reports/cmp-S9-%j.log
# sbatch --export=ALL,KSHA=<sha> compare_S9.sh g1: kryos-compare --restarts (cont, outputs and bundles); reports and summary.txt in $G/reports/S9.
#   g1: K(S9, vilma=1) against the earlier K (all four: S8 e4453dc, bit identical to S5b/S5) and sci, sci-const against K3-sci.
#   The restart bundles gain the Yelmo restart variables pc_taus, pc_tau, pc_tau_masked (new, N) in S9.
set -uo pipefail
module load python3 >/dev/null 2>&1 || true
export PATH=$HOME/.juliaup/bin:$HOME/.local/bin:$PATH
G=/work/ba1442/robinson/kryos-gate; T=$G/p4w7/kryos/tools; KSHA=${KSHA:?}; KS8=e4453dc
export JULIA_DEPOT_PATH=$G/julia-depot
OUT=$G/reports/S9; mkdir -p $OUT
kc () {   # kc <name> <ref> <run>
    if [ ! -d "$2" ] || [ ! -d "$3" ]; then echo "$1 | MISSING | $2 $3" >> $OUT/summary.txt; return; fi
    julia +1.13 --project=$T $T/bin/kryos-compare.jl "$2" "$3" --md $OUT/$1.md --restarts > $OUT/$1.log 2>&1
    echo "$1 | $(sed -n 1p $OUT/$1.md | sed 's/# kryos-compare: //') | $(grep -m1 '^- compared:' $OUT/$1.md | sed 's/- compared: //')" >> $OUT/summary.txt
}
case ${1:-} in
  g1)
    for g in grl ant ant-full nor; do
        kc $g-K-S9-vs-S8     $G/runs/$g/kx-$KS8/S8/cont  $G/runs/$g/kx-$KSHA/S9/cont
    done
    for g in sci sci-const; do
        kc $g-K-S9-vs-K3sci  $G/ref/K3-sci/kx-2d80b43/S5b/$g/cont $G/runs/$g/kx-$KSHA/S9/cont
    done ;;
  *) echo "usage: KSHA=<sha> compare_S9.sh g1"; exit 1 ;;
esac
echo "== $(date -Is) done $1" >> $OUT/summary.txt
