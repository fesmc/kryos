#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=ba1442
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --job-name=cmp-S8
#SBATCH --output=/work/ba1442/robinson/kryos-gate/reports/cmp-S8-%j.log
# sbatch --export=ALL,KSHA=<sha> compare_S8.sh g1: kryos-compare --restarts (cont, outputs and bundles); reports and summary.txt in $G/reports/S8.
#   g1: K(S8, vilma=1) against the earlier K (grl, ant: S5b; ant-full: S5 5a663a3; nor: S5 888884d) and sci, sci-const against K3-sci
set -uo pipefail
module load python3 >/dev/null 2>&1 || true
export PATH=$HOME/.juliaup/bin:$HOME/.local/bin:$PATH
G=/work/ba1442/robinson/kryos-gate; T=$G/p4w7/kryos/tools; KSHA=${KSHA:?}
export JULIA_DEPOT_PATH=$G/julia-depot
OUT=$G/reports/S8; mkdir -p $OUT
kc () {   # kc <name> <ref> <run>
    if [ ! -d "$2" ] || [ ! -d "$3" ]; then echo "$1 | MISSING | $2 $3" >> $OUT/summary.txt; return; fi
    julia +1.13 --project=$T $T/bin/kryos-compare.jl "$2" "$3" --md $OUT/$1.md --restarts > $OUT/$1.log 2>&1
    echo "$1 | $(sed -n 1p $OUT/$1.md | sed 's/# kryos-compare: //') | $(grep -m1 '^- compared:' $OUT/$1.md | sed 's/- compared: //')" >> $OUT/summary.txt
}
case ${1:-} in
  g1)
    kc grl-K-S8-vs-S5b       $G/runs/grl/kx-5a663a3/S5b/cont      $G/runs/grl/kx-$KSHA/S8/cont
    kc ant-K-S8-vs-S5b       $G/runs/ant/kx-5a663a3/S5b/cont      $G/runs/ant/kx-$KSHA/S8/cont
    kc ant-full-K-S8-vs-S5   $G/runs/ant-full/kx-5a663a3/S5/cont  $G/runs/ant-full/kx-$KSHA/S8/cont
    kc nor-K-S8-vs-S5        $G/runs/nor/kx-888884d/S5/cont       $G/runs/nor/kx-$KSHA/S8/cont
    for g in sci sci-const; do
        kc $g-K-S8-vs-K3sci  $G/ref/K3-sci/kx-2d80b43/S5b/$g/cont $G/runs/$g/kx-$KSHA/S8/cont
    done ;;
  *) echo "usage: KSHA=<sha> compare_S8.sh g1"; exit 1 ;;
esac
echo "== $(date -Is) done $1" >> $OUT/summary.txt
