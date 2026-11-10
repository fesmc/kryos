#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=ba1442
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --job-name=cmp-S5b
#SBATCH --output=/work/ba1442/robinson/kryos-gate/reports/cmp-S5b-%j.log
# sbatch compare_S5b.sh: K(S5b) cont against K(S5) cont, R2(S5b) cont against R2(S5) cont, grl and ant (outputs and
# restart bundles, kryos-compare --restarts); reports and summary.txt in $G/reports/S5b.
set -uo pipefail
module load python3 >/dev/null 2>&1 || true
export PATH=$HOME/.juliaup/bin:$HOME/.local/bin:$PATH
G=/work/ba1442/robinson/kryos-gate; T=$G/kryos/tools; KSHA=5a663a3   # kryos-compare of the shared checkout (its Manifest finds ../../FesmUtils.jl)
export JULIA_DEPOT_PATH=$G/julia-depot
OUT=$G/reports/S5b; mkdir -p $OUT
kc () {   # kc <name> <ref> <run>
    if [ ! -d "$2" ] || [ ! -d "$3" ]; then echo "$1 | MISSING | $2 $3" >> $OUT/summary.txt; return; fi
    julia +1.13 --project=$T $T/bin/kryos-compare.jl "$2" "$3" --md $OUT/$1.md --restarts > $OUT/$1.log 2>&1
    echo "$1 | $(sed -n 1p $OUT/$1.md | sed 's/# kryos-compare: //') | $(grep -m1 '^- compared:' $OUT/$1.md | sed 's/- compared: //')" >> $OUT/summary.txt
}
for g in grl ant; do
    kc $g-K-S5b-vs-S5  $G/runs/$g/kx-$KSHA/S5/cont $G/runs/$g/kx-$KSHA/S5b/cont
    kc $g-R2-S5b-vs-S5 $G/ref-S5/$g/yx-5ac89c8/S5/cont $G/ref-S5b/$g/yx-5ac89c8/S5b/cont
done
echo "== $(date -Is) done" >> $OUT/summary.txt
