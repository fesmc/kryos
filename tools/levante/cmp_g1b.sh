#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=ba1442
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --job-name=cmp-g1b
#SBATCH --output=/work/ba1442/robinson/kryos-gate/reports/cmp-g1b-%j.log
# G1b: (1) sci-hub cont against K3-sci sci cont (kryos-compare --restarts; the hub adds globe_ice.nc: a new file), (3) cont2 and the
# domain swap identical (kryos-gate-kk.jl check sci-hub), (4) globe_ice.nc against the offline merge (hub_offline.jl). Reports in
# $G/reports/S6; the hub closure of (2) is in the logs (grep "hub:").
set -uo pipefail
module load python3 >/dev/null 2>&1 || true
export PATH=$HOME/.juliaup/bin:$HOME/.local/bin:$PATH JULIA_DEPOT_PATH=/work/ba1442/robinson/kryos-gate/julia-depot
G=/work/ba1442/robinson/kryos-gate; T=$G/p4w7/kryos/tools; SHA=cd85c14; OUT=$G/reports/S6; mkdir -p $OUT
R=$G/runs/sci-hub/kx-$SHA/S6
echo "== (1) vs K3-sci"
julia +1.13 --project=$T $T/bin/kryos-compare.jl $G/ref/K3-sci/kx-2d80b43/S5b/sci/cont $R/cont --md $OUT/sci-hub-vs-K3sci.md --restarts > $OUT/sci-hub-vs-K3sci.log 2>&1
echo "exit $?"; sed -n 1,12p $OUT/sci-hub-vs-K3sci.md
echo "== (3) cont2, swap"
julia +1.13 --project=$T $T/bin/kryos-gate-kk.jl check sci-hub $R --md $OUT/sci-hub-kk.md > $OUT/sci-hub-kk.log 2>&1
echo "exit $?"; sed -n 1,25p $OUT/sci-hub-kk.md
echo "== (4) offline merge"
julia +1.13 --project=$T $G/p4w7/levante/hub_offline.jl $R/cont $G/p4w7/kryos-mg/libkryos/bin/hub_areas.x /work/ba1442/fesmdata/Global/GLOBAL-0.5DEG/GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc --tol-abs 3e-3 2>&1 | tail -3   # 3e-3 m: the ice inputs are the float32 of the Yelmo output (ulp 5e-4 m at 4.5 km)
echo "== (2) hub closure of the logs"
grep -h "hub:" $R/cont/out.out | sed 's/^INFO *//' | sort | uniq -c | sort -rn | head -8
