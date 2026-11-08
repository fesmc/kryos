#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=ba1442
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --job-name=jl-score
#SBATCH --output=score-%j.out
# score.sh <kryos run dir> <label>: the run's VILMA2 stream against GAPSLIP 2.1 with pm_compare.jl, beside B1 (B1_l170, rho_ice 931)
# and I910_it0 (B1 at 910 kg m-3), as in docs/dev/vilma2-paleomist.md. Output: runs/score/<label>/out.
set -e
G=/work/ba1442/robinson/kryos-gate; P=$G/runs/paleomist; K=/work/ba1442/robinson/kryos-ice-plastic
RUN=$1; LBL=$2; S=$K/runs/score/$LBL
export PATH=$HOME/.juliaup/bin:$PATH JULIA_DEPOT_PATH=$G/julia-depot
J="julia +1.13 --project=$G/p4w7/kryos/tools"
$J $K/tools/levante/pm_kryos_out.jl $RUN $S/$LBL
$J $K/tools/levante/pm_compare.jl /work/ba1442/fesmdata $S $S/out B1=$P/B1_l170 B1_170=$P/ref/I910_it0 B3=$S/$LBL
