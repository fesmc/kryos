#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=ba1442
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --job-name=cmp-g6
#SBATCH --output=/work/ba1442/robinson/kryos-gate/reports/WP7c/cmp-g6-%j.log
# G6: chunk (restart at the middle) == cont, bundles included, for the gate g6 (2 kyr) or g6-smb (smb_simple, 1 kyr): GATE=g6|g6-smb sbatch cmp_g6.sh
set -uo pipefail
export PATH=$HOME/.juliaup/bin:$HOME/.local/bin:$PATH JULIA_DEPOT_PATH=/work/ba1442/robinson/kryos-gate/julia-depot
G=/work/ba1442/robinson/kryos-gate; T=$G/p4w7/kryos/tools; GATE=${GATE:-g6}
julia +1.13 --project=$T $T/bin/kryos-gate-kk.jl check $G/p4w7c/gate/$GATE.toml $G/runs/$GATE/kx-e4453dc/S8 --md $G/reports/WP7c/$GATE-kk.md
echo "exit $?"; sed -n 1,40p $G/reports/WP7c/$GATE-kk.md
