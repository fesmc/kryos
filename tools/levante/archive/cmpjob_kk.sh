#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=ba1442
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --job-name=cmp-kk
#SBATCH --output=/work/ba1442/robinson/kryos-gate/reports/cmp-kk-%j.log
# sbatch cmpjob_kk.sh <gate> <kryos sha> [stack]
exec bash -l /work/ba1442/robinson/kryos-gate/tools/compare_kk.sh "$@"
