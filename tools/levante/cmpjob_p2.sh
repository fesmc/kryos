#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=ba1442
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --job-name=cmp-p2
#SBATCH --output=/work/ba1442/robinson/kryos-gate/reports/cmp-p2-%j.log
# sbatch cmpjob_p2.sh <gate> <kryos sha> [stack] [R2 root]
exec bash -l /work/ba1442/robinson/kryos-gate/tools/compare_p2.sh "$@"
