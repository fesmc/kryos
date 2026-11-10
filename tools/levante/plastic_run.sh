#!/bin/bash
# plastic_run.sh <run folder> <par file> [<domain factor file> [<scale file>]] [group.key=value ...]
#
# One kryos run of a plastic par file on Levante: the run folder with a copy of the driver, the data links and the par
# file (its taub_domain_factor_file and taub_scale_file keys set to the files given, "" = none; further keys set with
# set_nml.py), and a Slurm job, submitted. The run writes status ("<exit code> <seconds>") when it ends.
#
#   tools/levante/plastic_run.sh runs/sci-3d par/kryos_NH-ANT-PAT_plastic_vilma2.nml \
#       /work/ba1442/robinson/kryos-ice-plastic/runs/taub-fit-v3/report/taub_domain_factor.nc
#   tools/levante/plastic_run.sh runs/b1-short par/kryos_North_plastic_vilma2_b1.nml "" "" 'clock.time_end=-60e3'
#
# Environment (defaults): K (the kryos checkout: this script's), FESMDATA (/work/ba1442/fesmdata, read-only),
# KRYOS_DATA (/work/ba1442/kryos_data), MAPS (/work/ba1442/robinson/kryos-ice-plastic-maps: the grid descriptions and
# the remapping weights; the first run of a new grid writes its weights there), PARTITION (shared), ACCOUNT (ba1442),
# THREADS (64), MEM (60G), TIME (08:00:00), SUBMIT (1; 0: stage only). A 3-D Earth run of NH-ANT(-PAT) takes 3.5-4.2 h
# on 64 threads; the 1-D (b1) runs about 1 h.
set -euo pipefail
[ $# -ge 2 ] || { sed -n '2,3p' "$0"; exit 1; }
HERE=$(cd "$(dirname "$0")/../.." && pwd)
K=${K:-$HERE}; FESMDATA=${FESMDATA:-/work/ba1442/fesmdata}; KRYOS_DATA=${KRYOS_DATA:-/work/ba1442/kryos_data}
MAPS=${MAPS:-/work/ba1442/robinson/kryos-ice-plastic-maps}
PARTITION=${PARTITION:-shared}; ACCOUNT=${ACCOUNT:-ba1442}; THREADS=${THREADS:-64}; MEM=${MEM:-60G}; TIME=${TIME:-08:00:00}
SUBMIT=${SUBMIT:-1}
R=$1; P=$2; FAC=${3:-}; SC=${4:-}; shift $(( $# < 4 ? $# : 4 ))
case $P in /*) ;; *) P=$K/$P ;; esac
[ -x "$K/libkryos/bin/kryos.x" ] || { echo "no $K/libkryos/bin/kryos.x: make kryos plastic=1 vilma=1"; exit 1; }
mkdir -p "$R"; R=$(cd "$R" && pwd)
cp "$K/libkryos/bin/kryos.x" "$R/"
ln -sfn "$FESMDATA" "$R/fesmdata"; ln -sfn "$KRYOS_DATA" "$R/kryos_data"; ln -sfn "$K/input" "$R/input"; ln -sfn "$MAPS" "$R/maps"
sed -e "s|^\( *taub_domain_factor_file *= *\)\"[^\"]*\"|\1\"$FAC\"|" -e "s|^\( *taub_scale_file *= *\)\"[^\"]*\"|\1\"$SC\"|" \
    "$P" > "$R/par.nml"
[ $# -gt 0 ] && python3 "$K/tools/levante/set_nml.py" "$R/par.nml" "$@"
grep -E "^ *taub_(domain_factor|scale)_file" "$R/par.nml" || true
cat > "$R/job.submit" <<EOS
#!/bin/bash
#SBATCH --partition=$PARTITION
#SBATCH --account=$ACCOUNT
#SBATCH --time=$TIME
#SBATCH --job-name=$(basename "$R")
#SBATCH --output=./out.out
#SBATCH --error=./out.err
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$MEM
export OMP_PROC_BIND=true OMP_NUM_THREADS=$THREADS OMP_STACKSIZE=512M OMP_MAX_ACTIVE_LEVELS=1
ulimit -s unlimited
t0=\$(date +%s); ./kryos.x par.nml; echo "\$? \$(( \$(date +%s) - t0 ))" > status
EOS
echo "staged $R"
[ "$SUBMIT" = 1 ] && (cd "$R" && sbatch job.submit)
