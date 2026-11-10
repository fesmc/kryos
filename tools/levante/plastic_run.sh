#!/bin/bash
# plastic_run.sh <run folder> <par file> [<domain factor file> [<scale file>]] [group.key=value ...]
#
# One kryos run of a par file on Levante: the run folder with a copy of the driver, the data links (env.sh: fesmdata,
# kryos_data, input, maps) and the par file, its taub_domain_factor_file and taub_scale_file keys set to the files given
# on every instance ("" = none), further keys set with set_nml.py; then its Slurm job, submitted. The run writes
# `status` ("<exit code> <seconds>") when it ends; out.out, out.err its log.
#
#   tools/levante/plastic_run.sh $KRYOS_WORK/runs/sci-3d par/kryos_NH-ANT-PAT_plastic_vilma2.nml \
#       '{kryos_data}/plastic/taub_domain_factor_v3.nc'
#   tools/levante/plastic_run.sh $KRYOS_WORK/runs/b1-short par/kryos_NH-ANT_plastic_vilma2_b1.nml \
#       '{kryos_data}/plastic/taub_domain_factor_v2.nc' "" 'clock.time_end=-60e3'
#
# {kryos_data}, {fesmdata}, {kryos} in the file arguments are replaced by the paths of env.sh. A relative par file is
# taken from the checkout. Job settings (environment, defaults): PARTITION (shared), THREADS (64), MEM (60G),
# TIME (08:00:00), SUBMIT (1; 0: stage only). A 3-D Earth run of NH-ANT(-PAT) takes 3.5-4.2 h on 64 threads, a 1-D
# (b1) run about 15 min on 16 threads. The run folder must not exist. Build first: make kryos plastic=1 vilma=1.
set -euo pipefail
. "$(dirname "$(readlink -f "$0")")/env.sh"
[ $# -ge 2 ] || { sed -n '2,3p' "$0"; exit 1; }
PARTITION=${PARTITION:-shared}; THREADS=${THREADS:-64}; MEM=${MEM:-60G}; TIME=${TIME:-08:00:00}; SUBMIT=${SUBMIT:-1}
sub() { local s=${1//\{kryos_data\}/$KRYOS_DATA}; s=${s//\{fesmdata\}/$FESMDATA}; echo "${s//\{kryos\}/$KRYOS}"; }
R=$1; P=$(sub "$2"); FAC=$(sub "${3:-}"); SC=$(sub "${4:-}"); shift $(( $# < 4 ? $# : 4 ))
case $P in /*) ;; *) P=$KRYOS/$P ;; esac
[ -f "$P" ] || { echo "no par file $P"; exit 1; }
for f in "$FAC" "$SC"; do [ -z "$f" ] || [ -f "$f" ] || { echo "no file $f"; exit 1; }; done
[ -x "$KRYOS/libkryos/bin/kryos.x" ] || { echo "no $KRYOS/libkryos/bin/kryos.x: make kryos plastic=1 vilma=1"; exit 1; }
[ -e "$R" ] && { echo "exists: $R"; exit 1; }
mkdir -p "$R"; R=$(cd "$R" && pwd)
cp "$KRYOS/libkryos/bin/kryos.x" "$R/"
kryos_link_data "$R"
sed -e "s|^\( *taub_domain_factor_file *= *\)\"[^\"]*\"|\1\"$FAC\"|" -e "s|^\( *taub_scale_file *= *\)\"[^\"]*\"|\1\"$SC\"|" \
    "$P" > "$R/par.nml"
[ $# -gt 0 ] && python3 "$KRYOS/tools/levante/set_nml.py" "$R/par.nml" "$@"
echo "par file $P; $(grep -cE "^ *taub_domain_factor_file *= *\"$FAC\"" "$R/par.nml") instances with the factor file"
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
if [ "$SUBMIT" = 1 ]; then (cd "$R" && sbatch job.submit); fi
