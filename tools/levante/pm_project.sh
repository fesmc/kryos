#!/bin/bash
# pm_project.sh [<folder>] [<data folder>]: a folder of standalone VILMA2 runs as a runme project (default PM_DIR; safe to
# rerun): .runme/ (info.json for vilma.x of the VILMA2 checkout, the job templates and config.toml of the kryos checkout) and
# the links runme puts into every run: input (VILMA2's input/: vilma_defaults.nml), fesmdata, pm_data (the inputs of
# pm_prep.jl: <data folder>, default $PM_DIR/data if it exists, else the reference inputs $PM_REF/data). Then, in the folder,
#   runme -r -s -e vilma -n $KRYOS/tools/levante/paleomist/pm_1d.nml -o B1_l170 -q compute -w 08:00:00 --omp 64 \
#         -p vilma.rho_ice=931.0 vilma.equil_time_max=2.0e6
# (docs/experiment-gowan.md: the runs of Gowan et al. (2021) as cases, tools/levante/paleomist/pm_<earth>.nml).
set -euo pipefail
_h=$(dirname "$(readlink -f "$0")"); . "$_h/env.sh"
P=${1:-$PM_DIR}
D=${2:-}; if [ -z "$D" ]; then D=$PM_DIR/data; [ -d "$D" ] || D=$PM_REF/data; fi
[ -x "$VILMA2/bin/vilma.x" ] || { echo "no $VILMA2/bin/vilma.x (VILMA2: the vilma2 checkout; make vilma in it)"; exit 1; }
[ -d "$D" ] || { echo "no data folder $D (pm_prep.jl)"; exit 1; }
mkdir -p "$P/.runme"; P=$(cd "$P" && pwd)
cat > "$P/.runme/info.json" <<EOF
{
    "exe_default" : "vilma", "par_path_as_argument" : true,
    "exe_aliases" : {"vilma" : "$VILMA2/bin/vilma.x"}, "grp_aliases" : {},
    "par_paths" : {}, "files" : [], "dir-special" : {}, "par_defaults" : ["input/vilma_defaults.nml"], "par_library" : [],
    "links" : ["input", "fesmdata", "pm_data"]
}
EOF
for f in submit_slurm submit_slurm_omp; do cp "$KRYOS/.runme/$f" "$P/.runme/"; done
c=$KRYOS/.runme/config.toml; [ -f "$c" ] || c=$KRYOS/.runme/config.default.toml
cp "$c" "$P/.runme/config.toml"
ln -sfn "$VILMA2/input" "$P/input"; ln -sfn "$FESMDATA" "$P/fesmdata"; ln -sfn "$(cd "$D" && pwd)" "$P/pm_data"
echo "pm_project: $P (vilma.x of $VILMA2, pm_data -> $(readlink "$P/pm_data"))"
