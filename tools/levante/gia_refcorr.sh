#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=ba1442
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=128
#SBATCH --mem=60G
#SBATCH --job-name=gia_refcorr
# The present-day reference correction r of a glacial-cycle run (gia_refcorr.jl; pm-ref in docs/dev/vilma2-paleomist.md): standalone
# VILMA2 (vilma.x of a stack) under an ice history with the Earth of a kryos par file, iterated with i_eq = 3 until its 0 ka bed is
# today's (r_{n+1} = rsl_n(0 ka), r_0 = 0), then the kryos reference files (today's topography + r, H_ice present-day).
#   sbatch --export=ALL,CASE=PaleoMIST-sci,PAR=<kryos par file> gia_refcorr.sh
# Env:
#   CASE     name of the case: everything goes to $KD/refcorr/$CASE (KD: /work/ba1442/kryos_data, the kryos_data link of a run)
#   PAR      kryos par file: its &vilma group (overrides of the stack's input/vilma_defaults.nml, as the geo = vilma2 backend reads it)
#   STACK    the stack of vilma.x (S11)
#   HIST     the ice history (FesmData IceHistory product, lon-lat, years ago; default PaleoMIST on GLOBAL-0.25DEG)
#   INTERP   linear (default: the records linear in time, as a datasrc history) | steps (SELEN's Heaviside steps)
#   RHOI, RHOW  densities of the load (kryos sets VILMA2's from phys_const: rho_ice 910, rho_sw 1028)
#   TOL, NMAX   stop when max |r_n - rsl_n(0 ka)| < TOL [m] (1.5), or after NMAX runs (8)
#   XGRID    the globe_ice background grid (GLOBAL-0.5DEG) and its history file (XHIST, default HIST on that grid)
#   DOMAINS  "<domain>:<hub grid>:<topo file under fesmdata> ..." the domain reference topographies to correct
#   NTHR     OpenMP threads (64)
# Restartable: an iteration with a finished vilma_out.nc is not rerun.
set -uo pipefail
G=/work/ba1442/robinson/kryos-gate; F=/work/ba1442/fesmdata
T=$(dirname $(readlink -f $0)); [ -f $T/gia_refcorr.jl ] || T=${TOOLS:?TOOLS: the folder of gia_refcorr.jl (sbatch copies the script)}
export PATH=$HOME/.juliaup/bin:$PATH JULIA_DEPOT_PATH=$G/julia-depot
JL="julia +1.13 --project=$G/p4w7/kryos/tools"
CASE=${CASE:?}; PAR=${PAR:?}; STACK=${STACK:-S11}
HIST=${HIST:-$F/Global/GLOBAL-0.25DEG/GLOBAL-0.25DEG_IceHistory-PaleoMIST.nc}; INTERP=${INTERP:-linear}
RHOI=${RHOI:-910.0}; RHOW=${RHOW:-1028.0}; TOL=${TOL:-1.5}; NMAX=${NMAX:-8}; NTHR=${NTHR:-64}
XGRID=${XGRID:-GLOBAL-0.5DEG}; XHIST=${XHIST:-$F/Global/$XGRID/${XGRID}_IceHistory-$(basename $HIST | sed 's/.*_IceHistory-//')}
DOMAINS=${DOMAINS:-"North:NH-32KM:North/NH-32KM/NH-32KM_TOPO-GEBCO2026.nc Antarctica:ANT-16KM:Antarctica/ANT-16KM/ANT-16KM_TOPO-BedMachine-v4.nc"}
VX=$G/stacks/$STACK/vilma2/bin/vilma.x
KD=${KD:-/work/ba1442/kryos_data}; O=$KD/refcorr/$CASE; mkdir -p $O/data
cp $PAR $O/par.nml
# the &vilma group of the par file (the geo = vilma2 backend reads it over the stack's defaults)
awk '/^&vilma[ \t]*$/{f=1} f{print} f&&/^\//{exit}' $PAR | sed '1d;$d' > $O/vilma.nml
[ -s $O/vilma.nml ] || { echo "no &vilma group in $PAR"; exit 1; }
echo "== gia_refcorr $CASE: $PAR, $HIST ($INTERP), rho $RHOI/$RHOW, stack $STACK, tol $TOL m $(date -Is)" | tee -a $O/iter.log
[ -f $O/data/force.nc ] || $JL $T/gia_refcorr.jl prep $HIST $INTERP $O/data || exit 1

export OMP_NUM_THREADS=$NTHR OMP_PROC_BIND=spread OMP_PLACES=cores OMP_STACKSIZE=512M
ulimit -s unlimited
prev=""; conv=0
for n in $(seq 0 $((NMAX - 1))); do
  R=$O/it$n; r=$O/r_$((n + 1)).nc
  if ! grep -q "vilma.x exit 0" $R/run.out 2>/dev/null; then
    rm -rf $R; mkdir -p $R/input; cd $R
    for f in $G/stacks/$STACK/vilma2/input/*; do ln -s $f input/; done
    ln -s $F fesmdata
    cat > vilma_rc.nml <<EOF
&ctl
 file_forcing = "$O/data/force.nc"
 name_ice = "h_ice"
 name_time = "time"
 file_ref = ""
 name_zbed_eq = "z_bed_eq"
 name_hice_ref = "h_ice_ref"
 file_out = "vilma_out.nc"
 file_hor = ""
 write_lossless = .true.
 time_init = -1.0e30
 time_end = 1.0e30
 remap_input = .true.
 name_lon = "lon"
 name_lat = "lat"
 i_eq = $([ -n "$prev" ] && echo 3 || echo 1)
 z_bed_ref_file = "$O/data/ref.nc"
 h_ice_ref_file = "$O/data/ref.nc"
 name_z_bed_ref = "z_bed"
 name_h_ice_ref = "H_ice"
 z_bed_eq_file = ""
 h_ice_eq_file = ""
 rsl_restart_file = "$prev"
 name_rsl = "rsl"
 restart_in_file = ""
/
&vilma
$(cat $O/vilma.nml)
 rho_ice = $RHOI
 rho_water = $RHOW
/
EOF
    echo "== it$n $(date -Is) on $(hostname)" | tee -a $O/iter.log
    { /usr/bin/time -v $VX vilma_rc.nml > run.log 2>&1; echo "vilma.x exit $?"; } > run.out 2>&1
    grep -E "exit|Elapsed" run.out | tee -a $O/iter.log
    grep -q "vilma.x exit 0" run.out || { echo "it$n failed: $R/run.log" | tee -a $O/iter.log; exit 1; }
  fi
  $JL $T/gia_refcorr.jl next $R/vilma_out.nc $r $prev | tee -a $O/iter.log
  m=$(tail -1 $O/iter.log | sed -n 's/.*max |d| = \([0-9.]*\).*/\1/p')
  if [ -n "$prev" ] && awk -v m=$m -v t=$TOL 'BEGIN{exit !(m < t)}'; then conv=1; final=$prev; echo "converged: it$n (r_$n), max misfit $m m" | tee -a $O/iter.log; break; fi
  prev=$r
done
[ $conv = 1 ] || { echo "not converged after $NMAX runs" | tee -a $O/iter.log; exit 1; }
ln -sf $(basename $final) $O/r.nc
DOMS=""; for d in $DOMAINS; do IFS=: read dom hub topo <<< "$d"; DOMS="$DOMS $dom:$hub:$F/$topo"; done
RHO_ICE=$RHOI RHO_SW=$RHOW $JL $T/gia_refcorr.jl write $O/r.nc $O $XGRID:$XHIST $DOMS | tee -a $O/iter.log
echo "== done $(date -Is)" | tee -a $O/iter.log
