#!/bin/bash
# The present-day reference correction r of a glacial-cycle run (gia_refcorr.jl; pm-ref in docs/dev/vilma2-paleomist.md): standalone
# VILMA2 (vilma.x of a stack) under an ice history with the Earth of a kryos par file, iterated with i_eq = 3 until its 0 ka bed is
# today's (r_{n+1} = rsl_n(0 ka), r_0 = 0), then the kryos reference files (today's topography + r, H_ice present-day).
#   CASE=PaleoMIST-sci PAR=par/kryos_NH-ANT_plastic_vilma2.nml tools/levante/gia_refcorr.sh     (submits itself: shared,
#   2 days, 128 cpus, 60 GB)
# Env:
#   CASE     name of the case: everything goes to $KD/refcorr/$CASE (KD default $KRYOS_WORK/kryos_data: a run reads the case
#            through its kryos_data link, so point KRYOS_DATA there, or copy the case into the shared KRYOS_DATA)
#   PAR      kryos par file: its &vilma group (overrides of VILMA2's input/vilma_defaults.nml, as the geo = vilma2 backend reads it)
#   VILMA2   the vilma2 checkout of vilma.x (env.sh; the published cases used the gate stack S11)
#   HIST     the ice history (FesmData IceHistory product, lon-lat, years ago; default PaleoMIST on GLOBAL-0.25DEG)
#   INTERP   linear (default: the records linear in time, as a datasrc history) | steps (SELEN's Heaviside steps)
#   RHOI, RHOW  densities of the load (kryos sets VILMA2's from phys_const: rho_ice 910, rho_sw 1028)
#   TOL, NMAX   stop when max |r_n - rsl_n(0 ka)| < TOL [m] (1.5), or after NMAX runs (8)
#   XGRID    the globe_ice background grid (GLOBAL-0.5DEG) and its history file (XHIST, default HIST on that grid)
#   DOMAINS  "<domain>:<hub grid>:<topo file under fesmdata> ..." the domain reference topographies to correct (default North and
#            Antarctica; the case PaleoMIST-sci also has PatagoniaPaleo:PAT-PAL-8KM:PatagoniaPaleo/PAT-PAL-8KM/PAT-PAL-8KM_TOPO-GEBCO2026.nc)
#   NTHR     OpenMP threads (64)
#   R0       a present-day correction r to run with, once and without iterating (the members of an Earth-only ensemble: the r of a
#            nearby Earth); the run is it0, no reference files are written
# Restartable: an iteration whose status is 0 is not rerun. The folder is a PaleoMIST runme project (pm_project.sh): the
# iterations it<n> are runs of its case vilma_rc.nml, staged by runme.
set -uo pipefail
_h=$(dirname "$(readlink -f "$0")"); [ -f "$_h/env.sh" ] || _h=$KRYOS/tools/levante; . "$_h/env.sh"
T=$_h; F=$FESMDATA
JL="$KRYOS_JULIA --project=$_h"
CASE=${CASE:?}; PAR=${PAR:?}
case $PAR in /*) ;; *) PAR=$KRYOS/$PAR ;; esac
KD=${KD:-$KRYOS_WORK/kryos_data}
if [ "${KRYOS_JOB:-}" != 1 ]; then
    [ -f "$PAR" ] || { echo "no par file $PAR"; exit 1; }
    mkdir -p $KD/refcorr/$CASE
    sbatch --parsable -A $ACCOUNT $KRYOS_PART_LONG -t 2-00:00:00 -c 128 --mem=60G -J gia_refcorr-$CASE \
        -o $KD/refcorr/$CASE/job-%j.out --export=ALL,KRYOS_JOB=1,PAR=$PAR,KD=$KD "$_h/gia_refcorr.sh"
    exit
fi
HIST=${HIST:-$F/Global/GLOBAL-0.25DEG/GLOBAL-0.25DEG_IceHistory-PaleoMIST.nc}; INTERP=${INTERP:-linear}
RHOI=${RHOI:-910.0}; RHOW=${RHOW:-1028.0}; TOL=${TOL:-1.5}; NMAX=${NMAX:-8}; NTHR=${NTHR:-64}
XGRID=${XGRID:-GLOBAL-0.5DEG}; XHIST=${XHIST:-$F/Global/$XGRID/${XGRID}_IceHistory-$(basename $HIST | sed 's/.*_IceHistory-//')}
DOMAINS=${DOMAINS:-"North:NH-32KM:North/NH-32KM/NH-32KM_TOPO-GEBCO2026.nc Antarctica:ANT-16KM:Antarctica/ANT-16KM/ANT-16KM_TOPO-BedMachine-v4.nc"}
VX=$VILMA2/bin/vilma.x
[ -x $VX ] || { echo "no $VX (VILMA2)"; exit 1; }
O=$KD/refcorr/$CASE; mkdir -p $O/data
cp $PAR $O/par.nml
# the &vilma group of the par file (the geo = vilma2 backend reads it over the stack's defaults)
awk '/^&vilma[ \t]*$/{f=1} f{print} f&&/^\//{exit}' $PAR | sed '1d;$d' > $O/vilma.nml
[ -s $O/vilma.nml ] || { echo "no &vilma group in $PAR"; exit 1; }
echo "== gia_refcorr $CASE: $PAR, $HIST ($INTERP), rho $RHOI/$RHOW, vilma.x $VX, tol $TOL m $(date -Is)" | tee -a $O/iter.log
[ -f $O/data/force.nc ] || $JL $T/gia_refcorr.jl prep $HIST $INTERP $O/data || exit 1

# the folder as a PaleoMIST runme project (pm_project.sh: pm_data -> data/), the iterations its runs of the case vilma_rc.nml
bash $T/pm_project.sh $O $O/data > /dev/null || exit 1
cat > $O/vilma_rc.nml <<EOF
&ctl
 file_forcing = "pm_data/force.nc"
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
 i_eq = 1
 z_bed_ref_file = "pm_data/ref.nc"
 h_ice_ref_file = "pm_data/ref.nc"
 name_z_bed_ref = "z_bed"
 name_h_ice_ref = "H_ice"
 z_bed_eq_file = ""
 h_ice_eq_file = ""
 rsl_restart_file = ""
 name_rsl = "rsl"
 restart_in_file = ""
/
&vilma
$(cat $O/vilma.nml)
 rho_ice = $RHOI
 rho_water = $RHOW
/
EOF
export OMP_NUM_THREADS=$NTHR OMP_PROC_BIND=true OMP_STACKSIZE=512M OMP_MAX_ACTIVE_LEVELS=1
ulimit -s unlimited
prev=${R0:-}; conv=0
for n in $(seq 0 $((NMAX - 1))); do
  R=$O/it$n; r=$O/r_$((n + 1)).nc
  if [ "$(cut -d' ' -f1 $R/status 2> /dev/null)" != 0 ]; then
    ref=""; [ -n "$prev" ] && ref="-p ctl.i_eq=3 ctl.rsl_restart_file=$prev"
    (cd $O && $KRYOS_RUNME -e vilma -n $O/vilma_rc.nml -o it$n $ref > /dev/null) || { echo "it$n: runme failed" | tee -a $O/iter.log; exit 1; }
    echo "== it$n $(date -Is) on $(hostname)" | tee -a $O/iter.log
    (cd $R && t0=$(date +%s) && ./vilma.x vilma_rc.nml > out.out 2> out.err; echo "$? $(( $(date +%s) - t0 ))" > status)
    read rc sec < $R/status; echo "it$n: exit $rc, $sec s" | tee -a $O/iter.log
    [ "$rc" = 0 ] || { echo "it$n failed: $R/out.err" | tee -a $O/iter.log; exit 1; }
  fi
  $JL $T/gia_refcorr.jl next $R/vilma_out.nc $r $prev | tee -a $O/iter.log
  m=$(tail -1 $O/iter.log | sed -n 's/.*max |d| = \([0-9.]*\).*/\1/p')
  if [ -n "$prev" ] && awk -v m=$m -v t=$TOL 'BEGIN{exit !(m < t)}'; then conv=1; final=$prev; echo "converged: it$n (r_$n), max misfit $m m" | tee -a $O/iter.log; break; fi
  [ -n "${R0:-}" ] && { echo "== done (one run with r = $R0) $(date -Is)" | tee -a $O/iter.log; exit 0; }
  prev=$r
done
[ $conv = 1 ] || { echo "not converged after $NMAX runs" | tee -a $O/iter.log; exit 1; }
ln -sf $(basename $final) $O/r.nc
DOMS=""; for d in $DOMAINS; do IFS=: read dom hub topo <<< "$d"; DOMS="$DOMS $dom:$hub:$F/$topo"; done
RHO_ICE=$RHOI RHO_SW=$RHOW $JL $T/gia_refcorr.jl write $O/r.nc $O $XGRID:$XHIST $DOMS | tee -a $O/iter.log
echo "== done $(date -Is)" | tee -a $O/iter.log
