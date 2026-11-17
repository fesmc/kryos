#!/bin/bash
# The present-day reference correction r of a glacial-cycle run (gia_refcorr.jl; pm-ref in docs/dev/vilma2-paleomist.md): standalone
# VILMA2 (vilma.x of a stack) under an ice history with the Earth of a kryos par file, iterated with i_eq = 3 until its 0 ka bed is
# today's (r_{n+1} = rsl_n(0 ka), r_0 = 0), then the kryos reference files (today's topography + r, H_ice present-day).
#   CASE=PaleoMIST-sci PAR=par/kryos_NH-ANT_plastic_vilma2.nml tools/levante/gia_refcorr.sh     (submits itself: shared,
#   2 days, 128 cpus, 60 GB; CHAIN=3: and up to 3 successors, for passes longer than a job)
# Env:
#   CASE     name of the case: everything goes to $KD/refcorr/$CASE (KD default $KRYOS_WORK/kryos_data: a run reads the case
#            through its kryos_data link, so point KRYOS_DATA there, or copy the case into the shared KRYOS_DATA)
#   PAR      kryos par file: its &vilma group (overrides of VILMA2's input/vilma_defaults.nml, as the geo = vilma2 backend reads it)
#   VILMA2   the vilma2 checkout of vilma.x (env.sh; the published cases used the gate stack S11)
#   HIST     the ice history (FesmData IceHistory product, lon-lat, years ago; default PaleoMIST on GLOBAL-0.25DEG)
#   INTERP   linear (default: the records linear in time, as a datasrc history) | steps (SELEN's Heaviside steps)
#   ZREF     an IceHistory product on the grid of HIST whose 0 ka z_bed is today's bed, for a history without z_bed (GLAC1D_80k:
#            the PaleoMIST product)
#   RHOI, RHOW  densities of the load (kryos sets VILMA2's from phys_const: rho_ice 910, rho_sw 1028)
#   TOL, NMAX   stop when max |r_n - rsl_n(0 ka)| < TOL [m] (1.5), or after NMAX runs (8)
#   XGRID    the globe_ice background grid (GLOBAL-0.5DEG) and its history file (XHIST, default HIST on that grid)
#   XZREF    ZREF on XGRID, the 0 ka bed of the reference files for a history without z_bed (default ZREF on that grid; none
#            without ZREF)
#   DOMAINS  "<domain>:<hub grid>:<topo file under fesmdata> ..." the domain reference topographies to correct (default North on NH-16KM and
#            NH-32KM, and Antarctica; the case PaleoMIST-sci also has PatagoniaPaleo:PAT-PAL-8KM:PatagoniaPaleo/PAT-PAL-8KM/PAT-PAL-8KM_TOPO-GEBCO2026.nc)
#   NCPU, MEM, TIME  cpus, memory and time of the job (128, 60G, 2-00:00:00): a 3-D Earth; a 1-D Earth needs far fewer
#            (16-32 cpus, about 3 GB at degree 170; many 1-D cases: gia_refcorr_pool.sh)
#   NTHR     OpenMP threads (NCPU, at most 64)
#   KEEP     all (default) | last: once converged, delete the model output of the iterations (vilma_out.nc, spinup/, final/)
#            but the converged iteration's vilma_out.nc, what the scores read (a 1-D case: about 0.3 instead of 3 GB)
#   R0       a present-day correction r to run with, once and without iterating (the members of an Earth-only ensemble: the r of a
#            nearby Earth); the run is it0, no reference files are written
#   TMIN     [min] start a pass or a chunk of one only with at least this much of the job left (120; 240 for a 3-D Earth),
#            from SLURM_JOB_END_TIME (no job end: always)
#   CHAIN    m: at its start the job submits its successor (afterany, CHAIN = m - 1), which exits at once when the case is
#            done; a pass longer than a job goes on in the next
#   CHUNK_TEND  (tests) time_end [yr] of the chunks, e.g. "-60000 -30000": a pass in three chunks
# Provenance: <case>/provenance.txt, appended at every start (env, git commits, md5 of the par file and vilma.x, the inputs)
# and at every chunk (threads, node, vilma2 commit, md5 of vilma.x); a chunk that resumes one on another thread count or node
# type is warned of in iter.log (it may not be bit for bit).
# Restartable: an iteration whose status is 0 is not rerun. The folder is a PaleoMIST runme project (pm_project.sh): the
# iterations it<n> are runs of its case vilma_rc.nml, staged by runme. A pass runs in chunks it<n>/c<k> (vilma_out.nc,
# final/, out.out, status "exit seconds t_last"): vilma.x stops before the job's end (SLURM_JOB_END_TIME, exit 75) with
# final/vilma_restart.nc, and the next chunk resumes it, in this job if TMIN is left, else in the next ("it<n>: continues
# at t = ..." in iter.log; the script exits 0). A chunk without status died: it is deleted and redone from the chunk
# before (c0 from its spin-up restart if it got that far, kept as c0-spinup). When the last chunk has run to the end, the
# outputs are joined into it<n>/vilma_out.nc (gia_refcorr.jl join: bit for bit at the seams), the logs into it<n>/out.out,
# c0's spinup/ and the last chunk's final/ are moved to it<n>/, it<n>/status is written, and the chunks are deleted.
set -uo pipefail
_h=$(dirname "$(readlink -f "$0")"); [ -f "$_h/env.sh" ] || _h=$KRYOS/tools/levante; . "$_h/env.sh"
T=$_h; F=$FESMDATA
JL="$KRYOS_JULIA --project=$_h"
CASE=${CASE:?}; PAR=${PAR:?}
case $PAR in /*) ;; *) PAR=$KRYOS/$PAR ;; esac
KD=${KD:-$KRYOS_WORK/kryos_data}
submit() {   # the job of the case ($1: sbatch options before the script)
    sbatch --parsable -A $ACCOUNT ${PART:-$KRYOS_PART_LONG} -t ${TIME:-2-00:00:00} -c ${NCPU:-128} --mem=${MEM:-60G} \
        -J gia_refcorr-$CASE -o $KD/refcorr/$CASE/job-%j.out --export=ALL,KRYOS_JOB=1,PAR=$PAR,KD=$KD "$@" "$_h/gia_refcorr.sh"
}
if [ "${KRYOS_JOB:-}" != 1 ]; then
    [ -f "$PAR" ] || { echo "no par file $PAR"; exit 1; }
    mkdir -p $KD/refcorr/$CASE
    submit
    exit
fi
# CHAIN: a successor that waits for this job; a case that is done (or failed) has nothing left to do
if [ -n "${CHAIN:-}" ]; then
    tail -1 $KD/refcorr/$CASE/iter.log 2> /dev/null | grep -q "^== done\|^not converged\|failed" && {
        echo "gia_refcorr $CASE: nothing left to do ($(tail -1 $KD/refcorr/$CASE/iter.log))"; exit 0; }
    [ "$CHAIN" -gt 0 ] && echo "gia_refcorr $CASE: successor $(CHAIN=$((CHAIN - 1)) submit --dependency=afterany:$SLURM_JOB_ID)"
fi
HIST=${HIST:-$F/Global/GLOBAL-0.25DEG/GLOBAL-0.25DEG_IceHistory-PaleoMIST.nc}; INTERP=${INTERP:-linear}
RHOI=${RHOI:-910.0}; RHOW=${RHOW:-1028.0}; TOL=${TOL:-1.5}; NMAX=${NMAX:-8}; TMIN=${TMIN:-120}
NCPU=${NCPU:-128}; NTHR=${NTHR:-$(( NCPU < 64 ? NCPU : 64 ))}
XGRID=${XGRID:-GLOBAL-0.5DEG}; XHIST=${XHIST:-$F/Global/$XGRID/${XGRID}_IceHistory-$(basename $HIST | sed 's/.*_IceHistory-//')}
[ -n "${ZREF:-}" ] && XZREF=${XZREF:-$(dirname $(dirname $ZREF))/$XGRID/${XGRID}_IceHistory-$(basename $ZREF | sed 's/.*_IceHistory-//')}
DOMAINS=${DOMAINS:-"North:NH-16KM:North/NH-16KM/NH-16KM_TOPO-GEBCO2026.nc North:NH-32KM:North/NH-32KM/NH-32KM_TOPO-GEBCO2026.nc Antarctica:ANT-16KM:Antarctica/ANT-16KM/ANT-16KM_TOPO-BedMachine-v4.nc"}
VX=$VILMA2/bin/vilma.x
[ -x $VX ] || { echo "no $VX (VILMA2)"; exit 1; }
O=$KD/refcorr/$CASE; mkdir -p $O/data
gitrev() {   # the commit of checkout $1 (without git, as on the compute nodes: from .git)
    git -C $1 rev-parse HEAD 2> /dev/null && return
    local h; h=$(cat $1/.git/HEAD 2> /dev/null) || return 0
    case $h in
        "ref: "*) cat $1/.git/${h#ref: } 2> /dev/null || sed -n "s| ${h#ref: }\$||p" $1/.git/packed-refs 2> /dev/null ;;
        *) echo $h ;;
    esac
}
# provenance of the case, appended at every start (a restart too): the env, the checkouts (dirty: tracked files changed), the
# par file and vilma.x, the inputs with their version and history attributes
{
    echo "== $(date -Is) on $(hostname)"
    for v in CASE PAR HIST XHIST ZREF XZREF R0 DOMAINS NCPU NTHR TOL NMAX INTERP RHOI RHOW KEEP TMIN CHAIN CHUNK_TEND SLURM_JOB_ID; do
        echo "$v=${!v:-}"
    done
    for d in $KRYOS $VILMA2; do echo "git $d $(gitrev $d) dirty $(git -C $d status --porcelain -uno 2> /dev/null | wc -l)"; done
    md5sum $PAR $VX
    for v in HIST XHIST ZREF R0; do
        f=${!v:-}; [ -n "$f" ] || continue; echo "$v $f"
        command -v ncdump > /dev/null && ncdump -h $f 2> /dev/null | grep -E '^[[:space:]]*:(version|history) = ' | sed 's/^[[:space:]]*/    /; s/\\n.*/ (first line)/'
    done
} >> $O/provenance.txt
cp $PAR $O/par.nml
# the &vilma group of the par file (the geo = vilma2 backend reads it over the stack's defaults)
awk '/^&vilma[ \t]*$/{f=1} f{print} f&&/^\//{exit}' $PAR | sed '1d;$d' > $O/vilma.nml
[ -s $O/vilma.nml ] || { echo "no &vilma group in $PAR"; exit 1; }
echo "== gia_refcorr $CASE: $PAR, $HIST ($INTERP), rho $RHOI/$RHOW, vilma.x $VX, tol $TOL m $(date -Is)" | tee -a $O/iter.log
[ -f $O/data/force.nc ] || $JL $T/gia_refcorr.jl prep $HIST $INTERP $O/data ${ZREF:-} || exit 1

# the folder as a PaleoMIST runme project (pm_project.sh: pm_data -> data/), the iterations its runs of the case vilma_rc.nml
bash $T/pm_project.sh $O $O/data > /dev/null || exit 1
cat > $O/vilma_rc.nml <<EOF
&ctl
 file_forcing = "pm_data/force.nc"
 name_h_ice = "h_ice"
 name_time = "time"
 file_ref = ""
 name_z_bed_start = "z_bed_eq"
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
log() { echo "$*" | tee -a $O/iter.log; }
# the job's end [Unix s] (SLURM_JOB_END_TIME), and whether TMIN of it is left (always without one)
JEND=${SLURM_JOB_END_TIME:-}; [[ $JEND =~ ^[0-9]+$ ]] || JEND=""
enough() { [ -z "$JEND" ] || [ $(( JEND - $(date +%s) )) -ge $(( TMIN * 60 )) ]; }
info() { sed -n "s/^$2=//p" $1/chunk.txt 2> /dev/null; }
# the pass it<n> in chunks: 0 done (joined, status 0), 75 to be continued in a later job (TCONT: its model time), 1 failed
pass() {
    local n=$1 R=$O/it$1 k last rc sec t rst tend C ref was
    if [ ! -f $R/vilma_rc.nml ] || [ ! -x $R/vilma.x ]; then
        rm -rf $R; ref=""; [ -n "$prev" ] && ref="-p ctl.i_eq=3 ctl.rsl_restart_file=$prev"
        (cd $O && $KRYOS_RUNME -e vilma -n $O/vilma_rc.nml -o it$n $ref > /dev/null) || { log "it$n: runme failed"; return 1; }
        log "== it$n $(date -Is) on $(hostname)"
    fi
    while :; do
        # the last chunk with a status; one without died: deleted, and redone (c0 from its spin-up restart, if any)
        last=-1; k=0
        while [ -d $R/c$k ]; do
            [ -f $R/c$k/status ] && { last=$k; k=$((k + 1)); continue; }
            if [ $k = 0 ] && [ -f $R/c0/spinup/vilma_restart.nc ]; then rm -rf $R/c0-spinup; mv $R/c0/spinup $R/c0-spinup; fi
            log "it$n/c$k: died (no status), deleted; redone from $( [ $k = 0 ] && { [ -d $R/c0-spinup ] && echo the spin-up restart || echo the start; } || echo "c$((k - 1))/final")"
            while [ -d $R/c$k ]; do rm -rf $R/c$k; k=$((k + 1)); done
            break
        done
        t=""
        if [ $last -ge 0 ]; then
            read -r rc sec t < $R/c$last/status
            case $rc in
                0) [ -z "$(info $R/c$last tend)" ] && break ;;     # the pass's window is done
                75) ;;                                             # stopped before the job's end
                *) log "it$n/c$last failed: exit $rc, $R/c$last/out.err"; return 1 ;;
            esac
        fi
        k=$((last + 1)); C=$R/c$k
        enough || { TCONT=${t:-the start}; log "it$n/c$k: $(( (JEND - $(date +%s)) / 60 )) min left of the job, less than TMIN = $TMIN min: not started"; return 75; }
        rst=""
        if [ $last -ge 0 ]; then rst=c$last/final/vilma_restart.nc; elif [ -f $R/c0-spinup/vilma_restart.nc ]; then rst=c0-spinup/vilma_restart.nc; fi
        [ -z "$rst" ] || [ -f $R/$rst ] || { log "it$n/c$k: no restart $R/$rst"; return 1; }
        tend=""; for x in ${CHUNK_TEND:-}; do awk -v x=$x -v t=${t:--1e30} 'BEGIN{exit !(x > t)}' && { tend=$x; break; }; done
        mkdir -p $C
        sed -e "s|^ *file_out *=.*| file_out = 'c$k/vilma_out.nc'|" -e "s|^ *restart_in_file *=.*| restart_in_file = '$rst'|" \
            -e "s|^ *time_end *=.*| time_end = ${tend:-1e+30}|" $R/vilma_rc.nml > $C/vilma_rc.nml
        [ $(grep -c "^ *\(file_out\|restart_in_file\|time_end\) *=" $C/vilma_rc.nml) = 3 ] || { log "it$n/c$k: $C/vilma_rc.nml lacks a key"; return 1; }
        {   echo "nthr=$NTHR"; echo "cpu=$(sed -n 's/^model name[[:space:]]*: //p' /proc/cpuinfo | head -1)"
            echo "partition=${SLURM_JOB_PARTITION:-}"; echo "node=$(hostname)"; echo "job=${SLURM_JOB_ID:-}"
            echo "vilma2=$(gitrev $VILMA2)"; echo "vilma.x=$(md5sum < $R/vilma.x | cut -c1-32)"; echo "restart=$rst"; echo "tend=$tend"
        } > $C/chunk.txt
        echo "it$n/c$k $(date -Is): $(tr '\n' ' ' < $C/chunk.txt)" >> $O/provenance.txt
        if [ $last -ge 0 ]; then
            for v in nthr cpu; do
                was=$(info $R/c$last $v); [ "$was" = "$(info $C $v)" ] ||
                    log "WARNING it$n/c$k: resumes c$last with another $v ($was -> $(info $C $v)): may not be bit for bit"
            done
        fi
        log "it$n/c$k: from ${rst:-the start}${tend:+ to t = $tend yr} $(date -Is) on $(hostname), $NTHR threads"
        (cd $R && t0=$(date +%s) && ./vilma.x c$k/vilma_rc.nml > c$k/out.out 2> c$k/out.err; rc=$?
         t=$(sed -n 's/^ *t= *\([-0-9.]*\) yr.*/\1/p' c$k/out.out | tail -1)
         echo "$rc $(( $(date +%s) - t0 )) ${t:-?}" > c$k/status.tmp && mv c$k/status.tmp c$k/status)
        read -r rc sec t < $C/status; log "it$n/c$k: exit $rc, $sec s, t = $t yr"
    done
    # the window is done: join the chunks, then the status, then the clean-up (each step safe to redo)
    $JL $T/gia_refcorr.jl join $R | tee -a $O/iter.log; [ ${PIPESTATUS[0]} = 0 ] || { log "it$n: join failed"; return 1; }
    sec=0; for k in $(seq 0 $last); do read -r rc s t < $R/c$k/status; sec=$((sec + s)); done
    for f in out.out out.err; do for k in $(seq 0 $last); do cat $R/c$k/$f; done > $R/$f.tmp && mv $R/$f.tmp $R/$f; done
    [ -d $R/c$last/final ] && { rm -rf $R/final; mv $R/c$last/final $R/final; }    # the layout of an unsplit pass
    if [ -d $R/c0/spinup ]; then rm -rf $R/spinup; mv $R/c0/spinup $R/spinup
    elif [ -d $R/c0-spinup ]; then rm -rf $R/spinup; mv $R/c0-spinup $R/spinup; fi
    echo "0 $sec" > $R/status.tmp && mv $R/status.tmp $R/status
    log "it$n: exit 0, $sec s in $((last + 1)) chunks"
    rm -rf $R/c[0-9]*
}
prev=${R0:-}; conv=0
for n in $(seq 0 $((NMAX - 1))); do
  R=$O/it$n; r=$O/r_$((n + 1)).nc
  if [ "$(cut -d' ' -f1 $R/status 2> /dev/null)" != 0 ]; then
    pass $n; rc=$?
    [ $rc = 75 ] && { log "it$n: continues at t = $TCONT yr in a later job"; exit 0; }
    [ $rc = 0 ] || exit 1
  fi
  rm -rf $R/c[0-9]*     # the chunks of a finished pass (left by a job cut off while deleting them)
  $JL $T/gia_refcorr.jl next $R/vilma_out.nc $r $prev | tee -a $O/iter.log
  m=$(tail -1 $O/iter.log | sed -n 's/.*max |d| = \([0-9.]*\).*/\1/p')
  if [ -n "$prev" ] && awk -v m=$m -v t=$TOL 'BEGIN{exit !(m < t)}'; then conv=1; final=$prev; echo "converged: it$n (r_$n), max misfit $m m" | tee -a $O/iter.log; break; fi
  [ -n "${R0:-}" ] && { echo "== done (one run with r = $R0) $(date -Is)" | tee -a $O/iter.log; exit 0; }
  prev=$r
done
[ $conv = 1 ] || { echo "not converged after $NMAX runs" | tee -a $O/iter.log; exit 1; }
ln -sf $(basename $final) $O/r.nc
if [ "${KEEP:-all}" = last ]; then
    for d in $O/it*; do
        [ $d = $O/it$n ] || rm -f $d/vilma_out.nc; rm -rf $d/spinup $d/final
    done
    echo "KEEP=last: kept it$n/vilma_out.nc of the model output" | tee -a $O/iter.log
fi
DOMS=""; for d in $DOMAINS; do IFS=: read dom hub topo <<< "$d"; DOMS="$DOMS $dom:$hub:$F/$topo"; done
RHO_ICE=$RHOI RHO_SW=$RHOW $JL $T/gia_refcorr.jl write $O/r.nc $O $XGRID:$XHIST${XZREF:+:$XZREF} $DOMS | tee -a $O/iter.log || exit 1
echo "== done $(date -Is)" | tee -a $O/iter.log
