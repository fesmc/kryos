#!/bin/bash
# The joint ice + Earth ensemble (docs/science/joint-ensemble.md): C7's coupled run (tau-CBF-H300 | E1D-ES1-E117,
# experiments/coulomb-band.sh file; PlasticIceSheet online on VILMA2, the Coulomb law, 80 ka to 0) with the Coulomb law's
# H_N, lambda_min (north, Antarctica) and c_f per ice sheet, the LGM timing of the margin in North America and Eurasia, and
# the Earth sampled together (experiments/joint-ens.toml), Eurasia's margin from DATED-1 and Greenland's from PaleoGrIS
# throughout; scored on the 9 groups and emulated (history matching).
#
#   experiments/joint-ens.sh calib              the c_f calibration surface: lgm-calib-v4 (its plastic_coulomb_calib.x and
#                                               namelists) on a grid over the prior, $W/calib, checked against C7's c_f
#   experiments/joint-ens.sh earths             wave A's 1-D Earths not ruled out under either ice history ($W/earths.csv),
#                                               their par files ($W/earths) and r with the reference files of every domain
#                                               ($W/refcorr/ice-PMIST_E1D-ES1-<earth>, gia_refcorr.jl write): a job
#   experiments/joint-ens.sh margins            Greenland's PaleoGrIS patch (C4's), $W/margins: shared by every member
#   experiments/joint-ens.sh design <set> [<points.csv>]   the design, $W/<set>/design.csv: the maximin LHS of joint-ens.toml,
#                                               or the members of <points.csv> (the smoke test: experiments/joint-ens-smoke.csv)
#   experiments/joint-ens.sh refocus <set> <emulate dir> [n]   a next wave in what a wave did not rule out: n (120) points of the
#                                               emulator's sample with Imax <= 3, spread by maximin (plastic_joint_ens.jl nroy),
#                                               $W/<set>/points.csv, and their design
#   experiments/joint-ens.sh setup <set>        per member the par file, c_f file and margin patches (a job), the run folder
#   experiments/joint-ens.sh submit <set> [<member> ...]   runme jobs (default every member), my jobs kept below 50
#   experiments/joint-ens.sh score <set> [<folder>]   a score folder ($W/<set>/<folder>, default score, a new one): the groups,
#                                               data and error model of tools/levante/plastic_ensemble_v2.toml (the scores of the
#                                               checkout), and as references the bar (pub-PMIST, pub-GLAC1D), the published ice on
#                                               Gowan's Earth (I1, I5: the bar of the far field LGM), C7, C8 pDE and the members
#   experiments/joint-ens.sh post <set> [plots|-] [<folder>]   the emulator (plastic_joint_emulate.jl) of the scores of <folder>
#                                               (default score) in $W/<set>/emulate (emulate-<folder> for another); plots: the jumps
#                                               and volumes per ice sheet of the members (plastic_mass_jumps.jl,
#                                               plastic_coulomb_coupled.jl; a small set)
#
# The member <m> of <set> runs as tau-JE-<m> | E1D-ES1-<earth> in $W/<set>/runs/tau-JE-<m>_E1D-ES1-<earth>, its files in
# $W/<set>/members/<m>: cf.nc, margin_na.nc (region 1: PaleoMIST 27.5-20 ka, its 25 and 22.5 ka slices blended by tau_NA),
# margin_eu.nc (regions 3, 6: DATED-1 most-credible 27.5-7.5 ka, its 25-21 ka slices blended by tau_EU). The par file is C7's
# (TEMPLATE) with H_N, lambda_min, the c_f file, the Earth, its r and the three patches (plastic_joint_ens.jl files); C7's
# kryos.x, input/ and maps. About 25 min per run on 64 threads (C7).
#
# Env: W (default $KRYOS_WORK/joint-ens), C7 (C7's run folder, default $KRYOS_WORK/coulomb-taub/lambda-file/runs), TEMPLATE
# (default $C7/par_tau-CBF-H300_E1D-ES1-E117.nml), V4 (lgm-calib-v4), SEARCH (the Earth search of wave A: earths/, kryos_data/refcorr/;
# default /albedo/work/user/alrobi001/kryos-search), EMU (the wave-A emulator without ESL judging:
# $KRYOS_WORK/rescore-esl/earth-search/emulate-noesl), S (the FesmData products with DATED-1, default
# $KRYOS_WORK/fesmdata-staging/fesmdata), I1, I5 (PaleoMIST's and GLAC-1D's ice on Gowan's Earth, standalone VILMA2: default
# /albedo/work/user/alrobi001/kryos-gnss/kryos_data/refcorr/ice-PMIST_E1D-PMIST, .../kryos-work/kryos_data/refcorr/GLAC1D-1D),
# WALL (01:30:00), NJOBS (50).
set -euo pipefail
X=$(cd "$(dirname "$(readlink -f "$0")")" && pwd); T=$X/../tools/levante; . "$T/env.sh"
[ $# -ge 1 ] || { sed -n '2,42p' "$0"; exit 1; }
W=${W:-$KRYOS_WORK/joint-ens}; C7=${C7:-$KRYOS_WORK/coulomb-taub/lambda-file/runs}
TEMPLATE=${TEMPLATE:-$C7/par_tau-CBF-H300_E1D-ES1-E117.nml}; V4=${V4:-$KRYOS_WORK/coulomb-taub/lgm-calib-v4}
SEARCH=${SEARCH:-/albedo/work/user/alrobi001/kryos-search}; EMU=${EMU:-$KRYOS_WORK/rescore-esl/earth-search/emulate-noesl}
S=${S:-$KRYOS_WORK/fesmdata-staging/fesmdata}
I1=${I1:-/albedo/work/user/alrobi001/kryos-gnss/kryos_data/refcorr/ice-PMIST_E1D-PMIST}
I5=${I5:-/albedo/work/user/alrobi001/kryos-work/kryos_data/refcorr/GLAC1D-1D}
WALL=${WALL:-01:30:00}; NJOBS=${NJOBS:-50}
F=$FESMDATA; J="$KRYOS_JULIA --project=$T"
PM=$F/North/NH-32KM/NH-32KM_IceHistory-PaleoMIST.nc; DATED=$S/North/NH-32KM/NH-32KM_IceHistory-DATED1.nc
CALN=$W/calib/calib_NH-32KM.nc; CALA=$W/calib/calib_ANT-16KM.nc
setdir() { SD=$W/${1:?set}; }
# the members of a design and a column: "member value" lines
col() { python3 - "$1" "$2" << 'EOF'
import csv, sys
for r in csv.DictReader(open(sys.argv[1])): print(r["member"], r[sys.argv[2]])
EOF
}
names() { python3 - "$1" << 'EOF'
import csv, sys
for r in csv.DictReader(open(sys.argv[1])): print(r["member"], "tau-JE-%s_E1D-ES1-%s" % (r["member"], r["earth"]))
EOF
}
# a job that runs the lines of a command file, P at a time, and waits for it
pjob() {   # <command file> <log> <time> <P>
    sbatch --wait --parsable -A $ACCOUNT $KRYOS_PART -t $3 -c $4 --mem=$(( 4 * $4 ))G -J je-$(basename $1 .txt) -o $2 \
        --wrap "xargs -P $4 -I{} bash -c '{}' < $1" > /dev/null
}

case $1 in
calib)
    mkdir -p "$W/calib"; cd "$W/calib"; cp "$V4/plastic_coulomb_calib.x" .
    HN="300.0, 350.0, 420.0, 500.0, 600.0, 700.0, 850.0, 1000.0"
    sed -e "s|^ *HN *=.*|    HN            = $HN|" -e "s|^ *out *=.*|    out           = \"$CALN\"|" \
        -e "s|^ *lam_min *=.*|    lam_min       = 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90|" "$V4/calib_NH.nml" > calib_NH.nml
    sed -e "s|^ *HN *=.*|    HN            = $HN|" -e "s|^ *out *=.*|    out           = \"$CALA\"|" \
        -e "s|^ *lam_min *=.*|    lam_min       = 0.15, 0.20, 0.22, 0.25, 0.30, 0.34, 0.35, 0.40\n    zb0           = -600.0\n    zb1           = 0.0|" \
        "$V4/calib_ANT.nml" > calib_ANT.nml
    # C7's c_f (coulomb-band.sh file): H_N:lam_N:lam_A: NAIS, Greenland, Fennoscandia, BarentsKara, Antarctica
    sbatch --wait -A $ACCOUNT $KRYOS_PART -t 00:40:00 -c 28 --mem=60G -J je-calib -o "$W/calib/calib.out" --wrap \
        "cd $W/calib; module unload intel-oneapi-compilers; module load intel-oneapi-compilers/2024.0.0; ulimit -s unlimited
         export OMP_NUM_THREADS=28 OMP_STACKSIZE=512M; set -e
         for c in NH ANT; do ./plastic_coulomb_calib.x calib_\$c.nml > calib_\$c.log 2>&1; done
         $J $T/plastic_cf_surface.jl check $CALN $CALA $W/calib ref=300:0.5:0.22:0.019699:0.044263:0.021662:0.047985:0.042123 \
             ref=1000:0.75:0.34:0.008536:0.022060:0.010375:0.019428:0.019432"
    cat "$W/calib/cf_surface.md" ;;
earths)
    $J "$T/plastic_joint_ens.jl" earths "$EMU" "$SEARCH/earths/members.csv" "$W/earths.csv"
    mkdir -p "$W/earths" "$W/refcorr"; : > "$W/earths-write.txt"
    for e in $(tail -n +2 "$W/earths.csv" | cut -d, -f1); do
        c=$SEARCH/kryos_data/refcorr/ice-PMIST_E1D-ES1-$e; o=$W/refcorr/ice-PMIST_E1D-ES1-$e
        grep -q "^== done" "$c/iter.log" && grep -q "^converged" "$c/iter.log" || { echo "$c: not converged"; exit 1; }
        cp "$SEARCH/earths/$e.nml" "$W/earths/"
        [ -f "$o/North/NH-32KM/NH-32KM_TOPO-GEBCO2026-refcorr.nc" ] && continue
        mkdir -p "$o"; cp -L "$c/r.nc" "$o/r.nc"
        echo "r: copied from $c/r.nc ($(readlink $c/r.nc); $(grep '^converged' $c/iter.log)); reference files by gia_refcorr.jl write, $(date)" > "$o/provenance.txt"
        echo "RHO_ICE=910.0 RHO_SW=1028.0 $J $T/gia_refcorr.jl write $o/r.nc $o GLOBAL-0.5DEG:$F/Global/GLOBAL-0.5DEG/GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc North:NH-32KM:$F/North/NH-32KM/NH-32KM_TOPO-GEBCO2026.nc Antarctica:ANT-16KM:$F/Antarctica/ANT-16KM/ANT-16KM_TOPO-BedMachine-v4.nc PatagoniaPaleo:PAT-PAL-8KM:$F/PatagoniaPaleo/PAT-PAL-8KM/PAT-PAL-8KM_TOPO-GEBCO2026.nc > $o/write.log 2>&1" >> "$W/earths-write.txt"
    done
    [ -s "$W/earths-write.txt" ] && pjob "$W/earths-write.txt" "$W/earths-write.out" 02:00:00 8
    for o in "$W"/refcorr/*/; do [ -f "$o/North/NH-32KM/NH-32KM_TOPO-GEBCO2026-refcorr.nc" ] || echo "$o: no reference files"; done
    echo "$(ls -d "$W"/refcorr/*/ | wc -l) Earths in $W/refcorr" ;;
margins)
    mkdir -p "$W/margins"
    kryos_jl "$T/margin_patch.jl" "$F/GreenlandPaleo/GRL-PAL-1KM/GRL-PAL-1KM_IceHistory-PaleoGrIS.nc" age_deglac 2 15000 5000 \
        "$W/margins/NH-32KM_IceHistory-PaleoGrIS-patch.nc" ;;
design)
    setdir "$2"; mkdir -p "$SD"
    $J "$T/plastic_joint_ens.jl" design "$X/joint-ens.toml" "$W/earths.csv" "$CALN" "$CALA" "$SD/design.csv" ${3:+points=$3} ;;
refocus)
    setdir "$2"; mkdir -p "$SD"
    $J "$T/plastic_joint_ens.jl" nroy "$X/joint-ens.toml" "$3/sample.csv" "$SD/points.csv" n=${4:-120}
    $J "$T/plastic_joint_ens.jl" design "$X/joint-ens.toml" "$W/earths.csv" "$CALN" "$CALA" "$SD/design.csv" points="$SD/points.csv" ;;
setup)
    setdir "$2"
    $J "$T/plastic_joint_ens.jl" files "$SD/design.csv" "$TEMPLATE" "$W" "$SD"
    : > "$SD/patches.txt"
    paste -d' ' <(col "$SD/design.csv" tau_NA) <(col "$SD/design.csv" tau_EU | cut -d' ' -f2) | while read -r m na eu; do
        M=$SD/members/$m
        echo "$J $T/margin_patch.jl $PM f_ice 1 27500 20000 $M/margin_na.nc --blend=$na > $M/margin_na.log 2>&1" >> "$SD/patches.txt"
        echo "$J $T/margin_patch.jl $DATED f_ice 3,6 27500 7500 $M/margin_eu.nc --blend=$eu > $M/margin_eu.log 2>&1" >> "$SD/patches.txt"
    done
    pjob "$SD/patches.txt" "$SD/patches.out" 04:00:00 16
    for m in $(col "$SD/design.csv" member | cut -d' ' -f1); do
        for p in margin_na margin_eu; do [ -f "$SD/members/$m/$p.nc" ] || { echo "$m: no $p.nc ($SD/members/$m/$p.log)"; exit 1; }; done
    done
    R=$SD/runs; cd "$R"
    rm -rf .runme input maps; cp -r "$C7/.runme" .runme; cp -r "$C7/input" .; cp "$C7/kryos.x" .; cp -r "$C7/maps" maps
    ln -sfn "$FESMDATA" fesmdata; ln -sfn "$KRYOS_DATA" kryos_data
    echo "$(ls par_*.nml | wc -l) par files in $R; kryos.x, input/, maps of $C7" ;;
submit)
    setdir "$2"; shift 2; cd "$SD/runs"
    while read -r m n; do
        [ $# = 0 ] || [[ " $* " == *" $m "* ]] || continue
        [ -f "par_$n.nml" ] || { echo "no par_$n.nml: run setup"; exit 1; }
        [ -e "$n" ] && { echo "$n exists: skipped"; continue; }
        while [ $(squeue -u $USER -h | wc -l) -ge $NJOBS ]; do sleep 120; done
        $KRYOS_RUNME -r -s -e kryos -n "par_$n.nml" -o "$n" -q "$KRYOS_QUEUE" -w "$WALL" --omp 64 -m 30G --jobname "JE-$m" 2>&1 | tail -1
    done < <(names "$SD/design.csv") ;;
score)
    setdir "$2"; SC=$SD/${3:-score}
    [ -e "$SC/config.toml" ] && { echo "$SC/config.toml exists: give a new folder"; exit 1; }
    mkdir -p "$SC/iter_01"; cp "$C7/../score/iter_01/theta.csv" "$SC/iter_01/"
    C8=$KRYOS_WORK/coulomb-taub/dated-margins/runs/tau-CBF-H300-pDE_E1D-ES1-E117
    python3 - "$T/plastic_ensemble_v2.toml" "$SC/config.toml" "$SD/runs" <(names "$SD/design.csv") \
        "tau-CBF-H300 | E1D-ES1-E117=$C7/tau-CBF-H300_E1D-ES1-E117" "tau-CBF-H300-pDE | E1D-ES1-E117=$C8" \
        "ice-PMIST | E1D-PMIST=$I1" "ice-GLAC1D | E1D-PMIST=$I5" << 'EOF'
import os, re, sys
src, out, r, nf = sys.argv[1:5]; runs = [tuple(a.split("=", 1)) for a in sys.argv[5:]]
txt = open(src).read(); i = txt.index("[[references]]"); head, refs = txt[:i], txt[i:]
assert all(h == "[[references]]" for h in re.findall(r"^\s*\[\[?[^\]]+\]\]?", refs, re.M)), "tables after [[references]]"
# the paths of the config ({work}, {fesmdata}, {kryos_data}, {kryos}: plastic_esmda.jl expand_paths) from env.sh
for k, v in (("{work}", "KRYOS_WORK"), ("{fesmdata}", "FESMDATA"), ("{kryos_data}", "KRYOS_DATA"), ("{kryos}", "KRYOS")):
    head = head.replace(k, os.environ[v])
pubs = [b for b in refs.split("[[references]]") if re.search(r'name *= *"(pub-PMIST|pub-GLAC1D)"', b)]
assert len(pubs) == 2, "pub-PMIST, pub-GLAC1D not both in " + src
new = runs + [(n.replace("_", " | ", 1), "%s/%s" % (r, n)) for _, n in (l.split() for l in open(nf))]
open(out, "w").write("# joint-ens.sh score: %s, with the bar, I1, I5, C7, C8 pDE and the members as references\n" % src + head +
                     "".join("[[references]]" + b for b in pubs) + "".join('[[references]]\nname = "%s"\nrun  = "%s"\n' % nd for nd in new))
EOF
    "$X/scores.sh" "$SC" ;;
post)
    setdir "$2"; F=${4:-score}; SC=$SD/$F; E=$SD/emulate; [ "$F" = score ] || E=$SD/emulate-$F
    JL_TIME=02:00:00 JL_MEM=32G "$T/jl.sh" "$E.log" "$T/plastic_joint_emulate.jl" "$SD/design.csv" "$X/joint-ens.toml" \
        "$W/earths.csv" "$SC/iter_01/scores.csv" "$SC/config.toml" "$E"
    if [ "${3:-}" = plots ]; then
        L=("tau-CBF-H300=$C7/tau-CBF-H300_E1D-ES1-E117")
        while read -r m n; do [ -f "$SD/runs/$n/ice_nor_ts.nc" ] && L+=("tau-JE-$m=$SD/runs/$n"); done < <(names "$SD/design.csv")
        JL_TIME=00:30:00 "$T/jl.sh" "$SD/jumps.log" "$T/plastic_mass_jumps.jl" "$SD/jumps" "${L[@]}"
        JL_TIME=00:40:00 "$T/jl.sh" "$SD/coupled.log" "$T/plastic_coulomb_coupled.jl" "$SD/coupled" "${L[@]}"
    fi ;;
*)  sed -n '2,42p' "$0"; exit 1 ;;
esac
