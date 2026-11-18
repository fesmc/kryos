#!/bin/bash
# Dated ice-margin chronologies as margin patches (margin_patches; docs/backend-plastic.md#margin-patches) on C7's best run,
# tau-CBF-H300 | E1D-ES1-E117 (experiments/coulomb-band.sh file): North America from Dalton et al. (2020) and Eurasia from
# DATED-1 (Hughes et al. 2016, most-credible), the FesmData products IceHistory-Dalton2020 and IceHistory-DATED1 on NH-32KM.
#
#   experiments/coulomb-dated-margins.sh patches   the patches (tools/levante/margin_patch.jl) in $W/margins, with their checks
#   experiments/coulomb-dated-margins.sh setup     run folder and par files
#   experiments/coulomb-dated-margins.sh submit    runme jobs (about 25 min each on 64 threads)
#   experiments/coulomb-dated-margins.sh score     a score folder with C7's config and error model; references C7, esl-diag's
#                                                  ED-pNE25 and the runs here
#   experiments/coulomb-dated-margins.sh post      jumps (plastic_mass_jumps.jl), volumes per ice sheet (plastic_coulomb_coupled.jl)
#
#   tau-CBF-H300-pDA | E1D-ES1-E117    NAIS from Dalton et al. (2020) in region 1 (Laurentide) of the maps v2, NA_T0 to NA_T1
#   tau-CBF-H300-pDE | E1D-ES1-E117    Eurasia from DATED-1 (most-credible) in regions 3 and 6 (Fennoscandia, BarentsKara),
#                                      EU_T0 to EU_T1
#   tau-CBF-H300-pDAE | E1D-ES1-E117   both
#
# Each patch: PaleoMIST's slices at both ends of its window, the product's slices strictly inside it (margin_patch.jl).
# Greenland keeps PaleoMIST's margin (no PaleoGrIS patch, as C7). The par files are C7's unchanged with ice_nor's
# margin_patches and their groups added; C7's kryos.x, input/ and maps. Env: C7 (C7's run folder, default
# $KRYOS_WORK/coulomb-taub/lambda-file/runs), S (the FesmData products with the dated margins, default $FESMDATA:
# DATED-1 is released, IceHistory v2.7.0; Dalton2020 is not (no licence stated), so pDA and pDAE need S with it),
# NA_T0, NA_T1 (default 25000, 20000), EU_T0, EU_T1 (default 27500, 7500), W (default $KRYOS_WORK/coulomb-taub/dated-margins).
# Runs in $W/runs/<name>, scores in $W/score/iter_01/scores.csv.
set -euo pipefail
X=$(cd "$(dirname "$(readlink -f "$0")")" && pwd); T=$X/../tools/levante; . "$T/env.sh"
[ $# -ge 1 ] || { sed -n '2,26p' "$0"; exit 1; }
C7=${C7:-$KRYOS_WORK/coulomb-taub/lambda-file/runs}
S=${S:-$FESMDATA}
NA_T0=${NA_T0:-25000}; NA_T1=${NA_T1:-20000}; EU_T0=${EU_T0:-27500}; EU_T1=${EU_T1:-7500}
W=${W:-$KRYOS_WORK/coulomb-taub/dated-margins}; R=$W/runs; M=$W/margins
BASE=tau-CBF-H300_E1D-ES1-E117
NAMES=(tau-CBF-H300-pDA_E1D-ES1-E117 tau-CBF-H300-pDE_E1D-ES1-E117 tau-CBF-H300-pDAE_E1D-ES1-E117)
P_NA=$M/NH-32KM_IceHistory-Dalton2020-patch.nc
P_EU=$M/NH-32KM_IceHistory-DATED1-patch.nc

patches() {
    mkdir -p "$M"
    kryos_jl "$T/margin_patch.jl" "$S/North/NH-32KM/NH-32KM_IceHistory-Dalton2020.nc" f_ice 1 "$NA_T0" "$NA_T1" "$P_NA"
    kryos_jl "$T/margin_patch.jl" "$S/North/NH-32KM/NH-32KM_IceHistory-DATED1.nc" f_ice 3,6 "$EU_T0" "$EU_T1" "$P_EU"
}

setup() {
    for p in "$P_NA" "$P_EU"; do [ -f "$p" ] || { echo "no $p: run patches"; exit 1; }; done
    mkdir -p "$R"; cd "$R"
    rm -rf .runme; cp -r "$C7/.runme" .runme
    rm -rf input; cp -r "$C7/input" .; cp "$C7/kryos.x" .
    ln -sfn "$FESMDATA" fesmdata; ln -sfn "$KRYOS_DATA" kryos_data
    rm -rf maps; cp -r "$C7/maps" maps
    for n in "${NAMES[@]}"; do
        case $n in
            *-pDA_*)  g="margin_na=$P_NA=1" ;;
            *-pDE_*)  g="margin_eu=$P_EU=3,6" ;;
            *-pDAE_*) g="margin_na=$P_NA=1 margin_eu=$P_EU=3,6" ;;
        esac
        python3 - "$C7/par_$BASE.nml" "par_$n.nml" $g << 'EOF'
import re, sys
base, out, groups = sys.argv[1], sys.argv[2], [a.split("=") for a in sys.argv[3:]]
lines = open(base).read().split("\n")
i = next(k for k, l in enumerate(lines) if l.strip().lower() == "&ice_nor")
j = next(k for k in range(i + 1, len(lines)) if lines[k].strip() == "/")
assert not any(re.match(r"\s*margin_patches\s*=", l, re.I) for l in lines[i:j]), "margin_patches already set"
lines.insert(j, " %-15s = %s" % ("margin_patches", " ".join("'%s'" % g for g, _, _ in groups)))
for g, f, ids in groups:
    lines += ["&%s" % g, " file        = '%s'" % f, " var         = 'f_ice'", " H_min       = 0.5",
              " region_file = 'kryos_data/plastic/{grid_name}_taub_regions_v2.nc'", " region_var  = 'region'",
              " region_ids  = %s" % ", ".join("'%s'" % r for r in ids.split(",")), "/", ""]
open(out, "w").write("\n".join(lines))
EOF
        echo "par_$n.nml"; diff "$C7/par_$BASE.nml" "par_$n.nml" || true
    done
}

submit() {
    cd "$R"
    for n in "${NAMES[@]}"; do
        [ -f "par_$n.nml" ] || { echo "no par_$n.nml: run setup"; exit 1; }
        $KRYOS_RUNME -r -s -e kryos -n "par_$n.nml" -o "$n" -q "$KRYOS_QUEUE" -w 01:00:00 --omp 64 -m 30G --jobname "${n#tau-CBF-}" 2>&1 | tail -1
    done
}

# C7's score folder (config, error model) with the references pub-PMIST, pub-GLAC1D (the bar), ice-PMIST | E1D-ES1-E117, C7,
# esl-diag's ED-pNE25 (C7 H300 with PaleoMIST's 20 ka extent from 25 ka in NAIS and Eurasia) and the runs here
score() {
    local SC=$W/score C7S=$C7/../score ED=${ED:-$KRYOS_WORK/coulomb-taub/esl-diag/runs/tau-ED-pNE25-H300_E1D-ES1-E117}
    mkdir -p "$SC/iter_01"; cp "$C7S/errmodel.toml" "$SC/"; cp "$C7S/iter_01/theta.csv" "$SC/iter_01/"
    python3 - "$C7S/config.toml" "$SC/config.toml" "$C7/$BASE" "$ED" "$R" "${NAMES[@]}" << 'EOF'
import re, sys
src, out, c7, ed, r = sys.argv[1:6]; names = sys.argv[6:]
txt = open(src).read()
i = txt.index("[[references]]")
head, refs = txt[:i], txt[i:]
assert all(h == "[[references]]" for h in re.findall(r"^\s*\[\[?[^\]]+\]\]?", refs, re.M)), "tables after [[references]]"
keep = [b for b in refs.split("[[references]]") if re.search(r'name = "(pub-PMIST|pub-GLAC1D|ice-PMIST \| E1D-ES1-E117)"', b)]
new = [("tau-CBF-H300 | E1D-ES1-E117", c7), ("tau-ED-pNE25-H300 | E1D-ES1-E117", ed)]
new += [(n.replace("_", " | ", 1), "%s/%s" % (r, n)) for n in names]
open(out, "w").write(head + "".join("[[references]]" + b for b in keep) +
                     "".join('[[references]]\nname = "%s"\nrun  = "%s"\n' % nd for nd in new))
EOF
    "$X/scores.sh" "$SC"
}

post() {
    local L=("tau-CBF-H300=$C7/$BASE")
    for n in "${NAMES[@]}"; do L+=("${n%_E1D-ES1-E117}=$R/$n"); done
    JL_TIME=00:30:00 "$T/jl.sh" "$W/jumps.log" "$T/plastic_mass_jumps.jl" "$W/jumps" "${L[@]}"
    JL_TIME=00:40:00 "$T/jl.sh" "$W/coupled.log" "$T/plastic_coulomb_coupled.jl" "$W/coupled" "${L[@]}"
}

case $1 in
    patches) patches ;;
    setup)   setup ;;
    submit)  submit ;;
    score)   score ;;
    post)    post ;;
    *)       echo "unknown command $1 (patches, setup, submit, score, post)"; exit 1 ;;
esac
