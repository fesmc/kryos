#!/bin/bash
# Greenland's margin from PaleoGrIS 1.0 under the Coulomb law: an A/B pair of coupled runs, C2's tau-CB-H1000 | E1D-ES1-E117
# (experiments/coulomb-band.sh) with the margin chronology of PaleoMIST (A, the same run again) and with the merged chronology
# of tools/levante/margin_merge.jl (B: PaleoGrIS in Greenland 13.5-6.75 ka, PaleoMIST elsewhere):
#
#   experiments/coulomb-paleogris.sh setup     the merged margin (if missing), run folder and par files of A and B
#   experiments/coulomb-paleogris.sh submit    runme jobs of A and B (about 10 min each on 64 threads)
#   experiments/coulomb-paleogris.sh score     a score folder with C2's config and error model, references C2, A, B; scores.sh
#
#   tau-CB-H1000 | E1D-ES1-E117          A: C2's par file, only out_2D changed
#   tau-CB-H1000-mPG | E1D-ES1-E117      B: and ice_nor's margin_file = the merged chronology
#
# Both take C2's par file, kryos.x and input/ (the C2 run folder), so that only the margin differs between A and B, and A
# reproduces C2 at C2's records. The r (topography correction) is C2's, made for PaleoMIST's ice on E117: borrowed by B. The
# Coulomb regions are C2's ($KRYOS_DATA/plastic/NH-32KM_taub_regions.nc). out_2D: every 2.5 kyr to 30 ka, then every 500 yr
# (the scores need records at the 2.5 kyr ages; 500 yr resolves the 0.5 kyr isochrones).
#
# Env: CB (C2's run folder, default $KRYOS_WORK/coulomb-taub/coupled/runs), MARGIN (the merged chronology, default
# $KRYOS_WORK/margins/NH-32KM/NH-32KM_IceHistory-PaleoMIST-PaleoGrIS.nc), W (default $KRYOS_WORK/coulomb-taub/paleogris).
# Runs in $W/runs/<name>, scores in $W/score/iter_01/scores.csv.
set -euo pipefail
X=$(cd "$(dirname "$(readlink -f "$0")")" && pwd); T=$X/../tools/levante; . "$T/env.sh"
[ $# -ge 1 ] || { sed -n '2,23p' "$0"; exit 1; }
CB=${CB:-$KRYOS_WORK/coulomb-taub/coupled/runs}
MARGIN=${MARGIN:-$KRYOS_WORK/margins/NH-32KM/NH-32KM_IceHistory-PaleoMIST-PaleoGrIS.nc}
W=${W:-$KRYOS_WORK/coulomb-taub/paleogris}; R=$W/runs
BASE=tau-CB-H1000_E1D-ES1-E117
NAMES=(tau-CB-H1000_E1D-ES1-E117 tau-CB-H1000-mPG_E1D-ES1-E117)
OUT_2D="-80000:2500:-30000, -30000:500:0"

setup() {
    [ -f "$MARGIN" ] || kryos_jl "$T/margin_merge.jl" "$MARGIN"
    mkdir -p "$R"; cd "$R"
    rm -rf .runme; cp -r "$KRYOS/.runme" .runme; sed -e 's/^hpc = .*/hpc = "awi_albedo"/' -e "s/^account = .*/account = \"$ACCOUNT\"/" \
        "$KRYOS/.runme/config.default.toml" > .runme/config.toml
    cat > .runme/info.json << 'JSON'
{
    "exe_default" : "kryos", "par_path_as_argument" : true, "exe_aliases" : {"kryos" : "kryos.x"}, "grp_aliases" : {},
    "par_paths" : {}, "files" : [], "dir-special" : {}, "par_defaults" : ["input/kryos_defaults.nml", "input/vilma_defaults.nml"],
    "par_library" : [], "links" : ["input", "fesmdata", "kryos_data", "maps"]
}
JSON
    rm -rf input; cp -r "$CB/input" .; cp "$CB/kryos.x" .
    ln -sfn "$FESMDATA" fesmdata; ln -sfn "$KRYOS_DATA" kryos_data
    mkdir -p maps; cp -n "$KRYOS_DATA"/maps/grid_*.txt "$KRYOS_DATA"/maps/scrip-*.nc maps/ 2> /dev/null || true
    for n in "${NAMES[@]}"; do
        m=""; [ "$n" = "$BASE" ] || m=$MARGIN
        python3 - "$CB/par_$BASE.nml" "par_$n.nml" "$OUT_2D" "$m" << 'EOF'
import re, sys
base, out, times, margin = sys.argv[1:]
lines = open(base).read().split("\n")
def group(g):
    i = next(k for k, l in enumerate(lines) if l.strip().lower() == "&" + g.lower())
    j = next(k for k in range(i + 1, len(lines)) if lines[k].strip() == "/")
    return i, j
def put(g, key, val):
    i, j = group(g)
    for k in range(i + 1, j):
        if re.match(r"\s*" + key + r"\s*=", lines[k], re.I):
            lines[k] = " %-15s = %s" % (key, val); return
    lines.insert(j, " %-15s = %s" % (key, val))
put("out_2D", "method", "'times'")
put("out_2D", "times", "'%s'" % times)
if margin:
    put("ice_nor", "margin_file", "'%s'" % margin)
open(out, "w").write("\n".join(lines))
EOF
        echo "par_$n.nml"
    done
}

submit() {
    cd "$R"
    for n in "${NAMES[@]}"; do
        [ -f "par_$n.nml" ] || { echo "no par_$n.nml: run setup"; exit 1; }
        $KRYOS_RUNME -r -s -e kryos -n "par_$n.nml" -o "$n" -q "$KRYOS_QUEUE" -w 02:00:00 --omp 64 -m 30G --jobname "pg-${n#tau-CB-}" 2>&1 | tail -1
    done
}

# C2's score folder (config, error model, the wave's member) with the references pub-PMIST, pub-GLAC1D (the bar of the cuts),
# ice-PMIST | E1D-ES1-E117 (the r's case), C2 (as scored there), A and B
score() {
    local S=$W/score C2S=$CB/../score
    mkdir -p "$S/iter_01"; cp "$C2S/errmodel.toml" "$S/"; cp "$C2S/iter_01/theta.csv" "$S/iter_01/"
    ln -sfn "$(readlink -f "$C2S/iter_01/0")" "$S/iter_01/0"
    python3 - "$C2S/config.toml" "$S/config.toml" "$CB/$BASE" "$R" << 'EOF'
import re, sys
src, out, c2, R = sys.argv[1:]
txt = open(src).read()
i = txt.index("[[references]]")
head, refs = txt[:i], txt[i:]
assert all(h == "[[references]]" for h in re.findall(r"^\s*\[\[?[^\]]+\]\]?", refs, re.M)), "tables after [[references]]"
keep = [b for b in refs.split("[[references]]") if re.search(r'name = "(pub-PMIST|pub-GLAC1D|ice-PMIST \| E1D-ES1-E117)"', b)]
new = [("C2 tau-CB-H1000 | E1D-ES1-E117", c2), ("tau-CB-H1000 | E1D-ES1-E117", R + "/tau-CB-H1000_E1D-ES1-E117"),
       ("tau-CB-H1000-mPG | E1D-ES1-E117", R + "/tau-CB-H1000-mPG_E1D-ES1-E117")]
open(out, "w").write(head + "".join("[[references]]" + b for b in keep) +
                     "".join('[[references]]\nname = "%s"\nrun  = "%s"\n' % nd for nd in new))
EOF
    "$X/scores.sh" "$S"
}

case $1 in
    setup)  setup ;;
    submit) submit ;;
    score)  score ;;
    *)      echo "unknown command $1 (setup, submit, score)"; exit 1 ;;
esac
