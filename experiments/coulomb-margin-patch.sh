#!/bin/bash
# Greenland's margin from PaleoGrIS 1.0 as a margin patch (margin_patches; docs/backend-plastic.md) on C3's best run,
# tau-CBS-H1000 | E1D-ES1-E117 (experiments/coulomb-band.sh split), which is A:
#
#   experiments/coulomb-margin-patch.sh setup     the patch (tools/levante/margin_patch.jl, if missing), run folder, par files
#   experiments/coulomb-margin-patch.sh submit    runme jobs (about 25 min each on 64 threads)
#   experiments/coulomb-margin-patch.sh score     a score folder with C3's config and error model; references C3's, A and B
#
#   tau-CBS-H1000-pPM | E1D-ES1-E117   the identity check: a patch equal to the base (PaleoMIST's own file in Greenland's
#                                      region) must give A bit for bit
#   tau-CBS-H1000-pPG | E1D-ES1-E117   B: Greenland's margin from PaleoGrIS from 15 to 5 ka (the patch of margin_patch.jl in
#                                      region 2, Greenland, of the maps v2), PaleoMIST's elsewhere and at other times
#   tau-CBS-H1000-pFP | E1D-ES1-E117   the footprint control: B's patch with its slices moved 20 kyr later (never used, its
#                                      footprint only): equals B up to 15 ka if the patch's footprint is all that parts B
#                                      from A there
#
# Both take C3's par file unchanged (r, c_f, regions, output schedule) with ice_nor's margin_patches = 'margin_grl' and the
# group &margin_grl, and this checkout's kryos.x and input/ (the margin_patch template). Env: C3 (C3's run folder, default
# $KRYOS_WORK/coulomb-taub/split/runs), PATCH (default $KRYOS_WORK/margins/NH-32KM/NH-32KM_IceHistory-PaleoGrIS-patch.nc),
# W (default $KRYOS_WORK/coulomb-taub/margin-patch). Runs in $W/runs/<name>, scores in $W/score/iter_01/scores.csv.
set -euo pipefail
X=$(cd "$(dirname "$(readlink -f "$0")")" && pwd); T=$X/../tools/levante; . "$T/env.sh"
[ $# -ge 1 ] || { sed -n '2,22p' "$0"; exit 1; }
C3=${C3:-$KRYOS_WORK/coulomb-taub/split/runs}
PATCH=${PATCH:-$KRYOS_WORK/margins/NH-32KM/NH-32KM_IceHistory-PaleoGrIS-patch.nc}
W=${W:-$KRYOS_WORK/coulomb-taub/margin-patch}; R=$W/runs
BASE=tau-CBS-H1000_E1D-ES1-E117
NAMES=(tau-CBS-H1000-pPM_E1D-ES1-E117 tau-CBS-H1000-pPG_E1D-ES1-E117 tau-CBS-H1000-pFP_E1D-ES1-E117)
FUTURE=$W/NH-32KM_IceHistory-PaleoGrIS-patch-future.nc

setup() {
    [ -f "$PATCH" ] || kryos_jl "$T/margin_patch.jl" "$FESMDATA/GreenlandPaleo/GRL-PAL-1KM/GRL-PAL-1KM_IceHistory-PaleoGrIS.nc" \
        age_deglac 2 15000 5000 "$PATCH"
    mkdir -p "$W"; cp "$PATCH" "$FUTURE"; chmod u+w "$FUTURE"
    kryos_jl -e 'using NCDatasets; NCDataset(ARGS[1], "a") do ds; ds["time"][:] = ds["time"][:] .- 20000; end' "$FUTURE"
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
    rm -rf input; cp -r "$KRYOS/input" .; cp "$KRYOS/libkryos/bin/kryos.x" .
    ln -sfn "$FESMDATA" fesmdata; ln -sfn "$KRYOS_DATA" kryos_data
    rm -rf maps; cp -r "$C3/maps" maps
    for n in "${NAMES[@]}"; do
        case $n in
            *-pPM_*) f="fesmdata/{domain}/{grid_name}/{grid_name}_IceHistory-PaleoMIST.nc" ;;
            *-pFP_*) f=$FUTURE ;;
            *)       f=$PATCH ;;
        esac
        python3 - "$C3/par_$BASE.nml" "par_$n.nml" "$f" << 'EOF'
import re, sys
base, out, patch = sys.argv[1:]
lines = open(base).read().split("\n")
i = next(k for k, l in enumerate(lines) if l.strip().lower() == "&ice_nor")
j = next(k for k in range(i + 1, len(lines)) if lines[k].strip() == "/")
assert not any(re.match(r"\s*margin_patches\s*=", l, re.I) for l in lines[i:j]), "margin_patches already set"
lines.insert(j, " %-15s = 'margin_grl'" % "margin_patches")
lines += ["&margin_grl", " file        = '%s'" % patch, " var         = 'f_ice'", " H_min       = 0.5",
          " region_file = 'kryos_data/plastic/{grid_name}_taub_regions_v2.nc'", " region_var  = 'region'",
          " region_ids  = '2'", "/", ""]
open(out, "w").write("\n".join(lines))
EOF
        echo "par_$n.nml"
    done
}

submit() {
    cd "$R"
    for n in "${NAMES[@]}"; do
        [ -f "par_$n.nml" ] || { echo "no par_$n.nml: run setup"; exit 1; }
        $KRYOS_RUNME -r -s -e kryos -n "par_$n.nml" -o "$n" -q "$KRYOS_QUEUE" -w 01:00:00 --omp 64 -m 30G --jobname "${n#tau-CBS-}" 2>&1 | tail -1
    done
}

# C3's score folder (config, error model, the wave's member) with the references pub-PMIST, pub-GLAC1D (the bar of the cuts),
# ice-PMIST | E1D-ES1-E117 (the r's case), A (C3's run) and B
score() {
    local S=$W/score C3S=$C3/../score
    mkdir -p "$S/iter_01"; cp "$C3S/errmodel.toml" "$S/"; cp "$C3S/iter_01/theta.csv" "$S/iter_01/"
    ln -sfn "$(readlink -f "$C3S/iter_01/0")" "$S/iter_01/0"
    python3 - "$C3S/config.toml" "$S/config.toml" "$C3/$BASE" "$R/${NAMES[1]}" << 'EOF'
import re, sys
src, out, a, b = sys.argv[1:]
txt = open(src).read()
i = txt.index("[[references]]")
head, refs = txt[:i], txt[i:]
assert all(h == "[[references]]" for h in re.findall(r"^\s*\[\[?[^\]]+\]\]?", refs, re.M)), "tables after [[references]]"
keep = [b_ for b_ in refs.split("[[references]]") if re.search(r'name = "(pub-PMIST|pub-GLAC1D|ice-PMIST \| E1D-ES1-E117)"', b_)]
new = [("tau-CBS-H1000 | E1D-ES1-E117", a), ("tau-CBS-H1000-pPG | E1D-ES1-E117", b)]
open(out, "w").write(head + "".join("[[references]]" + b_ for b_ in keep) +
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
