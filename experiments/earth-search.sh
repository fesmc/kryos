#!/bin/bash
# The broad Earth search (docs/dev/earth-search-design.md; docs/simulations.md: ice-<ICE> | E1D-ES1-E<nnn>, E1D-ES2-E<nnn>,
# E3D-ES2-E<nnn>): PaleoMIST's and GLAC-1D's ice imposed on sampled Earths, wave A 200 1-D Earths of 5 radial parameters, wave B
# 50 3-D Earths of 13 parameters, each scored and emulated (plastic_earth_emulate.jl).
#
#   experiments/earth-search.sh sample-a         the 200 1-D Earths (seed 1) and Gowan's Earth in the parameters: $W/earths, $W/cases.txt
#   experiments/earth-search.sh run-a [njob]     wave A: the 402 cases in pools of 16 x 8 cores (default 9 jobs; resubmit to go on)
#   experiments/earth-search.sh score [b]        a job: the scoring folder, the GNSS error model refitted on wave A, scores.csv, the
#                                                emulators (all groups; without GNSS), with b the wave-B emulator too
#   experiments/earth-search.sh sample-b         the 50 3-D Earths from the space wave A left under either ice: $W/earthsB
#   experiments/earth-search.sh run-b 1d [njob]  wave B: the 1-D runs of the 50 profiles (pools of 4 x 8 cores, default 25 jobs)
#   experiments/earth-search.sh run-b 3d [njob]  wave B: the 3-D runs, one pass with the r of the 1-D run, 6-11 h on 32 cores
#                                                (pools of 2 x 32 cores for 36 h, default 35 jobs, my jobs kept below 50); after run-b 1d
#
# Env: W (default $KRYOS_WORK/earth-search; the search of 2026-11-13 is /albedo/work/user/alrobi001/kryos-search), KD (the
# cases: default $W/kryos_data, so $KD/refcorr/<case>), ENS (the ens-v2 folder whose members and config the scores use: linked
# read-only), ERRMODEL (the error model of the split RSL groups: refs-test's refit), GLAC1D (the folder of the GLAC1D_80k
# IceHistory products, until FesmData has them). The 1-D cases do not count against the job limit; the 3-D ones do.
# GLAC1D_80k has no z_bed: today's bed comes from PaleoMIST (ZREF; XZREF derived by gia_refcorr.sh).
set -euo pipefail
X=$(cd "$(dirname "$(readlink -f "$0")")" && pwd); T=$X/../tools/levante; . "$T/env.sh"
[ $# -ge 1 ] || { sed -n '2,22p' "$0"; exit 1; }
J="$KRYOS_JULIA --project=$T"
W=${W:-$KRYOS_WORK/earth-search}; export KD=${KD:-$W/kryos_data}
ENS=${ENS:-/albedo/work/user/alrobi001/kryos-plastic/kryos/runs/ens-v2}
ERRMODEL=${ERRMODEL:-/albedo/work/user/alrobi001/kryos-work/refs-test/ens/errmodel.toml}
GLAC1D=${GLAC1D:-/albedo/work/user/alrobi001/fesmdata-glac1d-80k/Global}
BASIS=${BASIS:-$KRYOS_WORK/kryos_data/plastic/earth_basis_v2.nc}; REGIONS=${REGIONS:-$KRYOS_DATA/plastic/earth_regions_v2.nc}
P1=$KRYOS/par/kryos_NH-ANT_plastic_vilma2_b1.nml; P3=$KRYOS/par/kryos_NH-ANT-PAT_plastic_vilma2.nml
GH="HIST=$GLAC1D/GLOBAL-0.25DEG/GLOBAL-0.25DEG_IceHistory-GLAC1D_80k.nc XHIST=$GLAC1D/GLOBAL-0.5DEG/GLOBAL-0.5DEG_IceHistory-GLAC1D_80k.nc ZREF=$FESMDATA/Global/GLOBAL-0.25DEG/GLOBAL-0.25DEG_IceHistory-PaleoMIST.nc"
export KEEP=last DOMAINS="North:NH-32KM:North/NH-32KM/NH-32KM_TOPO-GEBCO2026.nc"
mkdir -p "$W" "$KD/refcorr"

# the scoring folder: a linked copy of ens-v2 (members read-only) with its config (ens-v2's: the split RSL groups, the GNSS groups
# and ESL of plastic_ensemble_v2.toml, ESL keeping the bar) and errmodel ($ERRMODEL; GNSS refitted by score), and as references
# the bar, I1/I5, the 3-D Bagge references and every case of the lists; config_nognss.toml: the same without GNSS in the cut
score_folder() {
    W=$W KD=$KD ENS=$ENS ERRMODEL=$ERRMODEL V2=$T/plastic_ensemble_v2.toml python3 - <<'EOF'
import os, re
W, KD, E, V2 = os.environ["W"], os.environ["KD"] + "/refcorr", os.environ["ENS"], os.environ["V2"]
C = W + "/score/ens"; os.makedirs(C + "/iter_01", exist_ok=True)
for f in os.listdir(E):
    if f not in ("config.toml", "errmodel.toml", "iter_01") and not os.path.lexists(C + "/" + f): os.symlink(E + "/" + f, C + "/" + f)
skip = ("scores.csv", "weights.csv", "weights.md", "references.csv", "plausible.txt", "weights_union.md")
for f in os.listdir(E + "/iter_01"):
    if f not in skip and not os.path.lexists(C + "/iter_01/" + f): os.symlink(E + "/iter_01/" + f, C + "/iter_01/" + f)
s = open(E + "/config.toml").read(); s = s[:s.index("[[references]]")]
s = re.sub(r'^obs_groups = .*$', 'obs_groups = ["Laurentide", "Greenland", "Fennoscandia", "Barents", "Antarctica", "far field"]', s, flags=re.M)
s = s.replace("[ensemble]\n", "[ensemble]\nfixed_cut = []                           # ESL keeps the bar: a property of the imposed ice, not the Earth\n", 1)
s = s.replace('[likelihood."Fennoscandia/Barents"]', '[likelihood.Fennoscandia]', 1)
v2 = open(V2).read()
lk = re.findall(r'^\[likelihood\."(?:GNSS [^"]+|ESL|Barents)"\].*?(?=^\[|\Z)', v2, flags=re.M | re.S)
esl = re.search(r'^\[esl\].*?(?=^\[)', v2, flags=re.M | re.S).group(0)
gn = re.findall(r'^\[\[gnss\]\].*?(?=^\[|^#|\Z)', v2, flags=re.M | re.S)
s += "\n" + "".join(lk) + "\n" + esl + "\n" + "".join(gn) + "\n"
w = "/albedo/work/user/alrobi001/kryos-work/kryos_data/refcorr"
refs = [("pub-PMIST", 'esl     = "PaleoMIST"\nproduct = "Global/GLOBAL-0.5DEG/GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc"'),
        ("pub-GLAC1D", 'topo = "Global/GLOBAL-0.5DEG/GLOBAL-0.5DEG_IceHistory-GLAC1D.nc"'),
        ("ice-PMIST | E1D-PMIST", 'run  = "/albedo/work/user/alrobi001/kryos-gnss/kryos_data/refcorr/ice-PMIST_E1D-PMIST"'),
        ("ice-GLAC1D | E1D-PMIST", 'run  = "%s/GLAC1D-1D"' % w), ("ice-PMIST | E3D-B21-6", 'run  = "%s/PaleoMIST-3D"' % w),
        ("ice-GLAC1D | E3D-B21-6", 'run  = "%s/GLAC1D-3D"' % w)]
for f in ("cases.txt", "casesB1.txt", "casesB3.txt"):
    if os.path.exists(W + "/" + f):
        for l in open(W + "/" + f):
            c = l.split()[0]; ice, earth = c.split("_", 1)
            refs.append((ice + " | " + earth, 'run  = "%s/%s%s"' % (KD, c, "/it0" if "E3D-ES2" in c else "")))
s += "".join('[[references]]\nname = "%s"\n%s\n' % r for r in refs)
open(C + "/config.toml", "w").write(s)
og = '["Laurentide", "Greenland", "Fennoscandia", "Barents", "Antarctica", "far field", "ESL"]'
open(C + "/config_nognss.toml", "w").write(s.replace("[ensemble]\n", "[ensemble]\nuse_groups = " + og + "\n", 1))
open(C + "/errmodel.toml", "w").write("# copied from " + os.environ["ERRMODEL"] + "\n" + open(os.environ["ERRMODEL"]).read())
print("scoring folder", C, len(refs), "references")
EOF
}
# a pool of the cases of a list (gia_refcorr_pool.sh)
pools() { LIST=$1 NJOB=$2 bash "$T/gia_refcorr_pool.sh"; }

case $1 in
sample-a)
    [ -f "$W/earths/members.csv" ] || $J "$T/plastic_ensemble.jl" earth-search-sample "$P1" "$W/earths" n=200 seed=1
    mkdir -p "$W/check"
    g=$(python3 -c 'import math; print(math.log10(4e20), math.log10(4e20), math.log10(4e22), math.log10(4e22))')
    $J "$T/plastic_ensemble.jl" earth-search-par "$P1" "$W/check/E1D-ES1-PMIST.nml" 120 $g
    { echo "ice-PMIST_E1D-ES1-PMIST $W/check/E1D-ES1-PMIST.nml"; echo "ice-GLAC1D_E1D-ES1-PMIST $W/check/E1D-ES1-PMIST.nml $GH"
      for m in $(seq -f "%03g" 1 200); do
          echo "ice-PMIST_E1D-ES1-E$m $W/earths/E$m.nml"; echo "ice-GLAC1D_E1D-ES1-E$m $W/earths/E$m.nml $GH"
      done; } > "$W/cases.txt"
    wc -l "$W/cases.txt" ;;
run-a)
    NWORK=16 MEM=100G pools "$W/cases.txt" ${2:-9} ;;
score)
    score_folder
    sbatch --parsable -A $ACCOUNT $KRYOS_PART -t 06:00:00 -c 16 --mem=64G -J earth-search-score -o "$W/score/score-%j.out" --wrap \
        "export JULIA_NUM_THREADS=16; set -e
         $J $T/plastic_ensemble.jl errmodel $W/score/ens 1 6 'refs=E1D-ES1-E' 'groups=GNSS Fennoscandia;GNSS Hudson Bay'
         $J $T/plastic_ensemble.jl scores $W/score/ens 1
         $J $T/plastic_earth_emulate.jl $W/earths/members.csv $W/score/ens/iter_01/scores.csv $W/score/ens/config.toml $W/emulate
         $J $T/plastic_earth_emulate.jl $W/earths/members.csv $W/score/ens/iter_01/scores.csv $W/score/ens/config_nognss.toml $W/emulate-nognss
         [ '${2:-}' != b ] || $J $T/plastic_earth_emulate.jl wave-b $W/earths/members.csv $W/earthsB/members.csv $W/score/ens/iter_01/scores.csv \
             $W/score/ens/config.toml $BASIS $W/emulateB" ;;
sample-b)
    $J "$T/plastic_ensemble.jl" earth-search-sample3d "$W/emulate/sample.csv" "$BASIS" "$REGIONS" "$P3" "$P1" "$W/earthsB" n=50 seed=1 space=either
    : > "$W/casesB1.txt"; : > "$W/casesB3.txt"
    for m in $(seq -f "%03g" 1 50); do
        echo "ice-PMIST_E1D-ES2-E$m $W/earthsB/E${m}_1d.nml" >> "$W/casesB1.txt"
        echo "ice-GLAC1D_E1D-ES2-E$m $W/earthsB/E${m}_1d.nml $GH" >> "$W/casesB1.txt"
        echo "ice-PMIST_E3D-ES2-E$m $W/earthsB/E$m.nml R0=$KD/refcorr/ice-PMIST_E1D-ES2-E$m/r.nc" >> "$W/casesB3.txt"
        echo "ice-GLAC1D_E3D-ES2-E$m $W/earthsB/E$m.nml R0=$KD/refcorr/ice-GLAC1D_E1D-ES2-E$m/r.nc $GH" >> "$W/casesB3.txt"
    done ;;
run-b)
    case ${2:?1d or 3d} in
    1d) NCPU=32 NWORK=4 MEM=24G pools "$W/casesB1.txt" ${3:-25} ;;
    3d) for j in $(seq 1 ${3:-35}); do
            while [ $(squeue -u $USER -h | wc -l) -ge 49 ]; do sleep 120; done
            NCPU=64 NWORK=2 MEM=48G PART="$KRYOS_PART_LONG" TIME=36:00:00 pools "$W/casesB3.txt" 1
        done ;;
    esac ;;
*)  sed -n '2,22p' "$0"; exit 1 ;;
esac
