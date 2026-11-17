#!/bin/bash
# Coupled runs of the Coulomb law (taub_method = "coulomb"; docs/backend-plastic.md) along the feasible band of the LGM
# calibration (tests/plastic_coulomb_calib.f90): PlasticIceSheet online on VILMA2, 80 ka to 0, as D1 (tau-PMIST-PD |
# E1D-PMIST: the same par file, margin chronology, grids and output) but with the Coulomb law. Four sets:
#
#   experiments/coulomb-band.sh setup <set>                  run folder, c_f files and par files of the set
#   experiments/coulomb-band.sh submit <set> [<name> ...]    runme jobs (default all of the set)
#
#   band   tau-CB-H<H_N> | <Earth> (C1, C2; W default $KRYOS_WORK/coulomb-taub/coupled): the region maps v1 (Eurasia one region),
#          H_N 30, 100, 300, 1000 m with lambda_min 0.2, 0.3, 0.5, 0.75 (the cell nearest the band's centre line, an Antarctic
#          20 ka excess near 12.5 m) and c_f per ice sheet of the calibration lgm-calib (NAIS, EIS to 72.5, 19.5 m at 20 ka,
#          Greenland and Antarctica from the present-day fit; Patagonia: Eurasia's, it has none), on the Earths E1D-PMIST and
#          E1D-ES1-E117; out_2D every 2.5 kyr
#   split  tau-CBS-H<H_N> | E1D-ES1-E117 (C3; W default $KRYOS_WORK/coulomb-taub/split): the region maps v2 (Eurasia split
#          into Fennoscandia and BarentsKara, every cell with a region) and c_f of the calibration lgm-calib-v2 (NAIS 72.5,
#          Fennoscandia 11.20, BarentsKara 6.57 m: PaleoMIST's at 20 ka); H_N 1000 m, lambda_min 0.75 (the band's H1000, the
#          A/B of the split), H_N 2000 and 3000 m with lambda_min 1.0 (the centre line: the excess 12.5 m at 2000 m; at 3000 m
#          none reaches it, 1.0 is the largest); Patagonia keeps the c_f lgm-calib gave Eurasia at the point; out_2D every
#          2.5 kyr to 30 ka, then every 1 kyr, with the 2.5 kyr ages of the scores
#   delta  tau-CBD-H<H_N> | E1D-ES1-E117 (C5; W default $KRYOS_WORK/coulomb-taub/delta): split, with the floor of N,
#          taub_delta 0.05, and c_f of the calibration lgm-calib-v3 (lgm-calib-v2 with the floor); H_N 1000 m, lambda_min
#          0.75 (C3's H1000) and H_N 300 m, lambda_min 0.5 (C2's H300); Patagonia as split; out_2D every 2.5 kyr, every
#          500 yr from 30 ka
#   centre tau-CBDC-H<H_N> | E1D-ES1-E117 (C5b; W default $KRYOS_WORK/coulomb-taub/delta): delta on the centre line of the
#          band of lgm-calib-v3 (the Antarctic 20 ka excess near 12.5 m; lgm-calib-v3/fine: lambda_min in steps of 0.02):
#          H_N 1000 m with lambda_min 0.34, H_N 300 m with 0.24; Patagonia takes Fennoscandia's c_f (lgm-calib has no
#          Eurasia at these points)
#   file   tau-CBF-H<H_N> | E1D-ES1-E117 (C7; W default $KRYOS_WORK/coulomb-taub/lambda-file): C6 (tau-CBDAR: delta with
#          Antarctica of the calibration coulomb-taub/ant-calib) with lambda on the bed of the file in every instance
#          (taub_lambda_bed "file": {grid}_lambda_bed.nc, the calibration's bed, whatever the Earth) and c_f recalibrated
#          on it (lgm-calib-v4: the north as lgm-calib-v3; lgm-calib-v4/ant: Antarctica by ant-calib's pick rule); ice_nor
#          and ice_pat lambda_min 0.75 at H_N 1000 m and 0.5 at 300 m, zb -300 / 200 m (C5); ice_ant lambda_min 0.34 and
#          0.22, zb -600 / 0 m (ant-calib); Patagonia as delta
#   band and split were run before the floor (N = 0 at flotation); with this kryos they take the default taub_delta 0.05.
#
#   all    taub_p 1, taub_smooth_km 0 (as the calibration); taub_lambda_bed "relaxed" (the run's reference topography: on a
#          refcorr topography lambda depends on the Earth) except file. Earths: E1D-PMIST (Gowan's 1-D, r of
#          PaleoMIST-1D as D1), E1D-ES1-E117 (wave A's best 1-D under PaleoMIST ice; r of ice-PMIST | E1D-ES1-E117, its
#          reference files written for all domains)
#
# Env: BASE (D1's par file), R_PMIST (D1's r case), R_E117 (the E117 r case with reference files for North, Antarctica and
# PatagoniaPaleo: gia_refcorr.jl write of its r.nc), W. Runs in $W/runs/<name>, about 10 min (band) to 25 min (split, delta)
# each on 64 threads (D1: 7.4 min).
set -euo pipefail
X=$(cd "$(dirname "$(readlink -f "$0")")" && pwd); T=$X/../tools/levante; . "$T/env.sh"
[ $# -ge 2 ] || { sed -n '2,43p' "$0"; exit 1; }
SET=$2
BASE=${BASE:-/albedo/work/user/alrobi001/kryos-work/runs/pm-online/par_pdfac.nml}
R_PMIST=${R_PMIST:-/albedo/work/user/alrobi001/kryos-work/kryos_data/refcorr/PaleoMIST-1D}
case $SET in
    band)
        W=${W:-$KRYOS_WORK/coulomb-taub/coupled}; PREFIX=tau-CB; REGIONS=taub_regions; EARTHS=(E1D-PMIST E1D-ES1-E117)
        OUT2D="'const'|2500.0"; WALL=02:00:00; DELTA=; LAMBDA_BED=relaxed
        # H_N  lambda_min  c_f: Laurentide, Greenland, Eurasia, Antarctica, Patagonia (lgm-calib/calib_*.nc)
        POINTS=("30 0.2 0.187851 0.300839 0.306872 0.346899 0.306872"
                "100 0.3 0.056835 0.101503 0.091166 0.115349 0.091166"
                "300 0.5 0.020991 0.044263 0.033608 0.043862 0.033608"
                "1000 0.75 0.009206 0.022060 0.015105 0.019586 0.015105") ;;
    split)
        W=${W:-$KRYOS_WORK/coulomb-taub/split}; PREFIX=tau-CBS; REGIONS=taub_regions_v2; EARTHS=(E1D-ES1-E117)
        OUT2D="'times'|'-80000:2500:0, -30000:1000:0'"; WALL=01:00:00; DELTA=; LAMBDA_BED=relaxed
        # H_N  lambda_min  c_f: Laurentide, Greenland, Fennoscandia, Antarctica (lgm-calib-v2/calib_*.nc), Patagonia (lgm-calib's
        # Eurasia), BarentsKara (lgm-calib-v2)
        POINTS=("1000 0.75 0.009206 0.022060 0.010553 0.019586 0.015105 0.021621"
                "2000 1.0 0.006423 0.016654 0.007807 0.013282 0.010480 0.013762"
                "3000 1.0 0.005746 0.014740 0.007077 0.011415 0.009470 0.012578") ;;
    delta)
        W=${W:-$KRYOS_WORK/coulomb-taub/delta}; PREFIX=tau-CBD; REGIONS=taub_regions_v2; EARTHS=(E1D-ES1-E117)
        OUT2D="'times'|'-80000:2500:0, -30000:500:0'"; WALL=01:30:00; DELTA=0.05; LAMBDA_BED=relaxed
        # H_N  lambda_min  c_f: Laurentide, Greenland, Fennoscandia, Antarctica (lgm-calib-v3/calib_*.nc), Patagonia (lgm-calib's
        # Eurasia), BarentsKara (lgm-calib-v3)
        POINTS=("1000 0.75 0.008536 0.022060 0.010375 0.019212 0.015105 0.019428"
                "300 0.5 0.019699 0.044263 0.021662 0.043037 0.033608 0.047985") ;;
    centre)
        W=${W:-$KRYOS_WORK/coulomb-taub/delta}; PREFIX=tau-CBDC; REGIONS=taub_regions_v2; EARTHS=(E1D-ES1-E117)
        OUT2D="'times'|'-80000:2500:0, -30000:500:0'"; WALL=01:30:00; DELTA=0.05; LAMBDA_BED=relaxed
        # H_N  lambda_min  c_f: Laurentide, Greenland, Fennoscandia, Antarctica (lgm-calib-v3/fine/calib_*.nc), Patagonia
        # (Fennoscandia's), BarentsKara (lgm-calib-v3/fine)
        POINTS=("1000 0.34 0.011105 0.022204 0.012104 0.021030 0.012104 0.035994"
                "300 0.24 0.024691 0.044464 0.024552 0.045971 0.024552 0.082073") ;;
    file)
        W=${W:-$KRYOS_WORK/coulomb-taub/lambda-file}; PREFIX=tau-CBF; REGIONS=taub_regions_v2; EARTHS=(E1D-ES1-E117)
        OUT2D="'times'|'-80000:2500:0, -30000:500:0'"; WALL=01:30:00; DELTA=0.05; LAMBDA_BED=file
        # H_N  lambda_min (ice_nor, ice_pat)  c_f: Laurentide, Greenland, Fennoscandia (lgm-calib-v4/calib_NH-32KM.nc),
        # Antarctica (lgm-calib-v4/ant), Patagonia (C5's: lgm-calib's Eurasia), BarentsKara (lgm-calib-v4)  ice_ant:
        # lambda_min, zb0, zb1
        POINTS=("1000 0.75 0.008536 0.022060 0.010375 0.019432 0.015105 0.019428 0.34 -600.0 0.0"
                "300 0.5 0.019699 0.044263 0.021662 0.042123 0.033608 0.047985 0.22 -600.0 0.0") ;;
    *) echo "unknown set $SET (band, split, delta, centre, file)"; exit 1 ;;
esac
R=$W/runs
R_E117=${R_E117:-$KRYOS_WORK/coulomb-taub/coupled/refcorr/ice-PMIST_E1D-ES1-E117}
name() { printf "%s-H%03d_%s" "$PREFIX" "$1" "$2"; }      # folder of "<PREFIX>-H<H_N> | <Earth>"

setup() {
    mkdir -p "$R" "$W/cf"; cd "$R"
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
    mkdir -p maps; cp -n "$KRYOS_DATA"/maps/grid_*.txt "$KRYOS_DATA"/maps/scrip-*.nc maps/ 2> /dev/null || true
    for p in "${POINTS[@]}"; do
        read -r hn lam cl cg ce ca cp cb la za0 za1 <<< "$p"         # la, za0, za1: ice_ant's lambda_min, zb0, zb1
        la=${la:-$lam}; za0=${za0:--300.0}; za1=${za1:-200.0}
        cf=$W/cf/cf-${PREFIX#tau-}-H$(printf %03d $hn).nc       # per set: sets may share W
        if [ "$SET" = band ]; then           # v1 maps: Eurasia is region 3; region 6 (none on v1) takes Eurasia's c_f
            kryos_jl "$T/plastic_cf_regions.jl" cf "$cf" "$cl" "$cg" "$ce" "$ca" "$cp" "$ce"
        else
            kryos_jl "$T/plastic_cf_regions.jl" cf "$cf" "$cl" "$cg" "$ce" "$ca" "$cp" "$cb"
        fi
        for e in "${EARTHS[@]}"; do
            n=$(name "$hn" "$e"); rc=$R_PMIST; [ "$e" = E1D-ES1-E117 ] && rc=$R_E117
            python3 - "$BASE" "par_$n.nml" "$hn" "$lam" "$cf" "$rc" "$e" "$REGIONS" "$OUT2D" "$DELTA" "$la" "$za0" "$za1" "$LAMBDA_BED" << 'EOF'
import re, sys
base, out, hn, lam, cf, rc, earth, regions, out2d, delta, lam_ant, zb0_ant, zb1_ant, lambda_bed = sys.argv[1:]
lines = open(base).read().split("\n")
def group(g):
    i = next(k for k, l in enumerate(lines) if l.strip().lower() == "&" + g)
    j = next(k for k in range(i + 1, len(lines)) if lines[k].strip() == "/")
    return i, j
def put(g, key, val):
    i, j = group(g)
    for k in range(i + 1, j):
        if re.match(r"\s*" + key + r"\s*=", lines[k], re.I):
            lines[k] = " %-15s = %s" % (key, val); return
    lines.insert(j, " %-15s = %s" % (key, val))
for g in ("ice_nor", "ice_ant", "ice_pat"):
    for k, v in [("taub_method", "'coulomb'"), ("taub_domain_factor_file", "''"), ("taub_scale_file", "''"),
                 ("taub_region_file", "'kryos_data/plastic/{grid_name}_%s.nc'" % regions), ("taub_cf_file", "'%s'" % cf),
                 ("taub_HN", hn + ".0"), ("taub_p", "1.0"), ("taub_lam_min", lam_ant if g == "ice_ant" else lam),
                 ("taub_zb0", zb0_ant if g == "ice_ant" else "-300.0"),
                 ("taub_zb1", zb1_ant if g == "ice_ant" else "200.0"), ("taub_lambda_bed", "'%s'" % lambda_bed),
                 ("taub_smooth_km", "0.0"), ("out_plastic", "'out_2D'")]:
        put(g, k, v)
    if delta:
        put(g, "taub_delta", delta)
def drop(g, key):
    i, j = group(g)
    for k in range(i + 1, j):
        if re.match(r"\s*" + key + r"\s*=", lines[k], re.I):
            del lines[k]; return
method, when = out2d.split("|")
put("out_2d", "method", method)
if method == "'const'":
    drop("out_2d", "times"); put("out_2d", "dt", when)
else:
    drop("out_2d", "dt"); put("out_2d", "times", when)
topo = {"topo_nor": "North/{grid_name}/{grid_name}_TOPO-GEBCO2026-refcorr.nc",
        "topo_ant": "Antarctica/{grid_name}/{grid_name}_TOPO-BedMachine-v4-refcorr.nc",
        "topo_pat": "PatagoniaPaleo/{grid_name}/{grid_name}_TOPO-GEBCO2026-refcorr.nc"}
for g, f in topo.items():
    put(g, "path", "'%s/%s'" % (rc, f))
for g in ("x_z_bed", "x_h_ice"):
    put(g, "filename", "'%s/{grid_name}/{grid_name}_REFCORR.nc'" % rc)
if earth == "E1D-ES1-E117":     # wave A's E117 (earths/E117.nml of the Earth search)
    for k, v in [("n_layer", "15"),
                 ("r_top", "6371.000e3, 6291.000e3, 6260.110e3, 6251.000e3, 6151.000e3, 5971.000e3, 5961.000e3, 5771.000e3, "
                           "5701.000e3, 5600.000e3, 4943.000e3, 4871.000e3, 4287.000e3, 3630.000e3, 3480.000e3"),
                 ("r_bot", "6291.000e3, 6260.110e3, 6251.000e3, 6151.000e3, 5971.000e3, 5961.000e3, 5771.000e3, 5701.000e3, "
                           "5600.000e3, 4943.000e3, 4871.000e3, 4287.000e3, 3630.000e3, 3480.000e3, 0.000e3"),
                 ("eta", "1.0000e+40, 1.0000e+40, 1.3521e+20, 1.3521e+20, 1.3521e+20, 1.3521e+20, 5.5847e+21, 5.5847e+21, "
                         "1.4093e+23, 1.4093e+23, 1.4093e+23, 4.8641e+21, 4.8641e+21, 4.8641e+21, 0.0"),
                 ("rheology", "0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2")]:
        put("vilma", k, v)
open(out, "w").write("\n".join(lines))
EOF
            echo "par_$n.nml"
        done
    done
}

submit() {
    cd "$R"
    local names=("$@")
    if [ ${#names[@]} = 0 ]; then
        for p in "${POINTS[@]}"; do for e in "${EARTHS[@]}"; do names+=("$(name ${p%% *} $e)"); done; done
    fi
    for n in "${names[@]}"; do
        [ -f "par_$n.nml" ] || { echo "no par_$n.nml: run setup"; exit 1; }
        runme -r -s -e kryos -n "par_$n.nml" -o "$n" -q "$KRYOS_QUEUE" -w "$WALL" --omp 64 -m 30G --jobname "${n#tau-}" 2>&1 | tail -1
    done
}

case $1 in
    setup)  setup ;;
    submit) shift 2; submit "$@" ;;
    *)      echo "unknown command $1 (setup, submit)"; exit 1 ;;
esac
