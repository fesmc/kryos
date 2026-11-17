# The region maps of the Coulomb yield law of the plastic ice backend (taub_method = "coulomb", taub_region_file;
# docs/backend-plastic.md), one per grid, and c_f files (taub_cf_file):
#   julia --project=tools/levante tools/levante/plastic_cf_regions.jl regions [<out dir>]
#       writes <out dir>/<grid>_taub_regions_v2.nc for NH-32KM, ANT-16KM, PAT-PAL-8KM (default $KRYOS_DATA/plastic)
#   julia --project=tools/levante tools/levante/plastic_cf_regions.jl cf <out.nc> <c_f> [<c_f> x 6]
#       writes a c_f file: one c_f for every region, or one per region in the order of the ids
#
# Regions (ids 1-6): 1 Laurentide (North America with the Cordillera and Innuitian), 2 Greenland, 3 Fennoscandia
# (Fennoscandia with the British-Irish ice and the rest of Eurasia), 4 Antarctica, 5 Patagonia, 6 BarentsKara (the
# Barents and Kara seas with Svalbard, Franz Josef Land and Novaya Zemlya). A cell takes the region of its Gowan tau_b
# domain in the product {grid}_IceHistory-PaleoMIST_taub.nc (domain_id = 1000 region + Gowan's number; region 1 North
# America with Greenland, 2 Eurasia, 3 Antarctica, 4 Patagonia). North America and Eurasia are split as the scale groups
# of plastic_taub_scale.jl, by whole domains: a domain is Greenland when more than half of its cells on the grid lie in
# REGIONS region_2 code 103, BarentsKara when more than half lie in region_3 codes 10201 (Barents-Kara) and 10203
# (Svalbard); the "outside every polygon" domains (Gowan's number 0) are split cell by cell. Cells without a domain have
# region 0 (none). Only the geography of the product is used: none of its tau_b values.
#
# Versions: v1, {grid}_taub_regions.nc (2026-11-17, kryos 8ff43a2: ids 1-5, region 3 all of Eurasia; the runs C1, C2 of
# docs/science/results.md); v2, {grid}_taub_regions_v2.nc (this tool).
#
# Region file: region(xc, yc) (Int32, 0 = none; flag_values, flag_meanings), region_id(region), attribute region_names.
# c_f file: region_id(region), log_cf(region) (natural log of c_f).
using NCDatasets

const REGION_NAMES = ["Laurentide", "Greenland", "Fennoscandia", "Antarctica", "Patagonia", "BarentsKara"]
const VERSION = "v2"
const GRIDS = [("North", "NH-32KM"), ("Antarctica", "ANT-16KM"), ("PatagoniaPaleo", "PAT-PAL-8KM")]
const GREENLAND_CODE = 103               # REGIONS region_2
const BARENTS_CODES = (10201, 10203)     # REGIONS region_3: Barents-Kara, Svalbard

function region_map(taub_file, regions_file)
    dom, xc, yc = NCDataset(ds -> (Int.(nomissing(ds["taub_domain"][:, :], -9999)), ds["xc"][:], ds["yc"][:]), taub_file)
    r2, r3 = NCDataset(ds -> (Int.(nomissing(ds["region_2"][:, :], -9999)), Int.(nomissing(ds["region_3"][:, :], -9999))), regions_file)
    size(r2) == size(dom) || error("$regions_file is not on the grid of $taub_file")
    ingrl(c) = r2[c] == GREENLAND_CODE
    inbk(c) = r3[c] in BARENTS_CODES
    n = Dict{Int,Int}(); ngrl = Dict{Int,Int}(); nbk = Dict{Int,Int}()
    for c in eachindex(dom)
        d = dom[c]; d > 0 || continue
        n[d] = get(n, d, 0) + 1
        ngrl[d] = get(ngrl, d, 0) + ingrl(c)
        nbk[d] = get(nbk, d, 0) + inbk(c)
    end
    # a sub-region of a split region: the domain's majority, cell by cell in the "outside every polygon" domain
    sub(d, c, f, m) = mod(d, 1000) == 0 ? f(c) : m[d] > n[d] / 2
    reg = zeros(Int32, size(dom))
    for c in eachindex(dom)
        d = dom[c]; d > 0 || continue
        g = div(d, 1000)
        reg[c] = g == 1 ? (sub(d, c, ingrl, ngrl) ? 2 : 1) :
                 g == 2 ? (sub(d, c, inbk, nbk) ? 6 : 3) :
                 g in 3:4 ? g + 1 : error("$taub_file: domain $d of no known region")
    end
    reg, xc, yc
end

function write_regions(file, reg, xc, yc, src)
    isfile(file) && rm(file)
    NCDataset(file, "c") do ds
        defDim(ds, "region", length(REGION_NAMES))
        defVar(ds, "xc", xc, ("xc",); attrib = ["units" => "km", "axis" => "X"])
        defVar(ds, "yc", yc, ("yc",); attrib = ["units" => "km", "axis" => "Y"])
        defVar(ds, "region", reg, ("xc", "yc"); attrib = ["long_name" => "ice-sheet region of the Coulomb law's c_f (0 = none)",
               "flag_values" => Int32.(1:length(REGION_NAMES)), "flag_meanings" => join(REGION_NAMES, " ")])
        defVar(ds, "region_id", Int32.(1:length(REGION_NAMES)), ("region",); attrib = ["long_name" => "region id"])
        ds.attrib["region_names"] = join(REGION_NAMES, " ")
        ds.attrib["source"] = "kryos tools/levante/plastic_cf_regions.jl from $src (domain geography only)"
        ds.attrib["version"] = VERSION
    end
    file
end

function write_cf(file, cf::AbstractVector)
    length(cf) == length(REGION_NAMES) || error("one c_f per region ($(length(REGION_NAMES))) is needed")
    all(>(0), cf) || error("c_f must be positive")
    isfile(file) && rm(file)
    NCDataset(file, "c") do ds
        defDim(ds, "region", length(cf))
        defVar(ds, "region_id", Int32.(1:length(cf)), ("region",); attrib = ["long_name" => "region id"])
        defVar(ds, "log_cf", log.(Float64.(cf)), ("region",); attrib = ["units" => "1", "long_name" => "natural log of c_f"])
        ds.attrib["region_names"] = join(REGION_NAMES, " ")
    end
    file
end

if abspath(PROGRAM_FILE) == @__FILE__
    cmd = isempty(ARGS) ? "" : ARGS[1]
    if cmd == "regions"
        fesmdata = get(ENV, "FESMDATA", "fesmdata")
        out = length(ARGS) >= 2 ? ARGS[2] : joinpath(get(ENV, "KRYOS_DATA", "kryos_data"), "plastic")
        for (dom, g) in GRIDS
            tf = joinpath(fesmdata, dom, g, "$(g)_IceHistory-PaleoMIST_taub.nc")
            rf = joinpath(fesmdata, dom, g, "$(g)_REGIONS.nc")
            reg, xc, yc = region_map(tf, rf)
            f = write_regions(joinpath(out, "$(g)_taub_regions_$(VERSION).nc"), reg, xc, yc, basename(tf) * ", " * basename(rf))
            println(f, ": ", join(["$(REGION_NAMES[r]) $(count(==(r), reg))" for r in 1:length(REGION_NAMES) if any(==(r), reg)], ", "),
                    " cells; none ", count(==(0), reg))
        end
    elseif cmd == "cf" && length(ARGS) in (3, 2 + length(REGION_NAMES))
        v = parse.(Float64, ARGS[3:end])
        println("wrote ", write_cf(ARGS[2], length(v) == 1 ? fill(v[1], length(REGION_NAMES)) : v))
    else
        println("usage: plastic_cf_regions.jl regions [<out dir>] | cf <out.nc> <c_f> [<c_f> x $(length(REGION_NAMES))]")
        exit(1)
    end
end
