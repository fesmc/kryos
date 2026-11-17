# The region maps of the Coulomb yield law of the plastic ice backend (taub_method = "coulomb", taub_region_file;
# docs/backend-plastic.md), one per grid, and c_f files (taub_cf_file):
#   julia --project=tools/levante tools/levante/plastic_cf_regions.jl regions [<out dir>]
#       writes <out dir>/<grid>_taub_regions.nc for NH-32KM, ANT-16KM, PAT-PAL-8KM (default $KRYOS_DATA/plastic)
#   julia --project=tools/levante tools/levante/plastic_cf_regions.jl cf <out.nc> <c_f> [<c_f> x 5]
#       writes a c_f file: one c_f for every region, or one per region in the order of the ids
#
# Regions (the ice sheets; ids 1-5): 1 Laurentide (North America with the Cordillera and Innuitian), 2 Greenland,
# 3 Eurasia, 4 Antarctica, 5 Patagonia. A cell takes the region of its Gowan tau_b domain in the product
# {grid}_IceHistory-PaleoMIST_taub.nc (domain_id = 1000 region + Gowan's number; region 1 North America with
# Greenland, 2 Eurasia, 3 Antarctica, 4 Patagonia). North America is split as the scale groups of
# plastic_taub_scale.jl: a domain is Greenland when more than half of its cells on the grid lie in REGIONS
# region_2 code 103; the "outside every polygon" domains (Gowan's number 0) are split cell by cell. Cells without a
# domain have region 0 (none). Only the geography of the product is used: none of its tau_b values.
#
# Region file: region(xc, yc) (Int32, 0 = none; flag_values, flag_meanings), region_id(region), attribute region_names.
# c_f file: region_id(region), log_cf(region) (natural log of c_f).
using NCDatasets

const REGION_NAMES = ["Laurentide", "Greenland", "Eurasia", "Antarctica", "Patagonia"]
const GRIDS = [("North", "NH-32KM"), ("Antarctica", "ANT-16KM"), ("PatagoniaPaleo", "PAT-PAL-8KM")]
const GREENLAND_CODE = 103

function region_map(taub_file, regions_file)
    dom, xc, yc = NCDataset(ds -> (Int.(nomissing(ds["taub_domain"][:, :], -9999)), ds["xc"][:], ds["yc"][:]), taub_file)
    r2 = NCDataset(ds -> Int.(nomissing(ds["region_2"][:, :], -9999)), regions_file)
    size(r2) == size(dom) || error("$regions_file is not on the grid of $taub_file")
    ngrl = Dict{Int,Int}(); n = Dict{Int,Int}()
    for c in eachindex(dom)
        dom[c] > 0 || continue
        n[dom[c]] = get(n, dom[c], 0) + 1
        ngrl[dom[c]] = get(ngrl, dom[c], 0) + (r2[c] == GREENLAND_CODE)
    end
    reg = zeros(Int32, size(dom))
    for c in eachindex(dom)
        d = dom[c]; d > 0 || continue
        g = div(d, 1000)
        reg[c] = if g == 1
            grl = mod(d, 1000) == 0 ? r2[c] == GREENLAND_CODE : ngrl[d] > n[d] / 2
            grl ? 2 : 1
        elseif g in 2:4
            g + 1
        else
            error("$taub_file: domain $d of no known region")
        end
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
            f = write_regions(joinpath(out, "$(g)_taub_regions.nc"), reg, xc, yc, basename(tf) * ", " * basename(rf))
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
