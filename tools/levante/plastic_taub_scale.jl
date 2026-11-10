# The scale groups of Gowan's tau_b domains and the multiplier files of the plastic ice backend (taub_scale_file,
# docs/dev/plastic-backend.md), from the [scale] table of a calibration config (plastic_esmda.toml):
#   julia --project=tools/levante tools/levante/plastic_taub_scale.jl <config.toml> <out.nc> [<log factor>]
# writes a multiplier file with that log factor for every group and node (default 0: Gowan's tau_b), and <out>.groups.csv.
# Included by plastic_esmda.jl (scale_groups, write_scale).
#
# Grids: `taub_files` (the tau_b products of the grids the calibrated instances run on, all with the same domains; or one
# `taub_file`) and `regions_files` (their REGIONS maps; or `regions_file`): cells are counted over all of them.
# Groups: [[scale.groups]] in order, each with `name`, `region` (Gowan's region of the domain: 1 North America incl. Greenland,
# 2 Eurasia, 3 Antarctica, 4 Patagonia) and optionally a sub-region: `codes` (REGIONS codes of the variable `codes_var`, default
# "region_3") and `codes_in` (default true): the domain belongs to the group when more than half of its cells on the grids lie in
# those codes (codes_in = false: when not more than half do). FesmData Regions v2, a partition of land and shelf by territory:
# region_2 103 Greenland; region_3 10201 Russian Arctic, 10202 British Isles, 10203 Svalbard, 10204 Fennoscandia, 20101 West
# Antarctica, 20102 East Antarctica, 20103 Antarctic Peninsula. `greenland` (true/false) is short for codes = [greenland_code]
# (default 103) of `regions_var` (default "region_2") with codes_in = greenland. A domain matching no group is unscaled (group 0);
# one matching two stops. A group may set its own prior correlation time `corr_time_kyr` (plastic_esmda.jl prior_cov).
# Nodes: `nodes_ka` [ka BP], any order. `pinned_ka` (plastic_esmda.jl): nodes held at log factor 0 in every group.
#
# File: dimensions domain, scale_group, node; domain_id(domain), group(domain) (0 = unscaled), time(node) ["years ago"],
# log_factor(scale_group, node) (natural log of the multiplier); attribute group_names.
using NCDatasets, Printf, TOML

struct ScaleGroups
    names::Vector{String}
    nodes_ka::Vector{Float64}
    domain_id::Vector{Int}
    region::Vector{Int}
    group::Vector{Int}
    ncells::Vector{Int}
    f_sub::Vector{Float64}     # the fraction of the domain's cells in the sub-region of its group (1 without one)
end
ngroups(sg::ScaleGroups) = length(sg.names)
nnodes(sg::ScaleGroups) = length(sg.nodes_ka)

# the sub-region of a group spec: (variable, codes, inside) or nothing
function subregion(s, sc)
    haskey(s, "greenland") && return (get(sc, "regions_var", "region_2"), [Int(get(sc, "greenland_code", 103))], Bool(s["greenland"]))
    haskey(s, "codes") && return (String(get(s, "codes_var", "region_3")), Int.(s["codes"]), Bool(get(s, "codes_in", true)))
    nothing
end

function scale_groups(sc::AbstractDict)
    tfiles = haskey(sc, "taub_files") ? String.(sc["taub_files"]) : [String(sc["taub_file"])]
    rfiles = haskey(sc, "regions_files") ? String.(sc["regions_files"]) : [String(sc["regions_file"])]
    length(tfiles) == length(rfiles) || error("one REGIONS map per tau_b product is needed: $(tfiles), $(rfiles)")
    specs = sc["groups"]; names = [String(g["name"]) for g in specs]
    subs = [subregion(s, sc) for s in specs]
    ids, region = NCDataset(ds -> (Int.(ds["domain_id"][:]), Int.(ds["region"][:])), tfiles[1])
    n = zeros(Int, length(ids)); nin = zeros(Int, length(ids), length(specs))   # cells; cells in each group's sub-region
    pos = Dict(d => i for (i, d) in enumerate(ids))
    for (tf, rf) in zip(tfiles, rfiles)
        i2, dom = NCDataset(ds -> (Int.(ds["domain_id"][:]), Int.(nomissing(ds["taub_domain"][:, :], -9999))), tf)
        i2 == ids || error("$tf: the domain ids differ from those of $(tfiles[1])")
        maps = NCDataset(rf) do ds
            Dict(v => Int.(nomissing(ds[v][:, :], -9999)) for v in unique(first(x) for x in subs if x !== nothing))
        end
        for (v, m) in maps
            size(m) == size(dom) || error("the REGIONS map is not on the grid of the tau_b product: $rf $(size(m)), $tf $(size(dom))")
        end
        for c in eachindex(dom)
            i = get(pos, dom[c], 0); i == 0 && continue
            n[i] += 1
            for (g, x) in enumerate(subs)
                x === nothing && continue
                nin[i, g] += maps[x[1]][c] in x[2]
            end
        end
    end
    frac(i, g) = n[i] == 0 ? 0.0 : nin[i, g] / n[i]
    grp = zeros(Int, length(ids)); fsub = ones(length(ids))
    for i in eachindex(ids)
        hit = [g for (g, s) in enumerate(specs) if s["region"] == region[i] &&
               (subs[g] === nothing || subs[g][3] == (frac(i, g) > 0.5))]
        length(hit) > 1 && error("domain $(ids[i]) matches the groups $(names[hit])")
        grp[i] = isempty(hit) ? 0 : hit[1]
        grp[i] > 0 && subs[grp[i]] !== nothing && (fsub[i] = frac(i, grp[i]))
    end
    ScaleGroups(names, Float64.(sc["nodes_ka"]), ids, region, grp, n, fsub)
end

# log_factor (scale_group, node)
function write_scale(file, sg::ScaleGroups, lf::AbstractMatrix)
    size(lf) == (ngroups(sg), nnodes(sg)) || error("log_factor is $(size(lf)), not (group, node) = $((ngroups(sg), nnodes(sg)))")
    isfile(file) && rm(file)
    NCDataset(file, "c") do ds
        defDim(ds, "domain", length(sg.domain_id)); defDim(ds, "scale_group", ngroups(sg)); defDim(ds, "node", nnodes(sg))
        defVar(ds, "domain_id", Int32.(sg.domain_id), ("domain",); attrib = ["long_name" => "domain id of the tau_b product"])
        defVar(ds, "group", Int32.(sg.group), ("domain",);
               attrib = ["long_name" => "scale group of the domain (0 = unscaled)", "flag_meanings" => join(sg.names, " ")])
        defVar(ds, "time", 1000 .* sg.nodes_ka, ("node",); attrib = ["units" => "years ago", "long_name" => "time of the node"])
        defVar(ds, "log_factor", Float64.(lf), ("scale_group", "node");
               attrib = ["units" => "1", "long_name" => "natural log of the tau_b multiplier, linear in time between the nodes"])
        ds.attrib["group_names"] = join(sg.names, " ")
    end
    file
end

function write_groups_csv(file, sg::ScaleGroups)
    open(file, "w") do io
        println(io, "domain_id,region,ncells,f_subregion,group,group_name")
        for i in eachindex(sg.domain_id)
            g = sg.group[i]
            println(io, join((sg.domain_id[i], sg.region[i], sg.ncells[i], @sprintf("%.3f", sg.f_sub[i]), g,
                              g == 0 ? "unscaled" : sg.names[g]), ","))
        end
    end
    file
end

# one line per group: domains with cells on the grid and their cells; the domains split between Greenland and the rest
function groups_summary(sg::ScaleGroups)
    lines = String[]
    for g in 0:ngroups(sg)
        k = findall(i -> sg.group[i] == g && sg.ncells[i] > 0, eachindex(sg.domain_id))
        push!(lines, @sprintf("%-12s %4d domains on the grid, %6d cells", g == 0 ? "unscaled" : sg.names[g], length(k), sum(sg.ncells[k]; init = 0)))
    end
    mixed = [i for i in eachindex(sg.domain_id) if sg.ncells[i] > 0 && 0.2 < sg.f_sub[i] < 0.8]
    push!(lines, "domains split across their group's sub-region, 20-80 % inside (id: fraction inside -> group): " *
                 join([@sprintf("%d: %.2f -> %s", sg.domain_id[i], sg.f_sub[i], sg.group[i] == 0 ? "unscaled" : sg.names[sg.group[i]]) for i in mixed], ", "))
    lines
end

if abspath(PROGRAM_FILE) == @__FILE__
    sg = scale_groups(TOML.parsefile(ARGS[1])["scale"])
    v = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.0
    foreach(println, groups_summary(sg))
    println("wrote ", write_scale(ARGS[2], sg, fill(v, ngroups(sg), nnodes(sg))), ", ", write_groups_csv(ARGS[2] * ".groups.csv", sg))
end
