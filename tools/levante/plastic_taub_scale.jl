# The scale groups of Gowan's tau_b domains and the multiplier files of the plastic ice backend (taub_scale_file,
# docs/dev/plastic-backend.md), from the [scale] table of a calibration config (plastic_esmda.toml):
#   julia --project=tools tools/levante/plastic_taub_scale.jl <config.toml> <out.nc> [<log factor>]
# writes a multiplier file with that log factor for every group and node (default 0: Gowan's tau_b), and <out>.groups.csv.
# Included by plastic_esmda.jl (scale_groups, write_scale).
#
# Groups: [[scale.groups]] in order, each with `name`, `region` (Gowan's region of the domain: 1 North America incl. Greenland,
# 2 Eurasia, 3 Antarctica) and optionally `greenland` (true/false): a domain is in Greenland when more than half of its cells on
# the grid lie in the Greenland region of the REGIONS map (`regions_var` = `greenland_code`, FesmData Regions v2: region_2 = 103,
# a partition of land and shelf by territory). A domain matching no group is unscaled (group 0); one matching two stops.
# Nodes: `nodes_ka` [ka BP], any order.
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
    f_greenland::Vector{Float64}
end
ngroups(sg::ScaleGroups) = length(sg.names)
nnodes(sg::ScaleGroups) = length(sg.nodes_ka)

function scale_groups(sc::AbstractDict)
    ids, region, dom = NCDataset(sc["taub_file"]) do ds
        Int.(ds["domain_id"][:]), Int.(ds["region"][:]), Int.(nomissing(ds["taub_domain"][:, :], -9999))
    end
    reg = NCDataset(ds -> Int.(nomissing(ds[get(sc, "regions_var", "region_2")][:, :], -9999)), sc["regions_file"])
    size(reg) == size(dom) || error("the REGIONS map is not on the grid of the tau_b product: $(size(reg)), $(size(dom))")
    gcode = get(sc, "greenland_code", 103)
    n = [count(==(d), dom) for d in ids]
    fg = [n[i] == 0 ? 0.0 : count((dom .== d) .& (reg .== gcode)) / n[i] for (i, d) in enumerate(ids)]
    specs = sc["groups"]
    names = [String(g["name"]) for g in specs]
    grp = zeros(Int, length(ids))
    for i in eachindex(ids)
        hit = [g for (g, s) in enumerate(specs) if s["region"] == region[i] && (!haskey(s, "greenland") || s["greenland"] == (fg[i] > 0.5))]
        length(hit) > 1 && error("domain $(ids[i]) matches the groups $(names[hit])")
        grp[i] = isempty(hit) ? 0 : hit[1]
    end
    ScaleGroups(names, Float64.(sc["nodes_ka"]), ids, region, grp, n, fg)
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
        println(io, "domain_id,region,ncells,f_greenland,group,group_name")
        for i in eachindex(sg.domain_id)
            g = sg.group[i]
            println(io, join((sg.domain_id[i], sg.region[i], sg.ncells[i], @sprintf("%.3f", sg.f_greenland[i]), g,
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
    mixed = [i for i in eachindex(sg.domain_id) if 0 < sg.f_greenland[i] < 1]
    push!(lines, "domains partly in Greenland (id: fraction -> group): " *
                 join([@sprintf("%d: %.2f -> %s", sg.domain_id[i], sg.f_greenland[i], sg.group[i] == 0 ? "unscaled" : sg.names[sg.group[i]]) for i in mixed], ", "))
    lines
end

if abspath(PROGRAM_FILE) == @__FILE__
    sg = scale_groups(TOML.parsefile(ARGS[1])["scale"])
    v = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.0
    foreach(println, groups_summary(sg))
    println("wrote ", write_scale(ARGS[2], sg, fill(v, ngroups(sg), nnodes(sg))), ", ", write_groups_csv(ARGS[2] * ".groups.csv", sg))
end
