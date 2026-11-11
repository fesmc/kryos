# The Earth parameters of a weighted ensemble (plastic_ensemble.jl, [earth] of its config): shifts of log10(eta) of VILMA2's 3-D
# viscosity field (VILMA2 visc_shift_file), a global one in the upper mantle (surface to boundary_km) and one in the lower mantle
# (boundary_km to the core), and one per region in the upper mantle, weighted by the region's mask:
#   dlog10 eta(lon, lat, depth) = UM + sum_k w_k(lon, lat) R_k   above boundary_km,   LM below.
#   julia --project=tools/levante tools/levante/plastic_ensemble.jl earth-regions <config.toml> <out.nc>    the region masks
#   julia --project=tools/levante tools/levante/plastic_ensemble.jl earth-shift <config.toml> <out.nc> UM=0.5 LM=-0.5 <region>=1.0 ...
# Region masks (`regions`): [[earth.regions]] (name, groups: scale groups of [scale], plastic_taub_scale.jl), the core of a region
# the cells of the tau_b domains of its groups on the grids of [scale] taub_files (the ice of PaleoMIST), gathered on a global
# 1-degree grid; t_k(d) = 1 within buffer_km of the core, cos^2 to 0 over taper_km beyond (d the great-circle distance to the core);
# the masks w_k = t_k / max(1, sum_j t_j) (continuous, at most 1 in all; overlapping regions share). File: mask(lon, lat, region),
# attribute region_names (";"-separated). Shift files (VILMA2's generic form): dlog10_visc(lon, lat, r) [dex], r the radius [m]
# at the depths 0, boundary_km - 1, boundary_km + 1 and 2891 km (linear in r between them, so the step at the boundary is sharp).
# Included by plastic_esmda.jl (after plastic_taub_scale.jl).
using NCDatasets, Printf, TOML

const R_VILMA_KM = 6371.0          # VILMA2's R_PREM
const CMB_KM = 2891.0              # depth of the core-mantle boundary

# [earth] of a config: the parameter names (UM, LM, the regions), the key, the prior sd, theta_max, the boundary
function earth_spec(cfg)
    e = cfg["earth"]
    regs = [String(r["name"]) for r in get(e, "regions", Any[])]
    (key = String(get(e, "key", "vilma.visc_shift_file")), regions = regs, names = earth_names(regs),
     sigma = Float64(e["sigma"]), theta_max = Float64(get(e, "theta_max", 3 * e["sigma"])),
     boundary_km = Float64(get(e, "boundary_km", 670.0)), masks_file = String(e["regions_file"]))
end
earth_names(regs) = ["earth_UM"; "earth_LM"; ["earth_" * replace(r, r"[^A-Za-z0-9]+" => "_") for r in regs]]

function read_masks(file, regs)
    NCDataset(file) do ds
        names = split(ds.attrib["region_names"], ";")
        names == regs || error("$file: regions $(names), the config's $(regs)")
        Float64.(ds["lon"][:]), Float64.(ds["lat"][:]), Float64.(ds["mask"][:, :, :])
    end
end

# the shift file of one member: θe = [UM, LM, R_1, ...]
function write_earth(file, E, lon, lat, W, θe)
    length(θe) == 2 + length(E.regions) || error("earth: $(length(θe)) values for $(E.names)")
    um = fill(θe[1], length(lon), length(lat))
    for k in eachindex(E.regions); um .+= θe[2+k] .* W[:, :, k]; end
    depth = [CMB_KM, E.boundary_km + 1.0, E.boundary_km - 1.0, 0.0]          # radius ascending
    D = zeros(length(lon), length(lat), 4)
    D[:, :, 1] .= θe[2]; D[:, :, 2] .= θe[2]; D[:, :, 3] .= um; D[:, :, 4] .= um
    isfile(file) && rm(file)
    NCDataset(file, "c") do ds
        defVar(ds, "lon", lon, ("lon",); attrib = ["units" => "degrees_east"])
        defVar(ds, "lat", lat, ("lat",); attrib = ["units" => "degrees_north"])
        defVar(ds, "r", 1000.0 .* (R_VILMA_KM .- depth), ("r",); attrib = ["units" => "m", "long_name" => "radius"])
        defVar(ds, "dlog10_visc", D, ("lon", "lat", "r");
               attrib = ["units" => "dex", "long_name" => "shift of log10 of the viscosity (VILMA2 visc_shift_file)"])
        ds.attrib["parameters"] = join([@sprintf("%s=%.4f", n, x) for (n, x) in zip(E.names, θe)], " ")
    end
    file
end

# --- the region masks -------------------------------------------------------------------------------------------------------------
function hav(lo1, la1, lo2, la2)
    φ1, φ2 = deg2rad(la1), deg2rad(la2)
    2R_VILMA_KM * asin(min(1.0, sqrt(sin((φ2 - φ1) / 2)^2 + cos(φ1) * cos(φ2) * sin(deg2rad(lo2 - lo1) / 2)^2)))
end

function region_masks(cfg, out)
    e = cfg["earth"]; sc = cfg["scale"]; sg = scale_groups(sc)
    regs = e["regions"]; buf = Float64(get(e, "buffer_km", 300.0)); tap = Float64(get(e, "taper_km", 500.0))
    lon = collect(-179.5:1.0:179.5); lat = collect(-89.5:1.0:89.5)
    core = falses(length(lon), length(lat), length(regs))
    gi(x) = clamp(floor(Int, mod(x + 180, 360)) + 1, 1, 360); gj(y) = clamp(floor(Int, y + 90) + 1, 1, 180)
    gpos = Dict(n => g for (g, n) in enumerate(sg.names))
    for (k, r) in enumerate(regs), g in r["groups"]; haskey(gpos, g) || error("region $(r["name"]): no scale group $g"); end
    dgrp = Dict(d => sg.group[i] for (i, d) in enumerate(sg.domain_id))
    for tf in (haskey(sc, "taub_files") ? sc["taub_files"] : [sc["taub_file"]])
        grid = basename(dirname(tf)); gf = joinpath(dirname(tf), "$(grid)_grid.nc")
        dom = NCDataset(ds -> Int.(nomissing(ds["taub_domain"][:, :], -9999)), tf)
        lo, la = NCDataset(ds -> (Float64.(ds["lon2D"][:, :]), Float64.(ds["lat2D"][:, :])), gf)
        size(lo) == size(dom) || error("$gf is not the grid of $tf")
        for c in eachindex(dom)
            g = get(dgrp, dom[c], 0); g == 0 && continue
            for (k, r) in enumerate(regs)
                sg.names[g] in r["groups"] && (core[gi(lo[c]), gj(la[c]), k] = true)
            end
        end
    end
    T = zeros(length(lon), length(lat), length(regs))
    for k in eachindex(regs)
        cc = [(lon[I[1]], lat[I[2]]) for I in findall(core[:, :, k])]
        isempty(cc) && error("region $(regs[k]["name"]): no cells")
        Threads.@threads for j in eachindex(lat)
            for i in eachindex(lon)
                d = minimum(hav(lon[i], lat[j], a, b) for (a, b) in cc)
                T[i, j, k] = d <= buf ? 1.0 : d >= buf + tap ? 0.0 : cos(π / 2 * (d - buf) / tap)^2
            end
        end
        println(regs[k]["name"], ": ", length(cc), " core cells (1 deg)")
    end
    W = T ./ max.(1.0, sum(T; dims = 3))
    isfile(out) && rm(out)
    NCDataset(out, "c") do ds
        defVar(ds, "lon", lon, ("lon",); attrib = ["units" => "degrees_east"])
        defVar(ds, "lat", lat, ("lat",); attrib = ["units" => "degrees_north"])
        defVar(ds, "mask", W, ("lon", "lat", "region"); attrib = ["units" => "1",
               "long_name" => "weight of the region's upper-mantle viscosity shift"])
        ds.attrib["region_names"] = join([String(r["name"]) for r in regs], ";")
        ds.attrib["region_groups"] = join([String(r["name"]) * ": " * join(r["groups"], ", ") for r in regs], "; ")
        ds.attrib["buffer_km"] = buf; ds.attrib["taper_km"] = tap
        ds.attrib["source"] = "plastic_earth_shift.jl regions: the tau_b domains of the scale groups on " *
                              join(basename.(haskey(sc, "taub_files") ? sc["taub_files"] : [sc["taub_file"]]), ", ")
    end
    println("wrote ", out)
end
