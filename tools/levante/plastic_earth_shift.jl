# The Earth parameters of a weighted ensemble (plastic_ensemble.jl, [earth] of its config): shifts of log10(eta) of VILMA2's 3-D
# viscosity field (VILMA2 visc_shift_file), a global one in the upper mantle (surface to boundary_km) and one in the lower mantle
# (boundary_km to the core), and one per region in the upper mantle, weighted by the region's mask:
#   dlog10 eta(lon, lat, depth) = UM + sum_k w_k(lon, lat) R_k   above boundary_km,   LM below.
#   julia --project=tools/levante tools/levante/plastic_ensemble.jl earth-regions <config.toml> <out.nc>    the region masks
#   julia --project=tools/levante tools/levante/plastic_ensemble.jl earth-shift <config.toml> <out.nc> UM=0.5 LM=-0.5 <region>=1.0 ...
# Region masks (`regions`): [[earth.regions]] (name, groups: scale groups of [scale], plastic_taub_scale.jl), the core of a region
# the cells of the tau_b domains of its groups on the grids of [scale] taub_files (the ice of PaleoMIST), gathered on a global
# 1-degree grid; t_k(d) = 1 within buffer_km of the core, cos^2 to 0 over taper_km beyond (d the great-circle distance to the core);
# the masks w_k = t_k / max(1, sum_j t_j) (continuous, at most 1 in all; overlapping regions share). [earth] exclude_codes (REGIONS
# codes of exclude_codes_var, default region_2, of [scale] regions_files; Iceland: 105): no regional shift there, all masks times
# 1 - x(d), x = 1 in those areas and cos^2 to 0 over exclude_taper_km (200) around them. File: mask(lon, lat, region),
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
    tfiles = haskey(sc, "taub_files") ? sc["taub_files"] : [sc["taub_file"]]
    rfiles = haskey(sc, "regions_files") ? sc["regions_files"] : [sc["regions_file"]]
    xvar = String(get(e, "exclude_codes_var", "region_2")); xcodes = Int.(get(e, "exclude_codes", Int[]))
    xtap = Float64(get(e, "exclude_taper_km", 200.0))
    excl = falses(length(lon), length(lat))           # the areas of exclude_codes (all REGIONS cells there, not only ice)
    for (tf, rf) in zip(tfiles, rfiles)
        grid = basename(dirname(tf)); gf = joinpath(dirname(tf), "$(grid)_grid.nc")
        dom = NCDataset(ds -> Int.(nomissing(ds["taub_domain"][:, :], -9999)), tf)
        lo, la = NCDataset(ds -> (Float64.(ds["lon2D"][:, :]), Float64.(ds["lat2D"][:, :])), gf)
        size(lo) == size(dom) || error("$gf is not the grid of $tf")
        codes = isempty(xcodes) ? nothing : NCDataset(ds -> Int.(nomissing(ds[xvar][:, :], -9999)), rf)
        for c in eachindex(dom)
            if codes !== nothing && codes[c] in xcodes
                excl[gi(lo[c]), gj(la[c])] = true; continue
            end
            g = get(dgrp, dom[c], 0); g == 0 && continue
            for (k, r) in enumerate(regs)
                sg.names[g] in r["groups"] && (core[gi(lo[c]), gj(la[c]), k] = true)
            end
        end
    end
    taper(d, b, t) = d <= b ? 1.0 : d >= b + t ? 0.0 : cos(π / 2 * (d - b) / t)^2
    function near(mask, f)                            # f(distance to the cells of mask) on the 1-degree grid
        cc = [(lon[I[1]], lat[I[2]]) for I in findall(mask)]
        A = zeros(length(lon), length(lat))
        Threads.@threads for j in eachindex(lat)
            for i in eachindex(lon); A[i, j] = f(minimum(hav(lon[i], lat[j], a, b) for (a, b) in cc)); end
        end
        A, length(cc)
    end
    T = zeros(length(lon), length(lat), length(regs))
    for k in eachindex(regs)
        any(core[:, :, k]) || error("region $(regs[k]["name"]): no cells")
        T[:, :, k], n = near(core[:, :, k], d -> taper(d, buf, tap))
        println(regs[k]["name"], ": ", n, " core cells (1 deg)")
    end
    W = T ./ max.(1.0, sum(T; dims = 3))
    if any(excl)                                      # no regional shift in the excluded areas, tapered off around them
        X, n = near(excl, d -> taper(d, 0.0, xtap)); W .*= 1 .- X
        println("excluded ($xvar codes $(xcodes)): ", n, " cells (1 deg), taper ", xtap, " km")
    end
    isfile(out) && rm(out)
    NCDataset(out, "c") do ds
        defVar(ds, "lon", lon, ("lon",); attrib = ["units" => "degrees_east"])
        defVar(ds, "lat", lat, ("lat",); attrib = ["units" => "degrees_north"])
        defVar(ds, "mask", W, ("lon", "lat", "region"); attrib = ["units" => "1",
               "long_name" => "weight of the region's upper-mantle viscosity shift"])
        ds.attrib["region_names"] = join([String(r["name"]) for r in regs], ";")
        ds.attrib["region_groups"] = join([String(r["name"]) * ": " * join(r["groups"], ", ") for r in regs], "; ")
        ds.attrib["buffer_km"] = buf; ds.attrib["taper_km"] = tap
        isempty(xcodes) || (ds.attrib["excluded"] = "$xvar codes $(join(xcodes, ",")), taper $xtap km")
        ds.attrib["source"] = "plastic_earth_shift.jl regions: the tau_b domains of the scale groups on " *
                              join(basename.(haskey(sc, "taub_files") ? sc["taub_files"] : [sc["taub_file"]]), ", ")
    end
    println("wrote ", out)
end

# --- the radial mean of the 3-D field ------------------------------------------------------------------------------------------
# A 1-D Earth from a 3-D par file (l_visc_3d with the FesmData Bagge field, visc_3d_structure): the viscosity of each Maxwell layer
# (r_top, r_bot of its &vilma group) is the mean of log10 eta over the layer, weighted by area (cos lat) and by r^2 in depth, of
# the field clamped to [visc_log10_min, visc_log10_max] as VILMA2 clamps it (the field linear in depth between its levels, sampled
# every 1 km). Writes a copy of the par file with that eta, l_visc_3d = .false. and the spin-up of the 1-D Earths (pre_spinup_1d
# = .false., equil_time_max = 2e6 yr, as Gowan's B1): the reference that separates the 3-D field's radial structure from its
# lateral one.
#   julia --project=tools/levante tools/levante/plastic_ensemble.jl earth-radial <3-D par.nml> <out par.nml>     ($FESMDATA)
function nml_group(file, grp)
    s = read(file, String); m = match(Regex("^\\s*&" * grp * "\\s*\$(.*?)^\\s*/\\s*\$", "ims"), s); m === nothing && error("no &$grp in $file")
    d = Dict{String,String}()
    for l in split(m[1], "\n")
        l = strip(replace(l, r"!.*$" => "")); occursin("=", l) || continue
        k, v = strip.(split(l, "=", limit = 2)); d[lowercase(k)] = v
    end
    d
end
nml_reals(v) = parse.(Float64, replace.(strip.(split(v, ",")), r"[dD]" => "e"))
function earth_radial(par, out)
    v = nml_group(par, "vilma"); dflt = nml_group(joinpath(@__DIR__, "..", "..", "input", "vilma_defaults.nml"), "vilma")
    get(v, "l_visc_3d", ".false.") == ".true." || error("$par: no 3-D field (l_visc_3d)")
    g(k) = haskey(v, k) ? v[k] : dflt[k]
    rt = nml_reals(v["r_top"]) ./ 1e3; rb = nml_reals(v["r_bot"]) ./ 1e3; eta = nml_reals(v["eta"]); rh = parse.(Int, strip.(split(v["rheology"], ",")))
    lo, hi = nml_reals(g("visc_log10_min"))[1], nml_reals(g("visc_log10_max"))[1]
    f = replace(strip(g("visc_3d_file"), '"'), r"^fesmdata/" => ENV["FESMDATA"] * "/"); st = parse(Int, g("visc_3d_structure"))
    lat, dep, F = NCDataset(f) do ds
        Float64.(ds["lat"][:]), Float64.(ds["depth"][:]) ./ 1e3, clamp.(Float64.(ds["log10_visc"][:, :, :, st]), lo, hi)
    end
    w = cosd.(lat); A = [sum(F[:, j, k] .* w[j] for j in eachindex(lat)) / (size(F, 1) * sum(w)) for k in eachindex(dep)]   # area mean per level
    p = sortperm(dep); dep = dep[p]; A = A[p]
    at(d) = (d <= dep[1] ? A[1] : d >= dep[end] ? A[end] : (k = searchsortedlast(dep, d); x = (d - dep[k]) / (dep[k+1] - dep[k]); (1 - x) * A[k] + x * A[k+1]))
    new = copy(eta)
    for i in eachindex(eta)
        rh[i] == 1 || continue
        r = collect(range(rb[i], rt[i], length = max(2, round(Int, rt[i] - rb[i]) + 1))); wr = r .^ 2
        new[i] = 10.0^(sum(at(R_VILMA_KM - x) * y for (x, y) in zip(r, wr)) / sum(wr))
    end
    cp(par, out; force = true)
    fmt(x) = @sprintf("%.3e", x)
    Base.run(`python3 $(joinpath(@__DIR__, "set_nml.py")) $out vilma.eta=$(join(fmt.(new), ", ")) vilma.l_visc_3d=.false. vilma.pre_spinup_1d=.false. vilma.equil_time_max=2.0e6`)
    hdr = "! The radial mean of the 3-D Earth of $(basename(par)) (plastic_ensemble.jl earth-radial): eta per Maxwell layer = the mean of\n" *
          "! log10 eta of $(basename(f)) structure $st over the layer (area, r^2), clamped to [$lo, $hi]; l_visc_3d off, the 1-D spin-up.\n"
    write(out, hdr * read(out, String))
    for i in eachindex(eta); rh[i] == 1 && @printf("layer %2d  %6.0f-%6.0f km  log10 eta %.2f (par background %.2f)\n", i, R_VILMA_KM - rt[i], R_VILMA_KM - rb[i], log10(new[i]), log10(eta[i])); end
    println("wrote ", out)
end
