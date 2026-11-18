# A margin patch of the plastic ice backend (margin_patches; docs/backend-plastic.md#margin-patches) from a dated ice-margin
# product: a chronology on the grid of the base margin (PaleoMIST 1.0, Gowan et al. 2021, the instance's margin_file) that
# replaces the base's in some regions of the region map within a time window, the base's elsewhere and at other times untouched.
#
#   julia --project=tools/levante tools/levante/margin_patch.jl <product> <var> <region ids> <t0> <t1> <out> [<figure folder>]
#         [--base=<file>] [--regions=<file>] [--band=<cells>]
#
#   <product>     the product file (FesmData IceHistory), of one of two kinds, told by the dimensions of <var>:
#                 slices  var(xc, yc, time) on the base's grid, time in years ago (Dalton2020, DATED1, ...: var f_ice, or
#                         f_ice_min, f_ice_max): its slices within the window;
#                 ages    var(xc, yc) = age_deglac, the time of deglaciation of the areas between isochrones on a finer grid
#                         of the same projection (PaleoGrIS 1.0, Leger et al. 2024: GRL-PAL-1KM_IceHistory-PaleoGrIS.nc, with
#                         age_deglac_err and, in the same folder, GRL-PAL-1KM_TOPO-BedMachine-v5.nc): one slice per
#                         isochrone within the window (below);
#   <region ids>  ids of the region map, comma separated (e.g. 1 or 3,6): the patch's region, where kryos uses it online;
#   <t0> <t1>     the window [years ago], t0 > t1: the patch's slices are the base's at t0 and t1 (its slice or, between its
#                 slices, its mask interpolated by the signed distance as kryos does), so that the patch hands over to the base
#                 continuously at both ends, and the product's slices strictly inside (t1, t0) (at an end, the base's wins);
#   <out>         the patch (NetCDF), and next to it <out>_check.md and <out>_check.png (also copied to <figure folder>);
#   --base        the base chronology (default $FESMDATA/North/NH-32KM/NH-32KM_IceHistory-PaleoMIST.nc, f_ice > 0.5);
#   --regions     the region map on the base's grid (default $KRYOS_DATA/plastic/<grid>_taub_regions_v2.nc, variable region);
#   --band        the seam band, cells (default 2).
#
# The slices. At a product time a slice is the product's ice fraction in its domain, the base's (sdist-interpolated) mask
# elsewhere, and in a band of --band cells either side of the domain's boundary with other ice (the footprint of the base or
# of the product) the larger of the two, so that the patched ice and the ice beyond do not part at the seam. The domain is the
# patch's region for a slices product; for the ages product, Greenland (below). Online only the region is used: outside it the
# slices carry the base's ice, so that the signed distance of the patch near the region's boundary sees the ice beyond it.
#
# Use it in the ice instance on the base's grid (margin_file the base, margin_var 'f_ice', margin_H_min 0.5):
#   margin_patches = 'margin_x'
#   &margin_x  file = '<out>', var = 'f_ice', H_min = 0.5,
#              region_file = 'kryos_data/plastic/{grid_name}_taub_regions_v2.nc', region_var = 'region', region_ids = '<ids>'  /
#
# Examples ($S: a FesmData products folder with the dated margins):
#   Greenland from PaleoGrIS, 15-5 ka, region 2 (the patch of C4):
#     margin_patch.jl $FESMDATA/GreenlandPaleo/GRL-PAL-1KM/GRL-PAL-1KM_IceHistory-PaleoGrIS.nc age_deglac 2 15000 5000 <out>
#   North America from Dalton et al. (2020), Eurasia from DATED-1 (most-credible; Hughes et al. 2016):
#     margin_patch.jl $S/North/NH-32KM/NH-32KM_IceHistory-Dalton2020.nc f_ice 1 25000 20000 <out>
#     margin_patch.jl $S/North/NH-32KM/NH-32KM_IceHistory-DATED1.nc f_ice 3,6 27500 7500 <out>
#
# The ages product (PaleoGrIS). Each isochrone is a line with a time range (14-13, 13-12, 12-11.5, then every 500 years to
# 7-6.5 ka); each area between two isochrones has age_deglac = the middle of [younger bound of its inner isochrone, older bound
# of its outer one] and age_deglac_err = half that span (13000 +- 1000, 12250 +- 750, 11500 +- 500, ..., 7000 +- 500; checked
# here). The slice of an isochrone is taken at the middle of its range (13.5, 12.5, 11.75, 11.25, ..., 7.25, 6.75 ka): ice where
# the area deglaciated after it (age_deglac < t), under today's ice (BedMachine), and in the cells the reconstruction encloses
# without an age (inside its innermost isochrone, not today's ice); the 1 km fields are aggregated exactly by the overlap of
# the cells (same projection, checked). Greenland, its domain, is the cells nearer to Greenland's land (the 1 km land component
# with most of today's ice and the islands whose centroid lies within 120 km of it) than to any other land, within 400 km of it.
#
# Output: f_ice(xc, yc, time) (time in years ago), domain(xc, yc) (0 outside, 1 the product's domain, 2 the band), region(xc, yc)
# (1 the patch's region), source(time) (the domain from 0 the base, 1 the product). The check emulates the set online
# (margin_set_at: the base everywhere, the patch in the region within the window) every 100 yr from 1 kyr before t0 to 1 kyr
# after t1: outside the region it is the base exactly (by construction), and it reports the product's ice of the slices that
# falls outside the region (not used online) and, at each product slice, the product's and the base's ice in the region.
using CairoMakie, NCDatasets, Printf, Statistics, Dates

# --- arguments ----------------------------------------------------------------------------------------------------------------
const OPTS = ("base", "regions", "band")
opts = Dict{String,String}(); pos = String[]
for a in ARGS
    m = match(r"^--([a-z]+)=(.*)$", a)
    m === nothing ? push!(pos, a) : (opts[m[1]] = m[2])
end
issubset(keys(opts), OPTS) || error("unknown options $(setdiff(keys(opts), OPTS)) (available: $(join(OPTS, ", ")))")
length(pos) in (6, 7) || error("usage: margin_patch.jl <product> <var> <region ids> <t0> <t1> <out> [<figure folder>] " *
                               "[--base=<file>] [--regions=<file>] [--band=<cells>]")
f_prod, VAR = pos[1], pos[2]
const REGION_IDS = parse.(Int, split(pos[3], ','))
const T0, T1 = parse(Float64, pos[4]), parse(Float64, pos[5])
T0 > T1 || error("t0 = $T0 must be older (larger, years ago) than t1 = $T1")
out = pos[6]; figdir = length(pos) == 7 ? pos[7] : ""
const T_END = (T0, T1)                                                                                     # the base's slices
const BAND = parse(Int, get(opts, "band", "2"))                                                            # cells
f_pm = get(opts, "base", joinpath(ENV["FESMDATA"], "North", "NH-32KM", "NH-32KM_IceHistory-PaleoMIST.nc"))
const GRID = split(basename(f_pm), '_')[1]
f_rg = get(opts, "regions", joinpath(ENV["KRYOS_DATA"], "plastic", "$(GRID)_taub_regions_v2.nc"))
mkpath(dirname(abspath(out)))

# --- signed distance and its interpolation in time, as kryos_plastic_margin -----------------------------------------------------
function envelope_1d(f, s)                      # d(q) = min_p ((q - p) s)^2 + f(p) (Felzenszwalb and Huttenlocher 2012)
    n = length(f); d = fill(Inf, n); v = zeros(Int, n); z = zeros(n + 1); m = 0
    for q in 1:n
        isfinite(f[q]) || continue
        x = (q - 1) * s
        if m == 0; m = 1; v[1] = q; z[1] = -Inf; z[2] = Inf; continue; end
        xs = 0.0
        while true
            p = (v[m] - 1) * s; xs = ((f[q] + x^2) - (f[v[m]] + p^2)) / (2 * (x - p))
            xs > z[m] && break
            m -= 1
        end
        m += 1; v[m] = q; z[m] = xs; z[m + 1] = Inf
    end
    m == 0 && return d
    k = 1
    for q in 1:n
        x = (q - 1) * s
        while z[k + 1] < x; k += 1; end
        d[q] = (x - (v[k] - 1) * s)^2 + f[v[k]]
    end
    d
end
function edt2(fg, dx, dy)
    d = [fg[c] ? 0.0 : Inf for c in CartesianIndices(fg)]
    for j in axes(d, 2); d[:, j] = envelope_1d(d[:, j], dx); end
    for i in axes(d, 1); d[i, :] = envelope_1d(d[i, :], dy); end
    d
end
function sdist(ice, dx, dy)                     # [m] + inside, zero level midway between an ice cell and its ice-free neighbour
    h = min(dx, dy) / 2; L2 = (size(ice, 1) * dx)^2 + (size(ice, 2) * dy)^2
    dout = sqrt.(min.(edt2(.!ice, dx, dy), L2)); din = sqrt.(min.(edt2(ice, dx, dy), L2))
    ifelse.(ice, dout .- h, .-(din .- h))
end
# the mask at t of slices (times ts [yr ago], signed distances SD[k]): linear in time between the two that bracket t
function mask_at(ts, SD, t)
    k = findfirst(==(t), ts); k === nothing || return SD[k] .> 0
    p = sortperm(ts); tt = ts[p]
    (tt[1] < t < tt[end]) || error("t = $t outside the slices [$(tt[1]), $(tt[end])]")
    j = searchsortedlast(tt, t); a = (t - tt[j]) / (tt[j + 1] - tt[j])
    ((1 - a) .* SD[p[j]] .+ a .* SD[p[j + 1]]) .> 0
end

# --- grids and masks ------------------------------------------------------------------------------------------------------------
function wts(xt, dt, xs, ds)                    # W[i, j]: fraction of target cell i covered by source cell j, along one axis
    W = zeros(length(xt), length(xs))
    for (i, a) in enumerate(xt), (j, b) in enumerate(xs)
        o = min(a + dt / 2, b + ds / 2) - max(a - dt / 2, b - ds / 2); o > 0 && (W[i, j] = o / dt)
    end
    W
end
proj(f) = NCDataset(ds -> replace(ds["crs"].attrib["proj_params"], r"\s+" => " "), f)
function components(land)                       # 4-connected components of a mask
    lab = zeros(Int32, size(land)); n = 0; st = Tuple{Int,Int}[]
    for c in CartesianIndices(land)
        (land[c] && lab[c] == 0) || continue
        n += 1; lab[c] = n; push!(st, (c[1], c[2]))
        while !isempty(st)
            i, j = pop!(st)
            for (a, b) in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1))
                (1 <= a <= size(land, 1) && 1 <= b <= size(land, 2) && land[a, b] && lab[a, b] == 0) || continue
                lab[a, b] = n; push!(st, (a, b))
            end
        end
    end
    lab, n
end
function enclosed(free)                         # cells of a mask not reachable from the edge of the grid through it
    nx, ny = size(free); reach = falses(nx, ny); st = Tuple{Int,Int}[]
    for i in 1:nx, j in 1:ny
        (i in (1, nx) || j in (1, ny)) && free[i, j] && (reach[i, j] = true; push!(st, (i, j)))
    end
    while !isempty(st)
        i, j = pop!(st)
        for (a, b) in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1))
            (1 <= a <= nx && 1 <= b <= ny && free[a, b] && !reach[a, b]) || continue
            reach[a, b] = true; push!(st, (a, b))
        end
    end
    free .& .!reach
end

# --- the base on the grid, the region -----------------------------------------------------------------------------------------
ds = NCDataset(f_pm)
x, y = Float64.(ds["xc"][:]), Float64.(ds["yc"][:]); dx, dy = 1e3 * (x[2] - x[1]), 1e3 * (y[2] - y[1])
occursin("ago", ds["time"].attrib["units"]) || error("$f_pm: time is not in years ago")
t_pm = Float64.(ds["time"][:]); F_pm = Float64.(coalesce.(ds["f_ice"][:, :, :], 0.0))
close(ds)
all(t -> minimum(t_pm) <= t <= maximum(t_pm), T_END) || error("$f_pm: the window [$T1, $T0] is not within its slices")
SD_pm = [sdist(F_pm[:, :, k] .> 0.5, dx, dy) for k in eachindex(t_pm)]
reg = NCDataset(ds -> Int.(ds["region"][:, :]), f_rg)
size(reg) == size(F_pm)[1:2] || error("$f_rg is not on $GRID")
for id in REGION_IDS; any(==(id), reg) || error("$f_rg has no region $id"); end
R = in.(reg, Ref(REGION_IDS))
inwin(t) = T1 < t < T0

# --- the product: its slices in the window on the grid, its domain ------------------------------------------------------------
# slices: var(xc, yc, time) on the base's grid; the domain is the region
function product_slices()
    NCDataset(f_prod) do ds
        v = ds[VAR]
        dimnames(v) == ("xc", "yc", "time") || error("$f_prod: $VAR has dimensions $(dimnames(v)), not (xc, yc, time)")
        (Float64.(ds["xc"][:]) == x && Float64.(ds["yc"][:]) == y) || error("$f_prod is not on the grid of $f_pm")
        haskey(ds, "crs") && proj(f_prod) != proj(f_pm) && error("$f_prod and $f_pm are not on the same projection")
        occursin("ago", ds["time"].attrib["units"]) || error("$f_prod: time is not in years ago")
        tp = Float64.(ds["time"][:]); k = findall(inwin, tp)
        isempty(k) && error("$f_prod: no slice within ($T1, $T0)")
        P = Float64.(coalesce.(v[:, :, k], 0.0))
        lab = replace(basename(f_prod), r"^.*_IceHistory-" => "", r"\.nc$" => "") * (VAR == "f_ice" ? "" : " $VAR")
        (t = tp[k], P = P, D = copy(R), cov = ones(size(R)), label = lab,
         comment = "$(basename(f_prod)) ($(get(ds.attrib, "title", ""))): its $VAR at its slices within the window, in the region",
         notes = String[], refs = get(ds.attrib, "references", ""))
    end
end

# ages: PaleoGrIS's age_deglac on its 1 km grid, one slice per isochrone, aggregated; the domain is Greenland
const ISO = [(14000, 13000), (13000, 12000), (12000, 11500), [(t + 500, t) for t in 11000:-500:6500]...]   # isochrone ranges [yr ago]
function product_ages()
    VAR == "age_deglac" || error("$f_prod: a 2-D product must be the ages of deglaciation (age_deglac), not $VAR")
    proj(f_pm) == proj(f_prod) || error("$f_pm and $f_prod are not on the same projection")
    f_bm = joinpath(dirname(f_prod), replace(basename(f_prod), r"_IceHistory-.*$" => "_TOPO-BedMachine-v5.nc"))
    T_ISO = [(a + b) / 2 for (a, b) in ISO]                                                                # their slice times
    age, err, x1, y1 = NCDataset(f_prod) do ds
        Float64.(coalesce.(ds["age_deglac"][:, :], NaN)), Float64.(coalesce.(ds["age_deglac_err"][:, :], NaN)), Float64.(ds["xc"][:]),
        Float64.(ds["yc"][:])
    end
    bm = NCDataset(ds -> Int.(coalesce.(ds["mask"][:, :], 0)), f_bm)
    va = .!isnan.(age); pd = (bm .== 2) .| (bm .== 3)
    # the convention: each area's age is the middle of [the younger bound of its inner isochrone, the older bound of its outer
    # one], its uncertainty half that span; checked on the cells within 30 yr of such a value (the rest are blends of two
    # neighbouring areas by the remapping, about 10%), which must be most of them (ages are bit-rounded: 12248 for 12250)
    conv = let L = [((ISO[k + 1][2] + ISO[k][1]) / 2, (ISO[k][1] - ISO[k + 1][2]) / 2) for k in 1:length(ISO) - 1], n = 0, ok = 0
        for c in findall(va)
            n += 1; k = findfirst(l -> abs(age[c] - l[1]) < 30, L)
            k !== nothing && abs(err[c] - L[k][2]) < 5 && (ok += 1)
        end
        ok / n > 0.8 || error("$f_prod: only $(round(100ok / n))% of the cells have the age and uncertainty of an area between " *
                              "two isochrones of $ISO: the convention of this tool does not hold")
        ok / n
    end
    println("PaleoGrIS convention: $(round(100conv, digits = 1))% of the cells match an area between two isochrones")
    encl = enclosed(.!va .& .!pd)
    Wx = wts(x, dx / 1e3, x1, 1.0); Wy = wts(y, dy / 1e3, y1, 1.0)
    agg(F) = Wx * F * Wy'
    cov = agg(ones(size(age)))
    pg1(t) = (va .& (age .< t)) .| (.!va .& (pd .| encl))
    tw = filter(inwin, T_ISO)
    isempty(tw) && error("$f_prod: no isochrone within ($T1, $T0)")
    P = cat([agg(Float64.(pg1(t))) for t in tw]...; dims = 3)
    # Greenland on the grid
    lab, nl = components(bm .!= 0)
    gid = argmax([count(pd .& (lab .== k)) for k in 1:nl])
    sx = zeros(nl); sy = zeros(nl); sn = zeros(nl)
    for c in CartesianIndices(lab); k = lab[c]; k > 0 || continue; sx[k] += x1[c[1]]; sy[k] += y1[c[2]]; sn[k] += 1; end
    ML = [(x1[c[1]], y1[c[2]]) for c in CartesianIndices(lab) if lab[c] == gid && c[1] % 5 == 0 && c[2] % 5 == 0]
    isg = [k == gid || minimum((sx[k] / sn[k] - p[1])^2 + (sy[k] / sn[k] - p[2])^2 for p in ML) <= 120.0^2 for k in 1:nl]
    gl = agg(Float64.([l > 0 && isg[l] for l in lab])) .> 0.5; ol = agg(Float64.([l > 0 && !isg[l] for l in lab])) .> 0.5
    G = falses(length(x), length(y))
    let Gp = [(x[c[1]], y[c[2]]) for c in findall(gl)], Op = [(x[c[1]], y[c[2]]) for c in findall(ol)]
        for c in findall(cov .> 0.99)
            dg = minimum((x[c[1]] - p[1])^2 + (y[c[2]] - p[2])^2 for p in Gp)
            do_ = isempty(Op) ? Inf : minimum((x[c[1]] - p[1])^2 + (y[c[2]] - p[2])^2 for p in Op)
            G[c] = dg < do_ && dg <= 400.0^2
        end
    end
    (t = tw, P = P, D = G, cov = cov, label = "PaleoGrIS",
     comment = "PaleoGrIS 1.0 slices at the middle of each isochrone's range ($(join(Int.(tw), ", ")) years ago): ice where " *
               "age_deglac < t, under today's ice (BedMachine v5) and in the cells the reconstruction encloses without an age, " *
               "1 km aggregated; its domain Greenland: nearer to Greenland's land than to any other, within 400 km",
     notes = ["PaleoGrIS: enclosed cells without an age $(count(encl)) km2 (1 km); convention $(round(100conv, digits = 1))%"],
     refs = "Leger et al. 2024, Clim. Past 20, 701-755, doi:10.5194/cp-20-701-2024")
end
kind = NCDataset(ds -> ndims(ds[VAR]), f_prod) == 3 ? "slices" : "ages"
pr = kind == "slices" ? product_slices() : product_ages()
G = pr.D; cov = pr.cov

# the band: a cell with a cell of the other side within BAND cells, only where it separates the domain from other ice (the
# seam), not along an open-ocean edge
B = falses(size(G))
for c in CartesianIndices(G)
    cov[c] > 0.99 || continue
    for di in -BAND:BAND, dj in -BAND:BAND
        a, b = c[1] + di, c[2] + dj
        (1 <= a <= size(G, 1) && 1 <= b <= size(G, 2)) || continue
        G[a, b] != G[c] && (B[c] = true; break)
    end
end
fp_pm = dropdims(any(F_pm .> 0.5, dims = 3), dims = 3); fp_pr = dropdims(any(pr.P .> 0.5, dims = 3), dims = 3)
B .&= (fp_pm .| fp_pr)

# --- the patch's slices -------------------------------------------------------------------------------------------------------
t_p = sort([T_END...; pr.t], rev = true)
fpm(t) = (k = findfirst(==(t), t_pm); k === nothing ? Float64.(mask_at(t_pm, SD_pm, t)) : F_pm[:, :, k])
F = zeros(length(x), length(y), length(t_p)); src = zeros(Int8, length(t_p))
for (k, t) in enumerate(t_p)
    a = fpm(t)
    j = findfirst(==(t), pr.t)
    if j !== nothing
        b = pr.P[:, :, j]; a = ifelse.(G .& .!B, b, ifelse.(B, max.(a, b), a)); src[k] = 1
    end
    F[:, :, k] = a
end
dmap = Int8.(ifelse.(B, 2, ifelse.(G, 1, 0)))
ids = join(REGION_IDS, ",")
ids_nml = join(["'$i'" for i in REGION_IDS], ", ")

rm(out; force = true)
NCDataset(out, "c") do o
    NCDataset(f_pm) do s
        defVar(o, "xc", x, ("xc",); attrib = Dict(s["xc"].attrib)); defVar(o, "yc", y, ("yc",); attrib = Dict(s["yc"].attrib))
        defVar(o, "crs", Int32(0), (); attrib = Dict(s["crs"].attrib))
    end
    defVar(o, "time", t_p, ("time",); attrib = ["units" => "years ago", "long_name" => "time (years before 1950 CE)"])
    defVar(o, "f_ice", Float32.(F), ("xc", "yc", "time"); deflatelevel = 1, attrib = [
        "units" => "1", "long_name" => "ice area fraction", "grid_mapping" => "crs",
        "comment" => "$(Int(T0)) and $(Int(T1)) years ago: the base's f_ice (or its sdist-interpolated mask); the product's " *
                     "times: $(pr.label) in its domain, the base's sdist-interpolated mask (0 or 1) elsewhere, the larger of both " *
                     "in the band"])
    defVar(o, "domain", dmap, ("xc", "yc"); attrib = ["long_name" => "the product's domain in the patch's slices",
        "flag_values" => Int8[0, 1, 2], "flag_meanings" => "outside domain seam_band_max_of_both", "grid_mapping" => "crs"])
    defVar(o, "region", Int8.(R), ("xc", "yc"); attrib = ["long_name" => "the patch's region (region ids $ids of " *
        "$(basename(f_rg))): online the slices are used here only", "grid_mapping" => "crs"])
    defVar(o, "source", src, ("time",); attrib = ["long_name" => "source of the slice in the domain", "flag_values" => Int8[0, 1],
        "flag_meanings" => "base product"])
    o.attrib["title"] = "Margin patch on $GRID: $(pr.label) in regions $ids, $(T0 / 1000) to $(T1 / 1000) ka " *
                        "($(basename(f_pm)) at both ends)"
    o.attrib["usage"] = "kryos plastic backend, ice instance on $GRID with margin_file $(basename(f_pm)): margin_patches = " *
                        "'margin_x'; &margin_x file = '$out', var = 'f_ice', H_min = 0.5, region_file = " *
                        "'kryos_data/plastic/{grid_name}_taub_regions_v2.nc', region_var = 'region', region_ids = $ids_nml"
    o.attrib["comment"] = "$(pr.comment). Band of $BAND cells at the domain's boundary with other ice: max of both. Made by " *
                          "tools/levante/margin_patch.jl (kryos): margin_patch.jl $(join(pos[1:5], ' ')) <out>"
    o.attrib["references"] = "Gowan et al. 2021, Nat. Commun. 12, 1199, doi:10.1038/s41467-021-21469-w (base); $(pr.refs)"
    o.attrib["sources"] = join([f_pm, f_prod, f_rg], ", ")
    o.attrib["history"] = "$(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")): margin_patch.jl"
end
println("wrote $out: $(length(t_p)) slices, domain $(count(G)) cells, band $(count(B)) cells, region $(count(R)) cells")

# --- check: the set online (margin_set_at), areas, maps ------------------------------------------------------------------------
SD_p = [sdist(F[:, :, k] .> 0.5, dx, dy) for k in eachindex(t_p)]
SD_pr = [sdist(pr.P[:, :, k] .> 0.5, dx, dy) for k in eachindex(pr.t)]
set_at(t) = (m = mask_at(t_pm, SD_pm, t); T1 <= t <= T0 ? ifelse.(R, mask_at(t_p, SD_p, t), m) : m)
pr_at(t) = length(pr.t) == 1 ? SD_pr[1] .> 0 : mask_at(pr.t, SD_pr, t)
dA = dx * dy / 1e9                               # 10^3 km2
tt = collect(min(T0 + 1000, maximum(t_pm)):-100.0:max(T1 - 1000, minimum(t_pm)))
win(t) = minimum(pr.t) <= t <= maximum(pr.t)
A_pm = [count(mask_at(t_pm, SD_pm, t) .& R) * dA for t in tt]; A_s = [count(set_at(t) .& R) * dA for t in tt]
A_pr = [win(t) ? count(pr_at(t) .& R) * dA : NaN for t in tt]
d_out = [count(.!R .& (set_at(t) .!= mask_at(t_pm, SD_pm, t))) for t in tt]
fp_p = dropdims(any(F .> 0.5, dims = 3), dims = 3)
io = open(replace(out, r"\.nc$" => "_check.md"), "w")
println(io, "# $(basename(out))\n\n$(pr.label) in regions $ids of $(basename(f_rg)), window $(T0 / 1000) to $(T1 / 1000) ka\n")
println(io, "slices: $(join(Int.(t_p), ", ")) years ago (the domain from the product at ", join(Int.(t_p[src .== 1]), ", "), ")\n")
println(io, "regions $ids (the patch's): $(count(R)) cells; domain $(count(G)) cells ($(count(G .& R)) in the region), band ",
        "$(count(B)) cells ($(count(B .& R)) in the region)\n")
foreach(n -> println(io, n, "\n"), pr.notes)
println(io, "product ice of the slices outside the region (not used online): ", count(fp_pr .& G .& .!R), " cells\n")
println(io, "patch footprint in the region: $(count(fp_p .& R)) cells, the base's $(count(fp_pm .& R)); new: ",
        "$(count(fp_p .& R .& .!fp_pm))\n")
println(io, "the set (the base everywhere, the patch in the region $(T0 / 1000)-$(T1 / 1000) ka) against the base outside the ",
        "region, every 100 yr: $(sum(d_out)) cells differ (by construction 0)\n")
println(io, "At the product's slices in the window, ice in the region [10³ km²]: the product, the base (interpolated), their ",
        "symmetric difference and intersection over union\n")
println(io, "| t ka | product | base | sym. diff. | IoU |\n|---|---|---|---|---|")
for (k, t) in enumerate(pr.t)
    a = (pr.P[:, :, k] .> 0.5) .& R; b = mask_at(t_pm, SD_pm, t) .& R
    println(io, @sprintf("| %.2f | %.0f | %.0f | %.0f | %.3f |", t / 1000, count(a) * dA, count(b) * dA, count(a .⊻ b) * dA,
                         count(a .& b) / max(1, count(a .| b))))
end
println(io, "\nIce in the region through the window (kryos's interpolation) [10³ km²]\n")
println(io, "| t ka | base | product | the set |\n|---|---|---|---|")
for (i, t) in enumerate(tt)
    mod(t, 500) == 0 && println(io, @sprintf("| %.1f | %.0f | %.0f | %.0f |", t / 1000, A_pm[i], A_pr[i], A_s[i]))
end
close(io)
println(read(replace(out, r"\.nc$" => "_check.md"), String))
sum(d_out) == 0 || error("the set differs from the base outside regions $ids")

# the maps: the bounding box of the region and the domain, 200 km around
bb = findall(R .| G)
lim = ((x[minimum(c[1] for c in bb)] - 200, x[maximum(c[1] for c in bb)] + 200),
       (y[minimum(c[2] for c in bb)] - 200, y[maximum(c[2] for c in bb)] + 200))
fig = Figure(size = (1800, 1050))
ax = Axis(fig[1, 1:3], xreversed = true, xlabel = "ka", ylabel = "10³ km²",
          title = "ice area in regions $ids (kryos's interpolation, 100 yr)", limits = (extrema(tt ./ 1000), nothing))
lines!(ax, tt ./ 1000, A_pm, color = :black, linewidth = 2, label = "base ($(basename(f_pm)))")
lines!(ax, tt ./ 1000, A_pr, color = :red, linewidth = 2, label = pr.label)
lines!(ax, tt ./ 1000, A_s, color = :royalblue, linestyle = :dash, linewidth = 2.5, label = "base + patch")
vlines!(ax, t_p ./ 1000, color = (:red, 0.25)); vlines!(ax, t_pm ./ 1000, color = (:black, 0.25), linestyle = :dot)
axislegend(ax, position = :rt)
ax = Axis(fig[1, 4:5], title = "domain (1), band (2); regions $ids (green), patch footprint (black)", aspect = DataAspect(),
          limits = lim)
heatmap!(ax, x, y, ifelse.(fp_p .| (dmap .> 0), Float64.(dmap), NaN), colormap = [:gray80, :royalblue, :orange], colorrange = (0, 2))
contour!(ax, x, y, Float64.(fp_p .& R), levels = [0.5], color = :black)
contour!(ax, x, y, Float64.(R), levels = [0.5], color = :green, linewidth = 1.5)
tm = t_p[unique(round.(Int, range(1, length(t_p), length = 5)))]
for (p, t) in enumerate(tm)
    m = mask_at(t_pm, SD_pm, t); g = set_at(t)
    c = fill(NaN, size(m)); c[m .& g] .= 1; c[m .& .!g] .= 2; c[g .& .!m] .= 3
    local ax = Axis(fig[2, p], aspect = DataAspect(), title = "$(t / 1000) ka", limits = lim,
              xticklabelsvisible = false, yticklabelsvisible = false)
    heatmap!(ax, x, y, c, colormap = cgrad([:gray75, :royalblue, :darkorange], 3, categorical = true), colorrange = (1, 3))
    contour!(ax, x, y, Float64.(R), levels = [0.5], color = :green, linewidth = 0.8)
end
Legend(fig[3, 1:5], [PolyElement(color = :gray75), PolyElement(color = :royalblue), PolyElement(color = :darkorange), LineElement(color = :green)],
       ["both", "base only", "base + patch only", "regions $ids"], orientation = :horizontal, framevisible = false)
png = replace(out, r"\.nc$" => "_check.png")
save(png, fig)
isempty(figdir) || (mkpath(figdir); cp(png, joinpath(figdir, basename(png)); force = true))
