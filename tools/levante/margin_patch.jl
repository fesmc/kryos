# Greenland's margin from PaleoGrIS 1.0 (Leger et al. 2024) as a margin patch of the plastic ice backend (margin_patches;
# docs/backend-plastic.md): a chronology on NH-32KM that replaces PaleoMIST 1.0's (Gowan et al. 2021, margin_file) in
# Greenland's region from 15 to 5 ka, PaleoMIST's elsewhere and at other times untouched.
#
#   julia --project=tools/levante tools/levante/margin_patch.jl [<out file> [<figure folder>]]
#       default $KRYOS_WORK/margins/NH-32KM/NH-32KM_IceHistory-PaleoGrIS-patch.nc, and next to it <out>_check.png
#       (Greenland's ice area through time and maps; also copied to <figure folder>) and <out>_check.md (the numbers)
#
# Use it in the ice instance on NH-32KM (margin_file PaleoMIST's, margin_var 'f_ice', margin_H_min 0.5):
#   margin_patches = 'margin_grl'
#   &margin_grl  file = '<out file>', var = 'f_ice', H_min = 0.5,
#                region_file = 'kryos_data/plastic/{grid_name}_taub_regions_v2.nc', region_var = 'region', region_ids = '2'  /
#
# Inputs (FesmData, read only): {NH-32KM}_IceHistory-PaleoMIST.nc (f_ice, 33 slices 80-0 ka), GRL-PAL-1KM_IceHistory-PaleoGrIS.nc
# (age_deglac, age_deglac_err: the areas between two isochrones), GRL-PAL-1KM_TOPO-BedMachine-v5.nc (mask: today's ice), and the
# region map v2 ($KRYOS_DATA/plastic/NH-32KM_taub_regions_v2.nc, id 2 Greenland: the patch's region). The grids share one polar
# stereographic projection (checked), so the 1 km fields are aggregated exactly by the overlap of the cells.
#
# PaleoGrIS's ages. Each isochrone is a line with a time range (14-13, 13-12, 12-11.5, then every 500 years to 7-6.5 ka); each
# area between two isochrones has age_deglac = the middle of [younger bound of its inner isochrone, older bound of its outer one]
# and age_deglac_err = half that span (13000 +- 1000, 12250 +- 750, 11500 +- 500, ..., 7000 +- 500; checked here). The slice of an
# isochrone is taken at the middle of its range (13.5, 12.5, 11.75, 11.25, ..., 7.25, 6.75 ka): ice where the area deglaciated
# after it (age_deglac < t), under today's ice, and in the cells the reconstruction encloses without an age (inside its innermost
# isochrone, not today's ice).
#
# The patch. Its slices: PaleoMIST's at 15 and 5 ka (unchanged, so that the patch hands over to the base continuously at both
# ends of its window), and one at each isochrone time. Greenland is the cells nearer to Greenland's land (the 1 km land component
# with most of today's ice and the islands whose centroid lies within 120 km of it) than to any other land, within 400 km of it.
# At an isochrone time a slice is PaleoGrIS's ice fraction in Greenland, PaleoMIST's elsewhere (its mask at that time, the signed
# distance interpolated between its slices as kryos does); in a band of 2 cells either side of Greenland's boundary with other ice
# (Nares Strait, Kane Basin) the larger of the two, so that the Greenland and Innuitian ice do not part at the seam. Online only
# region 2 of the slices is used (the region map); outside it the slices carry PaleoMIST's ice, so that the signed distance of
# the patch near the region's boundary sees the ice beyond it.
#
# Output: f_ice(xc, yc, time) (time in years ago), greenland(xc, yc) (0 outside, 1 Greenland, 2 the band), region(xc, yc) (1 the
# patch's region), source(time) (Greenland from 0 PaleoMIST, 1 PaleoGrIS). The check emulates the set online (margin_set_at:
# PaleoMIST everywhere, the patch in region 2 from 15 to 5 ka), every 100 yr from 16 to 4 ka: outside region 2 it is PaleoMIST
# exactly (by construction), and it reports the PaleoGrIS ice of the slices that falls outside region 2 (not used online).
using CairoMakie, NCDatasets, Printf, Statistics, Dates

const GRID = "NH-32KM"
const FES = ENV["FESMDATA"]; const KD = ENV["KRYOS_DATA"]
const ISO = [(14000, 13000), (13000, 12000), (12000, 11500), [(t + 500, t) for t in 11000:-500:6500]...]   # isochrone ranges [yr ago]
const T_ISO = [(a + b) / 2 for (a, b) in ISO]                                                              # their slice times
const T_END = (15000.0, 5000.0)                                                                            # PaleoMIST's slices that end the window
const BAND = 2                                                                                             # cells
const REGION_ID = 2                                                                                        # Greenland, region maps v2

out = length(ARGS) >= 1 ? ARGS[1] : joinpath(ENV["KRYOS_WORK"], "margins", GRID, "$(GRID)_IceHistory-PaleoGrIS-patch.nc")
figdir = length(ARGS) >= 2 ? ARGS[2] : ""
mkpath(dirname(out))
f_pm = joinpath(FES, "North", GRID, "$(GRID)_IceHistory-PaleoMIST.nc")
f_pg = joinpath(FES, "GreenlandPaleo", "GRL-PAL-1KM", "GRL-PAL-1KM_IceHistory-PaleoGrIS.nc")
f_bm = joinpath(FES, "GreenlandPaleo", "GRL-PAL-1KM", "GRL-PAL-1KM_TOPO-BedMachine-v5.nc")
f_rg = joinpath(KD, "plastic", "$(GRID)_taub_regions_v2.nc")

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

# --- exact aggregation of the 1 km fields to the grid (same projection) ------------------------------------------------------
function wts(xt, dt, xs, ds)                    # W[i, j]: fraction of target cell i covered by source cell j, along one axis
    W = zeros(length(xt), length(xs))
    for (i, a) in enumerate(xt), (j, b) in enumerate(xs)
        o = min(a + dt / 2, b + ds / 2) - max(a - dt / 2, b - ds / 2); o > 0 && (W[i, j] = o / dt)
    end
    W
end
proj(f) = NCDataset(ds -> replace(ds["crs"].attrib["proj_params"], r"\s+" => " "), f)
proj(f_pm) == proj(f_pg) || error("$f_pm and $f_pg are not on the same projection")

# --- PaleoMIST on the grid --------------------------------------------------------------------------------------------------
ds = NCDataset(f_pm)
x, y = Float64.(ds["xc"][:]), Float64.(ds["yc"][:]); dx, dy = 1e3 * (x[2] - x[1]), 1e3 * (y[2] - y[1])
occursin("ago", ds["time"].attrib["units"]) || error("$f_pm: time is not in years ago")
t_pm = Float64.(ds["time"][:]); F_pm = Float64.(coalesce.(ds["f_ice"][:, :, :], 0.0))
close(ds)
all(t -> t in t_pm, T_END) || error("$f_pm has no slice at $(T_END)")
SD_pm = [sdist(F_pm[:, :, k] .> 0.5, dx, dy) for k in eachindex(t_pm)]
reg = NCDataset(ds -> Int.(ds["region"][:, :]), f_rg)
size(reg) == size(F_pm)[1:2] || error("$f_rg is not on $GRID")
R = reg .== REGION_ID

# --- PaleoGrIS on 1 km, its slices, aggregated ------------------------------------------------------------------------------
age, err, x1, y1 = NCDataset(f_pg) do ds
    Float64.(coalesce.(ds["age_deglac"][:, :], NaN)), Float64.(coalesce.(ds["age_deglac_err"][:, :], NaN)), Float64.(ds["xc"][:]),
    Float64.(ds["yc"][:])
end
bm = NCDataset(ds -> Int.(coalesce.(ds["mask"][:, :], 0)), f_bm)
va = .!isnan.(age); pd = (bm .== 2) .| (bm .== 3)
# the convention: each area's age is the middle of [the younger bound of its inner isochrone, the older bound of its outer one],
# its uncertainty half that span; checked on the cells within 30 yr of such a value (the rest are blends of two neighbouring
# areas by the remapping, about 10%), which must be most of them (ages are bit-rounded: 12248 for 12250)
let L = [((ISO[k + 1][2] + ISO[k][1]) / 2, (ISO[k][1] - ISO[k + 1][2]) / 2) for k in 1:length(ISO) - 1], n = 0, ok = 0
    for c in findall(va)
        n += 1; k = findfirst(l -> abs(age[c] - l[1]) < 30, L)
        k !== nothing && abs(err[c] - L[k][2]) < 5 && (ok += 1)
    end
    ok / n > 0.8 || error("PaleoGrIS: only $(round(100ok / n))% of the cells have the age and uncertainty of an area between two " *
                          "isochrones of $ISO: the convention of this tool does not hold")
    println("PaleoGrIS convention: $(round(100ok / n, digits = 1))% of the cells match an area between two isochrones")
end
# cells without an age, not today's ice, that the reconstruction encloses (not reachable from the edge of the 1 km domain)
function enclosed(free)
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
encl = enclosed(.!va .& .!pd)
Wx = wts(x, dx / 1e3, x1, 1.0); Wy = wts(y, dy / 1e3, y1, 1.0)
agg(F) = Wx * F * Wy'
cov = agg(ones(size(age)))
pg1(t) = (va .& (age .< t)) .| (.!va .& (pd .| encl))
F_pg = cat([agg(Float64.(pg1(t))) for t in T_ISO]...; dims = 3)
SD_pg = [sdist(F_pg[:, :, k] .> 0.5, dx, dy) for k in eachindex(T_ISO)]

# --- Greenland on the grid ----------------------------------------------------------------------------------------------------
function components(land)
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
B = falses(size(G))                              # the band: a cell with a cell of the other side within BAND cells
for c in CartesianIndices(G)
    cov[c] > 0.99 || continue
    for di in -BAND:BAND, dj in -BAND:BAND
        a, b = c[1] + di, c[2] + dj
        (1 <= a <= size(G, 1) && 1 <= b <= size(G, 2)) || continue
        G[a, b] != G[c] && (B[c] = true; break)
    end
end
# the band only where it separates Greenland from other ice (the seam), not along the open-ocean edge of the 400 km
fp_pm = dropdims(any(F_pm .> 0.5, dims = 3), dims = 3); fp_pg = dropdims(any(F_pg .> 0.5, dims = 3), dims = 3)
B .&= (fp_pm .| fp_pg)

# --- the patch's slices -------------------------------------------------------------------------------------------------------
t_p = sort([T_END...; T_ISO], rev = true)
fpm(t) = (k = findfirst(==(t), t_pm); k === nothing ? Float64.(mask_at(t_pm, SD_pm, t)) : F_pm[:, :, k])
F = zeros(length(x), length(y), length(t_p)); src = zeros(Int8, length(t_p))
for (k, t) in enumerate(t_p)
    a = fpm(t)
    j = findfirst(==(t), T_ISO)
    if j !== nothing
        b = F_pg[:, :, j]; a = ifelse.(G .& .!B, b, ifelse.(B, max.(a, b), a)); src[k] = 1
    end
    F[:, :, k] = a
end
gmap = Int8.(ifelse.(B, 2, ifelse.(G, 1, 0)))

rm(out; force = true)
NCDataset(out, "c") do o
    NCDataset(f_pm) do s
        defVar(o, "xc", x, ("xc",); attrib = Dict(s["xc"].attrib)); defVar(o, "yc", y, ("yc",); attrib = Dict(s["yc"].attrib))
        defVar(o, "crs", Int32(0), (); attrib = Dict(s["crs"].attrib))
    end
    defVar(o, "time", t_p, ("time",); attrib = ["units" => "years ago", "long_name" => "time (years before 1950 CE)"])
    defVar(o, "f_ice", Float32.(F), ("xc", "yc", "time"); deflatelevel = 1, attrib = [
        "units" => "1", "long_name" => "ice area fraction", "grid_mapping" => "crs",
        "comment" => "15 and 5 ka: PaleoMIST's f_ice; the isochrone times: PaleoGrIS's ice fraction (1 km, aggregated) in " *
                     "Greenland, PaleoMIST's sdist-interpolated mask (0 or 1) elsewhere, the larger of both in the band"])
    defVar(o, "greenland", gmap, ("xc", "yc"); attrib = ["long_name" => "Greenland of the patch's slices", "flag_values" => Int8[0, 1, 2],
        "flag_meanings" => "outside greenland seam_band_max_of_both", "grid_mapping" => "crs"])
    defVar(o, "region", Int8.(R), ("xc", "yc"); attrib = ["long_name" => "the patch's region (region id $REGION_ID of " *
        "$(basename(f_rg))): online the slices are used here only", "grid_mapping" => "crs"])
    defVar(o, "source", src, ("time",); attrib = ["long_name" => "source of Greenland's slice", "flag_values" => Int8[0, 1],
        "flag_meanings" => "PaleoMIST PaleoGrIS"])
    o.attrib["title"] = "Margin patch on $GRID: Greenland from PaleoGrIS 1.0 (13.5-6.75 ka), PaleoMIST 1.0 at 15 and 5 ka"
    o.attrib["usage"] = "kryos plastic backend, ice instance on $GRID with PaleoMIST's margin_file: margin_patches = 'margin_grl'; " *
                        "&margin_grl file = '$out', var = 'f_ice', H_min = 0.5, region_file = " *
                        "'kryos_data/plastic/{grid_name}_taub_regions_v2.nc', region_var = 'region', region_ids = '$REGION_ID'"
    o.attrib["comment"] = "PaleoGrIS slices at the middle of each isochrone's range ($(join(Int.(T_ISO), ", ")) years ago): ice " *
                          "where age_deglac < t, under today's ice (BedMachine v5) and in the cells the reconstruction encloses " *
                          "without an age. Greenland: nearer to Greenland's land than to any other, within 400 km; band of " *
                          "$BAND cells at its boundary with other ice: max of both. Made by tools/levante/margin_patch.jl (kryos)."
    o.attrib["references"] = "Gowan et al. 2021, Nat. Commun. 12, 1199, doi:10.1038/s41467-021-21469-w; Leger et al. 2024, " *
                             "Clim. Past 20, 701-755, doi:10.5194/cp-20-701-2024"
    o.attrib["sources"] = join([f_pm, f_pg, f_bm, f_rg], ", ")
    o.attrib["history"] = "$(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")): margin_patch.jl"
end
println("wrote $out: $(length(t_p)) slices, Greenland $(count(G)) cells, band $(count(B)) cells, region $(count(R)) cells")

# --- check: the set online (margin_set_at), areas, maps ------------------------------------------------------------------------
SD_p = [sdist(F[:, :, k] .> 0.5, dx, dy) for k in eachindex(t_p)]
set_at(t) = (m = mask_at(t_pm, SD_pm, t); T_END[2] <= t <= T_END[1] ? ifelse.(R, mask_at(t_p, SD_p, t), m) : m)
dA = dx * dy / 1e9                               # 10^3 km2
tt = collect(16000.0:-100.0:4000.0)
win(t) = T_ISO[end] <= t <= T_ISO[1]
A_pm = [count(mask_at(t_pm, SD_pm, t) .& R) * dA for t in tt]; A_s = [count(set_at(t) .& R) * dA for t in tt]
A_pg = [win(t) ? count(mask_at(T_ISO, SD_pg, t) .& R) * dA : NaN for t in tt]
d_out = [count(.!R .& (set_at(t) .!= mask_at(t_pm, SD_pm, t))) for t in tt]
fp_p = dropdims(any(F .> 0.5, dims = 3), dims = 3)
io = open(replace(out, r"\.nc$" => "_check.md"), "w")
println(io, "# $(basename(out))\n\nslices: $(join(Int.(t_p), ", ")) years ago (Greenland from PaleoGrIS at ",
        join(Int.(t_p[src .== 1]), ", "), ")\n")
println(io, "region $REGION_ID (the patch's): $(count(R)) cells; Greenland $(count(G)) cells ($(count(G .& R)) in the region), band ",
        "$(count(B)) cells ($(count(B .& R)) in the region); PaleoGrIS: enclosed cells without an age $(count(encl)) km2 (1 km)\n")
println(io, "PaleoGrIS ice of the slices outside the region (not used online): ", count(fp_pg .& G .& .!R), " cells (",
        join(["region $r: $(count(fp_pg .& G .& .!R .& (reg .== r)))" for r in sort(unique(reg[fp_pg .& G .& .!R]))], ", "), ")\n")
println(io, "patch footprint in the region: $(count(fp_p .& R)) cells, PaleoMIST's $(count(fp_pm .& R)); new: ",
        "$(count(fp_p .& R .& .!fp_pm))\n")
println(io, "the set (PaleoMIST everywhere, the patch in the region 15-5 ka) against PaleoMIST outside the region, 16-4 ka every ",
        "100 yr: $(sum(d_out)) cells differ (by construction 0)\n")
println(io, "| t ka | PaleoMIST | PaleoGrIS | the set |\n|---|---|---|---|")
for (i, t) in enumerate(tt)
    mod(t, 500) == 0 && println(io, @sprintf("| %.1f | %.0f | %.0f | %.0f |", t / 1000, A_pm[i], A_pg[i], A_s[i]))
end
close(io)
println(read(replace(out, r"\.nc$" => "_check.md"), String))
sum(d_out) == 0 || error("the set differs from PaleoMIST outside region $REGION_ID")

fig = Figure(size = (1800, 1050))
ax = Axis(fig[1, 1:3], xreversed = true, xlabel = "ka", ylabel = "10³ km²", title = "ice area in region $REGION_ID (kryos's interpolation, 100 yr)",
          limits = ((4, 16), nothing))
lines!(ax, tt ./ 1000, A_pm, color = :black, linewidth = 2, label = "PaleoMIST (base)")
lines!(ax, tt ./ 1000, A_pg, color = :red, linewidth = 2, label = "PaleoGrIS")
lines!(ax, tt ./ 1000, A_s, color = :royalblue, linestyle = :dash, linewidth = 2.5, label = "base + patch")
vlines!(ax, t_p ./ 1000, color = (:red, 0.25)); vlines!(ax, t_pm ./ 1000, color = (:black, 0.25), linestyle = :dot)
axislegend(ax, position = :rt)
ax = Axis(fig[1, 4:5], title = "greenland (1), band (2); region $REGION_ID (green), patch footprint (black)", aspect = DataAspect(),
          limits = ((-900, 1000), (-3500, -500)))
heatmap!(ax, x, y, ifelse.(fp_p .| (gmap .> 0), Float64.(gmap), NaN), colormap = [:gray80, :royalblue, :orange], colorrange = (0, 2))
contour!(ax, x, y, Float64.(fp_p .& R), levels = [0.5], color = :black)
contour!(ax, x, y, Float64.(R), levels = [0.5], color = :green, linewidth = 1.5)
for (p, t) in enumerate([14000.0, 12000.0, 10000.0, 8000.0, 7000.0])
    m = mask_at(t_pm, SD_pm, t); g = set_at(t)
    c = fill(NaN, size(m)); c[m .& g] .= 1; c[m .& .!g] .= 2; c[g .& .!m] .= 3
    local ax = Axis(fig[2, p], aspect = DataAspect(), title = "$(t / 1000) ka", limits = ((-900, 1000), (-3500, -500)),
              xticklabelsvisible = false, yticklabelsvisible = false)
    heatmap!(ax, x, y, c, colormap = cgrad([:gray75, :royalblue, :darkorange], 3, categorical = true), colorrange = (1, 3))
    contour!(ax, x, y, Float64.(R), levels = [0.5], color = :green, linewidth = 0.8)
end
Legend(fig[3, 1:5], [PolyElement(color = :gray75), PolyElement(color = :royalblue), PolyElement(color = :darkorange), LineElement(color = :green)],
       ["both", "PaleoMIST only", "base + patch only", "region $REGION_ID"], orientation = :horizontal, framevisible = false)
png = replace(out, r"\.nc$" => "_check.png")
save(png, fig)
isempty(figdir) || (mkpath(figdir); cp(png, joinpath(figdir, basename(png)); force = true))
