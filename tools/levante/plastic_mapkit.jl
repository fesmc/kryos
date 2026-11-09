# Global maps of kryos plastic runs (included as a module by plastic_lgm_map.jl and plastic_ensemble_sd_map.jl): the Equal Earth
# projection (Savric et al. 2018, written out here: no GeoMakie), fields on structured meshes drawn as triangles, contours, the VILMA2
# fields of a run at an age (geo_glob_vilma.nc), the 2-D fields of its plastic instances (<instance>.nc with the grid's lon2D, lat2D),
# coastlines, ice margins, and the GAPSLIP 2.1 index points scored by pm_compare.jl.
using CairoMakie, NCDatasets, Printf, Statistics

module PM
include("pm_compare.jl")
end
module KO
include("pm_kryos_out.jl")
end

const FES = "/work/ba1442/fesmdata"
const INSTANCES = [("ice_nor", "North", "NH-32KM"), ("ice_ant", "Antarctica", "ANT-16KM"), ("ice_pat", "PatagoniaPaleo", "PAT-PAL-8KM")]
const SIGMA_M = 5.0                                   # [m] the model error of z_m (as plastic_esmda_report.jl)
const NEAR = ["Laurentide", "Greenland", "Fennoscandia/Barents", "Antarctica"]
const LAND = RGBf(0.88, 0.86, 0.82)
const MARGIN = RGBf(0.0, 0.39, 0.0)                  # dark green: the target margin (PaleoMIST's, prescribed)
const ICEMAP = cgrad([RGBf(0.42, 0.42, 0.44), RGBf(1, 1, 1)])
const RSLMAP = cgrad([RGBf(0.52, 0.65, 0.80), RGBf(0.93, 0.96, 0.98)])   # toned down: light, low saturation
const RSL_LOW, RSL_HIGH = RGBf(0.42, 0.55, 0.71), RGBf(0.97, 0.98, 0.99)
kyr(age) = @sprintf("%g kyr", age / 1000)

# --- projection ---------------------------------------------------------------------------------------------------------------------
const LON0 = Ref(0.0)
const EA1, EA2, EA3, EA4 = 1.340264, -0.081106, 0.000893, 0.003796
relon(lon) = mod(lon - LON0[] + 180, 360) - 180       # longitude relative to the central meridian, [-180, 180)
function eqe(dl, lat)                                 # dl: relative longitude [deg], not wrapped
    θ = asin(sqrt(3) / 2 * sind(lat)); λ = deg2rad(dl)
    Point2f(2sqrt(3) * λ * cos(θ) / (3 * (9EA4 * θ^8 + 7EA3 * θ^6 + 3EA2 * θ^2 + EA1)), EA4 * θ^9 + EA3 * θ^7 + EA2 * θ^3 + EA1 * θ)
end
proj(lon, lat) = eqe(relon(lon), lat)
const SEAM = 1.0                                      # [map units] a step wider than this crosses the seam

# a field at the vertices of a structured mesh as triangles: quads with four finite corners where keep, not across the seam
function meshfield!(ax, P, V, keep; kw...)
    li = LinearIndices(V); F = Makie.GeometryBasics.GLTriangleFace[]
    for j in 1:size(V, 2)-1, i in 1:size(V, 1)-1
        c = (li[i, j], li[i+1, j], li[i+1, j+1], li[i, j+1])
        all(k -> keep[k] && isfinite(V[k]), c) || continue
        xs = [P[k][1] for k in c]; maximum(xs) - minimum(xs) < SEAM || continue
        push!(F, Makie.GeometryBasics.GLTriangleFace(c[1], c[2], c[3]), Makie.GeometryBasics.GLTriangleFace(c[1], c[3], c[4]))
    end
    isempty(F) && return nothing
    mesh!(ax, Makie.GeometryBasics.Mesh(vec(P), F); color = Float32.(vec(replace(V, NaN => 0.0))), shading = NoShading, kw...)
end
# polylines in map units, broken (NaN) where a step crosses the seam
function seamlines(pts)
    out = Point2f[]
    for (k, p) in enumerate(pts)
        k > 1 && abs(p[1] - pts[k-1][1]) > SEAM && push!(out, Point2f(NaN, NaN))
        push!(out, p)
    end
    push!(out, Point2f(NaN, NaN)); out
end
function split_nan(pts)
    r = Vector{Vector{Point2f}}(); c = Point2f[]
    for p in pts
        if isnan(p[1]); isempty(c) || push!(r, c); c = Point2f[]; else push!(c, p); end
    end
    isempty(c) || push!(r, c); r
end
# contour lines of Z (Z[i, j] at x[i], y[j]) at level, as (xs, ys) vertex vectors
clines(x, y, Z, level) = [Makie.Contours.coordinates(l) for l in Makie.Contours.lines(Makie.Contours.contour(x, y, Z, level))]

# --- the VILMA2 fields of a run at an age (Gauss grid, columns by relative longitude, padded to -180 and +180) ------------------------
function read_global(run, age)
    NCDataset(joinpath(run, "geo_glob_vilma.nc")) do gv
        tg = Float64.(gv["time"][:]); kt = findfirst(t -> abs(t + age) < 1, tg); k0 = argmin(abs.(tg))
        kt === nothing && error("$(run): no record at $(age) years ago in geo_glob_vilma.nc")
        glon = Float64.(gv["lon"][:]); glat = Float64.(gv["lat"][:])
        dl = relon.(glon); p = sortperm(dl); q = sortperm(glat)
        pad(F) = (G = F[p, q]; [G[1:1, :]; G; G[end:end, :]])
        f(v, k) = Float64.(gv[v][:, :, k])
        dlp = [-180.0; dl[p]; 180.0]; latq = glat[q]
        LON = [d for d in dlp, _ in latq]; LAT = [l for _ in dlp, l in latq]
        (dl = dlp, lat = latq, LON = LON, LAT = LAT, P = [eqe(LON[i], LAT[i]) for i in eachindex(LON)],
         rsl = pad(f("rsl", kt) .- f("rsl", k0)), ocn = pad(f("C_ocean", kt)), ocn0 = pad(f("C_ocean", k0)),
         h = pad(f("h_ice", kt)), zb = pad(f("z_bed", kt)), esl = Float64(gv["bsl_ive"][kt] - gv["bsl_ive"][k0]))
    end
end

# --- the plastic instances of a run at an age ---------------------------------------------------------------------------------------
const GRIDLL = Dict{String,Any}()
gridll(dom, grid) = get!(GRIDLL, grid) do
    NCDataset(joinpath(FES, dom, grid, "$(grid)_grid.nc")) do ds
        lo, la = Float64.(ds["lon2D"][:, :]), Float64.(ds["lat2D"][:, :])
        (lon = lo, lat = la, P = [proj(lo[i], la[i]) for i in eachindex(lo)])
    end
end
function read_regional(run, age; need = true)
    out = []
    for (inst, dom, grid) in INSTANCES
        f = joinpath(run, "$(inst).nc"); isfile(f) || continue
        r = NCDataset(f) do ds
            t = Float64.(ds["time"][:]); k = findfirst(x -> abs(x + age) < 1, t)
            k === nothing && (need ? error("$(f): no record at $(age) years ago (2-D output?)") : return nothing)
            g(v) = Float64.(nomissing(ds[v][:, :, k], NaN))
            (H = g("H_ice"), zs = g("z_srf") .- g("z_sl"))
        end
        r === nothing && continue
        ll = gridll(dom, grid)
        push!(out, (inst = inst, grid = grid, lon = ll.lon, lat = ll.lat, P = ll.P, H = r.H, zs = r.zs))
    end
    out
end
# the cells of the global grid inside the area the plastic grids cover with ice (0.5 degree bins)
function covered(G, regional)
    cov = falses(720, 360)
    binof(lon, lat) = (clamp(floor(Int, mod(lon, 360) / 0.5) + 1, 1, 720), clamp(floor(Int, (lat + 90) / 0.5) + 1, 1, 360))
    for r in regional, i in eachindex(r.H)
        r.H[i] > 1 || continue
        bi, bj = binof(r.lon[i], r.lat[i]); nl = ceil(Int, 0.6 / max(cosd(r.lat[i]), 0.05))
        for a in -nl:nl, b in -1:1; cov[mod1(bi + a, 720), clamp(bj + b, 1, 360)] = true; end
    end
    [cov[binof(G.LON[i, j] + LON0[], G.LAT[i, j])...] for i in axes(G.LON, 1), j in axes(G.LON, 2)]
end

# --- drawing ------------------------------------------------------------------------------------------------------------------------
const OUTLINE = [[eqe(-180.0, l) for l in -90:1.0:90]; [eqe(180.0, l) for l in 90:-1.0:-90]]
function mapaxis(pos; kw...)
    ax = Axis(pos; aspect = DataAspect(), kw...); hidedecorations!(ax); hidespines!(ax)
    poly!(ax, OUTLINE, color = LAND); ax
end
function graticule!(ax)
    for l in -60:30:60; lines!(ax, [eqe(d, l) for d in -180:2.0:180], color = (:black, 0.10), linewidth = 0.7); end
    for d in -180:30:180; lines!(ax, [eqe(Float64(d), l) for l in -90:2.0:90], color = (:black, 0.10), linewidth = 0.7); end
end
frame!(ax) = lines!(ax, [OUTLINE; OUTLINE[1:1]], color = :black, linewidth = 1.2)
# the ice surface elevation: the global load outside the plastic grids, then the plastic grids
function ice!(ax, G, regional)
    cov = covered(G, regional)
    m = meshfield!(ax, G.P, G.zb .+ G.h, (G.h .> 1) .& .!cov, colormap = ICEMAP, colorrange = (0, 4000))
    for r in regional
        mm = meshfield!(ax, r.P, r.zs, r.H .> 1, colormap = ICEMAP, colorrange = (0, 4000)); mm === nothing || (m = mm)
    end
    m
end
# contours of a global field over the ocean of the age, labelled (one label per level, between about 35 S and 35 N, apart)
function ocean_contours!(ax, G, F, levels; labelsize = 22, color = (:black, 0.45))
    labelled = []
    for lv in levels, (xs, ys) in clines(G.dl, G.lat, F, lv)
        pts = Point2f[]; seg = Point2f[]
        for (x, y) in zip(xs, ys)
            i = clamp(searchsortedfirst(G.dl, x), 1, length(G.dl)); j = clamp(searchsortedfirst(G.lat, y), 1, length(G.lat))
            if G.ocn[i, j] > 0.5
                push!(seg, eqe(x, y))
            elseif !isempty(seg)
                append!(pts, seamlines(seg)); seg = Point2f[]
            end
        end
        isempty(seg) || append!(pts, seamlines(seg)); isempty(pts) && continue
        lines!(ax, pts, color = color, linewidth = 1.0); push!(labelled, (lv, split_nan(pts)))
    end
    placed = Point2f[]
    for lv in levels
        cand = sort([p for (l, ps) in labelled if l == lv for p in ps if length(p) > 60], by = p -> -length(p))
        for p in cand, pt in (p[div(length(p), 2)], p[div(length(p), 4)], p[div(3length(p), 4)])
            abs(pt[2]) < 0.55 && all(q -> hypot((q - pt)...) > 0.5, placed) || continue
            text!(ax, pt, text = string(round(Int, lv)), fontsize = labelsize, align = (:center, :center), color = :gray10, font = :bold)
            push!(placed, pt); break
        end
    end
end
function coasts!(ax, G)
    for (C, col, lw) in ((G.ocn0, (:gray35, 0.8), 0.8), (G.ocn, :black, 1.3))
        for (xs, ys) in clines(G.dl, G.lat, C, 0.5)
            lines!(ax, seamlines([eqe(x, y) for (x, y) in zip(xs, ys)]), color = col, linewidth = lw)
        end
    end
end
# a point at fractional indices of a regional grid: its lon-lat (bilinear on the unit sphere, across the date line)
function grid_lonlat(r, fi, fj)
    i0 = clamp(floor(Int, fi), 1, size(r.lon, 1) - 1); j0 = clamp(floor(Int, fj), 1, size(r.lon, 2) - 1); a = fi - i0; b = fj - j0
    v = [0.0, 0.0, 0.0]
    for (ii, jj, w) in ((i0, j0, (1 - a) * (1 - b)), (i0 + 1, j0, a * (1 - b)), (i0, j0 + 1, (1 - a) * b), (i0 + 1, j0 + 1, a * b))
        lo, la = r.lon[ii, jj], r.lat[ii, jj]; v .+= w .* [cosd(la) * cosd(lo), cosd(la) * sind(lo), sind(la)]
    end
    atand(v[2], v[1]), atand(v[3], hypot(v[1], v[2]))
end
# the ice margin (H = 1 m) of the plastic grids
function margins!(ax, regional; color = MARGIN, linewidth = 1.8)
    for r in regional, (xs, ys) in clines(1.0:size(r.H, 1), 1.0:size(r.H, 2), Float64.(r.H .> 1), 0.5)
        lines!(ax, seamlines([proj(grid_lonlat(r, x, y)...) for (x, y) in zip(xs, ys)]), color = color, linewidth = linewidth)
    end
end

# --- the index points ---------------------------------------------------------------------------------------------------------------
const OBS = Ref{Any}(nothing)
obs() = (OBS[] === nothing && (OBS[] = PM.read_indicators(FES)); OBS[])
const AGES = Ref{Any}(nothing)
ages() = (AGES[] === nothing && (AGES[] = sort(NCDataset(ds -> Float64.(ds["time"][:]),
                                     joinpath(FES, "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc")))); AGES[])
zm(s) = isfinite(s.z) ? (s.z == 0 ? 0.0 : s.d / sqrt((s.d / s.z)^2 + SIGMA_M^2)) : NaN
sphere_mean(lo, la) = (v = sum([cosd(b) * cosd(a), cosd(b) * sind(a), sind(b)] for (a, b) in zip(lo, la));
                       (atand(v[2], v[1]), atand(v[3], hypot(v[1], v[2]))))
# the index points with an age in window (lo, hi) [yr]: per site, position and n; with a run, the RMS of z_m of those points
function sites(window; run = nothing)
    O = obs(); ip = [i for (i, o) in enumerate(O) if o.type == 0 && window[1] <= o.t <= window[2]]
    sc = nothing
    if run !== nothing
        isfile(joinpath(run, "vilma_out.nc")) || KO.kryos_vilma_out(run, run)
        S = PM.series_run(PM.read_run(run), [O[i].lon for i in ip], [O[i].lat for i in ip], ages())
        sc = Dict(i => PM.score(O[i], ages(), S[k, :]) for (k, i) in enumerate(ip))
        ip = [i for i in ip if isfinite(sc[i].z)]
    end
    names = unique(O[i].site_name for i in ip)
    idx = Dict(s => [i for i in ip if O[i].site_name == s] for s in names)
    map(names) do s
        lo, la = sphere_mean([O[i].lon for i in idx[s]], [O[i].lat for i in idx[s]])
        (name = s, pos = proj(lo, la), near = O[idx[s][1]].group in NEAR, n = length(idx[s]),
         rms = sc === nothing ? NaN : sqrt(mean(zm(sc[i])^2 for i in idx[s])))
    end, ip, sc
end
const SITEMAP = :YlOrRd
const SITE_MAX = 4.0
function sites!(ax, ss; color = nothing, markersize = 15)
    for (near, mk) in ((true, :utriangle), (false, :circle))
        k = [s for s in ss if s.near == near]; isempty(k) && continue
        c = color === nothing ? (colormap = SITEMAP, color = min.([s.rms for s in k], SITE_MAX), colorrange = (0, SITE_MAX)) : (color = color,)
        scatter!(ax, [s.pos for s in k]; marker = mk, markersize = markersize, strokewidth = 0.9, strokecolor = :black, c...)
    end
end
