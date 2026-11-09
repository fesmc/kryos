# Map of a kryos plastic run at one age, after Gowan et al. (2021, Fig. 2) (docs/dev/plastic-backend.md): the ice surface elevation
# above the sea level of that age (the plastic instances on their grids, the rest of the world's ice from VILMA2's load), the relative
# sea level of that age against today over the ocean of that age (filled and contoured), the coastlines of that age and of today, the
# ice margins, and the GAPSLIP 2.1 sites of index points coloured by the run's misfit there.
#   julia --project=<julia-plot> tools/levante/plastic_lgm_map.jl <run folder> <out.png> [age=20000] [label=<text>] [lon0=0]
# Equal Earth projection (Savric et al. 2018) written out here (no GeoMakie). Surface elevation of the plastic ice: z_srf - z_sl of
# <instance>.nc (needs its 2-D output) where H_ice > 1 m; elsewhere VILMA2's z_bed + h_ice of geo_glob_vilma.nc (bed relative to the
# sea level of the age) where h_ice > 1 m. RSL: rsl(age) - rsl(0 ka) of geo_glob_vilma.nc (as the scores), ocean where C_ocean > 0.5.
# Sites: RMS of z_m = d / sqrt(sigma^2 + sigma_m^2) (sigma_m 5 m) over the site's index points of all ages; triangles the ice-sheet
# groups (Laurentide, Greenland, Fennoscandia/Barents, Antarctica), circles the far field.
using CairoMakie, NCDatasets, Printf, Statistics

module PM
include("pm_compare.jl")
end
module KO
include("pm_kryos_out.jl")
end

const FES = "/work/ba1442/fesmdata"
const INSTANCES = [("ice_nor", "North", "NH-32KM"), ("ice_ant", "Antarctica", "ANT-16KM"), ("ice_pat", "PatagoniaPaleo", "PAT-PAL-8KM")]
const SIGMA_M = 5.0
const NEAR = ["Laurentide", "Greenland", "Fennoscandia/Barents", "Antarctica"]

run, outpng = abspath(ARGS[1]), ARGS[2]
opt = Dict(split(a, "=", limit = 2) for a in ARGS[3:end])
const AGE = parse(Float64, get(opt, "age", "20000"))
const LON0 = parse(Float64, get(opt, "lon0", "0"))
label = get(opt, "label", basename(run))

# --- projection ---------------------------------------------------------------------------------------------------------------------
const EA1, EA2, EA3, EA4 = 1.340264, -0.081106, 0.000893, 0.003796
relon(lon) = mod(lon - LON0 + 180, 360) - 180          # longitude relative to the central meridian, [-180, 180)
function eqe(dl, lat)                                  # dl: relative longitude [deg], not wrapped
    θ = asin(sqrt(3) / 2 * sind(lat)); λ = deg2rad(dl)
    Point2f(2sqrt(3) * λ * cos(θ) / (3 * (9EA4 * θ^8 + 7EA3 * θ^6 + 3EA2 * θ^2 + EA1)), EA4 * θ^9 + EA3 * θ^7 + EA2 * θ^3 + EA1 * θ)
end
proj(lon, lat) = eqe(relon(lon), lat)
const SEAM = 1.0                                       # [map units] a step wider than this crosses the seam

# a field at the vertices of a structured mesh (2-D lon, lat) as triangles: quads with four finite corners where keep, not across the seam
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

# polylines in map units from points, broken (NaN) where a step crosses the seam
function seamlines(pts)
    out = Point2f[]
    for (k, p) in enumerate(pts)
        k > 1 && abs(p[1] - pts[k-1][1]) > SEAM && push!(out, Point2f(NaN, NaN))
        push!(out, p)
    end
    push!(out, Point2f(NaN, NaN)); out
end
# the pieces of a NaN-separated polyline
function split_nan(pts)
    r = Vector{Vector{Point2f}}(); c = Point2f[]
    for p in pts
        if isnan(p[1]); isempty(c) || push!(r, c); c = Point2f[]; else push!(c, p); end
    end
    isempty(c) || push!(r, c); r
end
# contour lines of Z (Z[i, j] at x[i], y[j]) at level, as vectors of (x, y) vertices
clines(x, y, Z, level) = [Makie.Contours.coordinates(l) for l in Makie.Contours.lines(Makie.Contours.contour(x, y, Z, level))]

# --- the global fields (VILMA2) -----------------------------------------------------------------------------------------------------
gv = NCDataset(joinpath(run, "geo_glob_vilma.nc"))
tg = Float64.(gv["time"][:]); kt = findfirst(t -> abs(t + AGE) < 1, tg); k0 = argmin(abs.(tg))
kt === nothing && error("no record at $(AGE) years ago in geo_glob_vilma.nc")
glon = Float64.(gv["lon"][:]); glat = Float64.(gv["lat"][:])
rsl = Float64.(gv["rsl"][:, :, kt]) .- Float64.(gv["rsl"][:, :, k0])
ocn = Float64.(gv["C_ocean"][:, :, kt]); ocn0 = Float64.(gv["C_ocean"][:, :, k0])
hg = Float64.(gv["h_ice"][:, :, kt]); zbg = Float64.(gv["z_bed"][:, :, kt])
esl_age = Float64(gv["bsl_ive"][kt] - gv["bsl_ive"][k0])
close(gv)
# columns sorted by relative longitude, padded to -180 and +180 (the values of the nearest column), latitude ascending
dl = relon.(glon); p = sortperm(dl); q = sortperm(glat)
padl(F) = (G = F[p, q]; [G[1:1, :]; G; G[end:end, :]])
dlp = [-180.0; dl[p]; 180.0]; latq = glat[q]
GLON = [d for d in dlp, _ in latq]; GLAT = [l for _ in dlp, l in latq]
GP = [eqe(GLON[i], GLAT[i]) for i in eachindex(GLON)]
RSL, OCN, OCN0, HG, ZG = padl(rsl), padl(ocn), padl(ocn0), padl(hg), padl(zbg)

# --- the plastic instances ------------------------------------------------------------------------------------------------------------
regional = []
for (inst, dom, grid) in INSTANCES
    f = joinpath(run, "$(inst).nc"); isfile(f) || continue
    lon2, lat2 = NCDataset(joinpath(FES, dom, grid, "$(grid)_grid.nc")) do ds
        Float64.(ds["lon2D"][:, :]), Float64.(ds["lat2D"][:, :])
    end
    H, zs, zsl = NCDataset(f) do ds
        t = Float64.(ds["time"][:]); k = findfirst(x -> abs(x + AGE) < 1, t)
        k === nothing && error("$(inst).nc: no record at $(AGE) years ago (2-D output?)")
        g(v) = Float64.(nomissing(ds[v][:, :, k], NaN))
        g("H_ice"), g("z_srf"), g("z_sl")
    end
    push!(regional, (inst = inst, lon = lon2, lat = lat2, H = H, zs = zs .- zsl, P = [proj(lon2[i], lat2[i]) for i in eachindex(lon2)]))
end
isempty(regional) && @warn "no 2-D output of a plastic instance in $run: the ice is VILMA2's load only"
# the area the plastic grids cover, on 0.5 degree bins (the global ice is drawn outside it only)
covered = falses(720, 360)
binof(lon, lat) = (clamp(floor(Int, mod(lon, 360) / 0.5) + 1, 1, 720), clamp(floor(Int, (lat + 90) / 0.5) + 1, 1, 360))
for r in regional, i in eachindex(r.H)
    r.H[i] > 1 || continue
    bi, bj = binof(r.lon[i], r.lat[i]); nl = ceil(Int, 0.6 / max(cosd(r.lat[i]), 0.05))
    for a in -nl:nl, b in -1:1; covered[mod1(bi + a, 720), clamp(bj + b, 1, 360)] = true; end
end
gcov = [covered[binof(GLON[i, j] + LON0, GLAT[i, j])...] for i in axes(GLON, 1), j in axes(GLON, 2)]

# --- scores per site ----------------------------------------------------------------------------------------------------------------
isfile(joinpath(run, "vilma_out.nc")) || KO.kryos_vilma_out(run, run)
obs = PM.read_indicators(FES)
A = sort(NCDataset(ds -> Float64.(ds["time"][:]), joinpath(FES, "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc")))
S = PM.series_run(PM.read_run(run), [o.lon for o in obs], [o.lat for o in obs], A)
sc = [PM.score(o, A, S[i, :]) for (i, o) in enumerate(obs)]
zm(s) = isfinite(s.z) ? (s.z == 0 ? 0.0 : s.d / sqrt((s.d / s.z)^2 + SIGMA_M^2)) : NaN
ip = [i for (i, o) in enumerate(obs) if o.type == 0 && isfinite(sc[i].z)]
sites = unique(obs[i].site_name for i in ip)
sidx = Dict(s => [i for i in ip if obs[i].site_name == s] for s in sites)
srms = [sqrt(mean(zm(sc[i])^2 for i in sidx[s])) for s in sites]
spos = [proj(mean(obs[i].lon for i in sidx[s]), mean(obs[i].lat for i in sidx[s])) for s in sites]
snear = [obs[sidx[s][1]].group in NEAR for s in sites]
sn = [length(sidx[s]) for s in sites]
rmsset(f) = (k = [i for i in ip if f(obs[i])]; sqrt(mean(zm(sc[i])^2 for i in k)))
rms_near = rmsset(o -> o.group in NEAR); rms_far = rmsset(o -> !(o.group in NEAR))

# --- figure ------------------------------------------------------------------------------------------------------------------------
set_theme!(fontsize = 22)
fig = Figure(size = (2000, 1250))
ax = Axis(fig[1, 1:3], aspect = DataAspect(), title = @sprintf("%s: %.1f ka", label, AGE / 1000), titlesize = 32,
          subtitle = @sprintf("ESL %.1f m;  RMS z_m: ice sheets %.2f, far field %.2f (%d sites, %d index points)", esl_age, rms_near,
                              rms_far, length(sites), length(ip)), subtitlesize = 22)
hidedecorations!(ax); hidespines!(ax)
outline = [[eqe(-180.0, l) for l in -90:1.0:90]; [eqe(180.0, l) for l in 90:-1.0:-90]]
poly!(ax, outline, color = RGBf(0.86, 0.84, 0.79))                              # land of the age
# the colour range of the ocean: the far field (|lat| < 45), 10 m steps; beyond it clipped
ocmin, ocmax = let v = RSL[(OCN .> 0.5) .& (abs.(GLAT) .< 45)]; floor(quantile(v, 0.01) / 10) * 10, ceil(quantile(v, 0.99) / 10) * 10 end
const OCMAP = Reverse(:Blues)
ocm = meshfield!(ax, GP, RSL, OCN .> 0.5, colormap = OCMAP, colorrange = (ocmin, ocmax), lowclip = :midnightblue,
                 highclip = RGBf(0.93, 0.96, 1.0))
const ICEMAP = cgrad([RGBf(0.38, 0.38, 0.40), RGBf(1, 1, 1)])
# graticule
for l in -60:30:60; lines!(ax, [eqe(d, l) for d in -180:2.0:180], color = (:black, 0.12), linewidth = 0.7); end
for d in -180:30:180; lines!(ax, [eqe(Float64(d), l) for l in -90:2.0:90], color = (:black, 0.12), linewidth = 0.7); end
# RSL contours over the ocean, labelled
levels = ocmin:10:ocmax
labelled = []
for lv in levels
    for (xs, ys) in clines(dlp, latq, RSL, lv)
        pts = Point2f[]; seg = Point2f[]
        for (x, y) in zip(xs, ys)
            i = clamp(searchsortedfirst(dlp, x), 1, length(dlp)); j = clamp(searchsortedfirst(latq, y), 1, length(latq))
            if OCN[i, j] > 0.5
                push!(seg, eqe(x, y))
            elseif !isempty(seg)
                append!(pts, seamlines(seg)); seg = Point2f[]
            end
        end
        isempty(seg) || append!(pts, seamlines(seg))
        isempty(pts) && continue
        lines!(ax, pts, color = (:black, 0.4), linewidth = 1.0)
        push!(labelled, (lv, split_nan(pts)))
    end
end
# labels: one per level, on its longest pieces between about 35 S and 35 N (the far field), apart from the others
placed = Point2f[]
for (lv, pieces) in labelled
    cand = sort([p for p in pieces if length(p) > 60], by = p -> -length(p)); n = 0
    for p in cand, pt in (p[div(length(p), 2)], p[div(length(p), 4)], p[div(3length(p), 4)])
        abs(pt[2]) < 0.55 && all(q -> hypot((q - pt)...) > 0.45, placed) || continue
        text!(ax, pt, text = string(round(Int, lv)), fontsize = 16, align = (:center, :center), color = :gray10)
        push!(placed, pt); n += 1; break
    end
end
# ice: the global load outside the plastic grids, then the plastic grids
zg = ZG .+ HG
icem = meshfield!(ax, GP, zg, (HG .> 1) .& .!gcov, colormap = ICEMAP, colorrange = (0, 4000))
for r in regional
    m = meshfield!(ax, r.P, r.zs, r.H .> 1, colormap = ICEMAP, colorrange = (0, 4000))
    m === nothing || (icem = m)
end
# coastlines: of the age (black) and today (grey)
for (C, col, lw) in ((OCN0, (:gray30, 0.8), 0.8), (OCN, :black, 1.3))
    for (xs, ys) in clines(dlp, latq, C, 0.5)
        lines!(ax, seamlines([eqe(x, y) for (x, y) in zip(xs, ys)]), color = col, linewidth = lw)
    end
end
# ice margins of the plastic grids: H = 1 m, through the grid's lon-lat (on the unit sphere, across the date line)
function grid_lonlat(r, fi, fj)
    i0 = clamp(floor(Int, fi), 1, size(r.lon, 1) - 1); j0 = clamp(floor(Int, fj), 1, size(r.lon, 2) - 1); a = fi - i0; b = fj - j0
    v = [0.0, 0.0, 0.0]
    for (ii, jj, w) in ((i0, j0, (1 - a) * (1 - b)), (i0 + 1, j0, a * (1 - b)), (i0, j0 + 1, (1 - a) * b), (i0 + 1, j0 + 1, a * b))
        lo, la = r.lon[ii, jj], r.lat[ii, jj]; v .+= w .* [cosd(la) * cosd(lo), cosd(la) * sind(lo), sind(la)]
    end
    atand(v[2], v[1]), atand(v[3], hypot(v[1], v[2]))
end
for r in regional
    for (xs, ys) in clines(1.0:size(r.H, 1), 1.0:size(r.H, 2), Float64.(r.H .> 1), 0.5)
        lines!(ax, seamlines([proj(grid_lonlat(r, x, y)...) for (x, y) in zip(xs, ys)]), color = :firebrick, linewidth = 1.6)
    end
end
lines!(ax, [outline; outline[1:1]], color = :black, linewidth = 1.2)
# sites
smax = 4.0
for (near, mk) in ((true, :utriangle), (false, :circle))
    k = findall(==(near), snear)
    scatter!(ax, spos[k], color = min.(srms[k], smax), colormap = :YlOrRd, colorrange = (0, smax), marker = mk,
             markersize = 9 .+ 3 .* sqrt.(sn[k]), strokewidth = 0.8, strokecolor = :black)
end
cbw = Relative(0.8)
Colorbar(fig[2, 1], icem, label = "ice surface elevation [m]", vertical = false, width = cbw, labelsize = 22)
Colorbar(fig[2, 2], ocm, label = "RSL at $(round(Int, AGE / 1000)) ka against today [m]", vertical = false, width = cbw, labelsize = 22)
Colorbar(fig[2, 3], colormap = :YlOrRd, limits = (0, smax), highclip = :darkred, label = "site RMS of z_m", vertical = false,
         width = cbw, labelsize = 22)
rowgap!(fig.layout, 1, 5)
elems = [MarkerElement(marker = :utriangle, color = :white, strokecolor = :black, strokewidth = 0.8, markersize = 16),
         MarkerElement(marker = :circle, color = :white, strokecolor = :black, strokewidth = 0.8, markersize = 16),
         LineElement(color = :firebrick, linewidth = 2), LineElement(color = :black, linewidth = 1.5), LineElement(color = :gray30)]
axislegend(ax, elems, ["ice-sheet sites", "far-field sites", "plastic ice margin", "coast, $(round(Int, AGE / 1000)) ka",
                       "coast today"], position = :lb, labelsize = 18, framevisible = true, backgroundcolor = (:white, 0.85))
save(outpng, fig, px_per_unit = 1.5)
println("wrote ", outpng)
