# A margin chronology for the plastic ice backend (margin_file; docs/backend-plastic.md) that takes Greenland from PaleoGrIS
# 1.0 (Leger et al. 2024) inside its time window and PaleoMIST 1.0 (Gowan et al. 2021) everywhere else:
#
#   julia --project=tools/levante tools/levante/margin_merge.jl [<out file>]
#       default $KRYOS_WORK/margins/NH-32KM/NH-32KM_IceHistory-PaleoMIST-PaleoGrIS.nc, and next to it
#       <out>_check.png (Greenland's ice area through time and maps) and <out>_check.md (the numbers)
#
# Use it in an ice instance on NH-32KM in place of the PaleoMIST product (the other margin keys unchanged):
#   margin_file = '<out file>', margin_var = 'f_ice', margin_H_min = 0.5, margin_interp = 'sdist'
#
# Inputs (FesmData, read only): {NH-32KM}_IceHistory-PaleoMIST.nc (f_ice, 33 slices 80-0 ka), GRL-PAL-1KM_IceHistory-PaleoGrIS.nc
# (age_deglac, age_deglac_err: the areas between two isochrones) and GRL-PAL-1KM_TOPO-BedMachine-v5.nc (mask: today's ice). The
# grids share one polar stereographic projection (checked), so the 1 km fields are aggregated exactly by the overlap of the cells.
#
# PaleoGrIS's ages. Each isochrone is a line with a time range (14-13, 13-12, 12-11.5, then every 500 years to 7-6.5 ka); each
# area between two isochrones has age_deglac = the middle of [younger bound of its inner isochrone, older bound of its outer one]
# and age_deglac_err = half that span (13000 +- 1000, 12250 +- 750, 11500 +- 500, ..., 7000 +- 500; checked here). The slice of an
# isochrone is taken at the middle of its range (13.5, 12.5, 11.75, 11.25, ..., 7.25, 6.75 ka): ice where the area deglaciated
# after it (age_deglac < t), under today's ice, and in the cells the reconstruction encloses without an age (inside its innermost
# isochrone, not today's ice: ice until 6.75 ka, then deglaciated with the PaleoMIST slices that follow).
#
# The merge. Greenland is the cells nearer to Greenland's land (the 1 km land component with most of today's ice and the islands
# whose centroid lies within 120 km of it) than to any other land, within 400 km of it. The time axis is PaleoMIST's slices and
# the 13 isochrone times. Inside Greenland, from 13.5 to 6.75 ka, f_ice is PaleoGrIS's ice fraction (its slices; at PaleoMIST's
# 12.5, 10, 7.5 ka its own interpolation, below); outside that window (15 ka and older, 5, 2.5, 0 ka) and outside Greenland
# PaleoMIST's f_ice (its slices; at the isochrone times its own interpolation). kryos's interpolation then blends 15 -> 13.5 ka
# and 6.75 -> 5 ka by itself. In a band of 2 cells either side of Greenland's boundary (Nares Strait, Kane Basin) f_ice is the
# larger of the two within the window, so the Greenland and Innuitian ice do not part at the seam. A source's value at a time
# it has no slice at is the mask (0 or 1) of its signed distance interpolated in time, as kryos_plastic_margin does
# (margin_interp = "sdist"), so the margin of each source moves as it would online.
#
# Limitation: the masks are exact at every slice, but online kryos interpolates the signed distance of the merged slices, so
# between PaleoMIST's slices the margin outside Greenland also moves by about a cell (2.5% of its ice cells at most, at 14.2 ka;
# the check reports it). In the A/B of experiments/coulomb-paleogris.sh this costs NAIS up to 1.1 m SLE at 14 ka and moves its
# 10.4 ka drop by 100 yr. Exact for both sources only with a margin per region in kryos (one chronology each, a region map).
#
# Output: f_ice(xc, yc, time) (time in years ago), greenland(xc, yc) (0 outside, 1 Greenland, 2 the band), source(time)
# (Greenland from 0 PaleoMIST, 1 PaleoGrIS). The check reports the cells of the footprint without a Coulomb region or in a region
# other than Greenland's ($KRYOS_DATA/plastic/{grid}_taub_regions.nc).
using CairoMakie, NCDatasets, Printf, Statistics, Dates

const GRID = "NH-32KM"
const FES = ENV["FESMDATA"]; const KD = ENV["KRYOS_DATA"]
const ISO = [(14000, 13000), (13000, 12000), (12000, 11500), [(t + 500, t) for t in 11000:-500:6500]...]   # isochrone ranges [yr ago]
const T_ISO = [(a + b) / 2 for (a, b) in ISO]                                                              # their slice times
const BAND = 2                                                                                             # cells

out = length(ARGS) >= 1 ? ARGS[1] : joinpath(ENV["KRYOS_WORK"], "margins", GRID, "$(GRID)_IceHistory-PaleoMIST-PaleoGrIS.nc")
mkpath(dirname(out))
f_pm = joinpath(FES, "North", GRID, "$(GRID)_IceHistory-PaleoMIST.nc")
f_pg = joinpath(FES, "GreenlandPaleo", "GRL-PAL-1KM", "GRL-PAL-1KM_IceHistory-PaleoGrIS.nc")
f_bm = joinpath(FES, "GreenlandPaleo", "GRL-PAL-1KM", "GRL-PAL-1KM_TOPO-BedMachine-v5.nc")
f_rg = joinpath(KD, "plastic", "$(GRID)_taub_regions.nc")

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
SD_pm = [sdist(F_pm[:, :, k] .> 0.5, dx, dy) for k in eachindex(t_pm)]

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

# --- the merged slices --------------------------------------------------------------------------------------------------------
t_all = sort(unique([t_pm; T_ISO]), rev = true)
win(t) = T_ISO[end] <= t <= T_ISO[1]
fpm(t) = (k = findfirst(==(t), t_pm); k === nothing ? Float64.(mask_at(t_pm, SD_pm, t)) : F_pm[:, :, k])
fpg(t) = (k = findfirst(==(t), T_ISO); k === nothing ? Float64.(mask_at(T_ISO, SD_pg, t)) : F_pg[:, :, k])
F = zeros(length(x), length(y), length(t_all)); src = zeros(Int8, length(t_all))
for (k, t) in enumerate(t_all)
    a = fpm(t)
    if win(t)
        b = fpg(t); a = ifelse.(G .& .!B, b, ifelse.(B, max.(a, b), a)); src[k] = 1
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
    defVar(o, "time", t_all, ("time",); attrib = ["units" => "years ago", "long_name" => "time (years before 1950 CE)"])
    defVar(o, "f_ice", Float32.(F), ("xc", "yc", "time"); deflatelevel = 1, attrib = [
        "units" => "1", "long_name" => "ice area fraction", "grid_mapping" => "crs",
        "comment" => "PaleoMIST's f_ice outside Greenland and outside 13.5-6.75 ka; inside Greenland in that window PaleoGrIS's ice " *
                     "fraction (1 km, aggregated); at a time a source has no slice: its sdist-interpolated mask (0 or 1)"])
    defVar(o, "greenland", gmap, ("xc", "yc"); attrib = ["long_name" => "Greenland of the merge", "flag_values" => Int8[0, 1, 2],
        "flag_meanings" => "outside greenland seam_band_max_of_both", "grid_mapping" => "crs"])
    defVar(o, "source", src, ("time",); attrib = ["long_name" => "source of Greenland's slice", "flag_values" => Int8[0, 1],
        "flag_meanings" => "PaleoMIST PaleoGrIS"])
    o.attrib["title"] = "Margin chronology on $GRID: PaleoMIST 1.0 with Greenland from PaleoGrIS 1.0 (13.5-6.75 ka)"
    o.attrib["usage"] = "kryos plastic backend: margin_file = '$out', margin_var = 'f_ice', margin_H_min = 0.5, margin_interp = 'sdist'"
    o.attrib["comment"] = "PaleoGrIS slices at the middle of each isochrone's range ($(join(Int.(T_ISO), ", ")) years ago): ice " *
                          "where age_deglac < t, under today's ice (BedMachine v5) and in the cells the reconstruction encloses " *
                          "without an age. Greenland: nearer to Greenland's land than to any other, within 400 km; band of " *
                          "$BAND cells at its boundary with other ice: max of both. Made by tools/levante/margin_merge.jl (kryos)."
    o.attrib["references"] = "Gowan et al. 2021, Nat. Commun. 12, 1199, doi:10.1038/s41467-021-21469-w; Leger et al. 2024, " *
                             "Clim. Past 20, 701-755, doi:10.5194/cp-20-701-2024"
    o.attrib["sources"] = join([f_pm, f_pg, f_bm], ", ")
    o.attrib["history"] = "$(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")): margin_merge.jl"
end
println("wrote $out: $(length(t_all)) slices, Greenland $(count(G)) cells, band $(count(B)) cells")

# --- check: regions, areas, maps ---------------------------------------------------------------------------------------------
reg = NCDataset(ds -> Int.(ds["region"][:, :]), f_rg)
fp = dropdims(any(F .> 0.5, dims = 3), dims = 3)
SD_m = [sdist(F[:, :, k] .> 0.5, dx, dy) for k in eachindex(t_all)]
dA = dx * dy / 1e9                               # 10^3 km2
tt = collect(16000.0:-100.0:4000.0)
area(m) = count(m .& G) * dA
A_pm = [area(mask_at(t_pm, SD_pm, t)) for t in tt]; A_m = [area(mask_at(t_all, SD_m, t)) for t in tt]
A_pg = [win(t) ? area(mask_at(T_ISO, SD_pg, t)) : NaN for t in tt]
io = open(replace(out, r"\.nc$" => "_check.md"), "w")
println(io, "# $(basename(out))\n\nslices: $(join(Int.(t_all), ", ")) years ago (Greenland from PaleoGrIS at ",
        join(Int.(t_all[src .== 1]), ", "), ")\n")
println(io, "Greenland $(count(G)) cells, band $(count(B)) cells; PaleoGrIS: enclosed cells without an age $(count(encl)) km2 (1 km)\n")
println(io, "footprint (ice in any slice): $(count(fp)) cells; without a Coulomb region (would stop the run): ",
        count(fp .& (reg .<= 0)), "; Greenland or band cells of the footprint by region: ",
        join(["$r: $(count(fp .& (G .| B) .& (reg .== r)))" for r in sort(unique(reg[fp .& (G .| B)]))], ", "))
println(io, "new cells (merged footprint, not PaleoMIST's) by region: ",
        join(["$r: $(count(fp .& .!fp_pm .& (reg .== r)))" for r in sort(unique(reg[fp .& .!fp_pm]))], ", "))
# outside Greenland the merged slices give PaleoMIST back exactly at every slice; between them kryos interpolates the signed
# distance of the merged slices, not of PaleoMIST's, which moves the margin within a cell or so: the cells that differ
let O = .!(G .| B), d = [count(O .& (mask_at(t_pm, SD_pm, t) .!= mask_at(t_all, SD_m, t))) for t in tt], a = [count(O .& mask_at(t_pm, SD_pm, t)) for t in tt]
    k = argmax(d)
    println(io, "\noutside Greenland, 16-4 ka every 100 yr: cells where the merged margin differs from PaleoMIST's: max $(d[k]) at ",
            "$(tt[k] / 1000) ka (of $(a[k]) ice cells, $(round(100 * d[k] / a[k], digits = 2))%), mean $(round(mean(d), digits = 1))")
end
println(io, "\n| t ka | PaleoMIST | PaleoGrIS | merged (kryos interpolation) |\n|---|---|---|---|")
for (i, t) in enumerate(tt)
    mod(t, 500) == 0 && println(io, @sprintf("| %.1f | %.0f | %.0f | %.0f |", t / 1000, A_pm[i], A_pg[i], A_m[i]))
end
close(io)
println(read(replace(out, r"\.nc$" => "_check.md"), String))

fig = Figure(size = (1800, 1050))
ax = Axis(fig[1, 1:3], xreversed = true, xlabel = "ka", ylabel = "10³ km²", title = "Greenland's ice area (kryos's interpolation, 100 yr)",
          limits = ((4, 16), nothing))
lines!(ax, tt ./ 1000, A_pm, color = :black, linewidth = 2, label = "PaleoMIST")
lines!(ax, tt ./ 1000, A_pg, color = :red, linewidth = 2, label = "PaleoGrIS")
lines!(ax, tt ./ 1000, A_m, color = :royalblue, linestyle = :dash, linewidth = 2.5, label = "merged")
vlines!(ax, T_ISO ./ 1000, color = (:red, 0.25)); vlines!(ax, t_pm ./ 1000, color = (:black, 0.25), linestyle = :dot)
axislegend(ax, position = :rt)
ax = Axis(fig[1, 4:5], title = "greenland (1) and the band (2), merged footprint", aspect = DataAspect(), limits = ((-900, 1000), (-3500, -500)))
heatmap!(ax, x, y, ifelse.(fp .| (gmap .> 0), Float64.(gmap), NaN), colormap = [:gray80, :royalblue, :orange], colorrange = (0, 2))
contour!(ax, x, y, Float64.(fp), levels = [0.5], color = :black)
for (p, t) in enumerate([14000.0, 12000.0, 10000.0, 8000.0, 7000.0])
    m = mask_at(t_pm, SD_pm, t); g = mask_at(t_all, SD_m, t)
    c = fill(NaN, size(m)); c[m .& g] .= 1; c[m .& .!g] .= 2; c[g .& .!m] .= 3
    local ax = Axis(fig[2, p], aspect = DataAspect(), title = "$(t / 1000) ka", limits = ((-900, 1000), (-3500, -500)),
              xticklabelsvisible = false, yticklabelsvisible = false)
    heatmap!(ax, x, y, c, colormap = cgrad([:gray75, :royalblue, :darkorange], 3, categorical = true), colorrange = (1, 3))
    contour!(ax, x, y, Float64.(G), levels = [0.5], color = :green, linewidth = 0.8)
end
Legend(fig[3, 1:5], [PolyElement(color = :gray75), PolyElement(color = :royalblue), PolyElement(color = :darkorange), LineElement(color = :green)],
       ["both", "PaleoMIST only", "merged only", "Greenland of the merge"], orientation = :horizontal, framevisible = false)
save(replace(out, r"\.nc$" => "_check.png"), fig)
