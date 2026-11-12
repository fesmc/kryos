# Where and when the members of a wave miss the index points of a data group, against the references of the config (is the misfit
# the ice, which the members vary, or the Earth and the data, which they share?):
#   julia --project=tools/levante tools/levante/plastic_ensemble_misfit.jl <ensemble folder> <out folder> [wave=1]
#         [groups=Greenland,Fennoscandia/Barents,Laurentide] [members=36,28] [refs=PaleoMIST,PaleoMIST-3D,PaleoMIST-1D] [top=5]
# References: [[references]] of the config (plastic_ensemble.jl reference_series), those of refs= (default: the first three found);
# the time plot and the site maps need them over all ages (the references that end earlier, GLAC1D's own GIA, are left out there).
# Members: those of the wave that ran; the selection shown: the plausible ones (iter_kk/plausible.txt), or members= (e.g. the
# near misses when none is plausible). The error model: plastic_ensemble.jl's.
# Writes in <out folder>:
#   misfit.md          per group: the implausibility I_g of the references and the members (all points), the bias (mean d) and its spread over
#                      the members, and the sites that carry most of chi2 (d .* R^-1 d summed per site, mean over the selection)
#   misfit_time.png    per group: d = model - obs against age, the median in 1 kyr bins per selected member and per reference
#   misfit_sites.png   per group: the mean d per site, the selection's mean and each reference
#   misfit.md (end), misfit_shared.png   the shared misfit: per site what the good series (the references and the top= members)
#                      get wrong alike, and the selection's member part (the rest)
using CairoMakie, LinearAlgebra, NCDatasets, Printf, Statistics

module EN
include("plastic_ensemble.jl")
end
const ES = EN.ES
const PM = ES.PM

C, out = abspath(ARGS[1]), ARGS[2]
opt = Dict(split(a, "=", limit = 2) for a in ARGS[3:end] if occursin("=", a))
k = length(ARGS) >= 3 && !occursin("=", ARGS[3]) ? parse(Int, ARGS[3]) : 1
groups = split(get(opt, "groups", "Greenland,Fennoscandia/Barents,Laurentide"), ",")
mkpath(out)

cfg = ES.config(C)
ids, _, res, og, obs, A = EN.wave_residuals(C, k, cfg)
og = filter(in(EN.RSL_GROUPS), og)                     # the RSL groups (the GNSS groups have no index points)
em, emsrc = EN.error_model(cfg, C)
plaus = [m for m in sort(collect(keys(res))) if ES.mdir(C, k, m) in strip.(readlines(joinpath(ES.itdir(C, k), "plausible.txt")))]
# the selection shown: the plausible members, or those of members=
sel, sellbl = haskey(opt, "members") ? (parse.(Int, split(opt["members"], ",")), "members " * opt["members"]) : (plaus, "plausible members")
all(in(keys(res)), sel) || error("members not run: ", setdiff(sel, keys(res)))
isempty(sel) && error("no plausible members: select some with members=")
println("members: ", length(res), ", plausible: ", length(plaus), ", shown: ", sel)

# the references at the same index points, over all ages (refs=, default the first three)
fes = cfg["score"]["fesmdata"]
idx = [(g, findall(o -> o.group == g, obs)) for g in og]
refs = [r for r in EN.references(cfg, obs, A) if !any(k -> all(isnan, view(r[2], :, k)), axes(r[2], 2))]
haskey(opt, "refs") && (refs = [r for nm in split(opt["refs"], ",") for r in refs if r[1] == nm])
refs = first(refs, 3); isempty(refs) && error("no reference over all ages")
println("references: ", join(first.(refs), ", "))
lon, lat, zb0 = NCDataset(joinpath(fes, "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc")) do ds     # the coast
    Float64.(ds["lon"][:]), Float64.(ds["lat"][:]), Float64.(nomissing(ds["z_bed"][:, :, findfirst(==(0.0), ds["time"][:])], NaN))
end
rres = [nm => EN.group_residuals(obs, idx, S, A) for (nm, S) in refs]

# per point: chi2 contribution d .* R^-1 d; per site: sums
function contrib(r, p)
    R = ES.obs_cov(r.pts, r.σ, p...)
    r.d .* (Symmetric(R) \ r.d)
end
Ig(r, g) = (v = EN.loglik_group(r, em[g]...); EN.implausibility(v[2], v[3]))
sitekey(o) = string(o.site_name)

io = open(joinpath(out, "misfit.md"), "w")
println(io, "# Misfit by data group: wave $k of $(basename(C)) against the references\n")
println(io, "Error model: $(emsrc). d = model - obs [m]. I_g = (chi2 - n) / sqrt(2n), plausible at <= 3. Members: $(length(res)) ran, ",
        "$(length(plaus)) plausible; shown: $(sellbl). References: ", join(first.(refs), ", "), " (the config's [[references]]).\n")
println(io, "| group | n | ", join(["I " * nm for (nm, _) in refs], " | "), " | I members (min / median / max) | I selection (min / max) | ",
        "bias members [m] (sd over members) | ", join(["bias " * nm * " [m]" for (nm, _) in refs], " | "), " |")
println(io, "|---|---|", repeat("---|", 2length(refs) + 3))
for g in og
    Im = [Ig(res[m].rg[g], g) for m in keys(res)]; Ip = [Ig(res[m].rg[g], g) for m in sel]
    bm = [mean(res[m].rg[g].d) for m in keys(res)]
    println(io, "| ", g, " | ", length(first(values(res)).rg[g].d), " | ", join([@sprintf("%.1f", Ig(rr[g], g)) for (_, rr) in rres], " | "),
            " | ", @sprintf("%.1f / %.1f / %.1f", minimum(Im), median(Im), maximum(Im)), " | ",
            @sprintf("%.1f / %.1f", minimum(Ip), maximum(Ip)), " | ", @sprintf("%+.1f (%.1f)", mean(bm), std(bm)),
            " | ", join([@sprintf("%+.1f", mean(rr[g].d)) for (_, rr) in rres], " | "), " |")
end

# per group: the sites by their share of chi2 (the selection's mean)
sitetab = Dict{String,Any}()
for g in groups
    em_g = em[g]; ms = sel
    S = Dict{String,Any}()
    for m in ms
        r = res[m].rg[g]; c = contrib(r, em_g)
        for (i, o) in enumerate(r.pts)
            s = get!(S, sitekey(o), (lon = o.lon, lat = o.lat, c = 0.0, d = Float64[], t = Float64[]))
            S[sitekey(o)] = (lon = s.lon, lat = s.lat, c = s.c + c[i] / length(ms), d = push!(s.d, r.d[i]), t = push!(s.t, o.t))
        end
    end
    persite(r) = (D = Dict{String,Vector{Float64}}(); for (i, o) in enumerate(r.pts); push!(get!(D, sitekey(o), Float64[]), r.d[i]); end; D)
    dr = [nm => persite(rr[g]) for (nm, rr) in rres]; d3 = last(dr[1])
    sitetab[g] = (S = S, dr = dr)
    tot = sum(s.c for s in values(S))
    n = length(first(values(res)).rg[g].d)
    println(io, "\n## ", g, "\n")
    println(io, @sprintf("chi2 of the selection (mean) %.0f for n = %d; the sites with most of it:\n", tot, n))
    println(io, "| site | lon | lat | points | ages [kyr] | share of chi2 | mean d members [m] | mean d ", first(refs[1]), " [m] |\n|---|---|---|---|---|---|---|---|")
    for (nm, s) in first(sort(collect(S), by = x -> -x[2].c), 12)
        np = length(s.d) ÷ length(ms)
        println(io, "| ", nm, " | ", @sprintf("%.1f | %.1f | %d | %.1f-%.1f | %.0f%% | %+.1f | %+.1f", s.lon, s.lat, np, minimum(s.t) / 1000,
                maximum(s.t) / 1000, 100 * s.c / tot, mean(s.d), haskey(d3, nm) ? mean(d3[nm]) : NaN), " |")
    end
end
close(io)
println("wrote ", joinpath(out, "misfit.md"))

# --- d against age -----------------------------------------------------------------------------------------------------------------
set_theme!(fontsize = 20)
function binned(r; w = 1000.0)
    t = [o.t for o in r.pts]; e = 0:w:(ceil(maximum(t) / w) * w)
    xs = Float64[]; ys = Float64[]
    for a in e[1:end-1]
        j = findall(x -> a <= x < a + w, t); length(j) >= 3 || continue
        push!(xs, (a + w / 2) / 1000); push!(ys, median(r.d[j]))
    end
    xs, ys
end
RSTYLE = [(:black, :dash), (RGBf(0.85, 0.37, 0.0), :solid), (RGBf(0.13, 0.35, 0.65), :dot)]
RC = Dict(nm => RSTYLE[i] for (i, (nm, _)) in enumerate(rres))
fig = Figure(size = (1900, 420 * length(groups)))
for (row, g) in enumerate(groups)
    tmax = maximum(o.t for o in res[first(sel)].rg[g].pts) / 1000
    ax = Axis(fig[row, 1], title = g, xlabel = "age [kyr]", ylabel = "model - obs [m]", xreversed = true,
              limits = ((0, min(tmax, 25.0) + 0.5), nothing))
    σmed = median(res[first(sel)].rg[g].σ)
    band!(ax, [0, 30], [-σmed, -σmed], [σmed, σmed], color = (:gray, 0.15))
    hlines!(ax, 0, color = :black, linewidth = 0.8)
    mr = [res[m].rg[g] for m in sel]
    scatter!(ax, [o.t / 1000 for o in mr[1].pts], mean(r.d for r in mr if length(r.d) == length(mr[1].d)), color = (:gray40, 0.35), markersize = 5)
    for r in mr; lines!(ax, binned(r)..., color = (:gray30, 0.6), linewidth = 1.2); end
    for (nm, rr) in rres; c, ls = RC[nm]; lines!(ax, binned(rr[g])..., color = c, linestyle = ls, linewidth = 3); end
    row == 1 && Legend(fig[1, 2], [LineElement(color = (:gray30, 0.6)), [LineElement(color = RC[nm][1], linestyle = RC[nm][2], linewidth = 3) for (nm, _) in rres]...,
                                MarkerElement(marker = :circle, color = (:gray40, 0.5), markersize = 8), PolyElement(color = (:gray, 0.15))],
                           ["$sellbl (median per kyr)", first.(rres)..., "index points (members' mean)", "± median 1σ of the points"],
                           labelsize = 16, framevisible = false, valign = :top)
end
save(joinpath(out, "misfit_time.png"), fig, px_per_unit = 1.3)
println("wrote ", joinpath(out, "misfit_time.png"))

# --- mean d per site -------------------------------------------------------------------------------------------------------------
BOX = Dict("Greenland" => ((-75.0, -10.0), (59.0, 84.0)), "Fennoscandia/Barents" => ((-12.0, 65.0), (52.0, 82.0)),
           "Laurentide" => ((-140.0, -50.0), (40.0, 82.0)), "Antarctica" => ((-180.0, 180.0), (-90.0, -60.0)))
glon = [x > 180 ? x - 360 : x for x in lon]; p = sortperm(glon); q = sortperm(lat); glon = glon[p]; glat = lat[q]; zb = zb0[p, q]
DMAP = cgrad(:RdBu, rev = true); DLIM = 40.0
# one map panel of group g: values v [m] at the sites (lon, lat), marker sizes ms, today's coast (z_bed = 0)
function sitemap!(fig, row, col, g, title, lon, lat, v, ms)
    (x0, x1), (y0, y1) = get(BOX, g, ((-180.0, 180.0), (-90.0, 90.0))); c0 = cosd((y0 + y1) / 2)
    ax = Axis(fig[row, col], title = title, aspect = DataAspect(), limits = ((x0 * c0, x1 * c0), (y0, y1)),
              xticksvisible = false, xticklabelsvisible = false)
    for l in Makie.Contours.lines(Makie.Contours.contour(glon, glat, zb, 0.0))
        xs, ys = Makie.Contours.coordinates(l)
        lines!(ax, [i > 1 && abs(xs[i] - xs[i-1]) > 5 ? Point2f(NaN, NaN) : Point2f(xs[i] * c0, ys[i]) for i in eachindex(xs)],
               color = :gray40, linewidth = 0.8)
    end
    scatter!(ax, lon .* c0, lat, color = v, colormap = DMAP, colorrange = (-DLIM, DLIM), markersize = ms, strokecolor = :black, strokewidth = 0.6)
end
dcolorbar(fig, rows, col) = Colorbar(fig[rows, col], colormap = DMAP, limits = (-DLIM, DLIM), lowclip = DMAP[0.0], highclip = DMAP[1.0],
                                     label = "mean model - obs per site [m]", height = Relative(0.6))
fig = Figure(size = (900 * length(groups), 730 * (1 + length(rres))))
for (col, g) in enumerate(groups)
    st = sitetab[g]
    for (row, (lbl, D)) in enumerate([("$sellbl (mean)", Dict(nm => s.d for (nm, s) in st.S)); st.dr])
        ks = [nm for nm in keys(D) if haskey(st.S, nm)]
        sitemap!(fig, row, col, g, "$g: $lbl", [st.S[nm].lon for nm in ks], [st.S[nm].lat for nm in ks], [mean(D[nm]) for nm in ks],
                 [8 + 3sqrt(length(st.S[nm].d) / length(sel)) for nm in ks])
    end
end
dcolorbar(fig, 1:(1 + length(rres)), length(groups) + 1)
save(joinpath(out, "misfit_sites.png"), fig, px_per_unit = 1.2)
println("wrote ", joinpath(out, "misfit_sites.png"))

# --- the shared misfit -------------------------------------------------------------------------------------------------------------
# What every good series gets wrong alike (structural: the data, the Earth, what no member and no reconstruction changes) against what
# the members change. On the judging window (compare_max_age_ka): the good series are the references with records there (all of
# them, the published ones too) and the `top` members (default 5) of the highest log L over the groups. Per site the mean d of each
# good series; the shared misfit b = their median; the member part of the selection = its mean d - b. Per group the shared fraction
# of the selection's misfit: sum over sites of n b^2 / sum of n dbar^2 (n the points of the site).
W = EN.score_window(obs, A, og, 1000.0 * get(cfg["ensemble"], "compare_max_age_ka", 25.0))
wrefs = [r for r in EN.references(cfg, obs, A) if !any(k -> all(isnan, view(r[2], W.kp, k)), W.ka)]
ntop = parse(Int, get(opt, "top", "5"))
wsc = Dict(m => EN.score_series(res[m].S, W, em, og) for m in keys(res))
top = first(sort(collect(keys(wsc)), by = m -> -sum(wsc[m][g][1] for g in og)), ntop)
wres(S) = EN.group_residuals(W.obs, W.idx, S[W.kp, W.ka], W.A)
good = [[nm => wres(S) for (nm, S) in wrefs]; ["m" * lpad(m, 3, "0") => wres(res[m].S) for m in top]]
selr = [wres(res[m].S) for m in sel]
sitemean(r) = (D = Dict{String,Tuple{Float64,Float64,Vector{Float64},Vector{Float64}}}();
               for (i, o) in enumerate(r.pts); e = get!(D, sitekey(o), (o.lon, o.lat, Float64[], Float64[])); push!(e[3], r.d[i]); push!(e[4], o.t); end; D)
io = open(joinpath(out, "misfit.md"), "a")
println(io, "\n# The shared misfit (index points of age <= ", W.tmax / 1000, " ka)\n\nGood series: ", join(first.(good), ", "),
        ". Shared b = the median of their mean d per site; the selection's member part = its mean d - b. Shared fraction: ",
        "sum n b^2 / sum n dbar^2 over the sites (n points, dbar the selection's mean d).\n")
shared = Dict{String,Any}()
for g in groups
    G = [sitemean(r[g]) for (_, r) in good]; Ssel = [sitemean(r[g]) for r in selr]
    sites = [s for s in keys(G[1]) if all(haskey(D, s) for D in G) && all(haskey(D, s) for D in Ssel)]
    rows = map(sites) do s
        v = [mean(D[s][3]) for D in G]; b = median(v); ds = mean(mean(D[s][3]) for D in Ssel)
        (site = s, lon = G[1][s][1], lat = G[1][s][2], n = length(G[1][s][3]), t = G[1][s][4], b = b, lo = minimum(v), hi = maximum(v),
         same = all(>(0), v) || all(<(0), v), d = ds, dm = ds - b)
    end
    shared[g] = rows
    fs = sum(r.n * r.b^2 for r in rows) / sum(r.n * r.d^2 for r in rows)
    println(io, "## ", g, @sprintf("\n\nShared fraction of the selection's misfit: %.0f %%. Sites by n b^2:\n", 100fs))
    println(io, "| site | lon | lat | points | ages [kyr] | shared b [m] | range of the good series [m] | one sign | selection d [m] | member part [m] |")
    println(io, "|---|---|---|---|---|---|---|---|---|---|")
    for r in first(sort(rows, by = r -> -r.n * r.b^2), 15)
        println(io, "| ", r.site, @sprintf(" | %.1f | %.1f | %d | %.1f-%.1f | %+.1f | %+.1f to %+.1f | %s | %+.1f | %+.1f |", r.lon, r.lat, r.n,
                minimum(r.t) / 1000, maximum(r.t) / 1000, r.b, r.lo, r.hi, r.same ? "yes" : "no", r.d, r.dm))
    end
    println(io)
end
close(io)
fig = Figure(size = (900 * length(groups), 1460))
for (col, g) in enumerate(groups)
    rows = shared[g]; ms = [8 + 3sqrt(r.n) for r in rows]
    sitemap!(fig, 1, col, g, "$g: shared by the good series (median)", [r.lon for r in rows], [r.lat for r in rows], [r.b for r in rows], ms)
    sitemap!(fig, 2, col, g, "$g: $sellbl minus shared", [r.lon for r in rows], [r.lat for r in rows], [r.dm for r in rows], ms)
end
dcolorbar(fig, 1:2, length(groups) + 1)
save(joinpath(out, "misfit_shared.png"), fig, px_per_unit = 1.2)
println("wrote ", joinpath(out, "misfit_shared.png"))
