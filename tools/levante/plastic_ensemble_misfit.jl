# Where and when the members of a wave miss the index points of a data group, against the runs with PaleoMIST's ice (is the misfit
# the ice, which the members vary, or the Earth and the data, which they share?):
#   julia --project=tools/levante tools/levante/plastic_ensemble_misfit.jl <ensemble folder> <out folder> [wave=1]
#         [groups=Greenland,Fennoscandia/Barents,Laurentide]
# References (PM_REF, default /work/ba1442/kryos_data/paleomist): PaleoMIST (Gowan et al. 2021: their ice and SELEN, 1-D Earth),
# B1 (VILMA2, PaleoMIST's ice, the 1-D Earth of Gowan, lmax 170), B3 (VILMA2, PaleoMIST's ice, the 3-D Earth of Bagge, lmax 170).
# Members: those of the wave that ran, the plausible ones (iter_kk/plausible.txt) marked. The error model: plastic_ensemble.jl's.
# Writes in <out folder>:
#   misfit.md          per group: the implausibility I_g of the references and the members, the bias (mean d) and its spread over
#                      the members, and the sites that carry most of chi2 (d .* R^-1 d summed per site, mean over the plausible members)
#   misfit_time.png    per group: d = model - obs against age, the median in 1 kyr bins per plausible member and per reference
#   misfit_sites.png   per group: the mean d per site, the plausible members' mean, B3 and B1
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
pmref = get(ENV, "PM_REF", "/work/ba1442/kryos_data/paleomist")
mkpath(out)

cfg = ES.config(C)
ids, _, res, og = EN.wave_residuals(C, k, cfg)
em, emsrc = EN.error_model(cfg, C)
plaus = [m for m in sort(collect(keys(res))) if ES.mdir(C, k, m) in strip.(readlines(joinpath(ES.itdir(C, k), "plausible.txt")))]
println("members: ", length(res), ", plausible: ", length(plaus))

# the references at the same index points
fes = cfg["score"]["fesmdata"]
obs = [o for o in PM.read_indicators(fes) if o.type == 0 && o.group in og]
idx = [(g, findall(o -> o.group == g, obs)) for g in og]
pmf = joinpath(fes, "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc")
lon, lat, ages, sl, zb0 = NCDataset(pmf) do ds
    a = Float64.(ds["time"][:]); k0 = findfirst(==(0.0), a)
    Float64.(ds["lon"][:]), Float64.(ds["lat"][:]), a, Float64.(nomissing(ds["dz_sl"][:, :, :], NaN)),
    Float64.(nomissing(ds["z_bed"][:, :, k0], NaN))
end
A = sort(ages)
plon, plat = [o.lon for o in obs], [o.lat for o in obs]
refs = ["PaleoMIST" => PM.series_pm(lon, lat, ages, sl, plon, plat, A)]
for (nm, dir) in ["B1 (1-D)" => "B1_l170", "B3 (3-D)" => "B3_l170"]
    r = PM.read_run(joinpath(pmref, dir)); r === nothing && (println("no ", dir); continue)
    push!(refs, nm => PM.series_run(r, plon, plat, A))
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
println(io, "# Misfit by data group: wave $k of $(basename(C)) against PaleoMIST's ice\n")
println(io, "Error model: $(emsrc). d = model - obs [m]. I_g = (chi2 - n) / sqrt(2n), plausible at <= 3. Members: $(length(res)) ran, ",
        "$(length(plaus)) plausible. References: PaleoMIST (Gowan et al. 2021, SELEN, 1-D Earth), B1 and B3 (VILMA2 with PaleoMIST's ice, ",
        "1-D and 3-D Earth); the members have the 3-D Earth and their own ice.\n")
println(io, "| group | n | I PaleoMIST | I B1 | I B3 | I members (min / median / max) | I plausible (min / max) | bias members [m] (sd over members) | bias B3 [m] |")
println(io, "|---|---|---|---|---|---|---|---|---|")
for g in og
    Im = [Ig(res[m].rg[g], g) for m in keys(res)]; Ip = [Ig(res[m].rg[g], g) for m in plaus]
    bm = [mean(res[m].rg[g].d) for m in keys(res)]
    rI = Dict(nm => Ig(rr[g], g) for (nm, rr) in rres)
    b3 = haskey(Dict(rres), "B3 (3-D)") ? mean(Dict(rres)["B3 (3-D)"][g].d) : NaN
    println(io, "| ", g, " | ", length(first(values(res)).rg[g].d), " | ", join([@sprintf("%.1f", get(rI, nm, NaN)) for nm in ("PaleoMIST", "B1 (1-D)", "B3 (3-D)")], " | "),
            " | ", @sprintf("%.1f / %.1f / %.1f", minimum(Im), median(Im), maximum(Im)), " | ",
            isempty(Ip) ? "-" : @sprintf("%.1f / %.1f", minimum(Ip), maximum(Ip)), " | ", @sprintf("%+.1f (%.1f)", mean(bm), std(bm)),
            " | ", @sprintf("%+.1f", b3), " |")
end

# per group: the sites by their share of chi2 (plausible members' mean)
sitetab = Dict{String,Any}()
for g in groups
    em_g = em[g]; ms = isempty(plaus) ? collect(keys(res)) : plaus
    S = Dict{String,Any}()
    for m in ms
        r = res[m].rg[g]; c = contrib(r, em_g)
        for (i, o) in enumerate(r.pts)
            s = get!(S, sitekey(o), (lon = o.lon, lat = o.lat, c = 0.0, d = Float64[], t = Float64[]))
            S[sitekey(o)] = (lon = s.lon, lat = s.lat, c = s.c + c[i] / length(ms), d = push!(s.d, r.d[i]), t = push!(s.t, o.t))
        end
    end
    persite(r) = (D = Dict{String,Vector{Float64}}(); for (i, o) in enumerate(r.pts); push!(get!(D, sitekey(o), Float64[]), r.d[i]); end; D)
    d3 = persite(Dict(rres)["B3 (3-D)"][g])
    sitetab[g] = (S = S, d3 = d3, d1 = persite(Dict(rres)["B1 (1-D)"][g]))
    tot = sum(s.c for s in values(S))
    n = length(first(values(res)).rg[g].d)
    println(io, "\n## ", g, "\n")
    println(io, @sprintf("chi2 of the plausible members (mean) %.0f for n = %d; the sites with most of it:\n", tot, n))
    println(io, "| site | lon | lat | points | ages [kyr] | share of chi2 | mean d members [m] | mean d B3 [m] |\n|---|---|---|---|---|---|---|---|")
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
RC = Dict("PaleoMIST" => (:black, :dash), "B1 (1-D)" => (RGBf(0.13, 0.35, 0.65), :dot), "B3 (3-D)" => (RGBf(0.85, 0.37, 0.0), :solid))
fig = Figure(size = (1900, 420 * length(groups)))
for (row, g) in enumerate(groups)
    tmax = maximum(o.t for o in res[first(plaus)].rg[g].pts) / 1000
    ax = Axis(fig[row, 1], title = g, xlabel = "age [kyr]", ylabel = "model - obs [m]", xreversed = true,
              limits = ((0, min(tmax, 25.0) + 0.5), nothing))
    σmed = median(res[first(plaus)].rg[g].σ)
    band!(ax, [0, 30], [-σmed, -σmed], [σmed, σmed], color = (:gray, 0.15))
    hlines!(ax, 0, color = :black, linewidth = 0.8)
    mr = [res[m].rg[g] for m in plaus]
    scatter!(ax, [o.t / 1000 for o in mr[1].pts], mean(r.d for r in mr if length(r.d) == length(mr[1].d)), color = (:gray40, 0.35), markersize = 5)
    for r in mr; lines!(ax, binned(r)..., color = (:gray30, 0.6), linewidth = 1.2); end
    for (nm, rr) in rres; c, ls = RC[nm]; lines!(ax, binned(rr[g])..., color = c, linestyle = ls, linewidth = 3); end
    row == 1 && Legend(fig[1, 2], [LineElement(color = (:gray30, 0.6)), [LineElement(color = RC[nm][1], linestyle = RC[nm][2], linewidth = 3) for (nm, _) in rres]...,
                                MarkerElement(marker = :circle, color = (:gray40, 0.5), markersize = 8), PolyElement(color = (:gray, 0.15))],
                           ["plausible members (median per kyr)", first.(rres)..., "index points (members' mean)", "± median 1σ of the points"],
                           labelsize = 16, framevisible = false, valign = :top)
end
save(joinpath(out, "misfit_time.png"), fig, px_per_unit = 1.3)
println("wrote ", joinpath(out, "misfit_time.png"))

# --- mean d per site -------------------------------------------------------------------------------------------------------------
BOX = Dict("Greenland" => ((-75.0, -10.0), (59.0, 84.0)), "Fennoscandia/Barents" => ((-12.0, 65.0), (52.0, 82.0)),
           "Laurentide" => ((-140.0, -50.0), (40.0, 82.0)), "Antarctica" => ((-180.0, 180.0), (-90.0, -60.0)))
glon = [x > 180 ? x - 360 : x for x in lon]; p = sortperm(glon); q = sortperm(lat); glon = glon[p]; glat = lat[q]; zb = zb0[p, q]
DMAP = cgrad(:RdBu, rev = true); DLIM = 40.0
fig = Figure(size = (900 * length(groups), 2200))
for (col, g) in enumerate(groups)
    (x0, x1), (y0, y1) = get(BOX, g, ((-180.0, 180.0), (-90.0, 90.0))); c0 = cosd((y0 + y1) / 2)
    st = sitetab[g]
    for (row, (lbl, D)) in enumerate([("plausible members (mean)", Dict(nm => s.d for (nm, s) in st.S)), ("B3: PaleoMIST ice, 3-D Earth", st.d3),
                                         ("B1: PaleoMIST ice, 1-D Earth", st.d1)])
        ax = Axis(fig[row, col], title = "$g: $lbl", aspect = DataAspect(), limits = ((x0 * c0, x1 * c0), (y0, y1)),
                  xticksvisible = false, xticklabelsvisible = false)
        for l in Makie.Contours.lines(Makie.Contours.contour(glon, glat, zb, 0.0))      # the coast of today (z_bed = 0)
            xs, ys = Makie.Contours.coordinates(l)
            lines!(ax, [i > 1 && abs(xs[i] - xs[i-1]) > 5 ? Point2f(NaN, NaN) : Point2f(xs[i] * c0, ys[i]) for i in eachindex(xs)],
                   color = :gray40, linewidth = 0.8)
        end
        ks = [nm for nm in keys(D) if haskey(st.S, nm)]
        scatter!(ax, [st.S[nm].lon * c0 for nm in ks], [st.S[nm].lat for nm in ks], color = [mean(D[nm]) for nm in ks], colormap = DMAP,
                 colorrange = (-DLIM, DLIM), markersize = [8 + 3sqrt(length(st.S[nm].d) / max(length(plaus), 1)) for nm in ks],
                 strokecolor = :black, strokewidth = 0.6)
    end
end
Colorbar(fig[1:3, length(groups) + 1], colormap = DMAP, limits = (-DLIM, DLIM), lowclip = DMAP[0.0], highclip = DMAP[1.0],
         label = "mean model - obs per site [m]", height = Relative(0.6))
save(joinpath(out, "misfit_sites.png"), fig, px_per_unit = 1.2)
println("wrote ", joinpath(out, "misfit_sites.png"))
