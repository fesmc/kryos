# Diagnostic figures of a tau_b calibration (plastic_esmda.jl): where the posterior moves the ice, what the sea-level data ask for,
# and how strongly each group's parameters act on each set of data.
#   julia --project=<julia-plot> tools/levante/plastic_esmda_diag.jl <calibration folder> <out folder> <prior run> <posterior run>
# The runs: folders with ice_nor.nc and ice_ant.nc (2-D, the 2.5 kyr records) and geo_glob_vilma.nc (read in place by
# pm_compare.jl's read_run); for v2 the prior is runs/v2-baseline (= m000 of iteration 1, which has no 2-D output) and the
# posterior m000 of the last iteration. Writes to <out>:
#   volumes.png    ice volume [m SLE] by scale group (NH-32KM cells by the group of their tau_b domain; ANT-16KM), PaleoMIST, prior and
#                  posterior; ESL (VILMA2 bsl_ive) against PaleoMIST's and Lambeck et al. (2014); posterior - prior by group
#   sensitivity.png  iteration 1 (the prior ensemble): the mean misfit d of the index points of each data set against the factor of the
#                  parameters (Antarctica 80-20 ka, Laurentide 20 ka), with the posterior mean marked
#   residuals.png  d = model - observed RSL of the index points against age, per data set, prior and posterior (2.5 kyr bin means)
#   ant_sites.png  the Antarctic sites: index points (2 sigma) and limiting data, model RSL of the prior, the posterior and PaleoMIST
#   ant_maps.png   ANT-16KM ice thickness at 20 ka: PaleoMIST, prior, posterior, posterior - prior
using CairoMakie, NCDatasets, Printf, Statistics, DelimitedFiles, TOML

module PM
include("pm_compare.jl")
end

const RHO_I, RHO_SW, A_OCN = 910.0, 1028.0, 3.618e14
const GROUP_COLORS = Dict("Laurentide" => :steelblue, "Greenland" => :seagreen, "Eurasia" => :darkorange, "Antarctica" => :firebrick,
                          "other" => :gray50)

C, out, prior_dir, post_dir = ARGS
mkpath(out)
cfg = TOML.parsefile(joinpath(C, "config.toml"))
fes = cfg["score"]["fesmdata"]

# --- ice volume by scale group --------------------------------------------------------------------------------------------------------
gname = Dict(Int(r[1]) => String(r[6]) for r in eachrow(readdlm(joinpath(C, "groups.csv"), ',', Any; skipstart = 1)))
function group_mask(dom, grid)               # group name of each cell of the grid by its tau_b domain ("other": no domain)
    NCDataset(joinpath(fes, dom, grid, "$(grid)_IceHistory-PaleoMIST_taub.nc")) do ds
        D = ds["taub_domain"][:, :]
        [ismissing(d) ? "other" : get(gname, Int(d), "other") for d in D]
    end
end
cell_sle(x, y) = (x[2] - x[1]) * (y[2] - y[1]) * 1e6 * RHO_I / RHO_SW / A_OCN
function volumes(f, G, groups; tvar = "time", hvar = "H_ice", years_ago = false)   # m SLE per group at the 2.5 kyr ages (ka)
    NCDataset(f) do ds
        x = ds["xc"][:]; y = ds["yc"][:]; t = Float64.(ds[tvar][:]); a = years_ago ? t : -t
        w = cell_sle(x, y)
        V = Dict(g => zeros(length(t)) for g in groups)
        for k in eachindex(t)
            H = nomissing(ds[hvar][:, :, k], 0.0)
            for g in groups; V[g][k] = sum(H[G .== g]) * w; end
        end
        p = sortperm(a); (a[p] ./ 1000, Dict(g => v[p] for (g, v) in V))
    end
end
GN = group_mask("North", "NH-32KM"); GA = group_mask("Antarctica", "ANT-16KM")
nh_groups = ["Laurentide", "Greenland", "Eurasia", "other"]
vol = Dict{String,Any}()
vol["PaleoMIST"] = merge(volumes(joinpath(fes, "North", "NH-32KM", "NH-32KM_IceHistory-PaleoMIST.nc"), GN, nh_groups; years_ago = true)[2],
                         volumes(joinpath(fes, "Antarctica", "ANT-16KM", "ANT-16KM_IceHistory-PaleoMIST.nc"), fill("Antarctica", size(GA)),
                                 ["Antarctica"]; years_ago = true)[2])
ages_pm = volumes(joinpath(fes, "Antarctica", "ANT-16KM", "ANT-16KM_IceHistory-PaleoMIST.nc"), fill("Antarctica", size(GA)), ["Antarctica"];
                  years_ago = true)[1]
ages = Float64[]
for (key, D) in (("prior", prior_dir), ("posterior", post_dir))
    global ages
    a, vn = volumes(joinpath(D, "ice_nor.nc"), GN, nh_groups)
    _, va = volumes(joinpath(D, "ice_ant.nc"), fill("Antarctica", size(GA)), ["Antarctica"])
    ages = a; vol[key] = merge(vn, va)
end
runs = Dict("prior" => PM.read_run(prior_dir), "posterior" => PM.read_run(post_dir))
esl(r) = (ka = -r.t ./ 1000; p = sortperm(ka); (ka[p], (r.bsl_ive .- r.bsl_ive[end])[p]))
slrec(name, v) = NCDataset(joinpath(fes, "Records", "SeaLevel", "SeaLevel-$name.nc")) do ds
    t = Float64.(ds["time"][:]); M = Float64.(nomissing(ds[v][:, :], NaN)); s = size(M, 1) == length(t) ? M[:, 1] : M[1, :]
    a = abs.(t) ./ 1000; p = sortperm(a); (a[p], s[p])
end

fig = Figure(size = (1500, 900))
groups_all = ["Laurentide", "Greenland", "Eurasia", "Antarctica"]
for (i, g) in enumerate(groups_all)
    ax = Axis(fig[(i - 1) ÷ 2 + 1, (i - 1) % 2 + 1], title = "$g: ice volume", xlabel = "age [ka]", ylabel = "m SLE", xreversed = true)
    lines!(ax, ages_pm, vol["PaleoMIST"][g], color = :black, linestyle = :dash, label = "PaleoMIST")
    lines!(ax, ages, vol["prior"][g], color = :gray50, linewidth = 2, label = "prior mean (Stage 1)")
    lines!(ax, ages, vol["posterior"][g], color = GROUP_COLORS[g], linewidth = 2.5, label = "posterior mean")
    i == 1 && axislegend(ax, position = :lt, framevisible = false)
end
ax = Axis(fig[1, 3], title = "ESL (VILMA2, relative to 0 ka)", xlabel = "age [ka]", ylabel = "m", xreversed = true)
try; a, s = slrec("Lambeck2014", "sl"); lines!(ax, a, s, color = :purple, linestyle = :dot, label = "Lambeck et al. 2014")
catch e; @warn "Lambeck2014 not plotted" e; end
a, s = slrec("PaleoMIST", "sl"); lines!(ax, a, s, color = :black, linestyle = :dash, label = "PaleoMIST")
for (key, c) in (("prior", :gray50), ("posterior", :firebrick))
    a, s = esl(runs[key]); lines!(ax, a, s, color = c, linewidth = 2, label = "$key mean")
end
axislegend(ax, position = :rb, framevisible = false)
ax = Axis(fig[2, 3], title = "posterior - prior, ice volume", xlabel = "age [ka]", ylabel = "m SLE", xreversed = true)
for g in groups_all
    lines!(ax, ages, vol["posterior"][g] .- vol["prior"][g], color = GROUP_COLORS[g], linewidth = 2, label = g)
end
tot = sum(vol["posterior"][g] .- vol["prior"][g] for g in groups_all)
lines!(ax, ages, tot, color = :black, linewidth = 2, label = "sum")
hlines!(ax, [0.0], color = :gray70)
axislegend(ax, position = :lt, framevisible = false)
save(joinpath(out, "volumes.png"), fig)
k20 = findfirst(==(20.0), ages); k20p = findfirst(==(20.0), ages_pm)
println("20 ka [m SLE]: ", join(("$g PaleoMIST $(round(vol["PaleoMIST"][g][k20p]; digits=1)) prior $(round(vol["prior"][g][k20]; digits=1)) " *
                                 "posterior $(round(vol["posterior"][g][k20]; digits=1))" for g in groups_all), "; "))

# --- sensitivity in the prior ensemble ------------------------------------------------------------------------------------------------
th, hd = readdlm(joinpath(C, "iter_01", "theta.csv"), ','; header = true)
col(n) = findfirst(==(n), vec(hd))
sc, sh = readdlm(joinpath(C, "iter_01", "scores.csv"), ','; header = true)
setcol = findfirst(==("set"), vec(sh)); dcol = findfirst(==("mean_d"), vec(sh)); rcol = findfirst(==("rms_z"), vec(sh))
kl = maximum(parse(Int, d[6:7]) for d in readdir(C) if startswith(d, "iter_"))
thpost, _ = readdlm(joinpath(C, @sprintf("iter_%02d", kl), "theta.csv"), ','; header = true)
sets = ["far field", "Antarctica", "Laurentide", "Fennoscandia/Barents"]
pars = ["Antarctica_80-60-40-25-20", "Laurentide_20"]
fig = Figure(size = (1500, 650))
for (j, p) in enumerate(pars), (i, s) in enumerate(sets)
    rows = [r for r in eachrow(sc) if r[setcol] == s]
    mem = Int.([r[1] for r in rows]); d = Float64.([r[dcol] for r in rows])
    f = [exp(th[findfirst(==(m), Int.(th[:, 1])), col(p)]) for m in mem]
    ax = Axis(fig[j, i], title = "$s index points", xlabel = "factor $p (iteration 1 members)", ylabel = "mean d [m]",
              xscale = log10)
    scatter!(ax, f, d, color = (:gray30, 0.7))
    k = isfinite.(d); b = cov(log.(f[k]), d[k]) / var(log.(f[k]))
    vlines!(ax, [exp(thpost[1, col(p)])], color = :firebrick, label = "posterior mean")
    hlines!(ax, [0.0], color = :gray70)
    text!(ax, 0.02, 0.98, text = @sprintf("slope %.1f m per e-fold\ncorr %.2f", b, cor(log.(f[k]), d[k])), space = :relative,
          align = (:left, :top), fontsize = 12)
end
save(joinpath(out, "sensitivity.png"), fig)

# --- residuals by age -----------------------------------------------------------------------------------------------------------------
obs = [o for o in PM.read_indicators(fes) if o.group in cfg["esmda"]["obs_groups"]]
A = sort(NCDataset(ds -> Float64.(ds["time"][:]), joinpath(fes, "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc")))
S = Dict(k => PM.series_run(r, [o.lon for o in obs], [o.lat for o in obs], A) for (k, r) in runs)
pm = NCDataset(joinpath(fes, "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc")) do ds
    (lon = Float64.(ds["lon"][:]), lat = Float64.(ds["lat"][:]), t = Float64.(ds["time"][:]), sl = Float64.(nomissing(ds["dz_sl"][:, :, :], NaN)))
end
S["PaleoMIST"] = PM.series_pm(pm.lon, pm.lat, pm.t, pm.sl, [o.lon for o in obs], [o.lat for o in obs], A)
dval(k, i) = PM.tval(A, S[k][i, :], obs[i].t)[1] - obs[i].rsl
fig = Figure(size = (1500, 800))
for (j, g) in enumerate(cfg["esmda"]["obs_groups"])
    ip = [i for (i, o) in enumerate(obs) if o.group == g && o.type == 0]
    ax = Axis(fig[(j - 1) ÷ 3 + 1, (j - 1) % 3 + 1], title = "$g ($(length(ip)) index points)", xlabel = "age [ka]",
              ylabel = "d = model - observed [m]", xreversed = true)
    hlines!(ax, [0.0], color = :gray70)
    for (k, c) in (("prior", :gray50), ("posterior", :firebrick), ("PaleoMIST", :black))
        t = [obs[i].t / 1000 for i in ip]; d = [dval(k, i) for i in ip]
        k == "posterior" && scatter!(ax, t, d, color = (c, 0.25), markersize = 4)
        bins = 0:2.5:maximum(t; init = 0) + 2.5
        bc = Float64[]; bm = Float64[]
        for b in bins[1:end-1]
            m = (b .<= t .< b + 2.5) .& isfinite.(d); count(m) >= 3 || continue
            push!(bc, b + 1.25); push!(bm, mean(d[m]))
        end
        scatterlines!(ax, bc, bm, color = c, linestyle = k == "PaleoMIST" ? :dash : :solid, linewidth = 2, markersize = 6,
                      label = "$k: mean $(@sprintf("%.1f", mean(filter(isfinite, d)))) m")
    end
    axislegend(ax, position = :lb, framevisible = false, labelsize = 11)
end
save(joinpath(out, "residuals.png"), fig)

# --- Antarctic sites --------------------------------------------------------------------------------------------------------------------
ant = [i for (i, o) in enumerate(obs) if o.group == "Antarctica"]
sites = unique(obs[i].site_name for i in ant if obs[i].type == 0)
sites = sort(sites, by = s -> -count(i -> obs[i].site_name == s && obs[i].type == 0, ant))
n = length(sites); nc = 3
fig = Figure(size = (1500, 420 * cld(n, nc)))
for (j, s) in enumerate(sites)
    idx = [i for i in ant if obs[i].site_name == s]
    lo = mean(obs[i].lon for i in idx); la = mean(obs[i].lat for i in idx)
    ip = [i for i in idx if obs[i].type == 0]
    ax = Axis(fig[(j - 1) ÷ nc + 1, (j - 1) % nc + 1], title = @sprintf("%s (%.1f, %.1f): %d IP, %d limiting", s, la, lo, length(ip),
              length(idx) - length(ip)), xlabel = "age [ka]", ylabel = "RSL [m]", xreversed = true)
    xl = min(12.5, maximum(obs[i].tmax for i in ip) / 1000 + 1.5)
    Ap = A[A .<= 15000.0]
    for (k, c) in (("PaleoMIST", :black), ("prior", :gray50), ("posterior", :firebrick))
        r = k == "PaleoMIST" ? vec(PM.series_pm(pm.lon, pm.lat, pm.t, pm.sl, [lo], [la], Ap)) : vec(PM.series_run(runs[k], [lo], [la], Ap))
        lines!(ax, Ap ./ 1000, r, color = c, linestyle = k == "PaleoMIST" ? :dash : :solid, linewidth = 2, label = k)
    end
    for i in idx
        o = obs[i]; t = o.t / 1000
        if o.type == 0
            rangebars!(ax, [t], [o.rmin], [o.rmax], color = :royalblue); rangebars!(ax, [o.rsl], [o.tmin / 1000], [o.tmax / 1000],
                       direction = :x, color = :royalblue)
        else
            scatter!(ax, [t], [o.type == -1 ? o.rmin : o.rmax], marker = o.type == -1 ? :utriangle : :dtriangle,
                     color = o.type == -1 ? (:teal, 0.6) : (:sienna, 0.6), markersize = 8)
        end
    end
    xlims!(ax, 0, xl)
    j == 1 && axislegend(ax, position = :rt, framevisible = false)
end
Label(fig[0, :], "Antarctic sea-level data: index points (blue, 2 sigma), marine limiting (teal: sea level was above), terrestrial limiting " *
      "(brown: below)", fontsize = 14)
save(joinpath(out, "ant_sites.png"), fig)

# --- ANT thickness at 20 ka -------------------------------------------------------------------------------------------------------------
H20(f; years_ago = false) = NCDataset(f) do ds
    t = Float64.(ds["time"][:]); k = findfirst(==(years_ago ? 20000.0 : -20000.0), t)
    (Float64.(ds["xc"][:]), Float64.(ds["yc"][:]), Float64.(nomissing(ds["H_ice"][:, :, k], 0.0)))
end
x, y, Hp = H20(joinpath(fes, "Antarctica", "ANT-16KM", "ANT-16KM_IceHistory-PaleoMIST.nc"); years_ago = true)
_, _, H0 = H20(joinpath(prior_dir, "ice_ant.nc")); _, _, H1 = H20(joinpath(post_dir, "ice_ant.nc"))
fig = Figure(size = (1600, 450))
nanz(H) = replace(H, 0.0 => NaN)
for (j, (t, H, cm, cr)) in enumerate([("PaleoMIST", nanz(Hp), :ice, (0, 4500)), ("prior mean", nanz(H0), :ice, (0, 4500)),
                                      ("posterior mean", nanz(H1), :ice, (0, 4500)), ("posterior - prior", replace(H1 .- H0, 0.0 => NaN), :balance, (-1000, 1000))])
    ax = Axis(fig[1, 2j - 1], title = "ANT-16KM H at 20 ka: $t", aspect = DataAspect())
    hidedecorations!(ax)
    hm = heatmap!(ax, x, y, H, colormap = cm, colorrange = cr, nan_color = :gray92)
    Colorbar(fig[1, 2j], hm, label = "m")
end
save(joinpath(out, "ant_maps.png"), fig)
@printf("ANT 20 ka: mean H on the ice prior %.0f m, posterior %.0f m; area prior %.2f, posterior %.2f (1e6 km2)\n",
        mean(H0[H0 .> 0]), mean(H1[H1 .> 0]), count(H0 .> 0) * (x[2] - x[1])^2 / 1e6, count(H1 .> 0) * (x[2] - x[1])^2 / 1e6)
println("wrote $out")
