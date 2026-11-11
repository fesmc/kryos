# The results of a tau_b calibration (plastic_esmda.jl): per-iteration scores, the factor curves, the ice volume of the mean
# members (NH-32KM and, when run, ANT-16KM), the present-day thickness of the last mean member per region against the observations
# (when it wrote ice_nor.nc / ice_ant.nc: final_set of the config; plastic_pd_regions.jl), and an offline Greenland check.
#   julia --project=tools/levante tools/levante/plastic_esmda_report.jl <calibration folder> <out folder> [<prior slices.nc> <posterior slices.nc>]
# Writes <out>/report.md and <out>/theta_factors.png. The slices files: tests/plastic_slices.f90 with <out.nc>, without and with the
# posterior-mean multiplier file (the last iteration's scale/m000.nc): Greenland at 0 ka against PaleoMIST and BedMachine v6
# (FesmData NH-32KM TOPO-GEBCO2026, its H_ice) over the Greenland region (REGIONS region_2 = 103).
using CairoMakie, NCDatasets, Printf, Statistics, TOML
include(joinpath(@__DIR__, "plastic_pd_regions.jl"))
module PM
include(joinpath(@__DIR__, "pm_compare.jl"))
end

const RHO_I, RHO_SW, A_OCN = 910.0, 1028.0, 3.618e14
C, out = ARGS[1], ARGS[2]; mkpath(out)
cfg = TOML.parsefile(joinpath(C, "config.toml"))
iters = sort([parse(Int, d[6:end]) for d in readdir(C) if startswith(d, "iter_")])
itd(k) = joinpath(C, @sprintf("iter_%02d", k))
# the mean member 0 of iteration k: 0/ (staged by runme), or m000/ in a folder of before 2026-11-11
m0(k) = (d = joinpath(itd(k), "0"); isdir(d) ? d : joinpath(itd(k), "m000"))
readcsv(f) = (l = readlines(f); h = split(l[1], ","); [Dict(zip(h, split(r, ","))) for r in l[2:end]])
num(x) = parse(Float64, x)
io = open(joinpath(out, "report.md"), "w")

# --- per iteration ----------------------------------------------------------------------------------------------------------------
sets = ["obs groups"; cfg["esmda"]["obs_groups"]]
println(io, "### Per iteration (scores over all index points; mean member m000; ensemble median [min, max])\n")
println(io, "| iteration | alpha | members ok | set | mean member RMS z | median abs(z) | ML viol. | TL viol. | ensemble RMS z | ESL 20 ka mean member [m] | ESL 20 ka ensemble [m] | sd(theta) mean |")
println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|")
α = cfg["esmda"]["alpha"]
for k in iters
    isfile(joinpath(itd(k), "scores.csv")) || continue
    sc = readcsv(joinpath(itd(k), "scores.csv")); mb = readcsv(joinpath(itd(k), "members.csv"))
    ok = [r for r in mb if r["ok"] == "1"]; ens = [r["member"] for r in ok if r["member"] != "0"]
    th = readlines(joinpath(itd(k), "theta.csv"))[2:end]
    cols = [num.(split(l, ",")[2:end]) for l in th if split(l, ",")[1] in ens]
    sdθ = length(cols) < 2 ? NaN : mean(std(reduce(hcat, cols); dims = 2))
    e0 = [num(r["esl_20ka"]) for r in ok if r["member"] == "0"]; ee = [num(r["esl_20ka"]) for r in ok if r["member"] != "0"]
    for (j, s) in enumerate(sets)
        r0 = [r for r in sc if r["member"] == "0" && r["set"] == s]
        re = [num(r["rms_z"]) for r in sc if r["member"] in ens && r["set"] == s]
        f(r, key) = isempty(r) ? "-" : @sprintf("%.2f", num(r[1][key]))
        println(io, "| ", j == 1 ? "$k" : "", " | ", j == 1 ? (k <= length(α) ? "$(α[k])" : "eval") : "", " | ",
                j == 1 ? "$(length(ok)) of $(length(mb))" : "", " | ", s, " | ", f(r0, "rms_z"), " | ", f(r0, "median_abs_z"), " | ",
                f(r0, "ml_violated"), " | ", f(r0, "tl_violated"), " | ",
                isempty(re) ? "-" : @sprintf("%.2f [%.2f, %.2f]", median(re), extrema(re)...), " | ",
                j == 1 && !isempty(e0) ? @sprintf("%.1f", e0[1]) : "", " | ",
                j == 1 && !isempty(ee) ? @sprintf("%.1f [%.1f, %.1f]", median(ee), extrema(ee)...) : "", " | ",
                j == 1 ? @sprintf("%.3f", sdθ) : "", " |")
    end
end

# --- scores with the model error: z_m = d / sqrt(sigma_i^2 + sigma_m^2) (the diagonal of R), robust to an index point whose scorer sigma
# is ~0 (a one-sided elevation range on the side of d and a flat model: z of the scores unbounded) ----------------------------------------
σm = get(cfg["esmda"], "sigma_m", get(cfg["esmda"], "sigma_add", 0.0))
if σm > 0
    obs = [o for o in PM.read_indicators(cfg["score"]["fesmdata"]) if o.type == 0]
    A = sort(NCDataset(ds -> Float64.(ds["time"][:]),
                       joinpath(cfg["score"]["fesmdata"], "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc")))
    gs = cfg["esmda"]["obs_groups"]
    println(io, "\n### Mean members: RMS of z_m = d / sqrt(sigma^2 + sigma_m^2), sigma_m = $(σm) m (index points; n with scorer sigma < 0.1 m)\n")
    println(io, "| iteration | all | ", join(gs, " | "), " | sigma < 0.1 m |\n|---|---|", repeat("---|", length(gs) + 1))
    for k in iters
        D = m0(k); isfile(joinpath(D, "geo_glob_vilma.nc")) || continue
        r = PM.read_run(D); S = PM.series_run(r, [o.lon for o in obs], [o.lat for o in obs], A)
        zm = fill(NaN, length(obs)); small = 0
        for (i, o) in enumerate(obs)
            any(isnan, S[i, :]) && continue
            m, dmdt = PM.tval(A, S[i, :], o.t); d = m - o.rsl; s = PM.ip_sigma(o, d, dmdt)
            zm[i] = d / sqrt(s^2 + σm^2); small += s < 0.1
        end
        rz(g) = (z = [zm[i] for i in eachindex(obs) if (g === nothing || obs[i].group == g) && isfinite(zm[i])]; @sprintf("%.2f", sqrt(mean(z .^ 2))))
        println(io, "| ", k, " | ", rz(nothing), " | ", join([rz(g) for g in gs], " | "), " | ", small, " |")
    end
end

# --- factor curves --------------------------------------------------------------------------------------------------------------
hdr = split(readlines(joinpath(itd(1), "theta.csv"))[1], ",")[2:end]
pg = [String(rsplit(h, "_"; limit = 2)[1]) for h in hdr]
pa = [parse.(Float64, split(rsplit(h, "_"; limit = 2)[2], "-")) for h in hdr]
groups = [g["name"] for g in cfg["scale"]["groups"]]; nodes = Float64.(cfg["scale"]["nodes_ka"])
function curve(θ, g)          # factor at the nodes of group g, in the order of nodes (1 at a pinned node)
    [(q = findfirst(q -> pg[q] == g && a in pa[q], eachindex(hdr)); q === nothing ? 1.0 : exp(θ[q])) for a in nodes]
end
function ensemble(k)
    th = readlines(joinpath(itd(k), "theta.csv"))[2:end]
    rows = [num.(split(l, ",")) for l in th]
    m0 = rows[findfirst(r -> r[1] == 0, rows)][2:end]
    E = [r[2:end] for r in rows if r[1] != 0]
    m0, E
end
kl = maximum(iters); pr = cfg["prior"]
psd(t) = (r = get(pr, "ramp_ka", 0.0); r > 0 && t < r ? pr["sigma"] * sin(π * t / (2r))^2 : float(pr["sigma"]))   # plastic_esmda.jl
fig = Figure(size = (320 * length(groups) + 80, 400))
order = sortperm(nodes)
for (c, g) in enumerate(groups)
    ax = Axis(fig[1, c]; title = g, xlabel = "age [ka]", ylabel = c == 1 ? "tau_b factor" : "", yscale = log10,
              xreversed = true, yticks = [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4])
    tt = range(minimum(nodes), maximum(nodes); length = 200)
    band!(ax, tt, exp.(-psd.(tt)), exp.(psd.(tt)); color = (:grey, 0.2))
    hlines!(ax, [1.0]; color = :black, linewidth = 0.8)
    m0, E = ensemble(kl)
    if !isempty(E)
        F = reduce(hcat, [log.(curve(e, g)) for e in E])
        μ = vec(mean(F; dims = 2)); s = vec(std(F; dims = 2))
        band!(ax, nodes[order], exp.(μ .- s)[order], exp.(μ .+ s)[order]; color = (:tomato, 0.25))
    end
    for k in iters[2:end]
        m, _ = ensemble(k)
        lines!(ax, nodes[order], curve(m, g)[order]; color = k == kl ? :firebrick : (:firebrick, 0.25 + 0.5 * k / kl),
               linewidth = k == kl ? 2.5 : 1, label = k == kl ? "posterior mean" : "mean after update $(k - 1)")
    end
    ylims!(ax, 0.2, 4.6)
    c == 1 && axislegend(ax; position = :lt, framevisible = false, labelsize = 10)
end
save(joinpath(out, "theta_factors.png"), fig)
println(io, "\n### Posterior-mean factors (exp(theta) of m000 of iteration $kl)\n")
println(io, "| group | ", join([@sprintf("%g ka", a) for a in nodes], " | "), " |\n|---|", repeat("---|", length(nodes)))
m0, _ = ensemble(kl)
for g in groups; println(io, "| ", g, " | ", join([@sprintf("%.2f", f) for f in curve(m0, g)], " | "), " |"); end

# --- ice volume of the mean members ------------------------------------------------------------------------------------------------
fes = cfg["score"]["fesmdata"]
function pm_sle(dom, grid)             # PaleoMIST on the grid, m SLE at 20 and 0 ka; the cell area
    NCDataset(joinpath(fes, dom, grid, "$(grid)_IceHistory-PaleoMIST.nc")) do ds
        x = ds["xc"][:]; y = ds["yc"][:]; t = Float64.(ds["time"][:]); d = (x[2] - x[1]) * (y[2] - y[1]) * 1e6
        d, Dict(a => sum(nomissing(ds["H_ice"][:, :, findfirst(==(a), t)], 0.0)) * d * RHO_I / RHO_SW / A_OCN for a in (20000.0, 0.0))
    end
end
dxdy, slpm = pm_sle("North", "NH-32KM")
println(io, "\n### Ice volume (m SLE; M_ice of <instance>_ts.nc; ratio to PaleoMIST)\n\n| | grid | 20 ka | 0 ka |\n|---|---|---|---|")
for (inst, dom, grid) in (("ice_nor", "North", "NH-32KM"), ("ice_ant", "Antarctica", "ANT-16KM"))
    isfile(joinpath(m0(1), "$(inst)_ts.nc")) || continue
    _, sp = pm_sle(dom, grid)
    @printf(io, "| PaleoMIST | %s | %.2f | %.2f |\n", grid, sp[20000.0], sp[0.0])
    for k in (1, kl)
        f = joinpath(m0(k), "$(inst)_ts.nc"); isfile(f) || continue
        v = NCDataset(f) do ds
            t = Float64.(ds["time"][:]); M = Float64.(ds["M_ice"][:])
            Dict(a => M[argmin(abs.(t .+ a))] * 1e12 / RHO_SW / A_OCN for a in (20000.0, 0.0))
        end
        @printf(io, "| m000 of iteration %d%s | %s | %.2f (%.2f) | %.2f (%.2f) |\n", k, k == 1 ? " (prior mean)" : " (posterior mean)",
                grid, v[20000.0], v[20000.0] / sp[20000.0], v[0.0], v[0.0] / sp[0.0])
    end
end

# --- present day of the last mean member per region ----------------------------------------------------------------------------------
if isfile(joinpath(m0(kl), "ice_nor.nc"))
    println(io, "\n### Present-day thickness of m000 of iteration $kl against the observations\n")
    pd_run(io, m0(kl), fes; label = "posterior mean")
end

# --- cost -----------------------------------------------------------------------------------------------------------------------------
sec = Float64[]
for k in iters
    f = joinpath(itd(k), "members.csv"); isfile(f) || continue
    append!(sec, [num(r["seconds"]) for r in readcsv(f) if r["ok"] == "1"])
end
thr = cfg["run"]["threads"]; memgb = parse(Float64, replace(cfg["run"]["mem"], "G" => ""))
billed = max(thr, ceil(memgb / 0.9))
@printf(io, "\n### Cost\n\n%d runs, %.0f-%.0f min each (median %.0f), %.0f cpu-h billed for the members (%d cpus each)\n", length(sec),
        minimum(sec) / 60, maximum(sec) / 60, median(sec) / 60, sum(sec) / 3600 * billed, billed)

# --- Greenland 0 ka, offline -------------------------------------------------------------------------------------------------------
if length(ARGS) >= 4
    reg = NCDataset(ds -> Int.(nomissing(ds["region_2"][:, :], -1)), joinpath(fes, "North", "NH-32KM", "NH-32KM_REGIONS.nc"))
    bm = NCDataset(ds -> Float64.(nomissing(ds["H_ice"][:, :], 0.0)), joinpath(fes, "North", "NH-32KM", "NH-32KM_TOPO-GEBCO2026.nc"))
    sle(H, m) = sum(H[m]) * dxdy * RHO_I / RHO_SW / A_OCN
    println(io, "\n### Greenland at 0 ka (offline plastic_slices, PaleoMIST's 0 ka bed and margin; region_2 = 103)\n")
    println(io, "| | volume [m SLE] | mean H [m] on the plastic ice | dH vs PaleoMIST: mean, rms [m] | dH vs BedMachine: mean, rms [m] |\n|---|---|---|---|---|")
    first_ = Ref(true)
    for (lbl, f) in (("Gowan's tau_b", ARGS[3]), ("posterior mean", ARGS[4]))
        NCDataset(f) do ds
            t = Float64.(ds["time"][:]); k = argmin(abs.(t))
            H = Float64.(ds["H"][:, :, k]); P = Float64.(ds["H_pm"][:, :, k])
            g = reg .== 103; m = g .& (H .> 0)
            if first_[]
                @printf(io, "| PaleoMIST | %.2f | %.0f | | |\n", sle(P, g), mean(P[g .& (P .> 0)]))
                @printf(io, "| BedMachine v6 | %.2f | %.0f | | |\n", sle(bm, g), mean(bm[g .& (bm .> 0)]))
                first_[] = false
            end
            d1 = (H .- P)[m]; d2 = (H .- bm)[m]
            @printf(io, "| %s | %.2f | %.0f | %+.0f, %.0f | %+.0f, %.0f |\n", lbl, sle(H, g), mean(H[m]), mean(d1), sqrt(mean(d1 .^ 2)),
                    mean(d2), sqrt(mean(d2 .^ 2)))
        end
    end
end
close(io)
println("wrote ", out)
