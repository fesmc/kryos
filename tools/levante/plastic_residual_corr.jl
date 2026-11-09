# The structure of the misfit of the GAPSLIP 2.1 index points, to choose the model-error part of the observation error of the
# calibration (plastic_esmda.jl: R = diag(sigma_i^2) + sigma_m^2 C, C_ij = exp(-dist / L_x) exp(-|dt| / L_t); v2: sigma_m 5 m,
# L_x 150 km, L_t 2 kyr) (docs/dev/plastic-backend.md).
#   julia --project=<julia-plot> tools/levante/plastic_residual_corr.jl <out folder> <label>=<run folder> ...
# Per run, the misfit d = model - obs of the index points (pm_compare.jl) and their 1 sigma (the scorer's). Within each group of
# pm_compare.jl, each GAPSLIP region and each sector: the correlation of d over pairs of index points, rho = mean(d_i d_j) / sqrt(mean(d_i^2)
# mean(d_j^2)) (not centred: a common bias is a correlated error too), against the distance (pairs with |dt| <= 1 kyr) and against
# the age difference (pairs within 25 km); an exponential c exp(-h / L) fitted to the bins (weighted by the pairs); and the error
# beyond the observations, sigma_m^2 = mean(d^2) - mean(sigma_i^2). Writes to <out>:
#   corr_distance.png, corr_time.png   rho by group, the runs, the v2 assumption exp(-h / 150 km), exp(-dt / 2 kyr)
#   corr_map.png                       per GAPSLIP sector: fitted L_x (colour) and sigma_m (size) at its mean position
#   residual_corr.md                   the tables
using CairoMakie, NCDatasets, Printf, Statistics

module PM
include("pm_compare.jl")
end
module KO
include("pm_kryos_out.jl")
end

const FES = "/work/ba1442/fesmdata"
const DBINS = [0.0, 10, 25, 50, 100, 200, 400, 800, 1600, 3200, 6400]
const TBINS = [0.0, 0.5, 1, 2, 4, 8, 16]                 # kyr
const DT_MAX, DX_MAX = 1000.0, 25.0                     # yr for the distance bins, km for the time bins
const LX_V2, LT_V2 = 150.0, 2.0

out = ARGS[1]; mkpath(out)
runs = [(String(split(a, "=")[1]), abspath(String(split(a, "=", limit = 2)[2]))) for a in ARGS[2:end]]
obs = PM.read_indicators(FES)
A = sort(NCDataset(ds -> Float64.(ds["time"][:]), joinpath(FES, "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc")))

function residuals(D)
    isfile(joinpath(D, "vilma_out.nc")) || KO.kryos_vilma_out(D, D)
    S = PM.series_run(PM.read_run(D), [o.lon for o in obs], [o.lat for o in obs], A)
    sc = [PM.score(o, A, S[i, :]) for (i, o) in enumerate(obs)]
    ip = [i for (i, o) in enumerate(obs) if o.type == 0 && isfinite(sc[i].z) && sc[i].z != 0]
    (ip = ip, d = Dict(i => sc[i].d for i in ip), s = Dict(i => sc[i].d / sc[i].z for i in ip))
end
gcd(a, b) = 6371.0 * 2asin(min(1.0, sqrt(sind((b.lat - a.lat) / 2)^2 + cosd(a.lat) * cosd(b.lat) * sind((b.lon - a.lon) / 2)^2)))

# binned correlation of d over the pairs of idx: by distance (|dt| <= DT_MAX) and by age difference (dist <= DX_MAX)
function corr(idx, d)
    nb = length(DBINS) - 1; nt = length(TBINS) - 1
    sx = zeros(nb); sxx = zeros(nb); syy = zeros(nb); nx = zeros(Int, nb)
    st = zeros(nt); stt = zeros(nt); suu = zeros(nt); ntt = zeros(Int, nt)
    for a in eachindex(idx), b in a+1:length(idx)
        i, j = idx[a], idx[b]; h = gcd(obs[i], obs[j]); dt = abs(obs[i].t - obs[j].t)
        if dt <= DT_MAX && h < DBINS[end]
            k = searchsortedlast(DBINS, h); sx[k] += d[i] * d[j]; sxx[k] += d[i]^2; syy[k] += d[j]^2; nx[k] += 1
        end
        if h <= DX_MAX && dt / 1000 < TBINS[end]
            k = searchsortedlast(TBINS, dt / 1000); st[k] += d[i] * d[j]; stt[k] += d[i]^2; suu[k] += d[j]^2; ntt[k] += 1
        end
    end
    ρ(s, a, b, n) = [n[k] >= 10 ? s[k] / sqrt(a[k] * b[k]) : NaN for k in eachindex(n)]
    (rx = ρ(sx, sxx, syy, nx), nx = nx, rt = ρ(st, stt, suu, ntt), nt = ntt)
end
mid(b) = [k == 1 ? b[2] / 2 : sqrt(b[k] * b[k+1]) for k in 1:length(b)-1]   # bin centres (the first: half its upper edge)
# c exp(-h / L) by least squares over a grid of L (c in closed form), bins weighted by their pairs
function fitexp(h, r, n)
    k = findall(i -> isfinite(r[i]), eachindex(r)); isempty(k) && return (NaN, NaN)
    best = (Inf, NaN, NaN)
    for L in exp.(range(log(1.0), log(1e5), length = 400))
        e = exp.(-h[k] ./ L); w = n[k]; c = clamp(sum(w .* e .* r[k]) / sum(w .* e .^ 2), 0, 1)
        s = sum(w .* (r[k] .- c .* e) .^ 2); s < best[1] && (best = (s, c, L))
    end
    best[2], best[3]
end

R = [residuals(D) for (_, D) in runs]
groups = [first.(PM.GROUPS); "all"]
gidx(g, r) = g == "all" ? r.ip : [i for i in r.ip if obs[i].group == g]
hx = mid(DBINS); ht = mid(TBINS)
C = Dict((g, m) => corr(gidx(g, R[m]), R[m].d) for g in groups, m in eachindex(runs))

set_theme!(fontsize = 20)
cols = Makie.wong_colors()
for (kind, hs, xl, ref) in (("distance", hx, "distance [km] (pairs with |dt| <= 1 kyr)", h -> exp(-h / LX_V2)),
                           ("time", ht, "age difference [kyr] (pairs within 25 km)", h -> exp(-h / LT_V2)))
    fig = Figure(size = (2000, 420))
    for (c, g) in enumerate(groups)
        ax = Axis(fig[1, c], title = g, xlabel = xl, ylabel = c == 1 ? "correlation of d" : "", xscale = kind == "distance" ? log10 : identity,
                  xlabelsize = 15, limits = (nothing, (-0.3, 1.05)))
        hlines!(ax, [0.0], color = :gray70)
        hh = kind == "distance" ? exp.(range(log(5), log(6400), length = 100)) : range(0, 16, length = 100)
        lines!(ax, hh, ref.(hh), color = :black, linestyle = :dash, label = "v2 assumption")
        for m in eachindex(runs)
            cc = C[(g, m)]; r = kind == "distance" ? cc.rx : cc.rt
            scatterlines!(ax, hs, r, color = cols[m], label = runs[m][1], markersize = 9)
        end
        c == 1 && axislegend(ax, position = :rt, labelsize = 14)
    end
    save(joinpath(out, "corr_$(kind).png"), fig)
end

# per GAPSLIP region and per sector (at least 20 index points), the last run; the position: the mean on the sphere
function spheremean(idx)
    v = sum([cosd(obs[i].lat) * cosd(obs[i].lon), cosd(obs[i].lat) * sind(obs[i].lon), sind(obs[i].lat)] for i in idx)
    atand(v[2], v[1]), atand(v[3], hypot(v[1], v[2]))
end
m = length(runs); rr = R[m]
function unitrows(key)
    rows = []
    for u in sort(unique(getfield(obs[i], key) for i in rr.ip))
        idx = [i for i in rr.ip if getfield(obs[i], key) == u]; length(idx) < 20 && continue
        cc = corr(idx, rr.d); cx, Lx = fitexp(hx, cc.rx, cc.nx); ct, Lt = fitexp(ht, cc.rt, cc.nt)
        sm2 = mean(rr.d[i]^2 for i in idx) - mean(rr.s[i]^2 for i in idx); lo, la = spheremean(idx)
        push!(rows, (reg = u, n = length(idx), lon = lo, lat = la, cx = cx, Lx = Lx, ct = ct, Lt = Lt,
                     rmsd = sqrt(mean(rr.d[i]^2 for i in idx)), rmss = sqrt(mean(rr.s[i]^2 for i in idx)), sm = sm2 > 0 ? sqrt(sm2) : 0.0))
    end
    rows
end
rows_reg = unitrows(:region); rows = unitrows(:sector)
lon5, lat5, zb = NCDataset(joinpath(FES, "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc")) do ds
    t = Float64.(ds["time"][:]); k = argmin(abs.(t)); lo = Float64.(ds["lon"][:]); la = Float64.(ds["lat"][:])
    lo = ifelse.(lo .> 180, lo .- 360, lo); p = sortperm(lo); q = sortperm(la)
    lo[p], la[q], Float64.(nomissing(ds["z_bed"][:, :, k], NaN))[p, q]
end
fig = Figure(size = (1800, 950))
ax = Axis(fig[1, 1], title = "Misfit of the index points by GAPSLIP sector ($(runs[m][1])): correlation length L_x (colour), " *
          "error beyond the observations sigma_m (size)", aspect = DataAspect(), limits = (-180, 180, -80, 85))
contour!(ax, lon5, lat5, zb, levels = [0.0], color = :gray50, linewidth = 0.7)
sc = scatter!(ax, [r.lon for r in rows], [r.lat for r in rows], color = log10.(clamp.([r.Lx for r in rows], 1, 1e5)), colormap = :viridis,
              colorrange = (1, 4), markersize = [12 + 2.5 * min(r.sm, 30) for r in rows], strokewidth = 1, strokecolor = :black)
for r in rows
    text!(ax, r.lon, r.lat, text = @sprintf("%s\n%.0f km, %.0f m", replace(r.reg, "_" => " "), r.Lx, r.sm), fontsize = 10,
          align = (:center, :bottom), offset = (0, 14))
end
Colorbar(fig[1, 2], sc, label = "log10 L_x [km]")
save(joinpath(out, "corr_map.png"), fig)

open(joinpath(out, "residual_corr.md"), "w") do io
    println(io, "Runs: ", join(["$l = `$D`" for (l, D) in runs], ", "), "\n")
    println(io, "### By group: fitted c exp(-h / L) of the correlation of d, and the error beyond the observations\n")
    println(io, "| group | run | n IP | RMS d [m] | RMS sigma_i [m] | sigma_m [m] | c_x | L_x [km] | c_t | L_t [kyr] |\n|---|---|---|---|---|---|---|---|---|---|")
    for g in groups, (mi, (l, _)) in enumerate(runs)
        idx = gidx(g, R[mi]); cc = C[(g, mi)]; cx, Lx = fitexp(hx, cc.rx, cc.nx); ct, Lt = fitexp(ht, cc.rt, cc.nt)
        rd = sqrt(mean(R[mi].d[i]^2 for i in idx)); rs = sqrt(mean(R[mi].s[i]^2 for i in idx))
        @printf(io, "| %s | %s | %d | %.1f | %.1f | %.1f | %.2f | %.0f | %.2f | %.1f |\n", mi == 1 ? g : "", l, length(idx), rd, rs,
                sqrt(max(rd^2 - rs^2, 0)), cx, Lx, ct, Lt)
    end
    for (unit, rs) in (("region", rows_reg), ("sector", rows))
        println(io, "\n### By GAPSLIP $unit ($(runs[m][1]); at least 20 index points)\n")
        println(io, "| $unit | n IP | RMS d [m] | RMS sigma_i [m] | sigma_m [m] | c_x | L_x [km] | c_t | L_t [kyr] |\n|---|---|---|---|---|---|---|---|---|")
        for r in rs
            @printf(io, "| %s | %d | %.1f | %.1f | %.1f | %.2f | %.0f | %.2f | %.1f |\n", r.reg, r.n, r.rmsd, r.rmss, r.sm, r.cx, r.Lx, r.ct, r.Lt)
        end
    end
    println(io, "\nBins with fewer than 10 pairs are left out; L at the ends of the search range (1 km, 1e5 km) means no decay resolved.")
end
println("wrote ", out)
