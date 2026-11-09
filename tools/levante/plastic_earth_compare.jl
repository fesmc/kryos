# The plastic ice on different Earths (docs/dev/plastic-backend.md, "The 3-D Earth"): GAPSLIP 2.1 scores of kryos runs by group, by
# GAPSLIP region, by marine-bed area and by site, the change per site against the first run on maps, RSL curves at the best-sampled
# sites, ice volume per instance and ESL.
#   julia --project=<julia-plot> tools/levante/plastic_earth_compare.jl <out folder> <label>=<run folder> ...
# The runs: kryos folders with geo_glob_vilma.nc (vilma_out.nc is written from it, pm_kryos_out.jl) and <instance>_ts.nc. Index points
# as pm_compare.jl: d = model - obs, z = d / sigma; z_m = d / sqrt(sigma^2 + sigma_m^2) with the model error sigma_m (5 m, as
# plastic_esmda_report.jl). Writes to <out>:
#   compare.md       tables: volume and ESL; by group; by GAPSLIP region; by marine-bed area (pm_compare.jl FOCUS)
#   sites.csv        per site and run: n index points, RMS z, RMS z_m, mean d
#   sites_map.png    per site, RMS z_m of each run minus that of the first (blue: better than the first); Antarctica polar
#   curves.png       RSL at the sites with the most index points per group: index points (2 sigma) and the runs
using CairoMakie, NCDatasets, Printf, Statistics

module PM
include("pm_compare.jl")
end
module KO
include("pm_kryos_out.jl")
end

const FES = "/work/ba1442/fesmdata"
const SIGMA_M = 5.0
const RHO_SW, A_OCN = 1028.0, 3.618e14
const INSTANCES = ["ice_nor" => "NH", "ice_ant" => "ANT", "ice_pat" => "PAT"]
const NSITE = 4                     # curve panels per group
const MAP_BOXES = [("North America and Greenland", (-140.0, -10.0), (35.0, 84.0)), ("Europe", (-15.0, 70.0), (45.0, 82.0)),
                   ("Antarctica", (-180.0, 180.0), (-80.0, -60.0)), ("far field", (-180.0, 180.0), (-45.0, 40.0))]

out = ARGS[1]; mkpath(out)
runs = [(String(split(a, "=")[1]), abspath(String(split(a, "=", limit = 2)[2]))) for a in ARGS[2:end]]
length(runs) >= 2 || error("give at least two runs: <label>=<folder>")
for (_, D) in runs
    isfile(joinpath(D, "vilma_out.nc")) || KO.kryos_vilma_out(D, D)
end

obs = PM.read_indicators(FES)
A = sort(NCDataset(ds -> Float64.(ds["time"][:]), joinpath(FES, "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc")))
R = [PM.read_run(D) for (_, D) in runs]
S = [PM.series_run(r, [o.lon for o in obs], [o.lat for o in obs], A) for r in R]
SC = [[PM.score(o, A, s[i, :]) for (i, o) in enumerate(obs)] for s in S]
zm(sc) = isfinite(sc.z) ? (sc.z == 0 ? 0.0 : sc.d / sqrt((sc.d / sc.z)^2 + SIGMA_M^2)) : NaN
fmt(x) = isnan(x) ? "-" : @sprintf("%.2f", x)
rms(x) = isempty(x) ? NaN : sqrt(mean(x .^ 2))

function setrow(io, name, idx)
    for (m, (lbl, _)) in enumerate(runs)
        s = PM.summary(SC[m], obs, idx)
        ip = [i for i in idx if obs[i].type == 0 && isfinite(SC[m][i].z)]
        println(io, "| ", m == 1 ? name : "", " | ", lbl, " | ", s.nip, " | ", fmt(s.rmsz), " | ", fmt(rms([zm(SC[m][i]) for i in ip])), " | ",
                fmt(s.medz), " | ", isnan(s.bias) ? "-" : @sprintf("%.1f", s.bias), " | ", fmt(s.vml), " | ", fmt(s.vtl), " |")
    end
end
const HDR = "| set | run | n IP | RMS z | RMS z_m | median abs(z) | mean d [m] | ML violated | TL violated |\n|---|---|---|---|---|---|---|---|---|"

# --- volumes and ESL ----------------------------------------------------------------------------------------------------------------
sle(D, inst, a) = (f = joinpath(D, "$(inst)_ts.nc"); isfile(f) ? NCDataset(ds -> begin
    t = Float64.(ds["time"][:]); Float64(ds["M_ice"][argmin(abs.(t .+ a))]) * 1e12 / RHO_SW / A_OCN end, f) : NaN)
esl(r, a) = r.bsl_ive[PM.tindex(r, PM.trec(a))] - r.bsl_ive[end]

open(joinpath(out, "compare.md"), "w") do io
    println(io, "Runs: ", join(["$l = `$D`" for (l, D) in runs], ", "), "\n")
    println(io, "### Ice volume [m SLE] and ESL [m]\n")
    print(io, "| run |"); for (_, n) in INSTANCES, a in (20, 0); print(io, " $n $a ka |"); end
    println(io, " ESL 60 ka | ESL 20 ka |\n|---|", repeat("---|", 2 * length(INSTANCES) + 2))
    for (m, (l, D)) in enumerate(runs)
        print(io, "| ", l, " |")
        for (inst, _) in INSTANCES, a in (20000.0, 0.0)
            v = sle(D, inst, a); print(io, isnan(v) ? " - |" : @sprintf(" %.2f |", v))
        end
        @printf(io, " %.1f | %.1f |\n", esl(R[m], 60000.0), esl(R[m], 20000.0))
    end
    println(io, "\nGAPSLIP 2.1 (primary data, 0-80 ka): index points (IP), z = d / sigma, z_m = d / sqrt(sigma^2 + ", SIGMA_M,
            "^2), mean d = model - obs; limiting data violated (marine ML, terrestrial TL).")
    println(io, "\n### By group\n\n", HDR)
    setrow(io, "all", collect(eachindex(obs)))
    for (g, _) in PM.GROUPS; setrow(io, g, findall(o -> o.group == g, obs)); end
    println(io, "\n### By GAPSLIP region\n\n", HDR)
    for reg in sort(unique(o.region for o in obs)); setrow(io, reg, findall(o -> o.region == reg, obs)); end
    println(io, "\n### Marine-bed areas\n\n", HDR)
    for (f, (k, v)) in PM.FOCUS; setrow(io, f, findall(o -> getfield(o, k) in v, obs)); end
end

# --- per site -----------------------------------------------------------------------------------------------------------------------
sites = unique(o.site_name for o in obs)
sidx = Dict(s => findall(o -> o.site_name == s && o.type == 0, obs) for s in sites)
sites = [s for s in sites if !isempty(sidx[s])]
srms = [[rms(filter(isfinite, [zm(SC[m][i]) for i in sidx[s]])) for s in sites] for m in eachindex(runs)]
open(joinpath(out, "sites.csv"), "w") do io
    println(io, "site,group,region,lat,lon,n_ip", join([",rms_z_$l,rms_zm_$l,mean_d_$l" for (l, _) in runs]))
    for (k, s) in enumerate(sites)
        i = sidx[s]; o = obs[i[1]]
        print(io, s, ",", o.group, ",", o.region, ",", @sprintf("%.3f,%.3f", mean(obs[j].lat for j in i), mean(obs[j].lon for j in i)), ",", length(i))
        for m in eachindex(runs)
            z = filter(isfinite, [SC[m][j].z for j in i]); d = filter(isfinite, [SC[m][j].d for j in i])
            print(io, ",", fmt(rms(z)), ",", fmt(srms[m][k]), ",", isempty(d) ? "-" : @sprintf("%.1f", mean(d)))
        end
        println(io)
    end
end

# land outline: PaleoMIST's 0 ka bed on the 0.5 degree grid
lon5, lat5, zb = NCDataset(joinpath(FES, "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc")) do ds
    t = Float64.(ds["time"][:]); k = argmin(abs.(t))
    lo = Float64.(ds["lon"][:]); la = Float64.(ds["lat"][:]); Z = Float64.(nomissing(ds["z_bed"][:, :, k], NaN))
    lo = ifelse.(lo .> 180, lo .- 360, lo); p = sortperm(lo); q = sortperm(la)
    lo[p], la[q], Z[p, q]
end
slon = [(x = mean(obs[j].lon for j in sidx[s]); x > 180 ? x - 360 : x) for s in sites]
slat = [mean(obs[j].lat for j in sidx[s]) for s in sites]
nip = [length(sidx[s]) for s in sites]

# Antarctica in a polar view (azimuthal equidistant, degrees from the pole), the land as the 0.5 degree cells above sea level
polar(lo, la) = ((90 + la) * sind(lo), (90 + la) * cosd(lo))
land_s = [polar(x, y) for (i, x) in enumerate(lon5), (j, y) in enumerate(lat5) if y < -58 && zb[i, j] > 0]

fig = Figure(size = (1500, 820 * (length(runs) - 1)))
for m in 2:length(runs)
    dz = srms[m] .- srms[1]
    for (b, (title, (x0, x1), (y0, y1))) in enumerate(MAP_BOXES)
        pos = fig[2(m - 2) + (b - 1) ÷ 2 + 1, (b - 1) % 2 + 1]
        k = findall(i -> x0 <= slon[i] <= x1 && y0 <= slat[i] <= y1 && isfinite(dz[i]), eachindex(sites))
        if title == "Antarctica"
            ax = Axis(pos, title = "$(runs[m][1]) - $(runs[1][1]): $title", aspect = DataAspect(), titlesize = 12,
                      limits = (-32, 32, -32, 32))
            hidedecorations!(ax)
            scatter!(ax, first.(land_s), last.(land_s), color = :gray85, markersize = 3)
            px = [polar(slon[i], slat[i]) for i in k]; xs, ys = first.(px), last.(px)
        else
            ax = Axis(pos, title = "$(runs[m][1]) - $(runs[1][1]): $title", limits = (x0, x1, y0, y1), aspect = DataAspect(),
                      titlesize = 12)
            contour!(ax, lon5, lat5, zb, levels = [0.0], color = :gray40, linewidth = 0.6)
            xs, ys = slon[k], slat[k]
        end
        sc = scatter!(ax, xs, ys, color = dz[k], colormap = Reverse(:RdBu), colorrange = (-3, 3), markersize = 4 .+ 2 .* sqrt.(nip[k]),
                      strokewidth = 0.5, strokecolor = :black)
        b == length(MAP_BOXES) && Colorbar(fig[2(m - 2) + 1:2(m - 2) + 2, 3], sc,
                                           label = "RMS z_m of $(runs[m][1]) - $(runs[1][1]) (blue: better)")
        text!(ax, 0.02, 0.02, text = "$(count(<(0), dz[k])) better, $(count(>(0), dz[k])) worse", space = :relative, fontsize = 11)
    end
end
save(joinpath(out, "sites_map.png"), fig)

# --- RSL curves ---------------------------------------------------------------------------------------------------------------------
cols = Makie.wong_colors()
groups = first.(PM.GROUPS)
fig = Figure(size = (1600, 330 * length(groups)))
for (gi, g) in enumerate(groups)
    gs = sort([s for s in sites if obs[sidx[s][1]].group == g], by = s -> -length(sidx[s]))[1:min(NSITE, end)]
    for (si, s) in enumerate(gs)
        i = sidx[s]; la = mean(obs[j].lat for j in i); lo = mean(obs[j].lon for j in i)
        ax = Axis(fig[gi, si], title = "$g: $s ($(length(i)) IP)", xlabel = "age [ka]", ylabel = "RSL [m]", titlesize = 11,
                  xreversed = true)
        rangebars!(ax, [obs[j].t / 1000 for j in i], [obs[j].rmin for j in i], [obs[j].rmax for j in i], color = :gray50)
        rangebars!(ax, [obs[j].rsl for j in i], [obs[j].tmin / 1000 for j in i], [obs[j].tmax / 1000 for j in i], direction = :x, color = :gray70)
        amax = maximum(obs[j].tmax for j in i) * 1.3 + 2000
        k = findall(<=(amax), A)
        for (m, (l, _)) in enumerate(runs)
            ser = PM.series_run(R[m], [lo], [la], A)[1, :]
            lines!(ax, A[k] ./ 1000, ser[k], color = cols[m], linewidth = 2, label = l)
        end
        gi == 1 && si == 1 && axislegend(ax, position = :lt, labelsize = 10)
    end
end
save(joinpath(out, "curves.png"), fig)
println("wrote ", out)
