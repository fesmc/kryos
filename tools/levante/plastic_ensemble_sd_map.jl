# The spread of an ensemble of kryos plastic runs at one age on a global map (drawing in plastic_mapkit.jl): the standard deviation
# across the members of the relative sea level against today (VILMA2, over the ensemble-mean ocean of the age) and of the ice thickness
# (the plastic instances, on their grids, where every member has ice: the margin is prescribed), with the GAPSLIP 2.1 sites that have
# index points of ages within age +- window, the target margin and the coastlines of the age (the ensemble mean ocean).
#   julia --project=tools/levante tools/levante/plastic_ensemble_sd_map.jl <out.png> <member folder> ... [age=20000] [window=2500]
#         [label=<text>] [lon0=0]
# A member folder: a run with status 0, geo_glob_vilma.nc and (for the ice) <instance>.nc with 2-D output at the age.
include(joinpath(@__DIR__, "plastic_mapkit.jl"))

args = [a for a in ARGS[2:end] if !occursin("=", a)]
opt = Dict(split(a, "=", limit = 2) for a in ARGS[2:end] if occursin("=", a))
outpng = ARGS[1]
age = parse(Float64, get(opt, "age", "20000")); window = parse(Float64, get(opt, "window", "2500"))
LON0[] = parse(Float64, get(opt, "lon0", "0"))
label = get(opt, "label", "ensemble")
ok(D) = isfile(joinpath(D, "status")) && split(read(joinpath(D, "status"), String))[1] == "0"
mems = [abspath(D) for D in args if ok(D)]
length(mems) >= 3 || error("only $(length(mems)) finished members")

Gs = [read_global(D, age) for D in mems]
Rs = [read_regional(D, age; need = false) for D in mems]
sdof(fs) = (A = cat(fs..., dims = 3); dropdims(std(A, dims = 3), dims = 3))
meanof(fs) = (A = cat(fs..., dims = 3); dropdims(mean(A, dims = 3), dims = 3))
G1 = Gs[1]
G = merge(G1, (ocn = meanof([g.ocn for g in Gs]), h = meanof([g.h for g in Gs])))
rsl_sd = sdof([g.rsl for g in Gs])
esl = [g.esl for g in Gs]
# the ice thickness per instance, over the members that have it
regsd = []
for inst in unique(r.inst for R in Rs for r in R)
    rs = [r for R in Rs for r in R if r.inst == inst]
    length(rs) >= 3 || continue
    allice = reduce(.&, [r.H .> 1 for r in rs])
    sd = sdof([r.H for r in rs]); sd[.!allice] .= NaN
    push!(regsd, (r = rs[1], sd = sd, n = length(rs), allice = allice))
end
nice = isempty(regsd) ? 0 : minimum(x.n for x in regsd)

set_theme!(fontsize = 24)
fig = Figure(size = (2300, 1150))
ss, ip, _ = sites((age - window, age + window))
ax = mapaxis(fig[1, 1], title = "$label: spread at $(kyr(age))", titlesize = 34,
             subtitle = @sprintf("%d members (ice: %d);  ESL %.1f ± %.1f m;  %d index points of %s to %s at %d sites", length(mems), nice,
                                 mean(esl), std(esl), length(ip), kyr(age - window), kyr(age + window), length(ss)), subtitlesize = 24)
ocean = (G.h .<= 1) .& (G.ocn .> 0.5)
smax = ceil(quantile(filter(isfinite, rsl_sd[ocean]), 0.99))
SDMAP = cgrad([RGBf(0.97, 0.97, 0.99), RGBf(0.47, 0.40, 0.62)])
rm = meshfield!(ax, G.P, rsl_sd, ocean, colormap = SDMAP, colorrange = (0, smax), highclip = RGBf(0.38, 0.30, 0.52))
graticule!(ax)
ocean_contours!(ax, G, rsl_sd, 1:1:smax, color = (:black, 0.35))
hmax = isempty(regsd) ? 1.0 : ceil(quantile(filter(isfinite, vcat([vec(x.sd) for x in regsd]...)), 0.99) / 50) * 50
ims = [meshfield!(ax, x.r.P, x.sd, x.allice, colormap = :YlOrBr, colorrange = (0, hmax), highclip = RGBf(0.40, 0.15, 0.02))
       for x in regsd]
im = isempty(ims) ? nothing : something(ims...)
coasts!(ax, G); margins!(ax, [x.r for x in regsd]); frame!(ax)
sites!(ax, ss, color = :white, markersize = 17)
k = 1
if im !== nothing
    Colorbar(fig[1, 2][k, 1], im, label = "sd of ice thickness [m]", height = Relative(0.85), labelsize = 24); k += 1
end
Colorbar(fig[1, 2][k, 1], rm, label = "sd of RSL [m]", height = Relative(0.85), labelsize = 24)
elems = [MarkerElement(marker = :utriangle, color = :white, strokecolor = :black, strokewidth = 0.9, markersize = 18),
         MarkerElement(marker = :circle, color = :white, strokecolor = :black, strokewidth = 0.9, markersize = 18),
         LineElement(color = MARGIN, linewidth = 2.5), LineElement(color = :black, linewidth = 1.5)]
axislegend(ax, elems, ["ice-sheet sites", "far-field sites", "target ice margin",
                       "coast"], position = :lb, labelsize = 20, framevisible = true, backgroundcolor = (:white, 0.85))
colsize!(fig.layout, 2, Auto(0.08))
save(outpng, fig, px_per_unit = 1.3)
println("wrote ", outpng)
