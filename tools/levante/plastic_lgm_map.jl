# Maps of a kryos plastic run, after Gowan et al. (2021, Fig. 2) (docs/dev/plastic-backend.md; drawing in plastic_mapkit.jl):
#   julia --project=<julia-plot> tools/levante/plastic_lgm_map.jl <run folder> <out.png> [mode=age|margins] [age=20000] [window=2500]
#         [ages=25000,20000,15000,12500,10000,7500] [label=<text>] [lon0=0]
# mode = age (default): the run at one age: the ice surface elevation above the sea level of that age (the plastic instances on their
#   grids, needs their 2-D output; the rest of the world's ice from VILMA2's load), the relative sea level of that age against today over
#   the ocean of that age (filled and contoured), the coastlines of that age and of today, the target margin (PaleoMIST's, prescribed),
#   and the GAPSLIP 2.1 sites with index points of ages within age +- window, coloured by the RMS of z_m of those points.
# mode = margins: the target margins at the ages (isochrones), the coastline of today, and every site with index points coloured by the
#   RMS of z_m over all its points (0-80 kyr).
# z_m = d / sqrt(sigma^2 + sigma_m^2), sigma_m 5 m. Triangles: the ice-sheet groups (Laurentide, Greenland, Fennoscandia/Barents,
# Antarctica); circles: the far field.
include(joinpath(@__DIR__, "plastic_mapkit.jl"))

run, outpng = abspath(ARGS[1]), ARGS[2]
opt = Dict(split(a, "=", limit = 2) for a in ARGS[3:end])
mode = get(opt, "mode", "age")
age = parse(Float64, get(opt, "age", "20000")); window = parse(Float64, get(opt, "window", "2500"))
LON0[] = parse(Float64, get(opt, "lon0", "0"))
label = get(opt, "label", basename(run))
set_theme!(fontsize = 24)

function cbars!(fig, items)
    for (k, (src, lbl, kw)) in enumerate(items)
        Colorbar(fig[1, 2][k, 1], src...; label = lbl, height = Relative(0.85), labelsize = 24, kw...)
    end
end

if mode == "age"
    G = read_global(run, age); R = read_regional(run, age)
    ss, ip, sc = sites((age - window, age + window); run = run)
    rmsall = isempty(ip) ? NaN : sqrt(mean(zm(sc[i])^2 for i in ip))
    fig = Figure(size = (2300, 1150))
    ax = mapaxis(fig[1, 1], title = "$label: $(kyr(age))", titlesize = 34,
                 subtitle = @sprintf("ESL %.1f m;  index points of %s to %s: %d at %d sites, RMS z_m %.2f", G.esl, kyr(age - window),
                                     kyr(age + window), length(ip), length(ss), rmsall), subtitlesize = 24)
    v = G.rsl[(G.ocn .> 0.5) .& (abs.(G.LAT) .< 45)]
    lo, hi = floor(quantile(v, 0.01) / 10) * 10, ceil(quantile(v, 0.99) / 10) * 10
    om = meshfield!(ax, G.P, G.rsl, G.ocn .> 0.5, colormap = RSLMAP, colorrange = (lo, hi), lowclip = RSL_LOW, highclip = RSL_HIGH)
    graticule!(ax)
    ocean_contours!(ax, G, G.rsl, lo:10:hi)
    im = ice!(ax, G, R)
    coasts!(ax, G); margins!(ax, R); frame!(ax)
    sites!(ax, ss)
    cbars!(fig, [((im,), "ice surface elevation [m]", ()), ((om,), "RSL at $(kyr(age)) against today [m]", ()),
                 ((), "site RMS of z_m (± $(kyr(window)))",
                  (colormap = SITEMAP, limits = (0, SITE_MAX), highclip = :darkred))])
    elems = [MarkerElement(marker = :utriangle, color = :white, strokecolor = :black, strokewidth = 0.9, markersize = 18),
             MarkerElement(marker = :circle, color = :white, strokecolor = :black, strokewidth = 0.9, markersize = 18),
             LineElement(color = MARGIN, linewidth = 2.5), LineElement(color = :black, linewidth = 1.5), LineElement(color = :gray35)]
    axislegend(ax, elems, ["ice-sheet sites", "far-field sites", "target ice margin", "coast at $(kyr(age))", "coast today"],
               position = :lb, labelsize = 20, framevisible = true, backgroundcolor = (:white, 0.85))
elseif mode == "margins"
    iso = parse.(Float64, split(get(opt, "ages", "25000,20000,15000,12500,10000,7500"), ","))
    G0 = read_global(run, 0.0)
    ss, ip, sc = sites((0.0, 80000.0); run = run)
    fig = Figure(size = (2300, 1150))
    ax = mapaxis(fig[1, 1], title = "$label: target margins and the misfit over all ages", titlesize = 34,
                 subtitle = @sprintf("%d index points at %d sites, 0 to 80 kyr; RMS z_m %.2f", length(ip), length(ss),
                                     sqrt(mean(zm(sc[i])^2 for i in ip))), subtitlesize = 24)
    graticule!(ax)
    for (xs, ys) in clines(G0.dl, G0.lat, G0.ocn, 0.5)
        lines!(ax, seamlines([eqe(x, y) for (x, y) in zip(xs, ys)]), color = :gray35, linewidth = 0.9)
    end
    isocm = cgrad(cgrad(:viridis)[range(0, 0.8, length = length(iso))], categorical = true)
    for (k, a) in enumerate(iso)
        margins!(ax, read_regional(run, a), color = isocm[k], linewidth = 2.0)
    end
    frame!(ax)
    sites!(ax, ss)
    Colorbar(fig[1, 2][1, 1], colormap = isocm, limits = (0.5, length(iso) + 0.5), ticks = (1:length(iso), kyr.(iso)),
             label = "target margin", height = Relative(0.85), labelsize = 24)
    Colorbar(fig[1, 2][2, 1], colormap = SITEMAP, limits = (0, SITE_MAX), highclip = :darkred, label = "site RMS of z_m, all ages",
             height = Relative(0.85), labelsize = 24)
    elems = [MarkerElement(marker = :utriangle, color = :white, strokecolor = :black, strokewidth = 0.9, markersize = 18),
             MarkerElement(marker = :circle, color = :white, strokecolor = :black, strokewidth = 0.9, markersize = 18)]
    axislegend(ax, elems, ["ice-sheet sites", "far-field sites"], position = :lb, labelsize = 20, framevisible = true,
               backgroundcolor = (:white, 0.85))
else
    error("mode: age or margins")
end
colsize!(fig.layout, 2, Auto(0.08))
save(outpng, fig, px_per_unit = 1.3)
println("wrote ", outpng)
