# Maps of a kryos plastic run, after Gowan et al. (2021, Fig. 2) (docs/dev/plastic-backend.md; drawing in plastic_mapkit.jl):
#   julia --project=tools/levante tools/levante/plastic_lgm_map.jl <run folder> <out.png> [mode=age|margins|smb] [age=20000] [window=2500]
#         [ages=25000,20000,15000,12500,10000,7500] [label=<text>] [lon0=0]
# mode = age (default): the run at one age: the ice surface elevation above the sea level of that age (the plastic instances on their
#   grids, needs their 2-D output; the rest of the world's ice from VILMA2's load), the relative sea level of that age against today over
#   the ocean of that age (filled and contoured), the coastlines of that age and of today, the target margin (PaleoMIST's, prescribed),
#   and the GAPSLIP 2.1 sites with index points of ages within age +- window, coloured by the RMS of z_m of those points.
# mode = margins: the target margins at the ages (isochrones), the coastline of today, and every site with index points coloured by the
#   RMS of z_m over all its points (0-80 kyr).
# mode = smb: the implied mass balance at the age (mass_balance = "implied"): smb on the ice [m w.e. yr-1] (the land-margin
#   ablation in the margin cells), the totals of smb and cmb per instance [Gt yr-1], the target margin and the coastlines.
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
elseif mode == "smb"
    # the implied mass balance at the age (mass_balance = "implied"): smb on the ice [m w.e. yr-1], the margin cells' land ablation
    # included; the totals [Gt yr-1] of smb and cmb (calving) per instance
    G = read_global(run, age); R = read_regional(run, age)
    all(r -> all(isnan, r.smb), R) && error("$(run): no smb in the 2-D ice output (mass_balance = \"implied\", out_vars)")
    tot(r, F) = sum(F[isfinite.(F)] .* r.area[isfinite.(F)]) * 1e-12
    txt = join([@sprintf("%s: SMB %.0f, calving %.0f", r.inst, tot(r, r.smb), tot(r, r.cmb)) for r in R], ";  ")
    fig = Figure(size = (2300, 1150))
    ax = mapaxis(fig[1, 1], title = "$label: implied mass balance at $(kyr(age))", titlesize = 34,
                 subtitle = "[Gt yr-1] " * txt, subtitlesize = 22)
    graticule!(ax)
    SMBMAP = cgrad([RGBf(0.70, 0.13, 0.15), RGBf(0.97, 0.97, 0.97), RGBf(0.13, 0.35, 0.65)])
    sms = [meshfield!(ax, r.P, r.smb ./ 1000.0, r.H .> 1, colormap = SMBMAP, colorrange = (-1, 1),   # m w.e. yr-1
                      lowclip = RGBf(0.4, 0.0, 0.05), highclip = RGBf(0.03, 0.15, 0.40)) for r in R]
    sm = something(sms...)
    coasts!(ax, G); margins!(ax, R); frame!(ax)
    Colorbar(fig[1, 2][1, 1], sm, label = "implied SMB [m w.e. yr-1]", height = Relative(0.85), labelsize = 24)
    elems = [LineElement(color = MARGIN, linewidth = 2.5), LineElement(color = :black, linewidth = 1.5), LineElement(color = :gray35)]
    axislegend(ax, elems, ["target ice margin", "coast at $(kyr(age))", "coast today"], position = :lb, labelsize = 20,
               framevisible = true, backgroundcolor = (:white, 0.85))
else
    error("mode: age, margins or smb")
end
colsize!(fig.layout, 2, Auto(0.08))
save(outpng, fig, px_per_unit = 1.3)
println("wrote ", outpng)
