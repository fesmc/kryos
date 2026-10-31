# PaleoMIST test of the standalone VILMA2: figures from pm_compare.jl's output (CairoMakie).
#   julia --project=<julia-plot> tools/levante/pm_maps.jl <compare folder>
# pm_B1_vs_PM.png  : RSL of PaleoMIST, of B1 (lmax 256) and their difference at 22.5, 15, 10 ka
# pm_effects.png   : lmax 170 - 256 (B1) and 3-D - 1-D (B3 - B1, lmax 170) at the same records
# pm_decomp.png    : the 3-D effect split: lithosphere (B1_L80 - B1), radial profile (B1_bm - B1_L80), lateral (B3 - B1_bm), lmax 170
# pm_ind_sites.png : RSL relative to present at the best-constrained indicator sites: the indicators (index points: 2-sigma boxes,
#                    marine limiting: up triangles, terrestrial limiting: down triangles) and the models
# pm_ind_map.png   : per site, RMS of the normalized misfit of the index points, PaleoMIST and B1 (lmax 256)
# pm_A_minus_B.png : equilibrium - transient at 22.5 ka (A - B), RSL and u, 1-D (lmax 256) and 3-D
# pm_rsl0.png      : the present-day disequilibrium of B (rsl at 0 ka against the relaxed reference), 1-D and 3-D
# pm_sites.png     : RSL relative to present at the sites, B runs (every output step) and PaleoMIST (records)
# pm_gmsl.png      : VILMA2 bsl - bsl(0) and bsl_ive, PaleoMIST ocean mean of dz_sl and ESL
# RSL of the B runs relative to present (pm_compare.jl).
using NCDatasets, CairoMakie, DelimitedFiles, Printf

function gmap!(fig, pos, lon, lat, F, lim, title; cmap = :balance, ice = nothing)
    ax = Axis(fig[pos...], title = title, aspect = DataAspect(), titlesize = 13)
    hidedecorations!(ax)
    hm = heatmap!(ax, lon, lat, F; colormap = cmap, colorrange = (-lim, lim), rasterize = 2)
    ice === nothing || contour!(ax, lon, lat, ice; levels = [1.0], color = :black, linewidth = 0.5)
    limits!(ax, -180, 180, -90, 90)
    hm
end

const MODEL_STYLE = Dict("PaleoMIST" => (:black, :dash, 1.6), "B1_l256" => (:royalblue, :solid, 1.6), "B1_l170" => (:lightskyblue, :solid, 1.0),
                         "B1_L80" => (:darkorange, :solid, 1.4), "B1_bm" => (:seagreen, :solid, 1.4), "B3_l170" => (:firebrick, :solid, 1.6))

function indicator_figures(dir, lon, lat, land)
    c, _ = readdlm(joinpath(dir, "pm_ind_curves.csv"), ',', header = true)
    o, h = readdlm(joinpath(dir, "pm_ind_samples.csv"), ',', header = true)
    col(n) = o[:, findfirst(==(n), vec(h))]
    site = col("site"); t = Float64.(col("time")) ./ 1000; t0 = Float64.(col("time_min")) ./ 1000; t1 = Float64.(col("time_max")) ./ 1000
    r = Float64.(col("rsl")); r0 = Float64.(col("rsl_min")); r1 = Float64.(col("rsl_max")); ty = Int.(col("type"))
    sites = unique(c[:, 2]); srcs = unique(c[:, 1]); nc = 4
    fig = Figure(size = (1500, 300 * cld(length(sites), nc) + 60), fontsize = 11)
    local ax
    for (n, st) in enumerate(sites)
        m = site .== st
        grp = c[findfirst(==(st), c[:, 2]), 3]
        tmax = min(80.0, maximum(t1[m]) + 1.0)
        ax = Axis(fig[(n - 1) ÷ nc + 1, (n - 1) % nc + 1], title = "$(replace(st, '_' => ' ')) ($grp)", xlabel = "time [ka]", ylabel = "RSL [m]",
                  titlesize = 12)
        ip = findall(m .& (ty .== 0))
        isempty(ip) || poly!(ax, [Rect2f(-t1[i], r0[i], t1[i] - t0[i], r1[i] - r0[i]) for i in ip], color = (:gray, 0.35),
                             strokecolor = :gray30, strokewidth = 0.4)
        ml = findall(m .& (ty .== -1)); tl = findall(m .& (ty .== 1))
        isempty(ml) || scatter!(ax, -t[ml], r[ml], marker = :utriangle, color = :steelblue, markersize = 6)
        isempty(tl) || scatter!(ax, -t[tl], r[tl], marker = :dtriangle, color = :sienna, markersize = 6)
        lo, hi = extrema(vcat(r0[m], r1[m]))
        for src in srcs
            k = (c[:, 1] .== src) .& (c[:, 2] .== st)
            a = Float64.(c[k, 4]) ./ 1000; v = Float64.(c[k, 5]); w = a .<= tmax
            col_, ls, lw = get(MODEL_STYLE, src, (:gray, :solid, 1.0))
            lines!(ax, -a[w], v[w], color = col_, linestyle = ls, linewidth = lw, label = src)
            lo = min(lo, minimum(v[w])); hi = max(hi, maximum(v[w]))
        end
        xlims!(ax, -tmax, 0.5); ylims!(ax, lo - 0.05 * (hi - lo) - 1, hi + 0.05 * (hi - lo) + 1)
    end
    Legend(fig[cld(length(sites), nc) + 1, 1:nc], ax, orientation = :horizontal, framevisible = false, merge = true)
    save(joinpath(dir, "pm_ind_sites.png"), fig, px_per_unit = 1)

    s, hs = readdlm(joinpath(dir, "pm_ind_sites.csv"), ',', header = true)
    sc(n) = s[:, findfirst(==(n), vec(hs))]
    model = sc("model"); nip = Int.(sc("n_ip")); rz = Float64.(sc("rms_z")); la = Float64.(sc("lat")); lo_ = Float64.(sc("lon"))
    show = [mn for mn in ("PaleoMIST", "B1_l256") if mn in model]
    fig = Figure(size = (560 * length(show) + 120, 340), fontsize = 12)
    local hm
    for (j, mn) in enumerate(show)
        k = (model .== mn) .& (nip .> 0)
        ax = Axis(fig[1, j], title = "$mn: RMS of z per site (index points)", aspect = DataAspect())
        hidedecorations!(ax)
        land === nothing || contour!(ax, lon, lat, land; levels = [0.5], color = :gray60, linewidth = 0.5)
        hm = scatter!(ax, lo_[k], la[k], color = max.(rz[k], 0.1), colormap = :viridis, colorscale = log10, colorrange = (0.5, 100),
                      markersize = 4 .+ 1.5 .* sqrt.(nip[k]), strokecolor = :black, strokewidth = 0.3)
        limits!(ax, -180, 180, -90, 90)
    end
    Colorbar(fig[1, length(show) + 1], hm, label = "RMS z per site")
    save(joinpath(dir, "pm_ind_map.png"), fig, px_per_unit = 1)
end

function main(args)
    dir = args[1]
    ds = NCDataset(joinpath(dir, "pm_fields.nc"))
    lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:]); ages = Float64.(ds["age"][:])
    has(v) = haskey(ds, v)
    get3(v) = Float64.(ds[v][:, :, :]); get2(v) = Float64.(ds[v][:, :])
    H = get2("H_lgm")
    nr = length(ages)

    if has("rsl_B1")
        PM = get3("rsl_PM"); B1 = get3("rsl_B1")
        fig = Figure(size = (430 * nr + 100, 3 * 250), fontsize = 12)
        local h1, h2
        for (c, a) in enumerate(ages)
            h1 = gmap!(fig, (1, c), lon, lat, PM[:, :, c], 200, @sprintf("PaleoMIST, %.1f ka", a / 1000); ice = H)
            gmap!(fig, (2, c), lon, lat, B1[:, :, c], 200, @sprintf("VILMA2 B1 l256, %.1f ka", a / 1000); ice = H)
            h2 = gmap!(fig, (3, c), lon, lat, B1[:, :, c] .- PM[:, :, c], 30, "B1 - PaleoMIST"; ice = H)
        end
        Colorbar(fig[1:2, nr+1], h1, label = "RSL rel. to present [m]"); Colorbar(fig[3, nr+1], h2, label = "difference [m]")
        save(joinpath(dir, "pm_B1_vs_PM.png"), fig, px_per_unit = 1)
    end
    if has("rsl_B1") && (has("rsl_B1_170") || has("rsl_B3"))
        rows = [(k1, k2, lbl) for (k1, k2, lbl) in (("rsl_B1_170", "rsl_B1", "lmax 170 - 256 (B1)"), ("rsl_B3", "rsl_B1_170", "3-D - 1-D (B3 - B1, l170)"))
                if has(k1) && has(k2)]
        fig = Figure(size = (430 * nr + 100, 250 * length(rows)), fontsize = 12)
        for (r, (k1, k2, lbl)) in enumerate(rows)
            local h; lim = r == 1 && k1 == "rsl_B1_170" ? 5 : 30
            for (c, a) in enumerate(ages)
                h = gmap!(fig, (r, c), lon, lat, get3(k1)[:, :, c] .- get3(k2)[:, :, c], lim, @sprintf("%s, %.1f ka", lbl, a / 1000); ice = H)
            end
            Colorbar(fig[r, nr+1], h, label = "RSL difference [m]")
        end
        save(joinpath(dir, "pm_effects.png"), fig, px_per_unit = 1)
    end
    rows = [(k1, k2, lbl) for (k1, k2, lbl) in (("rsl_B1_L80", "rsl_B1_170", "lithosphere: B1_L80 - B1"),
                                                 ("rsl_B1_bm", "rsl_B1_L80", "radial profile: B1_bm - B1_L80"),
                                                 ("rsl_B3", "rsl_B1_bm", "lateral: B3 - B1_bm")) if has(k1) && has(k2)]
    if !isempty(rows)
        fig = Figure(size = (430 * nr + 100, 250 * length(rows)), fontsize = 12)
        for (r, (k1, k2, lbl)) in enumerate(rows)
            local h
            for (c, a) in enumerate(ages)
                h = gmap!(fig, (r, c), lon, lat, get3(k1)[:, :, c] .- get3(k2)[:, :, c], 30, @sprintf("%s, %.1f ka", lbl, a / 1000); ice = H)
            end
            Colorbar(fig[r, nr+1], h, label = "RSL difference [m]")
        end
        save(joinpath(dir, "pm_decomp.png"), fig, px_per_unit = 1)
    end
    pairs = [(b, lbl) for (b, lbl) in (("B1", "1-D, l256"), ("B3", "3-D, l170")) if has("AmB_rsl_" * b)]
    if !isempty(pairs)
        fig = Figure(size = (2 * 430 + 100, 250 * length(pairs)), fontsize = 12)
        for (r, (b, lbl)) in enumerate(pairs)
            h = gmap!(fig, (r, 1), lon, lat, get2("AmB_rsl_" * b), 300, "A - B, RSL, " * lbl; ice = H)
            gmap!(fig, (r, 2), lon, lat, get2("AmB_u_" * b), 300, "A - B, u (bed displacement), " * lbl; ice = H)
            Colorbar(fig[r, 3], h, label = "equilibrium - transient, 22.5 ka [m]")
        end
        save(joinpath(dir, "pm_A_minus_B.png"), fig, px_per_unit = 1)
    end
    zs = [(b, lbl) for (b, lbl) in (("B1", "1-D, l256"), ("B3", "3-D, l170")) if has("rsl0_" * b)]
    if !isempty(zs)
        fig = Figure(size = (430 * length(zs) + 100, 260), fontsize = 12)
        local h
        for (c, (b, lbl)) in enumerate(zs)
            h = gmap!(fig, (1, c), lon, lat, get2("rsl0_" * b), 300, "rsl at 0 ka vs relaxed reference, " * lbl; ice = H)
        end
        Colorbar(fig[1, length(zs)+1], h, label = "[m]")
        save(joinpath(dir, "pm_rsl0.png"), fig, px_per_unit = 1)
    end
    land = has("land") ? Float64.(ds["land"][:, :]) : nothing
    close(ds)

    isfile(joinpath(dir, "pm_ind_curves.csv")) && indicator_figures(dir, lon, lat, land)

    s, _ = readdlm(joinpath(dir, "pm_sites.csv"), ',', header = true)
    sites = unique(s[:, 2]); srcs = unique(s[:, 1])
    fig = Figure(size = (1200, 740), fontsize = 12)
    local ax
    for (n, site) in enumerate(sites)
        ax = Axis(fig[(n - 1) ÷ 3 + 1, (n - 1) % 3 + 1], title = site, xlabel = "time [ka]", ylabel = "RSL [m]")
        for src in srcs
            m = (s[:, 1] .== src) .& (s[:, 2] .== site)
            t = Float64.(s[m, 3]) ./ 1000; v = Float64.(s[m, 4])
            if src == "PaleoMIST"
                scatter!(ax, t, v, color = :black, markersize = 6, label = src)
            else
                lines!(ax, t, v, label = src)
            end
        end
    end
    Legend(fig[3, 1:3], ax, orientation = :horizontal, framevisible = false)
    save(joinpath(dir, "pm_sites.png"), fig, px_per_unit = 1)

    g, _ = readdlm(joinpath(dir, "pm_gmsl.csv"), ',', header = true)
    fig = Figure(size = (1000, 420), fontsize = 12)
    ax = Axis(fig[1, 1], xlabel = "time [ka]", ylabel = "sea level rel. to present [m]")
    cols = Makie.wong_colors()
    for (n, src) in enumerate(filter(x -> !startswith(x, "B1_l170"), unique(g[:, 1])))   # B1 lmax 170 is B1 lmax 256 to 0.3 m
        m = g[:, 1] .== src; t = Float64.(g[m, 2]) ./ 1000; v = Float64.(g[m, 3])
        c = cols[mod1(n, length(cols))]
        startswith(src, "PaleoMIST") ? scatterlines!(ax, t, v, label = src, markersize = 6, color = c, linestyle = :dash) :
                                       lines!(ax, t, v, label = src, color = c)
    end
    Legend(fig[1, 2], ax, framevisible = false)
    save(joinpath(dir, "pm_gmsl.png"), fig, px_per_unit = 1)
    println("wrote figures to ", dir)
end
main(ARGS)
