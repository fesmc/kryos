# The present-day test of the Coulomb yield law (tests/plastic_coulomb_pd.f90; docs/backend-plastic.md): the scan table
# (best c_f per H_N x lambda_min), the best case against the reference (Gowan's tau_b x the domain factors), figures.
#   julia --project=tools/levante tools/levante/plastic_coulomb_pd_report.jl <out folder> <pd_<grid>_<tag>.nc> ...
# Writes <out>/report.md and per file pd_<grid>_<tag>.png: H_obs, the thickness error of the best Coulomb case and of the
# reference over the fit region, chi2/n against H_N per lambda_min, and a profile (the row through the thickest observed
# ice of the fit region) of the bed, the observed surface and both modelled surfaces.
using CairoMakie, NCDatasets, Printf

out = ARGS[1]; mkpath(out)
io = open(joinpath(out, "report.md"), "w")
for f in ARGS[2:end]
    tag = replace(basename(f), r"\.nc$" => "")
    ds = NCDataset(f)
    g(v) = Float64.(ds[v][:, :])
    HN = Float64.(ds["HN"][:]); lam = Float64.(ds["lam_min"][:])
    lcf = g("lcf"); chi2n = g("chi2n"); rms = g("rms"); vol = g("volume"); wide = Int.(ds["widened"][:, :])
    Hb = g("H_best"); Hr = g("H_ref"); Ho = g("H_obs"); zb = g("z_bed"); fit = Int.(ds["fit"][:, :]) .== 1
    margin = Int.(ds["margin"][:, :]) .== 1; x = Float64.(ds["xc"][:]); y = Float64.(ds["yc"][:])
    vobs = ds["obs_volume"][1]; cref = ds["ref_chi2n"][1]; rref = ds["ref_rms"][1]; vref = ds["ref_volume"][1]
    sm = ds["smooth_km"][1]; nfit = ds["n_fit"][1]
    close(ds)
    ib = argmin(chi2n)

    println(io, "## ", tag, "\n\n", nfit, " fit cells, smoothing ", sm, " km, observed volume ", @sprintf("%.3f", vobs * 1e-15),
            " 10^6 km3.\n\n| H_N [m] | lambda_min | c_f | chi2/n | RMS [m] | volume / obs | widened |\n|---|---|---|---|---|---|---|")
    for il in eachindex(lam), ih in eachindex(HN)
        @printf(io, "| %g | %g | %.4f | %.3f | %.0f | %.3f | %d |\n", HN[ih], lam[il], exp(lcf[ih, il]), chi2n[ih, il],
                rms[ih, il], vol[ih, il] / vobs, wide[ih, il])
    end
    @printf(io, "\nBest: H_N %g m, lambda_min %g, c_f %.4f: chi2/n %.3f, RMS %.0f m, volume/obs %.3f. Reference (Gowan x factors): chi2/n %.3f, RMS %.0f m, volume/obs %.3f.\n\n",
            HN[ib[1]], lam[ib[2]], exp(lcf[ib]), chi2n[ib], rms[ib], vol[ib] / vobs, cref, rref, vref / vobs)

    # the fit region's bounding box (a margin of 4 cells)
    I = findall(fit); i0 = max(1, minimum(c[1] for c in I) - 4); i1 = min(length(x), maximum(c[1] for c in I) + 4)
    j0 = max(1, minimum(c[2] for c in I) - 4); j1 = min(length(y), maximum(c[2] for c in I) + 4)
    xs = x[i0:i1]; ys = y[j0:j1]; box(a) = a[i0:i1, j0:j1]
    err(H) = box(ifelse.(fit, H .- Ho, NaN))
    asp = (xs[end] - xs[1]) / (ys[end] - ys[1])
    fig = Figure(size = (1500, 900))
    ax1 = Axis(fig[1, 1], title = "H_obs [m]", aspect = DataAspect())
    hm1 = heatmap!(ax1, xs, ys, box(ifelse.(margin, Ho, NaN)), colormap = :ice, colorrange = (0, 4000))
    Colorbar(fig[2, 1], hm1, vertical = false, flipaxis = false)
    ttl = @sprintf("Coulomb best (H_N %g m, λ_min %g, c_f %.4f): H - H_obs [m]\nχ²/n %.2f, RMS %.0f m", HN[ib[1]], lam[ib[2]],
                   exp(lcf[ib]), chi2n[ib], rms[ib])
    ax2 = Axis(fig[1, 2], title = ttl, aspect = DataAspect())
    hm2 = heatmap!(ax2, xs, ys, err(Hb), colormap = :balance, colorrange = (-1000, 1000))
    ax3 = Axis(fig[1, 3], title = @sprintf("Reference (Gowan τ_b × factors v3): H - H_obs [m]\nχ²/n %.2f, RMS %.0f m", cref, rref),
               aspect = DataAspect())
    heatmap!(ax3, xs, ys, err(Hr), colormap = :balance, colorrange = (-1000, 1000))
    Colorbar(fig[2, 2:3], hm2, vertical = false, flipaxis = false)
    # the profile: the row through the thickest observed ice of the fit region
    jp = argmax(ifelse.(fit, Ho, -Inf))[2]
    for ax in (ax1, ax2, ax3); hlines!(ax, [y[jp]], color = :black, linewidth = 1); end
    ip = [i for i in i0:i1 if margin[i, jp] || zb[i, jp] > 0 || Ho[i, jp] > 0]
    ax4 = Axis(fig[3, 1:2], xlabel = "x [km]", ylabel = "elevation [m]", title = @sprintf("profile at y = %.0f km", y[jp]))
    zs(H) = [margin[i, jp] ? zb[i, jp] + H[i, jp] : NaN for i in ip]
    lines!(ax4, x[ip], zb[ip, jp], color = :saddlebrown, label = "bed")
    lines!(ax4, x[ip], [Ho[i, jp] > 0 ? zb[i, jp] + Ho[i, jp] : NaN for i in ip], color = :black, linewidth = 2, label = "observed")
    lines!(ax4, x[ip], zs(Hb), color = :dodgerblue, linewidth = 2, label = "Coulomb best")
    lines!(ax4, x[ip], zs(Hr), color = :orangered, linewidth = 2, linestyle = :dash, label = "reference")
    axislegend(ax4, position = :lb)
    ax5 = Axis(fig[3, 3], xlabel = "H_N [m]", ylabel = "χ²/n (best c_f)", xscale = log10, title = "χ²/n against H_N and λ_min")
    for il in eachindex(lam)
        scatterlines!(ax5, HN, chi2n[:, il], label = "λ_min $(lam[il])")
    end
    hlines!(ax5, [cref], color = :orangered, linestyle = :dash, label = "reference")
    axislegend(ax5, position = :rt)
    rowsize!(fig.layout, 1, Relative(0.55))
    save(joinpath(out, tag * ".png"), fig)
end
close(io)
println(read(joinpath(out, "report.md"), String))
