# The figure of a gia_refcorr.sh case (CairoMakie):  julia --project=<julia-plot> gia_refcorr_plot.jl <case folder> <out.png>
# Left: the converged correction r on the 0.5 degree grid (GLOBAL-0.5DEG_r.nc); right: the 0 ka bed misfit per run (iter.log).
using NCDatasets, CairoMakie

function main(dir, png)
    lon, lat, r = NCDataset(ds -> (Float64.(ds["lon"][:]), Float64.(ds["lat"][:]), Float64.(ds["r"][:, :])), joinpath(dir, "GLOBAL-0.5DEG_r.nc"))
    log = filter(l -> occursin("max |d| =", l), readlines(joinpath(dir, "iter.log")))
    mx = [parse(Float64, match(r"max \|d\| = ([0-9.]+)", l)[1]) for l in log]
    rms = [parse(Float64, match(r"RMS ([0-9.]+)", l)[1]) for l in log]
    fig = Figure(size = (1250, 380), fontsize = 12)
    ax = Axis(fig[1, 1], title = "present-day reference correction r [m], $(basename(dir))", aspect = DataAspect())
    hidedecorations!(ax)
    lim = ceil(maximum(abs.(r)) / 50) * 50
    hm = heatmap!(ax, lon, lat, r; colormap = :balance, colorrange = (-lim, lim), rasterize = 2)
    limits!(ax, -180, 180, -90, 90)
    Colorbar(fig[1, 2], hm, label = "[m]")
    ax2 = Axis(fig[1, 3], xlabel = "run", ylabel = "0 ka bed misfit [m]", yscale = log10, title = "convergence")
    scatterlines!(ax2, 0:length(mx)-1, mx, label = "max"); scatterlines!(ax2, 0:length(rms)-1, rms, label = "RMS")
    axislegend(ax2)
    colsize!(fig.layout, 1, Relative(0.62))
    save(png, fig, px_per_unit = 1)
    println("wrote ", png, ": r in ", extrema(r))
end
main(ARGS...)
