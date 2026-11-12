# The Earth basis of plastic_earth_basis.jl in a figure: the area-mean profile of each field and of the mean, and the mean and
# the components at three depths:
#   julia --project=tools/levante tools/levante/plastic_earth_basis_plot.jl <basis.nc> <out.png>
using CairoMakie, NCDatasets, Printf
function earth_basis_plot(basis, out; depths_km = [100.0, 250.0, 600.0])
    NCDataset(basis) do b
        lon = Float64.(b["lon"][:]); lat = Float64.(b["lat"][:]); dep = Float64.(b["depth"][:]) ./ 1e3
        names = split(b.attrib["fields"], "; "); P = Float64.(b["profile"][:, :]); K = size(b["pc"], 4); frac = b["variance_fraction"][:]
        lo = [x > 180 ? x - 360 : x for x in lon]; o = sortperm(lo)
        fig = Figure(size = (450 + 420 * length(depths_km), 330 * (K + 1)))
        ax = Axis(fig[1:K+1, 1], xlabel = "log10 eta (area mean)", ylabel = "depth [km]", yreversed = true, title = "the fields")
        for (m, nm) in enumerate(names)
            lines!(ax, P[m, :], dep, label = String(nm), linewidth = startswith(nm, "B21") ? 1.2 : 3,
                   color = startswith(nm, "B21") ? Makie.wong_colors()[1 + (m - 1) % 7] : (nm == "L24" ? :black : :red),
                   linestyle = startswith(nm, "B21") ? :solid : :dash)
        end
        μp = vec(sum(P, dims = 1) ./ size(P, 1)); lines!(ax, μp, dep, color = :gray30, linewidth = 4, label = "mean")
        axislegend(ax, position = :lb, labelsize = 11)
        for (c, dk) in enumerate(depths_km)
            kd = argmin(abs.(dep .- dk))
            for row in 0:K
                F = row == 0 ? Float64.(b["log10_visc_mean"][:, :, kd]) : Float64.(b["pc"][:, :, kd, row])
                lim = row == 0 ? (19.5, 24.0) : (-1.0, 1.0) .* max(maximum(abs, F), 1e-3)
                a = Axis(fig[row + 1, c + 1], title = row == 0 ? @sprintf("mean, %.0f km", dep[kd]) : @sprintf("pc %d (%.0f %%), %.0f km", row, 100frac[row], dep[kd]),
                         aspect = DataAspect(), xticksvisible = false, xticklabelsvisible = false, yticksvisible = false, yticklabelsvisible = false)
                heatmap!(a, lo[o], lat, F[o, :], colormap = row == 0 ? :viridis : cgrad(:RdBu, rev = true), colorrange = lim)
            end
        end
        save(out, fig, px_per_unit = 1)
    end
    println("wrote ", out)
end

earth_basis_plot(ARGS[1], ARGS[2])
