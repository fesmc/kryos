# Maps of the offline check of the ice backend plastic (tests/plastic_slices.f90 with an <out.nc>): H - PaleoMIST H_ice at given ages.
#   julia --project=<julia-plot> tools/levante/plastic_maps.jl <out.png> <ages ka, comma separated> <label>=<slices.nc> ...
# One row per file, one column per age; diverging colours +-1000 m; PaleoMIST's margin (H > 1 m) as a contour.
using CairoMakie, NCDatasets
png = ARGS[1]; ages = parse.(Float64, split(ARGS[2], ","))
runs = [(split(a, "=")[1], split(a, "=")[2]) for a in ARGS[3:end]]
fig = Figure(size = (420 * length(ages) + 80, 400 * length(runs)))
for (r, (lbl, f)) in enumerate(runs)
    NCDataset(f) do ds
        x = Float64.(ds["xc"][:]); y = Float64.(ds["yc"][:]); t = Float64.(ds["time"][:])
        for (c, a) in enumerate(ages)
            k = argmin(abs.(t .- a))
            H = Float64.(ds["H"][:, :, k]); P = Float64.(ds["H_pm"][:, :, k])
            d = H .- P; d[(H .<= 0) .& (P .<= 1)] .= NaN
            ax = Axis(fig[r, c]; aspect = DataAspect(), title = "$lbl, $(t[k]) ka: H - PaleoMIST", xticksvisible = false,
                      yticksvisible = false, xticklabelsvisible = false, yticklabelsvisible = false)
            heatmap!(ax, x, y, d; colormap = :balance, colorrange = (-1000, 1000), nan_color = :grey90)
            contour!(ax, x, y, Float64.(P .> 1); levels = [0.5], color = :black, linewidth = 0.6)
        end
    end
end
Colorbar(fig[:, length(ages) + 1]; colormap = :balance, limits = (-1000, 1000), label = "H - H_PaleoMIST [m]")
save(png, fig); println("wrote ", png)
