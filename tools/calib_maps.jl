#!/usr/bin/env julia
# Maps and series of the calib-cf report (Mac, CairoMakie): from the files of tools/levante/calib_metrics.jl.
#   julia tools/calib_maps.jl <dir with <Domain>-<case>_fields.nc and _ts.tsv> <out dir>
using CairoMakie, NCDatasets, DelimitedFiles, Printf
dir, out = ARGS; mkpath(out)
rdf(f, v) = Float64.(nomissing(NCDataset(f)[v][:, :], NaN))
function maps(dom, case, refs)
    f = joinpath(dir, "$dom-$(case)_fields.nc"); isfile(f) || return
    ds = NCDataset(f); x = Float64.(ds["xc"][:]) ./ 1e3; y = Float64.(ds["yc"][:]) ./ 1e3; close(ds)
    Hm = rdf(f, "H_model"); n = length(refs)
    fig = Figure(size = (430 * (n + 1), 460))
    ax = Axis(fig[1, 1]; title = "model grounded H [m]", aspect = DataAspect()); hidedecorations!(ax)
    hm = heatmap!(ax, x, y, replace(Hm, 0.0 => NaN); colormap = :ice, colorrange = (0, 4500)); Colorbar(fig[2, 1], hm; vertical = false, flipaxis = false)
    contour!(ax, x, y, Float64.(Hm .> 0); levels = [0.5], color = :black, linewidth = 0.8)
    for (k, (key, nm)) in enumerate(refs)
        Hr = rdf(f, key); d = Hm .- Hr; d[(Hm .== 0) .& (Hr .== 0)] .= NaN
        ax = Axis(fig[1, k + 1]; title = "model - $nm [m]", aspect = DataAspect()); hidedecorations!(ax)
        hd = heatmap!(ax, x, y, d; colormap = :balance, colorrange = (-1500, 1500)); Colorbar(fig[2, k + 1], hd; vertical = false, flipaxis = false)
        contour!(ax, x, y, Float64.(Hr .> 0); levels = [0.5], color = :magenta, linewidth = 0.8)
        contour!(ax, x, y, Float64.(Hm .> 0); levels = [0.5], color = :black, linewidth = 0.8)
    end
    Label(fig[0, 1:n+1], "$dom $case: black = model grounded margin, magenta = reference", fontsize = 14)
    save(joinpath(out, "$dom-$case-maps.png"), fig)
end
maps("Antarctica", "pd", [("H_BedMachine_v4", "BedMachine v4")])
maps("North", "pd", [("H_GEBCO_2026", "GEBCO 2026")])
for d in ("Antarctica", "North")
    maps(d, "lgm", [("H_PaleoMIST_20_ka", "PaleoMIST 20 ka"), ("H_GLAC_1D_21_ka", "GLAC-1D 21 ka"), ("H_ICE_6G_C_21_ka", "ICE-6G_C 21 ka")])
end
fig = Figure(size = (1000, 400))
for (j, dom) in enumerate(("Antarctica", "North"))
    ax = Axis(fig[1, j]; title = dom, xlabel = "time [kyr]", ylabel = "V_sle [m]")
    for (case, col) in (("pd", :steelblue), ("lgm", :firebrick))
        f = joinpath(dir, "$dom-$(case)_ts.tsv"); isfile(f) || continue
        d = readdlm(f, '\t', Float64; skipstart = 1); lines!(ax, d[:, 1] ./ 1e3, d[:, 2]; color = col, label = case)
    end
    axislegend(ax; position = :rb)
end
save(joinpath(out, "vsle-series.png"), fig)
