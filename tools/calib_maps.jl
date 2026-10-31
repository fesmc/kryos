#!/usr/bin/env julia
# Maps and series of the calib-cf report (CairoMakie): from the files of tools/levante/calib_metrics.jl. On Levante only, as a job
# (tools/levante/calib_maps.sh, the script copied to $G/calib/calib_maps.jl):
#   JULIA_DEPOT_PATH=$G/julia-depot julia +1.13 --project=$G/julia-plot calib_maps.jl <dir with <label>_fields.nc, _ts.tsv, .md (and _reg_ts.tsv)> <out dir>
# $G = /work/ba1442/robinson/kryos-gate; $G/julia-plot: the plot environment (CairoMakie, NCDatasets, DelimitedFiles), see calib_maps.sh.
# label: ant|nor - pd|lgm|lgmrecon [-tag] (the out prefix of calib_metrics.jl). One map figure per label (<label>-maps.png; the model and
# model - reference for each reference of the fields file, PaleoMIST 22.5 ka left out; ice fronts dashed on the model panel), V_sle series per domain (vsle-series.png) and the
# regional V_af series of the North LGM runs with the range of the references (nor-lgm-regions.png).
using CairoMakie, NCDatasets, DelimitedFiles, Printf
dir, out = ARGS; mkpath(out)
rdf(f, v) = Float64.(nomissing(NCDataset(f)[v][:, :], NaN))
labels = sort([replace(f, "_fields.nc" => "") for f in readdir(dir) if endswith(f, "_fields.nc") && !startswith(f, "t-")])
refname(k) = foldl(replace, ("_" => " ", "GLAC 1D" => "GLAC-1D", "ICE 6G C" => "ICE-6G_C"); init = k[3:end])

function maps(lab)
    f = joinpath(dir, "$(lab)_fields.nc")
    ds = NCDataset(f); x = Float64.(ds["xc"][:]) ./ 1e3; y = Float64.(ds["yc"][:]) ./ 1e3
    refs = sort([k for k in keys(ds) if startswith(k, "H_") && k != "H_model" && !occursin("22_5", k)]); close(ds)
    Hm = rdf(f, "H_model"); n = length(refs)
    fig = Figure(size = (430 * (n + 1), 460))
    ax = Axis(fig[1, 1]; title = "model grounded H [m]", aspect = DataAspect()); hidedecorations!(ax)
    hm = heatmap!(ax, x, y, replace(Hm, 0.0 => NaN); colormap = :ice, colorrange = (0, 4500)); Colorbar(fig[2, 1], hm; vertical = false, flipaxis = false, ticks = 0:1500:4500)
    contour!(ax, x, y, Float64.(Hm .> 0); levels = [0.5], color = :black, linewidth = 0.8)
    ds = NCDataset(f); fr = [k for k in ("ice_model", "ice_ref") if haskey(ds, k)]; close(ds)        # ice fronts (grounded and floating), dashed
    for (k, col) in zip(fr, (:black, :magenta)); contour!(ax, x, y, rdf(f, k); levels = [0.5], color = col, linewidth = 0.8, linestyle = :dash); end
    for (k, key) in enumerate(refs)
        Hr = rdf(f, key); d = Hm .- Hr; d[(Hm .== 0) .& (Hr .== 0)] .= NaN
        ax = Axis(fig[1, k + 1]; title = "model - $(refname(key)) [m]", aspect = DataAspect()); hidedecorations!(ax)
        hd = heatmap!(ax, x, y, d; colormap = :balance, colorrange = (-1500, 1500)); Colorbar(fig[2, k + 1], hd; vertical = false, flipaxis = false, ticks = -1500:750:1500)
        contour!(ax, x, y, Float64.(Hr .> 0); levels = [0.5], color = :magenta, linewidth = 0.8)
        contour!(ax, x, y, Float64.(Hm .> 0); levels = [0.5], color = :black, linewidth = 0.8)
    end
    Label(fig[0, 1:n+1], "$lab: black = model grounded margin, magenta = reference" * (isempty(fr) ? "" : "; dashed (left): ice front"), fontsize = 14)
    save(joinpath(out, "$lab-maps.png"), fig)
end
foreach(maps, labels)

fig = Figure(size = (1100, 420))
for (j, (dom, pfx)) in enumerate((("Antarctica", "ant-"), ("North", "nor-")))
    ax = Axis(fig[1, j]; title = dom, xlabel = "model time [kyr]", ylabel = "V_sle [m]")
    for lab in filter(startswith(pfx), labels)
        f = joinpath(dir, "$(lab)_ts.tsv"); isfile(f) || continue
        d = readdlm(f, '\t', Float64; skipstart = 1); lines!(ax, d[:, 1] ./ 1e3, d[:, 2]; label = lab)
    end
    axislegend(ax; position = :lt, labelsize = 10)
end
save(joinpath(out, "vsle-series.png"), fig)

function refrange(md, reg)          # min/max of the references in the region row of the metrics table
    for l in eachline(md)
        startswith(l, "| $reg |") || continue
        v = parse.(Float64, strip.(split(l, "|")[4:end-1])); return extrema(v)
    end
    (NaN, NaN)
end
regs = filter(l -> startswith(l, "nor-lgm") && isfile(joinpath(dir, "$(l)_reg_ts.tsv")), labels)
if !isempty(regs)
    fig = Figure(size = (1200, 380))
    for (j, (reg, col)) in enumerate((("NA", 3), ("EUR", 4), ("GRL", 2)))
        ax = Axis(fig[1, j]; title = "$reg (grey: references, LGM)", xlabel = "model time [kyr]", ylabel = "V_af [m SLE]")
        lo, hi = refrange(joinpath(dir, "$(regs[1]).md"), reg); hspan!(ax, lo, hi; color = (:grey, 0.3))
        for lab in regs
            d = readdlm(joinpath(dir, "$(lab)_reg_ts.tsv"), '\t', Float64; skipstart = 1); lines!(ax, d[:, 1] ./ 1e3, d[:, col]; label = lab)
        end
        j == 1 && axislegend(ax; position = :lt, labelsize = 10)
    end
    save(joinpath(out, "nor-lgm-regions.png"), fig)
end
