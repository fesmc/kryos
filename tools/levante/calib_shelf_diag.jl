#!/usr/bin/env julia
# calib-cf (docs/dev/calib-cf.md): ice-shelf budget of the final state of ANT calib runs: dTo, ice and grounded area, floating area and mean thickness,
# area-mean bmb_shlf/bmb/smb/cmb over floating ice, shelf T and thermal forcing. On Levante: jl.sh <log> calib_shelf_diag.jl <run dir> ...
using NCDatasets, Statistics, Printf
lastf(f, v) = (ds = NCDataset(f); n = length(ds["time"]); x = Float64.(nomissing(ds[v][:, :, n], NaN)); close(ds); x)
for r in ARGS
    println("=== ", r)
    ds = NCDataset(joinpath(r, "forcing_ts.nc")); println("coef_dTo first/last ", ds["coef_dTo"][1], " ", ds["coef_dTo"][end], "  alpha ", ds["coef_alpha"][end]); close(ds)
    H = lastf(joinpath(r, "ice_ant.nc"), "H_ice"); fg = lastf(joinpath(r, "ice_ant.nc"), "f_grnd")
    a = Float64.(nomissing(NCDataset(joinpath(r, "ice_ant_yelmo.nc"))["area"][:, :, 1], 0.0)); a ./= (maximum(a) < 1e6 ? 1.0 : 1e6)  # km2
    fl = (H .> 1) .& (fg .< 0.5)
    println(@sprintf("ice area (H > 1 m) %.3f M km2, grounded %.3f", sum(a .* (H .> 1)) / 1e6, sum(a .* ((H .> 1) .& (fg .>= 0.5))) / 1e6))
    println(@sprintf("floating area %.3f M km2, mean H float %.0f m", sum(a .* fl) / 1e6, mean(H[fl])))
    y = joinpath(r, "ice_ant_yelmo.nc")
    for v in ("bmb_shlf", "bmb", "smb", "cmb")
        x = lastf(y, v); ok = .!isnan.(x) .& fl
        println(@sprintf("%-8s floating: area-mean %.3f m/yr, median %.3f, total %.0f Gt/yr", v, sum(x[ok] .* a[ok]) / sum(a[ok]), median(x[ok]), sum(x[ok] .* a[ok]) * 1e6 * 0.917 / 1e9))
    end
    s = joinpath(r, "shelf_ant.nc")
    for v in ("T_shlf", "tf_shlf"); x = lastf(s, v); println(@sprintf("%-8s (shelf grid) mean %.3f median %.3f over non-NaN", v, mean(filter(!isnan, x)), median(filter(!isnan, x)))); end
end
