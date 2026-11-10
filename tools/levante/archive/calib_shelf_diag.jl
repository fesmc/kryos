#!/usr/bin/env julia
# calib-cf (docs/dev/calib-cf.md): ice-shelf budget of the final state of ANT calib runs: dTo, ice and grounded area, floating area and mean thickness,
# area-mean bmb_shlf/bmb/smb/cmb over floating ice, shelf T and thermal forcing; the ice front against BedMachine v4 (ice area, floating area,
# model-only / BedMachine-only ice area); the stress at the ice front (round 5, `front`). On Levante: jl.sh <log> calib_shelf_diag.jl <run dir> ...
using NCDatasets, Statistics, Printf
const BM = "/work/ba1442/fesmdata/Antarctica/ANT-32KM/ANT-32KM_TOPO-BedMachine-v4.nc"
rd2(f, v) = (ds = NCDataset(f); x = Float64.(nomissing(ds[v][:, :], 0.0)); close(ds); x)
lastf(f, v) = (ds = NCDataset(f); n = length(ds["time"]); x = Float64.(nomissing(ds[v][:, :, n], NaN)); close(ds); x)
# Ice front of the last restart bundle (tau_vm is only there): front cells are ice cells (H > 1 m) with an ice-free ocean edge neighbour, floating
# (f_grnd < 0.5) or marine-grounded (f_grnd >= 0.5, z_bed < z_sl). tau_vm [kPa] and H at the floating front, whole domain and the Ross / Filchner-Ronne
# boxes (x, y in km), with the floating area and the model-only / BedMachine-only ice area in each box; the number of marine-grounded front cells.
const BOX = (("Ross", (-700, 500), (-1500, -300)), ("Filchner-Ronne", (-1600, -300), (100, 1100)))
q(x) = isempty(x) ? "-" : join((@sprintf("%.0f", quantile(x, p)) for p in (0.1, 0.25, 0.5, 0.75, 0.9)), " / ")
function front(r, a)
    rs = filter(d -> startswith(d, "restart-"), readdir(r)); isempty(rs) && return
    f = joinpath(r, rs[argmax(parse.(Int, replace.(rs, "restart-" => "")))], "Antarctica", "ice_ant_restart.nc")
    ds = NCDataset(f); g(v) = Float64.(nomissing(ds[v][:, :, 1], 0.0))
    H, fg, zb, zs, tv = g("H_ice"), g("f_grnd"), g("z_bed"), g("z_sl"), g("tau_vm") ./ 1e3; x, y = ds["xc"][:], ds["yc"][:]; close(ds)
    maximum(abs, x) > 1e5 && (x = x ./ 1e3; y = y ./ 1e3)
    nx, ny = size(H); ocn = (H .<= 1) .& (zb .< zs)
    edge = falses(nx, ny)
    for j in 2:ny-1, i in 2:nx-1; edge[i, j] = ocn[i-1, j] | ocn[i+1, j] | ocn[i, j-1] | ocn[i, j+1]; end
    ff = edge .& (H .> 1) .& (fg .< 0.5); fgm = edge .& (H .> 1) .& (fg .>= 0.5) .& (zb .< zs)
    println("front (", basename(dirname(dirname(f))), "): tau_vm [kPa] p10 / p25 / p50 / p75 / p90, H [m] p50; marine-grounded front cells ", count(fgm))
    println(@sprintf("  %-15s n %4d  tau_vm %s  H %.0f", "all", count(ff), q(tv[ff]), median(H[ff])))
    Hb = rd2(BM, "H_ice"); ib = Hb .> 1; flb = ib .& (rd2(BM, "f_grnd") .< 0.5); ii = H .> 1
    for (nm, bx, by) in BOX
        m = [bx[1] <= x[i] <= bx[2] && by[1] <= y[j] <= by[2] for i in 1:nx, j in 1:ny]; fm = ff .& m
        println(@sprintf("  %-15s n %4d  tau_vm %s  H %s  floating model / BedMachine %.3f / %.3f, model-only / BedMachine-only %.3f / %.3f M km2",
            nm, count(fm), q(tv[fm]), any(fm) ? @sprintf("%.0f", median(H[fm])) : "-", sum(a .* (ii .& (fg .< 0.5) .& m)) / 1e6, sum(a .* (flb .& m)) / 1e6,
            sum(a .* (ii .& .!ib .& m)) / 1e6, sum(a .* (ib .& .!ii .& m)) / 1e6))
    end
end
for r in ARGS
    println("=== ", r)
    ds = NCDataset(joinpath(r, "forcing_ts.nc")); println("coef_dTo first/last ", ds["coef_dTo"][1], " ", ds["coef_dTo"][end], "  alpha ", ds["coef_alpha"][end]); close(ds)
    H = lastf(joinpath(r, "ice_ant.nc"), "H_ice"); fg = lastf(joinpath(r, "ice_ant.nc"), "f_grnd")
    a = Float64.(nomissing(NCDataset(joinpath(r, "ice_ant_yelmo.nc"))["area"][:, :, 1], 0.0)); a ./= (maximum(a) < 1e6 ? 1.0 : 1e6)  # km2
    fl = (H .> 1) .& (fg .< 0.5)
    println(@sprintf("ice area (H > 1 m) %.3f M km2, grounded %.3f", sum(a .* (H .> 1)) / 1e6, sum(a .* ((H .> 1) .& (fg .>= 0.5))) / 1e6))
    println(@sprintf("floating area %.3f M km2, mean H float %.0f m", sum(a .* fl) / 1e6, mean(H[fl])))
    Hb = rd2(BM, "H_ice"); ib = Hb .> 1; ii = H .> 1
    println(@sprintf("BedMachine: ice area %.3f, floating %.3f M km2; model-only / BedMachine-only ice area %.3f / %.3f M km2", sum(a .* ib) / 1e6,
        sum(a .* (ib .& (rd2(BM, "f_grnd") .< 0.5))) / 1e6, sum(a .* (ii .& .!ib)) / 1e6, sum(a .* (ib .& .!ii)) / 1e6))
    y = joinpath(r, "ice_ant_yelmo.nc")
    for v in ("bmb_shlf", "bmb", "smb", "cmb")
        x = lastf(y, v); ok = .!isnan.(x) .& fl
        println(@sprintf("%-8s floating: area-mean %.3f m/yr, median %.3f, total %.0f Gt/yr", v, sum(x[ok] .* a[ok]) / sum(a[ok]), median(x[ok]), sum(x[ok] .* a[ok]) * 1e6 * 0.917 / 1e9))
    end
    s = joinpath(r, "shelf_ant.nc")
    for v in ("T_shlf", "tf_shlf"); x = lastf(s, v); println(@sprintf("%-8s (shelf grid) mean %.3f median %.3f over non-NaN", v, mean(filter(!isnan, x)), median(filter(!isnan, x)))); end
    front(r, a)
end
