# Maps of the present-day implied mass balance of plastic_closure_fit sets (tests/plastic_closure_fit.f90, the NetCDF of one set):
# the implied SMB rho_i (div q - q_land) on the ice, smoothed conservatively as the backend's mb_smooth_km (scatter form, sums
# unchanged), and log10(model / observed) of the surface speed; optionally an observed SMB on the same grid as a first column.
#   julia --project=<julia-plot> tools/levante/plastic_closure_maps.jl <out.png> <sd_km> [ref=<smb.nc>] <label>=<set.nc>[:inv] ...
# <set.nc>:inv: the set's joint inversion (prior_file of plastic_closure_fit: smb_inv, speed_inv, q_marine_inv), unsmoothed.
# Prints per set the integrals [Gt yr-1] (net, accumulation, ablation, outflow to the ocean) and the ablation area (the fraction of
# the ice cells with smoothed SMB < 0), and the same for ref over the model's ice.
using CairoMakie, NCDatasets, Printf, Statistics

const RHO_I = 910.0
const KG_GT = 1e-12

out, sd = ARGS[1], 1000 * parse(Float64, ARGS[2])
opts = [split(a, "=", limit = 2) for a in ARGS[3:end]]
ref = [p[2] for p in opts if p[1] == "ref"]
sets = [(p[1], p[2]) for p in opts if p[1] != "ref"]

g2(ds, v) = Float64.(nomissing(ds[v][:, :], NaN))
function smooth_conservative(f, mask, dx, dy, sd)
    sd > 0 || return copy(f)
    nx, ny = size(f); ri, rj = ceil(Int, 3sd / dx), ceil(Int, 3sd / dy)
    kern(a, b) = exp(-((a * dx)^2 + (b * dy)^2) / (2sd^2))
    wsum = zeros(nx, ny)
    for j in 1:ny, i in 1:nx
        mask[i, j] || continue
        for q in max(j - rj, 1):min(j + rj, ny), p in max(i - ri, 1):min(i + ri, nx)
            mask[p, q] && (wsum[i, j] += kern(p - i, q - j))
        end
    end
    g = zeros(nx, ny)
    for j in 1:ny, i in 1:nx
        mask[i, j] || continue
        for q in max(j - rj, 1):min(j + rj, ny), p in max(i - ri, 1):min(i + ri, nx)
            mask[p, q] && (g[p, q] += f[i, j] * kern(p - i, q - j) / wsum[i, j])
        end
    end
    r = copy(f); r[mask] .= g[mask]; r
end

D = map(sets) do (lab, fs)
    inv = endswith(fs, ":inv"); f = inv ? fs[1:end-4] : fs
    NCDataset(f) do ds
        x, y = Float64.(ds["xc"][:]), Float64.(ds["yc"][:])
        dx, dy = 1000 * (x[2] - x[1]), 1000 * (y[2] - y[1])
        m = Int.(ds["mask"][:, :]) .== 1
        if inv
            smb, qm = RHO_I .* g2(ds, "smb_inv"), g2(ds, "q_marine_inv")
            sp = g2(ds, "speed_inv")
        else
            divq, ql, qm = g2(ds, "div_q"), g2(ds, "q_land"), g2(ds, "q_marine")
            smb = smooth_conservative(RHO_I .* (divq .- ql), m, dx, dy, sd)
            sp = g2(ds, "speed")
        end
        so = g2(ds, "speed_obs")
        cell = dx * dy
        tot = (net = sum(smb[m]) * cell * KG_GT, acc = sum(max.(smb[m], 0)) * cell * KG_GT,
               abl = sum(min.(smb[m], 0)) * cell * KG_GT, marine = RHO_I * sum(qm) * cell * KG_GT,
               farea = count(smb[m] .< 0) / count(m))
        rat = fill(NaN, size(sp)); k = m .& (so .> 0.1) .& (sp .> 0); rat[k] .= log10.(sp[k] ./ so[k])
        S = fill(NaN, size(smb)); S[m] .= smb[m]
        (lab = lab, x = x, y = y, m = m, smb = S, rat = rat, tot = tot, cell = cell)
    end
end
for d in D
    t = d.tot
    @printf("%-28s net %8.1f  acc %8.1f  abl %8.1f  ocean %8.1f Gt/yr  ablation area %.3f\n", d.lab, t.net, t.acc, t.abl, t.marine, t.farea)
end
R = nothing
if !isempty(ref)
    d = D[1]
    R = NCDataset(ds -> g2(ds, "smb"), ref[1])
    r = fill(NaN, size(R)); r[d.m] .= R[d.m]
    v = r[d.m]; v = v[isfinite.(v)]
    @printf("%-28s net %8.1f  acc %8.1f  abl %8.1f  (over the model's ice)  ablation area %.3f\n", "ref", sum(v) * d.cell * KG_GT,
            sum(max.(v, 0)) * d.cell * KG_GT, sum(min.(v, 0)) * d.cell * KG_GT, count(v .< 0) / length(v))
    R = r
end

nc = length(D) + (R === nothing ? 0 : 1)
fig = Figure(size = (330nc + 120, 700), fontsize = 15)
lim = 1000.0; grey = RGBf(0.92, 0.92, 0.92)
col = 0
function panel!(r, c, d, F, cm, cr, title)
    ax = Axis(fig[r, c], aspect = DataAspect(), title = title, titlesize = 14)
    hidedecorations!(ax)
    heatmap!(ax, d.x, d.y, F, colormap = cm, colorrange = cr, nan_color = grey, lowclip = :black, highclip = :black)
end
if R !== nothing
    col += 1; panel!(1, col, D[1], R, :RdBu, (-lim, lim), "reference SMB")
end
hs = hr = nothing
for d in D
    global col += 1
    global hs = panel!(1, col, d, d.smb, :RdBu, (-lim, lim), @sprintf("%s\nnet %.0f Gt/yr, abl. area %.0f%%", d.lab, d.tot.net, 100d.tot.farea))
    global hr = panel!(2, col, d, d.rat, :PuOr, (-1.5, 1.5), "log10 speed / observed")
end
Colorbar(fig[1, col + 1], hs, label = "implied SMB [kg m-2 yr-1]")
Colorbar(fig[2, col + 1], hr, label = "log10(model / observed)")
save(out, fig)
println("wrote ", out)
