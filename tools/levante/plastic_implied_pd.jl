# The implied mass balance of a kryos plastic run (mass_balance = "implied", docs/dev/plastic-backend.md) at the present day against
# the observations: SMB of MAR v3.14 (Greenland, 1981-2010) and RACMO2.4p1 (Antarctica, 1981-2010), surface speed of Joughin et al.
# (2018) and Mouginot et al. (2019) (FesmData SMB, Velocity).
#   julia --project=tools/levante tools/levante/plastic_implied_pd.jl <run folder> <out folder> [age=0]
# The run: ice_nor.nc, ice_ant.nc with the 2-D exports smb, cmb, ux_s, uy_s, H_ice (out_vars all) and <instance>_plastic.nc with
# q_land, q_marine (out_plastic). Greenland: the observations on GRL-32KM, the model (NH-32KM) at the nearest cell centre with ice
# (within 30 km); Antarctica: both on ANT-16KM. Compared over the model's ice cells (and there where the observation is ice).
# Writes to <out>: implied_pd.md (the integrals [Gt yr-1] and the speed statistics), implied_pd.png (maps and speed scatter).
using CairoMakie, NCDatasets, Printf, Statistics

const FES = get(ENV, "FESMDATA", "/work/ba1442/fesmdata")   # tools/levante/env.sh
const RHO_I = 910.0
const KG_GT = 1e-12

run, out = abspath(ARGS[1]), ARGS[2]; mkpath(out)
opt = Dict(split(a, "=", limit = 2) for a in ARGS[3:end]); age = parse(Float64, get(opt, "age", "0"))

g2(ds, v) = Float64.(nomissing(ds[v][:, :], NaN))
function model(inst, dom, grid)
    f = joinpath(run, "$(inst).nc"); fp = joinpath(run, "$(inst)_plastic.nc")
    M = NCDataset(f) do ds
        t = Float64.(ds["time"][:]); k = argmin(abs.(t .+ age))
        g(v) = Float64.(nomissing(ds[v][:, :, k], NaN))
        (H = g("H_ice"), smb = g("smb"), cmb = g("cmb"), ux = g("ux_s"), uy = g("uy_s"), t = t[k])
    end
    P = NCDataset(fp) do ds
        t = Float64.(ds["time"][:]); k = argmin(abs.(t .+ age))
        (ql = Float64.(nomissing(ds["q_land"][:, :, k], 0.0)), qm = Float64.(nomissing(ds["q_marine"][:, :, k], 0.0)), t = t[k])
    end
    lon, lat, area = NCDataset(ds -> (g2(ds, "lon2D"), g2(ds, "lat2D"), g2(ds, "area")), joinpath(FES, dom, grid, "$(grid)_grid.nc"))
    reg = NCDataset(ds -> Int.(nomissing(ds["region_3"][:, :], -1)), joinpath(FES, dom, grid, "$(grid)_REGIONS.nc"))
    reg2 = NCDataset(ds -> Int.(nomissing(ds["region_2"][:, :], -1)), joinpath(FES, dom, grid, "$(grid)_REGIONS.nc"))
    merge(M, (ql = P.ql, qm = P.qm, lon = lon, lat = lat, area = area, reg = reg, reg2 = reg2, speed = hypot.(M.ux, M.uy)))
end
obs(dom, grid, smbfile, velfile) = begin
    S = NCDataset(ds -> (smb = g2(ds, "smb"), fice = g2(ds, "f_ice")), joinpath(FES, dom, grid, "$(grid)_$(smbfile).nc"))
    V = NCDataset(ds -> g2(ds, "uxy_srf"), joinpath(FES, dom, grid, "$(grid)_$(velfile).nc"))
    lon, lat, area = NCDataset(ds -> (g2(ds, "lon2D"), g2(ds, "lat2D"), g2(ds, "area")), joinpath(FES, dom, grid, "$(grid)_grid.nc"))
    (smb = S.smb, fice = S.fice, speed = V, lon = lon, lat = lat, area = area)
end
xyz(lo, la) = (cosd(la) * cosd(lo), cosd(la) * sind(lo), sind(la))

# --- Greenland: model (NH-32KM, region_2 = 103) at the nearest ice cell of each GRL-32KM cell -------------------------------------------
mN = model("ice_nor", "North", "NH-32KM")
oG = obs("Greenland", "GRL-32KM", "SMB-MARv3.14-1981-2010-annual-mean", "Velocity-Joughin2018")
gi = findall((mN.H .> 1) .& (mN.reg2 .== 103))
P = [xyz(mN.lon[i], mN.lat[i]) for i in gi]
near = fill(0, size(oG.smb))
for j in eachindex(oG.smb)
    q = xyz(oG.lon[j], oG.lat[j]); d = [sum(abs2, p .- q) for p in P]; k = argmin(d)
    near[j] = 6371.0 * sqrt(d[k]) < 30.0 ? k : 0
end
G = (smb_m = reshape([near[j] > 0 ? mN.smb[gi[near[j]]] : NaN for j in eachindex(near)], size(near)),
     sp_m = reshape([near[j] > 0 ? mN.speed[gi[near[j]]] : NaN for j in eachindex(near)], size(near)),
     smb_o = oG.smb, sp_o = oG.speed, ice_o = oG.fice .> 0.5, has = near .> 0)
# integrals of the model on its own grid (Greenland cells), of MAR on its grid (its ice)
gsum(x, w) = sum(x .* w) * KG_GT
mod_g = (net = gsum(mN.smb[gi], mN.area[gi]), acc = gsum(max.(mN.smb[gi], 0), mN.area[gi]), abl = gsum(min.(mN.smb[gi], 0), mN.area[gi]),
         margin = -gsum(RHO_I .* mN.ql[gi], mN.area[gi]), disch = gsum(RHO_I .* mN.qm[gi], mN.area[gi]))
oi = findall(G.ice_o .& isfinite.(oG.smb))
obs_g = (net = gsum(oG.smb[oi], oG.area[oi]), acc = gsum(max.(oG.smb[oi], 0), oG.area[oi]), abl = gsum(min.(oG.smb[oi], 0), oG.area[oi]))

# --- Antarctica: both on ANT-16KM -----------------------------------------------------------------------------------------------------
mA = model("ice_ant", "Antarctica", "ANT-16KM")
oA = obs("Antarctica", "ANT-16KM", "SMB-RACMO2.4p1-1981-2010-annual-mean", "Velocity-Mouginot2019")
ai = mA.H .> 1
ant_regions = [("West Antarctica", 20101), ("East Antarctica", 20102), ("Antarctic Peninsula", 20103), ("all", 0)]
function ant_row(code)
    k = findall(ai .& (code == 0 ? trues(size(ai)) : mA.reg .== code))
    ko = [i for i in k if isfinite(oA.smb[i])]
    (net = gsum(mA.smb[k], mA.area[k]), acc = gsum(max.(mA.smb[k], 0), mA.area[k]), abl = gsum(min.(mA.smb[k], 0), mA.area[k]),
     margin = -gsum(RHO_I .* mA.ql[k], mA.area[k]), disch = gsum(RHO_I .* mA.qm[k], mA.area[k]),
     onet = gsum(oA.smb[ko], oA.area[ko]), oacc = gsum(max.(oA.smb[ko], 0), oA.area[ko]), oabl = gsum(min.(oA.smb[ko], 0), oA.area[ko]))
end

# --- speed statistics -------------------------------------------------------------------------------------------------------------
function speedstats(m, o)
    k = findall(isfinite.(m) .& isfinite.(o) .& (m .> 0.1) .& (o .> 0.1))
    r = log10.(m[k] ./ o[k])
    (n = length(k), med_m = median(m[k]), med_o = median(o[k]), ratio = 10^median(r), lrms = sqrt(mean(r .^ 2)))
end
sG = speedstats(G.sp_m[G.has .& G.ice_o], G.sp_o[G.has .& G.ice_o])
sA = speedstats(mA.speed[ai], oA.speed[ai])

open(joinpath(out, "implied_pd.md"), "w") do io
    println(io, "Run `$run`, model record ", @sprintf("%.0f", -mN.t), " years ago. Mass balance [Gt yr-1]: net SMB (with the margin ",
            "ablation), accumulation (smb > 0), ablation (smb < 0, the margin ablation included), the margin ablation alone, the discharge; ",
            "observed over its own ice (Greenland: MAR on GRL-32KM; Antarctica: RACMO where the model has ice).\n")
    println(io, "| ice sheet | model net | model acc | model abl | margin abl | discharge | obs net | obs acc | obs abl |\n|---|---|---|---|---|---|---|---|---|")
    @printf(io, "| Greenland | %.0f | %.0f | %.0f | %.0f | %.0f | %.0f | %.0f | %.0f |\n", mod_g.net, mod_g.acc, mod_g.abl, mod_g.margin,
            mod_g.disch, obs_g.net, obs_g.acc, obs_g.abl)
    for (n, c) in ant_regions
        r = ant_row(c)
        @printf(io, "| %s | %.0f | %.0f | %.0f | %.0f | %.0f | %.0f | %.0f | %.0f |\n", n, r.net, r.acc, r.abl, r.margin, r.disch, r.onet,
                r.oacc, r.oabl)
    end
    println(io, "\nSurface speed [m/yr] over the ice where both are > 0.1 m/yr: median model, median observed, median of the ratio, RMS of log10 ratio.\n")
    println(io, "| ice sheet | n | median model | median obs | ratio | RMS log10 |\n|---|---|---|---|---|---|")
    for (n, s) in (("Greenland (Joughin 2018)", sG), ("Antarctica (Mouginot 2019)", sA))
        @printf(io, "| %s | %d | %.1f | %.1f | %.2f | %.2f |\n", n, s.n, s.med_m, s.med_o, s.ratio, s.lrms)
    end
end

# --- figure -----------------------------------------------------------------------------------------------------------------------
set_theme!(fontsize = 18)
fig = Figure(size = (2600, 1100))
smbcm = :RdBu; lim = 1000.0                  # blue: accumulation
for (r, (title, so, sm, spo, spm, msk)) in enumerate((
        ("Greenland (GRL-32KM)", oG.smb, G.smb_m, oG.speed, G.sp_m, G.ice_o .& G.has),
        ("Antarctica (ANT-16KM)", oA.smb, mA.smb, oA.speed, mA.speed, ai)))
    f(x) = (y = copy(x); y[.!msk] .= NaN; y)
    for (c, (t, F, cm, cr, sc)) in enumerate((("observed SMB [kg m-2 yr-1]", f(so), smbcm, (-lim, lim), identity),
                                              ("implied SMB", f(sm), smbcm, (-lim, lim), identity),
                                              ("implied - observed", f(sm .- so), smbcm, (-lim, lim), identity),
                                              ("observed speed [m/yr]", f(log10.(max.(spo, 0.1))), :viridis, (-1, 3.5), identity),
                                              ("model surface speed", f(log10.(max.(spm, 0.1))), :viridis, (-1, 3.5), identity)))
        ax = Axis(fig[r, c][1, 1], title = c == 1 ? "$title: $t" : t, aspect = DataAspect()); hidedecorations!(ax)
        hm = heatmap!(ax, F, colormap = cm, colorrange = cr, nan_color = RGBf(0.92, 0.92, 0.92), lowclip = :black, highclip = :black)
        c in (3, 5) && Colorbar(fig[r, c][1, 2], hm, label = c == 5 ? "log10 speed" : "kg m-2 yr-1")
    end
    ax = Axis(fig[r, 6], xscale = log10, yscale = log10, xlabel = "observed speed [m/yr]", ylabel = "model [m/yr]",
              limits = ((0.1, 1e4), (0.1, 1e4)))
    k = findall(msk .& isfinite.(spm) .& isfinite.(spo) .& (spm .> 0.1) .& (spo .> 0.1))
    scatter!(ax, spo[k], spm[k], markersize = 3, color = (:black, 0.3)); lines!(ax, [0.1, 1e4], [0.1, 1e4], color = :red)
end
save(joinpath(out, "implied_pd.png"), fig)
println("wrote ", out)
