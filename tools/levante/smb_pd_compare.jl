# The present-day anchor of the SMB tuning: chion on the reference topography (par/kryos_GRL-ANT_data_paleoemu_pd.nml)
# against the regional climate models of the same period (1981-2010), as in smb_fit.jl: Greenland (GRL-16KM) with
# MAR v3.14 averaged over the ice of each cell only ($KRYOS_DATA/smb, tools/levante/mar_ice16.jl), weighted by its ice
# share f_ice; Antarctica (ANT-32KM) with RACMO2.3p2, which has values on the ice only, unweighted. ∫smb, snowfall, melt
# [Gt/yr], the ablation area, the T_srf bias (where the target has it) and the SMB rmse on the ice of both, and maps of
# chion's and the RCM's smb and their difference.
#
#   julia --project=tools/levante tools/levante/smb_pd_compare.jl <run dir> <out.png>
using NCDatasets, Statistics, CairoMakie
R, out = ARGS
F = get(ENV, "FESMDATA", "/work/ba1442/fesmdata"); K = get(ENV, "KRYOS_DATA", "/work/ba1442/kryos_data")
# instance suffix, name, target, its label, the target's variable of the ice share weighting it ("": none)
cases = [("grl", "Greenland", "$K/smb/GRL-16KM_SMB-MARv3.14-ice-1981-2010-annual-mean.nc", "MAR 1981-2010 (ice)", "f_ice"),
         ("ant", "Antarctica", "$F/Antarctica/ANT-32KM/ANT-32KM_SMB-RACMO2.3p2-1981-2010-annual-mean.nc", "RACMO2.3p2 1981-2010", "")]
last2(ds, v) = (a = Array(ds[v]); ndims(a) == 3 ? Float64.(coalesce.(a[:, :, end], NaN)) : Float64.(coalesce.(a, NaN)))
fig = Figure(size=(1500, 900))
for (r, (s, dom, rcm, lab, wv)) in enumerate(cases)
    m = NCDataset(joinpath(R, "surface_$s.nc")); i = NCDataset(joinpath(R, "ice_$s.nc")); o = NCDataset(rcm)
    x = Array(m["xc"]); dA = (x[2] - x[1])^2 * 1e6
    H = last2(i, "H_ice")
    smb, tm, mlt, sf = last2(m, "smb"), last2(m, "tsrf"), last2(m, "melt"), last2(m, "sf")
    so, sfo = last2(o, "smb"), last2(o, "sf")
    mo = replace(last2(o, "melt"), NaN => 0.0)                       # RACMO: no melt value where none
    w = isempty(wv) ? ones(size(so)) : replace(last2(o, wv), NaN => 0.0)
    ok = (H .> 0) .& isfinite.(so) .& isfinite.(smb) .& (w .> 0)
    W = w[ok]
    gt(a) = sum(W .* a[ok]) * dA * 1e-12
    wm(a) = sum(W .* a[ok]) / sum(W)
    println(dom, " (", count(ok), " ice cells; ", lab, isempty(wv) ? ")" : ", weighted by $wv)")
    println("  SMB   [Gt/yr]  chion ", round(gt(smb); digits=0), "   RCM ", round(gt(so); digits=0))
    println("  snow  [Gt/yr]  chion ", round(gt(sf); digits=0), "   RCM ", round(gt(sfo); digits=0))
    println("  melt  [Gt/yr]  chion ", round(gt(mlt); digits=0), "   RCM ", round(gt(mo); digits=0))
    println("  ablation area [%]  chion ", round(100 * wm(Float64.(smb .< 0)); digits=1),
            "   RCM ", round(100 * wm(Float64.(so .< 0)); digits=1))
    if haskey(o, "T_srf")
        to = last2(o, "T_srf"); d = tm .- to
        println("  T_srf bias [K]  mean ", round(wm(d); digits=2), "   rmse ", round(sqrt(wm(d .^ 2)); digits=2))
    end
    println("  SMB rmse [kg m-2 yr-1] ", round(sqrt(wm((smb .- so) .^ 2)); digits=0))
    msk(a) = ifelse.(ok, a, NaN)
    for (c, (a, t, cr, cm)) in enumerate([(msk(smb), "$dom chion smb", (-1500, 1500), :RdBu),
                                          (msk(so), "$dom $lab smb", (-1500, 1500), :RdBu),
                                          (msk(smb .- so), "chion - RCM smb", (-1000, 1000), :PuOr)])
        ax = Axis(fig[r, 2c - 1]; title=t, aspect=DataAspect()); hidedecorations!(ax)
        h = heatmap!(ax, a; colorrange=cr, colormap=cm); Colorbar(fig[r, 2c], h)
    end
    close(m); close(i); close(o)
end
save(out, fig)
