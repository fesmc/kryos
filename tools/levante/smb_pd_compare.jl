# The present-day anchor of the SMB tuning: chion on the reference topography (par/kryos_GRL-ANT_data_paleoemu_pd.nml)
# against the regional climate models, Greenland (GRL-16KM) with MAR v3.14 1981-2010 and Antarctica (ANT-32KM) with
# RACMO2.3p2 1981-2010 (the reference climate of the PD run is 1981-2010): ∫smb, snowfall, melt [Gt/yr], ablation area, T_srf bias and SMB rmse on the ice of both,
# and maps of chion's and the RCM's smb and their differences. chion's smb, melt and snowfall are weighted by the RCM's
# ice-covered fraction f_ice of the cell, as in smb_fit.jl (fice=true): the RCM's values on the grid are cell means, ice
# and ice-free (smb about 0), so f_ice times chion's ice column is their counterpart.
#
#   julia --project=tools/levante tools/levante/smb_pd_compare.jl <run dir> <out.png>
using NCDatasets, Statistics, CairoMakie
R, out = ARGS
F = get(ENV, "FESMDATA", "/work/ba1442/fesmdata")
cases = [("grl", "Greenland", "$F/Greenland/GRL-16KM/GRL-16KM_SMB-MARv3.14-1981-2010-annual-mean.nc", "MAR 1981-2010"),
         ("ant", "Antarctica", "$F/Antarctica/ANT-32KM/ANT-32KM_SMB-RACMO2.3p2-1981-2010-annual-mean.nc", "RACMO2.3p2 1981-2010")]
last2(ds, v) = (a = Array(ds[v]); ndims(a) == 3 ? Float64.(coalesce.(a[:, :, end], NaN)) : Float64.(coalesce.(a, NaN)))
fig = Figure(size=(1500, 900))
for (r, (s, dom, rcm, lab)) in enumerate(cases)
    m = NCDataset(joinpath(R, "surface_$s.nc")); i = NCDataset(joinpath(R, "ice_$s.nc")); o = NCDataset(rcm)
    x = Array(m["xc"]); dA = (x[2] - x[1])^2 * 1e6
    H = last2(i, "H_ice"); ice = H .> 0
    smb, tm, mlt, sf = last2(m, "smb"), last2(m, "tsrf"), last2(m, "melt"), last2(m, "sf")
    so, to, mo, sfo = last2(o, "smb"), last2(o, "T_srf"), last2(o, "melt"), last2(o, "sf")
    fi = replace(last2(o, "f_ice"), NaN => 0.0); smb, mlt, sf = smb .* fi, mlt .* fi, sf .* fi
    ok = ice .& isfinite.(so) .& isfinite.(smb)
    gt(a) = sum(a[ok]) * dA * 1e-12
    println(dom, " (", count(ok), " ice cells; ", lab, ")")
    println("  SMB   [Gt/yr]  chion ", round(gt(smb); digits=0), "   RCM ", round(gt(so); digits=0))
    println("  snow  [Gt/yr]  chion ", round(gt(sf); digits=0), "   RCM ", round(gt(sfo); digits=0))
    println("  melt  [Gt/yr]  chion ", round(gt(mlt); digits=0), "   RCM ", round(gt(mo); digits=0))
    println("  ablation area [%]  chion ", round(100 * count(smb[ok] .< 0) / count(ok); digits=1),
            "   RCM ", round(100 * count(so[ok] .< 0) / count(ok); digits=1))
    println("  T_srf bias [K]  mean ", round(mean((tm .- to)[ok]); digits=2), "   rmse ", round(sqrt(mean(((tm .- to)[ok]) .^ 2)); digits=2))
    println("  SMB rmse [kg m-2 yr-1] ", round(sqrt(mean(((smb .- so)[ok]) .^ 2)); digits=0))
    msk(a) = ifelse.(ok, a, NaN)
    for (c, (a, t, cr, cm)) in enumerate([(msk(smb), "$dom chion smb", (-1500, 1500), :RdBu),
                                          (msk(so), "$dom $lab smb", (-1500, 1500), :RdBu),
                                          (msk(smb .- so), "chion - RCM smb", (-1000, 1000), :PuOr),
                                          (msk(tm .- to), "chion - RCM T_srf [K]", (-8, 8), :balance)])
        ax = Axis(fig[r, 2c - 1]; title=t, aspect=DataAspect()); hidedecorations!(ax)
        h = heatmap!(ax, a; colorrange=cr, colormap=cm); Colorbar(fig[r, 2c], h)
    end
end
save(out, fig)
