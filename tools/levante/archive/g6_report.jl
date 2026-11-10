# G6 report (phase 4 WP7c): drift of the science case with geo = vilma2 against K3-sci (the fastiso `sci` run over its 1 kyr).
#   julia --project=tools tools/levante/g6_report.jl <vilma2 run> <K3-sci run> [--md FILE]
# Per domain the series of Yelmo's ice (ice_<dom>_yelmo_ts.nc, 10 yr): V_ice, A_ice (the units of the file), V_sle (m of sea level), grounded V_ice_g, at
# 0, 100, 250, 500, 750, 1000 yr (and the end of the vilma2 run), as change since the start; the drift rate (mean of d V_sle / dt over the
# first 1000 yr) of both runs, and the stop-free check: the run's last time against the planned end.
using NCDatasets, Printf, Statistics
rd(ds, v) = Float64.(replace(Array(ds[v]), missing => NaN))
function main(args)
    rv, rk = args[1:2]; md = nothing
    i = findfirst(==("--md"), args); i === nothing || (md = args[i+1])
    io = IOBuffer()
    println(io, "# G6: drift of ", rv, " against K3-sci (", rk, ")\n")
    for (dom, nm) in (("nor", "North"), ("ant", "Antarctica"))
        a = NCDataset(joinpath(rv, "ice_$(dom)_yelmo_ts.nc")); b = NCDataset(joinpath(rk, "ice_$(dom)_yelmo_ts.nc"))
        ta = rd(a, "time") .- rd(a, "time")[1]; tb = rd(b, "time") .- rd(b, "time")[1]
        println(io, "## ", nm, "\n")
        println(io, "Change since the start of each run (vilma2: from -80 ka, PaleoMIST start, transient forcing; K3-sci: from the ICE-6G/PD start, constant forcing):\n")
        println(io, "| t [yr] | V_ice vilma2 | d V_sle vilma2 [m] | d A_ice vilma2 | V_ice K3 | d V_sle K3 [m] | d A_ice K3 |\n|---|---|---|---|---|---|---|")
        for t in [0, 100, 250, 500, 750, 1000, 2000]
            ka = findfirst(x -> abs(x - t) < 1e-6, ta); kb = findfirst(x -> abs(x - t) < 1e-6, tb)
            ka === nothing && continue
            va = rd(a, "V_ice"); sa = rd(a, "V_sle"); aa = rd(a, "A_ice")
            vb = kb === nothing ? NaN : rd(b, "V_ice")[kb]; sb = kb === nothing ? NaN : rd(b, "V_sle")[kb] - rd(b, "V_sle")[1]; ab = kb === nothing ? NaN : rd(b, "A_ice")[kb] - rd(b, "A_ice")[1]
            @printf(io, "| %g | %.5g | %.3f | %.4g | %.5g | %.3f | %.4g |\n", t, va[ka], sa[ka] - sa[1], aa[ka] - aa[1], vb, sb, ab)
        end
        sa = rd(a, "V_sle"); sb = rd(b, "V_sle")
        k1a = findfirst(x -> abs(x - 1000) < 1e-6, ta); k1b = findfirst(x -> abs(x - 1000) < 1e-6, tb)
        @printf(io, "\nMean rate of V_sle over the first 1000 yr: vilma2 %.3e m/yr, K3-sci %.3e m/yr.\n\n", k1a === nothing ? NaN : (sa[k1a] - sa[1]) / 1000, k1b === nothing ? NaN : (sb[k1b] - sb[1]) / 1000)
    end
    # the Earth
    ts = NCDataset(joinpath(rv, "geo_glob_ts.nc")); t = rd(ts, "time"); bsl = rd(ts, "bsl"); mr = rd(ts, "mass_resid")
    @printf(io, "## Earth\n\nsteps %d from %g to %g yr; bsl %.3f m at the start, %.3f m at the end; max mass_resid %.2e.\n", length(t), t[1], t[end], bsl[1], bsl[end], maximum(abs.(mr)))
    s = String(take!(io)); print(s); md === nothing || write(md, s)
end
main(ARGS)
