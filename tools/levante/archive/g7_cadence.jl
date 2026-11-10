# G7: the error of a coarser geo cadence. julia --project=tools tools/levante/g7_cadence.jl <g7 root> [--md FILE]
# dz_bed of geo_glob.nc (Gauss grid, rounded, 2D every 100 yr) and bsl of geo_glob_ts.nc of the runs dt20, dt50, dt100 against dt10 (coupled
# runs: the ice responds to the bed, so the differences are those of the whole coupled system), at the common records.
using NCDatasets, Printf, Statistics
rd(ds, v) = Float64.(replace(Array(ds[v]), missing => NaN))
function main(args)
    root = args[1]; md = nothing
    i = findfirst(==("--md"), args); i === nothing || (md = args[i+1])
    io = IOBuffer()
    println(io, "# G7: coarser geo cadence against dt = 10 yr (coupled runs, 400 yr from -80 ka)\n")
    println(io, "| dt [yr] | t [yr] | max abs d(dz_bed) [m] | rms d(dz_bed) [m] | max abs dz_bed of dt10 [m] | d(bsl) [m] | bsl dt10 [m] |\n|---|---|---|---|---|---|---|")
    r = NCDataset(joinpath(root, "dt10", "geo_glob.nc")); tr = rd(r, "time"); zr = rd(r, "dz_bed")
    rs = NCDataset(joinpath(root, "dt10", "geo_glob_ts.nc")); trs = rd(rs, "time"); br = rd(rs, "bsl")
    for dt in (20, 50, 100)
        a = NCDataset(joinpath(root, "dt$dt", "geo_glob.nc")); ta = rd(a, "time"); za = rd(a, "dz_bed")
        as = NCDataset(joinpath(root, "dt$dt", "geo_glob_ts.nc")); tas = rd(as, "time"); ba = rd(as, "bsl")
        for (k, t) in enumerate(ta)
            kr = findfirst(x -> abs(x - t) < 1e-6, tr); kr === nothing && continue
            d = za[:, :, k] .- zr[:, :, kr]
            k1 = findfirst(x -> abs(x - t) < 1e-6, tas); k2 = findfirst(x -> abs(x - t) < 1e-6, trs)
            @printf(io, "| %d | %g | %.4f | %.5f | %.2f | %.4f | %.3f |\n", dt, t, maximum(abs.(d)), sqrt(mean(d .^ 2)), maximum(abs.(zr[:, :, kr])), ba[k1] - br[k2], br[k2])
        end
    end
    s = String(take!(io)); print(s); md === nothing || write(md, s)
end
main(ARGS)
