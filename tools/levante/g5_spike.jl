# G5: the steps around a given time of the (iv) series. julia ... g5_spike.jl <run> <t> [n steps each side]
using NCDatasets, Printf
rd(ds, v) = Float64.(replace(Array(ds[v]), missing => NaN))
r = ARGS[1]; t0 = parse(Float64, ARGS[2]); n = length(ARGS) > 2 ? parse(Int, ARGS[3]) : 5
gl = Dict{Float64,NamedTuple}()
for l in eachline(joinpath(r, "out.out"))
    m = match(r"t = (\S+) yr, sea-level equivalent of the ice change \[m\]:(.*), external (?:ice, all \(floating included\) )?(\S+);", l)
    m === nothing && continue
    p = split(strip(m[2])); d = Dict(p[j] => parse(Float64, p[j+1]) for j in 1:2:length(p)-1)
    gl[parse(Float64, m[1])] = (ant = get(d, "Antarctica", NaN), ext = parse(Float64, m[3]))
end
ts = NCDataset(joinpath(r, "geo_glob_ts.nc")); t = rd(ts, "time"); bi = rd(ts, "bsl_ive"); ao = rd(ts, "A_ocn"); bs = rd(ts, "bsl"); mr = rd(ts, "mass_resid")
k0 = argmin(abs.(t .- t0))
println("   t        ESL_ANT    ESL_ext    ESL_sum    bsl_ive      bsl        A_ocn       bsl_ive*A_ocn - ESL*A_sea [m3]")
for k in k0-n:k0+n
    g = gl[t[k] - 10.0]; esl = g.ant + g.ext
    @printf("%8.0f %10.4f %10.4f %10.4f %10.4f %10.4f %.5e %12.4e\n", t[k], g.ant, g.ext, esl, bi[k], bs[k], ao[k], bi[k] * ao[k] - esl * 3.618e14)
end
