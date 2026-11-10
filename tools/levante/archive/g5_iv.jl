# G5: G3 (iv) as a series. julia --project=tools tools/levante/g5_iv.jl <vilma2 run> [--md FILE]
# Relative difference of bsl_ive * A_ocn (VILMA2's solve) and ESL * A_sea (the constants), per step, with the sea-level equivalent of the
# gathered grounded ice (log) split into ANT and external; the time of the maximum; and the relative difference of bsl_ive to the ESL
# scaled by the ocean area of the step.
using NCDatasets, Printf, Statistics
rd(ds, v) = Float64.(replace(Array(ds[v]), missing => NaN))
function main(args)
    r = args[1]; md = nothing
    i = findfirst(==("--md"), args); i === nothing || (md = args[i+1])
    gl = Dict{Float64,NamedTuple}(); asea = NaN
    for l in eachline(joinpath(r, "out.out"))
        m = match(r"t = (\S+) yr, sea-level equivalent of the ice change \[m\]:(.*), external (?:ice, all \(floating included\) )?(\S+);", l)
        if m !== nothing
            p = split(strip(m[2])); d = Dict(p[j] => parse(Float64, p[j+1]) for j in 1:2:length(p)-1)
            gl[parse(Float64, m[1])] = (ant = get(d, "Antarctica", NaN), ext = parse(Float64, m[3]))
        end
        m = match(r"ocean area of the reference \S+ m2 \(\S+ % of the sea surface (\S+) m2", l); m === nothing || (asea = parse(Float64, m[1]))
    end
    ts = NCDataset(joinpath(r, "geo_glob_ts.nc")); t = rd(ts, "time"); bi = rd(ts, "bsl_ive"); ao = rd(ts, "A_ocn"); bs = rd(ts, "bsl")
    rows = []
    for n in eachindex(t)
        g = get(gl, t[n] - 10.0, nothing); g === nothing && continue
        esl = g.ant + g.ext
        push!(rows, (t[n], abs(bi[n] * ao[n] - esl * asea) / abs(esl * asea), abs(bi[n] - esl * asea / ao[n]) / abs(esl * asea / ao[n]), bi[n] * ao[n] - esl * asea, esl, bi[n], ao[n]))
    end
    io = IOBuffer()
    k = argmax([x[2] for x in rows])
    @printf(io, "# G5: (iv) as a series (%d steps)\n\nmax relative difference as volumes %.3e at t = %g yr (bsl_ive %.3f m, ESL %.3f m, A_ocn %.4e m2, difference %.3e m3); at the first step %.3e, at the last %.3e.\nWith the ESL scaled by the ocean area of the step (ESL * A_sea / A_ocn): max %.3e.\n\n", length(rows), rows[k][2], rows[k][1], rows[k][6], rows[k][5], rows[k][7], rows[k][4], rows[1][2], rows[end][2], maximum(x[3] for x in rows))
    println(io, "| t [yr] | rel. diff (volumes) | rel. diff (scaled ESL) | difference [m3] | ESL [m] | bsl_ive [m] |\n|---|---|---|---|---|---|")
    for j in unique(round.(Int, range(1, length(rows), length = 21)))
        x = rows[j]; @printf(io, "| %g | %.2e | %.2e | %.3e | %.3f | %.3f |\n", x[1], x[2], x[3], x[4], x[5], x[6])
    end
    s = String(take!(io)); print(s); md === nothing || write(md, s)
end
main(ARGS)
