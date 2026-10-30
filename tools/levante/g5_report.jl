# G5 report (phase 4 WP7c): ANT deglaciation with geo = vilma2 and the same run with FastIsostasy.
#   julia --project=tools tools/levante/g5_report.jl <vilma2 run> <fastiso run> [--md FILE] [--dt 10] [--nsel 11]
# (a) VILMA2 run: sign and magnitude of dz_bed on the ice grid (geo_ant.nc, ice_ant.nc) against the ice change since the start
#     (cells that lost > 100 m: uplift expected), (b) bsl against the sea-level equivalent of the ANT ice change (log of the gathers: ANT and
#     external PaleoMIST; Yelmo V_sle), (c) the documented differences against the FastIsostasy run.
using NCDatasets, Printf, Statistics

rd(ds, v) = Float64.(replace(Array(ds[v]), missing => NaN))
function gather_log(f)
    gl = Dict{Float64,NamedTuple}()
    for l in eachline(f)
        m = match(r"t = (\S+) yr, sea-level equivalent of the ice change \[m\]:(.*), external (\S+); clipped external ice (\S+) Gt; merge closure (\S+)", l)
        m === nothing && continue
        parts = split(strip(m[2]))
        doms = Dict{String,Float64}(parts[i] => parse(Float64, parts[i+1]) for i in 1:2:length(parts)-1)
        gl[parse(Float64, m[1])] = (ant = get(doms, "Antarctica", NaN), ext = parse(Float64, m[3]))
    end
    gl
end
function main(args)
    rv, rf = args[1:2]
    md = nothing; dt = 10.0; nsel = 11
    i = findfirst(==("--md"), args); i === nothing || (md = args[i+1])
    i = findfirst(==("--dt"), args); i === nothing || (dt = parse(Float64, args[i+1]))
    io = IOBuffer()
    println(io, "# G5: ANT deglaciation (", rv, " | ", rf, ")\n")
    # ---- (a) dz_bed against the ice change on the ice grid
    gv = NCDataset(joinpath(rv, "geo_ant.nc")); iv = NCDataset(joinpath(rv, "ice_ant.nc"))
    tg = rd(gv, "time"); ti = rd(iv, "time")
    println(io, "## (a) dz_bed (VILMA2, port on ANT-32KM) against the change of the ice thickness since the first record\n")
    println(io, "Cells that lost more than 100 m since the start: mean dz_bed change (uplift > 0), fraction with uplift; cells that gained more than 100 m: mean change, fraction with subsidence; the ice change of the cell is H_ice(t) - H_ice(t0) of ice_ant.nc.\n")
    println(io, "| t [yr] | lost cells | mean dH [m] | mean d(dz_bed) [m] | uplift fraction | gained cells | mean dH [m] | mean d(dz_bed) [m] | subsidence fraction | max d(dz_bed) [m] | min d(dz_bed) [m] | max abs(dz_ss) [m] |\n|---|---|---|---|---|---|---|---|---|---|---|---|")
    H0 = rd(iv, "H_ice")[:, :, 1]; Z0 = rd(gv, "dz_bed")[:, :, 1]
    for k in 2:length(tg)
        kk = findfirst(x -> abs(x - tg[k]) < 1e-6, ti); kk === nothing && continue
        dH = rd(iv, "H_ice")[:, :, kk] .- H0; dZ = rd(gv, "dz_bed")[:, :, k] .- Z0; dS = rd(gv, "dz_ss")[:, :, k]
        lo = dH .< -100; ga = dH .> 100
        @printf(io, "| %g | %d | %.0f | %.2f | %.3f | %d | %.0f | %.2f | %.3f | %.1f | %.1f | %.2f |\n", tg[k], count(lo), isempty(dH[lo]) ? NaN : mean(dH[lo]),
                isempty(dZ[lo]) ? NaN : mean(dZ[lo]), isempty(dZ[lo]) ? NaN : mean(dZ[lo] .> 0), count(ga), isempty(dH[ga]) ? NaN : mean(dH[ga]),
                isempty(dZ[ga]) ? NaN : mean(dZ[ga]), isempty(dZ[ga]) ? NaN : mean(dZ[ga] .< 0), maximum(dZ), minimum(dZ), maximum(abs.(dS)))
    end
    println(io)
    # ---- (b) bsl against the ESL of the ANT ice change
    ts = NCDataset(joinpath(rv, "geo_glob_ts.nc")); t = rd(ts, "time"); bsl = rd(ts, "bsl"); bi = rd(ts, "bsl_ive")
    gl = gather_log(joinpath(rv, "out.out"))
    yt = NCDataset(joinpath(rv, "ice_ant_yelmo_ts.nc")); tt = rd(yt, "time"); vsle = rd(yt, "V_sle")
    println(io, "## (b) bsl against the sea-level equivalent of the ANT ice change (VILMA2 run)\n")
    println(io, "ESL of the gather at t - dt (log: the change against the present-day reference, ice more than the reference < 0): ANT, external (PaleoMIST outside the footprint), sum; VILMA2 `bsl_ive` (solved at t, ice of t - dt) and `bsl` (ocean mean of the relative sea level change); Yelmo `V_sle` change since the start (ANT, > 0 for more ice).\n")
    println(io, "| t [yr] | ESL ANT [m] | ESL external [m] | sum [m] | bsl_ive [m] | bsl [m] | bsl - sum [m] | d V_sle ANT (Yelmo) [m] | -(ESL ANT(t) - ESL ANT(t0)) [m] |\n|---|---|---|---|---|---|---|---|---|")
    ks = unique(round.(Int, range(1, length(t), length = nsel)))
    k0 = findfirst(x -> haskey(gl, x - dt), t)
    a0 = k0 === nothing ? NaN : gl[t[k0]-dt].ant
    v0 = vsle[1]
    for n in ks
        g = get(gl, t[n] - dt, nothing); g === nothing && continue
        kv = findfirst(x -> abs(x - t[n]) < 1e-6, tt)
        @printf(io, "| %g | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |\n", t[n], g.ant, g.ext, g.ant + g.ext, bi[n], bsl[n], bsl[n] - (g.ant + g.ext),
                kv === nothing ? NaN : vsle[kv] - v0, -(g.ant - a0))
    end
    ant_esl = [(tn, gl[tn].ant) for tn in sort(collect(keys(gl)))]
    @printf(io, "\nANT ESL of the gathers: first %.3f m at %g yr, last %.3f m at %g yr (change %.3f m); external: first %.3f, last %.3f m; bsl at the start %.3f m, at the end %.3f m.\n\n",
            ant_esl[1][2], ant_esl[1][1], ant_esl[end][2], ant_esl[end][1], ant_esl[end][2] - ant_esl[1][2], gl[ant_esl[1][1]].ext, gl[ant_esl[end][1]].ext, bsl[1], bsl[end])
    # ---- (c) FastIsostasy
    if isdir(rf) && isfile(joinpath(rf, "geo_ant.nc"))
        gf = NCDataset(joinpath(rf, "geo_ant.nc")); ifi = NCDataset(joinpath(rf, "ice_ant.nc"))
        tf = rd(gf, "time"); tif = rd(ifi, "time")
        yf = NCDataset(joinpath(rf, "ice_ant_yelmo_ts.nc")); tyf = rd(yf, "time")
        bf = isfile(joinpath(rf, "geo_glob_bsl_ts.nc")) ? NCDataset(joinpath(rf, "geo_glob_bsl_ts.nc")) : NCDataset(joinpath(rf, "geo_glob_ts.nc"))
        tb = rd(bf, "time"); bslf = haskey(bf, "bsl") ? rd(bf, "bsl") : fill(NaN, length(tb))
        println(io, "## (c) VILMA2 against FastIsostasy (same run: `geo = fastisostasy`, `[barysealevel] method` as the base)\n")
        println(io, "Per record: ANT grounded ice volume V_ice [1e6 km3] and V_sle [m] of both runs; dz_bed over the cells of the start ice mask (H_ice(t0) > 0): mean and min/max; bsl [m]; uplift fraction in cells that lost > 100 m.\n")
        println(io, "| t [yr] | V_sle VILMA2 | V_sle FI | diff [m] | bsl VILMA2 | bsl FI | dz_bed mean V2 | dz_bed mean FI | dz_bed max V2 | dz_bed max FI | dz_bed min V2 | dz_bed min FI | uplift frac V2 | uplift frac FI |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        Hf0 = rd(ifi, "H_ice")[:, :, 1]; Zf0 = rd(gf, "dz_bed")[:, :, 1]; mask = H0 .> 0
        for k in 1:length(tg)
            kf = findfirst(x -> abs(x - tg[k]) < 1e-6, tf); kf === nothing && continue
            kvv = findfirst(x -> abs(x - tg[k]) < 1e-6, tt); kvf = findfirst(x -> abs(x - tg[k]) < 1e-6, tyf)
            kb = findfirst(x -> abs(x - tg[k]) < 1e-6, t); kbf = findfirst(x -> abs(x - tg[k]) < 1e-6, tb)
            kk = findfirst(x -> abs(x - tg[k]) < 1e-6, ti); kkf = findfirst(x -> abs(x - tg[k]) < 1e-6, tif)
            (kvv === nothing || kvf === nothing || kk === nothing || kkf === nothing) && continue
            zv = rd(gv, "dz_bed")[:, :, k] .- Z0; zf = rd(gf, "dz_bed")[:, :, kf] .- Zf0
            dHv = rd(iv, "H_ice")[:, :, kk] .- H0; dHf = rd(ifi, "H_ice")[:, :, kkf] .- Hf0
            lov = dHv .< -100; lof = dHf .< -100
            @printf(io, "| %g | %.3f | %.3f | %.3f | %.3f | %.3f | %.2f | %.2f | %.1f | %.1f | %.1f | %.1f | %.3f | %.3f |\n", tg[k], vsle[kvv], rd(yf, "V_sle")[kvf], vsle[kvv] - rd(yf, "V_sle")[kvf],
                    kb === nothing ? NaN : bsl[kb], kbf === nothing ? NaN : bslf[kbf], mean(zv[mask]), mean(zf[mask]), maximum(zv[mask]), maximum(zf[mask]), minimum(zv[mask]), minimum(zf[mask]),
                    isempty(zv[lov]) ? NaN : mean(zv[lov] .> 0), isempty(zf[lof]) ? NaN : mean(zf[lof] .> 0))
        end
        println(io)
    end
    s = String(take!(io)); print(s); md === nothing || write(md, s)
end
main(ARGS)
