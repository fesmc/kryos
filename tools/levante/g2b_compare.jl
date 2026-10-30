# G2b: kryos (tests/g2b_replay.f90: mock ice, host-resolved flotation, the globe) against the standalone vilma.x on the same history.
#   julia --project=tools tools/levante/g2b_compare.jl <run base: .../kx-<sha>/S9> <tag lmax a> <tag lmax b> ... [--md FILE] [--sign-thr 100] [--sign-dt 500]
# Each tag is a folder <base>/<tag>/{kryos,vilma}; lmax is read from the Gauss grid (nphi = 4 lmax).
# Per lmax: differences kryos - standalone of dz_bed, rsl and dz_ss in the host regions (kryos host_w > 0.5, NH and ANT by hemisphere), in the
# rim (0 < w <= 0.5) and outside (w = 0), and of bsl; then the sign test on the ice grids of the two domains.
using NCDatasets, Printf, Statistics

# ---- helpers -------------------------------------------------------------------------------------------------------------

function wmean(x, w, m)                      # area-weighted mean of x over the mask m
    s = sum(w[m]); s == 0 ? NaN : sum(w[m] .* x[m]) / s
end
wrms(d, w, m) = sqrt(wmean(d .^ 2, w, m))
amax(d, m) = any(m) ? maximum(abs.(d[m])) : NaN

# bilinear interpolation of a Gauss-grid field (lon ascending from lon[1], periodic; lat ascending) onto points (lonp, latp) [deg]
function bilin(field, lon, lat, lonp, latp)
    nx, ny = length(lon), length(lat); dlon = lon[2] - lon[1]
    out = similar(lonp, Float64)
    for k in eachindex(lonp)
        x = mod(lonp[k] - lon[1], 360.0) / dlon
        ii = floor(Int, x); fx = x - ii; i0 = mod(ii, nx) + 1; i1 = mod(ii + 1, nx) + 1
        j = searchsortedlast(lat, latp[k])
        if j < 1
            v0 = field[i0, 1]; v1 = field[i1, 1]; out[k] = (1 - fx) * v0 + fx * v1
        elseif j >= ny
            v0 = field[i0, ny]; v1 = field[i1, ny]; out[k] = (1 - fx) * v0 + fx * v1
        else
            fy = (latp[k] - lat[j]) / (lat[j+1] - lat[j])
            a = (1 - fx) * field[i0, j] + fx * field[i1, j]
            b = (1 - fx) * field[i0, j+1] + fx * field[i1, j+1]
            out[k] = (1 - fy) * a + fy * b
        end
    end
    out
end

function readrec(ds, name, k)
    v = ds[name]; nd = ndims(v)
    nd == 3 ? Float64.(nomissing(v[:, :, k], 0.0)) : error("not a 3-D variable: $name")
end

# ---- one lmax ------------------------------------------------------------------------------------------------------------

function analyse(base, tag, io; sthr, sdt)
    kd = joinpath(base, tag, "kryos"); vd = joinpath(base, tag, "vilma")
    k = NCDataset(joinpath(kd, "kryos_out.nc")); v = NCDataset(joinpath(vd, "vilma_out.nc"))
    lon = Float64.(k["lon"][:]); lat = Float64.(k["lat"][:]); nphi = length(lon); nlat = length(lat); lmax = nphi ÷ 4
    vlat = Float64.(v["lat"][:]); vlon = Float64.(v["lon"][:])
    flip = (vlat[1] > vlat[end])
    maximum(abs.(sort(vlat) .- lat)) < 1e-6 || error("the Gauss latitudes of the two runs differ")
    rd(name, r) = (a = Float64.(nomissing(v[name][:, :, r], 0.0)); flip ? a[:, end:-1:1] : a)
    tr = Float64.(k["rec"][:]); tstep = Float64.(k["step"][:]); nrec = length(tr)
    nv = size(v["dz_bed"], 3)
    w = Float64.(k["host_w"][:, :]); area = Float64.(k["area"][:, :]); LAT = repeat(reshape(lat, 1, :), nphi, 1)
    host = w .> 0.5; rim = (w .> 0.0) .& .!host; out = w .== 0.0
    regs = [("host NH", host .& (LAT .> 0)), ("host ANT", host .& (LAT .< 0)), ("rim", rim), ("outside", out)]
    # standalone record index of a kryos step n (record 0 of both = the seed)
    stepof(r) = round(Int, tr[r] - tr[1]) ÷ 10
    println(io, "\n## lmax ", lmax, " (", tag, ")\n")
    println(io, "records: kryos ", nrec, " (every ", stepof(2) - stepof(1), " steps), standalone ", nv, "; host region (w > 0.5): ",
            sum(host .* area) / sum(area) * 100 |> x -> @sprintf("%.2f", x), " % of the sphere, NH ", @sprintf("%.2f", sum((host .& (LAT .> 0)) .* area) / sum(area) * 100),
            " %, ANT ", @sprintf("%.2f", sum((host .& (LAT .< 0)) .* area) / sum(area) * 100), " %; rim ", @sprintf("%.2f", sum(rim .* area) / sum(area) * 100), " %\n")
    # series
    bk = Float64.(k["bsl"][:]); bv = Float64.(v["bsl"][:]); n = min(length(bk), length(bv))
    db = bk[1:n] .- bv[1:n]; rng = maximum(bv[1:n]) - minimum(bv[1:n])
    iw = argmax(abs.(db))
    println(io, "bsl [m]: standalone ", @sprintf("%.3f", bv[1]), " at the seed to ", @sprintf("%.3f", bv[n]), " at the end (range ", @sprintf("%.3f", rng), ");")
    println(io, @sprintf("  kryos - standalone: seed %.4f, end %.4f, max |d| %.4f at t = %.0f yr, rms %.4f; max |d| / range %.3e, at the seed / range %.3e, end / range %.3e",
        db[1], db[n], abs(db[iw]), tstep[iw], sqrt(mean(db .^ 2)), abs(db[iw]) / rng, abs(db[1]) / rng, abs(db[n]) / rng), "\n")
    dbs = db .- db[1]
    println(io, @sprintf("  of the change since the seed: max |d(bsl - bsl_0)| %.4f m = %.3e of the range\n", maximum(abs.(dbs)), maximum(abs.(dbs)) / rng))
    # fields
    stats = Dict{String,Any}()
    println(io, "Differences kryos - standalone, area-weighted rms [m] (max abs [m]) at selected times, and the maximum over all records:\n")
    sel = unique(sort([1, argmin(abs.(tr .- (-15000))), argmin(abs.(tr .- (-10000))), argmin(abs.(tr .- (-5000))), nrec]))
    for fld in ("dz_bed", "rsl", "dz_ss")
        println(io, "**", fld, "**\n")
        println(io, "| region | ", join([@sprintf("t = %.0f", tr[r]) for r in sel], " | "), " | max over records (rms) | max over records (abs) | max abs standalone | max rms / max abs standalone |")
        println(io, "|---|", "---|"^(length(sel) + 4))
        D0 = Dict{String,Matrix{Float64}}()
        series = Dict(nm => (rms = zeros(nrec), mx = zeros(nrec), ref = zeros(nrec)) for (nm, _) in regs)
        for r in 1:nrec
            a = readrec(k, fld, r); b = rd(fld, stepof(r) + 1)
            d = a .- b
            for (nm, m) in regs
                series[nm].rms[r] = wrms(d, area, m); series[nm].mx[r] = amax(d, m); series[nm].ref[r] = amax(b, m)
            end
        end
        for (nm, m) in regs
            s = series[nm]
            println(io, "| ", nm, " | ", join([@sprintf("%.3f (%.3f)", s.rms[r], s.mx[r]) for r in sel], " | "),
                    @sprintf(" | %.3f | %.3f | %.2f | %.2e |", maximum(s.rms), maximum(s.mx), maximum(s.ref), maximum(s.rms) / maximum(s.ref)))
            stats[fld*"|"*nm] = (maxrms = maximum(s.rms), maxabs = maximum(s.mx), ref = maximum(s.ref), endrms = s.rms[nrec], endabs = s.mx[nrec], seedrms = s.rms[1])
        end
        println(io)
    end
    # the changes since the seed (response), dz_bed
    println(io, "**dz_bed, change since the seed** (kryos - standalone of dz_bed(t) - dz_bed(seed)), maximum over records of the rms [m] (abs [m]); the largest change of the standalone [m]:\n")
    a0 = readrec(k, "dz_bed", 1); b0 = rd("dz_bed", 1)
    println(io, "| region | max rms | max abs | largest change (standalone) | ratio rms |")
    println(io, "|---|---|---|---|---|")
    rr = Dict(nm => (zeros(nrec), zeros(nrec), zeros(nrec)) for (nm, _) in regs)
    for r in 1:nrec
        d = (readrec(k, "dz_bed", r) .- a0) .- (rd("dz_bed", stepof(r) + 1) .- b0); c = rd("dz_bed", stepof(r) + 1) .- b0
        for (nm, m) in regs
            rr[nm][1][r] = wrms(d, area, m); rr[nm][2][r] = amax(d, m); rr[nm][3][r] = amax(c, m)
        end
    end
    for (nm, _) in regs
        println(io, @sprintf("| %s | %.3f | %.3f | %.2f | %.2e |", nm, maximum(rr[nm][1]), maximum(rr[nm][2]), maximum(rr[nm][3]), maximum(rr[nm][1]) / maximum(rr[nm][3])))
        stats["dchg|"*nm] = (maxrms = maximum(rr[nm][1]), maxabs = maximum(rr[nm][2]), ref = maximum(rr[nm][3]))
    end
    println(io)
    # sign test on the ice grids
    println(io, "**Sign test** (ice grid, footprint cells; interval ", sdt, " yr, loss > ", sthr, " m): fraction of the cells with loss whose bed rose, and mean change of dz_bed over them [m]; kryos port (coupler bilinear), kryos Gauss field and standalone, both bilinear in Julia\n")
    println(io, "| domain | intervals | cells with loss (median) | uplift fraction, kryos port: min / median | kryos Gauss: min / median | standalone: min / median | mean d(dz_bed) over lost cells: kryos / standalone (median) | gain: subsidence fraction kryos / standalone (median) |")
    println(io, "|---|---|---|---|---|---|---|---|")
    ni = round(Int, sdt / (tr[2] - tr[1]))
    for dn in ("North", "Antarctica")
        dd = NCDataset(joinpath(kd, "g2b_$dn.nc"))
        LON = Float64.(dd["lon"][:, :]); LATd = Float64.(dd["lat"][:, :]); fp = Float64.(dd["footprint"][:, :]) .> 0.5
        idx = findall(fp); lonp = LON[idx]; latp = LATd[idx]
        fk = Float64[]; fg = Float64[]; fs = Float64[]; nl = Int[]; mk = Float64[]; ms = Float64[]; gk = Float64[]; gs = Float64[]
        for r in 1:ni:nrec-ni
            r2 = r + ni
            H1 = Float64.(dd["H_ice"][:, :, r])[idx]; H2 = Float64.(dd["H_ice"][:, :, r2])[idx]; dH = H2 .- H1
            pk = (Float64.(dd["dz_bed"][:, :, r2]) .- Float64.(dd["dz_bed"][:, :, r]))[idx]
            gk_ = bilin(readrec(k, "dz_bed", r2), lon, lat, lonp, latp) .- bilin(readrec(k, "dz_bed", r), lon, lat, lonp, latp)
            gs_ = bilin(rd("dz_bed", stepof(r2) + 1), lon, lat, lonp, latp) .- bilin(rd("dz_bed", stepof(r) + 1), lon, lat, lonp, latp)
            lost = dH .< -sthr; gain = dH .> sthr
            sum(lost) < 20 && continue
            push!(nl, sum(lost)); push!(fk, mean(pk[lost] .> 0)); push!(fg, mean(gk_[lost] .> 0)); push!(fs, mean(gs_[lost] .> 0))
            push!(mk, mean(gk_[lost])); push!(ms, mean(gs_[lost]))
            if sum(gain) >= 20; push!(gk, mean(gk_[gain] .< 0)); push!(gs, mean(gs_[gain] .< 0)); end
        end
        close(dd)
        f(x) = isempty(x) ? "-" : @sprintf("%.3f / %.3f", minimum(x), median(x))
        println(io, @sprintf("| %s | %d | %s | %s | %s | %s | %.3f / %.3f | %s / %s |", dn, length(nl), isempty(nl) ? "-" : string(round(Int, median(nl))), f(fk), f(fg), f(fs),
                isempty(mk) ? NaN : median(mk), isempty(ms) ? NaN : median(ms), isempty(gk) ? "-" : @sprintf("%.3f", median(gk)), isempty(gs) ? "-" : @sprintf("%.3f", median(gs))))
        stats["sign|"*dn] = (kryos = isempty(fk) ? NaN : minimum(fk), gauss = isempty(fg) ? NaN : minimum(fg), standalone = isempty(fs) ? NaN : minimum(fs), n = length(nl))
    end
    println(io)
    close(k); close(v)
    return (lmax = lmax, bsl_max = abs(db[iw]), bsl_rel = abs(db[iw]) / rng, bsl_seed = abs(db[1]), bsl_end = abs(db[n]), rng = rng, stats = stats)
end

function main(args)
    md = nothing; sthr = 100.0; sdt = 500.0
    i = findfirst(==("--md"), args); i === nothing || (md = args[i+1]; deleteat!(args, i:i+1))
    i = findfirst(==("--sign-thr"), args); i === nothing || (sthr = parse(Float64, args[i+1]); deleteat!(args, i:i+1))
    i = findfirst(==("--sign-dt"), args); i === nothing || (sdt = parse(Float64, args[i+1]); deleteat!(args, i:i+1))
    base = args[1]; tags = args[2:end]
    io = IOBuffer()
    println(io, "# G2b: kryos (mock ice, host-resolved flotation) vs the standalone vilma.x\n\nbase: ", base)
    res = [analyse(base, t, io; sthr = sthr, sdt = sdt) for t in tags]
    sort!(res, by = r -> r.lmax)
    println(io, "\n## Convergence with lmax\n")
    println(io, "| lmax | bsl: max abs d [m] | max abs d / range | seed d [m] | end d [m] | dz_bed host NH: max rms [m] | host ANT | rim | outside | rsl host NH: max rms | host ANT | dz_bed change since seed, host NH: max rms | host ANT |")
    println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in res
        s = r.stats
        println(io, @sprintf("| %d | %.4f | %.2e | %.4f | %.4f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |", r.lmax, r.bsl_max, r.bsl_rel, r.bsl_seed, r.bsl_end,
                s["dz_bed|host NH"].maxrms, s["dz_bed|host ANT"].maxrms, s["dz_bed|rim"].maxrms, s["dz_bed|outside"].maxrms,
                s["rsl|host NH"].maxrms, s["rsl|host ANT"].maxrms, s["dchg|host NH"].maxrms, s["dchg|host ANT"].maxrms))
    end
    if length(res) >= 2
        println(io, "\nRatio of each lmax to the next finer one (maximum over records, rms; values > 1: the difference falls with lmax):\n")
        println(io, "| from -> to | bsl max abs | dz_bed host NH | dz_bed host ANT | rsl host NH | rsl host ANT |")
        println(io, "|---|---|---|---|---|---|")
        for j in 1:length(res)-1
            a, b = res[j], res[j+1]
            println(io, @sprintf("| %d -> %d | %.2f | %.2f | %.2f | %.2f | %.2f |", a.lmax, b.lmax, a.bsl_max / b.bsl_max,
                    a.stats["dz_bed|host NH"].maxrms / b.stats["dz_bed|host NH"].maxrms, a.stats["dz_bed|host ANT"].maxrms / b.stats["dz_bed|host ANT"].maxrms,
                    a.stats["rsl|host NH"].maxrms / b.stats["rsl|host NH"].maxrms, a.stats["rsl|host ANT"].maxrms / b.stats["rsl|host ANT"].maxrms))
        end
    end
    s = String(take!(io)); print(s)
    md === nothing || write(md, s)
    return 0
end
exit(main(copy(ARGS)))
