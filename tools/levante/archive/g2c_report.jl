# G2c report (phase 4 WP7c): the spin-up grid. julia --project=tools tools/levante/g2c_report.jl <runs dir> [--ref NAME] [--md FILE]
# For every run folder e<T>_t<tol> of <runs dir>: the spin-up from the log (passes, summed cpu seconds of VILMA2, last mean rate v, the
# cap warning), the job wall (time.txt), and from geo_glob_vilma.nc (lossless, every 10 yr, record 1 = the seed at t0):
#   dz_bed at the seed, the mean rate |dz_bed(t_k) - dz_bed(t_{k-1})| / dt of the first intervals, the drift |dz_bed(t_k) - dz_bed(t_0)|,
#   each as the mean over the globe and over the host region (w_ref > 0.5: the ice footprints), area weights cos(lat);
#   and the differences against the reference run (default: the run with the largest cap and the smallest tolerance).
using NCDatasets, Printf, Statistics

function parse_log(f)
    passes = NamedTuple[]; label = ""; warn = ""
    for l in eachline(f)
        m = match(r"spin-up \[(.*?)\]", l)
        m !== nothing && (label = strip(m[1]))
        m = match(r"^\s+pass\s+(\d+):\s+(\S+) s\s+<dt>=\s+(\S+) yr \((\d+) sub\)\s+d\|z_bed\|=\s+(\S+) m\s+max=\s+(\S+) m\s+v=\s+(\S+) m/yr\s+v/tol=\s*(\S+)", l)
        if m !== nothing
            push!(passes, (label = label, n = parse(Int, m[1]), cpu = (c = tryparse(Float64, m[2]); c === nothing ? 1e4 : c),   # f7.2 overflows above 9999.99 s: counted as 1e4 (a lower bound)
                            dt = parse(Float64, m[3]) * parse(Int, m[4]),
                           nsub = parse(Int, m[4]), d = parse(Float64, m[5]), mx = parse(Float64, m[6]), v = parse(Float64, m[7]), r = parse(Float64, m[8])))
        end
        occursin("WARNING: spin-up", l) && (warn = strip(l))
    end
    passes, warn
end
function wall(dir)
    for l in eachline(joinpath(dir, "time.txt"))
        m = match(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.*)", l)
        if m !== nothing
            p = parse.(Float64, split(m[1], ":")); return length(p) == 3 ? p[1]*3600 + p[2]*60 + p[3] : p[1]*60 + p[2]
        end
    end
    NaN
end
function main(args)
    dir = args[1]; md = nothing; ref = nothing
    i = 2
    while i <= length(args)
        args[i] == "--md" && (md = args[i+1]; i += 1)
        args[i] == "--ref" && (ref = args[i+1]; i += 1)
        i += 1
    end
    names = sort(filter(n -> occursin(r"^e.*_t", n) && isfile(joinpath(dir, n, "geo_glob_vilma.nc")), readdir(dir)))
    key(n) = (m = match(r"^e(.*)_t(.*)$", n); (parse(Float64, replace(m[1], "e" => "e")), parse(Float64, replace(m[2], "m" => "-"))))
    names = sort(names, by = n -> (key(n)[1], -key(n)[2]))
    isnothing(ref) && (ref = names[argmax([key(n)[1] * 1e6 - key(n)[2] * 1e12 for n in names])])   # the most converged: the largest cap, the smallest tolerance
    D = Dict{String,Any}()
    for n in names
        ds = NCDataset(joinpath(dir, n, "geo_glob_vilma.nc"))
        D[n] = (t = ds["time"][:], dz = Array(ds["dz_bed"][:, :, :]), w = Array(ds["w_ref"][:, :]), lat = ds["lat"][:]); close(ds)
    end
    w = D[ref].w; lat = D[ref].lat
    A = [cosd(lat[j]) for i in axes(w, 1), j in axes(w, 2)]
    host = w .> 0.5
    wm(f, mask) = sum(f[mask] .* A[mask]) / sum(A[mask])
    glob = trues(size(w))
    out = IOBuffer()
    println(out, "# G2c: spin-up grid (", dir, ")\n")
    println(out, "Reference for the differences: **", ref, "**. Host region: w_ref > 0.5 (", count(host), " of ", length(host), " cells of the Gauss grid, ",
            @sprintf("%.1f", 100 * sum(A[host]) / sum(A)), " % of the area); means are cos(lat)-weighted.\n")
    println(out, "## Spin-up (log)\n")
    println(out, "| run | 1-D passes | 1-D cpu s | 1-D last v [m/yr] | full passes | full sub-steps | full cpu s | full relaxed [yr] | full last v [m/yr] (v/tol) | cap hit | job wall [min] |")
    println(out, "|---|---|---|---|---|---|---|---|---|---|---|")
    for n in names
        P, warn = parse_log(joinpath(dir, n, "out.out"))
        p1 = filter(p -> startswith(p.label, "1-D"), P); pf = filter(p -> startswith(p.label, "full"), P)
        f(x) = isempty(x) ? "-" : @sprintf("%.2e", x[end].v) * @sprintf(" (%.2g)", x[end].r)
        @printf(out, "| %s | %d | %.0f | %.2e | %d | %d | %.0f | %.0f | %s | %s | %.1f |\n", n, length(p1), sum(p.cpu for p in p1; init = 0.0),
                isempty(p1) ? NaN : p1[end].v, length(pf), sum(p.nsub for p in pf; init = 0), sum(p.cpu for p in pf; init = 0.0),
                sum(p.dt for p in pf; init = 0.0), f(pf), isempty(warn) ? "no" : "yes", wall(joinpath(dir, n)) / 60)
    end
    println(out, "\n## Bed at the seed and the first intervals (geo_glob_vilma.nc, dz_bed)\n")
    println(out, "Mean rate k = mean |dz_bed(t_k) - dz_bed(t_{k-1})| / 10 yr [mm/yr]; drift 5 = mean |dz_bed(t_5) - dz_bed(t_0)| [m]; seed = mean |dz_bed(t_0)| [m].\n")
    println(out, "| run | region | seed mean [m] | seed max [m] | rate 1 | rate 2 | rate 3 | rate 4 | rate 5 | drift 5 mean [m] | drift 5 max [m] |")
    println(out, "|---|---|---|---|---|---|---|---|---|---|---|")
    for n in names, (rn, mask) in (("globe", glob), ("host", host))
        z = D[n].dz; dt = D[n].t[2] - D[n].t[1]
        rate(k) = wm(abs.(z[:, :, k+1] .- z[:, :, k]) ./ dt .* 1e3, mask)
        @printf(out, "| %s | %s | %.3f | %.2f | %s | %.3f | %.2f |\n", n, rn, wm(abs.(z[:, :, 1]), mask), maximum(abs.(z[:, :, 1][mask])),
                join([@sprintf("%.3f", rate(k)) for k in 1:5], " | "), wm(abs.(z[:, :, 6] .- z[:, :, 1]), mask), maximum(abs.((z[:, :, 6] .- z[:, :, 1])[mask])))
    end
    println(out, "\n## Difference to the reference run\n")
    println(out, "mean |dz_bed(run) - dz_bed(ref)| [m] at the seed and at t_5 (50 yr), and of the rate of interval 1 [mm/yr] (the leftover relaxation of the run against the reference under the same load history).\n")
    println(out, "| run | region | seed [m] | seed diff / ref seed mean | t5 [m] | seed max [m] | rate 1 diff [mm/yr] |")
    println(out, "|---|---|---|---|---|---|---|")
    zr = D[ref].dz
    for n in names, (rn, mask) in (("globe", glob), ("host", host))
        z = D[n].dz; dt = D[n].t[2] - D[n].t[1]
        r1 = (z[:, :, 2] .- z[:, :, 1]) ./ dt .* 1e3 .- (zr[:, :, 2] .- zr[:, :, 1]) ./ dt .* 1e3
        @printf(out, "| %s | %s | %.3f | %.2f | %.3f | %.2f | %.3f |\n", n, rn, wm(abs.(z[:, :, 1] .- zr[:, :, 1]), mask), wm(abs.(z[:, :, 1] .- zr[:, :, 1]), mask) / wm(abs.(zr[:, :, 1]), mask), wm(abs.(z[:, :, 6] .- zr[:, :, 6]), mask),
                maximum(abs.((z[:, :, 1] .- zr[:, :, 1])[mask])), wm(abs.(r1), mask))
    end
    s = String(take!(out)); print(s)
    md === nothing || write(md, s)
end
main(ARGS)
