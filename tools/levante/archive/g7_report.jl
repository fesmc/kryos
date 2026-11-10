# G7 report (phase 4 WP7c): cost table. julia --project=tools tools/levante/g7_report.jl <g7 root with dt10 dt20 dt50 dt100> <lossless run (G4, sci-vilma2)> <rounded run (G4)> [--md FILE]
# Per geo cadence: seconds per 10-yr step by component (timing.txt, cumulative columns over the 40 steps), the VILMA2 share
# (the column `forcing` holds the gather, the VILMA2 call and the forcing), fit of the cost per call against dt, job wall, node-hours
# (DKRZ job summary in out.out), peak RSS; restart bundle size and write time (jump of the output column at the restart step);
# output sizes of the rounded and the lossless run (400 yr, 2D every 50 yr).
using Printf, Statistics
function timing(f)
    L = readlines(f); h = split(strip(L[1])); rows = [parse.(Float64, split(strip(l))) for l in L[2:end] if !isempty(strip(l))]
    h, reduce(hcat, rows)'    # rows x columns
end
function summary(dir)
    el = NaN; nh = NaN; rss = NaN
    for l in eachline(joinpath(dir, "out.out"))
        m = match(r"Elapsed time\s*:\s*(\d+):(\d+):(\d+)", l); m === nothing || (el = parse(Int, m[1]) * 3600 + parse(Int, m[2]) * 60 + parse(Int, m[3]))
        m = match(r"^\* 0\s+\| kryos.x\s+(\S+)\s+(\d+)K", l); m === nothing || (nh = parse(Float64, m[1]); rss = parse(Float64, m[2]) / 1e6)
    end
    el, nh, rss
end
dirsize(d) = (s = 0; for (r, _, fs) in walkdir(d), f in fs; s += filesize(joinpath(r, f)); end; s)
function main(args)
    root, ll, rr = args[1:3]; md = nothing
    i = findfirst(==("--md"), args); i === nothing || (md = args[i+1])
    io = IOBuffer()
    println(io, "# G7: cost of the science case with geo = vilma2 (S8, 64 threads, exclusive node)\n")
    println(io, "## Per 10-yr step by component and geo cadence (400 yr = 40 steps from -80 ka, rounded output, bundle at -79800)\n")
    println(io, "| dt geo [yr] | calls | forcing (gather + VILMA2 + forcing) [s/step] | per call [s] | ice [s/step] | output [s/step] | climate + surface + shelf [s/step] | total [s/step] | job wall [min] | node-h | peak RSS [GB] |\n|---|---|---|---|---|---|---|---|---|---|---|")
    pts = Tuple{Float64,Float64}[]
    for dt in (10, 20, 50, 100)
        d = joinpath(root, "dt$dt"); isdir(d) || continue
        h, T = timing(joinpath(d, "timing.txt")); n = size(T, 1) - 1
        c(name) = T[end, findfirst(==(name), h)] / n
        calls = 400 / dt
        el, nh, rss = summary(d)
        @printf(io, "| %d | %d | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f | %.1f | %.2f | %.1f |\n", dt, calls, c("forcing"), c("forcing") * n / calls, c("ice"), c("output"),
                c("climate") + c("surface") + c("shelf"), c("total"), el / 60, nh, rss)
        push!(pts, (dt, c("forcing") * n / calls))
    end
    if length(pts) >= 3
        X = [ones(length(pts)) first.(pts)]; b = X \ last.(pts)
        @printf(io, "\nCost per call against dt: %.2f s + %.3f s/yr * dt (least squares over %d cadences): the integration cost is proportional to the simulated time (VILMA2's explicit step is about 10 yr in the 3D model: one sub-step per 10 yr), a coarser cadence saves only the fixed part of the call.\n", b[1], b[2], length(pts))
    end
    # restart
    d = joinpath(root, "dt10"); h, T = timing(joinpath(d, "timing.txt"))
    o = T[:, findfirst(==("output"), h)]; dO = diff(o); k = argmax(dO)
    rb = joinpath(d, "restart--79800")
    @printf(io, "\n## Restart\n\nBundle %s: %.0f MB; the output column jumps by %.1f s at the step after the restart (against %.2f s of a step without output, %.2f s with a 2D output).\n", basename(rb), dirsize(rb) / 1e6, dO[k], median(dO), maximum(dO[dO .< dO[k] * 0.5]))
    for f in sort(readdir(rb)); p = joinpath(rb, f); (isfile(p) && filesize(p) > 5e6) && @printf(io, "- %s %.1f MB\n", f, filesize(p) / 1e6); end
    # output sizes
    println(io, "\n## Output of 400 yr (G4 runs: 2D and streams every 50 yr), rounded against lossless\n")
    println(io, "| file | rounded [MB] | lossless [MB] | ratio |\n|---|---|---|---|")
    fs = sort(filter(f -> endswith(f, ".nc") && isfile(joinpath(rr, f)) && isfile(joinpath(ll, f)), readdir(rr)))
    tr = 0.0; tl = 0.0
    for f in fs
        a = filesize(joinpath(rr, f)) / 1e6; b = filesize(joinpath(ll, f)) / 1e6; tr += a; tl += b
        (a > 3 || b > 3) && @printf(io, "| %s | %.1f | %.1f | %.2f |\n", f, a, b, b / a)
    end
    @printf(io, "| all %d files | %.1f | %.1f | %.2f |\n", length(fs), tr, tl, tl / tr)
    s = String(take!(io)); print(s); md === nothing || write(md, s)
end
main(ARGS)
