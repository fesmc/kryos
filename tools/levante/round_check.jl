# G4 rounding check: a run with rounded output against the same run with [kryos] write_lossless.
#   julia --project=tools tools/levante/round_check.jl <rounded run> <lossless run> <kryos input folder> [--md FILE]
# (1) restart bundles: every variable of every file of every restart-* folder identical (kryos-compare --restarts on the bundles
#     is the stricter tool; here a bit comparison of the same files);
# (2) output files (top level of the run folder, not the bundles): every variable of a file present in both runs, tolerance from the
#     `prec` column of the variable tables of the input folder (|a-b| <= prec_abs + prec_rel |b|; no entry = lossless: identical);
#     reported: max error / tolerance per file, the worst variables, the number of variables identical, in tolerance, outside.
using NCDatasets, Printf

function load_prec(dir)
    tol = Dict{String,Tuple{Float64,Float64}}()
    for f in readdir(dir; join=true)
        (endswith(f, ".md") && occursin("variables", basename(f))) || continue
        for l in eachline(f)
            c = strip.(split(strip(l), "|"))
            filter!(!isempty, c)
            length(c) == 6 || continue
            c[1] == "id" && continue
            all(isdigit, c[1]) || continue
            name, p = c[2], c[6]
            pa, pr = 0.0, 0.0
            if p != "-"
                for part in split(p, ",")
                    part = strip(part)
                    if endswith(part, "r"); pr = parse(Float64, part[1:end-1]) else pa = parse(Float64, part) end
                end
            end
            old = get(tol, name, nothing)
            tol[name] = old === nothing ? (pa, pr) : (max(old[1], pa), max(old[2], pr))   # a name in two tables: the looser
        end
    end
    tol
end

function main(args)
    pos = filter(a -> !startswith(a, "--"), args)
    md = nothing
    i = findfirst(==("--md"), args); i === nothing || (md = args[i+1]; pos = filter(!=(md), pos))
    rdir, ldir, idir = pos
    tol = load_prec(idir)
    io = IOBuffer()
    println(io, "# G4 rounding check: $(basename(dirname(rdir)))/$(basename(rdir)) (rounded) vs $(basename(ldir)) (lossless)\n")
    # (1) bundles
    nb = 0; nbad = 0; nvar = 0
    for b in sort(filter(d -> startswith(d, "restart-") && isdir(joinpath(rdir, d)), readdir(rdir)))
        isdir(joinpath(ldir, b)) || (println(io, "- bundle $b missing in the lossless run"); nbad += 1; continue)
        for f in sort(filter(endswith(".nc"), readdir(joinpath(rdir, b))))
            fl = joinpath(ldir, b, f); isfile(fl) || (nbad += 1; println(io, "- $b/$f missing in the lossless run"); continue)
            NCDataset(joinpath(rdir, b, f)) do a; NCDataset(fl) do c
                for v in keys(a)
                    haskey(c, v) || (nbad += 1; println(io, "- $b/$f:$v absent in the lossless run"); continue)
                    x = Array(a[v]); y = Array(c[v]); nvar += 1
                    eltype(x) <: Number || continue
                    (size(x) == size(y) && all(((p, q),) -> p === q || (p isa AbstractFloat && isnan(p) && isnan(q)) || p == q, zip(x, y))) ||
                        (nbad += 1; println(io, "- $b/$f:$v DIFFERS"))
                end
            end; end
        end
        nb += 1
    end
    println(io, "## Restart bundles\n\n$(nb) bundles, $(nvar) variables compared, $(nbad) differing or missing -> ", nbad == 0 && nb > 0 ? "IDENTICAL" : "NOT IDENTICAL", "\n")
    # (2) outputs
    println(io, "## Output files\n\n| file | variables | identical | in tolerance (rounded) | outside | max error / tolerance | worst variable |\n|---|---|---|---|---|---|---|")
    worst_all = 0.0; nout_all = 0
    for f in sort(filter(endswith(".nc"), readdir(rdir)))
        fl = joinpath(ldir, f); isfile(fl) || continue
        nv = 0; nid = 0; nin = 0; nout = 0; wr = 0.0; wn = ""
        NCDataset(joinpath(rdir, f)) do a; NCDataset(fl) do c
            for v in keys(a)
                haskey(c, v) || continue
                x = a[v]; eltype(x) <: Union{Number,Missing} || continue
                xa = Array(x); ya = Array(c[v])
                (xa isa AbstractArray{<:Number} || eltype(xa) <: Union{Missing,Number}) || continue
                size(xa) == size(ya) || continue
                nv += 1
                pa, pr = get(tol, v, (0.0, 0.0))
                r = 0.0; same = true
                for (p, q) in zip(xa, ya)
                    (ismissing(p) || ismissing(q)) && continue
                    (isnan(p) && isnan(q)) && continue
                    p == q && continue
                    same = false
                    t = pa + pr * abs(q)
                    r = max(r, t > 0 ? abs(p - q) / t : Inf)
                end
                if same; nid += 1 elseif r <= 1 + 1e-9; nin += 1 else nout += 1 end
                if r > wr; wr = r; wn = v end
            end
        end; end
        nv == 0 && continue
        @printf(io, "| %s | %d | %d | %d | %d | %.3g | %s |\n", f, nv, nid, nin, nout, wr, wn)
        worst_all = max(worst_all, wr); nout_all += nout
    end
    @printf(io, "\nmax error / tolerance over all files: %.3g -> %s\n", worst_all, worst_all <= 1 + 1e-9 ? "PASS (<= 1)" : "FAIL")
    s = String(take!(io)); print(s); md === nothing || write(md, s)
    return 0
end
exit(main(ARGS))
