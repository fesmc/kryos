# G2a: kryos_out.nc (tests/g2a_cap.f90) against vilma_out.nc (vilma.x, same load); prints per record and variable
# max |kryos - standalone|, max |standalone| and their ratio, and the maximum over all records. The files are aligned by
# their lat axes (kryos: the host grid, ascending; vilma.x: VILMA2's layout, north first); other axes stop.
#   julia --project=tools tools/levante/g2a_compare.jl <run folder> [--md FILE] [--ll]
# --ll: the independent variant, vilma_out_ll.nc (vilma.x remapped the 0.5 degree lon-lat inputs itself). No bitwise
# expectation: the tolerance is the remap error of the load, the relative L2 difference of the Gauss-grid ice h_ice of
# vilma_out_ll.nc against vilma_out.nc (the exact Gauss load); dz_bed and dz_ss of kryos against vilma_out_ll.nc must differ
# by no more than it (relative L2, per record). The same against the standalone mirrored north-south shows the margin.
using NCDatasets, Printf

rel2(a, b) = sqrt(sum(abs2, a .- b) / max(sum(abs2, b), eps()))

function main(args)
    dir = args[1]
    md = nothing
    i = findfirst(==("--md"), args); i === nothing || (md = args[i+1])
    ll = "--ll" in args
    k = NCDataset(joinpath(dir, "kryos_out.nc")); v = NCDataset(joinpath(dir, ll ? "vilma_out_ll.nc" : "vilma_out.nc"))
    latk = Float64.(k["lat"][:]); latv = Float64.(v["lat"][:])
    maximum(abs.(Float64.(k["lon"][:]) .- Float64.(v["lon"][:]))) <= 1e-6 || error("lon axes differ")
    flip = if maximum(abs.(latk .- latv)) <= 1e-6
        false
    elseif maximum(abs.(latk .- reverse(latv))) <= 1e-6
        true
    else
        error("lat axes differ (not the same Gauss grid in either order)")
    end
    sa(x) = ndims(x) == 1 ? reshape(x, 1, 1, :) : x
    vrows(x) = (y = sa(Float64.(Array(x))); flip ? y[:, end:-1:1, :] : y)       # vilma.x field in the kryos row order
    pairs = [("dz_bed", "dz_bed"), ("dz_ss", "dz_ss"), ("C_ocn", "C_ocean"), ("dphi", "dphi"), ("bsl", "bsl"), ("bsl_ive", "bsl_ive")]
    nt = length(k["time"][:])
    io = IOBuffer()
    println(io, "# G2a: adapter (geo_glob) vs vilma.x", ll ? " on the lon-lat inputs (remap_input = .true.)" : "", "\n")
    println(io, "records: ", nt, " (0 = after the spin-up and seed), standalone records: ", length(v["time"][:]), "\n")
    println(io, "lat axes: kryos ", latk[1] < latk[end] ? "ascending" : "descending", ", vilma.x ", latv[1] < latv[end] ? "ascending" : "descending", flip ? " (vilma.x rows reversed for the comparison)" : "", "\n")
    println(io, "| variable | max over records of max abs diff | max abs value (standalone) | ratio | first record with a difference |")
    println(io, "|---|---|---|---|---|")
    for (kn, vn) in pairs
        haskey(v, vn) || (println(io, "| ", vn, " | absent in the standalone file | | | |"); continue)
        a = sa(Float64.(Array(k[kn]))); b = vrows(v[vn])
        n = min(size(a, 3), size(b, 3))
        dmax = 0.0; bmax = 0.0; first = 0
        for r in 1:n
            d = maximum(abs.(a[:, :, r] .- b[:, :, r])); dmax = max(dmax, d); bmax = max(bmax, maximum(abs.(b[:, :, r])))
            d > 0 && first == 0 && (first = r - 1)
        end
        @printf(io, "| %s | %.3e | %.3e | %.3e | %s |\n", kn, dmax, bmax, bmax == 0 ? NaN : dmax / bmax, first == 0 && dmax == 0 ? "none (bit identical)" : string(first))
    end
    println(io, "\nper record, dz_bed and dz_ss (max abs diff / max abs standalone):\n")
    for (kn, vn) in pairs[1:2]
        a = sa(Float64.(Array(k[kn]))); b = vrows(v[vn])
        println(io, "- ", kn, ": ", join([@sprintf("%.1e", maximum(abs.(a[:, :, r] .- b[:, :, r])) / max(maximum(abs.(b[:, :, r])), eps())) for r in 1:min(size(a, 3), size(b, 3))], " "))
    end
    pass = true
    if ll
        g = NCDataset(joinpath(dir, "vilma_out.nc"))
        hl = sa(Float64.(Array(v["h_ice"]))); hg = sa(Float64.(Array(g["h_ice"])))   # both in VILMA2's layout
        n = min(size(hl, 3), size(hg, 3), nt)
        tol = [rel2(hl[:, :, r], hg[:, :, r]) for r in 1:n]
        println(io, "\n## Tolerance from the remap error\n")
        println(io, "relative L2 difference of the load (h_ice, remapped lon-lat against the exact Gauss load), per record: ",
                join([@sprintf("%.2e", x) for x in tol], " "), "\n")
        println(io, "| variable | relative L2 kryos vs lon-lat standalone, per record | max | max tolerance | within | mirrored standalone, min |")
        println(io, "|---|---|---|---|---|---|")
        for (kn, vn) in pairs[1:2]
            a = sa(Float64.(Array(k[kn]))); b = vrows(v[vn]); bm = b[:, end:-1:1, :]
            e = [rel2(a[:, :, r], b[:, :, r]) for r in 1:n]; em = [rel2(a[:, :, r], bm[:, :, r]) for r in 1:n]
            ok = all(e[r] <= tol[r] for r in 1:n if tol[r] > 0); pass &= ok
            @printf(io, "| %s | %s | %.2e | %.2e | %s | %.2e |\n", kn, join([@sprintf("%.2e", x) for x in e], " "), maximum(e),
                    maximum(tol), ok ? "yes" : "NO", minimum(em))
        end
        println(io, "\nverdict: ", pass ? "PASS" : "FAIL", " (relative L2 of dz_bed and dz_ss within the load's remap error at every record)")
        close(g)
    end
    s = String(take!(io)); print(s)
    md === nothing || write(md, s)
    return pass ? 0 : 1
end
exit(main(ARGS))
