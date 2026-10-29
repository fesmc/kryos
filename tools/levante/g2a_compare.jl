# G2a: kryos_out.nc (tests/g2a_cap.f90) against vilma_out.nc (vilma.x, same load); prints per record and variable
# max |kryos - standalone|, max |standalone| and their ratio, and the maximum over all records.
#   julia --project=tools tools/levante/g2a_compare.jl <run folder> [--md FILE]
using NCDatasets, Printf

function main(args)
    dir = args[1]
    md = nothing
    i = findfirst(==("--md"), args); i === nothing || (md = args[i+1])
    k = NCDataset(joinpath(dir, "kryos_out.nc")); v = NCDataset(joinpath(dir, "vilma_out.nc"))
    pairs = [("dz_bed", "dz_bed"), ("dz_ss", "dz_ss"), ("C_ocn", "C_ocean"), ("dphi", "dphi"), ("bsl", "bsl"), ("bsl_ive", "bsl_ive")]
    nt = length(k["time"][:])
    io = IOBuffer()
    println(io, "# G2a: adapter (geo_glob) vs vilma.x\n")
    println(io, "records: ", nt, " (0 = after the spin-up and seed), standalone records: ", length(v["time"][:]), "\n")
    println(io, "| variable | max over records of max abs diff | max abs value (standalone) | ratio | first record with a difference |")
    println(io, "|---|---|---|---|---|")
    for (kn, vn) in pairs
        haskey(v, vn) || (println(io, "| ", vn, " | absent in vilma_out.nc | | | |"); continue)
        a = Array(k[kn][:, :, :] .* 1.0); b = Array(v[vn][:, :, :] .* 1.0)
        a = ndims(a) == 1 ? reshape(a, 1, 1, :) : a
        b = ndims(b) == 1 ? reshape(b, 1, 1, :) : b
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
        a = k[kn][:, :, :]; b = v[vn][:, :, :]
        println(io, "- ", kn, ": ", join([@sprintf("%.1e", maximum(abs.(a[:, :, r] .- b[:, :, r])) / max(maximum(abs.(b[:, :, r])), eps())) for r in 1:min(size(a, 3), size(b, 3))], " "))
    end
    s = String(take!(io)); print(s)
    md === nothing || write(md, s)
    return 0
end
exit(main(ARGS))
