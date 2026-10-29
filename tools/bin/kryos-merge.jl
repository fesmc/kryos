#!/usr/bin/env julia
# Single-domain kryos par files -> one par file for a run with all their domains (phase 3 WP1;
# KryosTools.Convert.merge_par).
#
#   julia --project=tools tools/bin/kryos-merge.jl par/kryos_North.nml par/kryos_Antarctica.nml \
#       --out par/kryos_NH-ANT.nml --report par/kryos_NH-ANT_merge.md \
#       --resolve clock.time_end=15e3 --resolve clock.time_const=0.0 --resolve 'barysealevel.sl_path="none"'
#
# The order of the files is the order of [kryos] domains. Exit 0 merged, 1 stopped (a run-wide key
# differs and is not resolved, ...), 2 usage.

using KryosTools
using KryosTools.Convert

const USAGE = """
usage: kryos-merge.jl FILE FILE [FILE ...] --out FILE [options]

  --out FILE        the merged par file
  --report FILE     the merge report (markdown): renamed groups, repointed keys, resolved keys
  --resolve g.k=v   the value of a run-wide key that differs between the files (repeatable; v as written
                    in a par file, strings quoted)
  --swap            domains in the reverse order of the files (the order-swap runs of the gates)
"""

function main(args)
    files = String[]
    res = Dict{String,String}()
    opts = Dict{String,Any}("swap" => false)
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("-h", "--help")
            print(USAGE); return 0
        elseif a == "--swap"
            opts["swap"] = true; i += 1
        elseif a == "--resolve"
            i < length(args) || (print(stderr, USAGE); return 2)
            kv = split(args[i+1], '='; limit=2)
            length(kv) == 2 || (println(stderr, "--resolve needs group.key=value"); return 2)
            res[String(kv[1])] = String(kv[2])
            i += 2
        elseif a in ("--out", "--report")
            i < length(args) || (print(stderr, USAGE); return 2)
            opts[a[3:end]] = args[i+1]; i += 2
        elseif startswith(a, "--")
            println(stderr, "unknown argument $a"); print(stderr, USAGE); return 2
        else
            push!(files, a); i += 1
        end
    end
    (length(files) >= 2 && haskey(opts, "out")) || (print(stderr, USAGE); return 2)
    opts["swap"] && reverse!(files)
    try
        m = merge_par(files; resolve=res)
        write_merged(opts["out"], m)
        haskey(opts, "report") && open(io -> merge_report(io, m, opts["out"]), opts["report"], "w")
        println("wrote ", opts["out"], ": domains ", join(m.domains, ", "), "; ", length(m.renamed), " groups renamed, ",
                length(m.repointed), " keys repointed, ", length(m.resolved), " resolved")
        return 0
    catch e
        e isa ConversionError || rethrow()
        println(stderr, sprint(showerror, e))
        return 1
    end
end

exit(main(ARGS))
