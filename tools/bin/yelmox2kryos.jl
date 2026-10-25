#!/usr/bin/env julia
# yelmox par file -> kryos par file, fully expanded, with its accounting report
# (regression-design.md §2; KryosTools.Convert).
#
#   julia --project=tools tools/bin/yelmox2kryos.jl --yelmox $YELMOXROOT \
#       --par yelmox/yelmox_Greenland.nml -p ctrl.time_end=1000 ... \
#       --out par/gate/grl.nml --report par/gate/grl_conversion.md
#
# Exit 0 converted, 1 refused (outside phase 1) or invalid, 2 usage.

using KryosTools
using KryosTools.Convert

const USAGE = """
usage: yelmox2kryos.jl --yelmox ROOT --par FILE --out FILE [options]

  --yelmox ROOT    yelmox checkout (holds .runme/info.json and input/)
  --par FILE       the par file: absolute, relative to the cwd, or to ROOT. Either the
                   file as staged by runme in a reference run directory (-p already
                   applied; give no -p) or a base file with the same -p list
  --out FILE       the kryos par file to write
  --report FILE    the accounting report (markdown)
  -p g.k=v         a runme override (repeatable; g.k==a,b for a vector); or
  -p "g.k=v g.k=v" several in one argument (a gate's plist_string)
  --kryos ROOT     the kryos checkout whose templates validate the file
                   (default: the checkout of this script)
  --keep-library   do not copy the library groups in: reference them in [kryos] library
  --spec FILE      the mapping table (default: tools/data/yelmox2kryos.toml)
"""

function main(args)
    opts = Dict{String,Any}("p" => String[], "keep" => false)
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("-h", "--help")
            print(USAGE); return 0
        elseif a == "--keep-library"
            opts["keep"] = true; i += 1
        elseif a == "-p"
            i < length(args) || (print(stderr, USAGE); return 2)
            append!(opts["p"], split(args[i+1]))
            i += 2
        elseif a in ("--yelmox", "--par", "--out", "--report", "--kryos", "--spec")
            i < length(args) || (print(stderr, USAGE); return 2)
            opts[a[3:end]] = args[i+1]
            i += 2
        else
            println(stderr, "unknown argument $a"); print(stderr, USAGE); return 2
        end
    end
    for k in ("yelmox", "par", "out")
        haskey(opts, k) || (println(stderr, "missing --$k"); print(stderr, USAGE); return 2)
    end
    kw = Dict{Symbol,Any}(:yelmox => opts["yelmox"], :par => opts["par"], :overrides => String.(opts["p"]),
                          :keep_library => opts["keep"])
    haskey(opts, "kryos") && (kw[:kryos] = opts["kryos"])
    haskey(opts, "spec") && (kw[:spec] = opts["spec"])
    try
        conv = convert_yelmox(; kw...)
        write_par(opts["out"], conv)
        haskey(opts, "report") && write_report(opts["report"], conv)
        s = account_summary(conv)
        println("wrote ", opts["out"], haskey(opts, "report") ? " and " * opts["report"] : "")
        println("keys: mapped ", s.mapped, ", default-equal ", s.default, ", dropped ", s.dropped, " (total ", s.total, ")")
        for w in unique(conv.warnings)
            println("warning: ", w)
        end
        return 0
    catch e
        if e isa ConversionError
            println(stderr, e.msg)
            return 1
        end
        rethrow()
    end
end

exit(main(ARGS))
