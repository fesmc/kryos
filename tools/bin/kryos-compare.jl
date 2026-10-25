#!/usr/bin/env julia
# Compare two run directories bit for bit (regression-design.md §3; KryosTools.Compare).
#
#   julia --project=tools tools/bin/kryos-compare.jl REF RUN \
#       --map tools/gate/map_phase1.toml --deviations tools/gate/deviations.toml \
#       --gate grl [--restarts] [--json out.json] [--md out.md] [--plot VAR ...]
#
# REF and RUN are run directories; `--map identity` (the default) compares two
# runs of one model. `--chunks` compares a chunk run against its continuous
# run (REF = cont, RUN = the folder holding c00, c01, ...; needs --gate).
# Exit 0 pass, 1 fail, 2 error (missing file or record, unreadable map).

using KryosTools
using KryosTools.Compare

const USAGE = """
usage: kryos-compare.jl REF RUN [options]

  --map FILE|identity   mapping table (default identity: same file names, class E)
  --deviations FILE     registry of accepted differences (mapped mode)
  --gate NAME|FILE      the gate (tools/gate/<NAME>.toml): instance names, restart times, chunks
  --restarts            also compare the restart bundles
  --chunks              REF = cont, RUN = chunk run (c00, c01, ...): the restart round trip
  --json FILE           the report as JSON
  --md FILE             the report as markdown (default: standard output)
  --plot VAR            figures of a variable (repeatable; needs CairoMakie and the plot extension)
  --plot-dir DIR        where the figures go (default .)
"""

function main(args)
    pos = String[]
    opts = Dict{String,Any}("plot" => String[], "restarts" => false, "chunks" => false)
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("-h", "--help")
            print(USAGE); return 0
        elseif a in ("--restarts", "--chunks")
            opts[a[3:end]] = true; i += 1
        elseif a == "--plot"
            i < length(args) || (print(stderr, USAGE); return 2)
            push!(opts["plot"], args[i+1]); i += 2
        elseif a in ("--map", "--deviations", "--gate", "--json", "--md", "--plot-dir")
            i < length(args) || (print(stderr, USAGE); return 2)
            opts[a[3:end]] = args[i+1]; i += 2
        elseif startswith(a, "--")
            println(stderr, "unknown option $a"); print(stderr, USAGE); return 2
        else
            push!(pos, a); i += 1
        end
    end
    length(pos) == 2 || (print(stderr, USAGE); return 2)
    ref, run = pos
    gate = get(opts, "gate", "")
    rep = if opts["chunks"]
        isempty(gate) && (println(stderr, "--chunks needs --gate"); return 2)
        compare_chunks(ref, run; gate=gate, restarts=true)
    else
        compare(ref, run; map=get(opts, "map", "identity"), gate=gate,
                deviations=get(opts, "deviations", ""), restarts=opts["restarts"])
    end
    if haskey(opts, "md")
        open(io -> report_markdown(io, rep), opts["md"], "w")
    else
        report_markdown(stdout, rep)
    end
    if haskey(opts, "json")
        open(opts["json"], "w") do io
            write(io, report_json(rep))
        end
    end
    if !isempty(opts["plot"])
        try
            @eval using CairoMakie
            for v in opts["plot"]
                for f in Base.invokelatest(plot_variable, rep, v; dir=get(opts, "plot-dir", "."))
                    println("figure: ", f)
                end
            end
        catch e
            println(stderr, "plot: ", sprint(showerror, e))
        end
    end
    return exit_code(rep)
end

exit(main(ARGS))
