#!/usr/bin/env julia
# Kryos-against-kryos gates (phase 3, regression-design.md §1.6; KryosTools.Gate).
#
#   julia --project=tools tools/bin/kryos-gate-kk.jl check GATE RUNS [--md FILE] [--only NAME ...]
#   julia --project=tools tools/bin/kryos-gate-kk.jl stage GATE RUN
#   julia --project=tools tools/bin/kryos-gate-kk.jl list GATE
#
# GATE is a name of tools/gate (kind = "kryos") or a path to a gate file. `check` runs every [[check]] of the
# gate over the run folders in RUNS (exit 0 all pass, 1 a failure, 2 an error). `stage` prints, for the run RUN,
# `PAR=<par/gate/NAME.nml name>` and one line `KP=<override>` per kryos -p override, and for a chunk run `CHUNK=<t0 t1 t2>` for the staging script.

using KryosTools
using KryosTools.Compare
using KryosTools.Gate

const USAGE = """
usage: kryos-gate-kk.jl check GATE RUNS [--md FILE] [--only NAME ...]
       kryos-gate-kk.jl stage GATE RUN
       kryos-gate-kk.jl list GATE
"""

function main(args)
    length(args) >= 2 || (print(stderr, USAGE); return 2)
    cmd, gate = args[1], args[2]
    g = load_gate(gate)
    if cmd == "list"
        for r in kk_runs(g)
            println("run   ", r["name"], ": par/gate/", r["par"], ".nml ", get(r, "kp", String[]), get(r, "swap", false) ? " (swap)" : "")
        end
        for c in kk_checks(g)
            println("check ", get(c, "type", "?"), ": ", join(["$k=$v" for (k, v) in c if k != "type"], " "))
        end
        return 0
    elseif cmd == "stage"
        length(args) == 3 || (print(stderr, USAGE); return 2)
        par, p = kk_stage(g, args[3])
        println("PAR=", par)
        foreach(x -> println("KP=", x), p)
        run = only(r for r in kk_runs(g) if r["name"] == args[3])
        if get(run, "chunk", false)
            ts = Float64.(g["chunks"])
            println("CHUNK=", join(ts, " "))
        end
        return 0
    elseif cmd == "check"
        length(args) >= 3 || (print(stderr, USAGE); return 2)
        runs = args[3]
        opts = Dict{String,Any}("only" => String[])
        i = 4
        while i <= length(args)
            a = args[i]
            if a == "--md"
                i < length(args) || (print(stderr, USAGE); return 2)
                opts[a[3:end]] = args[i+1]; i += 2
            elseif a == "--only"
                i < length(args) || (print(stderr, USAGE); return 2)
                push!(opts["only"], args[i+1]); i += 2
            else
                println(stderr, "unknown argument $a"); return 2
            end
        end
        rs = run_kk_checks(g, runs; only=isempty(opts["only"]) ? nothing : opts["only"])
        if haskey(opts, "md")
            open(io -> kk_markdown(io, g, rs; runs=runs), opts["md"], "w")
        end
        kk_markdown(stdout, g, rs; runs=runs)
        out = kk_outcome(rs)
        return out == "PASS" ? 0 : out == "FAIL" ? 1 : 2
    end
    print(stderr, USAGE)
    return 2
end

exit(main(ARGS))
