#!/usr/bin/env julia
# globe_ice.nc against the merge computed offline (phase4-plan 5.4, gate G1b (4); KryosTools.GlobeIce).
#
#   julia --project=tools tools/bin/kryos-globe-ice.jl GLOBE_ICE_NC INPUTS_NC [--tol-abs X] [--rtol X]
#
# INPUTS_NC holds the global attribute `domains` (comma separated), `map_<d>` (the coupler's
# scrip-con map of the domain's ice grid to the background grid), and the variables `area_bg`,
# `H_ext`, and per domain `H_<d>`, `area_<d>`, optionally `owner_<d>` (module docstring).
# --tol-abs: the rounding of the file (the `prec` of the variable table; 0 for a lossless run).
# Exit 0 equal, 1 different.

using KryosTools
using KryosTools.GlobeIce

function main(args)
    pos = String[]
    tol_abs = 0.0
    rtol = 1e-12
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--tol-abs"
            tol_abs = parse(Float64, args[i+1]); i += 2
        elseif a == "--rtol"
            rtol = parse(Float64, args[i+1]); i += 2
        elseif startswith(a, "--")
            println(stderr, "unknown option $a"); return 2
        else
            push!(pos, a); i += 1
        end
    end
    length(pos) == 2 || (println(stderr, "usage: kryos-globe-ice.jl GLOBE_ICE_NC INPUTS_NC [--tol-abs X] [--rtol X]"); return 2)
    r = check_globe_ice(pos[1], pos[2]; tol_abs = tol_abs, rtol = rtol)
    println("globe_ice.nc vs the offline merge of ", r.n_domains, " domains: max |H_ice| error ", r.max_err_H,
            ", max |cov| error ", r.max_err_cov, " -> ", r.ok ? "equal" : "DIFFERENT")
    return r.ok ? 0 : 1
end

exit(main(ARGS))
