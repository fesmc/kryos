#!/usr/bin/env julia
# Par files of the cf calibration (docs/dev/calib-cf.md), as a text edit by group (nothing else changes, every other key and comment is the base's).
#   julia tools/bin/kryos-calib-par.jl sci IN OUT
#       IN: par/kryos_NH-ANT.nml or par/kryos_NH-ANT_vilma2.nml. [ytill_ant]: scale_zb = 2 (exp, z0/z1 = -400/0, cf_min = 0.01 as in the file), cf_ref = 0.2;
#       [ytill_nor]: cf_ref = 0.2 (exp -100/100, cf_min = 0.001 as in the file).
#   julia tools/bin/kryos-calib-par.jl run IN OUT DOMAIN CASE
#       IN: par/kryos_NH-ANT_calib.nml (FastIsostasy). One active domain (Antarctica | North), 30 kyr, frozen forcing (clock "const").
#       CASE pd:  t 0 -> 30e3, alpha = 0 (present-day climate), barysealevel const 0, orbit at 0, start from the hub reference (no experiment).
#       CASE lgm: t 30e3 -> 60e3, restarted from the end of the pd run ([kryos] restart is set by the run script), alpha = 1 (the LGM terms),
#                 barysealevel from a constant -120 m series (input/sl_const_lgm120.dat), orbit fixed_time -21e3.
#       Antarctica: the paleo terms lgm + ocean (DJF as the second PMIP3 tas) as of WP7c kryos-ant-par.jl, in both cases (alpha = 0 is the present day);
#       North: its terms are the file's (lgm + ocean), the ICE-6G start of the base is removed.
function groups(lines)
    out = Tuple{String,Int,Int}[]
    i = 1
    while i <= length(lines)
        m = match(r"^&(\w+)", lines[i])
        if m !== nothing
            j = i
            while !startswith(strip(lines[j]), "/"); j += 1; end
            push!(out, (String(m.captures[1]), i, j)); i = j
        end
        i += 1
    end
    out
end
grp(L, name) = first(filter(x -> x[1] == name, groups(L)))
insert_after(L, name, body) = (g = grp(L, name); vcat(L[1:g[3]], [""], body, L[g[3]+1:end]))
function setkey(L, name, key, val; note="")      # keeps the comment, or sets one
    g = grp(L, name); L = copy(L)
    for i in g[2]:g[3]
        m = match(Regex("^(\\s*" * key * "\\s*=\\s*)([^!]*?)(\\s*!.*)?\$"), L[i])
        if m !== nothing
            c = isempty(note) ? (m[3] === nothing ? "" : "   " * strip(m[3])) : "   ! " * note
            L[i] = m[1] * val * c; return L
        end
    end
    error("no key $key in &$name")
end
body(s) = split(s, "\n")

function sci(L)
    L = setkey(L, "ytill_ant", "scale_zb", "2"; note="2: exp with z0/z1 below (calib-cf; the base has 0: none)")
    L = setkey(L, "ytill_ant", "cf_ref", "0.2"; note="[-- or deg] Reference/const/max value of cf (calib-cf; the base has 0.1)")
    L = setkey(L, "ytill_nor", "cf_ref", "0.2"; note="[-- or deg] Reference/const/max value of cf (calib-cf; the base has 0.1)")
    L
end

function run(L, dom, case)
    ant = dom == "Antarctica"; sfx = ant ? "ant" : "nor"
    pd = case == "pd"; pd || case == "lgm" || error("case pd|lgm")
    L = setkey(L, "kryos", "domains", "\"$dom\"")
    L = setkey(L, "clock", "time_init", pd ? "0.0" : "30e3"; note="calib-cf: lgm restarts from the end of the pd run")
    L = setkey(L, "clock", "time_end", pd ? "30e3" : "60e3"; note="calib-cf: 30 kyr per case")
    L = setkey(L, "clock", "time_const", "0.0"; note="frozen forcing (the orbit is set in [orbit_$sfx])")
    L = setkey(L, "out_1D", "dt", "10.0"); L = setkey(L, "out_2D", "dt", "500.0")
    L = setkey(L, "out_2Dsm", "dt", "500.0"); L = setkey(L, "out_rst", "dt", "10e3")
    if pd
        L = setkey(L, "barysealevel", "method", "\"const\""; note="calib-cf: present-day sea level")
        L = setkey(L, "barysealevel", "bsl_init", "0.0")
    else
        L = setkey(L, "barysealevel", "method", "\"file\""; note="calib-cf: constant -120 m series")
        L = setkey(L, "barysealevel", "sl_path", "\"input/sl_const_lgm120.dat\"")
    end
    for k in ("f_min", "f_max"); L = setkey(L, "coef_alpha", k, pd ? "0.0" : "1.0"); end
    L = setkey(L, "orbit_$sfx", "method", "\"fixed_time\""); L = setkey(L, "orbit_$sfx", "time", pd ? "0.0" : "-21e3")
    if ant
        L = setkey(L, "Antarctica", "experiments", "\"\"")
        L = setkey(L, "climate_ant", "terms", "\"lgm\" \"ocean\"")
        for nm in ("snapshot_picontrol", "snapshot_lgm")
            L = setkey(L, nm, "tas", "\"pmip3_mean_tas_ann\" \"pmip3_mean_tas_djf\"")
        end
        L = insert_after(L, "pmip3_mean_tas_jja", body("""
&pmip3_mean_tas_djf
    filename = "fesmdata/Global/{grid_global}/{grid_global}_Paleoclimate-PMIP3-{exp}-mean.nc"
    name     = "t2m_djf"
    units    = "degC" "K"
    scaling  = 1.0 273.15
    time     = 0 0 0 0 1
    remap    = "bilin"
    grid_src = "{grid_global}"
/"""))
    else
        L = setkey(L, "North", "experiments", "\"\"")
    end
    L
end

function main(args)
    mode = get(args, 1, "")
    (mode == "sci" && length(args) == 3) || (mode == "run" && length(args) == 5) || (println(stderr, "usage: see the header"); return 2)
    L = readlines(args[2]); k = findfirst(startswith("&kryos"), L)
    hdr = L[1:k-1]; L = L[k:end]
    if mode == "sci"
        L = sci(L)
        open(args[3], "w") do io
            println(io, "! ", basename(args[3]), ": calib-cf variant of ", basename(args[2]), " generated by tools/bin/kryos-calib-par.jl sci. Changes: [ytill_ant] scale_zb 0 -> 2, cf_ref 0.1 -> 0.2; [ytill_nor] cf_ref 0.1 -> 0.2.")
            foreach(l -> println(io, l), hdr); foreach(l -> println(io, l), L)
        end
    else
        L = run(L, args[4], args[5])
        open(args[3], "w") do io
            println(io, "! kryos par file of the cf calibration (docs/dev/calib-cf.md): ", args[4], " only, case ", args[5], ", generated by tools/bin/kryos-calib-par.jl from ", basename(args[2]), ".")
            println(io)
            foreach(l -> println(io, l), L)
        end
    end
    return 0
end
exit(main(ARGS))
