# The joint ice + Earth ensemble (experiments/joint-ens.sh, experiments/joint-ens.toml; docs/science/joint-ensemble.md): the
# Coulomb law's parameters, the LGM timing of the margin and the Earth sampled together, PlasticIceSheet online on VILMA2 as C7
# (docs/science/results.md). This tool: the Earths, the design and the members' files.
#
#   julia --project=tools/levante tools/levante/plastic_joint_ens.jl earths <emulate dir> <members.csv> <out earths.csv> [cut=3]
#       the Earths of the ensemble: wave A's 1-D Earths of the Earth search not ruled out under either ice history, by the
#       wave-A emulator (plastic_earth_emulate.jl; <emulate dir>: its loo.csv and gp.csv, the groups that judged there) in
#       leave-one-out at the Earth: max over the groups of (y_loo - y_cut) / sd_loo <= cut under PaleoMIST's or GLAC-1D's ice.
#       The Earth itself, its parameters (<members.csv>) and per ice history its max implausibility and the group of it.
#   julia ... plastic_joint_ens.jl design <joint-ens.toml> <earths.csv> <calib NH-32KM.nc> <calib ANT-16KM.nc> <out design.csv>
#         [n=<members>] [seed=<seed>] [points=<points.csv>]
#       the design: a maximin Latin hypercube over the ice parameters of <joint-ens.toml> ([[param]]), the Earths stratified
#       (each drawn n/nE times, give or take one); the maximin criterion is phi_p (p = 50; Morris and Mitchell 1995) of the
#       distances in the standardized parameters, the Earth by its 5 parameters (wave A's ranges), minimized by swaps within a
#       column (an LHS stays one, the Earths keep their counts) from several random starts. points=<csv> instead: the members of
#       the file (columns member, the parameters, earth). The c_f of each member (c_f* of the calibration surface,
#       plastic_cf_surface.jl, times exp of its deviation) are added. design.csv: member, the parameters, earth, its 5
#       parameters, c_f per region (cf_NAIS, cf_Greenland, cf_Fennoscandia, cf_BarentsKara, cf_Antarctica, cf_Patagonia).
#   julia ... plastic_joint_ens.jl files <design.csv> <template par.nml> <W> <S> [members=<m1,m2,...>]
#       per member <S>/members/<member>/cf.nc (taub_cf_file) and <S>/runs/par_<name>.nml (<S>: the set's folder, <W> the
#       ensemble's: its Earths, their r and the shared margin patch): the template (C7's par file) with
#       H_N, lambda_min (ice_nor and ice_pat: lam_N; ice_ant: lam_A), the c_f file, the Earth (&vilma of its par file of the
#       Earth search, <W>/earths/<earth>.nml) and its present-day correction r (<W>/refcorr/ice-PMIST_E1D-ES1-<earth>: the
#       reference topographies and the X background), and ice_nor's margin patches margin_na (region 1,
#       <S>/members/<member>/margin_na.nc), margin_eu (regions 3, 6, margin_eu.nc) and margin_grl (region 2,
#       <W>/margins/NH-32KM_IceHistory-PaleoGrIS-patch.nc); the patches are made by margin_patch.jl (joint-ens.sh). The
#       template must have the fixed settings of the design (the Coulomb law with taub_delta 0.05 and lambda from the file in
#       every instance, the region maps v2, ice_ant's zb -600 / 0) and no margin patches. Run name: tau-JE-<member>_E1D-ES1-<earth>
#       (the series "tau-JE-<member> | E1D-ES1-<earth>").
using DelimitedFiles, LinearAlgebra, Printf, Random, Statistics, TOML
include("plastic_cf_surface.jl")          # cf_surface, cf_star, cf_by_id, write_cf (plastic_cf_regions.jl)
include("plastic_earth_basis.jl")         # EARTH1D_PARAMS: the Earth's parameters and wave A's ranges

const ICE_PARAMS = ["H_N", "lam_N", "lam_A", "dlcf_NAIS", "dlcf_Greenland", "dlcf_Fennoscandia", "dlcf_BarentsKara",
                    "dlcf_Antarctica", "tau_NA", "tau_EU"]
const EARTH_COLS = [p.name for p in EARTH1D_PARAMS]
run_name(m, earth) = "tau-JE-$(m)_E1D-ES1-$(earth)"

kv(args) = Dict(String(a[1]) => String(a[2]) for a in (split(x, "=", limit = 2) for x in args))
function read_csv(f)
    M, h = readdlm(f, ',', Any; header = true); h = String.(vec(h))
    [Dict(h[j] => M[i, j] for j in eachindex(h)) for i in axes(M, 1)], h
end
str(x) = x isa AbstractString ? String(x) : string(x)

# --- the Earths ----------------------------------------------------------------------------------------------------------------
function earths(emu, membersf, out; cut = 3.0)
    gp, _ = read_csv(joinpath(emu, "gp.csv")); lo, _ = read_csv(joinpath(emu, "loo.csv"))
    ycut = Dict((str(r["ice"]), str(r["group"])) => Float64(r["ycut"]) for r in gp)
    I = Dict{Tuple{String,String},Vector{Tuple{Float64,String}}}()           # (ice, member) => [(I_loo, group)]
    for r in lo
        k = (str(r["ice"]), str(r["group"]))
        push!(get!(I, (k[1], str(r["member"])), Tuple{Float64,String}[]), ((Float64(r["y_loo"]) - ycut[k]) / Float64(r["sd_loo"]), k[2]))
    end
    ices = sort(unique(first.(collect(keys(I))))); groups = sort(unique(last.(collect(keys(ycut)))))
    mem, h = read_csv(membersf); h[2:end] == EARTH_COLS || error("$membersf: columns $(h[2:end]), not $EARTH_COLS")
    ng = Dict(i => length(unique([g for (ii, g) in keys(ycut) if ii == i])) for i in ices)
    rows = []
    for r in mem
        m = str(r["member"]); mx = Dict{String,Tuple{Float64,String}}()
        for ice in ices
            v = get(I, (ice, m), Tuple{Float64,String}[])
            length(v) == ng[ice] || error("$m under $ice: $(length(v)) of $(ng[ice]) groups in loo.csv")
            mx[ice] = maximum(v)
        end
        kept = [ice for ice in ices if mx[ice][1] <= cut]
        isempty(kept) || push!(rows, (m = m, θ = [Float64(r[c]) for c in EARTH_COLS], mx = mx, kept = kept))
    end
    open(out, "w") do io
        println(io, "earth,", join(EARTH_COLS, ","), ",", join(["Imax_$i,group_$i" for i in ices], ","), ",kept_under")
        for r in rows
            println(io, r.m, ",", join([@sprintf("%.6g", x) for x in r.θ], ","), ",",
                    join([@sprintf("%.2f,%s", r.mx[i]...) for i in ices], ","), ",", join(r.kept, " "))
        end
    end
    println("wrote $out: $(length(rows)) of $(length(mem)) Earths not ruled out under either ice history (",
            join(["$i $(count(r -> i in r.kept, rows))" for i in ices], ", "), ", both ", count(r -> length(r.kept) == length(ices), rows),
            "); groups: ", join(groups, ", "), "; leave-one-out of $emu, cut $cut")
    for (k, c) in enumerate(EARTH_COLS)
        q = [r.θ[k] for r in rows]
        @printf("  %-14s %8.2f - %8.2f (prior %g - %g)\n", c, minimum(q), maximum(q), EARTH1D_PARAMS[k].lo, EARTH1D_PARAMS[k].hi)
    end
    out
end

# --- the design ----------------------------------------------------------------------------------------------------------------
# the ice parameters of the config: name, lo, hi, scale; to [0, 1] and back
function ice_params(cfg)
    P = [(name = String(p["name"]), lo = Float64(p["lo"]), hi = Float64(p["hi"]), log = get(p, "scale", "lin") == "log")
         for p in cfg["param"]]
    [p.name for p in P] == ICE_PARAMS || error("[[param]]: $([p.name for p in P]), not $ICE_PARAMS (in this order)")
    all(p -> p.hi > p.lo && (!p.log || p.lo > 0), P) || error("[[param]]: hi > lo, and lo > 0 on a log scale")
    P
end
from01(p, u) = p.log ? exp(log(p.lo) + u * (log(p.hi) - log(p.lo))) : p.lo + u * (p.hi - p.lo)
to01(p, x) = p.log ? (log(x) - log(p.lo)) / (log(p.hi) - log(p.lo)) : (x - p.lo) / (p.hi - p.lo)
earth01(θ) = [(θ[k] - EARTH1D_PARAMS[k].lo) / (EARTH1D_PARAMS[k].hi - EARTH1D_PARAMS[k].lo) for k in eachindex(θ)]

# phi_p of the points (columns of X): (sum_{i<j} d_ij^-p)^(1/p), and the minimum distance
function phip(X, p)
    s = 0.0; dmin = Inf
    for i in axes(X, 2), j in i+1:size(X, 2)
        d = norm(view(X, :, i) .- view(X, :, j)); s += d^-p; dmin = min(dmin, d)
    end
    s^(1 / p), dmin
end
# the maximin LHS: U (d x n) the ice coordinates in [0, 1], e (n) the Earth of each member; X = [U; Earth coordinates]
function maximin(nice, n, E01, rng; p = 50.0, starts = 8, swaps = 40_000)
    nE = size(E01, 2); best = nothing
    for s in 1:starts
        U = permutedims(hcat([(randperm(rng, n) .- rand(rng, n)) ./ n for _ in 1:nice]...))
        o = randperm(rng, nE); e = shuffle(rng, [o[mod1(k, nE)] for k in 1:n])    # each Earth n / nE times, give or take one
        X = vcat(U, E01[:, e]); D = [i == j ? Inf : norm(X[:, i] .- X[:, j]) for i in 1:n, j in 1:n]
        S = sum(D[i, j]^-p for i in 1:n for j in i+1:n)
        for it in 1:swaps
            a, b = rand(rng, 1:n), rand(rng, 1:n); a == b && continue
            k = rand(rng, 0:nice)                                   # 0: swap the Earths of a and b, else ice column k
            k == 0 && e[a] == e[b] && continue
            Xa = copy(X[:, a]); Xb = copy(X[:, b])
            if k == 0; Xa[nice+1:end], Xb[nice+1:end] = X[nice+1:end, b], X[nice+1:end, a]
            else; Xa[k], Xb[k] = X[k, b], X[k, a]; end
            dS = 0.0; da = zeros(n); db = zeros(n)
            for j in 1:n
                (j == a || j == b) && continue
                da[j] = norm(Xa .- view(X, :, j)); db[j] = norm(Xb .- view(X, :, j))
                dS += da[j]^-p + db[j]^-p - D[a, j]^-p - D[b, j]^-p
            end
            if dS < 0
                X[:, a] = Xa; X[:, b] = Xb; S += dS
                for j in 1:n
                    (j == a || j == b) && continue
                    D[a, j] = D[j, a] = da[j]; D[b, j] = D[j, b] = db[j]
                end
                if k == 0; e[a], e[b] = e[b], e[a]; else; U[k, a], U[k, b] = U[k, b], U[k, a]; end
            end
        end
        φ, dmin = phip(X, p)
        @printf("  start %d: phi_p %.4f, min distance %.4f\n", s, φ, dmin)
        (best === nothing || φ < best.φ) && (best = (U = copy(U), e = copy(e), φ = φ, dmin = dmin))
    end
    best
end

function design(cfgf, earthsf, fN, fA, out; n = nothing, seed = nothing, points = nothing)
    cfg = TOML.parsefile(cfgf); P = ice_params(cfg); d = cfg["design"]
    n = n === nothing ? Int(d["n"]) : n; seed = seed === nothing ? Int(d["seed"]) : seed
    E, _ = read_csv(earthsf); enames = [str(r["earth"]) for r in E]; Eθ = Dict(str(r["earth"]) => [Float64(r[c]) for c in EARTH_COLS] for r in E)
    S = cf_surface(fN, fA)
    rows = []                                                       # (member, ice parameters, earth)
    if points === nothing
        rng = MersenneTwister(seed); E01 = hcat([earth01(Eθ[e]) for e in enames]...)
        φ0, d0 = phip(vcat(rand(rng, length(P), n), E01[:, rand(rng, 1:length(enames), n)]), 50.0)
        println("design: $n members, $(length(P)) ice parameters x $(length(enames)) Earths (random design: phi_p $(round(φ0, digits = 4)), ",
                "min distance $(round(d0, digits = 4)))")
        b = maximin(length(P), n, E01, rng; starts = get(d, "starts", 8), swaps = get(d, "swaps", 40_000))
        @printf("maximin: phi_p %.4f, min distance %.4f\n", b.φ, b.dmin)
        nd = ndigits(n)
        for i in 1:n
            push!(rows, ("m" * lpad(i, max(3, nd), '0'), [from01(P[k], b.U[k, i]) for k in eachindex(P)], enames[b.e[i]]))
        end
    else
        R, h = read_csv(points)
        issubset(["member"; ICE_PARAMS; "earth"], h) || error("$points: columns member, $(join(ICE_PARAMS, ", ")), earth")
        for r in R
            x = [Float64(r[c]) for c in ICE_PARAMS]; e = str(r["earth"])
            haskey(Eθ, e) || error("$points: $(r["member"]) has the Earth $e, not one of $earthsf")
            all(k -> P[k].lo - 1e-9 <= x[k] <= P[k].hi + 1e-9, eachindex(P)) || error("$points: $(r["member"]) outside the ranges of $cfgf")
            push!(rows, (str(r["member"]), x, e))
        end
    end
    open(out, "w") do io
        println(io, "member,", join(ICE_PARAMS, ","), ",earth,", join(EARTH_COLS, ","), ",",
                join(["cf_$r" for r in CF_REGIONS], ","), ",cf_Patagonia")
        for (m, x, e) in rows
            v = Dict(zip(ICE_PARAMS, x)); cs = cf_star(S, v["H_N"], v["lam_N"], v["lam_A"])
            cf = Dict(r => cs[r] * exp(v["dlcf_$r"]) for r in CF_REGIONS)
            println(io, m, ",", join([@sprintf("%.6g", y) for y in x], ","), ",", e, ",", join([@sprintf("%.6g", y) for y in Eθ[e]], ","), ",",
                    join([@sprintf("%.6g", cf[r]) for r in CF_REGIONS], ","), ",", @sprintf("%.6g", cs["Fennoscandia"]))
        end
    end
    cnt = Dict(e => count(r -> r[3] == e, rows) for e in enames)
    println("wrote $out: $(length(rows)) members; Earths drawn $(minimum(values(cnt)))-$(maximum(values(cnt))) times each")
    out
end

# --- the members' files ---------------------------------------------------------------------------------------------------------
# a namelist file as lines: set / read a key of a group (the group's name case-insensitive), append groups
function nml_group(lines, g)
    i = findfirst(l -> lowercase(strip(l)) == "&" * lowercase(g), lines); i === nothing && error("no &$g")
    j = findnext(l -> strip(l) == "/", lines, i + 1); (i, j)
end
keyrx(key) = Regex("^\\s*" * key * "\\s*=", "i")
function nml_get(lines, g, key)
    i, j = nml_group(lines, g)
    k = findfirst(l -> occursin(keyrx(key), l), lines[i+1:j-1]); k === nothing && return nothing
    v = strip(split(lines[i+k], "=", limit = 2)[2]); strip(replace(v, r"\s*!.*$" => ""))
end
function nml_put!(lines, g, key, val)
    i, j = nml_group(lines, g)
    for k in i+1:j-1
        occursin(keyrx(key), lines[k]) && (lines[k] = @sprintf(" %-15s = %s", key, val); return lines)
    end
    insert!(lines, j, @sprintf(" %-15s = %s", key, val)); lines
end
q(s) = "'$s'"

function files(designf, template, W, S; members = nothing)
    D, _ = read_csv(designf)
    members === nothing || (D = [r for r in D if str(r["member"]) in members]; length(D) == length(members) || error("members not in $designf"))
    T = readlines(template)
    # the fixed settings of the design (C7's): the template must have them
    for g in ("ice_nor", "ice_ant", "ice_pat")
        for (k, v) in [("taub_method", "'coulomb'"), ("taub_delta", "0.05"), ("taub_lambda_bed", "'file'"), ("taub_p", "1.0"),
                       ("taub_region_file", "'kryos_data/plastic/{grid_name}_taub_regions_v2.nc'"), ("taub_smooth_km", "0.0")]
            nml_get(T, g, k) == v || error("$template: &$g $k = $(nml_get(T, g, k)), not $v (the design's fixed setting)")
        end
        nml_get(T, g, "margin_patches") === nothing || error("$template: &$g has margin_patches already")
    end
    (nml_get(T, "ice_ant", "taub_zb0"), nml_get(T, "ice_ant", "taub_zb1")) == ("-600.0", "0.0") || error("$template: ice_ant's zb is not -600 / 0")
    (nml_get(T, "ice_nor", "taub_zb0"), nml_get(T, "ice_nor", "taub_zb1")) == ("-300.0", "200.0") || error("$template: ice_nor's zb is not -300 / 200")
    mkpath(joinpath(S, "runs")); pg = joinpath(W, "margins", "NH-32KM_IceHistory-PaleoGrIS-patch.nc")
    isfile(pg) || error("no $pg (joint-ens.sh margins)")
    for r in D
        m = str(r["member"]); e = str(r["earth"]); M = joinpath(S, "members", m); mkpath(M)
        rc = joinpath(W, "refcorr", "ice-PMIST_E1D-ES1-$e"); isfile(joinpath(rc, "r.nc")) || error("no $rc/r.nc (joint-ens.sh earths)")
        en = readlines(joinpath(W, "earths", "$e.nml"))
        cf = [Float64(r["cf_$k"]) for k in ("NAIS", "Greenland", "Fennoscandia", "Antarctica", "Patagonia", "BarentsKara")]
        write_cf(joinpath(M, "cf.nc"), cf)                         # the region ids of plastic_cf_regions.jl
        L = copy(T)
        for g in ("ice_nor", "ice_ant", "ice_pat")
            nml_put!(L, g, "taub_cf_file", q(joinpath(M, "cf.nc"))); nml_put!(L, g, "taub_HN", @sprintf("%.4f", r["H_N"]))
            nml_put!(L, g, "taub_lam_min", @sprintf("%.6f", g == "ice_ant" ? r["lam_A"] : r["lam_N"]))
        end
        for (g, f) in [("topo_nor", "North/{grid_name}/{grid_name}_TOPO-GEBCO2026-refcorr.nc"),
                       ("topo_ant", "Antarctica/{grid_name}/{grid_name}_TOPO-BedMachine-v4-refcorr.nc"),
                       ("topo_pat", "PatagoniaPaleo/{grid_name}/{grid_name}_TOPO-GEBCO2026-refcorr.nc")]
            nml_put!(L, g, "path", q(joinpath(rc, f)))
        end
        for g in ("X_z_bed", "X_H_ice"); nml_put!(L, g, "filename", q(joinpath(rc, "{grid_name}", "{grid_name}_REFCORR.nc"))); end
        for k in ("n_layer", "r_top", "r_bot", "eta", "rheology")
            v = nml_get(en, "vilma", k); v === nothing && error("$e.nml: no $k in &vilma"); nml_put!(L, "vilma", k, v)
        end
        nml_put!(L, "ice_nor", "margin_patches", "'margin_na' 'margin_eu' 'margin_grl'")
        for (g, f, ids) in [("margin_na", joinpath(M, "margin_na.nc"), "'1'"), ("margin_eu", joinpath(M, "margin_eu.nc"), "'3', '6'"),
                            ("margin_grl", pg, "'2'")]
            append!(L, ["", "&$g", " file        = '$f'", " var         = 'f_ice'", " H_min       = 0.5",
                        " region_file = 'kryos_data/plastic/{grid_name}_taub_regions_v2.nc'", " region_var  = 'region'",
                        " region_ids  = $ids", "/"])
        end
        n = run_name(m, e)
        open(joinpath(S, "runs", "par_$n.nml"), "w") do io; foreach(l -> println(io, l), L); end
        println("par_$n.nml")
    end
end

function joint_main(args)
    isempty(args) && error("plastic_joint_ens.jl earths | design | files (see the header)")
    o = kv([a for a in args if occursin("=", a)]); pos = [a for a in args if !occursin("=", a)]
    if pos[1] == "earths" && length(pos) == 4
        earths(pos[2], pos[3], pos[4]; cut = parse(Float64, get(o, "cut", "3")))
    elseif pos[1] == "design" && length(pos) == 6
        design(pos[2:6]...; n = haskey(o, "n") ? parse(Int, o["n"]) : nothing, seed = haskey(o, "seed") ? parse(Int, o["seed"]) : nothing,
               points = get(o, "points", nothing))
    elseif pos[1] == "files" && length(pos) == 5
        files(pos[2:5]...; members = haskey(o, "members") ? split(o["members"], ",") : nothing)
    else
        error("plastic_joint_ens.jl earths | design | files: wrong arguments $args (see the header)")
    end
end
abspath(PROGRAM_FILE) == @__FILE__() && joint_main(ARGS)
