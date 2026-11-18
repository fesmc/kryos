# The emulator of the joint ice + Earth ensemble (experiments/joint-ens.sh; docs/science/joint-ensemble.md): per data group a
# Gaussian process of the score over the 10 ice parameters of joint-ens.toml and the 5 parameters of the member's Earth, and the
# history-matching cut of a dense sample of the prior (plastic_hm.jl: the GP, its leave-one-out and the cut, as the Earth search).
#   julia --project=tools/levante tools/levante/plastic_joint_emulate.jl <design.csv> <joint-ens.toml> <earths.csv> <scores.csv> <config.toml> <outdir>
#         [npts=100000] [seed=2] [min=10]
# The cases: the references "tau-JE-<member> | E1D-ES1-<earth>" of scores.csv (member and Earth as in design.csv). The groups that
# judge and their cut as the scoring config says (plastic_hm.jl cut_setup: the bar's best I_g where the publications fail, the
# plain cut 3 for ESL). Inputs to [0, 1]: the ice parameters over their prior ranges (H_N in log), the Earth's over wave A's.
# The dense sample is the prior of the design: a Latin hypercube over the ice parameters, the Earth drawn uniformly from
# <earths.csv> (stratified). min: the fewest cases for a GP of a group (10; a smoke test of 3 members needs min=3, and its GPs mean
# nothing).
# Writes in outdir: summary.md, gp.csv, loo.csv, cases.csv (per member its parameters, I per group, log L, groups passed), sample.csv
# (the dense sample and its max I_emu over the judging groups), nroi_ice.png and nroi_earth.png (per parameter pair the fraction
# not ruled out).
include("plastic_hm.jl")
include("plastic_joint_ens.jl")          # ICE_PARAMS, EARTH_COLS, ice_params, to01, earth01, read_csv

const JOINT_RX = r"^tau-JE-(\S+) \| E1D-ES1-(E\d+)$"
const SET = "joint"

function joint_emulate(designf, cfgf, earthsf, scoresf, scfgf, outdir; npts = 100_000, seed = 2, cutI = 3.0, min_cases = 10)
    mkpath(outdir)
    P = ice_params(TOML.parsefile(cfgf)); pn = [ICE_PARAMS; EARTH_COLS]; ni = length(P); d = length(pn)
    D, _ = read_csv(designf)
    θ = Dict(str(r["member"]) => [[Float64(r[c]) for c in ICE_PARAMS]; [Float64(r[c]) for c in EARTH_COLS]] for r in D)
    earth = Dict(str(r["member"]) => str(r["earth"]) for r in D)
    std01(t) = [[to01(P[k], t[k]) for k in 1:ni]; earth01(t[ni+1:end])]
    E, _ = read_csv(earthsf); E01 = hcat([earth01([Float64(r[c]) for c in EARTH_COLS]) for r in E]...)
    sc = read_scores(scoresf); C = cut_setup(sc, TOML.parsefile(scfgf)); ug = C.ug
    cs = Tuple{String,String}[]
    for s in unique(first.(collect(keys(sc))))
        m = match(JOINT_RX, s); m === nothing && continue
        haskey(θ, m[1]) || error("$s: no member $(m[1]) in $designf")
        earth[m[1]] == m[2] || error("$s: member $(m[1]) has the Earth $(earth[m[1]]) in $designf")
        push!(cs, (String(m[1]), s))
    end
    isempty(cs) && error("no case \"tau-JE-<member> | E1D-ES1-<earth>\" in $scoresf")
    cases = Dict(SET => sort(cs)); println("cases: $(length(cs)) of $(length(D)) members")
    # the GPs and the cut on the dense sample: the ice parameters a Latin hypercube, the Earth stratified over the Earths
    rng = MersenneTwister(seed); U = lhs01(ni + 1, npts, rng)
    ie = clamp.(floor.(Int, U[end, :] .* size(E01, 2)) .+ 1, 1, size(E01, 2))
    Xs01 = vcat(U[1:ni, :], E01[:, ie])
    Imax = Dict(SET => fill(-Inf, npts)); Igrp = Dict((SET, g) => zeros(npts) for g in ug)
    gp, loorows, EMU = hm_fit(cases, θ, std01, sc, C, Xs01, Imax, Igrp; cutI = cutI, min_cases = min_cases)
    write_gp(joinpath(outdir, "gp.csv"), gp, pn); write_loo(joinpath(outdir, "loo.csv"), loorows)
    fitted = [g for g in ug if haskey(EMU, (SET, g))]
    keep = Imax[SET] .<= cutI
    from(k, u) = k <= ni ? from01(P[k], u) : EARTH1D_PARAMS[k-ni].lo + u * (EARTH1D_PARAMS[k-ni].hi - EARTH1D_PARAMS[k-ni].lo)
    Θs = [from(k, Xs01[k, i]) for k in 1:d, i in 1:npts]
    open(joinpath(outdir, "sample.csv"), "w") do io
        println(io, join(pn, ","), ",earth,Imax")
        for i in 1:npts
            println(io, join([@sprintf("%.4g", x) for x in Θs[:, i]], ","), ",", str(E[ie[i]]["earth"]), ",", @sprintf("%.2f", Imax[SET][i]))
        end
    end
    lo01 = zeros(d); hi01 = ones(d)
    if !isempty(fitted)
        plot_nroi(Xs01, keep, pn, lo01, hi01, 1:ni, joinpath(outdir, "nroi_ice.png"), "joint ensemble: not ruled out (ice parameters on [0, 1] of their ranges)")
        plot_nroi(Xs01, keep, pn, lo01, hi01, ni+1:d, joinpath(outdir, "nroi_earth.png"), "joint ensemble: not ruled out (Earth parameters on [0, 1] of wave A's ranges)")
    end
    # the cases
    LL(s) = sum(sc[(s, g)].logL for g in ug if haskey(sc, (s, g)); init = 0.0)
    np(s) = count(g -> haskey(sc, (s, g)) && sc[(s, g)].I <= C.cut[g], ug)
    fI(s, g) = haskey(sc, (s, g)) ? @sprintf("%.2f", sc[(s, g)].I) : "NaN"
    open(joinpath(outdir, "cases.csv"), "w") do io
        println(io, "member,", join(pn, ","), ",earth,", join(["I_$g" for g in ug], ","), ",logL,passed")
        for (m, s) in cases[SET]
            println(io, m, ",", join([@sprintf("%.6g", x) for x in θ[m]], ","), ",", earth[m], ",", join([fI(s, g) for g in ug], ","), ",",
                    @sprintf("%.1f", LL(s)), ",", np(s))
        end
    end
    esc(s) = replace(s, "|" => "\\|")
    refs = [s for s in unique(first.(collect(keys(sc)))) if match(JOINT_RX, s) === nothing]
    open(joinpath(outdir, "summary.md"), "w") do io
        println(io, "# The joint ice + Earth ensemble: scores and emulator\n\nCases: $(length(cs)) of $(length(D)) members ($designf; $scoresf). ",
                "Groups that judge: ", join([@sprintf("%s (cut %.2f)", g, C.cut[g]) for g in ug], ", "), ".\n")
        println(io, "## The emulators\n\nPer group a GP of y = log(chi2/n) in the ", d, " parameters (", join(pn, ", "), "; on [0, 1] of their ranges); ",
                "leave-one-out: RMSE and R² of y, the fraction of cases within 2 LOO sd, the cases the LOO prediction would rule out though they ",
                "pass (wrong) or keep though they fail (kept). Groups with fewer than $min_cases cases: no GP.\n")
        println(io, "| group | n | pass | y_cut | LOO RMSE | R² | within 2 sd | wrong | kept | length scales |\n|---|---|---|---|---|---|---|---|---|---|")
        for r in gp
            println(io, "| ", r.g, " | ", r.n, " | ", r.pass, " | ", @sprintf("%.3f | %.3f | %.3f | %.2f", r.ycut, r.rmse, r.r2, r.cov2), " | ",
                    r.wrong, " | ", r.missed, " | ", join([@sprintf("%.2f", x) for x in r.h.ℓ], ", "), " |")
        end
        println(io, "\n## Not ruled out\n\nOf $npts points of the prior (the ice parameters a Latin hypercube, the Earth one of ",
                size(E01, 2), " drawn uniformly), with the GPs of ", isempty(fitted) ? "no group" : join(fitted, ", "),
                @sprintf(": %.4f not ruled out (I_emu <= %.0f in every group).\n", mean(keep), cutI))
        if any(keep) && !isempty(fitted)
            println(io, "| parameter | prior | not ruled out: min-max (5-95 %) |\n|---|---|---|")
            for k in 1:d
                q = Θs[k, keep]
                println(io, "| ", pn[k], " | ", @sprintf("%.4g-%.4g", from(k, 0.0), from(k, 1.0)), " | ",
                        @sprintf("%.4g-%.4g (%.4g-%.4g)", minimum(q), maximum(q), quantile(q, 0.05), quantile(q, 0.95)), " |")
            end
            println(io, "\nRuled out by each group alone (fraction of the prior with I_emu > $cutI there): ",
                    join([@sprintf("%s %.3f", g, mean(Igrp[(SET, g)] .> cutI)) for g in fitted], ", "), ".\n")
        end
        # the far-field RMSE of the index points (sites.csv of the score folder, plastic_ensemble.jl write_sites), for context
        ffr = far_rmse(joinpath(dirname(scoresf), "sites.csv")); ffg = sort(unique(last.(collect(keys(ffr)))))
        rep = [g for g in unique(last.(collect(keys(sc)))) if !(g in ug)]
        println(io, "\n## The runs\n\nRanked by the groups passed, then log L over the groups that judge; the references first. ",
                isempty(rep) ? "" : "Reported, not judged: " * join(rep, ", ") * ". ",
                isempty(ffg) ? "" : "RMSE [m]: of the far-field index points (" * join(ffg, ", ") * "; sites.csv).", "\n")
        println(io, "| series | log L | passed | ", join(["I " * g for g in [ug; rep]], " | "), isempty(ffg) ? "" : " | " * join(["RMSE " * g for g in ffg], " | "),
                " |\n|---|---|---|", repeat("---|", length(ug) + length(rep) + length(ffg)))
        for s in [refs; [s for (_, s) in sort(cases[SET], by = x -> (-np(x[2]), -LL(x[2])))]]
            println(io, "| ", esc(s), " | ", @sprintf("%.1f", LL(s)), " | ", np(s), "/", length(ug), " | ", join([fI(s, g) for g in [ug; rep]], " | "),
                    isempty(ffg) ? "" : " | " * join([haskey(ffr, (s, g)) ? @sprintf("%.1f", ffr[(s, g)]) : "-" for g in ffg], " | "), " |")
        end
    end
    println("wrote ", joinpath(outdir, "summary.md"), @sprintf("; not ruled out %.4f", mean(keep)))
end

# the far-field index-point RMSE per (series, group) of a sites.csv (its "all" rows), empty without the file
function far_rmse(f)
    out = Dict{Tuple{String,String},Float64}()
    isfile(f) || return out
    M, h = readdlm(f, ',', Any; header = true); h = String.(vec(h)); c(n) = findfirst(==(n), h)
    for i in axes(M, 1)
        M[i, c("site")] == "all" || continue
        x = M[i, c("rmse_m")]; out[(strip(String(M[i, c("series")]), '"'), String(M[i, c("group")]))] = x isa Number ? Float64(x) : NaN
    end
    out
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 6 || error("plastic_joint_emulate.jl <design.csv> <joint-ens.toml> <earths.csv> <scores.csv> <config.toml> <outdir> [npts=] [seed=] [min=]")
    o = kv(ARGS[7:end])
    joint_emulate(ARGS[1:6]...; npts = parse(Int, get(o, "npts", "100000")), seed = parse(Int, get(o, "seed", "2")),
                  min_cases = parse(Int, get(o, "min", "10")))
end
