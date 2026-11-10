# A weighted ensemble of plastic ice histories (docs/dev/plastic-ensemble-design.md): a space-filling sample of the prior of
# the tau_b multipliers, each member a kryos run, and its weight from the GAPSLIP 2.1 index points with a correlated error per
# data group; members the data rule out (implausibility) get no weight. The jobs: plastic_ensemble.sh.
#   julia --project=tools/levante plastic_ensemble.jl init <config.toml> <ensemble folder>    wave 1: the members (iter_01/mNNN)
#   julia --project=tools/levante plastic_ensemble.jl weights <ensemble folder> [<wave>]     score the members, the weights
#   julia --project=tools/levante plastic_ensemble.jl errmodel <ensemble folder> [<wave>] [<top>]   the error model from a wave
# The folder, the members and the scale groups are those of plastic_esmda.jl (prepare, write_iteration, member_series). Config:
# [scale], [prior] (sigma, ramp_ka, corr_time_kyr, theta_max), [run], [score] as plastic_esmda_v2.toml; [ensemble]: n (members
# besides the prior mean, m000), seed, obs_groups; [likelihood.<obs group>]: sigma_m [m], corr_dist_km, corr_time_kyr.
# Sample: a Latin hypercube in the standard normal space of the parameters (n strata per parameter, one point in each, random
# order), correlated by the Cholesky factor of the prior covariance, clipped to |theta| <= theta_max; m000 = the prior mean.
# Likelihood: per data group g the index points of the group (pm_compare.jl GROUPS), d = model - obs, R_g = diag(sigma_i^2) +
# sigma_m,g^2 C_g, C_ij = exp(-dist_ij / L_x,g - |t_i - t_j| / L_t,g), sigma_i the scorer's ip_sigma of the member (its slope and
# the side of d); log L = sum_g -1/2 (d' R_g^-1 d + log det R_g) (+ const). Groups are independent. Weights: w_m ∝ L_m (the
# sample follows the prior, so the prior and the sampling density cancel), normalised over the members that ran; the effective
# sample size ESS = 1 / sum w^2. The error model: [likelihood] of the config, or <folder>/errmodel.toml once errmodel has
# estimated it from a wave (maximum likelihood from the residuals of the members that fit each group best). Plausibility
# (history matching): the implausibility of a group I_g = (chi2_g - n_g) / sqrt(2 n_g), chi2_g = d' R_g^-1 d; a member is
# plausible when max_g I_g <= [ensemble] implausibility (3); the weights are over the plausible members only. Writes
# iter_kk/weights.csv (per member and group log L, chi2/n, I_g; plausible, ESL at 20 ka, weight), weights.md and plausible.txt.
using LinearAlgebra, Printf, Random, Statistics, TOML, NCDatasets

module ES
include("plastic_esmda.jl")
end
const PM = ES.PM

# inverse of the standard normal CDF (Acklam 2003; relative error < 1.2e-9)
function norminv(p)
    a = (-3.969683028665376e1, 2.209460984245205e2, -2.759285104469687e2, 1.383577518672690e2, -3.066479806614716e1, 2.506628277459239)
    b = (-5.447609879822406e1, 1.615858368580409e2, -1.556989798598866e2, 6.680131188771972e1, -1.328068155288572e1)
    c = (-7.784894002430293e-3, -3.223964580411365e-1, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783)
    d = (7.784695709041462e-3, 3.224671290700398e-1, 2.445134137142996, 3.754408661907416)
    pl = 0.02425
    if p < pl
        q = sqrt(-2log(p)); return (((((c[1]q + c[2])q + c[3])q + c[4])q + c[5])q + c[6]) / ((((d[1]q + d[2])q + d[3])q + d[4])q + 1)
    elseif p > 1 - pl
        q = sqrt(-2log(1 - p)); return -(((((c[1]q + c[2])q + c[3])q + c[4])q + c[5])q + c[6]) / ((((d[1]q + d[2])q + d[3])q + d[4])q + 1)
    end
    q = p - 0.5; r = q * q
    (((((a[1]r + a[2])r + a[3])r + a[4])r + a[5])r + a[6])q / (((((b[1]r + b[2])r + b[3])r + b[4])r + b[5])r + 1)
end
# Latin hypercube of n points in np standard normal dimensions
function lhs_normal(np, n, rng)
    Z = zeros(np, n)
    for i in 1:np
        u = (randperm(rng, n) .- rand(rng, n)) ./ n
        Z[i, :] = norminv.(clamp.(u, 1e-12, 1 - 1e-12))
    end
    Z
end

function init(cfgfile, C)
    C = abspath(C)
    en = TOML.parsefile(cfgfile)["ensemble"]
    cfg, sg, obs = ES.prepare(cfgfile, C, en["obs_groups"])
    P = ES.params(cfg, sg); pr = cfg["prior"]
    Cθ = ES.prior_cov(P, pr, cfg["scale"]["groups"])
    Z = lhs_normal(length(P), en["n"], MersenneTwister(en["seed"]))
    Θ = clamp.(cholesky(Symmetric(Cθ)).L * Z, -pr["theta_max"], pr["theta_max"])
    ids = collect(0:en["n"])
    ES.write_iteration(C, 1, cfg, sg, P, ids, hcat(zeros(length(P)), Θ), ids)
    open(joinpath(C, "ensemble.md"), "w") do io
        println(io, "# Weighted ensemble of tau_b (plastic backend)\n\nconfig: ", cfgfile, "; ", ES.ngroups(sg), " groups x ", ES.nnodes(sg),
                " nodes, ", length(P), " parameters; prior sigma ", pr["sigma"], ", wave 1: ", en["n"], " members (Latin hypercube) and m000 ",
                "(the prior mean)\n\n```\n", join(ES.groups_summary(sg), "\n"), "\n```\n\nparameters (prior sd): ",
                join([@sprintf("%s (%.2f)", p.name, sqrt(Cθ[q, q])) for (q, p) in enumerate(P)], ", "))
    end
    println("init: ", C, ", wave 1 with members 0-", en["n"])
end

# --- likelihood ---------------------------------------------------------------------------------------------------------------------
# the misfit of a member per data group: the index points with a finite score, d = model - obs and the scorer's 1 sigma
function group_residuals(obs, idx, S, A)
    out = Dict{String,Any}()
    for (g, ii) in idx
        sc = [PM.score(obs[i], A, S[i, :]) for i in ii]
        k = [j for j in eachindex(ii) if isfinite(sc[j].z) && sc[j].z != 0]
        out[g] = (pts = [obs[ii[j]] for j in k], d = [sc[j].d for j in k], σ = [sc[j].d / sc[j].z for j in k])
    end
    out
end
# log L (+ const), chi2 = d' R^-1 d and n of one group's misfit r under the model error (sigma_m [m], L_x [km], L_t [yr])
function loglik_group(r, σm, Lx, Lt; DT = nothing)
    R = DT === nothing ? ES.obs_cov(r.pts, r.σ, σm, Lx, Lt) : σm^2 .* exp.(.-DT.D ./ Lx .- DT.T ./ Lt) .+ Diagonal(r.σ .^ 2)
    F = cholesky(Symmetric(R)); x = F.L \ r.d; χ2 = sum(abs2, x)
    (-0.5 * (χ2 + logdet(F)), χ2, length(r.d))
end
# the error model per data group: <ensemble folder>/errmodel.toml (plastic_ensemble.jl errmodel) when it exists, else the config's
function error_model(cfg, C)
    f = joinpath(C, "errmodel.toml")
    lk = isfile(f) ? TOML.parsefile(f)["likelihood"] : cfg["likelihood"]
    Dict(g => (p["sigma_m"], p["corr_dist_km"], 1000.0 * p["corr_time_kyr"]) for (g, p) in lk), (isfile(f) ? f : "config")
end
# the implausibility of a group's misfit: (chi2 - n) / sqrt(2 n), the standardized chi2 (N(0, 1) for large n when R is right)
implausibility(χ2, n) = (χ2 - n) / sqrt(2n)

# the members of wave k that ran, their residuals per group and ESL at 20 ka
function wave_residuals(C, k, cfg)
    og = cfg["ensemble"]["obs_groups"]
    obs = [o for o in PM.read_indicators(cfg["score"]["fesmdata"]) if o.type == 0 && o.group in og]
    idx = [(g, findall(o -> o.group == g, obs)) for g in og]
    A = sort(NCDataset(ds -> Float64.(ds["time"][:]),
                       joinpath(cfg["score"]["fesmdata"], "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc")))
    ids, Θ = ES.read_theta(C, k)
    res = Dict{Int,Any}()
    for m in ids
        r, why = ES.member_series(ES.mdir(C, k, m), obs, A)
        r === nothing && (println("member $m: $why"); continue)
        res[m] = (esl = r.esl, rg = group_residuals(obs, idx, r.S, A))
    end
    ids, Θ, res, og
end

# The error model from a wave: per data group, sigma_m, L_x, L_t by maximum likelihood over a grid, from the residuals of the
# `top` members that fit that group best under the current error model (the members the data do not rule out; the rest would
# inflate the error by their own misfit). Writes <folder>/errmodel.toml (read by weights from then on) and errmodel.md.
const GRID_SIGMA = [3.0, 5.0, 8.0, 12.0, 16.0, 20.0, 25.0, 30.0, 40.0, 50.0, 65.0, 80.0]   # [m]
const GRID_LX = [25.0, 50.0, 100.0, 175.0, 250.0, 350.0, 500.0, 750.0, 1000.0]             # [km]
const GRID_LT = [1.0, 2.0, 5.0, 10.0, 20.0]                                                # [kyr]
function errmodel(C, k; top = 6)
    C = abspath(C); cfg = ES.config(C); ids, _, res, og = wave_residuals(C, k, cfg)
    em, src = error_model(cfg, C); ok = sort(collect(keys(res)))
    lk = Dict{String,Any}(); rows = []
    for g in og
        ll0 = Dict(m => loglik_group(res[m].rg[g], em[g]...) for m in ok)
        best = sort(ok, by = m -> -ll0[m][1])[1:min(top, length(ok))]
        # the distances [km] and age differences [yr] of each member's points, once
        DT = Dict(m => (p = res[m].rg[g].pts; (D = [ES.gcdist(a, b) for a in p, b in p], T = [abs(a.t - b.t) for a in p, b in p]))
                  for m in best)
        fit = (-Inf, 0.0, 0.0, 0.0)
        for σm in GRID_SIGMA, Lx in GRID_LX, Lt in GRID_LT
            L = sum(loglik_group(res[m].rg[g], σm, Lx, 1000Lt; DT = DT[m])[1] for m in best)
            L > fit[1] && (fit = (L, σm, Lx, Lt))
        end
        _, σm, Lx, Lt = fit
        edge = σm in (first(GRID_SIGMA), last(GRID_SIGMA)) || Lx in (first(GRID_LX), last(GRID_LX)) || Lt in (first(GRID_LT), last(GRID_LT))
        χn(m, p...) = (v = loglik_group(res[m].rg[g], p...); v[2] / v[3])
        lk[g] = Dict("sigma_m" => σm, "corr_dist_km" => Lx, "corr_time_kyr" => Lt)
        println(g, ": sigma_m ", σm, " m, L_x ", Lx, " km, L_t ", Lt, " kyr (members ", best, ")")
        push!(rows, (g = g, old = em[g], new = (σm, Lx, Lt), edge = edge, n = res[best[1]].rg[g].d |> length, best = best,
                     χold = mean(χn(m, em[g]...) for m in best), χnew = mean(χn(m, σm, Lx, 1000Lt) for m in best)))
    end
    open(joinpath(C, "errmodel.toml"), "w") do io
        println(io, "# the model error per data group from wave $k (plastic_ensemble.jl errmodel, top $top members per group, ",
                "maximum likelihood on a grid; previous: $src)")
        TOML.print(io, Dict("likelihood" => lk))
    end
    open(joinpath(C, "errmodel.md"), "w") do io
        println(io, "# Error model from wave $k

Per data group: sigma_m, L_x, L_t by maximum likelihood on a grid, from the residuals ",
                "of the $top members that fit the group best under the previous error model ($src). chi2/n: the mean over those members ",
                "(1 when the error model matches their misfit).
")
        println(io, "| group | n | sigma_m [m] | L_x [km] | L_t [kyr] | chi2/n | members |
|---|---|---|---|---|---|---|")
        for r in rows
            println(io, "| ", r.g, " | ", r.n, " | ", @sprintf("%g → %g", r.old[1], r.new[1]), " | ", @sprintf("%g → %g", r.old[2], r.new[2]),
                    " | ", @sprintf("%g → %g", r.old[3] / 1000, r.new[3]), " | ", @sprintf("%.2f → %.2f", r.χold, r.χnew), " | ",
                    join(["m" * lpad(m, 3, "0") for m in r.best], " "), r.edge ? " (on the edge of the grid)" : "", " |")
        end
    end
    println("wrote ", joinpath(C, "errmodel.toml"), ", errmodel.md")
end

# The weights of wave k: per member and data group log L, chi2/n and the implausibility I_g under the error model (error_model);
# a member is plausible when max_g I_g <= [ensemble] implausibility (default 3). Weights w ∝ L over the plausible members (0
# for the rest). Writes iter_kk/weights.csv, weights.md and plausible.txt (the folders of the plausible members, one per line).
function weights(C, k)
    C = abspath(C); cfg = ES.config(C); ids, Θ, res, og = wave_residuals(C, k, cfg)
    em, src = error_model(cfg, C); cut = get(cfg["ensemble"], "implausibility", 3.0)
    sg = ES.scale_groups(cfg["scale"]); P = ES.params(cfg, sg)
    ll = Dict(m => Dict(g => loglik_group(res[m].rg[g], em[g]...) for g in og) for m in keys(res))
    I = Dict(m => Dict(g => implausibility(ll[m][g][2], ll[m][g][3]) for g in og) for m in keys(res))
    ok = sort(collect(keys(res)))
    plaus = [m for m in ok if maximum(values(I[m])) <= cut]
    LL = Dict(m => sum(ll[m][g][1] for g in og) for m in ok)
    w = Dict(m => 0.0 for m in ok)
    if !isempty(plaus)
        mx = maximum(LL[m] for m in plaus); for m in plaus; w[m] = exp(LL[m] - mx); end
        t = sum(values(w)); for m in plaus; w[m] /= t; end
    end
    ess = isempty(plaus) ? 0.0 : 1 / sum(w[m]^2 for m in plaus)
    It = ES.itdir(C, k)
    open(joinpath(It, "weights.csv"), "w") do io
        println(io, "member,ok,", join(["logL_$g" for g in og], ","), ",logL,", join(["chi2n_$g" for g in og], ","), ",",
                join(["I_$g" for g in og], ","), ",plausible,esl_20ka,weight")
        for m in ids
            haskey(res, m) || (println(io, m, ",0", repeat(",NaN", 3length(og) + 4)); continue)
            println(io, m, ",1,", join([@sprintf("%.2f", ll[m][g][1]) for g in og], ","), ",", @sprintf("%.2f", LL[m]), ",",
                    join([@sprintf("%.3f", ll[m][g][2] / ll[m][g][3]) for g in og], ","), ",", join([@sprintf("%.2f", I[m][g]) for g in og], ","),
                    ",", m in plaus ? 1 : 0, ",", @sprintf("%.2f", res[m].esl), ",", @sprintf("%.6g", w[m]))
        end
    end
    open(joinpath(It, "plausible.txt"), "w") do io
        for m in plaus; println(io, ES.mdir(C, k, m)); end
    end
    th = Θ[:, [findfirst(==(m), ids) for m in ok]]; wv = [w[m] for m in ok]
    μ = th * wv; sd = sqrt.(max.(((th .- μ) .^ 2) * wv, 0))
    eslp = [res[m].esl for m in plaus]
    open(joinpath(It, "weights.md"), "w") do io
        println(io, "Wave $k: $(length(ok)) of $(length(ids)) members ran; error model: $src; implausibility cut $cut.\n")
        println(io, "Plausible (max_g I_g <= $cut): $(length(plaus)) members", isempty(plaus) ? "" :
                @sprintf(", ESL at 20 ka %.1f to %.1f m; effective sample size of their weights %.2f", minimum(eslp), maximum(eslp), ess), ".\n")
        println(io, "Ruled out by group (I_g > $cut): ", join([@sprintf("%s %d", g, count(m -> I[m][g] > cut, ok)) for g in og], ", "), "\n")
        println(io, "| member | plausible | weight | log L | ", join(["I $g" for g in og], " | "), " | ESL 20 ka [m] |\n|---|---|---|---|",
                repeat("---|", length(og) + 1))
        for m in sort(ok, by = m -> (m in plaus ? 0 : 1, -LL[m]))
            println(io, "| m", lpad(m, 3, "0"), " | ", m in plaus ? "yes" : "no", " | ", @sprintf("%.3g", w[m]), " | ", @sprintf("%.1f", LL[m]),
                    " | ", join([@sprintf("%.1f", I[m][g]) for g in og], " | "), " | ", @sprintf("%.1f", res[m].esl), " |")
        end
        if !isempty(plaus)
            println(io, "\nWeighted parameters of the plausible members (log factor; prior mean 0): mean (sd)\n")
            for (q, p) in enumerate(P); println(io, "- ", p.name, ": ", @sprintf("%+.2f (%.2f)", μ[q], sd[q])); end
        end
    end
    println("wrote ", joinpath(It, "weights.md"), "; ", length(plaus), " plausible, ESS ", @sprintf("%.2f", ess))
end

function main(args)
    args[1] == "init" && return init(args[2], args[3])
    args[1] == "weights" && return weights(args[2], length(args) >= 3 ? parse(Int, args[3]) : 1)
    args[1] == "errmodel" && return errmodel(args[2], length(args) >= 3 ? parse(Int, args[3]) : 1;
                                             top = length(args) >= 4 ? parse(Int, args[4]) : 6)
    error("plastic_ensemble.jl init <config> <folder> | weights <folder> [<wave>] | errmodel <folder> [<wave>] [<top>]")
end
abspath(PROGRAM_FILE) == @__FILE__() && main(ARGS)
