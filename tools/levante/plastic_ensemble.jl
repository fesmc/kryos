# A weighted ensemble of plastic ice histories (docs/dev/plastic-ensemble-design.md): a space-filling sample of the prior of
# the tau_b multipliers, each member a kryos run, and its weight from the GAPSLIP 2.1 index points with a correlated error per
# data group; no member is discarded. The jobs: plastic_ensemble.sh.
#   julia --project=tools/levante plastic_ensemble.jl init <config.toml> <ensemble folder>    wave 1: the members (iter_01/mNNN)
#   julia --project=tools/levante plastic_ensemble.jl weights <ensemble folder> [<wave>]     score the members, the weights
# The folder, the members and the scale groups are those of plastic_esmda.jl (prepare, write_iteration, member_series). Config:
# [scale], [prior] (sigma, ramp_ka, corr_time_kyr, theta_max), [run], [score] as plastic_esmda_v2.toml; [ensemble]: n (members
# besides the prior mean, m000), seed, obs_groups; [likelihood.<obs group>]: sigma_m [m], corr_dist_km, corr_time_kyr.
# Sample: a Latin hypercube in the standard normal space of the parameters (n strata per parameter, one point in each, random
# order), correlated by the Cholesky factor of the prior covariance, clipped to |theta| <= theta_max; m000 = the prior mean.
# Likelihood: per data group g the index points of the group (pm_compare.jl GROUPS), d = model - obs, R_g = diag(sigma_i^2) +
# sigma_m,g^2 C_g, C_ij = exp(-dist_ij / L_x,g - |t_i - t_j| / L_t,g), sigma_i the scorer's ip_sigma of the member (its slope and
# the side of d); log L = sum_g -1/2 (d' R_g^-1 d + log det R_g) (+ const). Groups are independent. Weights: w_m ∝ L_m (the
# sample follows the prior, so the prior and the sampling density cancel), normalised over the members that ran; the effective
# sample size ESS = 1 / sum w^2. Writes iter_kk/weights.csv (member, ok, log L per group and total, chi2/n per group, ESL at 20 ka,
# weight) and weights.md (ESS overall and per group, the members by weight, the weighted parameter means and sd).
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
function loglik(cfg, obs, idx, S, A)
    lk = cfg["likelihood"]; out = Dict{String,Tuple{Float64,Float64,Int}}()
    for (g, ii) in idx
        p = lk[g]; σm = p["sigma_m"]; Lx = p["corr_dist_km"]; Lt = 1000.0 * p["corr_time_kyr"]
        sc = [PM.score(obs[i], A, S[i, :]) for i in ii]
        k = [j for j in eachindex(ii) if isfinite(sc[j].z) && sc[j].z != 0]
        d = [sc[j].d for j in k]; σ = [sc[j].d / sc[j].z for j in k]
        R = ES.obs_cov([obs[ii[j]] for j in k], σ, σm, Lx, Lt)
        F = cholesky(Symmetric(R)); x = F.L \ d
        out[g] = (-0.5 * (sum(abs2, x) + logdet(F)), sum(abs2, x) / length(k), length(k))
    end
    out
end

function weights(C, k)
    C = abspath(C); cfg = ES.config(C); en = cfg["ensemble"]; og = en["obs_groups"]
    obs = [o for o in PM.read_indicators(cfg["score"]["fesmdata"]) if o.type == 0 && o.group in og]
    idx = [(g, findall(o -> o.group == g, obs)) for g in og]
    A = sort(NCDataset(ds -> Float64.(ds["time"][:]),
                       joinpath(cfg["score"]["fesmdata"], "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc")))
    ids, Θ = ES.read_theta(C, k)
    sg = ES.scale_groups(cfg["scale"]); P = ES.params(cfg, sg)
    res = Dict{Int,Any}()
    for m in ids
        r, why = ES.member_series(ES.mdir(C, k, m), obs, A)
        r === nothing && (println("member $m: $why"); continue)
        res[m] = merge(r, (ll = loglik(cfg, obs, idx, r.S, A),))
        println("member $m: log L ", @sprintf("%.1f", sum(first(v) for v in values(res[m].ll))))
    end
    ok = sort(collect(keys(res)))
    LL = [sum(first(res[m].ll[g]) for g in og) for m in ok]
    w = exp.(LL .- maximum(LL)); w ./= sum(w)
    ess(w) = 1 / sum(w .^ 2)
    wg = Dict(g => (x = [first(res[m].ll[g]) for m in ok]; v = exp.(x .- maximum(x)); v ./ sum(v)) for g in og)
    It = ES.itdir(C, k)
    open(joinpath(It, "weights.csv"), "w") do io
        println(io, "member,ok,", join(["logL_$g" for g in og], ","), ",logL,", join(["chi2n_$g" for g in og], ","), ",esl_20ka,weight")
        for m in ids
            j = findfirst(==(m), ok)
            if j === nothing; println(io, m, ",0", repeat(",NaN", 2length(og) + 3)); continue; end
            r = res[m]
            println(io, m, ",1,", join([@sprintf("%.2f", r.ll[g][1]) for g in og], ","), ",", @sprintf("%.2f", LL[j]), ",",
                    join([@sprintf("%.3f", r.ll[g][2]) for g in og], ","), ",", @sprintf("%.2f", r.esl), ",", @sprintf("%.6g", w[j]))
        end
    end
    th = Θ[:, [findfirst(==(m), ids) for m in ok]]
    μ = th * w; sd = sqrt.(max.(((th .- μ) .^ 2) * w, 0))
    open(joinpath(It, "weights.md"), "w") do io
        println(io, "Wave $k: $(length(ok)) of $(length(ids)) members ran. Effective sample size ", @sprintf("%.2f", ess(w)),
                "; per data group (weights from that group alone): ", join([@sprintf("%s %.1f", g, ess(wg[g])) for g in og], ", "), "\n")
        println(io, "| member | weight | log L | ", join(["chi2/n $g" for g in og], " | "), " | ESL 20 ka [m] |\n|---|---|---|",
                repeat("---|", length(og) + 1))
        for j in sortperm(w, rev = true)
            m = ok[j]; r = res[m]
            println(io, "| m", lpad(m, 3, "0"), " | ", @sprintf("%.3g", w[j]), " | ", @sprintf("%.1f", LL[j]), " | ",
                    join([@sprintf("%.2f", r.ll[g][2]) for g in og], " | "), " | ", @sprintf("%.1f", r.esl), " |")
        end
        println(io, "\nWeighted parameters (log factor; prior mean 0): mean (sd)\n")
        for (q, p) in enumerate(P); println(io, "- ", p.name, ": ", @sprintf("%+.2f (%.2f)", μ[q], sd[q])); end
    end
    println("wrote ", joinpath(It, "weights.md"), "; ESS ", @sprintf("%.2f", ess(w)))
end

function main(args)
    args[1] == "init" && return init(args[2], args[3])
    args[1] == "weights" && return weights(args[2], length(args) >= 3 ? parse(Int, args[3]) : 1)
    error("plastic_ensemble.jl init <config> <folder> | weights <folder> [<wave>]")
end
abspath(PROGRAM_FILE) == @__FILE__() && main(ARGS)
