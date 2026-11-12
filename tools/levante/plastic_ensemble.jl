# A weighted ensemble of plastic ice histories (docs/dev/plastic-ensemble-design.md): a space-filling sample of the prior of
# the tau_b multipliers, each member a kryos run, and its weight from the GAPSLIP 2.1 index points with a correlated error per
# data group; members the data rule out (implausibility) get no weight. The jobs: plastic_ensemble.sh.
#   julia --project=tools/levante plastic_ensemble.jl init <config.toml> <ensemble folder>    wave 1: the members (iter_01/<m>)
#   julia --project=tools/levante plastic_ensemble.jl weights <ensemble folder> [<wave>]     score the members, the weights
#   julia --project=tools/levante plastic_ensemble.jl errmodel <ensemble folder> [<wave>] [<top>]   the error model from a wave
#   julia --project=tools/levante plastic_ensemble.jl earth-regions <config.toml> <out.nc>   the region masks of [earth]
#   julia --project=tools/levante plastic_ensemble.jl earth-shift <config.toml> <out.nc> UM=0.5 LM=0 <region>=1.0   one shift file
#   julia --project=tools/levante plastic_ensemble.jl earth-radial <3-D par.nml> <out par.nml>   the 1-D radial mean of the 3-D Earth
# The folder, the members and the scale groups are those of plastic_esmda.jl (prepare, write_iteration, member_series). Config:
# [scale], [prior] (sigma, ramp_ka, corr_time_kyr, theta_max), [run], [score] as plastic_esmda_v2.toml; [ensemble]: n (members
# besides the prior mean, m000), seed, obs_groups, compare_max_age_ka (the window of the judging and of the comparison with the
# references, 25), bar (the references whose best I_g per group sets the cut, weights), implausibility, use_groups;
# [[references]] (name and product | topo | run: the reconstructions the members are compared with, reference_series); [likelihood.<obs group>]: sigma_m [m], corr_dist_km, corr_time_kyr.
# [earth] (optional): the Earth as parameters, sampled with tau_b: shifts of log10(eta) of VILMA2's 3-D field, global in the upper
# and lower mantle and per region in the upper mantle (plastic_earth_shift.jl): sigma (prior sd [dex]), theta_max, boundary_km,
# regions_file (the masks, earth-regions), [[earth.regions]] (name, groups: scale groups), exclude_codes (REGIONS areas without a regional shift),
# key (vilma.visc_shift_file).
# Sample: a Latin hypercube in the standard normal space of the parameters (n strata per parameter, one point in each, random
# order), correlated by the Cholesky factor of the prior covariance, clipped to |theta| <= theta_max; m000 = the prior mean.
# Likelihood: per data group g the index points of the group (pm_compare.jl GROUPS), d = model - obs, R_g = diag(sigma_i^2) +
# sigma_m,g^2 C_g, C_ij = exp(-dist_ij / L_x,g - |t_i - t_j| / L_t,g), sigma_i the scorer's ip_sigma of the member (its slope and
# the side of d); log L = sum_g -1/2 (d' R_g^-1 d + log det R_g) (+ const). Groups are independent. Weights: w_m ∝ L_m (the
# sample follows the prior, so the prior and the sampling density cancel), normalised over the members that ran; the effective
# sample size ESS = 1 / sum w^2. The error model: [likelihood] of the config, or <folder>/errmodel.toml once errmodel has
# estimated it from a wave (maximum likelihood from the residuals of the members that fit each group best). Plausibility
# (history matching): the implausibility of a group I_g = (chi2_g - n_g) / sqrt(2 n_g), chi2_g = d' R_g^-1 d; a member is
# plausible when max_g I_g <= [ensemble] implausibility (3) over [ensemble] use_groups (default all); the weights (of those
# groups) are over the plausible members only. Writes
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
    P = ES.params(cfg, sg); pr = cfg["prior"]; nP = length(P)
    E = haskey(cfg, "earth") ? ES.earth_spec(cfg) : nothing; nE = E === nothing ? 0 : length(E.names)
    Cθ = ES.prior_cov(P, pr, cfg["scale"]["groups"])
    nE > 0 && (Cθ = [Cθ zeros(nP, nE); zeros(nE, nP) E.sigma^2 * I(nE)])     # the Earth shifts: independent, sd sigma
    Z = lhs_normal(nP + nE, en["n"], MersenneTwister(en["seed"]))
    Θ = cholesky(Symmetric(Cθ)).L * Z
    Θ[1:nP, :] .= clamp.(Θ[1:nP, :], -pr["theta_max"], pr["theta_max"])
    nE > 0 && (Θ[nP+1:end, :] .= clamp.(Θ[nP+1:end, :], -E.theta_max, E.theta_max))
    ids = collect(0:en["n"])
    ES.write_iteration(C, 1, cfg, sg, P, ids, hcat(zeros(nP + nE), Θ), ids)
    open(joinpath(C, "ensemble.md"), "w") do io
        println(io, "# Weighted ensemble of tau_b (plastic backend)\n\nconfig: ", cfgfile, "; ", ES.ngroups(sg), " groups x ", ES.nnodes(sg),
                " nodes, ", length(P), " parameters; prior sigma ", pr["sigma"], ", wave 1: ", en["n"], " members (Latin hypercube) and m000 ",
                "(the prior mean)\n\n```\n", join(ES.groups_summary(sg), "\n"), "\n```\n\nparameters (prior sd): ",
                join([@sprintf("%s (%.2f)", p.name, sqrt(Cθ[q, q])) for (q, p) in enumerate(P)], ", "),
                nE == 0 ? "" : "\n\nEarth (log10 viscosity shifts [dex], prior sd $(E.sigma), |shift| <= $(E.theta_max); " *
                "upper mantle above $(E.boundary_km) km, regions weighted by $(E.masks_file)): " * join(E.names, ", "))
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
        res[m] = (esl = r.esl, rg = group_residuals(obs, idx, r.S, A), S = r.S)
    end
    ids, Θ, res, og, obs, A
end

# --- the references --------------------------------------------------------------------------------------------------------------
# The reconstructions the members are compared with ([[references]] of the config): per entry a name and one of
#   product  a FesmData IceHistory file (under [score] fesmdata) with dz_sl, its own GIA (PaleoMIST: SELEN, Gowan's 1-D Earth)
#   topo     a FesmData IceHistory file with the topography relative to the sea level of the time (GLAC1D: its own GIA, PM.rsl_topo)
#   run      a standalone VILMA2 run under a reconstruction's ice: a gia_refcorr.sh case (its converged iteration, r.nc -> r_<n>.nc:
#            it<n>) or a folder with vilma_out.nc ($VAR in the path from the environment)
# The RSL at the index points (npoints x nages), NaN at the ages the reference does not have; nothing (and why) when it is missing.
expandenv(p) = replace(p, r"\$\{?(\w+)\}?" => m -> ENV[strip(m, ['$', '{', '}'])])
function reference_series(ref, fes, obs, A)
    plon = [o.lon for o in obs]; plat = [o.lat for o in obs]
    if haskey(ref, "run")
        d = expandenv(ref["run"]); rc = joinpath(d, "r.nc")
        islink(rc) && (d = joinpath(d, "it" * match(r"r_(\d+)\.nc", readlink(rc))[1]))
        r = PM.read_run(d); r === nothing && return nothing, "no run in $d"
        all(a -> any(x -> abs(x - PM.trun(r, a)) < 0.01, r.t), A) || return nothing, "$d: records missing"
        return PM.series_run(r, plon, plat, A; rec = a -> PM.trun(r, a)), ""
    end
    f = joinpath(fes, haskey(ref, "product") ? ref["product"] : ref["topo"]); isfile(f) || return nothing, "no $f"
    lon, lat, ages, sl = haskey(ref, "product") ? NCDataset(ds -> (Float64.(ds["lon"][:]), Float64.(ds["lat"][:]), Float64.(ds["time"][:]),
                                                                   Float64.(nomissing(ds["dz_sl"][:, :, :], NaN))), f) : PM.rsl_topo(f)
    S = fill(NaN, length(obs), length(A)); k = findall(in(ages), A)
    S[:, k] = PM.series_pm(lon, lat, ages, sl, plon, plat, A[k])
    S, ""
end
function references(cfg, obs, A)
    out = Pair{String,Matrix{Float64}}[]
    for ref in get(cfg, "references", [])
        S, why = reference_series(ref, cfg["score"]["fesmdata"], obs, A)
        S === nothing ? println("reference ", ref["name"], ": ", why) : push!(out, ref["name"] => S)
    end
    out
end

# A scoring window: the index points of age <= tmax [yr] and the records of age <= tmax (the references that end there: GLAC1D's own
# GIA at 26 ka, so tmax = 25 ka, its last record of the 2.5 kyr ages of the scores). The series of all, members and references, are
# scored on the same window.
function score_window(obs, A, og, tmax)
    kp = findall(o -> o.t <= tmax, obs); ka = findall(<=(tmax), A)
    (tmax = tmax, kp = kp, ka = ka, A = A[ka], obs = obs[kp], idx = [(g, findall(o -> o.group == g, obs[kp])) for g in og])
end
# per data group (log L, chi2, n, I) of a series S (all points x all ages) on the window W; nothing when S misses a record there (the
# points where S has gaps are dropped, as in group_residuals)
function score_series(S, W, em, og)
    Sw = S[W.kp, W.ka]; any(k -> all(isnan, view(Sw, :, k)), axes(Sw, 2)) && return nothing
    rg = group_residuals(W.obs, W.idx, Sw, W.A)
    Dict(g => (v = loglik_group(rg[g], em[g]...); (v..., implausibility(v[2], v[3]))) for g in og)
end

# The error model from a wave: per data group, sigma_m, L_x, L_t by maximum likelihood over a grid, from the residuals of the
# `top` members that fit that group best under the current error model (the members the data do not rule out; the rest would
# inflate the error by their own misfit). Writes <folder>/errmodel.toml (read by weights from then on) and errmodel.md.
const GRID_SIGMA = [3.0, 5.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 23.0, 26.0, 30.0, 35.0, 40.0, 50.0, 65.0, 80.0]   # [m]
const GRID_LX = [10.0, 25.0, 50.0, 100.0, 175.0, 250.0, 350.0, 500.0, 750.0, 1000.0, 1500.0, 2500.0]   # [km]
const GRID_LT = [1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 80.0, 160.0]                            # [kyr]
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

# The weights of wave k, on the judging window (score_window: the index points of age <= [ensemble] compare_max_age_ka, all when it
# is not set), per member and data group log L, chi2/n and the implausibility I_g under the error model (error_model).
# Plausibility against the published reconstructions ([ensemble] bar, names of [[references]]): the cut of group g is the larger
# of [ensemble] implausibility (3) and the lowest I_g of the bar's references, so a member is plausible when it is in every group
# that judges ([ensemble] use_groups, default all) at least as good as the best of the reconstructions whose misfit was accepted
# (without bar: the plain cut). Weights w ∝ L over the plausible members (0 for the rest), tempered: per group the error covariance
# is scaled by s_g = max(1, the lowest chi2/n of the members and the references), so that the best series fits its error model
# (chi2/n = 1) and log L differences are those the data resolve (log L_t = log L + chi2 (1 - 1/s) / 2 - n log(s) / 2).
# Writes iter_kk/weights.csv, weights.md and plausible.txt (the folders of the plausible members, one per line).
function weights(C, k)
    C = abspath(C); cfg = ES.config(C); ids, Θ, res, og, obs, A = wave_residuals(C, k, cfg); en = cfg["ensemble"]
    em, src = error_model(cfg, C); cut0 = get(en, "implausibility", 3.0)
    ug = String.(get(en, "use_groups", og))     # the groups that judge: the cut and the weights; the rest reported
    all(in(og), ug) || error("use_groups $(ug) not all in obs_groups $(og)")
    sg = ES.scale_groups(cfg["scale"]); P = ES.params(cfg, sg)
    W = score_window(obs, A, og, 1000.0 * get(en, "compare_max_age_ka", Inf))
    ok = sort(collect(keys(res)))
    sc = Dict(m => score_series(res[m].S, W, em, og) for m in ok)                  # g => (log L, chi2, n, I)
    refs = references(cfg, obs, A)
    rsc = [(nm, s) for (nm, s) in ((nm, score_series(S, W, em, og)) for (nm, S) in refs) if s !== nothing]
    # the bar: per group the lowest I_g of the bar's references, and which one
    bar = String.(get(en, "bar", String[])); nob = setdiff(bar, first.(rsc))
    isempty(nob) || error("bar references not scored on the window: ", join(nob, ", "))
    best = Dict(g => (isempty(bar) ? (-Inf, "") : minimum((s[g][4], nm) for (nm, s) in rsc if nm in bar)) for g in og)
    cut = Dict(g => max(cut0, best[g][1]) for g in og)
    plaus = [m for m in ok if all(sc[m][g][4] <= cut[g] for g in ug)]
    # tempering: the scale factor of each group's error covariance
    sf = Dict(g => max(1.0, minimum(s[g][2] / s[g][3] for s in [[sc[m] for m in ok]; last.(rsc)])) for g in og)
    LL = Dict(m => sum(sc[m][g][1] for g in ug) for m in ok)
    LLt = Dict(m => sum(sc[m][g][1] + 0.5 * sc[m][g][2] * (1 - 1 / sf[g]) - 0.5 * sc[m][g][3] * log(sf[g]) for g in ug) for m in ok)
    w = Dict(m => 0.0 for m in ok)
    if !isempty(plaus)
        mx = maximum(LLt[m] for m in plaus); for m in plaus; w[m] = exp(LLt[m] - mx); end
        t = sum(values(w)); for m in plaus; w[m] /= t; end
    end
    ess = isempty(plaus) ? 0.0 : 1 / sum(w[m]^2 for m in plaus)
    It = ES.itdir(C, k)
    I(m, g) = sc[m][g][4]
    open(joinpath(It, "weights.csv"), "w") do io
        println(io, "member,ok,", join(["logL_$g" for g in og], ","), ",logL,logL_tempered,", join(["chi2n_$g" for g in og], ","), ",",
                join(["I_$g" for g in og], ","), ",plausible,esl_20ka,weight")
        for m in ids
            haskey(res, m) || (println(io, m, ",0", repeat(",NaN", 3length(og) + 5)); continue)
            println(io, m, ",1,", join([@sprintf("%.2f", sc[m][g][1]) for g in og], ","), ",", @sprintf("%.2f,%.2f", LL[m], LLt[m]), ",",
                    join([@sprintf("%.3f", sc[m][g][2] / sc[m][g][3]) for g in og], ","), ",", join([@sprintf("%.2f", I(m, g)) for g in og], ","),
                    ",", m in plaus ? 1 : 0, ",", @sprintf("%.2f", res[m].esl), ",", @sprintf("%.6g", w[m]))
        end
    end
    open(joinpath(It, "plausible.txt"), "w") do io
        for m in plaus; println(io, ES.mdir(C, k, m)); end
    end
    th = Θ[:, [findfirst(==(m), ids) for m in ok]]; wv = [w[m] for m in ok]
    names = haskey(cfg, "earth") ? [ES.pnames(P); ES.earth_spec(cfg).names] : ES.pnames(P)
    μ = th * wv; sd = sqrt.(max.(((th .- μ) .^ 2) * wv, 0))
    eslp = [res[m].esl for m in plaus]
    open(joinpath(It, "weights.md"), "w") do io
        println(io, "Wave $k: $(length(ok)) of $(length(ids)) members ran; error model: $src; judged on the index points of age <= ",
                isfinite(W.tmax) ? @sprintf("%g ka", W.tmax / 1000) : "80 ka (all)", "; groups that judge (the cut and the weights): ",
                join(ug, ", "), length(ug) < length(og) ? "; reported only: " * join(setdiff(og, ug), ", ") : "", ".\n")
        println(io, "Cut per group (the larger of $cut0 and the best I_g of the bar: ", isempty(bar) ? "none" : join(bar, ", "),
                "); scale factor s_g of the tempered weights:\n")
        println(io, "| group | n | cut | set by | s_g |\n|---|---|---|---|---|")
        for g in og
            println(io, "| ", g, " | ", sc[ok[1]][g][3], " | ", @sprintf("%.2f", cut[g]), " | ", cut[g] > cut0 ? best[g][2] : "-",
                    " | ", @sprintf("%.2f", sf[g]), " |")
        end
        println(io, "\nPlausible (I_g <= cut in every group that judges): $(length(plaus)) members", isempty(plaus) ? "" :
                @sprintf(", ESL at 20 ka %.1f to %.1f m; effective sample size of their weights %.2f", minimum(eslp), maximum(eslp), ess), ".\n")
        println(io, "Ruled out by group (I_g > cut): ", join([@sprintf("%s %d", g, count(m -> I(m, g) > cut[g], ok)) for g in og], ", "), "\n")
        println(io, "| member | plausible | weight | log L | log L tempered | ", join(["I $g" for g in og], " | "), " | ESL 20 ka [m] |\n|---|---|---|---|---|",
                repeat("---|", length(og) + 1))
        for m in sort(ok, by = m -> (m in plaus ? 0 : 1, -LLt[m]))
            println(io, "| m", lpad(m, 3, "0"), " | ", m in plaus ? "yes" : "no", " | ", @sprintf("%.3g", w[m]), " | ", @sprintf("%.1f", LL[m]),
                    " | ", @sprintf("%.1f", LLt[m]), " | ", join([@sprintf("%.1f", I(m, g)) for g in og], " | "), " | ", @sprintf("%.1f", res[m].esl), " |")
        end
        if !isempty(plaus)
            println(io, "\nWeighted parameters of the plausible members (tau_b: log factor; Earth: dex; prior mean 0): mean (sd)\n")
            for (q, n) in enumerate(names); println(io, "- ", n, ": ", @sprintf("%+.2f (%.2f)", μ[q], sd[q])); end
        end
    end
    compare_references(It, cfg, res, og, obs, A, em, ug, cut, refs)
    println("wrote ", joinpath(It, "weights.md"), "; ", length(plaus), " plausible, ESS ", @sprintf("%.2f", ess))
end

# The members against the references (references), on the window of [ensemble] compare_max_age_ka (25 ka: GLAC1D's own GIA ends at
# 26 ka) and on all points (the references that cover them): per series the log L of the groups that judge, I_g per group and the
# groups passed (I_g <= cut[g], weights). The references, the 5 members of the highest log L and the members' median. Appends to
# weights.md, writes references.csv (every series).
function compare_references(It, cfg, res, og, obs, A, em, ug, cut, refs)
    tw = 1000.0 * get(cfg["ensemble"], "compare_max_age_ka", 25.0)
    wins = [(@sprintf("index points of age <= %g ka", tw / 1000), score_window(obs, A, og, tw)), ("all index points", score_window(obs, A, og, Inf))]
    LL(s) = sum(s[g][1] for g in ug); np(s) = count(g -> s[g][4] <= cut[g], ug)
    io = open(joinpath(It, "weights.md"), "a"); csv = open(joinpath(It, "references.csv"), "w")
    println(csv, "window,series,reference,", join(["n_$g" for g in og], ","), ",", join(["I_$g" for g in og], ","), ",logL,passed")
    for (wn, (lbl, W)) in enumerate(wins)
        rs = [(nm, s) for (nm, s) in ((nm, score_series(S, W, em, og)) for (nm, S) in refs) if s !== nothing]
        ms = Dict(m => score_series(res[m].S, W, em, og) for m in keys(res))
        for (nm, s, isref) in [[(nm, s, true) for (nm, s) in rs]; [("m" * lpad(m, 3, "0"), ms[m], false) for m in sort(collect(keys(ms)))]]
            println(csv, wn == 1 ? "window" : "all", ",", nm, ",", isref ? 1 : 0, ",", join([s[g][3] for g in og], ","), ",",
                    join([@sprintf("%.2f", s[g][4]) for g in og], ","), ",", @sprintf("%.2f", LL(s)), ",", np(s))
        end
        best = first(sort(collect(keys(ms)), by = m -> -LL(ms[m])), 5)
        println(io, "\n## Against the references: ", lbl, "\n\nlog L and groups passed over the groups that judge (", join(ug, ", "),
                "); I_g per group (passed at <= the cut of the group, above). References missing here: ",
                isempty(setdiff(first.(refs), first.(rs))) ? "none" : join(setdiff(first.(refs), first.(rs)), ", "), ".\n")
        println(io, "| series | log L | ", join(["I $g" for g in og], " | "), " | passed |\n|---|---|", repeat("---|", length(og) + 1))
        row(nm, s) = println(io, "| ", nm, " | ", @sprintf("%.1f", LL(s)), " | ", join([@sprintf("%.1f", s[g][4]) for g in og], " | "), " | ",
                             np(s), "/", length(ug), " |")
        for (nm, s) in rs; row("**" * nm * "**", s); end
        for m in best; row("m" * lpad(m, 3, "0"), ms[m]); end
        println(io, "| members' median | ", @sprintf("%.1f", median(LL(s) for s in values(ms))), " | ",
                join([@sprintf("%.1f", median(s[g][4] for s in values(ms))) for g in og], " | "), " | |")
        println(io, "| points (members) | | ", join([string(first(values(ms))[g][3]) for g in og], " | "), " | |")
        println(io, "\nMembers passing each group: ", join([@sprintf("%s %d", g, count(s -> s[g][4] <= cut[g], values(ms))) for g in og], ", "),
                "; all groups that judge: ", count(s -> np(s) == length(ug), values(ms)), " of ", length(ms), ".")
    end
    close(io); close(csv)
end

# the Earth parameters: the region masks of a config ([earth]), and a shift file of given values (a test, or a check of r)
function earth_regions(cfgfile, out)
    ES.region_masks(ES.load_config(cfgfile), out)
end
function earth_shift(cfgfile, out, kv)
    cfg = ES.load_config(cfgfile); E = ES.earth_spec(cfg); lon, lat, W = ES.read_masks(E.masks_file, E.regions)
    v = Dict(split(a, "=")[1] => parse(Float64, split(a, "=")[2]) for a in kv)
    unknown = setdiff(keys(v), ["UM"; "LM"; E.regions]); isempty(unknown) || error("earth-shift: no parameter $(unknown); UM, LM, $(E.regions)")
    θe = [get(v, "UM", 0.0); get(v, "LM", 0.0); [get(v, r, 0.0) for r in E.regions]]
    ES.write_earth(out, E, lon, lat, W, θe); println("wrote ", out, ": ", join([n * "=" * string(x) for (n, x) in zip(E.names, θe)], " "))
end

function main(args)
    args[1] == "init" && return init(args[2], args[3])
    args[1] == "earth-regions" && return earth_regions(args[2], args[3])
    args[1] == "earth-shift" && return earth_shift(args[2], args[3], args[4:end])
    args[1] == "earth-radial" && return ES.earth_radial(args[2], args[3])
    args[1] == "weights" && return weights(args[2], length(args) >= 3 ? parse(Int, args[3]) : 1)
    args[1] == "errmodel" && return errmodel(args[2], length(args) >= 3 ? parse(Int, args[3]) : 1;
                                             top = length(args) >= 4 ? parse(Int, args[4]) : 6)
    error("plastic_ensemble.jl init <config> <folder> | weights <folder> [<wave>] | errmodel <folder> [<wave>] [<top>] | " *
          "earth-regions <config> <out.nc> | earth-shift <config> <out.nc> UM=.. LM=.. <region>=.. | earth-radial <par> <out par>")
end
abspath(PROGRAM_FILE) == @__FILE__() && main(ARGS)
