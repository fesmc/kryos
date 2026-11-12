# A weighted ensemble of plastic ice histories (docs/dev/plastic-ensemble-design.md): a space-filling sample of the prior of
# the tau_b multipliers, each member a kryos run, and its weight from the GAPSLIP 2.1 index points with a correlated error per
# data group; members the data rule out (implausibility) get no weight. The jobs: plastic_ensemble.sh.
#   julia --project=tools/levante plastic_ensemble.jl init <config.toml> <ensemble folder>    wave 1: the members (iter_01/<m>)
#   julia --project=tools/levante plastic_ensemble.jl weights <ensemble folder> [<wave>]     score the members, the weights
#   julia --project=tools/levante plastic_ensemble.jl errmodel <ensemble folder> [<wave>] [<top>]   the error model from a wave
#   julia --project=tools/levante plastic_ensemble.jl scores <ensemble folder> [<wave>]     RMSE, RMS z, I_g, log L of every series
#   julia --project=tools/levante plastic_ensemble.jl earth-regions <config.toml> <out.nc>   the region masks of [earth]
#   julia --project=tools/levante plastic_ensemble.jl earth-shift <config.toml> <out.nc> UM=0.5 LM=0 <region>=1.0   one shift file
#   julia --project=tools/levante plastic_ensemble.jl earth-basis <out.nc> <pan visc.nc> <pan litho.nc> | earth-field <basis.nc> <out.nc> c1 ...   (plastic_earth_basis.jl)
#   julia --project=tools/levante plastic_ensemble.jl earth-radial <3-D par.nml> <out par.nml>   the 1-D radial mean of the 3-D Earth
#   julia --project=tools/levante plastic_ensemble.jl earth-search-sample <1-D par.nml> <outdir> [n=200] [seed=1]   wave A of the
#       Earth search: 1-D Earths of the 5 radial parameters (plastic_earth_basis.jl; docs/dev/earth-search-design.md)
#   julia --project=tools/levante plastic_ensemble.jl earth-search-par <1-D par.nml> <out par.nml> L um tz lm1 lm2   one of them
# The folder, the members and the scale groups are those of plastic_esmda.jl (prepare, write_iteration, member_series). Config:
# [scale], [prior] (sigma, ramp_ka, corr_time_kyr, theta_max), [run], [score] as plastic_esmda_v2.toml; [ensemble]: n (members
# besides the prior mean, m000), seed, obs_groups, compare_max_age_ka (the window of the judging and of the comparison with the
# references, 25), bar (the references whose best I_g per group sets the cut, weights), implausibility, use_groups;
# [[references]] (name and product | topo | run: the reconstructions the members are compared with, reference_series); [likelihood.<obs group>]: sigma_m [m], corr_dist_km, corr_time_kyr.
# [esl] (optional): the ESL group, global sea level relative to the present against an independent curve: name ("ESL"), product
# (FesmData Records/SeaLevel/SeaLevel-<product>.nc, e.g. Spratt2016), var (sl_short), ages_ka (25 to 10 every 2.5), judge (true),
# reference_curves (SeaLevel products shown in weights.md, not scored, e.g. Lambeck2014); a reference's esl = "<product>" gives a
# published reconstruction its ESL (pub-PMIST: PaleoMIST). Model ESL: bsl_ive of the run's records (pm_compare.jl esl_at).
# [[gnss]] (optional): GNSS data groups, the present-day uplift rate at GNSS stations: name (the group), product
# (FesmData Records/GPS/GPS-<product>.nc under [score] gnss_dir, default <fesmdata>/Records/GPS), rate (gia | corrected | observed),
# box, exclude, min_years (pm_compare.jl read_gnss), judge (default true; false: reported only, not in the default use_groups);
# several entries of one name: one group of their stations, each once (gnss_data); every [[gnss]] group is a data group (after
# obs_groups: data_groups); [likelihood."<name>"]: sigma_m [mm/yr], corr_dist_km (no time term).
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
# plausible when max_g I_g <= [ensemble] implausibility (3) over [ensemble] use_groups (default all but [[gnss]] judge = false); the weights (of those
# groups) are over the plausible members only. Writes
# iter_kk/weights.csv (per member and group log L, chi2/n, I_g; plausible, ESL at 20 ka, weight), weights.md and plausible.txt.
# GNSS groups ([[gnss]]): d = model - obs at the stations of the group, the model dz_bed_rate of the run's last VILMA2 record [m/yr]
# (pm_compare.jl read_rate: the rate over VILMA2's last update, 100 yr in a kryos run) bilinear at the station, in mm/yr; R as above
# with sigma_i the station's sd and C_ij = exp(-dist_ij / L_x) (no time term). A series (member or reference) without dz_bed_rate (a
# published reconstruction, a run of VILMA2 before it) has no GNSS score: it is left out of the GNSS groups, the tables show "-",
# and a GNSS group judges only when every member that ran has a score in it. Frame: the model rate is VILMA2's (deg1_frame, centre of
# mass), the data ITRF/IGS (Schumacher 2018: CM): they differ by a few tenths of a mm/yr, not corrected.
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
    cfg, sg, obs = ES.prepare(cfgfile, C, [g for g in en["obs_groups"] if g in RSL_GROUPS])
    gd = gnss_data(cfg)
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
    isempty(gd.idx) || write_gnss(C, gd, cfg)
    println("init: ", C, ", wave 1 with members 0-", en["n"])
end

# --- GNSS -------------------------------------------------------------------------------------------------------------------------
# The data groups are of two kinds: RSL (the index points of pm_compare.jl's GROUPS) and GNSS ([[gnss]] of the config, by name).
const RSL_GROUPS = [first.(PM.GROUPS); collect(keys(PM.GROUP_UNIONS))]   # the RSL groups, and their unions (older configs)
# The data groups of a config: [ensemble] obs_groups and the GNSS groups of [[gnss]] (in the order of their entries); and those that
# judge by default (use_groups absent): all but the GNSS groups of judge = false.
esl_name(cfg) = haskey(cfg, "esl") ? String(get(cfg["esl"], "name", "ESL")) : nothing
data_groups(cfg) = unique([String.(cfg["ensemble"]["obs_groups"]); [String(sp["name"]) for sp in get(cfg, "gnss", Any[])];
                           esl_name(cfg) === nothing ? String[] : [esl_name(cfg)]])
judging_groups(cfg) = [g for g in data_groups(cfg) if !any(sp -> sp["name"] == g && !get(sp, "judge", true), get(cfg, "gnss", Any[])) &&
                                                     !(g == esl_name(cfg) && !get(cfg["esl"], "judge", true))]
# the kind of a data group: :rsl (index points), :gnss ([[gnss]]) or :esl ([esl])
group_kind(cfg, g) = g in RSL_GROUPS ? :rsl : g == esl_name(cfg) ? :esl : :gnss
# The point data of a config, the data groups that are not RSL index points: the stations of its GNSS groups (pm_compare.jl
# read_gnss) and the ages of its ESL group (read_esl), as points (group, kind, site_name, lat, lon, t, value, sd), and per group their
# indices (idx); every group of obs_groups must be an RSL group, a [[gnss]] entry or the [esl] group. A group of several entries (products) takes their
# stations in the order of the entries, without those of a later entry that an earlier one has (the same site name, or within
# MERGE_KM): a station counts once.
const MERGE_KM = 1.0
function gnss_data(cfg)
    og = data_groups(cfg); specs = get(cfg, "gnss", Any[])
    bad = [g for g in og if !(g in RSL_GROUPS) && !any(sp -> sp["name"] == g, specs) && g != esl_name(cfg)]
    isempty(bad) || error("obs_groups: neither an RSL group (pm_compare.jl GROUPS) nor a [[gnss]] entry nor [esl]: ", join(bad, ", "))
    dir = get(cfg["score"], "gnss_dir", joinpath(cfg["score"]["fesmdata"], "Records", "GPS"))
    sts = NamedTuple[]
    for sp in specs
        sp["name"] in og || continue
        have = [s for s in sts if s.group == sp["name"]]
        append!(sts, [st for st in PM.read_gnss(dir, sp)
                      if !any(s -> s.site_name == st.site_name || ES.gcdist(s, st) < MERGE_KM, have)])
    end
    esl_name(cfg) in og && append!(sts, PM.read_esl(cfg["score"]["fesmdata"], cfg["esl"]))
    idx = [(g, findall(s -> s.group == g, sts)) for g in og if !(g in RSL_GROUPS)]
    none = [g for (g, ii) in idx if isempty(ii)]
    isempty(none) || error("[[gnss]]: no station selected in ", join(none, ", "))
    (obs = sts, idx = idx, dir = dir)
end
# the GNSS stations of an ensemble: <folder>/gnss.csv (group, site, lat, lon, rate, sd [mm/yr]), and per group a line in ensemble.md
function write_gnss(C, gd, cfg)
    open(joinpath(C, "gnss.csv"), "w") do io
        println(io, "group,site,lat,lon,rate,sd")
        for s in gd.obs; s.kind == :gnss && println(io, join((s.group, s.site_name, s.lat, s.lon, @sprintf("%.3f", s.value), @sprintf("%.3f", s.sd)), ",")); end
    end
    open(joinpath(C, "ensemble.md"), "a") do io
        println(io, "\n\nGNSS groups (", gd.dir, "; the stations in gnss.csv):\n")
        for (g, ii) in gd.idx
            gd.obs[ii[1]].kind == :gnss || continue
            println(io, "- ", g, ": ", join(["GPS-" * x["product"] * " (rate " * get(x, "rate", "corrected") * ")"
                                             for x in cfg["gnss"] if x["name"] == g], " + "), ", ", length(ii), " stations")
        end
    end
end
# the model's values at the points of gd for the run folder dir: the rates [mm/yr] at the GNSS stations (pm_compare.jl read_rate;
# NaN without dz_bed_rate) and ESL [m] at the ages of the ESL group (esl_at); nothing when it has none
function gnss_model(dir, gd)
    isempty(gd.obs) && return nothing
    v = fill(NaN, length(gd.obs)); ig = findall(s -> s.kind == :gnss, gd.obs); ie = findall(s -> s.kind == :esl, gd.obs)
    if !isempty(ig)
        R = PM.run_done(dir) ? PM.read_rate(dir) : nothing                   # the present-day rate: of a finished run only
        R === nothing || (v[ig] = PM.rate_at(R, [gd.obs[i].lon for i in ig], [gd.obs[i].lat for i in ig]))
    end
    isempty(ie) || (v[ie] = PM.esl_at(dir, [gd.obs[i].t for i in ie]))
    all(isnan, v) ? nothing : v
end
# the misfit per GNSS group of the model rates G (nothing: no group): d = model - obs, sigma_i the station's sd
function gnss_residuals(gd, G)
    out = Dict{String,Any}()
    G === nothing && return out
    for (g, ii) in gd.idx
        k = [i for i in ii if isfinite(G[i])]
        isempty(k) || (out[g] = (pts = gd.obs[k], d = G[k] .- [gd.obs[i].value for i in k], σ = [gd.obs[i].sd for i in k]))
    end
    out
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
# the error model per data group: <ensemble folder>/errmodel.toml (plastic_ensemble.jl errmodel) when it exists, else the config's;
# the config's too for a group the file does not have (a GNSS group added after the fit). No corr_time_kyr: no time term (GNSS).
function error_model(cfg, C)
    f = joinpath(C, "errmodel.toml")
    lk = isfile(f) ? TOML.parsefile(f)["likelihood"] : cfg["likelihood"]
    add = isfile(f) ? sort([g for g in keys(get(cfg, "likelihood", Dict())) if !haskey(lk, g)]) : String[]
    for g in add; lk[g] = cfg["likelihood"][g]; end
    Dict(g => (p["sigma_m"], p["corr_dist_km"], 1000.0 * get(p, "corr_time_kyr", Inf)) for (g, p) in lk),
    (isfile(f) ? f : "config") * (isempty(add) ? "" : "; " * join(add, ", ") * ": config")
end
# the implausibility of a group's misfit: (chi2 - n) / sqrt(2 n), the standardized chi2 (N(0, 1) for large n when R is right)
implausibility(χ2, n) = (χ2 - n) / sqrt(2n)

# the members of wave k that ran, their residuals per group (RSL from S; GNSS from G, the model rates at the stations, nothing without
# dz_bed_rate) and ESL at 20 ka; gd the GNSS data (gnss_data)
function wave_residuals(C, k, cfg)
    og = data_groups(cfg); gd = gnss_data(cfg)
    obs = [o for o in PM.read_indicators(cfg["score"]["fesmdata"]) if o.type == 0 && any(g -> PM.in_group(o, g), og)]
    idx = [(g, findall(o -> PM.in_group(o, g), obs)) for g in og if g in RSL_GROUPS]
    A = sort(NCDataset(ds -> Float64.(ds["time"][:]),
                       joinpath(cfg["score"]["fesmdata"], "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc")))
    ids, Θ = ES.read_theta(C, k)
    res = Dict{Int,Any}()
    for m in ids
        r, why = ES.member_series(ES.mdir(C, k, m), obs, A)
        r === nothing && (println("member $m: $why"); continue)
        G = gnss_model(ES.mdir(C, k, m), gd)
        res[m] = (esl = r.esl, rg = merge(group_residuals(obs, idx, r.S, A), gnss_residuals(gd, G)), S = r.S, G = G)
    end
    isempty(gd.idx) || println("GNSS: ", count(m -> res[m].G !== nothing, keys(res)), " of ", length(res), " members with dz_bed_rate")
    ids, Θ, res, og, obs, A, gd
end

# --- the references --------------------------------------------------------------------------------------------------------------
# The reconstructions the members are compared with ([[references]] of the config): per entry a name and one of
#   product  a FesmData IceHistory file (under [score] fesmdata) with dz_sl, its own GIA (PaleoMIST: SELEN, Gowan's 1-D Earth)
#   topo     a FesmData IceHistory file with the topography relative to the sea level of the time (GLAC1D: its own GIA, PM.rsl_topo)
#   run      a standalone VILMA2 run under a reconstruction's ice: a gia_refcorr.sh case (its converged iteration, r.nc -> r_<n>.nc:
#            it<n>) or a folder with vilma_out.nc ($VAR in the path from the environment)
# The RSL at the index points (npoints x nages), NaN at the ages the reference does not have; nothing (and why) when it is missing.
expandenv(p) = replace(p, r"\$\{?(\w+)\}?" => m -> ENV[strip(m, ['$', '{', '}'])])
# the folder of a run reference: a gia_refcorr.sh case's converged iteration, or the folder itself
function ref_dir(ref)
    d = expandenv(ref["run"]); rc = joinpath(d, "r.nc")
    islink(rc) && (d = joinpath(d, "it" * match(r"r_(\d+)\.nc", readlink(rc))[1]))
    d
end
function reference_series(ref, fes, obs, A)
    plon = [o.lon for o in obs]; plat = [o.lat for o in obs]
    if haskey(ref, "run")
        d = ref_dir(ref)
        PM.run_done(d) || return nothing, "$d: not finished (no status 0)"
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
# the ESL of a published reconstruction at the points of the ESL group: from the reference's esl = "<SeaLevel product>" (pm_compare.jl
# sealevel_at; SeaLevel-PaleoMIST for pub-PMIST); nothing without one
function esl_reference(cfg, ref, gd)
    (haskey(ref, "esl") && any(s -> s.kind == :esl, gd.obs)) || return nothing
    v = fill(NaN, length(gd.obs)); ie = findall(s -> s.kind == :esl, gd.obs)
    v[ie] = PM.sealevel_at(cfg["score"]["fesmdata"], ref["esl"], [gd.obs[i].t for i in ie])[1]
    v
end
# the model rates at the GNSS stations of every reference (gnss_model): nothing for a published reconstruction (product, topo) and a
# run without dz_bed_rate
function reference_rates(cfg, gd)
    G = Dict{String,Any}(ref["name"] => (haskey(ref, "run") ? gnss_model(ref_dir(ref), gd) : esl_reference(cfg, ref, gd))
                         for ref in get(cfg, "references", []))
    isempty(gd.idx) || println("GNSS: references with dz_bed_rate: ", join(sort([n for (n, g) in G if g !== nothing]), ", "))
    G
end

# A scoring window: the index points of age <= tmax [yr] and the records of age <= tmax (the references that end there: GLAC1D's own
# GIA at 26 ka, so tmax = 25 ka, its last record of the 2.5 kyr ages of the scores). The series of all, members and references, are
# scored on the same window.
function score_window(obs, A, og, tmax)
    kp = findall(o -> o.t <= tmax, obs); ka = findall(<=(tmax), A)
    (tmax = tmax, kp = kp, ka = ka, A = A[ka], obs = obs[kp], idx = [(g, findall(o -> PM.in_group(o, g), obs[kp])) for g in og if g in RSL_GROUPS])
end
# per data group (log L, chi2, n, I) of a series S (all points x all ages) on the window W; nothing when S misses a record there (the
# points where S has gaps are dropped, as in group_residuals). The GNSS groups (present day, in every window) from the series' model
# rates G at the stations of gd: none when G is nothing (no dz_bed_rate), so a series may lack the GNSS groups.
function score_series(S, W, em, og; G = nothing, gd = nothing)
    Sw = S[W.kp, W.ka]; any(k -> all(isnan, view(Sw, :, k)), axes(Sw, 2)) && return nothing
    rg = group_residuals(W.obs, W.idx, Sw, W.A)
    gd === nothing || merge!(rg, gnss_residuals(gd, G))
    Dict(g => (v = loglik_group(rg[g], em[g]...); (v..., implausibility(v[2], v[3]))) for g in og if haskey(rg, g))
end
# a value of a series' score s in group g (i: 1 log L, 2 chi2, 3 n, 4 I), NaN (n: 0) when it has no score there
sv(s, g, i) = haskey(s, g) ? s[g][i] : (i == 3 ? 0 : NaN)
fmtg(s, g, i, f = "%.1f") = haskey(s, g) ? Printf.format(Printf.Format(f), s[g][i]) : "-"

# The error model from a wave: per data group, sigma_m, L_x, L_t by maximum likelihood over a grid, from the residuals of the
# `top` members that fit that group best under the current error model (the members the data do not rule out; the rest would
# inflate the error by their own misfit). A GNSS group: sigma_m [mm/yr] and L_x only (no time term), from the members with a score
# there (none: not fitted, the previous error model kept). Writes <folder>/errmodel.toml (read by weights from then on) and errmodel.md.
const GRID_SIGMA = [3.0, 5.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 23.0, 26.0, 30.0, 35.0, 40.0, 50.0, 65.0, 80.0]   # [m]
const GRID_LX = [10.0, 25.0, 50.0, 100.0, 175.0, 250.0, 350.0, 500.0, 750.0, 1000.0, 1500.0, 2500.0]   # [km]
const GRID_LT = [1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 80.0, 160.0]                            # [kyr]
const GRID_SIGMA_GNSS = [0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0]                # [mm/yr]
const GRID_SIGMA_ESL = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 30.0]                       # [m]
function errmodel(C, k; top = 6)
    C = abspath(C); cfg = ES.config(C); ids, _, res, og = wave_residuals(C, k, cfg)
    em, src = error_model(cfg, C); ok = sort(collect(keys(res)))
    lk = Dict{String,Any}(); rows = []
    for g in og
        okg = [m for m in ok if haskey(res[m].rg, g)]
        if isempty(okg)
            println(g, ": no member has a score (no dz_bed_rate); not fitted")
            haskey(em, g) && (lk[g] = Dict{String,Any}("sigma_m" => em[g][1], "corr_dist_km" => em[g][2]);
                              isfinite(em[g][3]) && (lk[g]["corr_time_kyr"] = em[g][3] / 1000))
            continue
        end
        ll0 = Dict(m => loglik_group(res[m].rg[g], em[g]...) for m in okg)
        best = sort(okg, by = m -> -ll0[m][1])[1:min(top, length(okg))]
        # the distances [km] and age differences [yr] of each member's points, once
        DT = Dict(m => (p = res[m].rg[g].pts; (D = [ES.gcdist(a, b) for a in p, b in p], T = [abs(a.t - b.t) for a in p, b in p]))
                  for m in best)
        kd = group_kind(cfg, g)
        gs, glx, glt = kd == :gnss ? (GRID_SIGMA_GNSS, GRID_LX, [Inf]) : kd == :esl ? (GRID_SIGMA_ESL, [1.0], GRID_LT) : (GRID_SIGMA, GRID_LX, GRID_LT)
        fit = (-Inf, 0.0, 0.0, 0.0)
        for σm in gs, Lx in glx, Lt in glt
            L = sum(loglik_group(res[m].rg[g], σm, Lx, 1000Lt; DT = DT[m])[1] for m in best)
            L > fit[1] && (fit = (L, σm, Lx, Lt))
        end
        _, σm, Lx, Lt = fit
        edge = σm in (first(gs), last(gs)) || (length(glx) > 1 && Lx in (first(glx), last(glx))) || (length(glt) > 1 && Lt in (first(glt), last(glt)))
        χn(m, p...) = (v = loglik_group(res[m].rg[g], p...); v[2] / v[3])
        lk[g] = kd == :gnss ? Dict("sigma_m" => σm, "corr_dist_km" => Lx) : Dict("sigma_m" => σm, "corr_dist_km" => Lx, "corr_time_kyr" => Lt)
        println(g, ": sigma_m ", σm, kd == :gnss ? " mm/yr" : " m", ", L_x ", Lx, " km", kd == :gnss ? "" : ", L_t $(Lt) kyr", " (members ", best, ")")
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
                "(1 when the error model matches their misfit).", any(r -> !(r.g in RSL_GROUPS), rows) ?
                " GNSS groups: sigma_m in mm/yr, no L_t." : "", "
")
        println(io, "| group | n | sigma_m [m] | L_x [km] | L_t [kyr] | chi2/n | members |
|---|---|---|---|---|---|---|")
        for r in rows
            println(io, "| ", r.g, " | ", r.n, " | ", @sprintf("%g → %g", r.old[1], r.new[1]), " | ", @sprintf("%g → %g", r.old[2], r.new[2]),
                    " | ", isfinite(r.old[3]) || isfinite(r.new[3]) ? @sprintf("%g → %g", r.old[3] / 1000, r.new[3]) : "-", " | ", @sprintf("%.2f → %.2f", r.χold, r.χnew), " | ",
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
# GNSS groups: a group that judges loses its judging role when a member that ran has no score there (no dz_bed_rate); a group where no
# bar reference has a score keeps the plain cut. Writes iter_kk/weights.csv, weights.md and plausible.txt (the folders of the
# plausible members, one per line).
function weights(C, k)
    C = abspath(C); cfg = ES.config(C); ids, Θ, res, og, obs, A, gd = wave_residuals(C, k, cfg); en = cfg["ensemble"]
    em, src = error_model(cfg, C); cut0 = get(en, "implausibility", 3.0)
    ug0 = String.(get(en, "use_groups", judging_groups(cfg)))     # the groups that judge: the cut and the weights; the rest reported
    all(in(og), ug0) || error("use_groups $(ug0) not all in obs_groups $(og)")
    sg = ES.scale_groups(cfg["scale"]); P = ES.params(cfg, sg)
    W = score_window(obs, A, og, 1000.0 * get(en, "compare_max_age_ka", Inf))
    ok = sort(collect(keys(res)))
    sc = Dict(m => score_series(res[m].S, W, em, og; G = res[m].G, gd = gd) for m in ok)   # g => (log L, chi2, n, I)
    nos = [g for g in ug0 if !all(m -> haskey(sc[m], g), ok)]      # judging groups without a score in every member
    ug = [g for g in ug0 if !(g in nos)]
    refs = references(cfg, obs, A); rG = reference_rates(cfg, gd)
    rsc = [(nm, s) for (nm, s) in ((nm, score_series(S, W, em, og; G = rG[nm], gd = gd)) for (nm, S) in refs) if s !== nothing]
    # the bar: per group the lowest I_g of the bar's references that have a score there, and which one
    bar = String.(get(en, "bar", String[])); nob = setdiff(bar, first.(rsc))
    isempty(nob) || error("bar references not scored on the window: ", join(nob, ", "))
    best = Dict(g => (b = [(s[g][4], nm) for (nm, s) in rsc if nm in bar && haskey(s, g)]; isempty(b) ? (-Inf, "") : minimum(b)) for g in og)
    cut = Dict(g => max(cut0, best[g][1]) for g in og)
    plaus = [m for m in ok if all(sc[m][g][4] <= cut[g] for g in ug)]
    # tempering: the scale factor of each group's error covariance (the groups with a score)
    sf = Dict(g => max(1.0, minimum(s[g][2] / s[g][3] for s in [[sc[m] for m in ok]; last.(rsc)] if haskey(s, g)))
              for g in og if any(m -> haskey(sc[m], g), ok) || any(r -> haskey(r[2], g), rsc))
    LL = Dict(m => sum(sc[m][g][1] for g in ug) for m in ok)
    LLt = Dict(m => sum(sc[m][g][1] + 0.5 * sc[m][g][2] * (1 - 1 / sf[g]) - 0.5 * sc[m][g][3] * log(sf[g]) for g in ug) for m in ok)
    w = Dict(m => 0.0 for m in ok)
    if !isempty(plaus)
        mx = maximum(LLt[m] for m in plaus); for m in plaus; w[m] = exp(LLt[m] - mx); end
        t = sum(values(w)); for m in plaus; w[m] /= t; end
    end
    ess = isempty(plaus) ? 0.0 : 1 / sum(w[m]^2 for m in plaus)
    It = ES.itdir(C, k)
    I(m, g) = sv(sc[m], g, 4)
    open(joinpath(It, "weights.csv"), "w") do io
        println(io, "member,ok,", join(["logL_$g" for g in og], ","), ",logL,logL_tempered,", join(["chi2n_$g" for g in og], ","), ",",
                join(["I_$g" for g in og], ","), ",plausible,esl_20ka,weight")
        for m in ids
            haskey(res, m) || (println(io, m, ",0", repeat(",NaN", 3length(og) + 5)); continue)
            println(io, m, ",1,", join([@sprintf("%.2f", sv(sc[m], g, 1)) for g in og], ","), ",", @sprintf("%.2f,%.2f", LL[m], LLt[m]), ",",
                    join([@sprintf("%.3f", sv(sc[m], g, 2) / sv(sc[m], g, 3)) for g in og], ","), ",", join([@sprintf("%.2f", I(m, g)) for g in og], ","),
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
                join(ug, ", "), length(ug0) < length(og) ? "; reported only: " * join(setdiff(og, ug0), ", ") : "",
                isempty(nos) ? "" : "; not judging, without a score in every member (no dz_bed_rate: no GNSS score): " * join(nos, ", "), ".\n")
        println(io, "Cut per group (the larger of $cut0 and the best I_g of the bar: ", isempty(bar) ? "none" : join(bar, ", "),
                "); scale factor s_g of the tempered weights:\n")
        println(io, "| group | n | cut | set by | s_g |\n|---|---|---|---|---|")
        for g in og
            m1 = findfirst(m -> haskey(sc[m], g), ok)
            println(io, "| ", g, " | ", m1 === nothing ? "-" : sc[ok[m1]][g][3], " | ", @sprintf("%.2f", cut[g]), " | ", cut[g] > cut0 ? best[g][2] : "-",
                    " | ", haskey(sf, g) ? @sprintf("%.2f", sf[g]) : "-", " |")
        end
        println(io, "\nPlausible (I_g <= cut in every group that judges): $(length(plaus)) members", isempty(plaus) ? "" :
                @sprintf(", ESL at 20 ka %.1f to %.1f m; effective sample size of their weights %.2f", minimum(eslp), maximum(eslp), ess), ".\n")
        println(io, "Ruled out by group (I_g > cut): ", join([@sprintf("%s %d", g, count(m -> haskey(sc[m], g) && I(m, g) > cut[g], ok)) for g in og], ", "), "\n")
        println(io, "| member | plausible | weight | log L | log L tempered | ", join(["I $g" for g in og], " | "), " | ESL 20 ka [m] |\n|---|---|---|---|---|",
                repeat("---|", length(og) + 1))
        for m in sort(ok, by = m -> (m in plaus ? 0 : 1, -LLt[m]))
            println(io, "| m", lpad(m, 3, "0"), " | ", m in plaus ? "yes" : "no", " | ", @sprintf("%.3g", w[m]), " | ", @sprintf("%.1f", LL[m]),
                    " | ", @sprintf("%.1f", LLt[m]), " | ", join([fmtg(sc[m], g, 4) for g in og], " | "), " | ", @sprintf("%.1f", res[m].esl), " |")
        end
        if !isempty(plaus)
            println(io, "\nWeighted parameters of the plausible members (tau_b: log factor; Earth: dex; prior mean 0): mean (sd)\n")
            for (q, n) in enumerate(names); println(io, "- ", n, ": ", @sprintf("%+.2f (%.2f)", μ[q], sd[q])); end
        end
    end
    compare_references(It, cfg, res, og, obs, A, em, ug, cut, refs; gd = gd, rG = rG)
    println("wrote ", joinpath(It, "weights.md"), "; ", length(plaus), " plausible, ESS ", @sprintf("%.2f", ess))
end

# The members against the references (references), on the window of [ensemble] compare_max_age_ka (25 ka: GLAC1D's own GIA ends at
# 26 ka) and on all points (the references that cover them): per series the log L of the groups that judge, I_g per group and the
# groups passed (I_g <= cut[g], weights). The references, the 5 members of the highest log L and the members' median. A series without
# a score in a group (no dz_bed_rate: no GNSS score) shows "-" there (references.csv: n 0, I NaN), and its log L and passed are over
# the groups that judge it has a score in. Appends to weights.md, writes references.csv (every series).
function compare_references(It, cfg, res, og, obs, A, em, ug, cut, refs; gd = nothing, rG = Dict{String,Any}())
    best1 = Int[]
    tw = 1000.0 * get(cfg["ensemble"], "compare_max_age_ka", 25.0)
    wins = [(@sprintf("index points of age <= %g ka", tw / 1000), score_window(obs, A, og, tw)), ("all index points", score_window(obs, A, og, Inf))]
    LL(s) = sum(s[g][1] for g in ug if haskey(s, g)); np(s) = count(g -> haskey(s, g) && s[g][4] <= cut[g], ug)
    nj(s) = count(g -> haskey(s, g), ug)
    io = open(joinpath(It, "weights.md"), "a"); csv = open(joinpath(It, "references.csv"), "w")
    println(csv, "window,series,reference,", join(["n_$g" for g in og], ","), ",", join(["I_$g" for g in og], ","), ",logL,passed")
    for (wn, (lbl, W)) in enumerate(wins)
        rs = [(nm, s) for (nm, s) in ((nm, score_series(S, W, em, og; G = get(rG, nm, nothing), gd = gd)) for (nm, S) in refs) if s !== nothing]
        ms = Dict(m => score_series(res[m].S, W, em, og; G = get(res[m], :G, nothing), gd = gd) for m in keys(res))
        for (nm, s, isref) in [[(nm, s, true) for (nm, s) in rs]; [("m" * lpad(m, 3, "0"), ms[m], false) for m in sort(collect(keys(ms)))]]
            println(csv, wn == 1 ? "window" : "all", ",", nm, ",", isref ? 1 : 0, ",", join([sv(s, g, 3) for g in og], ","), ",",
                    join([@sprintf("%.2f", sv(s, g, 4)) for g in og], ","), ",", @sprintf("%.2f", LL(s)), ",", np(s))
        end
        best = first(sort(collect(keys(ms)), by = m -> -LL(ms[m])), 5); wn == 1 && (best1 = best)
        lack = any(s -> any(g -> !haskey(s, g), og), [last.(rs); collect(values(ms))])
        println(io, "\n## Against the references: ", lbl, "\n\nlog L and groups passed over the groups that judge (", join(ug, ", "),
                "); I_g per group (passed at <= the cut of the group, above). References missing here: ",
                isempty(setdiff(first.(refs), first.(rs))) ? "none" : join(setdiff(first.(refs), first.(rs)), ", "), ".",
                lack ? " -: no score in the group (no dz_bed_rate: no GNSS score); log L and passed are over the groups the series has." : "", "\n")
        println(io, "| series | log L | ", join(["I $g" for g in og], " | "), " | passed |\n|---|---|", repeat("---|", length(og) + 1))
        row(nm, s) = println(io, "| ", nm, " | ", @sprintf("%.1f", LL(s)), " | ", join([fmtg(s, g, 4) for g in og], " | "), " | ",
                             np(s), "/", nj(s), " |")
        for (nm, s) in rs; row("**" * nm * "**", s); end
        for m in best; row("m" * lpad(m, 3, "0"), ms[m]); end
        med(g) = (v = [s[g][4] for s in values(ms) if haskey(s, g)]; isempty(v) ? "-" : @sprintf("%.1f", median(v)))
        npts(g) = (v = [s[g][3] for s in values(ms) if haskey(s, g)]; isempty(v) ? "-" : string(first(v)))
        println(io, "| members' median | ", @sprintf("%.1f", median(LL(s) for s in values(ms))), " | ",
                join([med(g) for g in og], " | "), " | |")
        println(io, "| points (members) | | ", join([npts(g) for g in og], " | "), " | |")
        println(io, "\nMembers passing each group: ", join([@sprintf("%s %d", g, count(s -> haskey(s, g) && s[g][4] <= cut[g], values(ms))) for g in og], ", "),
                "; all groups that judge: ", count(s -> np(s) == length(ug), values(ms)), " of ", length(ms), ".")
    end
    # ESL: the data, the reference curves of [esl] reference_curves (SeaLevel products, not scored: Lambeck et al. 2014 uses far-field
    # RSL that the far-field group already has), the references and the 5 members of the highest log L on the window
    ie = gd === nothing ? Int[] : findall(s -> s.kind == :esl, gd.obs)
    if !isempty(ie)
        ages = [gd.obs[i].t for i in ie]; f0(v) = join([isfinite(x) ? @sprintf("%.0f", x) : "-" for x in v], " | ")
        println(io, "\n## ESL [m], relative to the present: the data (", cfg["esl"]["product"], ", 1 sigma), reference curves (not scored), ",
                "the references and the 5 best members\n")
        println(io, "| series | ", join([@sprintf("%g ka", a / 1000) for a in ages], " | "), " |\n|---|", repeat("---|", length(ages)))
        println(io, "| **data** | ", join([@sprintf("%.0f ± %.0f", gd.obs[i].value, gd.obs[i].sd) for i in ie], " | "), " |")
        for c in get(cfg["esl"], "reference_curves", String[])
            println(io, "| ", c, " (curve) | ", f0(PM.sealevel_at(cfg["score"]["fesmdata"], c, ages)[1]), " |")
        end
        for (nm, G) in sort(collect(rG), by = first)
            (G === nothing || all(isnan, G[ie])) || println(io, "| **", nm, "** | ", f0(G[ie]), " |")
        end
        for m in best1
            G = get(res[m], :G, nothing); G === nothing || println(io, "| m", lpad(m, 3, "0"), " | ", f0(G[ie]), " |")
        end
    end
    close(io); close(csv)
end

# The scores of every series, references and members, for tables (docs/simulations.md): per window (the index points of age <=
# compare_max_age_ka; all points) and data group, and over all groups: n, the RMSE of d [m], the RMS of z = d / sigma_i (the classic
# score), chi2/n, I_g and log L under the error model. Writes iter_kk/scores.csv. A reference without records in a window is left out.
# GNSS groups (the same in both windows): their rows after the RSL groups' row "all" (the RSL groups only), rmse_m then in mm/yr, and
# their own total "all GNSS"; a series without dz_bed_rate has no GNSS rows.
#   julia --project=tools/levante plastic_ensemble.jl scores <ensemble folder> [<wave>]
function scores(C, k)
    C = abspath(C); cfg = ES.config(C); ids, Θ, res, og, obs, A, gd = wave_residuals(C, k, cfg)
    em, src = error_model(cfg, C); refs = references(cfg, obs, A); rG = reference_rates(cfg, gd)
    tw = 1000.0 * get(cfg["ensemble"], "compare_max_age_ka", 25.0)
    series = [[(nm, S, 1, rG[nm]) for (nm, S) in refs]; [("m" * lpad(m, 3, "0"), res[m].S, 0, res[m].G) for m in sort(collect(keys(res)))]]
    f = joinpath(ES.itdir(C, k), "scores.csv")
    open(f, "w") do io
        println(io, "window,series,reference,group,n,rmse_m,rmsz,chi2n,I,logL")
        for (wn, W) in (("window", score_window(obs, A, og, tw)), ("all", score_window(obs, A, og, Inf)))
            for (nm, S, isref, G) in series
                Sw = S[W.kp, W.ka]; any(c -> all(isnan, view(Sw, :, c)), axes(Sw, 2)) && continue
                for (rg, tname) in ((group_residuals(W.obs, W.idx, Sw, W.A), "all"), (gnss_residuals(gd, G), "all points"))
                    tot = (n = 0, d2 = 0.0, z2 = 0.0, chi2 = 0.0, ll = 0.0)
                    for g in og
                        haskey(rg, g) || continue
                        r = rg[g]; ll, χ2, n = loglik_group(r, em[g]...); z = r.d ./ r.σ
                        @printf(io, "%s,%s,%d,%s,%d,%.2f,%.3f,%.3f,%.2f,%.2f\n", wn, nm, isref, g, n, sqrt(mean(abs2, r.d)), sqrt(mean(abs2, z)),
                                χ2 / n, implausibility(χ2, n), ll)
                        tot = (n = tot.n + n, d2 = tot.d2 + sum(abs2, r.d), z2 = tot.z2 + sum(abs2, z), chi2 = tot.chi2 + χ2, ll = tot.ll + ll)
                    end
                    tot.n > 0 && @printf(io, "%s,%s,%d,%s,%d,%.2f,%.3f,%.3f,,%.2f\n", wn, nm, isref, tname, tot.n, sqrt(tot.d2 / tot.n),
                                         sqrt(tot.z2 / tot.n), tot.chi2 / tot.n, tot.ll)
                end
            end
        end
    end
    println("wrote ", f, " (error model: ", src, ")")
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
    if args[1] == "earth-basis"
        kv = Dict(split(x, "=") for x in args[5:end])
        return ES.earth_basis(args[2], args[3], args[4]; kmax = parse(Int, get(kv, "kmax", "6")), lo = parse(Float64, get(kv, "lo", "19.5")),
                              hi = parse(Float64, get(kv, "hi", "24")), weights = String(get(kv, "weights", "field")),
                              lid = haskey(kv, "lid") ? parse(Float64, kv["lid"]) : nothing)
    end
    args[1] == "earth-field" && return ES.earth_field(args[2], args[3], parse.(Float64, args[4:end]))
    if args[1] == "earth-sample"
        kv = Dict(split(x, "=") for x in args[5:end])
        return ES.earth_sample(args[2], args[3], args[4]; n = parse(Int, get(kv, "n", "30")), K = parse(Int, get(kv, "K", "4")),
                               seed = parse(Int, get(kv, "seed", "1")))
    end
    args[1] == "earth-radial" && return ES.earth_radial(args[2], args[3])
    if args[1] == "earth-search-sample"
        kv = Dict(split(x, "=") for x in args[4:end])
        return ES.earth_search_sample(args[2], args[3]; n = parse(Int, get(kv, "n", "200")), seed = parse(Int, get(kv, "seed", "1")))
    end
    args[1] == "earth-search-par" && return println("wrote ", ES.earth1d_par(args[2], args[3], parse.(Float64, args[4:8])))
    args[1] == "weights" && return weights(args[2], length(args) >= 3 ? parse(Int, args[3]) : 1)
    args[1] == "scores" && return scores(args[2], length(args) >= 3 ? parse(Int, args[3]) : 1)
    args[1] == "errmodel" && return errmodel(args[2], length(args) >= 3 ? parse(Int, args[3]) : 1;
                                             top = length(args) >= 4 ? parse(Int, args[4]) : 6)
    error("plastic_ensemble.jl init <config> <folder> | weights <folder> [<wave>] | errmodel <folder> [<wave>] [<top>] | " *
          "earth-regions <config> <out.nc> | earth-shift <config> <out.nc> UM=.. LM=.. <region>=.. | earth-radial <par> <out par> | " *
          "earth-search-sample <1-D par> <outdir> [n=200] [seed=1] | earth-search-par <1-D par> <out par> L um tz lm1 lm2")
end
abspath(PROGRAM_FILE) == @__FILE__() && main(ARGS)
