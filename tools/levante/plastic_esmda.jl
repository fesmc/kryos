# ES-MDA (Emerick and Reynolds 2013) calibration of tau_b of the plastic ice backend against the GAPSLIP 2.1 sea-level index points,
# with VILMA2 in the loop (docs/dev/plastic-backend.md, "Calibration"; the jobs: plastic_esmda.sh):
#   julia --project=<tools> plastic_esmda.jl init <config.toml> <calibration folder>   prior ensemble = iteration 1
#   julia --project=<tools> plastic_esmda.jl update <calibration folder> <k>            score iteration k, update, write k + 1
#   julia --project=<tools> plastic_esmda.jl selftest                                   the update against a linear-Gaussian posterior
# Parameters theta: the log multipliers of tau_b per scale group and node (plastic_taub_scale.jl), group-major; nodes of a group may be
# tied to one parameter (`tie` of [[scale.groups]], lists of node ages: the nodes without direct data), and the nodes of `pinned_ka`
# ([scale]) are no parameter: log factor 0 in every group (the present day, fitted in Stage 1); prior N(0, C), C of a group
# s(t_p) s(t_q) exp(-|t_p - t_q| / L) with t_p the mean age of a parameter's nodes and s(t) = sigma sin^2(pi t / (2 ramp_ka)) below
# ramp_ka, sigma above (0 and flat at 0 ka, smooth), groups independent; every member clipped to |theta| <= theta_max after an
# update. Observations: the index points of pm_compare.jl's GROUPS in obs_groups, the model RSL at a point the scorer's (series_run:
# the run's VILMA2 stream geo_glob_vilma.nc read by pm_compare.jl's read_run, bilinear, linear in time), so d = model - rsl is the d of the scores.
# Their error covariance R = diag(sigma_i^2) + sigma_m^2 C, sigma_i the scorer's ip_sigma (elevation on the side of d, and the age
# through the model slope) at the ensemble-mean RSL of the iteration, fixed within the iteration; sigma_m the model error the indicator
# lacks (GIA and Earth model, resolution, RSL interpolation), correlated between index points: C_ij = exp(-dist_ij / L_s)
# exp(-|t_i - t_j| / L_t), dist the great-circle distance, t the indicator age (corr_dist_km, corr_time_kyr; 0: independent, the
# sigma_add of v1). A gross-error cut, fixed in iteration 1 at the prior mean member (m000): index points with
# |d| / sqrt(sigma_i^2 + sigma_m^2) > z_cut there (the diagonal of R) are not assimilated (obs_kept.csv); every score is over all
# points. The update with the full R: R = L L' (Cholesky, once per iteration), predictions and observations whitened by L^-1, then
# perturbed observations N(0, alpha I) in the whitened space (= N(0, alpha R)) and the update in the ensemble subspace (Woodbury: an
# Ne x Ne solve). Limiting data are not assimilated: their violation fractions are diagnostics.
# Folder: config.toml, kryos.x, input/, par.nml, groups.csv, obs.csv, obs_kept.csv, esmda.md (the report), next (the iteration to run, or done);
# iter_kk/: theta.csv (all members of the iteration), members.txt (those to run), mNNN/ (taub_scale.nc, par.nml, the run, status:
# exit code and seconds), scores.csv (per member and set), members.csv (per member: ok, chi2, ESL at 20 ka, run time).
# Every member writes the product (PRODUCT_VARS): the 2-D ice of each instance of scale_keys (<instance>.nc on its out_2D schedule,
# 2.5 kyr) and its 0-D series, besides the VILMA2 stream (RSL, bed, ocean, load) of the par file; [run] set and final_set come after.
# Member 0 is the mean member: the prior mean (0, Gowan's tau_b) in iteration 1, the ensemble mean after an update.
using LinearAlgebra, Printf, Random, Statistics, TOML, NCDatasets

module PM
include("pm_compare.jl")
end
include("plastic_taub_scale.jl")

const ESL_AGE = 20000.0
# the exports every member writes per plastic instance: the 2-D ice (thickness, surface, bed and sea level as the ice saw them, ice
# fraction) and the ice mass (out_vars filters the 0-D series too)
const PRODUCT_VARS = ["H_ice", "z_srf", "z_bed", "z_sl", "f_ice", "M_ice"]
product_set(keys) = [s for inst in unique(first.(split.(keys, "."))) for s in
                     ("$(inst).out_2D=\"out_2D\"", "$(inst).out_vars=" * join(("\"$v\"" for v in PRODUCT_VARS), " "))]

itdir(C, k) = joinpath(C, @sprintf("iter_%02d", k))
mdir(C, k, m) = joinpath(itdir(C, k), @sprintf("m%03d", m))
config(C) = TOML.parsefile(joinpath(C, "config.toml"))

# the parameters: (group, its node indices, name, mean age [ka]), group-major, nodes in the order of nodes_ka; `tie` joins nodes
function params(cfg, sg)
    P = NamedTuple{(:g, :ks, :name, :t),Tuple{Int,Vector{Int},String,Float64}}[]
    pinned = Float64.(get(cfg["scale"], "pinned_ka", Float64[]))
    all(a -> a in sg.nodes_ka, pinned) || error("a pinned node is not in nodes_ka: $(pinned)")
    for (g, spec) in enumerate(cfg["scale"]["groups"])
        ties = [Float64.(t) for t in get(spec, "tie", Any[])]
        any(t -> any(in(pinned), t), ties) && error("group $(sg.names[g]): a tied node is pinned")
        seen = Set{Int}()
        for (k, a) in enumerate(sg.nodes_ka)
            (k in seen || a in pinned) && continue
            tset = findfirst(t -> a in t, ties)
            ks = tset === nothing ? [k] : [findfirst(==(x), sg.nodes_ka) for x in ties[tset]]
            any(isnothing, ks) && error("group $(sg.names[g]): a tied node is not in nodes_ka: $(ties[tset])")
            union!(seen, ks)
            push!(P, (g = g, ks = ks, name = sg.names[g] * "_" * join([@sprintf("%g", sg.nodes_ka[i]) for i in ks], "-"),
                      t = mean(sg.nodes_ka[ks])))
        end
    end
    P
end
pnames(P) = [p.name for p in P]
# log_factor (group, node) of a parameter vector (0 at the pinned nodes)
function log_factor(sg, P, θ)
    lf = zeros(ngroups(sg), nnodes(sg))
    for (q, p) in enumerate(P), k in p.ks; lf[p.g, k] = θ[q]; end
    lf
end
logline(C, s...) = open(io -> println(io, s...), joinpath(C, "esmda.md"), "a")

# --- prior --------------------------------------------------------------------------------------------------------------------------
# the prior sd at age t [ka]: sigma, ramped smoothly from 0 at 0 ka to sigma at ramp_ka (sin^2; no ramp_ka: sigma)
prior_sd(pr, t) = (r = get(pr, "ramp_ka", 0.0); r > 0 && t < r ? pr["sigma"] * sin(π * t / (2r))^2 : float(pr["sigma"]))
# groups independent; within a group exponential in time, the group's corr_time_kyr ([[scale.groups]]) or the prior's
function prior_cov(P, pr, specs)
    L = [Float64(get(s, "corr_time_kyr", pr["corr_time_kyr"])) for s in specs]
    [p.g == q.g ? prior_sd(pr, p.t) * prior_sd(pr, q.t) * exp(-abs(p.t - q.t) / L[p.g]) : 0.0 for p in P, q in P]
end

# --- observation error ------------------------------------------------------------------------------------------------------------
const R_EARTH_KM = 6371.0
function gcdist(a, b)            # great-circle distance [km] (haversine)
    φ1, φ2 = deg2rad(a.lat), deg2rad(b.lat); dφ = φ2 - φ1; dλ = deg2rad(b.lon - a.lon)
    2R_EARTH_KM * asin(min(1.0, sqrt(sin(dφ / 2)^2 + cos(φ1) * cos(φ2) * sin(dλ / 2)^2)))
end
# R = diag(σ^2) + σm^2 C, C_ij = exp(-dist_ij / Ls - |t_i - t_j| / Lt) (Ls, Lt [km, yr]; Ls = 0: C = I)
function obs_cov(pts, σ, σm, Ls, Lt)
    n = length(pts); R = zeros(n, n)
    for j in 1:n, i in j:n
        c = i == j ? 1.0 : Ls > 0 ? exp(-gcdist(pts[i], pts[j]) / Ls - abs(pts[i].t - pts[j].t) / Lt) : 0.0
        R[i, j] = R[j, i] = σm^2 * c + (i == j ? σ[i]^2 : 0.0)
    end
    R
end
# the model error of a config: sigma_m (v1: sigma_add, independent), its correlation lengths [km, yr]
obs_err(es) = (get(es, "sigma_m", get(es, "sigma_add", 0.0)), get(es, "corr_dist_km", 0.0), 1000.0 * get(es, "corr_time_kyr", 0.0))

# --- members ------------------------------------------------------------------------------------------------------------------------
# theta (P x n) of the members ids; run: the ids to run
function write_iteration(C, k, cfg, sg, P, ids, Θ, run; final = false)
    It = itdir(C, k); mkpath(It)
    open(joinpath(It, "theta.csv"), "w") do io
        println(io, "member,", join(pnames(P), ","))
        for (j, m) in enumerate(ids); println(io, m, ",", join([@sprintf("%.6f", x) for x in Θ[:, j]], ",")); end
    end
    for (j, m) in enumerate(ids)
        m in run || continue
        D = mdir(C, k, m); mkpath(D)
        f = write_scale(joinpath(D, "taub_scale.nc"), sg, log_factor(sg, P, Θ[:, j]))
        cp(joinpath(C, "par.nml"), joinpath(D, "par.nml"); force = true)
        keys = get(cfg["run"], "scale_keys", ["ice_nor.taub_scale_file"])
        set = [["$(key)=\"$f\"" for key in keys]; product_set(keys)]
        final && append!(set, String.(get(cfg["run"], "final_set", String[])))
        Base.run(`python3 $(cfg["run"]["set_nml"]) $(joinpath(D, "par.nml")) $set`)
    end
    open(io -> foreach(m -> println(io, m), run), joinpath(It, "members.txt"), "w")
end

function read_theta(C, k)
    lines = readlines(joinpath(itdir(C, k), "theta.csv"))
    rows = [parse.(Float64, split(l, ",")) for l in lines[2:end]]
    Int.(first.(rows)), reduce(hcat, [r[2:end] for r in rows])
end

# --- init ---------------------------------------------------------------------------------------------------------------------------
# The folder of a calibration or an ensemble (plastic_ensemble.jl): the config, kryos.x, input/, the par file with [run] set and
# the frozen domain factors, slurm.env, groups.csv, obs.csv of the observation groups og. Returns (cfg, sg, obs).
function prepare(cfgfile, C, og)
    isfile(joinpath(C, "config.toml")) && error("$C is a calibration already")
    mkpath(C); cp(cfgfile, joinpath(C, "config.toml"))
    cfg = config(C); rn = cfg["run"]; K = rn["kryos_dir"]
    cp(joinpath(K, "libkryos", "bin", "kryos.x"), joinpath(C, "kryos.x"))
    mkpath(joinpath(C, "input"))
    for f in readdir(joinpath(K, "input")); cp(joinpath(K, "input", f), joinpath(C, "input", f); follow_symlinks = true); end
    cp(joinpath(K, rn["par"]), joinpath(C, "par.nml"))
    isempty(rn["set"]) || Base.run(`python3 $(rn["set_nml"]) $(joinpath(C, "par.nml")) $(rn["set"])`)
    if !isempty(get(rn, "domain_factor_file", ""))          # Stage 1: frozen in the folder, set for every instance
        f = joinpath(C, "taub_domain_factor.nc"); cp(rn["domain_factor_file"], f; follow_symlinks = true)
        Base.run(`python3 $(rn["set_nml"]) $(joinpath(C, "par.nml")) $(["$(key)=\"$f\"" for key in rn["domain_factor_keys"]])`)
    end
    open(joinpath(C, "slurm.env"), "w") do io
        for (k, v) in (("THREADS", rn["threads"]), ("MEM", rn["mem"]), ("TIME", rn["time"]), ("PART", rn["partition"]),
                       ("ACCT", rn["account"]), ("THROTTLE", rn["throttle"]), ("MAPS", rn["maps"]))
            println(io, k, "=", v)
        end
    end

    sg = scale_groups(cfg["scale"]); write_groups_csv(joinpath(C, "groups.csv"), sg)
    obs = PM.read_indicators(cfg["score"]["fesmdata"])
    open(joinpath(C, "obs.csv"), "w") do io
        println(io, "group,n_index_points,n_marine_limiting,n_terrestrial_limiting")
        for g in og; println(io, g, ",", join((count(o -> o.group == g && o.type == t, obs) for t in (0, -1, 1)), ",")); end
    end
    cfg, sg, obs
end

function init(cfgfile, C)
    C = abspath(C)
    og = TOML.parsefile(cfgfile)["esmda"]["obs_groups"]
    cfg, sg, obs = prepare(cfgfile, C, og)
    rn = cfg["run"]
    es = cfg["esmda"]; ne = es["ne"]; P = params(cfg, sg)
    Cθ = prior_cov(P, cfg["prior"], cfg["scale"]["groups"])
    Θ = clamp.(cholesky(Symmetric(Cθ)).L * randn(MersenneTwister(cfg["prior"]["seed"]), size(Cθ, 1), ne), -cfg["prior"]["theta_max"],
               cfg["prior"]["theta_max"])
    ids = collect(0:ne)
    write_iteration(C, 1, cfg, sg, P, ids, hcat(zeros(size(Cθ, 1)), Θ), ids)
    write(joinpath(C, "next"), "1\n")

    σm, Ls, Lt = obs_err(es)
    logline(C, "# ES-MDA calibration of tau_b (plastic backend)\n\nconfig: ", cfgfile, "; ", length(sg.domain_id), " domains, ",
            ngroups(sg), " groups x ", nnodes(sg), " nodes (pinned at 0: ", get(cfg["scale"], "pinned_ka", "none"), " ka), ",
            size(Cθ, 1), " parameters (", join(pnames(P), ", "), "); prior sigma ", cfg["prior"]["sigma"], " (ramp from 0 at 0 ka to ",
            get(cfg["prior"], "ramp_ka", 0), " ka; prior sd of the parameters ", join([@sprintf("%.3f", sqrt(Cθ[q, q])) for q in axes(Cθ, 1)], ", "),
            "), correlation time ", cfg["prior"]["corr_time_kyr"], " kyr, |theta| <= ", cfg["prior"]["theta_max"],
            "; Ne = ", ne, ", alpha = ", es["alpha"], ", final: ", es["final"], "; model error sigma_m ", σm, " m, correlated over ", Ls,
            " km and ", Lt / 1000, " kyr; gross-error cut |z| > ", es["z_cut"], " at the prior mean member", 
            isempty(get(rn, "domain_factor_file", "")) ? "" : "; domain factors " * rn["domain_factor_file"], "\n")
    logline(C, "```\n", join(groups_summary(sg), "\n"), "\n```\n")
    logline(C, "Index points assimilated (GAPSLIP 2.1 primary, 0-80 ka) and limiting data (diagnostics) by group:\n\n",
            "| group | IP | ML | TL |\n|---|---|---|---|")
    for g in og; logline(C, "| ", g, " | ", join((count(o -> o.group == g && o.type == t, obs) for t in (0, -1, 1)), " | "), " |"); end
    println("init: ", C, ", iteration 1 with members 0-", ne)
end

# --- scoring ------------------------------------------------------------------------------------------------------------------------
# the model of member m: RSL series at the samples (npoints x nages), the ESL at ESL_AGE, or nothing if the run failed
function member_series(D, obs, A)
    st = joinpath(D, "status")
    isfile(st) || return nothing, "no status"
    s = split(read(st, String)); s[1] == "0" || return nothing, "exit $(s[1])"
    isfile(joinpath(D, "geo_glob_vilma.nc")) || return nothing, "no geo_glob_vilma.nc"
    r = PM.read_run(D)
    all(a -> any(x -> abs(x - PM.trec(a)) < 0.01, r.t), A) || return nothing, "records missing"
    S = PM.series_run(r, [o.lon for o in obs], [o.lat for o in obs], A)
    esl = r.bsl_ive[PM.tindex(r, PM.trec(ESL_AGE))] - r.bsl_ive[end]
    (S = S, esl = esl, sec = parse(Float64, s[2])), ""
end

function summaries(obs, sets, S, A)
    sc = [PM.score(o, A, S[i, :]) for (i, o) in enumerate(obs)]
    [PM.summary(sc, obs, idx) for (_, idx) in sets]
end

fmt(x) = isnan(x) ? "-" : @sprintf("%.2f", x)

# --- update -------------------------------------------------------------------------------------------------------------------------
# One ES-MDA step: parameters Θ (P x N), predictions D (m x N) and observations y (m) in units of sigma (R = I), inflation α:
# Θ + dΘ S' (S S' + α I)^-1 (y + sqrt(α) e - D), with dΘ, S the anomalies over sqrt(N - 1), as dΘ (I + G S)^-1 G (...),
# G = S'/α (Woodbury: an N x N solve).
function esmda_step(Θ, D, y, α, rng)
    N = size(Θ, 2)
    Res = y .+ sqrt(α) .* randn(rng, size(D)...) .- D
    dΘ = (Θ .- mean(Θ; dims = 2)) ./ sqrt(N - 1); S = (D .- mean(D; dims = 2)) ./ sqrt(N - 1)
    G = S' ./ α
    Θ .+ dΘ * ((I + G * S) \ (G * Res))
end

# a linear-Gaussian check: ES-MDA with alpha (9.333, 7, 4, 2) on a large ensemble against the exact posterior, with R = I and with a
# correlated R (diagonal plus an exponential covariance along a line, the update whitened by its Cholesky factor as in update); and
# one whitened step against the update written with R itself, Θ + dΘ S'(S S' + α R)^-1 (y + sqrt(α) L e - D), the same e
function esmda_direct(Θ, D, y, R, α, e)
    N = size(Θ, 2); L = cholesky(Symmetric(R)).L
    dΘ = (Θ .- mean(Θ; dims = 2)) ./ sqrt(N - 1); S = (D .- mean(D; dims = 2)) ./ sqrt(N - 1)
    Θ .+ dΘ * (S' * ((S * S' .+ α .* R) \ (y .+ sqrt(α) .* (L * e) .- D)))
end
function selftest()
    for corr in (false, true)
        rng = MersenneTwister(7); P, m, N = 6, 40, 4000
        H = randn(rng, m, P); C0 = Matrix(Diagonal(fill(0.09, P))); θt = 0.3 .* randn(rng, P)
        pos = sort(10 .* rand(rng, m)); s = 0.5 .+ rand(rng, m)
        R = corr ? Matrix(Diagonal(s .^ 2)) .+ exp.(-abs.(pos .- pos') ./ 1.5) : Matrix(1.0I, m, m)
        L = cholesky(Symmetric(R)).L
        y = H * θt .+ L * randn(rng, m)
        Cp = inv(inv(C0) + H' * (R \ H)); μp = Cp * (H' * (R \ y))
        Θ = sqrt(0.09) .* randn(rng, P, N)
        for α in (9.333, 7.0, 4.0, 2.0); Θ = esmda_step(Θ, L \ (H * Θ), L \ y, α, rng); end
        em = maximum(abs.(vec(mean(Θ; dims = 2)) .- μp) ./ sqrt.(diag(Cp))); es = maximum(abs.(sqrt.(diag(cov(Θ; dims = 2))) ./ sqrt.(diag(Cp)) .- 1))
        @printf("selftest (%s R): posterior mean error %.3f sd, sd ratio error %.3f (N = %d)\n", corr ? "correlated" : "identity", em, es, N)
        em < 0.1 && es < 0.1 || error("selftest failed")
        if corr                    # one step, whitened against direct, N = 30 < m
            Θ0 = sqrt(0.09) .* randn(rng, P, 30); D0 = H * Θ0
            e = randn(MersenneTwister(3), m, 30)
            Θw = esmda_step(Θ0, L \ D0, L \ y, 4.0, MersenneTwister(3))
            Θd = esmda_direct(Θ0, D0, y, R, 4.0, e)
            dmax = maximum(abs.(Θw .- Θd))
            @printf("selftest: whitened step against the update with R: max abs difference %.2e\n", dmax)
            dmax < 1e-10 || error("selftest failed: the whitened update differs from the update with R")
        end
    end
end

# The index points kept for the assimilation: in iteration 1, |d| / sqrt(sigma^2 + sigma_m^2) <= z_cut at the prior mean member
# (m000), written to obs_kept.csv (every index point, its z there and whether kept); later iterations read it.
function gross_error_cut(C, k, obs, A, res, σm, zcut)
    f = joinpath(C, "obs_kept.csv")
    if k > 1
        lines = readlines(f)
        parse(Int, split(lines[1], "=")[2]) == length(obs) || error("obs_kept.csv is for another list of samples")
        return [parse(Int, split(l, ",")[1]) for l in lines[3:end] if split(l, ",")[end] == "1"]
    end
    haskey(res, 0) || error("iteration 1: the prior mean member m000 failed; it sets the gross-error cut")
    S0 = res[0].S; kept = Int[]
    open(f, "w") do io
        println(io, "# n_samples=", length(obs))
        println(io, "index,group,site,lat,lon,time,rsl,z_prior_mean,kept")
        for (i, o) in enumerate(obs)
            o.type == 0 || continue
            any(isnan, S0[i, :]) && continue
            mv = PM.tval(A, S0[i, :], o.t); d = mv[1] - o.rsl
            z = d / sqrt(PM.ip_sigma(o, d, mv[2])^2 + σm^2)
            ok = abs(z) <= zcut; ok && push!(kept, i)
            println(io, join((i, o.group, o.site_name, o.lat, o.lon, o.t, o.rsl, @sprintf("%.3f", z), Int(ok)), ","))
        end
    end
    og = unique(o.group for o in obs)
    logline(C, "\nGross-error cut at the prior mean member (|z| > ", zcut, " with sigma_m ", σm, " m): kept ", length(kept), " of ",
            count(o -> o.type == 0, obs), " index points (", join(["$g $(count(i -> obs[i].group == g, kept)) of $(count(o -> o.group == g && o.type == 0, obs))" for g in og], ", "), ")")
    kept
end

function update(C, k)
    C = abspath(C)
    cfg = config(C); es = cfg["esmda"]; α = Float64.(es["alpha"]); Na = length(α)
    sg = scale_groups(cfg["scale"]); P = params(cfg, sg)
    og = es["obs_groups"]; σm, Ls, Lt = obs_err(es)
    obs = [o for o in PM.read_indicators(cfg["score"]["fesmdata"]) if o.group in og]
    A = sort(NCDataset(ds -> Float64.(ds["time"][:]),
                       joinpath(cfg["score"]["fesmdata"], "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc")))
    sets = [("obs groups", collect(eachindex(obs))); [(g, findall(o -> o.group == g, obs)) for g in og]]
    ids, Θall = read_theta(C, k)
    run = parse.(Int, readlines(joinpath(itdir(C, k), "members.txt")))

    res = Dict{Int,Any}(); failed = String[]
    for m in run
        r, why = member_series(mdir(C, k, m), obs, A)
        r === nothing ? push!(failed, "m$(m): $(why)") : (res[m] = r)
        println("member ", m, r === nothing ? ": failed, $why" : ": ok")
    end
    ens = [m for m in run if m != 0 && haskey(res, m)]
    analysis = k <= Na
    if isempty(res)
        logline(C, "\n## iteration $k: no member succeeded; stopped\n\nfailed: ", join(failed, "; "))
        write(joinpath(C, "next"), "done\n"); error("iteration $k: no member succeeded")
    end
    if analysis && length(ens) < max(3, ceil(Int, es["min_members"] * count(!=(0), run)))
        logline(C, "\n## iteration $k: only $(length(ens)) members succeeded; stopped\n\nfailed: ", join(failed, "; "))
        write(joinpath(C, "next"), "done\n"); error("iteration $k: too few members")
    end

    # the gross-error cut, fixed in iteration 1 at the prior mean member
    kept = gross_error_cut(C, k, obs, A, res, σm, es["z_cut"])
    # the assimilated index points: kept, finite in every member used; their sigma at the mean RSL series of the ensemble (of the
    # mean member in an evaluation)
    used = isempty(ens) ? collect(keys(res)) : ens
    ip = [i for i in kept if all(m -> !any(isnan, res[m].S[i, :]), used)]
    Sbar = mean(res[m].S for m in used)
    σ = [(mv = PM.tval(A, Sbar[i, :], obs[i].t); PM.ip_sigma(obs[i], mv[1] - obs[i].rsl, mv[2])) for i in ip]
    bad = findall(s -> !(s^2 + σm^2 > 0), σ); deleteat!(ip, bad); deleteat!(σ, bad)
    y = [obs[i].rsl for i in ip]
    Lc = cholesky(Symmetric(obs_cov(obs[ip], σ, σm, Ls, Lt))).L      # R = Lc Lc'
    pred(m) = [PM.tval(A, res[m].S[i, :], obs[i].t)[1] for i in ip]
    chi2(m) = isempty(ip) ? NaN : sum(abs2, Lc \ (pred(m) .- y)) / length(ip)

    # scores
    It = itdir(C, k)
    summ = Dict(m => summaries(obs, sets, res[m].S, A) for m in keys(res))
    open(joinpath(It, "scores.csv"), "w") do io
        println(io, "member,set,n_ip,rms_z,median_abs_z,in2,mean_d,n_ml,ml_violated,n_tl,tl_violated")
        for m in sort(collect(keys(res))), (j, (sn, _)) in enumerate(sets)
            s = summ[m][j]
            println(io, join((m, sn, s.nip, s.rmsz, s.medz, s.in2, s.bias, s.nml, s.vml, s.ntl, s.vtl), ","))
        end
    end
    open(joinpath(It, "members.csv"), "w") do io
        println(io, "member,ok,chi2,esl_20ka,seconds")
        for m in run
            haskey(res, m) ? println(io, join((m, 1, chi2(m), res[m].esl, res[m].sec), ",")) : println(io, m, ",0,NaN,NaN,NaN")
        end
    end

    # report
    logline(C, "\n## iteration $k", analysis ? " (alpha = $(α[k]))" : " (evaluation)", "\n")
    logline(C, "members run ", length(run), ", succeeded ", length(res), isempty(failed) ? "" : "; failed: " * join(failed, "; "),
            ". Assimilated index points ", length(ip), " (", join(["$g $(count(i -> obs[i].group == g, ip))" for g in og], ", "), ")",
            isempty(bad) ? "" : " ($(length(bad)) with sigma 0 left out)", " (sigma at the ensemble-mean RSL", σm > 0 ? "; model error $(σm) m" * (Ls > 0 ? ", correlated over $(Ls) km and $(Lt / 1000) kyr" : "") : "",
            "); scores over all index points. Run time per member ",
            isempty(res) ? "-" : @sprintf("%.0f-%.0f s", extrema(r.sec for r in values(res))...), ".\n")
    logline(C, "| set | n IP | mean member: RMS z | median abs(z) | abs(z) <= 2 | mean d [m] | ML viol. | TL viol. | ensemble RMS z: median [min, max] | ML viol. | TL viol. |")
    logline(C, "|---|---|---|---|---|---|---|---|---|---|---|")
    for (j, (sn, _)) in enumerate(sets)
        s0 = haskey(summ, 0) ? summ[0][j] : nothing
        e(f) = isempty(ens) ? "-" : (v = [getfield(summ[m][j], f) for m in ens]; @sprintf("%.2f [%.2f, %.2f]", median(v), extrema(v)...))
        logline(C, "| ", sn, " | ", summ[first(keys(summ))][j].nip, " | ",
                s0 === nothing ? "- | - | - | - | - | -" : join((fmt(s0.rmsz), fmt(s0.medz), fmt(s0.in2), fmt(s0.bias), fmt(s0.vml), fmt(s0.vtl)), " | "),
                " | ", e(:rmsz), " | ", e(:vml), " | ", e(:vtl), " |")
    end
    if !isempty(ens)
        c = [chi2(m) for m in ens]; esl = [res[m].esl for m in ens]
        logline(C, "\nchi2 / n (assimilated, d' R^-1 d / n with this iteration's R): mean member ", haskey(res, 0) ? fmt(chi2(0)) : "-", ", ensemble ",
                @sprintf("%.2f [%.2f, %.2f]", median(c), extrema(c)...), "; ESL at 20 ka [m]: mean member ",
                haskey(res, 0) ? @sprintf("%.1f", res[0].esl) : "-", ", ensemble ", @sprintf("%.1f [%.1f, %.1f]", median(esl), extrema(esl)...))
    elseif haskey(res, 0)
        logline(C, "\nESL at 20 ka [m]: mean member ", @sprintf("%.1f", res[0].esl))
    end

    if !analysis
        write(joinpath(C, "next"), "done\n"); println("iteration $k: evaluated"); return
    end

    # ES-MDA update of the members that succeeded, in units of sigma
    col = Dict(m => j for (j, m) in enumerate(ids))
    Θ = Θall[:, [col[m] for m in ens]]
    D = Lc \ reduce(hcat, [pred(m) for m in ens])                     # whitened: R = I
    Θa = esmda_step(Θ, D, Lc \ y, α[k], MersenneTwister(1000 * cfg["prior"]["seed"] + k))
    tmax = cfg["prior"]["theta_max"]; nclip = count(abs.(Θa) .> tmax)
    Θa = clamp.(Θa, -tmax, tmax)
    θm = vec(mean(Θa; dims = 2))
    nxt = k + 1
    run2 = nxt <= Na || es["final"] == "ensemble" ? [0; ens] : es["final"] == "mean" ? [0] : Int[]
    write_iteration(C, nxt, cfg, sg, P, [0; ens], hcat(θm, Θa), run2; final = nxt > Na)
    write(joinpath(C, "next"), isempty(run2) ? "done\n" : "$nxt\n")

    # theta before and after, as factors exp(theta)
    θb = vec(mean(Θ; dims = 2)); sb = vec(std(Θ; dims = 2)); sa = vec(std(Θa; dims = 2))
    logline(C, "\nUpdate -> iteration $nxt", isempty(run2) ? " (not run)" : "", ": ", nclip, " of ", length(Θa),
            " values clipped to |theta| <= ", tmax, "; the factor exp(mean theta) before -> after (sd of theta before -> after)\n")
    logline(C, "| group | ", join([@sprintf("%g ka", a) for a in sg.nodes_ka], " | "), " |")
    logline(C, "|---|", repeat("---|", nnodes(sg)))
    pidx = [findfirst(p -> p.g == g && k2 in p.ks, P) for g in 1:ngroups(sg), k2 in 1:nnodes(sg)]
    for (g, gn) in enumerate(sg.names)
        logline(C, "| ", gn, " | ", join([q === nothing ? "1 (pinned)" : @sprintf("%.2f -> %.2f (%.2f -> %.2f)", exp(θb[q]), exp(θm[q]), sb[q], sa[q])
                                          for q in pidx[g, :]], " | "), " |")
    end
    println("iteration $k: updated, next ", isempty(run2) ? "done" : nxt)
end

function main(args)
    args[1] == "init" && return init(args[2], args[3])
    args[1] == "update" && return update(args[2], parse(Int, args[3]))
    args[1] == "selftest" && return selftest()
    error("plastic_esmda.jl init <config> <folder> | update <folder> <k>")
end
abspath(PROGRAM_FILE) == @__FILE__() && main(ARGS)
