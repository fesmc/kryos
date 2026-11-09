# ES-MDA (Emerick and Reynolds 2013) calibration of tau_b of the plastic ice backend against the GAPSLIP 2.1 sea-level index points,
# with VILMA2 in the loop (docs/dev/plastic-backend.md, "Calibration"; the jobs: plastic_esmda.sh):
#   julia --project=<tools> plastic_esmda.jl init <config.toml> <calibration folder>   prior ensemble = iteration 1
#   julia --project=<tools> plastic_esmda.jl update <calibration folder> <k>            score iteration k, update, write k + 1
#   julia --project=<tools> plastic_esmda.jl selftest                                   the update against a linear-Gaussian posterior
# Parameters theta: the log multipliers of tau_b per scale group and node (plastic_taub_scale.jl), group-major; nodes of a group may be
# tied to one parameter (`tie` of [[scale.groups]], lists of node ages: the nodes without direct data); prior N(0, C), C of a group
# sigma^2 exp(-|t_p - t_q| / L) with t_p the mean age of a parameter's nodes, groups independent; every member clipped to
# |theta| <= theta_max after an update. Observations: the index points of pm_compare.jl's GROUPS in obs_groups,
# the model RSL at a point the scorer's (series_run: the run's VILMA2 stream as vilma_out.nc by pm_kryos_out.jl, bilinear, linear in
# time), so d = model - rsl is the d of the scores; sigma the scorer's ip_sigma (elevation on the side of d, and the age through the
# model slope) at the ensemble-mean RSL of the iteration, fixed within the iteration, plus sigma_add in quadrature (the model error the
# indicator lacks: GIA and Earth model, resolution, RSL interpolation). A gross-error cut, fixed in iteration 1 at the prior mean member
# (m000, Gowan's tau_b): index points with |d| / sqrt(sigma^2 + sigma_add^2) > z_cut there are not assimilated (obs_kept.csv); every
# score is over all points. Perturbed
# observations, update in the ensemble subspace (Woodbury: an Ne x Ne solve). Limiting data are not assimilated: their violation
# fractions are diagnostics.
# Folder: config.toml, kryos.x, input/, par.nml, groups.csv, obs.csv, obs_kept.csv, esmda.md (the report), next (the iteration to run, or done);
# iter_kk/: theta.csv (all members of the iteration), members.txt (those to run), mNNN/ (taub_scale.nc, par.nml, the run, status:
# exit code and seconds), scores.csv (per member and set), members.csv (per member: ok, chi2, ESL at 20 ka, run time).
# Member 0 is the mean member: the prior mean (0, Gowan's tau_b) in iteration 1, the ensemble mean after an update.
using LinearAlgebra, Printf, Random, Statistics, TOML, NCDatasets

module PM
include("pm_compare.jl")
end
module KO
include("pm_kryos_out.jl")
end
include("plastic_taub_scale.jl")

const ESL_AGE = 20000.0

itdir(C, k) = joinpath(C, @sprintf("iter_%02d", k))
mdir(C, k, m) = joinpath(itdir(C, k), @sprintf("m%03d", m))
config(C) = TOML.parsefile(joinpath(C, "config.toml"))

# the parameters: (group, its node indices, name, mean age [ka]), group-major, nodes in the order of nodes_ka; `tie` joins nodes
function params(cfg, sg)
    P = NamedTuple{(:g, :ks, :name, :t),Tuple{Int,Vector{Int},String,Float64}}[]
    for (g, spec) in enumerate(cfg["scale"]["groups"])
        ties = [Float64.(t) for t in get(spec, "tie", Any[])]
        seen = Set{Int}()
        for (k, a) in enumerate(sg.nodes_ka)
            k in seen && continue
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
# log_factor (group, node) of a parameter vector
function log_factor(sg, P, θ)
    lf = zeros(ngroups(sg), nnodes(sg))
    for (q, p) in enumerate(P), k in p.ks; lf[p.g, k] = θ[q]; end
    lf
end
logline(C, s...) = open(io -> println(io, s...), joinpath(C, "esmda.md"), "a")

# --- prior --------------------------------------------------------------------------------------------------------------------------
function prior_cov(P, pr)
    s = pr["sigma"]; L = pr["corr_time_kyr"]
    [p.g == q.g ? s^2 * exp(-abs(p.t - q.t) / L) : 0.0 for p in P, q in P]
end

# --- members ------------------------------------------------------------------------------------------------------------------------
# theta (P x n) of the members ids; run: the ids to run
function write_iteration(C, k, cfg, sg, P, ids, Θ, run)
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
        Base.run(`python3 $(cfg["run"]["set_nml"]) $(joinpath(D, "par.nml")) $("ice_nor.taub_scale_file=\"$f\"")`)
    end
    open(io -> foreach(m -> println(io, m), run), joinpath(It, "members.txt"), "w")
end

function read_theta(C, k)
    lines = readlines(joinpath(itdir(C, k), "theta.csv"))
    rows = [parse.(Float64, split(l, ",")) for l in lines[2:end]]
    Int.(first.(rows)), reduce(hcat, [r[2:end] for r in rows])
end

# --- init ---------------------------------------------------------------------------------------------------------------------------
function init(cfgfile, C)
    C = abspath(C)
    isfile(joinpath(C, "config.toml")) && error("$C is a calibration already")
    mkpath(C); cp(cfgfile, joinpath(C, "config.toml"))
    cfg = config(C); rn = cfg["run"]; K = rn["kryos_dir"]
    cp(joinpath(K, "libkryos", "bin", "kryos.x"), joinpath(C, "kryos.x"))
    mkpath(joinpath(C, "input"))
    for f in readdir(joinpath(K, "input")); cp(joinpath(K, "input", f), joinpath(C, "input", f); follow_symlinks = true); end
    cp(joinpath(K, rn["par"]), joinpath(C, "par.nml"))
    isempty(rn["set"]) || Base.run(`python3 $(rn["set_nml"]) $(joinpath(C, "par.nml")) $(rn["set"])`)
    open(joinpath(C, "slurm.env"), "w") do io
        for (k, v) in (("THREADS", rn["threads"]), ("MEM", rn["mem"]), ("TIME", rn["time"]), ("PART", rn["partition"]),
                       ("ACCT", rn["account"]), ("THROTTLE", rn["throttle"]), ("MAPS", rn["maps"]))
            println(io, k, "=", v)
        end
    end

    sg = scale_groups(cfg["scale"]); write_groups_csv(joinpath(C, "groups.csv"), sg)
    obs = PM.read_indicators(cfg["score"]["fesmdata"])
    og = cfg["esmda"]["obs_groups"]
    open(joinpath(C, "obs.csv"), "w") do io
        println(io, "group,n_index_points,n_marine_limiting,n_terrestrial_limiting")
        for g in og; println(io, g, ",", (count(o -> o.group == g && o.type == t, obs) for t in (0, -1, 1))...); end
    end

    es = cfg["esmda"]; ne = es["ne"]; P = params(cfg, sg)
    Cθ = prior_cov(P, cfg["prior"])
    Θ = clamp.(cholesky(Symmetric(Cθ)).L * randn(MersenneTwister(cfg["prior"]["seed"]), size(Cθ, 1), ne), -cfg["prior"]["theta_max"],
               cfg["prior"]["theta_max"])
    ids = collect(0:ne)
    write_iteration(C, 1, cfg, sg, P, ids, hcat(zeros(size(Cθ, 1)), Θ), ids)
    write(joinpath(C, "next"), "1\n")

    logline(C, "# ES-MDA calibration of tau_b (plastic backend)\n\nconfig: ", cfgfile, "; ", length(sg.domain_id), " domains, ",
            ngroups(sg), " groups x ", nnodes(sg), " nodes, ", size(Cθ, 1), " parameters (", join(pnames(P), ", "), "); prior sigma ",
            cfg["prior"]["sigma"], ", correlation time ", cfg["prior"]["corr_time_kyr"], " kyr, |theta| <= ", cfg["prior"]["theta_max"],
            "; Ne = ", ne, ", alpha = ", es["alpha"], ", final: ", es["final"], "; sigma_add ", es["sigma_add"], " m, gross-error cut |z| > ",
            es["z_cut"], " at the prior mean member\n")
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
    KO.kryos_vilma_out(D, D)
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

# a linear-Gaussian check: ES-MDA with alpha (9.333, 7, 4, 2) on a large ensemble against the exact posterior
function selftest()
    rng = MersenneTwister(7); P, m, N = 6, 40, 4000
    H = randn(rng, m, P); C0 = Matrix(Diagonal(fill(0.09, P))); θt = 0.3 .* randn(rng, P)
    y = H * θt .+ randn(rng, m)
    Cp = inv(inv(C0) + H' * H); μp = Cp * (H' * y)
    Θ = sqrt(0.09) .* randn(rng, P, N)
    for α in (9.333, 7.0, 4.0, 2.0); Θ = esmda_step(Θ, H * Θ, y, α, rng); end
    em = maximum(abs.(vec(mean(Θ; dims = 2)) .- μp) ./ sqrt.(diag(Cp))); es = maximum(abs.(sqrt.(diag(cov(Θ; dims = 2))) ./ sqrt.(diag(Cp)) .- 1))
    @printf("selftest: posterior mean error %.3f sd, sd ratio error %.3f (N = %d)\n", em, es, N)
    em < 0.1 && es < 0.1 || error("selftest failed")
end

# The index points kept for the assimilation: in iteration 1, |d| / sqrt(sigma^2 + sigma_add^2) <= z_cut at the prior mean member
# (m000), written to obs_kept.csv (every index point, its z there and whether kept); later iterations read it.
function gross_error_cut(C, k, obs, A, res, sadd, zcut)
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
            z = d / sqrt(PM.ip_sigma(o, d, mv[2])^2 + sadd^2)
            ok = abs(z) <= zcut; ok && push!(kept, i)
            println(io, join((i, o.group, o.site_name, o.lat, o.lon, o.t, o.rsl, @sprintf("%.3f", z), Int(ok)), ","))
        end
    end
    og = unique(o.group for o in obs)
    logline(C, "\nGross-error cut at the prior mean member (|z| > ", zcut, " with sigma_add ", sadd, " m): kept ", length(kept), " of ",
            count(o -> o.type == 0, obs), " index points (", join(["$g $(count(i -> obs[i].group == g, kept)) of $(count(o -> o.group == g && o.type == 0, obs))" for g in og], ", "), ")")
    kept
end

function update(C, k)
    C = abspath(C)
    cfg = config(C); es = cfg["esmda"]; α = Float64.(es["alpha"]); Na = length(α)
    sg = scale_groups(cfg["scale"]); P = params(cfg, sg)
    og = es["obs_groups"]; sadd = es["sigma_add"]
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
    kept = gross_error_cut(C, k, obs, A, res, sadd, es["z_cut"])
    # the assimilated index points: kept, finite in every member used; their sigma at the mean RSL series of the ensemble (of the
    # mean member in an evaluation)
    used = isempty(ens) ? collect(keys(res)) : ens
    ip = [i for i in kept if all(m -> !any(isnan, res[m].S[i, :]), used)]
    Sbar = mean(res[m].S for m in used)
    σ = [(mv = PM.tval(A, Sbar[i, :], obs[i].t); sqrt(PM.ip_sigma(obs[i], mv[1] - obs[i].rsl, mv[2])^2 + sadd^2)) for i in ip]
    bad = findall(s -> !(s > 0), σ); deleteat!(ip, bad); deleteat!(σ, bad)
    y = [obs[i].rsl for i in ip]
    pred(m) = [PM.tval(A, res[m].S[i, :], obs[i].t)[1] for i in ip]
    chi2(m) = isempty(ip) ? NaN : mean(((pred(m) .- y) ./ σ) .^ 2)

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
            isempty(bad) ? "" : " ($(length(bad)) with sigma 0 left out)", " (sigma at the ensemble-mean RSL", sadd > 0 ? ", + $(sadd) m" : "",
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
        logline(C, "\nchi2 / n (assimilated, this iteration's sigma): mean member ", haskey(res, 0) ? fmt(chi2(0)) : "-", ", ensemble ",
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
    D = reduce(hcat, [pred(m) for m in ens]) ./ σ
    Θa = esmda_step(Θ, D, y ./ σ, α[k], MersenneTwister(1000 * cfg["prior"]["seed"] + k))
    tmax = cfg["prior"]["theta_max"]; nclip = count(abs.(Θa) .> tmax)
    Θa = clamp.(Θa, -tmax, tmax)
    θm = vec(mean(Θa; dims = 2))
    nxt = k + 1
    run2 = nxt <= Na || es["final"] == "ensemble" ? [0; ens] : es["final"] == "mean" ? [0] : Int[]
    write_iteration(C, nxt, cfg, sg, P, [0; ens], hcat(θm, Θa), run2)
    write(joinpath(C, "next"), isempty(run2) ? "done\n" : "$nxt\n")

    # theta before and after, as factors exp(theta)
    θb = vec(mean(Θ; dims = 2)); sb = vec(std(Θ; dims = 2)); sa = vec(std(Θa; dims = 2))
    logline(C, "\nUpdate -> iteration $nxt", isempty(run2) ? " (not run)" : "", ": ", nclip, " of ", length(Θa),
            " values clipped to |theta| <= ", tmax, "; the factor exp(mean theta) before -> after (sd of theta before -> after)\n")
    logline(C, "| group | ", join([@sprintf("%g ka", a) for a in sg.nodes_ka], " | "), " |")
    logline(C, "|---|", repeat("---|", nnodes(sg)))
    pidx = [findfirst(p -> p.g == g && k2 in p.ks, P) for g in 1:ngroups(sg), k2 in 1:nnodes(sg)]
    for (g, gn) in enumerate(sg.names)
        logline(C, "| ", gn, " | ", join([@sprintf("%.2f -> %.2f (%.2f -> %.2f)", exp(θb[q]), exp(θm[q]), sb[q], sa[q]) for q in pidx[g, :]], " | "), " |")
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
