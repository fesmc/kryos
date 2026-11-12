# The emulator of the Earth search (docs/dev/earth-search-design.md, wave A): from the scores of the 1-D Earths of
# earth-search-sample (members.csv: the parameters) under each imposed ice history (gia_refcorr.sh cases scored as references,
# plastic_ensemble.jl scores: scores.csv), per ice history and data group a Gaussian process of the score in the standardized
# parameters, then the history-matching cut of a dense sample of the parameter space.
#   julia --project=tools/levante tools/levante/plastic_earth_emulate.jl <members.csv> <scores.csv> <config.toml> <outdir> [npts=100000] [seed=2]
# The cases: the references named "ice-<ICE> | E1D-ES1-E<nnn>" (member E<nnn> of members.csv under the ice <ICE>), scored on the
# judging window (scores.csv window = "window": the index points of age <= compare_max_age_ka). The groups that judge and the cut
# of the config as plastic_ensemble.jl weights: [ensemble] use_groups (default all but [[gnss]]/[esl] judge = false); cut_g the
# larger of [ensemble] implausibility (3) and the lowest I_g of the [ensemble] bar references in scores.csv.
# Emulated: y = log(chi2/n) per group, not I_g = (chi2 - n)/sqrt(2n). Both are monotone in chi2, so the cut is the same
# (y_cut = log(1 + cut sqrt(2/n))), but I_g spans orders of magnitude over the space (a few near the cut, thousands for a stiff
# lithosphere on a soft mantle): a stationary GP of I_g follows the worst Earths and smooths the region of the cut away; log(chi2/n)
# varies on one scale everywhere.
# GP: AbstractGPs.jl, constant mean (the sample mean of y, y standardized), Matérn 5/2 with one length scale per parameter (each
# parameter scaled to [0, 1] over its range), signal variance and a nugget; the hyperparameters maximize the marginal likelihood
# by a coordinate search on log grids (no optimizer). Cross-validation: leave-one-out, in closed form.
# Cut (history matching): a point is ruled out for an ice history when in a group that judges the emulated y exceeds y_cut by
# more than 3 emulator sd, I_emu = (E[y] - y_cut) / sd > 3. Not ruled out: per ice history, for both (the Earths the data allow
# under either ice: an Earth pattern counts when it holds for both) and for either.
# Writes in outdir: summary.md, gp.csv (hyperparameters and LOO per GP), loo.csv, cases.csv (per case and group the score),
# sample.csv (the dense sample: parameters and per ice the max I_emu over the judging groups) and the figures nroi_<set>.png
# (per parameter pair the fraction of the space not ruled out).
using AbstractGPs, KernelFunctions, CairoMakie, DelimitedFiles, LinearAlgebra, Printf, Random, Statistics, TOML
include("plastic_ensemble.jl")      # judging_groups; ES.EARTH1D_PARAMS

const CASE_RX = r"^ice-(\S+) \| E1D-ES1-(E\d+)$"

# --- data --------------------------------------------------------------------------------------------------------------------
function read_members(f)
    M, h = readdlm(f, ',', Any; header = true)
    names = String.(vec(h)[2:end]); [p.name for p in ES.EARTH1D_PARAMS] == names || error("$f: columns $names")
    Dict(String(M[i, 1]) => Float64.(M[i, 2:end]) for i in axes(M, 1)), names
end
# the window rows of scores.csv: (series, group) => (n, chi2n, I, logL)
function read_scores(f)
    M, h = readdlm(f, ',', Any; header = true); h = String.(vec(h)); c(n) = findfirst(==(n), h)
    out = Dict{Tuple{String,String},NamedTuple}()
    for i in axes(M, 1)
        M[i, c("window")] == "window" || continue
        g = String(M[i, c("group")]); g in ("all", "all points") && continue
        out[(String(M[i, c("series")]), g)] = (n = Int(M[i, c("n")]), chi2n = Float64(M[i, c("chi2n")]), I = Float64(M[i, c("I")]),
                                               logL = Float64(M[i, c("logL")]))
    end
    out
end

# --- GP ----------------------------------------------------------------------------------------------------------------------
kernel(ℓ, sf2) = sf2 * (Matern52Kernel() ∘ ARDTransform(1 ./ ℓ))
nlml(X, y, ℓ, sf2, sn2) = -logpdf(GP(kernel(ℓ, sf2))(X, sn2), y)
const GRID_L = exp.(range(log(0.03), log(50.0), length = 40))
const GRID_SF2 = exp.(range(log(0.05), log(100.0), length = 34))
const GRID_SN2 = exp.(range(log(1e-6), log(0.3), length = 28))
# maximum likelihood by coordinate search: each hyperparameter in turn over its grid, the others fixed, until a sweep changes none
function fit_gp(X, y)
    d = size(X.X, 1); ℓ = fill(0.5, d); sf2 = 1.0; sn2 = 1e-3; best = nlml(X, y, ℓ, sf2, sn2)
    for sweep in 1:20
        changed = false
        for k in 1:d+2
            grid = k <= d ? GRID_L : k == d + 1 ? GRID_SF2 : GRID_SN2
            for v in grid
                ℓt = copy(ℓ); sft = sf2; snt = sn2
                k <= d ? (ℓt[k] = v) : k == d + 1 ? (sft = v) : (snt = v)
                L = try nlml(X, y, ℓt, sft, snt) catch; Inf end
                L < best - 1e-9 && (best = L; ℓ = ℓt; sf2 = sft; sn2 = snt; changed = true)
            end
        end
        changed || break
    end
    (ℓ = ℓ, sf2 = sf2, sn2 = sn2, nlml = best)
end
# leave-one-out mean and sd (closed form: Rasmussen and Williams 2006, eq. 5.12)
function loo(X, y, h)
    K = kernelmatrix(kernel(h.ℓ, h.sf2), X) + h.sn2 * I
    Ki = inv(cholesky(Symmetric(K))); α = Ki * y
    y .- α ./ diag(Ki), sqrt.(1 ./ diag(Ki))
end
# an emulator of y over the parameters: the GP on the standardized y
struct Emu
    h::NamedTuple; post::Any; m::Float64; s::Float64
end
function emulator(X, y)
    m = mean(y); s = std(y); ys = (y .- m) ./ s; h = fit_gp(X, ys)
    Emu(h, posterior(GP(kernel(h.ℓ, h.sf2))(X, h.sn2), ys), m, s), ys
end
function predict(e::Emu, Xs; chunk = 10_000)
    μ = zeros(size(Xs, 2)); σ = similar(μ)
    for i0 in 1:chunk:size(Xs, 2)
        r = i0:min(i0 + chunk - 1, size(Xs, 2)); a, v = mean_and_var(e.post(ColVecs(Xs[:, r])))
        μ[r] = e.m .+ e.s .* a; σ[r] = e.s .* sqrt.(max.(v, 0))
    end
    μ, σ
end

# --- the Latin hypercube of the dense sample, in [0, 1]^d -----------------------------------------------------------------------
lhs01(d, n, rng) = permutedims(hcat([(randperm(rng, n) .- rand(rng, n)) ./ n for _ in 1:d]...))

function emulate(membersf, scoresf, cfgf, outdir; npts = 100_000, seed = 2, cutI = 3.0)
    mkpath(outdir)
    mem, pn = read_members(membersf); P = ES.EARTH1D_PARAMS; d = length(P)
    lo = [p.lo for p in P]; hi = [p.hi for p in P]; std01(θ) = (θ .- lo) ./ (hi .- lo)
    sc = read_scores(scoresf); cfg = TOML.parsefile(cfgf); en = cfg["ensemble"]
    ug = String.(get(en, "use_groups", judging_groups(cfg)))
    # the cut per group: the larger of the plain cut and the bar's best I_g
    cut0 = get(en, "implausibility", 3.0); bar = String.(get(en, "bar", String[]))
    cut = Dict(g => max(cut0, minimum([sc[(b, g)].I for b in bar if haskey(sc, (b, g))]; init = Inf)) for g in ug)
    for g in ug; isfinite(cut[g]) || (cut[g] = cut0); end
    # the cases per ice
    series = unique(first.(collect(keys(sc))))
    cases = Dict{String,Vector{Tuple{String,String}}}()       # ice => [(member, series)]
    for s in series
        m = match(CASE_RX, s); m === nothing && continue
        haskey(mem, m[2]) || error("$s: no member $(m[2]) in $membersf")
        push!(get!(cases, m[1], Tuple{String,String}[]), (String(m[2]), s))
    end
    ices = sort(collect(keys(cases))); isempty(ices) || println("cases: ", join(["$i $(length(cases[i]))" for i in ices], ", "))
    isempty(ices) && error("no case \"ice-<ICE> | E1D-ES1-E<nnn>\" in $scoresf")
    ngrp = Dict(g => first([v.n for ((s, gg), v) in sc if gg == g]) for g in ug if any(k -> k[2] == g, keys(sc)))
    ycut = Dict(g => log(1 + cut[g] * sqrt(2 / ngrp[g])) for g in keys(ngrp))
    # the per-case table
    open(joinpath(outdir, "cases.csv"), "w") do io
        println(io, "ice,member,", join(pn, ","), ",", join(["I_$g" for g in ug], ","), ",logL,passed")
        for ice in ices, (m, s) in sort(cases[ice])
            v = [get(sc, (s, g), nothing) for g in ug]
            println(io, ice, ",", m, ",", join(mem[m], ","), ",", join([x === nothing ? "NaN" : @sprintf("%.2f", x.I) for x in v], ","), ",",
                    @sprintf("%.1f", sum(x.logL for x in v if x !== nothing)), ",", count(k -> v[k] !== nothing && v[k].I <= cut[ug[k]], eachindex(ug)))
        end
    end
    # the GPs and their cross-validation
    rng = MersenneTwister(seed); Xs01 = lhs01(d, npts, rng)
    Imax = Dict(ice => fill(-Inf, npts) for ice in ices)           # per point the max I_emu over the judging groups
    Igrp = Dict((ice, g) => zeros(npts) for ice in ices for g in ug)
    gp = []; loorows = []
    for ice in ices, g in ug
        haskey(ngrp, g) || continue
        ok = [(m, s) for (m, s) in cases[ice] if haskey(sc, (s, g)) && sc[(s, g)].chi2n > 0]
        length(ok) < 10 && (println("$ice $g: $(length(ok)) cases, no GP"); continue)
        X = ColVecs(hcat([std01(mem[m]) for (m, _) in ok]...)); y = [log(sc[(s, g)].chi2n) for (_, s) in ok]
        e, ys = emulator(X, y)
        μl, σl = loo(X, ys, e.h); μl = e.m .+ e.s .* μl; σl = e.s .* σl
        z = (y .- μl) ./ σl; out = y .> ycut[g]; rule = (μl .- ycut[g]) ./ σl .> cutI
        r2 = 1 - sum(abs2, y .- μl) / sum(abs2, y .- mean(y))
        push!(gp, (ice = ice, g = g, n = length(y), h = e.h, rmse = sqrt(mean(abs2, y .- μl)), r2 = r2, cov2 = mean(abs.(z) .< 2),
                   pass = count(.!out), wrong = count(rule .& .!out), missed = count(.!rule .& out), ycut = ycut[g]))
        for (k, (m, _)) in enumerate(ok); push!(loorows, (ice, g, m, y[k], μl[k], σl[k])); end
        μ, σ = predict(e, Xs01); Ie = (μ .- ycut[g]) ./ σ
        Igrp[(ice, g)] = Ie; Imax[ice] = max.(Imax[ice], Ie)
        @printf("%s %-18s n %3d  LOO RMSE %.3f R2 %.3f |z|<2 %.2f  ℓ %s  sf2 %.3g sn2 %.2g  not ruled out %.3f\n", ice, g, length(y),
                gp[end].rmse, r2, gp[end].cov2, join([@sprintf("%.2f", x) for x in e.h.ℓ], " "), e.h.sf2, e.h.sn2, mean(Ie .<= cutI))
    end
    open(joinpath(outdir, "gp.csv"), "w") do io
        println(io, "ice,group,n,", join(["l_$p" for p in pn], ","), ",sf2,sn2,loo_rmse,loo_r2,loo_cov2,ycut,cases_pass,loo_wrong_out,loo_missed_out")
        for r in gp
            println(io, r.ice, ",", r.g, ",", r.n, ",", join([@sprintf("%.4f", x) for x in r.h.ℓ], ","), ",", @sprintf("%.4g,%.3g,%.4f,%.4f,%.3f,%.4f", r.h.sf2, r.h.sn2, r.rmse, r.r2, r.cov2, r.ycut),
                    ",", r.pass, ",", r.wrong, ",", r.missed)
        end
    end
    open(joinpath(outdir, "loo.csv"), "w") do io
        println(io, "ice,group,member,y,y_loo,sd_loo"); for r in loorows; println(io, join(r[1:3], ","), ",", join([@sprintf("%.4f", x) for x in r[4:6]], ",")); end
    end
    # the cut
    Θs = lo .+ Xs01 .* (hi .- lo)
    nroi = Dict(ice => Imax[ice] .<= cutI for ice in ices)
    sets = [(ice, nroi[ice]) for ice in ices]
    length(ices) > 1 && (push!(sets, ("both", reduce((a, b) -> a .& b, values(nroi)))); push!(sets, ("either", reduce((a, b) -> a .| b, values(nroi)))))
    open(joinpath(outdir, "sample.csv"), "w") do io
        println(io, join(pn, ","), ",", join(["Imax_$ice" for ice in ices], ","))
        for i in 1:npts; println(io, join([@sprintf("%.4f", x) for x in Θs[:, i]], ","), ",", join([@sprintf("%.2f", Imax[ice][i]) for ice in ices], ",")); end
    end
    # figures: per pair the fraction not ruled out (20 x 20 bins), and the cases that pass every judging group
    nb = 20
    for (nm, keep) in sets
        fig = Figure(size = (950, 900)); hm = nothing
        for i in 2:d, j in 1:i-1
            ax = Axis(fig[i-1, j], xlabel = i == d ? pn[j] : "", ylabel = j == 1 ? pn[i] : "",
                      xticklabelsvisible = i == d, yticklabelsvisible = j == 1, xticks = LinearTicks(3), yticks = LinearTicks(4))
            F = zeros(nb, nb); N = zeros(nb, nb)
            for k in 1:npts
                a = clamp(floor(Int, Xs01[j, k] * nb) + 1, 1, nb); b = clamp(floor(Int, Xs01[i, k] * nb) + 1, 1, nb)
                N[a, b] += 1; keep[k] && (F[a, b] += 1)
            end
            hm = heatmap!(ax, range(lo[j], hi[j], length = nb + 1), range(lo[i], hi[i], length = nb + 1), F ./ max.(N, 1),
                          colormap = :viridis, colorrange = (0, 1))
        end
        Colorbar(fig[1, 3:d-1], hm, label = "fraction not ruled out", vertical = false, width = Relative(0.8),
                 tellheight = false, valign = :center)
        Label(fig[0, :], "Earth search wave A: not ruled out, $nm ($(@sprintf("%.1f", 100mean(keep))) % of the space)", fontsize = 18)
        save(joinpath(outdir, "nroi_$nm.png"), fig)
    end
    # the best cases: per ice by the groups passed, then log L over the groups that judge; both: summed over the ices
    rows = Dict(ice => [(m = m, s = s, v = [get(sc, (s, g), nothing) for g in ug]) for (m, s) in cases[ice]] for ice in ices)
    LL(r) = sum(x.logL for x in r.v if x !== nothing); np(r) = count(k -> r.v[k] !== nothing && r.v[k].I <= cut[ug[k]], eachindex(ug))
    fmtθ(m) = @sprintf("%.0f km, %.2f / %.2f / %.2f / %.2f", mem[m]...)
    refrow(s) = (v = [get(sc, (s, g), nothing) for g in ug]; any(!isnothing, v) ? (m = "", s = s, v = v) : nothing)
    open(joinpath(outdir, "summary.md"), "w") do io
        println(io, "# Earth search, wave A: 1-D Earths and their emulator\n")
        println(io, "Cases: ", join(["$(length(cases[i])) under ice-$i" for i in ices], ", "), " (", membersf, "; ", scoresf, "). Groups that ",
                "judge: ", join(ug, ", "), ". Cut per group (the larger of $cut0 and the best I_g of ", join(bar, ", "), "): ",
                join([@sprintf("%s %.2f", g, cut[g]) for g in ug], ", "), ".\n")
        println(io, "## The emulators\n\nPer ice history and group a GP of y = log(chi2/n) (Matérn 5/2, a length scale per parameter on [0, 1] of its range); ",
                "leave-one-out: RMSE and R² of y, the fraction of cases within 2 LOO sd (0.95 if the sd is right), and the cases the LOO prediction would ",
                "rule out (I_emu > $cutI) though they pass the cut (wrong) or keep though they fail (kept).\n")
        println(io, "| ice | group | n | pass | y_cut | LOO RMSE | R² | within 2 sd | wrong | kept | length scales (", join(pn, ", "), ") |\n|---|---|---|---|---|---|---|---|---|---|---|")
        for r in gp
            println(io, "| ", r.ice, " | ", r.g, " | ", r.n, " | ", r.pass, " | ", @sprintf("%.3f | %.3f | %.3f | %.2f", r.ycut, r.rmse, r.r2, r.cov2), " | ",
                    r.wrong, " | ", r.missed, " | ", join([@sprintf("%.2f", x) for x in r.h.ℓ], ", "), " |")
        end
        println(io, "\n## Not ruled out\n\nOf $npts points (Latin hypercube over the ranges): the fraction where no judging group has I_emu > $cutI, and the ",
                "ranges of the parameters there (min-max; 5-95 %).\n")
        println(io, "| set | fraction | ", join(pn, " | "), " |\n|---|---|", repeat("---|", d))
        for (nm, keep) in sets
            if any(keep)
                rg = [(q = Θs[k, keep]; @sprintf("%.2f-%.2f (%.2f-%.2f)", minimum(q), maximum(q), quantile(q, 0.05), quantile(q, 0.95))) for k in 1:d]
            else
                rg = fill("-", d)
            end
            println(io, "| ", nm, " | ", @sprintf("%.4f", mean(keep)), " | ", join(rg, " | "), " |")
        end
        println(io, "\nRuled out by each group alone (fraction of the space with I_emu > $cutI there):\n")
        println(io, "| group | ", join(ices, " | "), " |\n|---|", repeat("---|", length(ices)))
        for g in ug; println(io, "| ", g, " | ", join([haskey(Igrp, (i, g)) ? @sprintf("%.3f", mean(Igrp[(i, g)] .> cutI)) : "-" for i in ices], " | "), " |"); end
        println(io, "\nFigures: ", join(["nroi_$nm.png" for (nm, _) in sets], ", "), " (per parameter pair the fraction of the sample not ruled out).\n")
        println(io, "## The best cases\n\nRanked by the groups passed (I_g <= cut), then log L over the groups that judge. Parameters: lithosphere, ",
                "log10 eta of UM / TZ / LM1 / LM2.\n")
        hdr = "| series | Earth | log L | passed | " * join(["I " * g for g in ug], " | ") * " |\n|---|---|---|---|" * repeat("---|", length(ug))
        row(r) = println(io, "| ", replace(r.s, "|" => "\\|"), " | ", isempty(r.m) ? "" : fmtθ(r.m), " | ", @sprintf("%.1f", LL(r)), " | ", np(r), "/", length(ug), " | ",
                         join([x === nothing ? "-" : @sprintf("%.1f", x.I) for x in r.v], " | "), " |")
        refs = filter(!isnothing, [refrow(s) for s in series if match(CASE_RX, s) === nothing && occursin("E1D-PMIST", s) || s in bar])
        for ice in ices
            println(io, "\n### ice-$ice\n\n", hdr)
            for r in refs; row(r); end
            for r in first(sort(rows[ice], by = r -> (-np(r), -LL(r))), 10); row(r); end
        end
        if length(ices) > 1
            both = [(m = m, rs = [only(filter(r -> r.m == m, rows[i])) for i in ices]) for m in intersect([Set(r.m for r in rows[i]) for i in ices]...)]
            println(io, "\n### Both ice histories\n\nRanked by the groups passed summed over the ice histories, then the summed log L.\n")
            println(io, "| Earth | parameters | log L (", join(ices, " + "), ") | passed (", join(ices, " + "), ") |\n|---|---|---|---|")
            for b in first(sort(both, by = b -> (-sum(np, b.rs), -sum(LL, b.rs))), 10)
                println(io, "| ", b.m, " | ", fmtθ(b.m), " | ", @sprintf("%.1f", sum(LL, b.rs)), " (", join([@sprintf("%.1f", LL(r)) for r in b.rs], " + "), ") | ",
                        join([string(np(r)) for r in b.rs], " + "), " of ", length(ug), " each |")
            end
        end
    end
    println("wrote ", joinpath(outdir, "summary.md"), "; not ruled out: ", join([@sprintf("%s %.4f", nm, mean(k)) for (nm, k) in sets], ", "))
end

function emulate_main(args)
    length(args) >= 4 || error("plastic_earth_emulate.jl <members.csv> <scores.csv> <config.toml> <outdir> [npts=100000] [seed=2]")
    kv = Dict(split(x, "=") for x in args[5:end])
    emulate(args[1], args[2], args[3], args[4]; npts = parse(Int, get(kv, "npts", "100000")), seed = parse(Int, get(kv, "seed", "2")))
end
abspath(PROGRAM_FILE) == @__FILE__() && emulate_main(ARGS)
