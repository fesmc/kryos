# The emulator of the Earth search (docs/dev/earth-search-design.md, wave A): from the scores of the 1-D Earths of
# earth-search-sample (members.csv: the parameters) under each imposed ice history (gia_refcorr.sh cases scored as references,
# plastic_ensemble.jl scores: scores.csv), per ice history and data group a Gaussian process of the score in the standardized
# parameters, then the history-matching cut of a dense sample of the parameter space.
#   julia --project=tools/levante tools/levante/plastic_earth_emulate.jl <members.csv> <scores.csv> <config.toml> <outdir> [npts=100000] [seed=2]
#   julia --project=tools/levante tools/levante/plastic_earth_emulate.jl wave-b <members A.csv> <members B.csv> <scores.csv> <config.toml> <basis.nc> <outdir> [npts=100000] [seed=3]
#       wave B (3-D Earths, earth-search-sample3d): see emulate_b
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
function read_members(f, want = [p.name for p in ES.EARTH1D_PARAMS])
    M, h = readdlm(f, ',', Any; header = true)
    names = String.(vec(h)[2:end]); want == names || error("$f: columns $names, not $want")
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

# per pair of the parameters idx the fraction of the sample not ruled out (nb x nb bins), a corner plot
function plot_nroi(Xs01, keep, pn, lo, hi, idx, file, title; nb = 20)
    d = length(idx); fig = Figure(size = (190 * d, 180 * d)); hm = nothing
    for (ii, i) in enumerate(idx), (jj, j) in enumerate(idx)
        jj < ii || continue
        ax = Axis(fig[ii-1, jj], xlabel = ii == d ? pn[j] : "", ylabel = jj == 1 ? pn[i] : "",
                  xticklabelsvisible = ii == d, yticklabelsvisible = jj == 1, xticks = LinearTicks(3), yticks = LinearTicks(4))
        F = zeros(nb, nb); N = zeros(nb, nb)
        for k in axes(Xs01, 2)
            a = clamp(floor(Int, Xs01[j, k] * nb) + 1, 1, nb); b = clamp(floor(Int, Xs01[i, k] * nb) + 1, 1, nb)
            N[a, b] += 1; keep[k] && (F[a, b] += 1)
        end
        hm = heatmap!(ax, range(lo[j], hi[j], length = nb + 1), range(lo[i], hi[i], length = nb + 1), F ./ max.(N, 1),
                      colormap = :viridis, colorrange = (0, 1))
    end
    Colorbar(fig[1, max(2, d - 2):d-1], hm, label = "fraction not ruled out", vertical = false, width = Relative(0.8), tellheight = false, valign = :center)
    Label(fig[0, :], title * @sprintf(" (%.2f %% of the sample)", 100mean(keep)), fontsize = 18)
    save(file, fig)
end

# --- the Latin hypercube of the dense sample, in [0, 1]^d -----------------------------------------------------------------------
lhs01(d, n, rng) = permutedims(hcat([(randperm(rng, n) .- rand(rng, n)) ./ n for _ in 1:d]...))

# the groups that judge, their cut (the larger of [ensemble] implausibility and the bar's best I_g), n and y_cut = log(1 + cut sqrt(2/n))
function cut_setup(sc, cfg)
    en = cfg["ensemble"]; ug = String.(get(en, "use_groups", judging_groups(cfg)))
    cut0 = get(en, "implausibility", 3.0); bar = String.(get(en, "bar", String[]))
    cut = Dict(g => max(cut0, minimum([sc[(b, g)].I for b in bar if haskey(sc, (b, g))]; init = Inf)) for g in ug)
    for g in ug; isfinite(cut[g]) || (cut[g] = cut0); end
    ngrp = Dict(g => first([v.n for ((s, gg), v) in sc if gg == g]) for g in ug if any(k -> k[2] == g, keys(sc)))
    (ug = ug, cut0 = cut0, bar = bar, cut = cut, ngrp = ngrp, ycut = Dict(g => log(1 + cut[g] * sqrt(2 / ngrp[g])) for g in keys(ngrp)))
end
# the series of scores.csv that match rx ("ice-<ICE> | <earth>-E<nnn>"): ice => [(member, series)]
function cases_by(sc, rx, mem, membersf)
    cases = Dict{String,Vector{Tuple{String,String}}}()
    for s in unique(first.(collect(keys(sc))))
        m = match(rx, s); m === nothing && continue
        haskey(mem, m[2]) || error("$s: no member $(m[2]) in $membersf")
        push!(get!(cases, m[1], Tuple{String,String}[]), (String(m[2]), s))
    end
    cases
end
# a GP of y on the points X (d x n in [0, 1]) and its leave-one-out against the cut ycut: the emulator and the LOO statistics
function gp_fit(X, y, ycut, cutI)
    Xc = ColVecs(X); e, ys = emulator(Xc, y)
    μl, σl = loo(Xc, ys, e.h); μl = e.m .+ e.s .* μl; σl = e.s .* σl
    z = (y .- μl) ./ σl; out = y .> ycut; rule = (μl .- ycut) ./ σl .> cutI
    e, (n = length(y), h = e.h, rmse = sqrt(mean(abs2, y .- μl)), r2 = 1 - sum(abs2, y .- μl) / sum(abs2, y .- mean(y)),
        cov2 = mean(abs.(z) .< 2), pass = count(.!out), wrong = count(rule .& .!out), missed = count(.!rule .& out), ycut = ycut, μl = μl, σl = σl)
end

function emulate(membersf, scoresf, cfgf, outdir; npts = 100_000, seed = 2, cutI = 3.0)
    mkpath(outdir)
    mem, pn = read_members(membersf); P = ES.EARTH1D_PARAMS; d = length(P)
    lo = [p.lo for p in P]; hi = [p.hi for p in P]; std01(θ) = (θ .- lo) ./ (hi .- lo)
    sc = read_scores(scoresf); cfg = TOML.parsefile(cfgf)
    C = cut_setup(sc, cfg); ug = C.ug; cut = C.cut; cut0 = C.cut0; bar = C.bar; ngrp = C.ngrp; ycut = C.ycut
    series = unique(first.(collect(keys(sc))))
    cases = cases_by(sc, CASE_RX, mem, membersf)
    ices = sort(collect(keys(cases))); isempty(ices) || println("cases: ", join(["$i $(length(cases[i]))" for i in ices], ", "))
    isempty(ices) && error("no case \"ice-<ICE> | E1D-ES1-E<nnn>\" in $scoresf")
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
    gp = []; loorows = []; EMU = Dict{Tuple{String,String},Emu}()
    for ice in ices, g in ug
        haskey(ngrp, g) || continue
        ok = [(m, s) for (m, s) in cases[ice] if haskey(sc, (s, g)) && sc[(s, g)].chi2n > 0]
        length(ok) < 10 && (println("$ice $g: $(length(ok)) cases, no GP"); continue)
        y = [log(sc[(s, g)].chi2n) for (_, s) in ok]
        e, st = gp_fit(hcat([std01(mem[m]) for (m, _) in ok]...), y, ycut[g], cutI); μl, σl = st.μl, st.σl; r2 = st.r2
        push!(gp, (ice = ice, g = g, st...)); EMU[(ice, g)] = e
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
    # figures: per pair the fraction not ruled out
    for (nm, keep) in sets; plot_nroi(Xs01, keep, pn, lo, hi, 1:d, joinpath(outdir, "nroi_$nm.png"), "Earth search wave A: not ruled out, $nm"); end
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
    (EMU = EMU, ices = ices, C = C, sc = sc)
end

# --- wave B: 3-D Earths ------------------------------------------------------------------------------------------------------
# The cases "ice-<ICE> | E3D-ES2-E<nnn>" (the 3-D Earth E<nnn> of earth-search-sample3d, one pass with the r of the 1-D run of its
# radial profile) and "ice-<ICE> | E1D-ES2-E<nnn>" (that 1-D run). 100 runs are too few for a GP in 13 dimensions, so the score of a
# 3-D Earth is split: y_3D(θ) = y_1D(θ_r) + Δ(θ_l), the 1-D part from the wave-A emulator of the radial parameters θ_r, the
# 3-D effect Δ = y_3D - y_1D of the same radial profile (known at the design points: both were run) from a GP in the 8 lateral and
# regional parameters θ_l. Variances add. The split assumes the lateral effect does not depend on the radial profile within the
# space wave A left: the check is the prediction of y_3D at the design points (wave-A emulator + LOO of Δ) against the runs.
const CASE3_RX = r"^ice-(\S+) \| E3D-ES2-(E\d+)$"
const CASE1B_RX = r"^ice-(\S+) \| E1D-ES2-(E\d+)$"
# the distance of a 3-D Earth from the expert fields (the second axis of the design): radial, the RMS difference of log10 eta of P
# (clamped as the basis, [19.5, 24]) from each expert's radial mean (the basis' profile), weighted by volume, the nearest expert and
# Gowan's 1-D Earth; lateral, the distance of (α, c2, c3) from the experts' (1, coef_2, coef_3), the nearest; the RMS of the shifts
function expert_distance(θ, basis)
    NCDataset(basis) do b
        dep = Float64.(b["depth"][:]) ./ 1e3; prof = Float64.(b["profile"][:, :]); coef = Float64.(b["coef"][:, :])
        lo, hi = parse.(Float64, split(b.attrib["clamp"])); names = strip.(split(b.attrib["fields"], ";"))
        rr = ES.R_VILMA_KM .- dep; dr = [abs(i == 1 ? dep[2] - dep[1] : i == length(dep) ? dep[end] - dep[end-1] : (dep[i+1] - dep[i-1]) / 2) for i in eachindex(dep)]
        w = rr .^ 2 .* dr; w ./= sum(w)
        pf(t) = [clamp(ES.earth1d_profile(t, d), lo, hi) for d in dep]
        p = pf(θ); rad = [sqrt(sum(w .* (p .- prof[m, :]) .^ 2)) for m in axes(prof, 1)]
        g = log10(4e20); gl = log10(4e22); dg = sqrt(sum(w .* (p .- pf([120.0, g, g, gl, gl])) .^ 2))
        lat = [sqrt((θ[6] - 1)^2 + (θ[7] - coef[m, 2])^2 + (θ[8] - coef[m, 3])^2) for m in axes(coef, 1)]
        mr = argmin(rad); ml = argmin(lat)
        (radial = rad[mr], radial_to = first(split(names[mr])), gowan = dg, lateral = lat[ml], lateral_to = first(split(names[ml])),
         shifts = sqrt(mean(abs2, θ[9:13])))
    end
end
function emulate_b(membersA, membersB, scoresf, cfgf, basis, outdir; npts = 100_000, seed = 3, cutI = 3.0)
    mkpath(outdir)
    A = emulate(membersA, scoresf, cfgf, joinpath(outdir, "waveA"); npts = 20_000, cutI = cutI)    # the wave-A emulators, refitted
    P = ES.earth3d_params(basis); pn = [p.name for p in P]; d = length(P); lo = [p.lo for p in P]; hi = [p.hi for p in P]
    std01(θ) = (θ .- lo) ./ (hi .- lo); il = 6:d; ir = 1:5
    mem, _ = read_members(membersB, pn); sc = A.sc; C = A.C; ug = C.ug
    c3 = cases_by(sc, CASE3_RX, mem, membersB); c1 = cases_by(sc, CASE1B_RX, mem, membersB)
    ices = sort([i for i in keys(c3) if haskey(c1, i) && i in A.ices]); isempty(ices) && error("no 3-D cases with their 1-D run")
    println("wave B cases: ", join(["$i $(length(c3[i])) 3-D, $(length(c1[i])) 1-D" for i in ices], "; "))
    y(s, g) = haskey(sc, (s, g)) && sc[(s, g)].chi2n > 0 ? log(sc[(s, g)].chi2n) : NaN
    rng = MersenneTwister(seed); Xs01 = lhs01(d, npts, rng)
    Imax = Dict(ice => fill(-Inf, npts) for ice in ices); Imax1 = Dict(ice => fill(-Inf, npts) for ice in ices)
    gp = []; chk = []
    for ice in ices, g in ug
        haskey(A.EMU, (ice, g)) || continue
        s1 = Dict(c1[ice]); pairs = [(m, s, s1[m]) for (m, s) in c3[ice] if haskey(s1, m) && isfinite(y(s, g)) && isfinite(y(s1[m], g))]
        length(pairs) < 10 && (println("$ice $g: $(length(pairs)) pairs, no GP"); continue)
        Δ = [y(s, g) - y(s1m, g) for (_, s, s1m) in pairs]
        e, st = gp_fit(hcat([std01(mem[m])[il] for (m, _, _) in pairs]...), Δ, Inf, cutI)
        push!(gp, (ice = ice, g = g, st..., meanΔ = mean(Δ), sdΔ = std(Δ)))
        # the check: y_3D at the design points from the wave-A emulator and the LOO of Δ
        μa, σa = predict(A.EMU[(ice, g)], hcat([std01(mem[m])[ir] for (m, _, _) in pairs]...))
        y3 = [y(s, g) for (_, s, _) in pairs]; μ3 = μa .+ st.μl; σ3 = sqrt.(σa .^ 2 .+ st.σl .^ 2)
        push!(chk, (ice = ice, g = g, rmse = sqrt(mean(abs2, y3 .- μ3)), cov2 = mean(abs.(y3 .- μ3) ./ σ3 .< 2),
                    out = count(y3 .> C.ycut[g]), wrong = count(((μ3 .- C.ycut[g]) ./ σ3 .> cutI) .& (y3 .<= C.ycut[g]))))
        μr, σr = predict(A.EMU[(ice, g)], Xs01[ir, :]); μl, σl = predict(e, Xs01[il, :])
        Imax[ice] = max.(Imax[ice], (μr .+ μl .- C.ycut[g]) ./ sqrt.(σr .^ 2 .+ σl .^ 2))
        Imax1[ice] = max.(Imax1[ice], (μr .- C.ycut[g]) ./ σr)
        @printf("%s %-18s n %2d  Δ mean %+.3f sd %.3f  LOO RMSE %.3f R2 %.3f  ℓ %s  check: y3 RMSE %.3f within 2 sd %.2f\n", ice, g, st.n, mean(Δ), std(Δ),
                st.rmse, st.r2, join([@sprintf("%.2f", x) for x in st.h.ℓ], " "), chk[end].rmse, chk[end].cov2)
    end
    nroi = Dict(ice => Imax[ice] .<= cutI for ice in ices); radok = reduce((a, b) -> a .& b, [Imax1[i] .<= cutI for i in ices])
    sets = [(ice, nroi[ice]) for ice in ices]
    length(ices) > 1 && (push!(sets, ("both", reduce((a, b) -> a .& b, values(nroi)))); push!(sets, ("either", reduce((a, b) -> a .| b, values(nroi)))))
    Θs = lo .+ Xs01 .* (hi .- lo)
    for (nm, keep) in sets
        plot_nroi(Xs01, keep, pn, lo, hi, il, joinpath(outdir, "nroi_lateral_$nm.png"), "Earth search wave B: not ruled out, $nm")
        plot_nroi(Xs01, keep, pn, lo, hi, ir, joinpath(outdir, "nroi_radial_$nm.png"), "Earth search wave B: not ruled out, $nm")
    end
    # the runs: per ice by the groups passed, then log L; both: summed
    LLs(s) = sum(sc[(s, g)].logL for g in ug if haskey(sc, (s, g))); nps(s) = count(g -> haskey(sc, (s, g)) && sc[(s, g)].I <= C.cut[g], ug)
    fI(s, g) = haskey(sc, (s, g)) ? @sprintf("%.1f", sc[(s, g)].I) : "-"
    esc(s) = replace(s, "|" => "\\|")
    refs = [s for s in unique(first.(collect(keys(sc)))) if s in C.bar || s in ("ice-PMIST | E1D-PMIST", "ice-PMIST | E3D-B21-6", "ice-GLAC1D | E1D-PMIST", "ice-GLAC1D | E3D-B21-6")]
    fmtθ(m) = (t = mem[m]; @sprintf("%.0f km, %.2f/%.2f/%.2f/%.2f; α %.2f, c2 %+.2f, c3 %+.2f; s %+.1f/%+.1f/%+.1f/%+.1f/%+.1f", t...))
    open(joinpath(outdir, "cases.csv"), "w") do io
        println(io, "ice,member,", join(pn, ","), ",", join(["I_$g" for g in ug], ","), ",logL,passed,logL_1d,passed_1d")
        for ice in ices, (m, s) in sort(c3[ice])
            s1 = get(Dict(c1[ice]), m, "")
            println(io, ice, ",", m, ",", join(mem[m], ","), ",", join([haskey(sc, (s, g)) ? @sprintf("%.2f", sc[(s, g)].I) : "NaN" for g in ug], ","), ",",
                    @sprintf("%.1f", LLs(s)), ",", nps(s), ",", isempty(s1) ? "NaN" : @sprintf("%.1f", LLs(s1)), ",", isempty(s1) ? "NaN" : nps(s1))
        end
    end
    open(joinpath(outdir, "summary.md"), "w") do io
        println(io, "# Earth search, wave B: 3-D Earths\n\nCases: ", join(["$(length(c3[i])) 3-D and $(length(c1[i])) 1-D under ice-$i" for i in ices], ", "),
                ". Groups that judge and their cut as wave A (", joinpath(outdir, "waveA", "summary.md"), "): ", join([@sprintf("%s %.2f", g, C.cut[g]) for g in ug], ", "), ".\n")
        println(io, "## The emulator\n\ny_3D = y_1D(radial; the wave-A GP) + Δ(lateral; a GP of Δ = y_3D - y_1D of the same profile in the 8 lateral and ",
                "regional parameters). Δ: its mean and sd over the runs, the LOO RMSE and R² of its GP. Check: y_3D of the runs against the wave-A ",
                "emulator plus the LOO of Δ (RMSE, within 2 sd, runs failing the cut, runs the prediction would rule out though they pass).\n")
        println(io, "| ice | group | n | Δ mean | Δ sd | LOO RMSE | R² | y_3D RMSE | within 2 sd | runs out | wrong | length scales (", join(pn[il], ", "), ") |\n|---|---|---|---|---|---|---|---|---|---|---|---|")
        for (r, c) in zip(gp, chk)
            println(io, "| ", r.ice, " | ", r.g, " | ", r.n, " | ", @sprintf("%+.3f | %.3f | %.3f | %.2f | %.3f | %.2f", r.meanΔ, r.sdΔ, r.rmse, r.r2, c.rmse, c.cov2),
                    " | ", c.out, " | ", c.wrong, " | ", join([@sprintf("%.2f", x) for x in r.h.ℓ], ", "), " |")
        end
        println(io, "\n## Not ruled out\n\nOf $npts points (Latin hypercube over the 13 ranges); radial only: the fraction the wave-A emulator alone keeps ",
                @sprintf("(both ice histories: %.4f); the fraction of those the 3-D emulator keeps too in brackets. Ranges: min-max (5-95 %%).\n", mean(radok)))
        println(io, "| set | fraction | ", join(pn, " | "), " |\n|---|---|", repeat("---|", d))
        for (nm, keep) in sets
            rg = any(keep) ? [(q = Θs[k, keep]; @sprintf("%.2f-%.2f (%.2f-%.2f)", minimum(q), maximum(q), quantile(q, 0.05), quantile(q, 0.95))) for k in 1:d] : fill("-", d)
            println(io, "| ", nm, " | ", @sprintf("%.4f (%.3f)", mean(keep), count(keep .& radok) / max(1, count(radok))), " | ", join(rg, " | "), " |")
        end
        println(io, "\nFigures: nroi_lateral_<set>.png, nroi_radial_<set>.png.\n\n## The runs\n\nRanked by the groups passed, then log L over the groups that judge; ",
                "the 1-D run of the same profile after it. Earth: lithosphere, log10 eta UM/TZ/LM1/LM2; α, c2, c3; shifts Fennoscandia/Barents/Laurentide/Greenland/WAIS.\n")
        hdr = "| series | Earth | log L | passed | log L 1-D | " * join(["I " * g for g in ug], " | ") * " |\n|---|---|---|---|---|" * repeat("---|", length(ug))
        for ice in ices
            println(io, "\n### ice-$ice\n\n", hdr)
            for s in refs; println(io, "| ", esc(s), " | | ", @sprintf("%.1f", LLs(s)), " | ", nps(s), "/", length(ug), " | | ", join([fI(s, g) for g in ug], " | "), " |"); end
            s1 = Dict(c1[ice])
            for (m, s) in first(sort(c3[ice], by = x -> (-nps(x[2]), -LLs(x[2]))), 10)
                println(io, "| ", esc(s), " | ", fmtθ(m), " | ", @sprintf("%.1f", LLs(s)), " | ", nps(s), "/", length(ug), " | ",
                        haskey(s1, m) ? @sprintf("%.1f", LLs(s1[m])) : "-", " | ", join([fI(s, g) for g in ug], " | "), " |")
            end
        end
        if length(ices) > 1
            ms = intersect([Set(first.(c3[i])) for i in ices]...); S3 = Dict(i => Dict(c3[i]) for i in ices)
            best = first(sort(collect(ms), by = m -> (-sum(nps(S3[i][m]) for i in ices), -sum(LLs(S3[i][m]) for i in ices))), 10)
            println(io, "\n### Both ice histories, and the distance from the experts\n\nRanked by the groups passed summed over the ice histories, then the ",
                    "summed log L. Distance: radial, the volume-weighted RMS of log10 eta of the profile from the nearest expert's radial mean (and from ",
                    "Gowan's 1-D Earth) [dex]; lateral, the distance of (α, c2, c3) from the nearest expert's (1, c2, c3); shifts, their RMS [dex].\n")
            println(io, "| Earth | parameters | log L (", join(ices, " + "), ") | passed | radial (nearest) | from Gowan | lateral (nearest) | shifts |\n|---|---|---|---|---|---|---|---|")
            for m in best
                x = expert_distance(mem[m], basis)
                println(io, "| ", m, " | ", fmtθ(m), " | ", @sprintf("%.1f", sum(LLs(S3[i][m]) for i in ices)), " (", join([@sprintf("%.1f", LLs(S3[i][m])) for i in ices], " + "), ") | ",
                        join([string(nps(S3[i][m])) for i in ices], " + "), " | ", @sprintf("%.2f (%s) | %.2f | %.2f (%s) | %.2f", x.radial, x.radial_to, x.gowan, x.lateral, x.lateral_to, x.shifts), " |")
            end
        end
    end
    println("wrote ", joinpath(outdir, "summary.md"), "; not ruled out: ", join([@sprintf("%s %.4f", nm, mean(k)) for (nm, k) in sets], ", "))
end

function emulate_main(args)
    if args[1] == "wave-b"
        kv = Dict(split(x, "=") for x in args[8:end])
        return emulate_b(args[2:7]...; npts = parse(Int, get(kv, "npts", "100000")), seed = parse(Int, get(kv, "seed", "3")))
    end
    length(args) >= 4 || error("plastic_earth_emulate.jl <members.csv> <scores.csv> <config.toml> <outdir> [npts=100000] [seed=2]")
    kv = Dict(split(x, "=") for x in args[5:end])
    emulate(args[1], args[2], args[3], args[4]; npts = parse(Int, get(kv, "npts", "100000")), seed = parse(Int, get(kv, "seed", "2")))
end
abspath(PROGRAM_FILE) == @__FILE__() && emulate_main(ARGS)
