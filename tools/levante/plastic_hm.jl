# History matching (Vernon et al. 2010; Williamson et al. 2013) of kryos's simulations: per data group a Gaussian-process emulator of the
# score over the parameters of an ensemble and the cut of the space whose emulated score fails a judging group with confidence.
# The library of plastic_earth_emulate.jl (the Earth search) and plastic_joint_emulate.jl (the joint ice + Earth ensemble):
#   read_scores        the window rows of scores.csv (plastic_ensemble.jl scores)
#   fit_gp, loo, emulator, predict, gp_fit   the GP (below) and its leave-one-out
#   cut_setup          the judging groups of a scoring config, their cut and y_cut
#   cases_by           the series of scores.csv that match a regex, by set
#   lhs01, hm_fit      the dense sample, the GPs of a set of cases and the implausibility on the sample
#   write_gp, write_loo, plot_nroi   gp.csv, loo.csv, the fraction not ruled out per parameter pair
# Emulated: y = log(chi2/n) per group, not I_g = (chi2 - n)/sqrt(2n). Both are monotone in chi2, so the cut is the same
# (y_cut = log(1 + cut sqrt(2/n))), but I_g spans orders of magnitude over the space (a few near the cut, thousands for a stiff
# lithosphere on a soft mantle): a stationary GP of I_g follows the worst runs and smooths the region of the cut away; log(chi2/n)
# varies on one scale everywhere.
# GP: AbstractGPs.jl, constant mean (the sample mean of y, y standardized), Matérn 5/2 with one length scale per parameter (each
# parameter scaled to [0, 1] over its range), signal variance and a nugget; the hyperparameters maximize the marginal likelihood
# by a coordinate search on log grids (no optimizer). Cross-validation: leave-one-out, in closed form.
# Cut: a point is ruled out when in a group that judges the emulated y exceeds y_cut by more than cutI (3) emulator sd,
# I_emu = (E[y] - y_cut) / sd > cutI.
using AbstractGPs, KernelFunctions, CairoMakie, DelimitedFiles, LinearAlgebra, Printf, Random, Statistics, TOML
include("plastic_ensemble.jl")      # judging_groups, plain_cut_groups; ES

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
    Colorbar(fig[1:d-1, d], hm, label = "fraction not ruled out", height = Relative(0.6))
    Label(fig[0, :], title * @sprintf(" (%.2f %% of the sample)", 100mean(keep)), fontsize = 18)
    save(file, fig; px_per_unit = 1)          # about 100-200 kB: the docs figures
end

# --- the Latin hypercube of the dense sample, in [0, 1]^d -----------------------------------------------------------------------
lhs01(d, n, rng) = permutedims(hcat([(randperm(rng, n) .- rand(rng, n)) ./ n for _ in 1:d]...))

# the groups that judge, their cut (the larger of [ensemble] implausibility and the bar's best I_g), n and y_cut = log(1 + cut sqrt(2/n))
function cut_setup(sc, cfg)
    en = cfg["ensemble"]; ug = String.(get(en, "use_groups", judging_groups(cfg)))
    cut0 = get(en, "implausibility", 3.0); bar = String.(get(en, "bar", String[]))
    fixed = plain_cut_groups(cfg)                             # groups that keep the plain cut (fixed_cut and ESL; plastic_ensemble.jl)
    cut = Dict(g => g in fixed ? cut0 : max(cut0, minimum([sc[(b, g)].I for b in bar if haskey(sc, (b, g))]; init = Inf)) for g in ug)
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

# The GPs of a set of cases and the history-matching cut on a dense sample. cases: set (an ice history, or one set) =>
# [(member, series)]; θ: member => its parameters, std01 maps them to [0, 1]; sc: read_scores; C: cut_setup; Xs01: the dense
# sample (d x npts in [0, 1]). Per set and judging group with at least min_cases cases a GP of y = log(chi2/n) (gp_fit) and on the
# sample I_emu = (E[y] - y_cut) / sd: Igrp[(set, group)] and their max over the groups, Imax[set] (filled in). Returns the GP rows
# (gp.csv), the leave-one-out rows (loo.csv) and the emulators.
function hm_fit(cases, θ, std01, sc, C, Xs01, Imax, Igrp; cutI = 3.0, min_cases = 10)
    gp = []; loorows = []; EMU = Dict{Tuple{String,String},Emu}()
    for ice in sort(collect(keys(cases))), g in C.ug
        haskey(C.ngrp, g) || continue
        ok = [(m, s) for (m, s) in cases[ice] if haskey(sc, (s, g)) && sc[(s, g)].chi2n > 0]
        length(ok) < min_cases && (println("$ice $g: $(length(ok)) cases, no GP"); continue)
        y = [log(sc[(s, g)].chi2n) for (_, s) in ok]
        e, st = gp_fit(hcat([std01(θ[m]) for (m, _) in ok]...), y, C.ycut[g], cutI); μl, σl = st.μl, st.σl; r2 = st.r2
        push!(gp, (ice = ice, g = g, st...)); EMU[(ice, g)] = e
        for (k, (m, _)) in enumerate(ok); push!(loorows, (ice, g, m, y[k], μl[k], σl[k])); end
        μ, σ = predict(e, Xs01); Ie = (μ .- C.ycut[g]) ./ σ
        Igrp[(ice, g)] = Ie; Imax[ice] = max.(Imax[ice], Ie)
        @printf("%s %-18s n %3d  LOO RMSE %.3f R2 %.3f |z|<2 %.2f  ℓ %s  sf2 %.3g sn2 %.2g  not ruled out %.3f\n", ice, g, length(y),
                gp[end].rmse, r2, gp[end].cov2, join([@sprintf("%.2f", x) for x in e.h.ℓ], " "), e.h.sf2, e.h.sn2, mean(Ie .<= cutI))
    end
    gp, loorows, EMU
end
function write_gp(f, gp, pn)
    open(f, "w") do io
        println(io, "ice,group,n,", join(["l_$p" for p in pn], ","), ",sf2,sn2,loo_rmse,loo_r2,loo_cov2,ycut,cases_pass,loo_wrong_out,loo_missed_out")
        for r in gp
            println(io, r.ice, ",", r.g, ",", r.n, ",", join([@sprintf("%.4f", x) for x in r.h.ℓ], ","), ",", @sprintf("%.4g,%.3g,%.4f,%.4f,%.3f,%.4f", r.h.sf2, r.h.sn2, r.rmse, r.r2, r.cov2, r.ycut),
                    ",", r.pass, ",", r.wrong, ",", r.missed)
        end
    end
end
function write_loo(f, loorows)
    open(f, "w") do io
        println(io, "ice,group,member,y,y_loo,sd_loo"); for r in loorows; println(io, join(r[1:3], ","), ",", join([@sprintf("%.4f", x) for x in r[4:6]], ",")); end
    end
end
