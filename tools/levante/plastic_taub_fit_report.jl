# Stage 1 of the tau_b calibration (tests/plastic_taub_fit.f90; docs/dev/plastic-backend.md, "Calibration v2"): the L-curve and the
# choice of lambda per grid, the backend's taub_domain_factor_file of all grids, the present-day thickness per region before
# (Gowan's tau_b) and after, maps.
#   julia --project=tools/levante tools/levante/plastic_taub_fit_report.jl <out folder> <fit_<grid>.nc> ... [<grid>=<lambda> ...]
# lambda, unless given: by the discrepancy principle, the largest lambda with chi2 / n <= 1 (sigma_H the error of the fit); when no
# lambda reaches it, the corner of the L-curve (log chi2 against log roughness, lambda > 0: the largest Menger curvature of three
# consecutive points). Writes <out>/taub_domain_factor.nc (domain_id, log_factor of the free domains of every grid; the others keep
# factor 1), report.md, lcurve.png and maps_<grid>.png.
using CairoMakie, NCDatasets, Printf, Statistics
include(joinpath(@__DIR__, "plastic_pd_regions.jl"))

const FES = get(ENV, "FESMDATA", "/work/ba1442/fesmdata")   # tools/levante/env.sh
out = ARGS[1]; mkpath(out)
fits = [a for a in ARGS[2:end] if endswith(a, ".nc")]
given = Dict(split(a, "=")[1] => parse(Float64, split(a, "=")[2]) for a in ARGS[2:end] if occursin("=", a))
gridof(f) = (m = match(r"fit_(.+)\.nc$", basename(f)); m === nothing ? error("$f: not fit_<grid>.nc") : String(m[1]))

# Menger curvature of consecutive points of the L-curve (log10 chi2, log10 roughness), lambda > 0
function corner(chi2, rough, lam)
    k = findall(l -> l > 0, lam); x = log10.(chi2[k]); y = log10.(max.(rough[k], 1e-12))
    c = fill(-Inf, length(k))
    for i in 2:length(k)-1
        a = hypot(x[i] - x[i-1], y[i] - y[i-1]); b = hypot(x[i+1] - x[i], y[i+1] - y[i]); d = hypot(x[i+1] - x[i-1], y[i+1] - y[i-1])
        area2 = (x[i] - x[i-1]) * (y[i+1] - y[i-1]) - (x[i+1] - x[i-1]) * (y[i] - y[i-1])
        c[i] = 2abs(area2) / (a * b * d)
    end
    k[argmax(c)], Dict(zip(k, c))
end

io = open(joinpath(out, "report.md"), "w")
ids = Int[]; lfs = Float64[]; chosen = Dict{String,Float64}()
figL = Figure(size = (450 * length(fits), 380))
for (gi, f) in enumerate(fits)
    grid = gridof(f)
    ds = NCDataset(f)
    lam = Float64.(ds["lambda"][:]); chi2 = Float64.(ds["chi2"][:]); rough = Float64.(ds["roughness"][:])
    it = Int.(ds["iterations"][:]); ncl = Int.(ds["n_clip"][:]); U = Float64.(ds["u"][:, :]); fid = Int.(ds["free_id"][:])
    fit = Int.(ds["fit"][:, :]) .== 1; margin = Int.(ds["margin"][:, :]) .== 1; nfit = count(fit)
    ℓc, curv = corner(chi2, rough, lam)
    disc = findall(l -> chi2[l] / nfit <= 1, eachindex(lam))
    ℓd = isempty(disc) ? nothing : disc[argmax(lam[disc])]
    ℓ = haskey(given, grid) ? findfirst(==(given[grid]), lam) : something(ℓd, ℓc)
    ℓ === nothing && error("$grid: lambda $(given[grid]) is not in the fit")
    chosen[grid] = lam[ℓ]
    Hp = Float64.(ds["H_prior"][:, :]); Hf = Float64.(ds["H"][:, :, ℓ]); Hobs = Float64.(ds["H_obs"][:, :])
    dom = Int.(nomissing(ds["taub_domain"][:, :], 0)); x = Float64.(ds["xc"][:]); y = Float64.(ds["yc"][:])
    rp = ((Hp .- Hobs) ./ Float64.(ds["sigma"][:, :]))[fit]
    close(ds)

    println(io, "## ", grid, "\n\n", length(fid), " free domains, ", nfit, " fit cells; chi2 / n with Gowan's tau_b ",
            @sprintf("%.3f", sum(rp .^ 2) / nfit), ". L-curve (discrepancy chi2 / n <= 1: lambda = ", ℓd === nothing ? "none reaches it" : lam[ℓd],
            "; the corner: lambda = ", lam[ℓc], "; chosen: ", lam[ℓ], "):\n")
    println(io, "| lambda | chi2 / n | roughness [km] | curvature | iterations | at the bound | factor range | median abs(u) |\n|---|---|---|---|---|---|---|---|")
    for l in eachindex(lam)
        @printf(io, "| %g%s | %.3f | %.4g | %s | %d | %d | %.2f-%.2f | %.3f |\n", lam[l], l == ℓ ? " (chosen)" : "", chi2[l] / nfit, rough[l],
                haskey(curv, l) && isfinite(curv[l]) ? @sprintf("%.2f", curv[l]) : "-", it[l], ncl[l], exp(minimum(U[:, l])),
                exp(maximum(U[:, l])), median(abs.(U[:, l])))
    end
    println(io)
    ob = pd_obs(FES, grid)
    pd_table(io, ob, margin, ["Gowan's tau_b" => Hp, "Stage 1" => Hf]; title = "Present day per region (H_obs: $(PD_GRIDS[grid].topo))")
    append!(ids, fid); append!(lfs, U[:, ℓ])

    # L-curve
    ax = Axis(figL[1, gi]; title = grid, xlabel = "chi2 / n", ylabel = gi == 1 ? "roughness sum w (u_d - u_d')^2 [km]" : "",
              xscale = log10, yscale = log10)
    k = findall(l -> l > 0, lam)
    scatterlines!(ax, chi2[k] ./ nfit, max.(rough[k], 1e-3); color = :grey30)
    for l in k; text!(ax, chi2[l] / nfit, max(rough[l], 1e-3); text = @sprintf(" %g", lam[l]), fontsize = 9); end
    scatter!(ax, [chi2[ℓ] / nfit], [max(rough[ℓ], 1e-3)]; color = :firebrick, markersize = 14)

    # maps: H_obs, H - H_obs before and after, the factor
    uc = fill(NaN, size(dom)); pos = Dict(d => i for (i, d) in enumerate(fid))
    for i in eachindex(dom); (margin[i] && haskey(pos, dom[i])) && (uc[i] = U[pos[dom[i]], ℓ]); end
    nan(A) = (B = copy(A); B[.!margin] .= NaN; B)
    fig = Figure(size = (1500, 420))
    xs, ys = (findfirst(vec(any(margin; dims = 2))), findlast(vec(any(margin; dims = 2)))), (findfirst(vec(any(margin; dims = 1))), findlast(vec(any(margin; dims = 1))))
    lims = (x[max(xs[1] - 3, 1)], x[min(xs[2] + 3, length(x))], y[max(ys[1] - 3, 1)], y[min(ys[2] + 3, length(y))])
    for (c, (ttl, F, cm, cr)) in enumerate((("H_obs [m]", nan(Hobs), :viridis, (0, 4000)),
                                           ("H - H_obs, Gowan's tau_b [m]", nan(Hp .- Hobs), :RdBu, (-800, 800)),
                                           ("H - H_obs, Stage 1 [m]", nan(Hf .- Hobs), :RdBu, (-800, 800)),
                                           ("tau_b factor exp(u), lambda = $(lam[ℓ])", exp.(uc), :PuOr, (0.1, 10))))
        ax = Axis(fig[1, 2c-1]; title = ttl, aspect = DataAspect(), limits = lims)
        hidedecorations!(ax)
        hm = c == 4 ? heatmap!(ax, x, y, F; colormap = cm, colorrange = cr, colorscale = log10, nan_color = :grey90) :
                      heatmap!(ax, x, y, F; colormap = c == 1 ? cm : Reverse(cm), colorrange = cr, nan_color = :grey90)
        Colorbar(fig[1, 2c], hm; height = Relative(0.8))
    end
    save(joinpath(out, "maps_$(grid).png"), fig)
end
save(joinpath(out, "lcurve.png"), figL)

f = joinpath(out, "taub_domain_factor.nc"); isfile(f) && rm(f)
NCDataset(f, "c") do ds
    defDim(ds, "domain", length(ids))
    defVar(ds, "domain_id", Int32.(ids), ("domain",); attrib = ["long_name" => "domain id of the tau_b product"])
    defVar(ds, "log_factor", lfs, ("domain",); attrib = ["units" => "1", "long_name" => "natural log of the tau_b factor of the domain, constant in time"])
    ds.attrib["source"] = "tests/plastic_taub_fit.f90, tools/levante/plastic_taub_fit_report.jl: present-day fit (Stage 1)"
    ds.attrib["lambda"] = join(["$g $(chosen[g])" for g in sort(collect(keys(chosen)))], ", ")
end
close(io)
println("wrote ", out, ": ", length(ids), " domains, lambda ", chosen)
