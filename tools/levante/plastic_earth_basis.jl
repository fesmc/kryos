# A basis of 3-D Earths from expert viscosity fields (docs/simulations.md: E3D-PCA): the mean and the principal components of
# log10 eta over the fields of FesmData's SolidEarth-Bagge2021 (its 8 structures), SolidEarth-Lloyd2024 and Pan et al. (2022),
# on Bagge's grid (GLOBAL-T170, its 164 depths to the core-mantle boundary). Each field is clamped to [lo, hi] first (the
# cratonic roots near 30 dex would dominate the spread; above about 23 dex the mantle is elastic over a glacial cycle, which the
# lid threshold of VILMA2 decides, visc_3d_lid_log10max). The fields count equally; the components are weighted by area (cos lat)
# and volume (r^2 dr). A component is stored as the pattern of one standard deviation over the fields, so a field of the basis is
#   log10 eta = mean + sum_k c_k pc_k,   c_k ~ N(0, 1),
# and the experts sit at c_k = coef(field, k). Fields of the basis: earth-field (a VILMA2 lon-lat-r file, visc_3d_structure = 0).
#   julia --project=tools/levante tools/levante/plastic_ensemble.jl earth-basis <out.nc> <pan2022.nc> [kmax=6] [lo=19.5] [hi=24]
#   julia --project=tools/levante tools/levante/plastic_ensemble.jl earth-field <basis.nc> <out.nc> c1 c2 ...
# Pan et al. (2022) is not released in FesmData yet: its file is vilma2's input/pan2022.nc (log10 eta on lon, lat, radius).
# Included by plastic_esmda.jl (after plastic_earth_shift.jl).
using LinearAlgebra, NCDatasets, Printf

# linear weights of x in an ascending axis a (clamped at the ends)
function linw(a, x)
    x <= a[1] && return (1, 1, 0.0); x >= a[end] && return (length(a), length(a), 0.0)
    k = searchsortedlast(a, x); (k, k + 1, (x - a[k]) / (a[k+1] - a[k]))
end
# a field F(lon, lat, depth) on ascending axes (lon periodic in [0, 360)) to the target axes, bilinear horizontally, linear in depth
function regrid(slon, slat, sdep, F, tlon, tlat, tdep)
    p = sortperm(mod.(slon, 360)); sl = mod.(slon[p], 360); q = sortperm(slat); r = sortperm(sdep)
    F = F[p, q, r]; slat = slat[q]; sdep = sdep[r]; n = length(sl)
    wx = map(tlon) do x
        xm = mod(x, 360); i = searchsortedlast(sl, xm)
        i0 = i == 0 ? n : i; i1 = i == n || i == 0 ? 1 : i + 1
        x0 = i == 0 ? sl[n] - 360 : sl[i0]; x1 = i == n ? sl[1] + 360 : sl[i1]
        (i0, i1, (xm - x0) / (x1 - x0))
    end
    wy = [linw(slat, y) for y in tlat]; wd = [linw(sdep, d) for d in tdep]
    G = zeros(Float32, length(tlon), length(tlat), length(tdep))
    Threads.@threads for kd in eachindex(tdep)
        d0, d1, fd = wd[kd]
        L = (1 - fd) .* view(F, :, :, d0) .+ fd .* view(F, :, :, d1)
        for (j, (j0, j1, fy)) in enumerate(wy), (i, (i0, i1, fx)) in enumerate(wx)
            G[i, j, kd] = (1 - fy) * ((1 - fx) * L[i0, j0] + fx * L[i1, j0]) + fy * ((1 - fx) * L[i0, j1] + fx * L[i1, j1])
        end
    end
    G
end

function earth_basis(out, panfile; kmax = 6, lo = 19.5, hi = 24.0)
    fes = ENV["FESMDATA"]
    bf = joinpath(fes, "Global", "GLOBAL-T170", "GLOBAL-T170_SolidEarth-Bagge2021.nc")
    lon, lat, dep, B, snames = NCDataset(bf) do ds
        Float64.(ds["lon"][:]), Float64.(ds["lat"][:]), Float64.(ds["depth"][:]), Float32.(ds["log10_visc"][:, :, :, :]),
        split(ds["structure"].attrib["structure_names"])
    end
    fields = [("B21-$s ($(snames[s]))", B[:, :, :, s]) for s in axes(B, 4)]
    lf = joinpath(fes, "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_SolidEarth-Lloyd2024.nc")
    push!(fields, ("L24", NCDataset(ds -> regrid(Float64.(ds["lon"][:]), Float64.(ds["lat"][:]), Float64.(ds["depth"][:]),
                                                  Float32.(nomissing(ds["log10_visc"][:, :, :], NaN)), lon, lat, dep), lf)))
    push!(fields, ("P22", NCDataset(panfile) do ds
        regrid(Float64.(ds["lon"][:]), Float64.(ds["lat"][:]), R_VILMA_KM * 1e3 .- Float64.(ds["r"][:]),
               Float32.(nomissing(ds["eta"][:, :, :], NaN)), lon, lat, dep)
    end))
    for (nm, F) in fields
        any(!isfinite, F) && error("$nm: non-finite values")
        F .= clamp.(F, Float32(lo), Float32(hi))
    end
    M = length(fields); N = length(lon) * length(lat) * length(dep)
    μ = sum(F for (_, F) in fields) ./ M
    A = zeros(Float32, M, N); for (m, (_, F)) in enumerate(fields); A[m, :] .= vec(F .- μ); end
    # weights: area and volume, normalized
    rr = R_VILMA_KM .- dep ./ 1e3; dr = [abs(i == 1 ? (dep[2] - dep[1]) : i == length(dep) ? (dep[end] - dep[end-1]) : (dep[i+1] - dep[i-1]) / 2) for i in eachindex(dep)]
    w = vec([cosd(y) * r^2 * d for _ in lon, y in lat, (r, d) in zip(rr, dr)]); w ./= sum(w)
    G = Symmetric(Float64.((A .* w') * A'))
    ev = eigen(G); o = sortperm(ev.values, rev = true); λ = max.(ev.values[o], 0); U = ev.vectors[:, o]
    K = min(kmax, M - 1)
    pcs = [reshape(A' * U[:, k], length(lon), length(lat), length(dep)) ./ Float32(sqrt(M - 1)) for k in 1:K]
    coef = sqrt(M - 1) .* U[:, 1:K]                             # the experts in the coefficients
    frac = λ ./ sum(λ)
    wa = cosd.(lat) ./ (length(lon) * sum(cosd.(lat)))
    prof = [sum(F[i, j, d] * wa[j] for i in eachindex(lon), j in eachindex(lat)) for (_, F) in fields, d in eachindex(dep)]   # area means
    isfile(out) && rm(out)
    NCDataset(out, "c") do ds
        defVar(ds, "lon", lon, ("lon",); attrib = ["units" => "degrees_east"])
        defVar(ds, "lat", lat, ("lat",); attrib = ["units" => "degrees_north"])
        defVar(ds, "depth", dep, ("depth",); attrib = ["units" => "m", "positive" => "down"])
        defVar(ds, "log10_visc_mean", μ, ("lon", "lat", "depth"); attrib = ["units" => "log10(Pa s)", "long_name" => "mean of the clamped fields"])
        P = zeros(Float32, length(lon), length(lat), length(dep), K); for k in 1:K; P[:, :, :, k] .= pcs[k]; end
        defVar(ds, "pc", P, ("lon", "lat", "depth", "component"); attrib = ["units" => "dex",
               "long_name" => "principal component: the pattern of one standard deviation over the fields"])
        defVar(ds, "variance_fraction", frac[1:K], ("component",))
        defVar(ds, "coef", Float64.(coef), ("field", "component"); attrib = ["long_name" => "the coefficients of the fields"])
        defVar(ds, "profile", Float64.(prof), ("field", "depth"); attrib = ["units" => "log10(Pa s)", "long_name" => "area mean of each clamped field"])
        ds.attrib["fields"] = join(first.(fields), "; "); ds.attrib["clamp"] = "$lo $hi"
        ds.attrib["title"] = "Basis of 3-D Earths: mean and principal components of log10 viscosity of expert fields (plastic_ensemble.jl earth-basis)"
    end
    println("wrote ", out, ": ", M, " fields, variance fractions ", join([@sprintf("%.3f", f) for f in frac[1:K]], " "),
            " (", @sprintf("%.3f", sum(frac[1:K])), " in ", K, ")")
    for (m, (nm, _)) in enumerate(fields); println("  ", rpad(nm, 26), join([@sprintf("%+6.2f", coef[m, k]) for k in 1:K], " ")); end
end

# a field of the basis, log10 eta = mean + sum_k c_k pc_k (clamped as the basis), as a VILMA2 lon-lat-r file (visc_3d_structure = 0,
# name_visc = "eta", name_visc_r = "r"); c may be shorter than the basis (the rest 0)
function earth_field(basis, out, c)
    NCDataset(basis) do b
        lo, hi = parse.(Float32, split(b.attrib["clamp"])); K = size(b["pc"], 4); length(c) <= K || error("$(length(c)) coefficients, $K components")
        F = Float32.(b["log10_visc_mean"][:, :, :])
        for (k, ck) in enumerate(c); F .+= Float32(ck) .* Float32.(b["pc"][:, :, :, k]); end
        F .= clamp.(F, lo, hi)
        dep = Float64.(b["depth"][:]); o = sortperm(dep, rev = true)            # radius ascending
        isfile(out) && rm(out)
        NCDataset(out, "c") do ds
            defVar(ds, "lon", Float64.(b["lon"][:]), ("lon",); attrib = ["units" => "degrees_east"])
            defVar(ds, "lat", Float64.(b["lat"][:]), ("lat",); attrib = ["units" => "degrees_north"])
            defVar(ds, "r", R_VILMA_KM * 1e3 .- dep[o], ("r",); attrib = ["units" => "m", "long_name" => "radius"])
            defVar(ds, "eta", F[:, :, o], ("lon", "lat", "r"); attrib = ["units" => "log10(Pa s)"], deflatelevel = 4)
            ds.attrib["coefficients"] = join([@sprintf("%.4f", x) for x in c], " "); ds.attrib["basis"] = abspath(basis)
        end
    end
    out
end

# the basis in a figure: the area-mean profile of each field and of the mean, and the components at three depths
#   julia --project=tools/levante tools/levante/plastic_ensemble.jl earth-basis-plot <basis.nc> <out.png> [depths_km=100,250,600]
function earth_basis_plot(basis, out; depths_km = [100.0, 250.0, 600.0])
    @eval using CairoMakie
    NCDataset(basis) do b
        lon = Float64.(b["lon"][:]); lat = Float64.(b["lat"][:]); dep = Float64.(b["depth"][:]) ./ 1e3
        names = split(b.attrib["fields"], "; "); P = Float64.(b["profile"][:, :]); K = size(b["pc"], 4); frac = b["variance_fraction"][:]
        lo = [x > 180 ? x - 360 : x for x in lon]; o = sortperm(lo)
        fig = Figure(size = (450 + 420 * length(depths_km), 330 * (K + 1)))
        ax = Axis(fig[1:K+1, 1], xlabel = "log10 eta (area mean)", ylabel = "depth [km]", yreversed = true, title = "the fields")
        for (m, nm) in enumerate(names)
            lines!(ax, P[m, :], dep, label = String(nm), linewidth = startswith(nm, "B21") ? 1.2 : 3,
                   color = startswith(nm, "B21") ? Makie.wong_colors()[1 + (m - 1) % 7] : (nm == "L24" ? :black : :red),
                   linestyle = startswith(nm, "B21") ? :solid : :dash)
        end
        μp = vec(sum(P, dims = 1) ./ size(P, 1)); lines!(ax, μp, dep, color = :gray30, linewidth = 4, label = "mean")
        axislegend(ax, position = :lb, labelsize = 11)
        for (c, dk) in enumerate(depths_km)
            kd = argmin(abs.(dep .- dk))
            for row in 0:K
                F = row == 0 ? Float64.(b["log10_visc_mean"][:, :, kd]) : Float64.(b["pc"][:, :, kd, row])
                lim = row == 0 ? (19.5, 24.0) : (-1.0, 1.0) .* max(maximum(abs, F), 1e-3)
                a = Axis(fig[row + 1, c + 1], title = row == 0 ? @sprintf("mean, %.0f km", dep[kd]) : @sprintf("pc %d (%.0f %%), %.0f km", row, 100frac[row], dep[kd]),
                         aspect = DataAspect(), xticksvisible = false, xticklabelsvisible = false, yticksvisible = false, yticklabelsvisible = false)
                heatmap!(a, lo[o], lat, F[o, :], colormap = row == 0 ? :viridis : cgrad(:RdBu, rev = true), colorrange = lim)
            end
        end
        save(out, fig, px_per_unit = 1)
    end
    println("wrote ", out)
end
