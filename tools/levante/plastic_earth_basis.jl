# A basis of 3-D Earths from expert viscosity fields (docs/simulations.md: E3D-PCA): the mean and the principal components of
# log10 eta over the fields of FesmData's SolidEarth-Bagge2021 (its 8 structures), SolidEarth-Lloyd2024 and Pan et al. (2022),
# on Bagge's grid (GLOBAL-T170, its 164 depths to the core-mantle boundary). Each field is clamped to [lo, hi] first (the
# cratonic roots near 30 dex would dominate the spread; above about 23 dex the mantle is elastic over a glacial cycle, which the
# lid threshold of VILMA2 decides, visc_3d_lid_log10max). The fields count equally; the components are weighted by area (cos lat)
# and volume (r^2 dr). A component is stored as the pattern of one standard deviation over the fields, so a field of the basis is
#   log10 eta = mean + sum_k c_k pc_k,   c_k ~ N(0, 1),
# and the experts sit at c_k = coef(field, k). Fields of the basis: earth-field (a VILMA2 lon-lat-r file, visc_3d_structure = 0).
#   julia --project=tools/levante tools/levante/plastic_ensemble.jl earth-basis <out.nc> <pan visc.nc> <pan litho.nc> [kmax=6] [lo=19.5] [hi=24]
#   julia --project=tools/levante tools/levante/plastic_ensemble.jl earth-field <basis.nc> <out.nc> c1 c2 ...
# Pan et al. (2022) is not released in FesmData yet: its viscosity is vilma2's input/pan2022.nc (log10 eta on lon, lat, radius), its
# lithosphere thickness isostasy_data's earth_structure/lithothickness/pan2022.nc (T_litho [km]). Pan's viscosity field starts below
# the lithosphere (Seakon carries the lithosphere as a thickness): above T_litho the field is set to hi (elastic), as in the others.
# Included by plastic_esmda.jl (after plastic_earth_shift.jl).
using LinearAlgebra, NCDatasets, Printf, Random

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

function earth_basis(out, panfile, panlitho; kmax = 6, lo = 19.5, hi = 24.0)
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
    P = NCDataset(panfile) do ds
        regrid(Float64.(ds["lon"][:]), Float64.(ds["lat"][:]), R_VILMA_KM * 1e3 .- Float64.(ds["r"][:]),
               Float32.(nomissing(ds["eta"][:, :, :], NaN)), lon, lat, dep)
    end
    T = NCDataset(panlitho) do ds                                           # [km] on Bagge's grid (one level)
        tl = Float64.(ds["lon"][:]); keep = findall(<(360.0), tl .- tl[1])      # the last column repeats the first
        regrid(tl[keep], Float64.(ds["lat"][:]), [0.0, 1.0], repeat(Float32.(ds["T_litho"][keep, :]), 1, 1, 2), lon, lat, [0.0])[:, :, 1]
    end
    for k in eachindex(dep), j in eachindex(lat), i in eachindex(lon)
        dep[k] < 1e3 * T[i, j] && (P[i, j, k] = Float32(hi))
    end
    push!(fields, ("P22", P))
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

# An Earth-only ensemble of the basis (docs/simulations.md: E3D-PCA-E<nnn>): n members, the first K coefficients in a Latin hypercube
# uniform over the range of the experts' coefficients extended by `ext` of it on each side (so the members span the experts, not a
# prior centred on the more numerous Bagge structures), and the lid threshold (VILMA2 visc_3d_lid_log10max: Maxwell elements wholly
# above lid_depth km stiffer than it are elastic) uniform in `lid`. Writes in outdir: E<nnn>.nc (earth-field), E<nnn>.nml (the 3-D
# par file with visc_3d_file = E<nnn>.nc, visc_3d_structure = 0 and the lid in its &vilma group, for gia_refcorr.sh) and members.csv.
#   julia --project=tools/levante tools/levante/plastic_ensemble.jl earth-sample <basis.nc> <3-D par.nml> <outdir> [n=30] [K=4] [seed=1]
function earth_sample(basis, par, outdir; n = 30, K = 4, seed = 1, ext = 0.1, lid = (21.5, 23.5), lid_depth = 400.0)
    coef = NCDataset(b -> Float64.(b["coef"][:, 1:K]), basis)
    lo = [minimum(coef[:, k]) for k in 1:K]; hi = [maximum(coef[:, k]) for k in 1:K]; w = hi .- lo
    lo .-= ext .* w; hi .+= ext .* w
    rng = MersenneTwister(seed)
    U = hcat([(randperm(rng, n) .- rand(rng, n)) ./ n for _ in 1:K+1]...)     # n x (K+1) in (0, 1)
    mkpath(outdir); txt = read(par, String)
    occursin(r"^\s*l_visc_3d\s*=\s*\.true\."mi, txt) || error("$par: not a 3-D par file")
    open(joinpath(outdir, "members.csv"), "w") do io
        println(io, "member,", join(["c$k" for k in 1:K], ","), ",lid_log10max")
        for m in 1:n
            c = lo .+ U[m, 1:K] .* (hi .- lo); λ = lid[1] + U[m, K+1] * (lid[2] - lid[1]); nm = @sprintf("E%03d", m)
            f = earth_field(basis, joinpath(outdir, nm * ".nc"), c)
            t = replace(txt, r"^(\s*visc_3d_file\s*=\s*).*$"m => SubstitutionString("\\1\"$(abspath(f))\""),
                             r"^(\s*visc_3d_structure\s*=\s*).*$"m => s"\g<1>0")
            t = replace(t, r"(^&vilma\b.*?)(^\s*/\s*$)"ms => SubstitutionString(
                "\\1    visc_3d_lid_depth    = $(lid_depth * 1e3)   ! earth-sample: the lid rule above this depth [m]\n" *
                "    visc_3d_lid_log10max = $(round(λ, digits = 3))   ! earth-sample: stiffer is elastic [dex]\n\\2"))
            write(joinpath(outdir, nm * ".nml"), "! earth-sample $(nm): $(abspath(basis)), c = $(join(round.(c, digits = 3), " ")), lid $(round(λ, digits = 3))\n" * t)
            println(io, nm, ",", join([@sprintf("%.4f", x) for x in c], ","), ",", @sprintf("%.3f", λ))
        end
    end
    println("wrote ", n, " members in ", outdir, " (coefficients in ", join([@sprintf("[%.2f, %.2f]", a, b) for (a, b) in zip(lo, hi)], " "), ")")
end
