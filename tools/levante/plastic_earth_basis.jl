# A basis of 3-D Earths from expert viscosity fields (docs/simulations.md: E3D-PCA): the mean and the principal components of
# log10 eta over the fields of FesmData's SolidEarth-Bagge2021 (its 8 structures), SolidEarth-Lloyd2024 and Pan et al. (2022),
# on Bagge's grid (GLOBAL-T170, its 164 depths to the core-mantle boundary). Each field is clamped to [lo, hi] first (the
# cratonic roots near 30 dex would dominate the spread; above about 23 dex the mantle is elastic over a glacial cycle, which the
# lid threshold of VILMA2 decides, visc_3d_lid_log10max). The fields count equally; the components are weighted by area (cos lat)
# and volume (r^2 dr). weights = "field" counts the fields equally (v1), "study" the three studies equally (Bagge's structures share
# one weight). lid (v2: 22.5) gives all fields one lithosphere: elastic (hi) from the surface down to where log10 eta first drops
# below lid. A component is stored as the pattern of one (weighted) standard deviation over the fields, so a field of the basis is
#   log10 eta = mean + sum_k c_k pc_k,   c_k ~ N(0, 1),
# and the experts sit at c_k = coef(field, k). Fields of the basis: earth-field (a VILMA2 lon-lat-r file, visc_3d_structure = 0).
#   julia --project=tools/levante tools/levante/plastic_ensemble.jl earth-basis <out.nc> <pan visc.nc> <pan litho.nc> [kmax=6] [lo=19.5] [hi=24] [weights=field|study] [lid=22.5]
#   julia --project=tools/levante tools/levante/plastic_ensemble.jl earth-field <basis.nc> <out.nc> c1 c2 ...
# Pan et al. (2022) is not released in FesmData yet: its viscosity is vilma2's input/pan2022.nc (log10 eta on lon, lat, radius), its
# lithosphere thickness isostasy_data's earth_structure/lithothickness/pan2022.nc (T_litho [km]). Pan's viscosity field starts below
# the lithosphere (Seakon carries the lithosphere as a thickness): above T_litho the field is set to hi (elastic), as in the others.
# The 1-D Earths of the broad Earth search (docs/dev/earth-search-design.md, wave A) from 5 radial parameters: earth1d_par,
# earth-search-sample (at the end of this file).
# Included by plastic_esmda.jl (after plastic_earth_shift.jl).
using DelimitedFiles, LinearAlgebra, NCDatasets, Printf, Random

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

function earth_basis(out, panfile, panlitho; kmax = 6, lo = 19.5, hi = 24.0, weights = "field", lid = nothing)
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
    issorted(dep) || error("$bf: depth not ascending")
    for (nm, F) in fields
        any(!isfinite, F) && error("$nm: non-finite values")
        F .= clamp.(F, Float32(lo), Float32(hi))
        # one lithosphere definition for all: elastic (hi) down to the first depth where log10 eta drops below lid
        lid === nothing || for j in eachindex(lat), i in eachindex(lon)
            k = findfirst(<(Float32(lid)), view(F, i, j, :)); k === nothing && (k = length(dep) + 1)
            F[i, j, 1:k-1] .= Float32(hi)
        end
    end
    M = length(fields); N = length(lon) * length(lat) * length(dep)
    # the weight of each field: equal ("field"), or equal per study with Bagge's 8 structures sharing one ("study")
    nb = count(f -> startswith(f[1], "B21"), fields)
    wm = weights == "study" ? [startswith(nm, "B21") ? 1 / (3nb) : 1 / 3 for (nm, _) in fields] :
         weights == "field" ? fill(1 / M, M) : error("weights = field | study, not $weights")
    μ = sum(wm[m] .* F for (m, (_, F)) in enumerate(fields))
    A = zeros(Float32, M, N); for (m, (_, F)) in enumerate(fields); A[m, :] .= vec(F .- μ); end
    # weights: area and volume, normalized
    rr = R_VILMA_KM .- dep ./ 1e3; dr = [abs(i == 1 ? (dep[2] - dep[1]) : i == length(dep) ? (dep[end] - dep[end-1]) : (dep[i+1] - dep[i-1]) / 2) for i in eachindex(dep)]
    w = vec([cosd(y) * r^2 * d for _ in lon, y in lat, (r, d) in zip(rr, dr)]); w ./= sum(w)
    sw = sqrt.(wm)
    G = Symmetric(sw .* Float64.((A .* w') * A') .* sw')            # the weighted covariance of the fields
    ev = eigen(G); o = sortperm(ev.values, rev = true); λ = max.(ev.values[o], 0); U = ev.vectors[:, o]
    K = min(kmax, M - 1)
    pcs = [reshape(A' * Float32.(sw .* U[:, k]), length(lon), length(lat), length(dep)) for k in 1:K]   # one weighted sd
    coef = U[:, 1:K] ./ sw                                      # the experts in the coefficients: field m = mean + sum coef pc
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
        ds.attrib["weights"] = weights * ": " * join([@sprintf("%.4f", x) for x in wm], " "); ds.attrib["lid"] = lid === nothing ? "none" : string(lid)
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

# --- 1-D Earths of the broad search (docs/dev/earth-search-design.md, wave A) ----------------------------------------------------
# A 1-D Earth from the 5 radial parameters: the elastic lithosphere thickness [km] and log10 eta of the upper mantle (lithosphere
# base to 410 km), the transition zone (410-670 km) and the lower mantle above and below 1500 km (to the core, 2891 km).
const EARTH1D_PARAMS = [(name = "litho_km", lo = 50.0, hi = 200.0), (name = "log10_eta_um", lo = 19.5, hi = 21.5),
                        (name = "log10_eta_tz", lo = 20.0, hi = 22.5), (name = "log10_eta_lm1", lo = 21.0, hi = 23.5),
                        (name = "log10_eta_lm2", lo = 21.0, hi = 23.5)]
const EARTH1D_BASE_KM = [410.0, 670.0, 1500.0, 2891.0]   # the bases of UM, TZ, LM1 and LM2 [km depth]; below: the fluid core
# The layers [km depth, surface first] of the 1-D Earth θ = (L, um, tz, lm1, lm2) on the layering of a template: VILMA2 takes
# density and rigidity per layer as PREM's shell average (earth = "PREM"), so the template's boundaries (Gowan's 12 layers: the
# PREM discontinuities and the lower mantle in parts) are kept, with the boundaries of the viscosity (L, 410, 670, 1500 km)
# added; a template boundary within merge_km of one of those is dropped (no sliver layers). Elastic (rheology 0, eta 1e40) above
# L, Maxwell (1) below, the fluid core (2, eta 0) below 2891 km. (top, bot) [km depth], eta [Pa s], rheology per layer.
function earth1d_layers(θ, tmpl_depths_km; merge_km = 5.0)
    L = θ[1]; vb = [L; EARTH1D_BASE_KM]
    keep = [d for d in tmpl_depths_km if 0 < d < R_VILMA_KM && all(abs(d - b) > merge_km for b in vb)]
    b = sort(unique([0.0; keep; vb; R_VILMA_KM]))
    top = b[1:end-1]; bot = b[2:end]; eta = similar(top); rh = zeros(Int, length(top))
    for (i, (t, d)) in enumerate(zip(top, bot))
        mid = (t + d) / 2
        eta[i], rh[i] = mid < L ? (1e40, 0) : mid < 410 ? (10.0^θ[2], 1) : mid < 670 ? (10.0^θ[3], 1) :
                        mid < 1500 ? (10.0^θ[4], 1) : mid < 2891 ? (10.0^θ[5], 1) : (0.0, 2)
    end
    top, bot, eta, rh
end
# A VILMA2 1-D par file of θ: a copy of the template (a 1-D kryos par file: its spin-up and the rest as they are) with the &vilma
# layer arrays of earth1d_layers (r in m, surface first) and l_visc_3d = .false. (set_nml.py: only keys the template has).
function earth1d_par(tmpl, out, θ; note = "")
    v = nml_group(tmpl, "vilma")
    tmplr = nml_reals(v["r_top"]) ./ 1e3
    top, bot, eta, rh = earth1d_layers(θ, R_VILMA_KM .- tmplr)
    length(top) <= 16 || error("$(length(top)) layers: VILMA2 takes at most 16 (MAX_LAYER)")
    fr(d) = @sprintf("%.3fe3", R_VILMA_KM - d); fe(x) = x == 0 ? "0.0" : @sprintf("%.4e", x)
    cp(tmpl, out; force = true); chmod(out, 0o644)
    kv = ["vilma.n_layer=$(length(top))", "vilma.r_top=" * join(fr.(top), ", "), "vilma.r_bot=" * join(fr.(bot), ", "),
          "vilma.eta=" * join(fe.(eta), ", "), "vilma.rheology=" * join(rh, ", "), "vilma.l_visc_3d=.false."]
    Base.run(`python3 $(joinpath(@__DIR__, "set_nml.py")) $out $kv`)
    hdr = "! A 1-D Earth of the broad search (plastic_ensemble.jl earth-search-*$(isempty(note) ? "" : ", " * note)): lithosphere " *
          @sprintf("%.2f km, log10 eta UM %.3f, TZ %.3f, LM1 %.3f, LM2 %.3f", θ...) * "; the layering of $(basename(tmpl)) with the\n" *
          "! viscosity boundaries (L, 410, 670, 1500 km) added (earth1d_layers).\n"
    write(out, hdr * read(out, String))
    out
end
# Wave A of the search: n 1-D Earths in a Latin hypercube, uniform over EARTH1D_PARAMS. Writes in outdir E<nnn>.nml (earth1d_par
# of the template, a 1-D par file such as par/kryos_NH-ANT_plastic_vilma2_b1.nml) and members.csv (the parameters, rounded as in
# the par files: L to 0.01 km, log10 eta to 0.001).
#   julia --project=tools/levante tools/levante/plastic_ensemble.jl earth-search-sample <template par> <outdir> [n=200] [seed=1]
#   julia --project=tools/levante tools/levante/plastic_ensemble.jl earth-search-par <template par> <out.nml> L um tz lm1 lm2
function earth_search_sample(tmpl, outdir; n = 200, seed = 1)
    rng = MersenneTwister(seed); np = length(EARTH1D_PARAMS)
    U = hcat([(randperm(rng, n) .- rand(rng, n)) ./ n for _ in 1:np]...)     # n x np in (0, 1)
    mkpath(outdir)
    open(joinpath(outdir, "members.csv"), "w") do io
        println(io, "member,", join([p.name for p in EARTH1D_PARAMS], ","))
        for m in 1:n
            θ = [p.lo + U[m, k] * (p.hi - p.lo) for (k, p) in enumerate(EARTH1D_PARAMS)]
            θ = [round(θ[1], digits = 2); round.(θ[2:end], digits = 3)]
            nm = @sprintf("E%03d", m)
            earth1d_par(tmpl, joinpath(outdir, nm * ".nml"), θ; note = "$nm of $n, seed $seed")
            println(io, nm, ",", @sprintf("%.2f", θ[1]), ",", join([@sprintf("%.3f", x) for x in θ[2:end]], ","))
        end
    end
    println("wrote ", n, " 1-D Earths in ", outdir, " (", join([@sprintf("%s [%g, %g]", p.name, p.lo, p.hi) for p in EARTH1D_PARAMS], ", "), ")")
end

# --- 3-D Earths of the broad search (wave B) -----------------------------------------------------------------------------------
# 13 parameters: the 5 radial ones of EARTH1D_PARAMS (the profile P), the amplitude α of the experts' lateral anomaly A (the basis
# mean minus its radial mean), c2 and c3 (the lateral part of the basis components 2 and 3, over the experts' range of their
# coefficients) and upper-mantle shifts [dex] in 5 regions (the masks of earth_regions_v2.nc; its Fennoscandia/Barents region is
# split, below). The field
#   log10 eta = P(z) + α A + c2 pc2' + c3 pc3' + Σ_j s_j w_j(x) [z < 670 km],   clamped to [19.5, 30], 30 above the lithosphere,
# with ' the lateral part (minus the area mean per depth): the radial profile is P's alone.
const EARTH3D_LATERAL = ["alpha", "c2", "c3"]
const EARTH3D_REGIONS = ["s_fenno", "s_barents", "s_laurentide", "s_greenland", "s_wais"]
const EARTH3D_SHIFT = 1.0             # [dex] the range of the regional shifts, ±
const EARTH3D_CLAMP = (19.5, 30.0)
# Fennoscandia/Barents of the masks split: the Barents share s(lon, lat) = max(clamp((lat - 69) / 4), clamp((lon - 40) / 10)) (north of
# 71-73°N or east of 40-50°E: Svalbard, Franz Josef Land, Novaya Zemlya, the Kara Sea; the mainland and the British Isles below), the
# rest Fennoscandia; the two add to the region's mask.
barents_share(lon, lat) = max(clamp((lat - 69.0) / 4.0, 0.0, 1.0), clamp((mod(lon + 180, 360) - 180 - 40.0) / 10.0, 0.0, 1.0))
# the parameter ranges of wave B: the radial ones, α in [0, 1.5], c2 and c3 over the experts' coefficients of the basis
function earth3d_params(basis)
    coef = NCDataset(b -> Float64.(b["coef"][:, :]), basis)              # field x component
    [EARTH1D_PARAMS; (name = "alpha", lo = 0.0, hi = 1.5); (name = "c2", lo = minimum(coef[:, 2]), hi = maximum(coef[:, 2]));
     (name = "c3", lo = minimum(coef[:, 3]), hi = maximum(coef[:, 3])); [(name = r, lo = -EARTH3D_SHIFT, hi = EARTH3D_SHIFT) for r in EARTH3D_REGIONS]]
end
# the parts of the field on the basis grid: A, pc2', pc3' (lon, lat, depth), the regional weights (lon, lat, 5) and the depths [km]
function earth3d_parts(basis, regions)
    NCDataset(basis) do b
        lon = Float64.(b["lon"][:]); lat = Float64.(b["lat"][:]); dep = Float64.(b["depth"][:]) ./ 1e3
        w = cosd.(lat) ./ sum(cosd.(lat))
        lateral(F) = F .- sum(F .* w', dims = (1, 2)) ./ size(F, 1)      # minus the area mean per depth
        A = lateral(Float64.(b["log10_visc_mean"][:, :, :]))
        p2 = lateral(Float64.(b["pc"][:, :, :, 2])); p3 = lateral(Float64.(b["pc"][:, :, :, 3]))
        mlon, mlat, M = NCDataset(r -> (Float64.(r["lon"][:]), Float64.(r["lat"][:]), Float64.(r["mask"][:, :, :])), regions)
        names = split(NCDataset(r -> r.attrib["region_names"], regions), ";")
        k(n) = (i = findfirst(==(n), names); i === nothing ? error("$regions: no region $n") : i)
        near(a, x) = argmin(abs.(a .- x))
        W = zeros(length(lon), length(lat), 5)
        for (j, y) in enumerate(lat), (i, x) in enumerate(lon)
            ii = near(mlon, mod(x + 180, 360) - 180); jj = near(mlat, y); s = barents_share(x, y); fb = M[ii, jj, k("Fennoscandia/Barents")]
            W[i, j, :] = [fb * (1 - s), fb * s, M[ii, jj, k("Laurentide")], M[ii, jj, k("Greenland")], M[ii, jj, k("West Antarctica")]]
        end
        (lon = lon, lat = lat, dep = dep, A = A, p2 = p2, p3 = p3, W = W)
    end
end
# the radial profile P at depth d [km]: 30 (elastic) above the lithosphere
earth1d_profile(θ, d) = d < θ[1] ? EARTH3D_CLAMP[2] : d < 410 ? θ[2] : d < 670 ? θ[3] : d < 1500 ? θ[4] : θ[5]
function earth3d_field(parts, out, θ)
    P = parts; F = similar(P.A)
    for (kd, d) in enumerate(P.dep)
        rad = earth1d_profile(θ, d)
        sh = d < 670 ? sum(θ[8+j] .* view(P.W, :, :, j) for j in 1:5) : 0.0
        F[:, :, kd] .= d < θ[1] ? EARTH3D_CLAMP[2] :
                       clamp.(rad .+ θ[6] .* view(P.A, :, :, kd) .+ θ[7] .* view(P.p2, :, :, kd) .+ θ[8] .* view(P.p3, :, :, kd) .+ sh, EARTH3D_CLAMP...)
    end
    o = sortperm(P.dep, rev = true)                                     # radius ascending
    isfile(out) && rm(out)
    NCDataset(out, "c") do ds
        defVar(ds, "lon", P.lon, ("lon",); attrib = ["units" => "degrees_east"])
        defVar(ds, "lat", P.lat, ("lat",); attrib = ["units" => "degrees_north"])
        defVar(ds, "r", R_VILMA_KM * 1e3 .- 1e3 .* P.dep[o], ("r",); attrib = ["units" => "m", "long_name" => "radius"])
        defVar(ds, "eta", Float32.(F[:, :, o]), ("lon", "lat", "r"); attrib = ["units" => "log10(Pa s)"], deflatelevel = 4)
        ds.attrib["parameters"] = join([@sprintf("%.4f", x) for x in θ], " ")
    end
    out
end
# the 3-D par file of θ: the template (a 3-D par file: visc_3d_file, visc_3d_structure, l_visc_3d and its spin-up) with the layers
# of P on the layering of the 1-D template (earth1d_layers: the same layers as the 1-D run of P whose r the 3-D run takes)
function earth3d_par(tmpl3, tmpl1, field, out, θ; note = "")
    v = nml_group(tmpl3, "vilma"); get(v, "l_visc_3d", ".false.") == ".true." || error("$tmpl3: not a 3-D par file")
    top, bot, eta, rh = earth1d_layers(θ[1:5], R_VILMA_KM .- nml_reals(nml_group(tmpl1, "vilma")["r_top"]) ./ 1e3)
    fr(d) = @sprintf("%.3fe3", R_VILMA_KM - d); fe(x) = x == 0 ? "0.0" : @sprintf("%.4e", x)
    cp(tmpl3, out; force = true); chmod(out, 0o644)
    kv = ["vilma.n_layer=$(length(top))", "vilma.r_top=" * join(fr.(top), ", "), "vilma.r_bot=" * join(fr.(bot), ", "),
          "vilma.eta=" * join(fe.(eta), ", "), "vilma.rheology=" * join(rh, ", "), "vilma.visc_3d_file=\"$(abspath(field))\"", "vilma.visc_3d_structure=0"]
    Base.run(`python3 $(joinpath(@__DIR__, "set_nml.py")) $out $kv`)
    write(out, "! A 3-D Earth of the broad search (plastic_ensemble.jl earth-search-sample3d$(isempty(note) ? "" : ", " * note)): " *
               join([@sprintf("%.3f", x) for x in θ], " ") * "\n! (earth3d_field: P + α A + c2 pc2' + c3 pc3' + regional shifts; the layers of P)\n" * read(out, String))
    out
end
# Wave B: n 3-D Earths. The radial parameters from the points of the wave-A emulator's dense sample (plastic_earth_emulate.jl
# sample.csv) not ruled out under every ice history (Imax_<ice> <= cut for all), spread by greedy maximin on [0, 1] of their ranges
# (the first: the point nearest the centre of the kept ones); the other 8 in a Latin hypercube over their ranges. Writes in outdir
# per member E<nnn>.nc (the field), E<nnn>.nml (the 3-D par file), E<nnn>_1d.nml (the 1-D par file of its P) and members.csv.
#   julia --project=tools/levante tools/levante/plastic_ensemble.jl earth-search-sample3d <sample.csv> <basis.nc> <regions.nc> <3-D par> <1-D par> <outdir> [n=50] [seed=1] [cut=3]
function earth_search_sample3d(samplef, basis, regions, tmpl3, tmpl1, outdir; n = 50, seed = 1, cut = 3.0)
    S, h = readdlm(samplef, ',', Float64; header = true); h = String.(vec(h))
    ir = [findfirst(==(p.name), h) for p in EARTH1D_PARAMS]; any(isnothing, ir) && error("$samplef: not the radial parameters")
    keep = [all(S[i, k] <= cut for k in findall(startswith("Imax_"), h)) for i in axes(S, 1)]
    count(keep) >= n || error("$samplef: $(count(keep)) points not ruled out, fewer than $n")
    R = S[keep, ir]; lo = [p.lo for p in EARTH1D_PARAMS]; hi = [p.hi for p in EARTH1D_PARAMS]; U = (R .- lo') ./ (hi .- lo)'
    c = vec(sum(U, dims = 1)) ./ size(U, 1); pick = [argmin(vec(sum(abs2, U .- c', dims = 2)))]
    dmin = vec(sum(abs2, U .- U[pick[1], :]', dims = 2))
    while length(pick) < n
        k = argmax(dmin); push!(pick, k); dmin = min.(dmin, vec(sum(abs2, U .- U[k, :]', dims = 2)))
    end
    P = earth3d_params(basis); rng = MersenneTwister(seed); nl = length(P) - 5
    L = hcat([(randperm(rng, n) .- rand(rng, n)) ./ n for _ in 1:nl]...)
    parts = earth3d_parts(basis, regions); mkpath(outdir)
    open(joinpath(outdir, "members.csv"), "w") do io
        println(io, "member,", join([p.name for p in P], ","))
        for m in 1:n
            θ = [R[pick[m], :]; [P[5+k].lo + L[m, k] * (P[5+k].hi - P[5+k].lo) for k in 1:nl]]
            θ = [round(θ[1], digits = 2); round.(θ[2:end], digits = 3)]; nm = @sprintf("E%03d", m)
            f = earth3d_field(parts, joinpath(outdir, nm * ".nc"), θ)
            earth3d_par(tmpl3, tmpl1, f, joinpath(outdir, nm * ".nml"), θ; note = "$nm of $n, seed $seed")
            earth1d_par(tmpl1, joinpath(outdir, nm * "_1d.nml"), θ[1:5]; note = "the radial profile of wave-B $nm")
            println(io, nm, ",", @sprintf("%.2f", θ[1]), ",", join([@sprintf("%.3f", x) for x in θ[2:end]], ","))
        end
    end
    println("wrote ", n, " 3-D Earths in ", outdir, " (radial: ", count(keep), " of ", size(S, 1), " points not ruled out; ",
            join([@sprintf("%s [%.2f, %.2f]", p.name, p.lo, p.hi) for p in P[6:end]], ", "), ")")
end
