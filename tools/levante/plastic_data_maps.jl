# Data-model maps of imposed-ice runs on VILMA2 (docs/science/maps.md): the GNSS uplift rates (the data, each run's present-day
# rate, the residual), the RSL residuals at the GAPSLIP 2.1 index points in time windows, and RSL at the best-sampled sites.
#   julia --project=tools/levante tools/levante/plastic_data_maps.jl <out folder> <label>=<run folder> ...
# The runs: gia_refcorr.sh cases (the converged iteration, r.nc -> it<n>: plastic_ensemble.jl ref_dir) or a folder with
# vilma_out.nc (a one-pass case: its it0); each a row of the figures, the label its name (`ice-GLAC1D | E1D-PMIST`).
# The data, the groups and the error model are those of the scores: the scoring folder KRYOS_SCORE (default the Earth search's,
# /albedo/work/user/alrobi001/kryos-search/score/ens), its config.toml ([[gnss]] groups) and errmodel.toml. The model values
# as scored: the rate at the stations and on the maps, dz_bed_rate of the run's last record (pm_compare.jl read_rate); RSL at the
# index points of age <= 25 ka at the 2.5 kyr ages of the scores (pm_compare.jl series_run, score: d = model - obs and its sigma_i).
# The maps show obs - model.
# Gaussian processes (GP) on the sphere: Matérn-3/2 on the chord distance (lon-lat as 3-D vectors of radius 6371 km), the noise
# of each point its own sigma_i. Data alone (GNSS): a constant mean (GLS) and sigma_f, l by maximum marginal likelihood on log
# grids per domain. Residuals: zero mean, sigma_m and l of the error model (errmodel.toml; RSL: times exp(-|t - t_c| / L_t), t_c
# the window's mid age); the GNSS groups not fitted there (Greenland, Antarctica: GNSS_FIT) get sigma_m, l by maximum marginal
# likelihood over the runs together. The fields are drawn where the GP's sd is at most MASK of its prior sd.
# Maps: stereographic about the data (Antarctica: about the South Pole), the far field in longitude-latitude; the coast is
# C_ocean = 0.5 of the first run's last record. Writes to <out>:
#   gnss-<domain>.png   per GNSS group: the data over their GP mean (top), per run its rate field and the residuals
#   rsl-<group>.png     per RSL group: the residuals at the index points over their GP per window (columns) and run (rows)
#   sites-<group>.png   per RSL group: RSL at the NSITE sites with the most index points, the data and the runs
#   maps.md             the GP hyperparameters, and n, mean and RMS of the residuals per panel
using CairoMakie, KernelFunctions, LinearAlgebra, NCDatasets, Printf, Statistics, TOML

module EN
include("plastic_ensemble.jl")
end
const PM = EN.PM

const FES = get(ENV, "FESMDATA", "/work/ba1442/fesmdata")   # tools/levante/env.sh
const SCORE = get(ENV, "KRYOS_SCORE", "/albedo/work/user/alrobi001/kryos-search/score/ens")
const RE = 6371.0                    # [km]
const TMAX = 25000.0                 # [yr] the scored index points
const MASK = 0.8                     # a GP field is drawn where its sd <= MASK * the prior sd
const WINDOWS = [("0-6 ka", 0.0, 6000.0), ("8 ± 1 ka", 7000.0, 9000.0), ("12 ± 1 ka", 11000.0, 13000.0),
                 ("16 ± 1 ka", 15000.0, 17000.0)]
const MIN_WINDOW = 5                 # index points of a group for its window column
const NSITE = 3
const GNSS_DOMAINS = ["GNSS Hudson Bay" => "hudson-bay", "GNSS Fennoscandia" => "fennoscandia",
                      "GNSS Greenland" => "greenland", "GNSS Antarctica" => "antarctica"]
const GNSS_FIT = ["GNSS Greenland", "GNSS Antarctica"]      # their model error is not fitted in errmodel.toml (data-scores.md)
const LGRID = exp.(range(log(25.0), log(5000.0), length = 47))      # [km]
const SGRID = exp.(range(log(0.05), log(50.0), length = 61))        # [mm/yr]
const CMAP = Reverse(:RdBu)
const LAND, OCEAN = RGBf(0.86, 0.85, 0.82), RGBf(0.97, 0.97, 0.98)
const COLS = Makie.wong_colors()
const IPCOL, IPSTROKE = (:gray, 0.3), (:gray30, 0.6)
const LIMITING = [(-1, :utriangle, :steelblue), (1, :dtriangle, :sienna)]   # GAPSLIP type: marine (RSL above), terrestrial (below)

# --- Gaussian processes ---------------------------------------------------------------------------------------------------------
xyz(lon, lat) = reduce(hcat, [RE .* [cosd(b) * cosd(a), cosd(b) * sind(a), sind(b)] for (a, b) in zip(lon, lat)])
k32(ℓ) = Matern32Kernel() ∘ ScaleTransform(1 / ℓ)
tfac(t1, t2, Lt) = exp.(-abs.(t1 .- t2') ./ Lt)
# a GP conditioned on the points X (3 x n) at times t, values y with sd σ: prior sd s, length ℓ [km], L_t [yr], mean μ
struct GPc
    X::Matrix{Float64}; t::Vector{Float64}; ℓ::Float64; s::Float64; Lt::Float64; μ::Float64
    L::LowerTriangular{Float64,Matrix{Float64}}; α::Vector{Float64}
end
function gp_cond(X, t, y, σ, ℓ, s, Lt; μ = 0.0)
    F = cholesky(Symmetric(s^2 .* kernelmatrix(k32(ℓ), ColVecs(X)) .* tfac(t, t, Lt) + Diagonal(σ .^ 2)))
    GPc(X, t, ℓ, s, Lt, μ, F.L, F \ (y .- μ))
end
# mean and sd at the points Xg (3 x N) at time tg
function gp_pred(g::GPc, Xg, tg)
    N = size(Xg, 2); m = zeros(N); sd = zeros(N)
    for c in Iterators.partition(1:N, 4000)
        Ks = g.s^2 .* kernelmatrix(k32(g.ℓ), ColVecs(g.X), ColVecs(Xg[:, c])) .* exp.(-abs.(g.t .- tg) ./ g.Lt)
        m[c] = g.μ .+ Ks' * g.α
        sd[c] = sqrt.(max.(g.s^2 .- vec(sum(abs2, g.L \ Ks; dims = 1)), 0.0))
    end
    m, sd
end
# the log marginal likelihood of y under prior sd s on K0 (the unit kernel matrix), and the GLS mean (fitmean) or 0
function lml(K0, y, σ, s; fitmean)
    F = cholesky(Symmetric(s^2 .* K0 + Diagonal(σ .^ 2)))
    μ = fitmean ? (a = F \ ones(length(y)); dot(a, y) / sum(a)) : 0.0
    r = y .- μ
    -0.5 * (dot(r, F \ r) + logdet(F) + length(y) * log(2π)), μ
end
# s and ℓ by maximum marginal likelihood on LGRID x SGRID over the sets P (each (X, y, σ)), summed; on the edge of a grid: a warning
function gp_fit(P; fitmean)
    best = (ll = -Inf, ℓ = NaN, s = NaN, μ = NaN)
    for ℓ in LGRID
        K0 = [kernelmatrix(k32(ℓ), ColVecs(X)) for (X, _, _) in P]
        for s in SGRID
            v = [lml(K0[i], P[i][2], P[i][3], s; fitmean) for i in eachindex(P)]
            ll = sum(first, v); ll > best.ll && (best = (ll = ll, ℓ = ℓ, s = s, μ = v[1][2]))
        end
    end
    (best.ℓ in extrema(LGRID) || best.s in extrema(SGRID)) && println("gp_fit: on the edge of the grid: ", best)
    best
end
masked(m, sd, s) = [sd[i] <= MASK * s ? m[i] : NaN for i in eachindex(m)]

# --- maps -------------------------------------------------------------------------------------------------------------------------
# A map: stereographic about (lon0, lat0) [km], or longitude-latitude (plate) [degrees]; its limits and a regular grid over them
struct Domain
    plate::Bool; lon0::Float64; lat0::Float64; xlim::Tuple{Float64,Float64}; ylim::Tuple{Float64,Float64}
    xs::Vector{Float64}; ys::Vector{Float64}; glon::Matrix{Float64}; glat::Matrix{Float64}; X::Matrix{Float64}
end
function stereo(lon0, lat0, lon, lat)
    dl = lon - lon0; k = 2RE / (1 + sind(lat0) * sind(lat) + cosd(lat0) * cosd(lat) * cosd(dl))
    k * cosd(lat) * sind(dl), k * (cosd(lat0) * sind(lat) - sind(lat0) * cosd(lat) * cosd(dl))
end
function istereo(lon0, lat0, x, y)
    ρ = hypot(x, y); ρ == 0 && return lon0, lat0
    c = 2atan(ρ / 2RE)
    mod(lon0 + atand(x * sin(c), ρ * cosd(lat0) * cos(c) - y * sind(lat0) * sin(c)) + 180, 360) - 180,
    asind(cos(c) * sind(lat0) + y * sin(c) * cosd(lat0) / ρ)
end
proj(D, lon, lat) = D.plate ? (mod(lon, 360), lat) : stereo(D.lon0, D.lat0, lon, lat)
sphere_mean(lo, la) = (v = sum([cosd(b) * cosd(a), cosd(b) * sind(a), sind(b)] for (a, b) in zip(lo, la));
                       (atand(v[2], v[1]), atand(v[3], hypot(v[1], v[2]))))
function domain(lon, lat; pole = false, margin = 400.0, n = 160)
    lon0, lat0 = pole ? (0.0, -90.0) : sphere_mean(lon, lat)
    P = [stereo(lon0, lat0, a, b) for (a, b) in zip(lon, lat)]
    x0, x1 = extrema(first.(P)) .+ (-margin, margin); y0, y1 = extrema(last.(P)) .+ (-margin, margin)
    pole && (r = max(abs(x0), abs(x1), abs(y0), abs(y1)); (x0, x1, y0, y1) = (-r, r, -r, r))
    d = max(x1 - x0, y1 - y0) / n; xs = collect(x0:d:x1); ys = collect(y0:d:y1)
    G = [istereo(lon0, lat0, x, y) for x in xs, y in ys]
    Domain(false, lon0, lat0, (x0, x1), (y0, y1), xs, ys, first.(G), last.(G), xyz(vec(first.(G)), vec(last.(G))))
end
# longitude-latitude over the points (longitude in [0, 360)) and a margin [degrees]
function plate_domain(lon, lat; margin = 8.0, d = 0.5)
    x0, x1 = extrema(mod.(lon, 360)) .+ (-margin, margin); y0, y1 = extrema(lat) .+ (-margin, margin)
    xs = collect(x0:d:x1); ys = collect(y0:d:y1)
    G = [(x, y) for x in xs, y in ys]
    Domain(true, 0.0, 0.0, (x0, x1), (y0, y1), xs, ys, first.(G), last.(G), xyz(vec(first.(G)), vec(last.(G))))
end
grid(D, v) = reshape(v, size(D.glon))
aspect(D) = (D.xlim[2] - D.xlim[1]) / (D.ylim[2] - D.ylim[1])
# the land of the first run's last record (C_ocean < 0.5) on a domain's grid
land(D, oc) = grid(D, [1.0 - PM.interp(oc.C, PM.xweight(oc.lon, x), PM.yweight(oc.lat, y)) for (x, y) in zip(vec(D.glon), vec(D.glat))])
function mapax(pos, D; title = "")
    ax = Axis(pos; aspect = DataAspect(), limits = (D.xlim..., D.ylim...), title, titlesize = 12, titlefont = :regular)
    hidedecorations!(ax); ax
end
base!(ax, D, Lm) = heatmap!(ax, D.xs, D.ys, Float64.(Lm .> 0.5), colormap = [OCEAN, LAND], colorrange = (0, 1))
function field!(ax, D, F, lim)
    heatmap!(ax, D.xs, D.ys, F, colormap = CMAP, colorrange = (-lim, lim), nan_color = :transparent)
end
function coast!(ax, D, Lm)
    contour!(ax, D.xs, D.ys, Lm, levels = [0.5], color = :black, linewidth = 0.6)
    D.plate && return
    for la in -80:10:80
        p = [stereo(D.lon0, D.lat0, lo, la) for lo in -180:1.0:180 if cosd(angdist(D, lo, la)) > 0]
        length(p) > 1 && lines!(ax, first.(p), last.(p), color = (:black, 0.12), linewidth = 0.5)
    end
    for lo in -180:30:150
        p = [stereo(D.lon0, D.lat0, lo, la) for la in -89:1.0:89 if cosd(angdist(D, lo, la)) > 0]
        length(p) > 1 && lines!(ax, first.(p), last.(p), color = (:black, 0.12), linewidth = 0.5)
    end
end
angdist(D, lo, la) = acosd(clamp(sind(D.lat0) * sind(la) + cosd(D.lat0) * cosd(la) * cosd(lo - D.lon0), -1, 1))
function points!(ax, D, lon, lat, v, lim; markersize = 7)
    P = [proj(D, a, b) for (a, b) in zip(lon, lat)]
    scatter!(ax, first.(P), last.(P), color = v, colormap = CMAP, colorrange = (-lim, lim), markersize = markersize,
             strokewidth = 0.5, strokecolor = :black)
end
note!(ax, s) = text!(ax, 0.02, 0.02, text = s, space = :relative, fontsize = 10, align = (:left, :bottom))
# a symmetric colour limit: the q quantile of |v| rounded up to a step of NICE
const NICE = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0]
nicelim(v; q = 0.98) = (x = quantile(abs.(filter(isfinite, v)), q); NICE[min(searchsortedfirst(NICE, x), length(NICE))])
rowlabel!(fig, r, s) = Label(fig[r, 1], replace(s, " | " => " |\n"), rotation = π / 2, tellheight = false, fontsize = 12, font = :bold)
fmt1(x) = @sprintf("%.1f", x)
rms(v) = sqrt(mean(abs2, v))

# --- inputs -------------------------------------------------------------------------------------------------------------------------
out = ARGS[1]; mkpath(out)
runs = [(String(strip(split(a, "=", limit = 2)[1])), String(split(a, "=", limit = 2)[2])) for a in ARGS[2:end]]
isempty(runs) && error("give the runs: <label>=<run folder> ...")
cfg = TOML.parsefile(joinpath(SCORE, "config.toml"))
em, emsrc = EN.error_model(cfg, SCORE)
gd = EN.gnss_data(cfg)
const RSLG = first.(PM.GROUPS)
obs = [o for o in PM.read_indicators(FES) if o.t <= TMAX && o.group in RSLG]
ip = findall(o -> o.type == 0, obs)
A = sort(NCDataset(ds -> Float64.(ds["time"][:]), joinpath(FES, "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc")))
A = A[A .<= TMAX]
gidx = Dict(g => [i for (g2, ii) in gd.idx if g2 == g for i in ii] for (g, _) in GNSS_DOMAINS)
gdom = Dict(g => domain([gd.obs[i].lon for i in gidx[g]], [gd.obs[i].lat for i in gidx[g]]; pole = g == "GNSS Antarctica")
            for (g, _) in GNSS_DOMAINS)
rdom = Dict(g => (lo = [obs[i].lon for i in ip if obs[i].group == g]; la = [obs[i].lat for i in ip if obs[i].group == g];
                 g == "far field" ? plate_domain(lo, la) : domain(lo, la; pole = g == "Antarctica", margin = 300.0)) for g in RSLG)
# the sites per group: those with the most index points, at the mean position of their points
sites = Dict(g => begin
    names = unique(obs[i].site_name for i in ip if obs[i].group == g)
    n = Dict(s => count(i -> obs[i].site_name == s, ip) for s in names)
    top = sort(names, by = s -> -n[s])[1:min(NSITE, length(names))]
    [(name = s, n = n[s], pos = sphere_mean([obs[i].lon for i in ip if obs[i].site_name == s], [obs[i].lat for i in ip if obs[i].site_name == s]))
     for s in top]
end for g in RSLG)

# --- the runs: everything needed from each, read once -------------------------------------------------------------------------------
function read_ocean(dir)
    NCDataset(joinpath(dir, "vilma_out.nc")) do ds
        lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:]); jl = sortperm(lat); il = sortperm(mod.(lon, 360))
        (lon = mod.(lon[il], 360), lat = lat[jl], C = Float64.(ds["C_ocean"][:, :, size(ds["C_ocean"], 3)])[il, jl])
    end
end
sims = map(runs) do (label, folder)
    dir = EN.ref_dir(Dict("run" => folder))
    PM.run_done(dir) || error("$label: $dir has not finished (no status 0)")
    println(label, ": ", dir)
    R = PM.read_rate(dir); R === nothing && error("$label: no dz_bed_rate in $dir")
    rate = PM.rate_at(R, [s.lon for s in gd.obs], [s.lat for s in gd.obs])
    fields = Dict(g => grid(D, PM.rate_at(R, vec(D.glon), vec(D.glat))) for (g, D) in gdom)
    r = PM.read_run(dir)
    all(a -> any(x -> abs(x - PM.trun(r, a)) < 0.01, r.t), A) || error("$label: records missing at the ages of the scores")
    S = PM.series_run(r, [obs[i].lon for i in ip], [obs[i].lat for i in ip], A; rec = a -> PM.trun(r, a))
    sc = [PM.score(obs[i], A, S[j, :]) for (j, i) in enumerate(ip)]
    res = Dict(i => (-s.d, s.d / s.z) for (i, s) in zip(ip, sc) if isfinite(s.z) && s.z != 0)   # obs - model, sigma_i
    k = findall(t -> -t <= TMAX + 5000, r.t)
    curves = Dict(s.name => (W = (PM.xweight(r.lon, s.pos[1]), PM.yweight(r.lat, s.pos[2]));
                             (-r.t[k] ./ 1000, [PM.interp(view(r.rsl, :, :, j), W...) - PM.interp(view(r.rsl, :, :, length(r.t)), W...) for j in k]))
                  for g in RSLG for s in sites[g])
    (label = label, dir = dir, rate = rate, fields = fields, res = res, curves = curves)
end
oc = read_ocean(sims[1].dir)
GC.gc()

io = open(joinpath(out, "maps.md"), "w")
println(io, "# plastic_data_maps.jl\n\nRuns:\n")
for s in sims; println(io, "- ", s.label, ": `", s.dir, "`"); end
println(io, "\nError model: ", emsrc, ". Mask: sd > ", MASK, " of the prior sd.\n")

# --- GNSS ---------------------------------------------------------------------------------------------------------------------------
println(io, "## GNSS\n\n| group | n | data GP: l [km] | sigma_f | mean [mm/yr] | residual GP: sigma_m | l [km] | source |\n|---|---|---|---|---|---|---|---|")
gnss_rows = String[]
for (g, file) in GNSS_DOMAINS
    D = gdom[g]; ii = gidx[g]; st = gd.obs[ii]
    lon = [s.lon for s in st]; lat = [s.lat for s in st]; y = [s.value for s in st]; σ = [s.sd for s in st]
    X = xyz(lon, lat); t0 = zeros(length(ii))
    fd = gp_fit([(X, y, σ)]; fitmean = true)
    m, sd = gp_pred(gp_cond(X, t0, y, σ, fd.ℓ, fd.s, Inf; μ = fd.μ), D.X, 0.0)
    Fdata = grid(D, masked(m, sd, fd.s))
    resid = [y .- s.rate[ii] for s in sims]
    if g in GNSS_FIT
        fr = gp_fit([(X, r, σ) for r in resid]; fitmean = false); σm, ℓr, src = fr.s, fr.ℓ, "fitted here, the runs together"
    else
        σm, ℓr, _ = em[g]; src = "errmodel.toml (data-scores.md)"
    end
    @printf(io, "| %s | %d | %.0f | %.2f | %.2f | %.2f | %.0f | %s |\n", g, length(ii), fd.ℓ, fd.s, fd.μ, σm, ℓr, src)
    for (s, r) in zip(sims, resid)
        push!(gnss_rows, @sprintf("| %s | %s | %d | %+.2f | %.2f |", g, s.label, length(r), mean(r), rms(r)))
    end
    lim = nicelim(y)
    Lm = land(D, oc); w = 330.0; h = w / aspect(D)
    fig = Figure(size = (2w + 120, (length(sims) + 1) * (h + 28) + 40), fontsize = 12)
    ax = mapax(fig[1, 2], D; title = "Data: $(length(ii)) stations over their GP mean")
    base!(ax, D, Lm); field!(ax, D, Fdata, lim); coast!(ax, D, Lm); points!(ax, D, lon, lat, y, lim)
    Label(fig[1, 3], "GP of the data: Matérn-3/2, l = $(round(Int, fd.ℓ)) km,\nσ_f = $(fmt1(fd.s)) mm/yr, mean $(fmt1(fd.μ)) mm/yr.\n" *
          "GP of the residual: σ_m = $(@sprintf("%.2f", σm)) mm/yr,\nl = $(round(Int, ℓr)) km ($(g in GNSS_FIT ? "fitted" : "error model")).\n" *
          "Fields shown where the GP sd\n≤ $(round(Int, 100MASK)) % of its prior sd.", tellwidth = false, tellheight = false,
          justification = :left, halign = :left, fontsize = 11)
    for (k, s) in enumerate(sims)
        rowlabel!(fig, k + 1, s.label)
        ax = mapax(fig[k+1, 2], D; title = k == 1 ? "Model: present-day uplift rate" : "")
        field!(ax, D, s.fields[g], lim); coast!(ax, D, Lm); points!(ax, D, lon, lat, y, lim)
        mr, sdr = gp_pred(gp_cond(X, t0, resid[k], σ, ℓr, σm, Inf), D.X, 0.0)
        ax = mapax(fig[k+1, 3], D; title = k == 1 ? "Residual: observed − model" : "")
        base!(ax, D, Lm); field!(ax, D, grid(D, masked(mr, sdr, σm)), lim); coast!(ax, D, Lm)
        points!(ax, D, lon, lat, resid[k], lim)
        note!(ax, "mean $(@sprintf("%+.1f", mean(resid[k]))), RMS $(fmt1(rms(resid[k]))) mm/yr")
    end
    Colorbar(fig[2:length(sims)+1, 4], limits = (-lim, lim), colormap = CMAP, label = "uplift rate, residual [mm/yr]", height = Relative(0.6))
    colsize!(fig.layout, 1, 36); rowgap!(fig.layout, 1, 4)
    save(joinpath(out, "gnss-$file.png"), fig, px_per_unit = 1.25)
    println("wrote gnss-$file.png")
end
println(io, "\n| group | run | n | mean obs - model [mm/yr] | RMS |\n|---|---|---|---|---|")
foreach(r -> println(io, r), gnss_rows)

# --- RSL maps -----------------------------------------------------------------------------------------------------------------------
println(io, "\n## RSL residuals (obs - model) at the index points\n\n| group | window | run | n | mean [m] | RMS [m] | max abs, points [m] | max abs, field shown [m] |\n|---|---|---|---|---|---|---|---|")
for g in RSLG
    D = rdom[g]; σm, ℓ, Lt = em[g]
    ig = [i for i in ip if obs[i].group == g]
    wins = [(w, a, b) for (w, a, b) in WINDOWS if count(i -> a <= obs[i].t <= b, ig) >= MIN_WINDOW]
    lim = nicelim([s.res[i][1] for s in sims for i in ig if haskey(s.res, i) && any(((_, a, b),) -> a <= obs[i].t <= b, wins)]; q = 0.95)
    Lm = land(D, oc); w = D.plate ? 440.0 : 300.0; h = w / aspect(D)
    fig = Figure(size = (length(wins) * w + 120, length(sims) * (h + 24) + 60), fontsize = 12)
    for (k, s) in enumerate(sims)
        rowlabel!(fig, k, s.label)
        for (c, (wn, a, b)) in enumerate(wins)
            ii = [i for i in ig if a <= obs[i].t <= b && haskey(s.res, i)]
            y = [s.res[i][1] for i in ii]; σ = [s.res[i][2] for i in ii]
            lon = [obs[i].lon for i in ii]; lat = [obs[i].lat for i in ii]; t = [obs[i].t for i in ii]
            m, sd = gp_pred(gp_cond(xyz(lon, lat), t, y, σ, ℓ, σm, Lt), D.X, (a + b) / 2)
            ax = mapax(fig[k, c+1], D; title = k == 1 ? "$wn ($(length(ii)) index points)" : "")
            F = masked(m, sd, σm)
            base!(ax, D, Lm); field!(ax, D, grid(D, F), lim); coast!(ax, D, Lm)
            lm = [o for o in obs if o.group == g && o.type != 0 && a <= o.t <= b]
            isempty(lm) || (P = [proj(D, o.lon, o.lat) for o in lm]; scatter!(ax, first.(P), last.(P), color = :gray45, markersize = 3.5))
            points!(ax, D, lon, lat, y, lim; markersize = 6)
            note!(ax, "mean $(@sprintf("%+.1f", mean(y))), RMS $(fmt1(rms(y))) m")
            @printf(io, "| %s | %s | %s | %d | %+.1f | %.1f | %.1f | %.1f |\n", g, wn, s.label, length(y), mean(y), rms(y),
                    maximum(abs, y), maximum(abs, filter(isfinite, F); init = 0.0))
        end
    end
    Colorbar(fig[1:length(sims), length(wins)+2], limits = (-lim, lim), colormap = CMAP, label = "RSL observed − model [m]",
             height = Relative(0.6))
    Label(fig[length(sims)+1, 2:length(wins)+1], "$g: GP of the residuals, σ_m = $(round(Int, σm)) m, l = $(round(Int, ℓ)) km, " *
          "L_t = $(round(Int, Lt / 1000)) kyr; drawn where its sd ≤ $(round(Int, 100MASK)) % of σ_m. Grey: limiting points.",
          fontsize = 11, tellwidth = false)
    colsize!(fig.layout, 1, 36)
    f = "rsl-" * replace(lowercase(g), " " => "-") * ".png"
    save(joinpath(out, f), fig, px_per_unit = 1.3); println("wrote ", f)
end

# --- RSL at sites -------------------------------------------------------------------------------------------------------------------
for g in RSLG
    ss = sites[g]
    fig = Figure(size = (430 * length(ss), 400), fontsize = 12)
    for (c, st) in enumerate(ss)
        so = [o for o in obs if o.site_name == st.name]
        ipo = [o for o in so if o.type == 0]
        tmax = min(TMAX, maximum(o.tmax for o in ipo)) / 1000 + 1.0
        ax = Axis(fig[1, c], title = "$(replace(st.name, '_' => ' ')) ($(st.n) index points)", xlabel = "age [ka]",
                  ylabel = c == 1 ? "RSL [m]" : "", xreversed = true, titlesize = 12, titlefont = :regular)
        poly!(ax, [Rect2f(o.tmin / 1000, o.rmin, (o.tmax - o.tmin) / 1000, o.rmax - o.rmin) for o in ipo], color = IPCOL,
              strokecolor = IPSTROKE, strokewidth = 0.4)
        lo, hi = extrema([[o.rmin for o in ipo]; [o.rmax for o in ipo]])
        for (ty, mk, col) in LIMITING
            l = [o for o in so if o.type == ty && o.t / 1000 <= tmax]
            isempty(l) && continue
            scatter!(ax, [o.t / 1000 for o in l], [o.rsl for o in l], marker = mk, color = :transparent, strokecolor = col,
                     strokewidth = 1.0, markersize = 7)
            lo = min(lo, minimum(o.rsl for o in l)); hi = max(hi, maximum(o.rsl for o in l))
        end
        for (k, s) in enumerate(sims)
            a, v = s.curves[st.name]; w = a .<= tmax + 2.5   # one record beyond the frame (PaleoMIST: 2.5 kyr)
            lines!(ax, a[w], v[w], color = COLS[k], linewidth = 1.8)
            lo = min(lo, minimum(v[a .<= tmax])); hi = max(hi, maximum(v[a .<= tmax]))
        end
        pad = 0.05 * (hi - lo) + 1
        xlims!(ax, tmax, 0); ylims!(ax, lo - pad, hi + pad)   # age decreasing to the right
    end
    Legend(fig[2, 1:length(ss)],
           [[PolyElement(color = IPCOL, strokecolor = IPSTROKE, strokewidth = 0.4)];
            [MarkerElement(marker = mk, color = :transparent, strokecolor = col, strokewidth = 1.0, markersize = 8) for (_, mk, col) in LIMITING];
            [LineElement(color = COLS[k], linewidth = 1.8) for k in eachindex(sims)]],
           ["index point (age and RSL range)", "marine limiting", "terrestrial limiting", [s.label for s in sims]...],
           orientation = :horizontal, framevisible = false, nbanks = 2)
    f = "sites-" * replace(lowercase(g), " " => "-") * ".png"
    save(joinpath(out, f), fig, px_per_unit = 1.5); println("wrote ", f)
end
println(io, "\n## Sites\n")
for g in RSLG; println(io, "- ", g, ": ", join(["$(s.name) ($(s.n))" for s in sites[g]], ", ")); end
close(io)
println("wrote ", out)
