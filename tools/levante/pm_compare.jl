# PaleoMIST test of the standalone VILMA2: compare the runs of pm_run.sh with PaleoMIST (Gowan et al. 2021) and with each other.
#   julia --project=tools tools/levante/pm_compare.jl <fesmdata> <runs folder> <out folder> [<key>=<folder> ...]
# The runs by key (RUNS: B1 = 1-D Gowan lmax 256, B1_170 = lmax 170, B1_L80, B1_bm, B3, A1, A1_170, A3), each the folder of RUNS unless
# given as key=folder (smoke: B1=S_1d B3=S_3d); a missing run is skipped. Everything on the 0.5 degree PaleoMIST grid: the VILMA2 fields
# (Gauss grid) bilinear in lon-lat.
# Definitions (PaleoMIST, checked): dz_sl = sea-level change from present (SELEN, positive = sea higher than the land), z_bed relative
# to the sea level of the time, so z_bed(t) - z_bed(0) = -dz_sl(t) (to 0.01 m). VILMA2: rsl = dz_ss - dz_bed (same sign as dz_sl),
# z_bed = z_bed_eq - rsl, so the bed change of Gowan's definition is -rsl of both and its metrics are those of rsl with the bias
# negated. dz_bed (u, the solid-Earth radial displacement) has no PaleoMIST counterpart; it is kept for the A - B comparison.
# Record t_k of B is the VILMA2 state at t_k + 1 yr (after the step of the record; pm_prep.jl), the 80 ka record the spun-up start.
# Relative to present: VILMA2's rsl is relative to the relaxed present-day reference (i_eq = 1), and the Earth is not relaxed at 0 ka
# (rsl(0) != 0), while PaleoMIST's dz_sl is S(t) - S(0) by construction; the B fields are therefore taken as f(t) - f(0 ka), as Gowan's.
# A - B (equilibrium - transient at 22.5 ka) compares the two states against the same relaxed reference, without that subtraction.
# The 3-D effect B3 - B1 split (EFFECTS): lithosphere B1_L80 - B1 (80 against 120 km, Gowan's viscosities), radial profile B1_bm - B1_L80
# (the science Earth collapsed to the 1-D profile of its pre_spinup_1d), lateral structure B3 - B1_bm.
# Sea-level indicators (FesmData RSL-GAPSLIP_global, GAPSLIP 2.1, the primary data, alternative = 0, of age <= 80 ka): RSL relative to
# present of PaleoMIST (dz_sl) and of the B runs (rsl(t) - rsl(0 ka)) at the samples, bilinear in space, linear in time between the
# records (2.5 kyr). Index points: d = model - rsl, sigma^2 = sigma_z^2 + (dRSL/dt sigma_t)^2 with sigma_z the 1-sigma (half the 2-sigma)
# elevation uncertainty on the side of d and sigma_t that of the age, dRSL/dt the slope of the model at the age; reported: RMS of d/sigma
# and the fraction with |d/sigma| <= 2. Limiting data: violated if no model RSL over the 2-sigma age range reaches the 2-sigma range of
# the sample (marine: max RSL < rsl_min; terrestrial: min RSL > rsl_max). By the groups of Gowan's regions (GROUPS) and the marine-bed
# areas (FOCUS).
# Writes: pm_metrics.md (tables), pm_metrics_B1_PM.csv (all records), pm_sites.csv, pm_gmsl.csv, pm_fields.nc (for pm_maps.jl),
# pm_indicators.md (scores), pm_ind_samples.csv (per sample), pm_ind_sites.csv (per site), pm_ind_curves.csv (model RSL at the sites).
using NCDatasets, Printf, Statistics

const MAP_AGES = [22500.0, 15000.0, 10000.0]
const SITES = [("Hudson Bay", -85.0, 60.0), ("Gulf of Bothnia", 20.0, 63.0), ("Barents Sea", 40.0, 75.0),
               ("Ross embayment", -175.0, -80.0), ("Weddell embayment", -50.0, -77.0), ("Barbados", -59.6, 13.1)]
const RUNS = ["B1" => "B1_l256", "B1_170" => "B1_l170", "B1_L80" => "B1_L80", "B1_bm" => "B1_bm", "B3" => "B3_l170",
              "A1" => "A1_l256", "A1_170" => "A1_l170", "A3" => "A3_l170"]
const BKEYS = ["B1", "B1_170", "B1_L80", "B1_bm", "B3"]
const EFFECTS = [("B1_170", "B1", "lmax 170 - 256 (B1)"), ("B1_L80", "B1_170", "lithosphere: B1_L80 - B1"),
                 ("B1_bm", "B1_L80", "radial profile: B1_bm - B1_L80"), ("B3", "B1_bm", "lateral: B3 - B1_bm"),
                 ("B3", "B1_170", "3-D - 1-D: B3 - B1 (l170)")]
const GROUPS = ["Laurentide" => ["North_America_Arctic", "North_America_Atlantic"], "Greenland" => ["Greenland"],
                "Fennoscandia/Barents" => ["Europe", "Eurasian_Arctic"], "Antarctica" => ["Antarctica"],
                "far field" => ["Caribbean", "East_Asia", "Southeast_Asia", "South_Asia", "Australia", "Pacific_Islands"]]
# the marine-bed areas where B1 and PaleoMIST differ most: by sector or site of GAPSLIP (no data in the Weddell embayment itself)
const FOCUS = ["Hudson Bay" => (:sector, ["Hudson_Bay"]), "Hudson Strait" => (:sector, ["Hudson_Strait"]),
               "Ross Sea (Scott Coast, Terra Nova Bay)" => (:site_name, ["Southern_Scott_Coast", "Terra_Nova_Bay"]),
               "Weddell side (James Ross Island)" => (:site_name, ["James_Ross_Island"])]
const NCURVE = 3          # sites per group in pm_ind_curves.csv: those with the most index points

# --- VILMA2 output ---------------------------------------------------------------------------------------------------------------
struct Run
    name::String
    lon::Vector{Float64}; lat::Vector{Float64}; t::Vector{Float64}
    rsl::Array{Float64,3}; u::Array{Float64,3}
    bsl::Vector{Float64}; bsl_ive::Vector{Float64}
end

function read_run(dir)
    f = joinpath(dir, "vilma_out.nc")
    isfile(f) || return nothing
    NCDataset(f) do ds
        lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:]); t = Float64.(ds["time"][:])
        jl = sortperm(lat); il = sortperm(mod.(lon, 360))
        g(v) = Float64.(nomissing(ds[v][:, :, :], NaN))[il, jl, :]
        Run(basename(dir), mod.(lon[il], 360), lat[jl], t, g("rsl"), g("dz_bed"),
            Float64.(ds["bsl"][:]), Float64.(ds["bsl_ive"][:]))
    end
end

tindex(r::Run, t) = (k = findfirst(x -> abs(x - t) < 0.01, r.t); k === nothing && error("$(r.name): no output at t = $t"); k)
trec(age) = age == 80000.0 ? -80000.0 : -age + 1.0          # the B state of record age
field(r::Run, v, t) = getfield(r, v)[:, :, tindex(r, t)]
relp(r::Run, v, t) = getfield(r, v)[:, :, tindex(r, t)] .- getfield(r, v)[:, :, end]   # relative to the last step (0 ka + 1 yr)

# bilinear from a grid periodic in lon (ascending in [0, 360)) to target lon-lat; latitude clamped beyond the outermost rows
function xweight(glon, x)
    n = length(glon); xm = mod(x, 360); i = searchsortedlast(glon, xm)
    i0 = i == 0 ? n : i; i1 = i == n || i == 0 ? 1 : i + 1
    x0 = i == 0 ? glon[n] - 360 : glon[i0]; x1 = i == n ? glon[1] + 360 : glon[i1]
    (i0, i1, (xm - x0) / (x1 - x0))
end
function yweight(glat, y)
    j = searchsortedlast(glat, y)
    j == 0 ? (1, 1, 0.0) : j == length(glat) ? (j, j, 0.0) : (j, j + 1, (y - glat[j]) / (glat[j+1] - glat[j]))
end
interp(F, (i0, i1, wx), (j0, j1, wy)) = (1 - wy) * ((1 - wx) * F[i0, j0] + wx * F[i1, j0]) + wy * ((1 - wx) * F[i0, j1] + wx * F[i1, j1])
function bilin(F, glon, glat, tlon, tlat)
    ix = [xweight(glon, x) for x in tlon]; jy = [yweight(glat, y) for y in tlat]
    [interp(F, a, b) for a in ix, b in jy]
end

# --- metrics ----------------------------------------------------------------------------------------------------------------------
function stats(d, w, m)
    k = m .& .!isnan.(d)
    ww = w[k]; dd = d[k]
    (rms = sqrt(sum(ww .* dd .^ 2) / sum(ww)), max = maximum(abs.(dd)), bias = sum(ww .* dd) / sum(ww))
end

# --- sea-level indicators ---------------------------------------------------------------------------------------------------------
# the primary data (alternative = 0) of age 0-80 ka of FesmData RSL-GAPSLIP_global, with the group of their region
function read_indicators(fesmdata)
    grp = Dict(r => g for (g, rs) in GROUPS for r in rs)
    NCDataset(joinpath(fesmdata, "Records", "RSL", "RSL-GAPSLIP_global.nc")) do ds
        v(n) = ds[n][:]
        t = Float64.(v("time")); keep = findall((v("alternative") .== 0) .& (0 .<= t .<= 80000))
        reg = v("region"); sec = v("sector"); site = v("site_name"); la = Float64.(v("lat")); lo = Float64.(v("lon"))
        t0 = Float64.(v("time_min")); t1 = Float64.(v("time_max")); r = Float64.(v("rsl")); r0 = Float64.(v("rsl_min"))
        r1 = Float64.(v("rsl_max")); ty = Int.(v("type"))
        [(region = reg[i], group = grp[reg[i]], sector = sec[i], site_name = site[i], lat = la[i], lon = lo[i], t = t[i],
          tmin = max(t0[i], 0.0), tmax = min(t1[i], 80000.0), st = (t1[i] - t0[i]) / 4, rsl = r[i], rmin = r0[i], rmax = r1[i], type = ty[i]) for i in keep]
    end
end

# RSL relative to present at points (plon, plat), at the ages A (ascending), as an (npoints, nages) matrix: PaleoMIST dz_sl, a B run
function series_pm(lon, lat, ages, sl, plon, plat, A)
    gl = mod.(lon, 360); p = sortperm(gl); q = sortperm(lat); gl = gl[p]; gt = lat[q]
    W = [(xweight(gl, x), yweight(gt, y)) for (x, y) in zip(plon, plat)]
    S = zeros(length(plon), length(A))
    for (k, a) in enumerate(A)
        F = sl[p, q, findfirst(==(a), ages)]
        for (i, (wx, wy)) in enumerate(W); S[i, k] = interp(F, wx, wy); end
    end
    S
end
function series_run(r::Run, plon, plat, A)
    W = [(xweight(r.lon, x), yweight(r.lat, y)) for (x, y) in zip(plon, plat)]
    S = zeros(length(plon), length(A))
    for (k, a) in enumerate(A)
        F = relp(r, :rsl, trec(a))
        for (i, (wx, wy)) in enumerate(W); S[i, k] = interp(F, wx, wy); end
    end
    S
end

# linear in time: the value and the slope at t; the range over [t0, t1]
function tval(A, s, t)
    k = clamp(searchsortedlast(A, t), 1, length(A) - 1); w = (t - A[k]) / (A[k+1] - A[k])
    (1 - w) * s[k] + w * s[k+1], (s[k+1] - s[k]) / (A[k+1] - A[k])
end
function trange(A, s, t0, t1)
    v = [tval(A, s, t0)[1], tval(A, s, t1)[1], (s[k] for k in eachindex(A) if t0 < A[k] < t1)...]
    extrema(v)
end

# score of one sample: index point -> (d, z), limiting -> violated (as 1.0 / 0.0)
function score(o, A, s)
    any(isnan, s) && return (m = NaN, d = NaN, z = NaN, viol = NaN)
    if o.type == 0
        m, dmdt = tval(A, s, o.t)
        d = m - o.rsl; sz = (d > 0 ? o.rmax - o.rsl : o.rsl - o.rmin) / 2          # 1 sigma: half the 2-sigma range, on the side of d
        return (m = m, d = d, z = d / sqrt(sz^2 + (dmdt * o.st)^2), viol = NaN)
    end
    lo, hi = trange(A, s, o.tmin, o.tmax)
    (m = tval(A, s, o.t)[1], d = NaN, z = NaN, viol = o.type == -1 ? Float64(hi < o.rmin) : Float64(lo > o.rmax))
end

function summary(sc, obs, idx)
    ip = [i for i in idx if obs[i].type == 0 && isfinite(sc[i].z)]
    ml = [i for i in idx if obs[i].type == -1 && isfinite(sc[i].viol)]; tl = [i for i in idx if obs[i].type == 1 && isfinite(sc[i].viol)]
    z = [sc[i].z for i in ip]; d = [sc[i].d for i in ip]
    (nip = length(ip), rmsz = isempty(z) ? NaN : sqrt(mean(z .^ 2)), medz = isempty(z) ? NaN : median(abs.(z)),
     in2 = isempty(z) ? NaN : mean(abs.(z) .<= 2),
     bias = isempty(d) ? NaN : mean(d), nml = length(ml), vml = isempty(ml) ? NaN : mean(sc[i].viol for i in ml),
     ntl = length(tl), vtl = isempty(tl) ? NaN : mean(sc[i].viol for i in tl))
end

function indicators(fesmdata, out, R, lon, lat, ages, sl)
    obs = read_indicators(fesmdata)
    A = sort(ages)
    models = [("PaleoMIST", series_pm(lon, lat, ages, sl, [o.lon for o in obs], [o.lat for o in obs], A))]
    for key in BKEYS
        haskey(R, key) && all(a -> any(x -> abs(x - trec(a)) < 0.01, R[key].t), A) || continue
        push!(models, (R[key].name, series_run(R[key], [o.lon for o in obs], [o.lat for o in obs], A)))
    end
    scores = [[score(o, A, S[i, :]) for (i, o) in enumerate(obs)] for (_, S) in models]
    println("indicators: ", length(obs), " samples (", count(o -> o.type == 0, obs), " index points), models ", join(first.(models), ", "))

    sets = [("all", collect(eachindex(obs))); [(g, findall(o -> o.group == g, obs)) for (g, _) in GROUPS]]
    focus = [(f, findall(o -> getfield(o, k) in v, obs)) for (f, (k, v)) in FOCUS]
    fmt(x) = isnan(x) ? "-" : @sprintf("%.2f", x)
    open(joinpath(out, "pm_indicators.md"), "w") do io
        println(io, "Sea-level indicators (GAPSLIP 2.1, primary data, 0-80 ka): index points (IP), RMS of the normalized misfit z = (model - obs)/sigma,",
                " median |z|, fraction |z| <= 2, mean misfit [m]; limiting data, fraction violated (marine ML, terrestrial TL).\n")
        for (title, S) in (("by group", sets), ("marine-bed areas", focus))
            println(io, "\n### ", title, "\n")
            println(io, "| set | model | n IP | RMS z | median abs(z) | abs(z) <= 2 | mean d [m] | n ML | ML violated | n TL | TL violated | ages [ka] |")
            println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|")
            for (name, idx) in S
                ar = isempty(idx) ? "-" : @sprintf("%.1f-%.1f", minimum(obs[i].t for i in idx) / 1000, maximum(obs[i].t for i in idx) / 1000)
                for (m, (mn, _)) in enumerate(models)
                    s = summary(scores[m], obs, idx)
                    println(io, "| ", m == 1 ? name : "", " | ", mn, " | ", s.nip, " | ", fmt(s.rmsz), " | ", fmt(s.medz), " | ", fmt(s.in2), " | ",
                            isnan(s.bias) ? "-" : @sprintf("%.1f", s.bias), " | ", s.nml, " | ", fmt(s.vml), " | ", s.ntl, " | ", fmt(s.vtl),
                            " | ", m == 1 ? ar : "", " |")
                end
            end
        end
    end
    open(joinpath(out, "pm_ind_samples.csv"), "w") do c
        println(c, "region,group,sector,site,lat,lon,time,time_min,time_max,rsl,rsl_min,rsl_max,type,",
                join(["$(m)_rsl,$(m)_score" for (m, _) in models], ","))
        for (i, o) in enumerate(obs)
            print(c, join((o.region, o.group, o.sector, o.site_name, o.lat, o.lon, o.t, o.tmin, o.tmax, o.rsl, o.rmin, o.rmax, o.type), ","))
            for m in eachindex(models)
                s = scores[m][i]; print(c, ",", @sprintf("%.3f", s.m), ",", @sprintf("%.3f", o.type == 0 ? s.z : s.viol))
            end
            println(c)
        end
    end
    sites = unique(o.site_name for o in obs)
    idxs = Dict(st => findall(o -> o.site_name == st, obs) for st in sites)
    pos(st) = (i = idxs[st]; (mean(obs[k].lat for k in i), rad2deg(atan(mean(sind(obs[k].lon) for k in i), mean(cosd(obs[k].lon) for k in i)))))
    open(joinpath(out, "pm_ind_sites.csv"), "w") do c
        println(c, "model,site,group,lat,lon,n_ip,rms_z,in2,n_lim,lim_violated")
        for (m, (mn, _)) in enumerate(models), st in sites
            s = summary(scores[m], obs, idxs[st]); (la, lo) = pos(st)
            println(c, join((mn, st, obs[idxs[st][1]].group, la, lo, s.nip, s.rmsz, s.in2, s.nml + s.ntl,
                             (s.nml * nan0(s.vml) + s.ntl * nan0(s.vtl)) / max(s.nml + s.ntl, 1)), ","))
        end
    end
    # model RSL at the best-constrained sites (most index points) of each group, and the focus sites with data
    nip(st) = count(k -> obs[k].type == 0, idxs[st])
    pick = String[]
    for (g, _) in GROUPS
        append!(pick, first(sort(filter(st -> obs[idxs[st][1]].group == g, sites); by = st -> -nip(st)), NCURVE))
    end
    for (_, idx) in focus
        isempty(idx) || push!(pick, first(sort(unique(obs[i].site_name for i in idx); by = st -> -length(idxs[st]))))
    end
    pick = unique(pick)
    plon = [pos(st)[2] for st in pick]; plat = [pos(st)[1] for st in pick]
    curves = [("PaleoMIST", series_pm(lon, lat, ages, sl, plon, plat, A));
              [(R[k].name, series_run(R[k], plon, plat, A)) for k in BKEYS if haskey(R, k) && any(m -> m[1] == R[k].name, models)]]
    open(joinpath(out, "pm_ind_curves.csv"), "w") do c
        println(c, "source,site,group,age,rsl")
        for (mn, S) in curves, (i, st) in enumerate(pick), (k, a) in enumerate(A)
            println(c, mn, ",", st, ",", obs[idxs[st][1]].group, ",", a, ",", S[i, k])
        end
    end
end
nan0(x) = isnan(x) ? 0.0 : x

function main(args)
    fesmdata, runs, out = args[1], args[2], args[3]
    folders = Dict(RUNS); foreach(a -> (kv = split(a, '='); folders[kv[1]] = kv[2]), args[4:end])
    mkpath(out)
    R = Dict{String,Run}()
    for (key, _) in RUNS
        nm = folders[key]
        r = read_run(joinpath(runs, nm)); r === nothing ? println("skip ", nm) : (R[key] = r; println("read ", nm, ": ", length(r.t), " steps"))
    end

    pm = NCDataset(joinpath(fesmdata, "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc"))
    lon = Float64.(pm["lon"][:]); lat = Float64.(pm["lat"][:]); ages = Float64.(pm["time"][:])
    sl = Float64.(nomissing(pm["dz_sl"][:, :, :], NaN)); H = Float64.(nomissing(pm["H_ice"][:, :, :], 0.0))
    zb = Float64.(nomissing(pm["z_bed"][:, :, :], NaN))
    close(pm)
    w = [cosd(lat[j]) for i in eachindex(lon), j in eachindex(lat)]
    glac = (dropdims(maximum(H; dims = 3); dims = 3) .> 0) .& (H[:, :, end] .== 0)
    ocean = (zb[:, :, end] .< 0) .& (H[:, :, end] .== 0) .& .!glac
    allc = trues(size(w))
    regions = [("glaciated", glac), ("ocean", ocean), ("global", allc)]
    kage(a) = findfirst(==(a), ages)
    on05(r, v, t) = bilin(field(r, v, t), r.lon, r.lat, lon, lat)        # relative to the relaxed reference
    rel05(r, v, t) = bilin(relp(r, v, t), r.lon, r.lat, lon, lat)        # relative to present

    io = open(joinpath(out, "pm_metrics.md"), "w")
    println(io, "Area-weighted metrics on the 0.5 degree grid [m]: RMS, max |d|, bias = mean(d); d = first - second. Regions: glaciated = ice at",
            " some record, ice-free today; ocean = present ocean outside it; global.\n")
    row(io, lbl, a, d) = (print(io, "| ", lbl, " | ", @sprintf("%.1f", a / 1000), " |");
                          for (_, m) in regions; s = stats(d, w, m); print(io, @sprintf(" %.1f | %.1f | %.1f |", s.rms, s.max, s.bias)); end; println(io))
    hdr(io, title) = (println(io, "\n### ", title, "\n"),
                      println(io, "| comparison | ka | glac RMS | glac max | glac bias | ocean RMS | ocean max | ocean bias | global RMS | global max | global bias |"),
                      println(io, "|---|---|---|---|---|---|---|---|---|---|---|"))

    # 1. B vs PaleoMIST, rsl (Gowan's bed change = -rsl: same RMS/max, bias negated)
    hdr(io, "RSL relative to present, B runs - PaleoMIST")
    for key in BKEYS, a in MAP_AGES
        haskey(R, key) || continue
        any(x -> abs(x - trec(a)) < 0.01, R[key].t) || continue
        row(io, "$(R[key].name) - PM", a, rel05(R[key], :rsl, trec(a)) .- sl[:, :, kage(a)])
    end
    if haskey(R, "B1")
        open(joinpath(out, "pm_metrics_B1_PM.csv"), "w") do c
            println(c, "age,", join(["$(n)_$(s)" for (n, _) in regions for s in ("rms", "max", "bias")], ","))
            for a in ages
                any(x -> abs(x - trec(a)) < 0.01, R["B1"].t) || continue
                d = rel05(R["B1"], :rsl, trec(a)) .- sl[:, :, kage(a)]
                println(c, a, ",", join([@sprintf("%.2f", getfield(stats(d, w, m), s)) for (_, m) in regions for s in (:rms, :max, :bias)], ","))
            end
        end
    end
    # 2. resolution and 3-D effects
    hdr(io, "RSL: resolution (lmax 170 - 256), and the 3-D effect (B3 - B1, lmax 170) by lithosphere, radial profile and lateral structure")
    for (k1, k2, _) in EFFECTS, a in MAP_AGES
        (haskey(R, k1) && haskey(R, k2)) || continue
        any(x -> abs(x - trec(a)) < 0.01, R[k1].t) || continue
        row(io, "$(R[k1].name) - $(R[k2].name)", a, rel05(R[k1], :rsl, trec(a)) .- rel05(R[k2], :rsl, trec(a)))
    end
    # 3. steady state: A (equilibrium under the 22.5 ka ice) - B (transient) at 22.5 ka
    hdr(io, "Equilibrium - transient at 22.5 ka (A - B): RSL, and u = dz_bed (solid-Earth displacement)")
    for (ka, kb) in (("A1", "B1"), ("A1_170", "B1_170"), ("A3", "B3"))
        (haskey(R, ka) && haskey(R, kb)) || continue
        tA = R[ka].t[end]
        row(io, "$(R[ka].name) - $(R[kb].name), rsl", 22500.0, on05(R[ka], :rsl, tA) .- on05(R[kb], :rsl, trec(22500.0)))
        row(io, "$(R[ka].name) - $(R[kb].name), u", 22500.0, on05(R[ka], :u, tA) .- on05(R[kb], :u, trec(22500.0)))
    end
    # the present-day disequilibrium of B: rsl and u at 0 ka against the relaxed reference
    hdr(io, "Present-day disequilibrium of B (0 ka state against the relaxed reference): RSL and u")
    for key in BKEYS
        haskey(R, key) || continue
        row(io, "$(R[key].name), rsl(0)", 0.0, on05(R[key], :rsl, R[key].t[end]))
        row(io, "$(R[key].name), u(0)", 0.0, on05(R[key], :u, R[key].t[end]))
    end
    close(io)

    # sites: rsl of every output step (B runs) and of the records (PaleoMIST); A at 22.5 ka
    open(joinpath(out, "pm_sites.csv"), "w") do c
        println(c, "source,site,time,rsl")
        for (nm, slon, slat) in SITES
            i = argmin(abs.(mod.(lon .- slon .+ 180, 360) .- 180)); j = argmin(abs.(lat .- slat))
            for (k, a) in enumerate(ages); println(c, "PaleoMIST,", nm, ",", -a, ",", sl[i, j, k]); end
            for key in BKEYS
                haskey(R, key) || continue
                r = R[key]
                for (k, t) in enumerate(r.t)
                    v = bilin(r.rsl[:, :, k] .- r.rsl[:, :, end], r.lon, r.lat, [lon[i]], [lat[j]])[1]
                    println(c, r.name, ",", nm, ",", t, ",", v)
                end
            end
        end
    end
    # global mean sea level: VILMA2 bsl (ocean mean of C rsl) and bsl_ive; PaleoMIST ocean mean of dz_sl over the ocean of the time
    # (z_bed < 0, no grounded ice) and the ESL series of the ice volumes
    esl = NCDataset(joinpath(fesmdata, "Records", "SeaLevel", "SeaLevel-PaleoMIST.nc")) do ds
        (Float64.(ds["time"][:]), Float64.(ds["sl"][:, 1]))
    end
    open(joinpath(out, "pm_gmsl.csv"), "w") do c
        println(c, "source,time,value")
        for (k, a) in enumerate(ages)
            m = (zb[:, :, k] .< 0) .& (H[:, :, k] .== 0)
            println(c, "PaleoMIST ocean mean dz_sl,", -a, ",", sum(w[m] .* sl[:, :, k][m]) / sum(w[m]))
        end
        for (a, v) in zip(esl...); println(c, "PaleoMIST ESL,", -a, ",", v); end
        for key in BKEYS
            haskey(R, key) || continue
            b0 = R[key].bsl[end]
            for (t, b, bi) in zip(R[key].t, R[key].bsl, R[key].bsl_ive)
                println(c, R[key].name, " bsl - bsl(0),", t, ",", b - b0); println(c, R[key].name, " bsl_ive,", t, ",", bi)
            end
        end
    end

    # fields for the maps
    f = joinpath(out, "pm_fields.nc"); isfile(f) && rm(f)
    NCDataset(f, "c") do nc
        defDim(nc, "lon", length(lon)); defDim(nc, "lat", length(lat)); defDim(nc, "rec", length(MAP_AGES))
        defVar(nc, "lon", lon, ("lon",)); defVar(nc, "lat", lat, ("lat",)); defVar(nc, "age", MAP_AGES, ("rec",))
        defVar(nc, "glac", Int8.(glac), ("lon", "lat")); defVar(nc, "ocean", Int8.(ocean), ("lon", "lat"))
        defVar(nc, "land", Int8.(zb[:, :, end] .>= 0), ("lon", "lat"))
        defVar(nc, "H_lgm", H[:, :, kage(22500.0)], ("lon", "lat"))
        defVar(nc, "rsl_PM", cat([sl[:, :, kage(a)] for a in MAP_AGES]...; dims = 3), ("lon", "lat", "rec"))
        for key in BKEYS
            haskey(R, key) || continue
            all(a -> any(x -> abs(x - trec(a)) < 0.01, R[key].t), MAP_AGES) || continue
            defVar(nc, "rsl_" * key, cat([rel05(R[key], :rsl, trec(a)) for a in MAP_AGES]...; dims = 3), ("lon", "lat", "rec"))
            defVar(nc, "rsl0_" * key, on05(R[key], :rsl, R[key].t[end]), ("lon", "lat"))
        end
        for (ka, kb) in (("A1", "B1"), ("A1_170", "B1_170"), ("A3", "B3"))
            (haskey(R, ka) && haskey(R, kb)) || continue
            for v in (:rsl, :u)
                defVar(nc, "AmB_$(v)_" * kb, on05(R[ka], v, R[ka].t[end]) .- on05(R[kb], v, trec(22500.0)), ("lon", "lat"))
            end
        end
        nc.attrib["runs"] = join(["$k=$(R[k].name)" for k in sort(collect(keys(R)))], " ")
    end
    indicators(fesmdata, out, R, lon, lat, ages, sl)
    println("wrote ", out)
end
main(ARGS)
