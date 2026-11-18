# PaleoMIST test of the standalone VILMA2: compare the standalone runs (the cases of paleomist/, pm_project.sh) with PaleoMIST (Gowan et al. 2021) and with each other.
#   julia --project=tools/levante tools/levante/pm_compare.jl <fesmdata> <runs folder> <out folder> [<key>=<folder> ...]
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
# GNSS uplift rates (the GNSS groups of plastic_ensemble.jl): read_gnss (the stations of a group, FesmData Records/GPS), read_rate and
# rate_at (the model's present-day uplift rate, dz_bed_rate of the run's last VILMA2 record, at the stations). Ice-core elevations (the
# ice-core group): read_cores (the data, FesmData Records/IceCores) and core_at (the model's surface elevation change at the sites).
# Earth rotation (the rotation group): read_rotation (the data of the config) and rotation_rates (J2dot and the true polar wander
# velocity).
# Included by plastic_esmda.jl (in a module) for the scoring functions: main runs only as a script.
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
                "Fennoscandia" => ["Europe"], "Barents" => ["Eurasian_Arctic"], "Antarctica" => ["Antarctica"],
                "far field" => ["Caribbean", "East_Asia", "Southeast_Asia", "South_Asia", "Australia", "Pacific_Islands"]]
# groups made of groups (the configs before the split of Fennoscandia/Barents, 2026-11-12: the Fennoscandian shield and the softer,
# thinner Barents shelf want different Earths); in_group: an index point is in a group or in a union that holds its group
const GROUP_UNIONS = Dict("Fennoscandia/Barents" => ["Fennoscandia", "Barents"])
in_group(o, g) = o.group == g || o.group in get(GROUP_UNIONS, g, String[])
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

# A standalone run (a case of paleomist/: vilma_out.nc) or a kryos run (geo = vilma2: the VILMA2 stream <inst>_vilma.nc of its global
# instance, read in place). The time of a kryos record is moved by +1 yr except the first (trec: record t_k of a standalone run is
# the state at t_k + 1 yr, its 80 ka record the spun-up start; a kryos record at t_k is the state after the step ending at t_k).
function read_run(dir; inst = "geo_glob")
    f = joinpath(dir, "vilma_out.nc"); kryos = !isfile(f)
    kryos && (f = joinpath(dir, inst * "_vilma.nc"))
    isfile(f) || return nothing
    NCDataset(f) do ds
        lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:]); t = Float64.(ds["time"][:])
        kryos && (t[2:end] .+= 1.0)
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
# rec: the record time of age a (trec: the B state of a step history; trun: the record at the age itself where the run has no B
# state, as gia_refcorr's runs of a linear history)
function series_run(r::Run, plon, plat, A; rec = trec)
    W = [(xweight(r.lon, x), yweight(r.lat, y)) for (x, y) in zip(plon, plat)]
    S = zeros(length(plon), length(A))
    for (k, a) in enumerate(A)
        F = relp(r, :rsl, rec(a))
        for (i, (wx, wy)) in enumerate(W); S[i, k] = interp(F, wx, wy); end
    end
    S
end

# RSL relative to present from a topography relative to the sea level of the time (FesmData IceHistory-GLAC1D: z_topo, the surface of
# land and grounded ice, the sea floor and the base of floating ice; z_srf, the surface; H_ice), the reconstruction's own GIA: the
# bed z_b = z_srf - H_ice under grounded ice (z_topo = z_srf), z_topo elsewhere, and RSL(t) = z_b(0 ka) - z_b(t). Returns lon, lat,
# ages and the RSL (lon, lat, age), as the dz_sl of PaleoMIST (series_pm).
function rsl_topo(file)
    NCDataset(file) do ds
        ages = Float64.(ds["time"][:]); k0 = findfirst(==(0.0), ages)
        k0 === nothing && error("$(file): no 0 ka record")
        g(v) = Float64.(nomissing(ds[v][:, :, :], NaN))
        zt = g("z_topo"); zs = g("z_srf"); H = g("H_ice")
        zb = @. ifelse((H > 0) & (abs(zs - zt) < 1.0), zs - H, zt)
        Float64.(ds["lon"][:]), Float64.(ds["lat"][:]), ages, zb[:, :, k0:k0] .- zb
    end
end

# --- GNSS uplift rates ----------------------------------------------------------------------------------------------------------
# The stations of a GNSS group (a [[gnss]] entry of an ensemble config: plastic_ensemble.jl), from FesmData Records/GPS/GPS-<product>.nc
# in dir (uplift_rate, uplift_rate_sd [mm/yr], optional uplift_rate_elastic and its sd; the attribute elastic_correction of uplift_rate
# says whether the product removed the elastic response to present-day ice change). rate (of the entry):
#   gia        uplift_rate - uplift_rate_elastic, the GIA estimate of a product that did not remove it (sd combined with the elastic
#              part's where the product has one)
#   corrected  uplift_rate of a product that removed it (default)
#   observed   uplift_rate as it is
# box, exclude (optional): lists of [lon_min, lon_max, lat_min, lat_max] [degrees east, north; lon_min to lon_max eastward]: the
# stations in a box (all without box) and in no exclude box. min_years (optional): stations whose series (time_start to time_end)
# spans at least that many years. Stations without a finite rate and a positive sd are left out.
# Returns per station (group, kind = :gnss, site_name, lat, lon, t = 0 (the present: no time term in the error covariance), value (the
# rate), sd) [mm/yr].
function read_gnss(dir, spec)
    nm = spec["name"]; f = joinpath(dir, "GPS-" * spec["product"] * ".nc")
    isfile(f) || error("[[gnss]] $nm: no $f ([score] gnss_dir)")
    kind = get(spec, "rate", "corrected")
    inbox(b, x, y) = b[3] <= y <= b[4] && (b[2] - b[1] >= 360 || mod(x - b[1], 360) <= mod(b[2] - b[1], 360))
    boxes = get(spec, "box", Any[]); excl = get(spec, "exclude", Any[])
    NCDataset(f) do ds
        ec = get(ds["uplift_rate"].attrib, "elastic_correction", "")
        v(n) = Float64.(nomissing(ds[n][:], NaN))
        r = v("uplift_rate"); sd = v("uplift_rate_sd")
        if kind == "gia"
            ec == "not removed" && haskey(ds, "uplift_rate_elastic") ||
                error("[[gnss]] $nm: rate = \"gia\" needs a product with uplift_rate_elastic, the elastic part not removed ($f: $ec)")
            r = r .- v("uplift_rate_elastic")
            haskey(ds, "uplift_rate_elastic_sd") && (sd = sqrt.(sd .^ 2 .+ v("uplift_rate_elastic_sd") .^ 2))
        elseif kind == "corrected"
            ec == "removed" || error("[[gnss]] $nm: rate = \"corrected\" needs a product with the elastic part removed ($f: $ec)")
        elseif kind != "observed"
            error("[[gnss]] $nm: rate = \"$kind\": gia | corrected | observed")
        end
        lat = v("lat"); lon = v("lon"); site = string.(ds["site_name"][:])
        span = haskey(spec, "min_years") ? (haskey(ds, "time_start") ? v("time_end") .- v("time_start") :
                                             error("[[gnss]] $nm: min_years, but $f has no time_start")) : fill(Inf, length(r))
        sel = [i for i in eachindex(r) if isfinite(r[i]) && isfinite(sd[i]) && sd[i] > 0 && span[i] >= get(spec, "min_years", 0.0) &&
               (isempty(boxes) || any(b -> inbox(b, lon[i], lat[i]), boxes)) && !any(b -> inbox(b, lon[i], lat[i]), excl)]
        [(group = nm, kind = :gnss, site_name = site[i], lat = lat[i], lon = lon[i], t = 0.0, value = r[i], sd = sd[i]) for i in sel]
    end
end

# The model's present-day uplift rate of a run: dz_bed_rate [m/yr] of the last record of its VILMA2 stream (vilma_out.nc of a standalone
# run, <inst>_vilma.nc of a kryos run), the rate of the bed over VILMA2's last update; nothing without the file or the variable (VILMA2
# before dz_bed_rate). Grid as read_run's.
function read_rate(dir; inst = "geo_glob")
    f = joinpath(dir, "vilma_out.nc"); isfile(f) || (f = joinpath(dir, inst * "_vilma.nc"))
    isfile(f) || return nothing
    NCDataset(f) do ds
        haskey(ds, "dz_bed_rate") || return nothing
        lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:]); jl = sortperm(lat); il = sortperm(mod.(lon, 360))
        n = size(ds["dz_bed_rate"], 3)
        (lon = mod.(lon[il], 360), lat = lat[jl], rate = Float64.(nomissing(ds["dz_bed_rate"][:, :, n], NaN))[il, jl])
    end
end
# the rate at points (plon, plat) [mm/yr], bilinear
rate_at(R, plon, plat) = [1000.0 * interp(R.rate, xweight(R.lon, x), yweight(R.lat, y)) for (x, y) in zip(plon, plat)]

# --- ESL --------------------------------------------------------------------------------------------------------------------------
# linear in time: the series x on the ascending axis t at a, its record where a is one; NaN outside [t[1], t[end]]
function tlin(t, x, a)
    (a < t[1] || a > t[end]) && return NaN
    k = clamp(searchsortedlast(t, a), 1, length(t) - 1); w = (a - t[k]) / (t[k+1] - t[k])
    w == 0 ? x[k] : w == 1 ? x[k+1] : (1 - w) * x[k] + w * x[k+1]
end
# A global sea-level curve of FesmData Records/SeaLevel/SeaLevel-<product>.nc (sl or the variable var, negative when lower; region 1),
# relative to its own present (the value at 0 ka subtracted: Spratt and Lisiecki 2016 are +8.5 m there), at the ages [yr], linear
# in time between its records (Spratt and Lisiecki: every 1 kyr, so 22.5 ka is the mean of 22 and 23 ka); with its 1 sigma (var_sd,
# linear the same way; NaN without one). The time of the file in years ago.
function sealevel_at(fesmdata, product, ages; var = "sl")
    NCDataset(joinpath(fesmdata, "Records", "SeaLevel", "SeaLevel-" * product * ".nc")) do ds
        t = Float64.(ds["time"][:]); x = Float64.(nomissing(ds[var][:], NaN)); x = ndims(x) > 1 ? vec(x[1, :]) : x
        sd = haskey(ds, var * "_sd") ? (y = Float64.(nomissing(ds[var * "_sd"][:], NaN)); ndims(y) > 1 ? vec(y[1, :]) : y) : fill(NaN, length(t))
        p = sortperm(t); t = t[p]; x = x[p]; sd = sd[p]
        x0 = tlin(t, x, 0.0)
        [tlin(t, x, a) - x0 for a in ages], [tlin(t, sd, a) for a in ages]
    end
end
# The points of the ESL group ([esl] of an ensemble config: name, product, var, ages_ka): the curve at the ages, each a point at
# (0, 0) with t its age, so that the error covariance has only its time term. Returns (group, kind = :esl, site_name, lat, lon, t,
# value [m], sd [m]).
function read_esl(fesmdata, spec)
    ages = 1000.0 .* Float64.(get(spec, "ages_ka", [25.0, 22.5, 20.0, 17.5, 15.0, 12.5, 10.0]))
    v, sd = sealevel_at(fesmdata, spec["product"], ages; var = get(spec, "var", "sl"))
    any(x -> !(isfinite(x) && x > 0), sd) && error("[esl]: $(spec["product"]) has no positive sd at every age")
    [(group = get(spec, "name", "ESL"), kind = :esl, site_name = @sprintf("ESL %g ka", a / 1000), lat = 0.0, lon = 0.0, t = a,
      value = v[i], sd = sd[i]) for (i, a) in enumerate(ages)]
end
# whether the run of a folder has finished without error: its status file (written at the end: the exit code and the run time, by
# runme's job and by gia_refcorr.sh) says 0. VILMA2 writes the whole time axis of its output when it starts, so the output alone does
# not tell a finished run from one still going (or stopped), whose later records are fill values.
run_done(dir) = (f = joinpath(dir, "status"); isfile(f) && (w = split(read(f, String)); !isempty(w) && w[1] == "0"))
# the age [yr] of each record of a run: its time rounded to 10 yr, without the 1 yr of trec / trun (records 100 yr apart or more)
rec_age(r::Run) = -10.0 .* round.(r.t ./ 10.0)
# The model's ESL: a scalar of VILMA2's own diagnostics (vilma_coupling.f90 update_bsl; in a kryos run geo_glob_vilma.nc of the
# global instance), var of esl_at, [esl] model_var of an ensemble config:
#   bsl      the barystatic sea level, the ocean mean of C rsl, ∫C rsl dA / ∫C dA over the ocean of the time; by the mass balance of
#            the sea-level equation the melt water (−ρ_i/ρ_w ∫ΔI_g dA, ΔI_g the grounded-ice change) plus the flooded-area fill
#            (vilma_coupling.f90 fill_int, with the host's F in a kryos run). The ocean that fills newly flooded area (the
#            cavity of lost marine ice, flooding shelves) does not raise it, so in effect the ice above flotation counts, with the
#            change of the ocean area: the ice-volume-equivalent sea level of a stack. The default.
#   bsl_ive  −ρ_i/ρ_w ∫ΔI_g dA / ∫C dA, the whole grounded column, below flotation too (about 10 m more ice at the LGM; bsl −
#            bsl_ive = fill_int / ∫C dA): the comparator before 2026-11-18, kept to reproduce the old scores
const ESL_VARS = (:bsl, :bsl_ive)
# a run's ESL [m] at the ages [yr]: var (ESL_VARS) relative to its last record (0 ka), linear in time between the records (rec_age; the
# record itself at an age it has); NaN for a run that has not finished
function esl_at(dir, ages; var = :bsl)
    var in ESL_VARS || error("esl_at: var $(var), not one of $(ESL_VARS)")
    run_done(dir) || return fill(NaN, length(ages))
    r = read_run(dir); r === nothing && return fill(NaN, length(ages))
    b = getfield(r, var); a = rec_age(r); p = sortperm(a); x = b .- b[end]
    [tlin(a[p], x[p], age) for age in ages]
end

# --- ice-core elevations ----------------------------------------------------------------------------------------------------------
# The points of the ice-core group ([cores] of an ensemble config: plastic_ensemble.jl): the surface elevation change relative to the
# present at the sites of FesmData Records/IceCores/IceCores-<product>.nc (dz_srf, dz_srf_sd [m] by station and time, time in years
# ago; Vinther2009: the elevation above sea level, the bed's isostatic change included) at the ages (ages_ka), each the mean of the
# records within window_kyr / 2 of the age (0.5 kyr: 25 of the 20-yr values of Vinther et al. 2009), with the mean of their sd; sites
# (optional): the site names to use (default all). Returns per site and age (group, kind = :core, site_name "<site> <age> ka", lat,
# lon, t = the age, value, sd) [m].
function read_cores(fesmdata, spec)
    nm = String(get(spec, "name", "ice cores")); f = joinpath(fesmdata, "Records", "IceCores", "IceCores-" * spec["product"] * ".nc")
    isfile(f) || error("[cores]: no $f")
    ages = 1000.0 .* Float64.(spec["ages_ka"]); hw = 500.0 * Float64(get(spec, "window_kyr", 0.5))
    NCDataset(f) do ds
        site = string.(ds["site_name"][:]); lat = Float64.(ds["lat"][:]); lon = Float64.(ds["lon"][:]); t = Float64.(ds["time"][:])
        g(v) = (x = Float64.(nomissing(ds[v][:, :], NaN)); dimnames(ds[v])[1] == "station" ? x : permutedims(x))    # (station, time)
        dz = g("dz_srf"); sd = g("dz_srf_sd")
        use = haskey(spec, "sites") ? [findfirst(==(String(s)), site) for s in spec["sites"]] : collect(eachindex(site))
        any(isnothing, use) && error("[cores] sites $(spec["sites"]): not all in $f ($(site))")
        pts = NamedTuple[]
        for s in use, a in ages
            k = [i for i in findall(x -> abs(x - a) <= hw, t) if isfinite(dz[s, i]) && isfinite(sd[s, i])]
            isempty(k) && error("[cores]: $(site[s]) has no record within $(hw) yr of $(a) years ago ($f)")
            push!(pts, (group = nm, kind = :core, site_name = @sprintf("%s %g ka", site[s], a / 1000), lat = lat[s], lon = lon[s], t = a,
                        value = mean(dz[s, k]), sd = mean(sd[s, k])))
        end
        pts
    end
end
# The projection of a FesmData grid file (its crs: polar_stereographic, latitude_of_projection_origin ±90, standard_parallel,
# straight_vertical_longitude_from_pole, semi_major_axis, inverse_flattening; Snyder 1987, eqs. 21-33 and 21-34 of the ellipsoid) and
# its axes xc, yc [km]: (xc, yc, xy), xy(lon, lat) the projected point [km].
function stereo_grid(file)
    NCDataset(file) do ds
        c = ds["crs"].attrib; get(c, "grid_mapping_name", "") == "polar_stereographic" || error("$file: crs not polar_stereographic")
        south = c["latitude_of_projection_origin"] < 0; φc = Float64(c["standard_parallel"]); λ0 = Float64(c["straight_vertical_longitude_from_pole"])
        a = Float64(c["semi_major_axis"]); e = sqrt(1 - (1 - 1 / Float64(c["inverse_flattening"]))^2)
        tf(φ) = tand(45 - φ / 2) / ((1 - e * sind(φ)) / (1 + e * sind(φ)))^(e / 2)
        sg = south ? -1.0 : 1.0; mc = cosd(sg * φc) / sqrt(1 - e^2 * sind(sg * φc)^2)
        u = get(ds["xc"].attrib, "units", "m") == "km" ? 1.0 : 1000.0
        function xy(lon, lat)
            ρ = a * mc * tf(sg * lat) / tf(sg * φc) / 1000; λ = sg * (lon - λ0)
            (sg * ρ * sind(λ), -sg * ρ * cosd(λ))
        end
        (xc = Float64.(ds["xc"][:]) ./ u, yc = Float64.(ds["yc"][:]) ./ u, xy = xy)
    end
end
# bilinear on a regular grid (ascending axes x, y) at the point (px, py)
function bilin_xy(F, x, y, px, py)
    i = clamp(searchsortedlast(x, px), 1, length(x) - 1); j = clamp(searchsortedlast(y, py), 1, length(y) - 1)
    wx = (px - x[i]) / (x[i+1] - x[i]); wy = (py - y[j]) / (y[j+1] - y[j])
    (1 - wy) * ((1 - wx) * F[i, j] + wx * F[i+1, j]) + wy * ((1 - wx) * F[i, j+1] + wx * F[i+1, j+1])
end
# The model's surface elevation change [m] at the points pts (lon, lat, t the age [yr]), relative to the present: e(t) - e(0 ka), e the
# surface above the sea level of the time, linear in time between the records (an age outside them: NaN); NaN for a run that has not
# finished (run_done). Of a kryos run with the ice output <gr.ice>.nc (the plastic or Yelmo instance on the grid of gr, stereo_grid):
# z_srf - z_sl, bilinear in the projected x, y; otherwise (a standalone VILMA2 run under a reconstruction's ice: vilma_out.nc) of the
# VILMA2 stream: z_bed + h_ice (z_bed relative to the sea level of the time), bilinear on the Gauss grid.
function core_at(dir, gr, pts; inst = "geo_glob")
    run_done(dir) || return fill(NaN, length(pts))
    fi = joinpath(dir, gr.ice * ".nc"); fv = joinpath(dir, "vilma_out.nc"); isfile(fv) || (fv = joinpath(dir, inst * "_vilma.nc"))
    ice = isfile(fi); f = ice ? fi : fv; isfile(f) || return fill(NaN, length(pts))
    NCDataset(f) do ds
        age = -10.0 .* round.(Float64.(ds["time"][:]) ./ 10.0); p = sortperm(age)      # records' ages (kryos: t = -age; VILMA2: -age + 1)
        if ice
            xy = [gr.xy(o.lon, o.lat) for o in pts]
            surf = k -> (F = Float64.(nomissing(ds["z_srf"][:, :, k], NaN)) .- Float64.(nomissing(ds["z_sl"][:, :, k], NaN));
                         [bilin_xy(F, gr.xc, gr.yc, x, y) for (x, y) in xy])
        else
            lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:]); jl = sortperm(lat); il = sortperm(mod.(lon, 360)); glon = mod.(lon[il], 360)
            w = [(xweight(glon, o.lon), yweight(lat[jl], o.lat)) for o in pts]
            surf = k -> (F = (Float64.(nomissing(ds["z_bed"][:, :, k], NaN)) .+ Float64.(nomissing(ds["h_ice"][:, :, k], NaN)))[il, jl];
                         [interp(F, a, b) for (a, b) in w])
        end
        series_at(age[p], k -> surf(p[k]), pts)
    end
end
# The surface elevation change of a published reconstruction (FesmData IceHistory file on a lon-lat grid: z_topo, the topography
# relative to the sea level of the time, the surface of grounded ice on the ice sheets; time in years ago) at pts, as core_at.
function core_product(file, pts)
    NCDataset(file) do ds
        lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:]); jl = sortperm(lat); il = sortperm(mod.(lon, 360)); glon = mod.(lon[il], 360)
        age = Float64.(ds["time"][:]); p = sortperm(age)
        w = [(xweight(glon, o.lon), yweight(lat[jl], o.lat)) for o in pts]
        series_at(age[p], k -> (F = Float64.(nomissing(ds["z_topo"][:, :, p[k]], NaN))[il, jl]; [interp(F, a, b) for (a, b) in w]), pts)
    end
end
# per point i the value of the series at its age pts[i].t minus that at 0 ka, linear between the records of the ascending ages A,
# srf(k) the values of record k at all points (read once per record needed)
function series_at(A, srf, pts)
    cache = Dict{Int,Vector{Float64}}(); rec(k) = haskey(cache, k) ? cache[k] : (cache[k] = srf(k))
    at(i, a) = (a < A[1] - 1 || a > A[end] + 1) ? NaN : (k = clamp(searchsortedlast(A, a), 1, length(A) - 1);
                w = clamp((a - A[k]) / (A[k+1] - A[k]), 0.0, 1.0); w == 0 ? rec(k)[i] : w == 1 ? rec(k + 1)[i] : (1 - w) * rec(k)[i] + w * rec(k + 1)[i])
    [at(i, o.t) - at(i, 0.0) for (i, o) in enumerate(pts)]
end

# --- Earth rotation ---------------------------------------------------------------------------------------------------------------
# The points of the rotation group ([rotation] of an ensemble config: plastic_ensemble.jl): per [[rotation.data]] entry the quantity
# (ROT_QUANTITIES), its value and 1 sigma, the data's sd and the model error sd_model combined (sqrt(sd^2 + sd_model^2)); each a point
# at (0, 0), t = 0 (the present). The true polar wander as a velocity, its two components along and across tpw_axis_deg ([rotation]:
# the observed direction [deg east]): a rate near zero has no direction, so the components, not rate and direction, are the data.
# Returns (group, kind = :rot, site_name = the quantity, lat, lon, t, value, sd).
const ROT_QUANTITIES = ("J2dot", "TPW along", "TPW across")   # [1e-11 / yr], [deg / Myr] along and across the axis (90 deg east of it)
function read_rotation(spec)
    nm = String(get(spec, "name", "rotation"))
    map(get(spec, "data", Any[])) do d
        q = String(d["quantity"])
        q in ROT_QUANTITIES || error("[[rotation.data]] quantity \"$q\": not one of ", join(ROT_QUANTITIES, ", "),
                                     " (TPW rate and direction became the components along and across tpw_axis_deg)")
        (group = nm, kind = :rot, site_name = q, lat = 0.0, lon = 0.0, t = 0.0, value = Float64(d["value"]),
         sd = hypot(Float64(d["sd"]), Float64(get(d, "sd_model", 0.0))))
    end
end
# the Gauss-Legendre weights of the latitudes of a Gauss grid (the nodes their sines; sum 2)
function gauss_weights(lat)
    x = sind.(lat); n = length(x)
    map(x) do xi
        p0, p1 = 1.0, xi
        for k in 2:n; p0, p1 = p1, ((2k - 1) * xi * p1 - (k - 1) * p0) / k; end
        dp = n * (xi * p1 - p0) / (xi^2 - 1)
        2 / ((1 - xi^2) * dp^2)
    end
end
const R_EARTH_M = 6.371e6
# J2 of the sea surface change N (VILMA2 dz_ss = N + Δφ [m], lon x lat on the Gauss grid of the latitudes lat, weights w): the
# degree-2 zonal Stokes coefficient of the potential perturbation, J2 = -√5 C̄20, C̄20 = (1 / 4π a) ∫ N P̄20 dΩ (P̄20 = √5 P2, fully
# normalized), so J2 = -(5 / 4π a) ∫ N P2 dΩ. Δφ (degree 0) and VILMA2's rotational potential (degree 2, order 1) do not enter.
function j2_of(N, lat, w)
    P2 = (3 .* sind.(lat) .^ 2 .- 1) ./ 2
    -5 / (4π * R_EARTH_M) * sum(w[j] * P2[j] * sum(N[:, j]) * 2π / size(N, 1) for j in eachindex(lat))
end
# The model's present-day rotation of a run, over the last interval of the records of its VILMA2 stream (vilma_out.nc of a standalone
# run, <inst>_vilma.nc of a kryos run): J2dot [1e-11 / yr] from dz_ss (j2_of), the true polar wander [deg / Myr] from the polar motion
# m = m1 + i m2 (rot_m_re, rot_m_im; m1 toward Greenwich, m2 toward 90 E), its components along the axis [deg east] and across it (90
# deg east of the axis). NaN for what the stream lacks (rot_m: a kryos run before 2026-11-18, or one without the rotation solver),
# and all for a run that has not finished. The last interval: 100 yr in the standalone runs, 500 yr (the record interval) in the kryos
# runs; the model's own ice change over it is in the rates.
function rotation_rates(dir; axis = 0.0, inst = "geo_glob")
    out = Dict(q => NaN for q in ROT_QUANTITIES)
    run_done(dir) || return out
    f = joinpath(dir, "vilma_out.nc"); isfile(f) || (f = joinpath(dir, inst * "_vilma.nc")); isfile(f) || return out
    NCDataset(f) do ds
        t = Float64.(ds["time"][:]); n = length(t); n >= 2 || return out
        dt = t[n] - t[n-1]
        if haskey(ds, "dz_ss")
            lat = Float64.(ds["lat"][:]); w = gauss_weights(lat)
            N(k) = Float64.(nomissing(ds["dz_ss"][:, :, k], NaN))
            out["J2dot"] = 1e11 * (j2_of(N(n), lat, w) - j2_of(N(n - 1), lat, w)) / dt
        end
        if haskey(ds, "rot_m_re") && haskey(ds, "rot_m_im")
            m1 = Float64.(ds["rot_m_re"][:]); m2 = Float64.(ds["rot_m_im"][:])
            v1 = rad2deg((m1[n] - m1[n-1]) / dt) * 1e6; v2 = rad2deg((m2[n] - m2[n-1]) / dt) * 1e6
            out["TPW along"] = v1 * cosd(axis) + v2 * sind(axis); out["TPW across"] = -v1 * sind(axis) + v2 * cosd(axis)
        end
        out
    end
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

# 1 sigma of an index point for the misfit d at the model slope dmdt: the elevation's (half the 2-sigma range, on the side of d) and the
# age's (dmdt times the 1-sigma age uncertainty)
function ip_sigma(o, d, dmdt)
    sz = (d > 0 ? o.rmax - o.rsl : o.rsl - o.rmin) / 2
    sqrt(sz^2 + (dmdt * o.st)^2)
end

# score of one sample: index point -> (d, z), limiting -> violated (as 1.0 / 0.0)
function score(o, A, s)
    any(isnan, s) && return (m = NaN, d = NaN, z = NaN, viol = NaN)
    if o.type == 0
        m, dmdt = tval(A, s, o.t)
        d = m - o.rsl
        return (m = m, d = d, z = d / ip_sigma(o, d, dmdt), viol = NaN)
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

# --- pm-diag: one-factor variants of B1 against PaleoMIST ------------------------------------------------------------------------
#   julia --project=tools/levante tools/levante/pm_compare.jl diag <fesmdata> <out> <ref run folder> [<target>=<nc> ...] <run folder> ...
# Per run, at MAP_AGES and over all records 2.5-80 ka: RSL relative to present minus the target (PM: FesmData dz_sl; SELEN: PaleoMIST's
# own SELEN sea level from its Stokes coefficients, pm_prep.jl selen: e.g. SELEN = the final ice, SELEN_h = its 500-yr version) by region
# (RMS, bias), the ocean
# mean over PaleoMIST's ocean of the time (z_bed(t) < 0, no ice), the sites, and the change against the reference run (B1_l170).
# The record state of a run: t_k + 1 yr for the step forcings, t_k for the linear one (the first of the two present).
const DIAG_AGES = [80000.0, 60000.0, 22500.0, 10000.0]
trun(r::Run, a) = a == 80000.0 ? -80000.0 : (any(x -> abs(x - (-a + 1.0)) < 0.01, r.t) ? -a + 1.0 : -a)
function diag(args)
    fesmdata, out, refdir = args[1], args[2], args[3]
    opts = [split(a, '=')[1] => split(a, '=')[2] for a in args[4:end] if occursin('=', a)]
    dirs = [refdir; [a for a in args[4:end] if !occursin('=', a)]]
    R = [r for r in read_run.(dirs) if r !== nothing]
    mkpath(out)
    pm = NCDataset(joinpath(fesmdata, "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc"))
    lon = Float64.(pm["lon"][:]); lat = Float64.(pm["lat"][:]); ages = Float64.(pm["time"][:])
    sl = Float64.(nomissing(pm["dz_sl"][:, :, :], NaN)); H = Float64.(nomissing(pm["H_ice"][:, :, :], 0.0))
    zb = Float64.(nomissing(pm["z_bed"][:, :, :], NaN))
    close(pm)
    targets = [("PM", sl)]
    for (tn, file) in opts
        S = fill(NaN, size(sl))
        NCDataset(file) do ds
            slon = Float64.(ds["lon"][:]); slat = Float64.(ds["lat"][:]); sa = Float64.(ds["time"][:]); v = ds["dz_sl"][:, :, :]
            ii = [findfirst(x -> abs(mod(x - l + 180, 360) - 180) < 0.01, slon) for l in lon]; jj = [findfirst(x -> abs(x - l) < 0.01, slat) for l in lat]
            for (k, a) in enumerate(ages); S[:, :, k] = Float64.(v[ii, jj, findfirst(==(a), sa)]); end
        end
        push!(targets, (String(tn), S))
    end
    w = [cosd(lat[j]) for i in eachindex(lon), j in eachindex(lat)]
    glac = (dropdims(maximum(H; dims = 3); dims = 3) .> 0) .& (H[:, :, end] .== 0)
    ocean = (zb[:, :, end] .< 0) .& (H[:, :, end] .== 0) .& .!glac
    regions = [glac, ocean, trues(size(w))]
    kage(a) = findfirst(==(a), ages)
    rel(r, a) = bilin(r.rsl[:, :, tindex(r, trun(r, a))] .- r.rsl[:, :, end], r.lon, r.lat, lon, lat)
    pocean(f, a) = (m = (zb[:, :, kage(a)] .< 0) .& (H[:, :, kage(a)] .== 0) .& .!isnan.(f); sum(w[m] .* f[m]) / sum(w[m]))
    site = [(nm, argmin(abs.(mod.(lon .- so .+ 180, 360) .- 180)), argmin(abs.(lat .- sa))) for (nm, so, sa) in SITES if nm != "Gulf of Bothnia" && nm != "Barents Sea"]
    F = Dict((r.name, a) => rel(r, a) for r in R for a in ages if a > 0 && any(x -> abs(x - trun(r, a)) < 0.01, r.t))
    ref = R[1].name
    f1(x) = @sprintf("%.1f", x)
    io = open(joinpath(out, "pm_diag.md"), "w")
    println(io, "RSL relative to present, run - target [m] (area weighted, 0.5 degree): glaciated / ocean / global RMS, bias; the ocean mean over ",
            "PaleoMIST's ocean of the time (target in brackets); sites; change against $(ref) (glaciated RMS, ocean bias of run - $(ref)). ",
            "'all' = mean over the records 2.5-80 ka.\n")
    for (tn, T) in targets
        println(io, "\n### against $(tn)\n")
        println(io, "| run | ka | glac RMS | glac bias | ocean RMS | ocean bias | global RMS | ocean mean (target) | ",
                join([s[1] for s in site], " | "), " | vs $(ref): glac RMS | ocean bias |")
        println(io, "|---|---|---|---|---|---|---|---|", repeat("---|", length(site)), "---|---|")
        tt = [(tn2, T2) for (tn2, T2) in targets if tn2 != tn]
        for (rn, A) in [[(r.name, a -> get(F, (r.name, a), nothing)) for r in R]; [(tn2 * " (target)", a -> T2[:, :, kage(a)]) for (tn2, T2) in tt]]
            for a in [MAP_AGES; 0.0]
                if a > 0
                    f = A(a); f === nothing && continue
                    d = f .- T[:, :, kage(a)]; s = [stats(d, w, m) for m in regions]
                    dr = haskey(F, (ref, a)) ? f .- F[(ref, a)] : fill(NaN, size(f))
                    sr = [stats(dr, w, m) for m in regions[1:2]]
                    println(io, "| ", a == MAP_AGES[1] ? rn : "", " | ", a / 1000, " | ", f1(s[1].rms), " | ", f1(s[1].bias), " | ", f1(s[2].rms), " | ",
                            f1(s[2].bias), " | ", f1(s[3].rms), " | ", f1(pocean(f, a)), " (", f1(pocean(T[:, :, kage(a)], a)), ") | ",
                            join([f1(f[i, j] - T[i, j, kage(a)]) for (_, i, j) in site], " | "), " | ", f1(sr[1].rms), " | ", f1(sr[2].bias), " |")
                else          # mean over all records 2.5-80 ka
                    as = [x for x in ages if x > 0 && A(x) !== nothing]
                    isempty(as) && continue
                    m = [let d = A(x) .- T[:, :, kage(x)]; [stats(d, w, regions[1]).rms, stats(d, w, regions[1]).bias, stats(d, w, regions[2]).rms,
                         stats(d, w, regions[2]).bias, stats(d, w, regions[3]).rms] end for x in as]
                    mm = sum(m) ./ length(m)
                    println(io, "| | all | ", join(f1.(mm[1:5]), " | "), " | | ", repeat(" | ", length(site)), " | |")
                end
            end
        end
    end
    close(io)
    # per record: the metrics against each target, the ocean mean and the site values (RSL relative to present)
    open(joinpath(out, "pm_diag.csv"), "w") do c
        println(c, "run,target,age,glac_rms,glac_bias,ocean_rms,ocean_bias,global_rms,ocean_mean,", join([replace(s[1], " " => "_") for s in site], ","))
        for (tn, T) in targets, r in R, a in ages
            haskey(F, (r.name, a)) || continue
            f = F[(r.name, a)]; d = f .- T[:, :, kage(a)]; s = [stats(d, w, m) for m in regions]
            println(c, join((r.name, tn, a, s[1].rms, s[1].bias, s[2].rms, s[2].bias, s[3].rms, pocean(f, a), (f[i, j] for (_, i, j) in site)...), ","))
        end
        for (tn, T) in targets, a in ages
            a > 0 || continue
            f = T[:, :, kage(a)]
            println(c, join(("target", tn, a, NaN, NaN, NaN, NaN, NaN, pocean(f, a), (f[i, j] for (_, i, j) in site)...), ","))
        end
    end
    # fields at DIAG_AGES for pm_maps.jl diag
    f = joinpath(out, "pm_diag_fields.nc"); isfile(f) && rm(f)
    NCDataset(f, "c") do nc
        defDim(nc, "lon", length(lon)); defDim(nc, "lat", length(lat)); defDim(nc, "rec", length(DIAG_AGES))
        defVar(nc, "lon", lon, ("lon",)); defVar(nc, "lat", lat, ("lat",)); defVar(nc, "age", DIAG_AGES, ("rec",))
        defVar(nc, "H_lgm", H[:, :, kage(22500.0)], ("lon", "lat"))
        for (tn, T) in targets
            defVar(nc, "rsl_" * tn, cat([T[:, :, kage(a)] for a in DIAG_AGES]...; dims = 3), ("lon", "lat", "rec"))
        end
        for r in R
            all(a -> haskey(F, (r.name, a)), DIAG_AGES) || continue
            defVar(nc, "rsl_" * r.name, cat([F[(r.name, a)] for a in DIAG_AGES]...; dims = 3), ("lon", "lat", "rec"))
        end
    end
    println("wrote ", out)
end

function main(args)
    args[1] == "diag" && return diag(args[2:end])
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
abspath(PROGRAM_FILE) == @__FILE__() && main(ARGS)
