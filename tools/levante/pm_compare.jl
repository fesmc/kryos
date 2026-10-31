# PaleoMIST test of the standalone VILMA2: compare the runs of pm_run.sh with PaleoMIST (Gowan et al. 2021) and with each other.
#   julia --project=tools tools/levante/pm_compare.jl <fesmdata> <runs folder> <out folder> [B1_l256 B1_l170 B3_l170 A1_l256 A1_l170 A3_l170]
# The six run names are the folders of B1 (lmax 256), B1 (lmax 170), B3, A1 (256), A1 (170), A3; a missing run is skipped (smoke:
# pass the smoke folders). Everything on the 0.5 degree PaleoMIST grid: the VILMA2 fields (Gauss grid) bilinear in lon-lat.
# Definitions (PaleoMIST, checked): dz_sl = sea-level change from present (SELEN, positive = sea higher than the land), z_bed relative
# to the sea level of the time, so z_bed(t) - z_bed(0) = -dz_sl(t) (to 0.01 m). VILMA2: rsl = dz_ss - dz_bed (same sign as dz_sl),
# z_bed = z_bed_eq - rsl, so the bed change of Gowan's definition is -rsl of both and its metrics are those of rsl with the bias
# negated. dz_bed (u, the solid-Earth radial displacement) has no PaleoMIST counterpart; it is kept for the A - B comparison.
# Record t_k of B is the VILMA2 state at t_k + 1 yr (after the step of the record; pm_prep.jl), the 80 ka record the spun-up start.
# Relative to present: VILMA2's rsl is relative to the relaxed present-day reference (i_eq = 1), and the Earth is not relaxed at 0 ka
# (rsl(0) != 0), while PaleoMIST's dz_sl is S(t) - S(0) by construction; the B fields are therefore taken as f(t) - f(0 ka), as Gowan's.
# A - B (equilibrium - transient at 22.5 ka) compares the two states against the same relaxed reference, without that subtraction.
# Writes: pm_metrics.md (tables), pm_metrics_B1_PM.csv (all records), pm_sites.csv, pm_gmsl.csv, pm_fields.nc (for pm_maps.jl).
using NCDatasets, Printf, Statistics

const MAP_AGES = [22500.0, 15000.0, 10000.0]
const SITES = [("Hudson Bay", -85.0, 60.0), ("Gulf of Bothnia", 20.0, 63.0), ("Barents Sea", 40.0, 75.0),
               ("Ross embayment", -175.0, -80.0), ("Weddell embayment", -50.0, -77.0), ("Barbados", -59.6, 13.1)]

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

# bilinear from the (periodic in lon) Gauss grid to target lon-lat; latitude clamped beyond the outermost Gauss rows
function bilin(F, glon, glat, tlon, tlat)
    n = length(glon); out = similar(F, length(tlon), length(tlat))
    ix = map(tlon) do x
        xm = mod(x, 360); i = searchsortedlast(glon, xm)
        i0 = i == 0 ? n : i; i1 = i == n || i == 0 ? 1 : i + 1
        x0 = i == 0 ? glon[n] - 360 : glon[i0]; x1 = i == n ? glon[1] + 360 : glon[i1]
        (i0, i1, (xm - x0) / (x1 - x0))
    end
    jy = map(tlat) do y
        j = searchsortedlast(glat, y)
        j == 0 ? (1, 1, 0.0) : j == length(glat) ? (j, j, 0.0) : (j, j + 1, (y - glat[j]) / (glat[j+1] - glat[j]))
    end
    for (b, (j0, j1, wy)) in enumerate(jy), (a, (i0, i1, wx)) in enumerate(ix)
        out[a, b] = (1 - wy) * ((1 - wx) * F[i0, j0] + wx * F[i1, j0]) + wy * ((1 - wx) * F[i0, j1] + wx * F[i1, j1])
    end
    out
end

# --- metrics ----------------------------------------------------------------------------------------------------------------------
function stats(d, w, m)
    k = m .& .!isnan.(d)
    ww = w[k]; dd = d[k]
    (rms = sqrt(sum(ww .* dd .^ 2) / sum(ww)), max = maximum(abs.(dd)), bias = sum(ww .* dd) / sum(ww))
end

function main(args)
    fesmdata, runs, out = args[1], args[2], args[3]
    names = length(args) >= 9 ? args[4:9] : ["B1_l256", "B1_l170", "B3_l170", "A1_l256", "A1_l170", "A3_l170"]
    mkpath(out)
    R = Dict{String,Run}()
    for (key, nm) in zip(["B1", "B1_170", "B3", "A1", "A1_170", "A3"], names)
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
    for key in ("B1", "B1_170", "B3"), a in MAP_AGES
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
    hdr(io, "RSL: resolution (lmax 170 - 256) and 3-D (B3 - B1, lmax 170) effects")
    for (k1, k2) in (("B1_170", "B1"), ("B3", "B1_170")), a in MAP_AGES
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
    for key in ("B1", "B1_170", "B3")
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
            for key in ("B1", "B1_170", "B3")
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
        for key in ("B1", "B1_170", "B3")
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
        defVar(nc, "H_lgm", H[:, :, kage(22500.0)], ("lon", "lat"))
        defVar(nc, "rsl_PM", cat([sl[:, :, kage(a)] for a in MAP_AGES]...; dims = 3), ("lon", "lat", "rec"))
        for key in ("B1", "B1_170", "B3")
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
    println("wrote ", out)
end
main(ARGS)
