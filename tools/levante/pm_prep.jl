# PaleoMIST test of the standalone VILMA2: the inputs of vilma.x from the FesmData PaleoMIST product (minimal MIS 3 scenario).
#   julia --project=tools tools/levante/pm_prep.jl <fesmdata> <out folder> [grid = GLOBAL-0.25DEG]
# pm_ref.nc     : z_bed, H_ice of the 0 ka record (the present-day reference, i_eq = 1).
# pm_ref80.nc   : z_bed, H_ice of the 80 ka record (pm-diag: SELEN's reference, the initial state; run with TMAX=0).
# pm_force_B.nc : h_ice(lon, lat, time), the records as steps (SELEN's Heaviside history: the load of record t_k applies from t_k
#                 to t_{k+1}). vilma.x ramps the ice linearly between two slices (solid_earth_update), so each record enters with
#                 a ramp of DT_STEP years: slices t_k (ice of t_{k-1}) and t_k + DT_STEP (ice of t_k). Starts at 80 ka (the spin-up
#                 load), ends at 0 ka + DT_STEP. The state at t_k + DT_STEP is the one to compare with record t_k.
# pm_force_A.nc : the 22.5 ka ice held, two slices (-22500, -22500 + DT_STEP): the spin-up under it is the A run.
# pm_force_S.nc : smoke, the records 25, 22.5, 20 ka as for B.
# Variants of B for the diagnosis of the differences to PaleoMIST (pm-diag), one factor each; combined by the letters:
# pm_force_Bg.nc : as in Gowan's SELEN input (ICESHEET pixel/create_selen_input_2.sh of evangowan/selen, branch old): SELEN applies
#                  the column of time t_k over (t_k + 2.5 ka, t_k], so Gowan shifted the columns one step and dropped the 2.5 ka
#                  ice; the record ages keep their ice, except that the 0 ka ice enters at 2.5 ka.
# pm_force_Bl.nc : the records linear in time (the slices at the record ages, vilma.x ramps between them): the 500-yr version.
# pm_force_Bm.nc : SELEN without the displaced water (hypothesis): where the ice is grounded at the time (rho_i H >= rho_w times the
#                  water depth of the time, PaleoMIST z_bed) on a cell that is ocean today (z_bed(0) < 0, no ice), the present water
#                  column as ice equivalent, (rho_w/rho_i) (-z_bed(0)), is added to H. vilma.x removes that column again with the
#                  water it displaces (its ocean load s = C rsl - z0 (C - C0)) and counts the added ice as melt water, so the net
#                  load is rho_i H and the ocean volume -(rho_i/rho_w) H, as in an SLE that ignores the displaced water.
#                  rho_i/rho_w = RHO_I/RHO_W of the run (931/1000, or 910/1000 for the Bm910 variants). Bm80*: the same against
#                  the 80 ka state (ocean at 80 ka, its water depth), for the runs on pm_ref80.
# Model time t = -age (years before 1950).
#   julia --project=tools tools/levante/pm_prep.jl budget <fesmdata> <PaleoMIST bundle> <out folder>
# pm_budget.md, pm_budget.csv: ice volume change against 0 ka per region and record (FesmData 0.25 and 0.5 degree, total and above
#                  flotation, the displaced water), against Gowan's volume_sorted.txt (the 5 km ICESHEET grids) and his
#                  Supplementary Table 1.
#   julia --project=tools tools/levante/pm_prep.jl margins <fesmdata> <PaleoMIST bundle>
#                  registration of FesmData H_ice against Gowan's margin polygons (areas, overlap, centroid offset).
#   julia --project=tools tools/levante/pm_prep.jl taboo <prem200.dat>
#                  the layer means of TABOO (NV = 2, CDE = 0) for 120 km lithosphere, as SELEN's Earth of PaleoMIST (pm_run.sh t3).
#   julia --project=tools tools/levante/pm_prep.jl selen <shs.bin> <NN> <DELTA [kyr]> <out.nc> [<rsl_spreadsheet.dat>]
#                  PaleoMIST's SELEN sea level from its Stokes coefficients at the records on the 0.5 degree grid (checked against
#                  Gowan's calc_sl at points, pm_selen.sh).
#   julia --project=tools tools/levante/pm_prep.jl visc1d <profile.csv> <out.nc>
# pm_visc_bm.nc : the 1-D profile of pre_spinup_1d (pm_visc1d.x: per FE node, the lateral mean of the 3-D log10 eta) as a
#                 laterally uniform field in the form of the FesmData Bagge2021 product (log10_visc(lon, lat, depth,
#                 structure), one structure), on the node depths, so that vilma.x reads it back at its nodes unchanged.
using NCDatasets

const DT_STEP = 1.0     # [yr] the ramp of a record step; << the shortest Maxwell time (~15 yr at log10 eta = 19.5)

function write_force(f, lon, lat, t, slices)
    isfile(f) && rm(f)
    NCDataset(f, "c") do nc
        defDim(nc, "lon", length(lon)); defDim(nc, "lat", length(lat)); defDim(nc, "time", length(t))
        defVar(nc, "lon", Float64, ("lon",); attrib = ["units" => "degrees_east"])[:] = lon
        defVar(nc, "lat", Float64, ("lat",); attrib = ["units" => "degrees_north"])[:] = lat
        defVar(nc, "time", Float64, ("time",); attrib = ["units" => "years", "comment" => "model time = -age (years before 1950)"])[:] = t
        vh = defVar(nc, "h_ice", Float32, ("lon", "lat", "time"); chunksizes = (length(lon), length(lat), 1), deflatelevel = 2,
                    attrib = ["units" => "m", "comment" => "PaleoMIST H_ice; record t_k enters over [t_k, t_k + $(DT_STEP)] (step history)"])
        for (k, s) in enumerate(slices)
            vh[:, :, k] = s
        end
    end
    println("wrote ", f, ": ", length(t), " slices, t = ", t[1], " to ", t[end])
end

# the step history of the records ages (oldest first): t_0 with the first ice, then per later record t_k (previous ice), t_k + DT_STEP
function steps(ages, ice)
    t = [-ages[1]]; s = [ice[1]]
    for k in 2:length(ages)
        push!(t, -ages[k]); push!(s, ice[k-1])
        push!(t, -ages[k] + DT_STEP); push!(s, ice[k])
    end
    t, s
end

function write_visc_uniform(csv, f)
    l = filter(!isempty ∘ strip, readlines(csv)); h = split(l[1], ',')
    col(name) = [parse(Float64, split(r, ',')[findfirst(==(name), h)]) for r in l[2:end]]
    depth = col("depth"); v = col("log10_visc_mean")
    p = sortperm(depth); depth = depth[p]; v = v[p]
    depth[1] = max(depth[1], 0.0)                                      # the surface node (6371 km - r rounding)
    lon = collect(0.0:30.0:330.0); lat = collect(-90.0:30.0:90.0)
    isfile(f) && rm(f)
    NCDataset(f, "c") do nc
        defDim(nc, "lon", length(lon)); defDim(nc, "lat", length(lat)); defDim(nc, "depth", length(depth)); defDim(nc, "structure", 1)
        defVar(nc, "lon", lon, ("lon",); attrib = ["units" => "degrees_east"])
        defVar(nc, "lat", lat, ("lat",); attrib = ["units" => "degrees_north"])
        defVar(nc, "depth", depth, ("depth",); attrib = ["units" => "m", "positive" => "down", "comment" => "the FE node depths of vilma.x"])
        defVar(nc, "structure", Int32[1], ("structure",))
        defVar(nc, "log10_visc", [v[k] for i in lon, j in lat, k in eachindex(depth), s in 1:1], ("lon", "lat", "depth", "structure");
               attrib = ["units" => "log10(Pa s)", "long_name" => "lateral mean of log10 viscosity (pre_spinup_1d profile), laterally uniform"])
        nc.attrib["source"] = csv
    end
    println("wrote ", f, ": ", length(depth), " depths, log10 eta ", round(minimum(v), digits = 3), " to ", round(maximum(v), digits = 3))
end

# --- pm-diag: variants of the B forcing -----------------------------------------------------------------------------------------
const RHO_W = 1000.0
# the ice of the records with the marine-ice emulation (pm_force_Bm*): + (rho_w/rho_i) (-z_bed(ref)) where grounded at the time on a
# cell that is ocean in the reference state of the run (kref: the 0 ka record for pm_ref, the 80 ka record for pm_ref80)
function marine_ice(ice, zb, rho_i, kref = length(ice))
    ocean0 = (zb[kref] .< 0) .& (ice[kref] .== 0)
    map(eachindex(ice)) do k
        grounded = (ice[k] .> 0) .& (rho_i .* ice[k] .>= RHO_W .* max.(0.0f0, .-zb[k]))
        add = ocean0 .& grounded
        s = copy(ice[k]); s[add] .+= Float32(RHO_W / rho_i) .* .-zb[kref][add]; s
    end
end
# Gowan's SELEN input: the 0 ka ice from 2.5 ka on (the 2.5 ka ice dropped)
function gowan_ice(ages, ice)
    s = copy(ice); s[findfirst(==(2500.0), ages)] = ice[end]; s
end
function write_variants(out, lon, lat, ages, ice, zb)
    lin(ice) = (-ages, ice)
    for (nm, (t, s)) in ("Bg" => steps(ages, gowan_ice(ages, ice)), "Bl" => lin(ice),
                         "Bm" => steps(ages, marine_ice(ice, zb, 931.0)),
                         "Bm910" => steps(ages, marine_ice(ice, zb, 910.0)),
                         "Bgm910" => steps(ages, gowan_ice(ages, marine_ice(ice, zb, 910.0))),
                         "Bm80" => steps(ages, marine_ice(ice, zb, 931.0, 1)),
                         "Bgm80_910" => steps(ages, gowan_ice(ages, marine_ice(ice, zb, 910.0, 1))))
        write_force(joinpath(out, "pm_force_$(nm).nc"), lon, lat, t, s)
    end
end

# --- pm-diag: ice budgets -------------------------------------------------------------------------------------------------------
# Gowan et al. (2021), Supplementary Table 1 (minimal scenario): total volume change [1e6 km3] and the "calculated ESL" [m], by age
const ST1_AGES = collect(0.0:2500.0:80000.0)
const ST1_VOL = [0.0, 0.5, 1.2, 2.5, 10.2, 20.0, 29.6, 39.2, 42.2, 40.4, 35.2, 28.8, 22.0, 18.5, 15.8, 14.2, 10.5, 9.4, 14.2, 15.4,
                 13.7, 17.3, 21.3, 24.8, 29.1, 25.9, 18.6, 14.0, 10.0, 7.1, 4.5, 2.9, 3.3]
const ST1_ESL = [0.0, 0.0, -4.6, -8.4, -28.2, -54.1, -80.3, -107.1, -116.0, -111.4, -97.3, -79.7, -60.9, -51.2, -43.4, -38.5,
                 -28.4, -25.3, -37.7, -40.6, -36.2, -45.6, -56.2, -66.4, -78.9, -70.3, -50.7, -37.9, -26.8, -19.3, -12.6, -8.5, -9.5]
const REGIONS = ["North_America", "Eurasia", "Antarctica", "Patagonia"]
# Gowan's four ICESHEET domains on the lon-lat grid (North America with Greenland and Iceland)
function region_of(lo, la)
    lo = mod(lo + 180, 360) - 180
    la < -55 && return 3
    la < -35 && return (-80 < lo < -60 ? 4 : 0)
    la < 35 && return 0
    (lo < (la < 62 ? -14.0 : -10.0) && lo > -180 && !(lo > 170)) ? 1 : 2
end
cellarea(lat, dlat, dlon) = 6371.0^2 * deg2rad(dlon) * (sind(min(lat + dlat / 2, 90)) - sind(max(lat - dlat / 2, -90)))   # [km2]

function budget(fesmdata, bundle, out)
    mkpath(out)
    rows = Dict{Any,Any}()
    for grid in ("GLOBAL-0.25DEG", "GLOBAL-0.5DEG")
        ds = NCDataset(joinpath(fesmdata, "Global", grid, "$(grid)_IceHistory-PaleoMIST.nc"))
        ages = Float64.(ds["time"][:]); lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:])
        H = Float64.(nomissing(ds["H_ice"][:, :, :], 0.0)); H[H .< 0] .= 0
        zb = Float64.(nomissing(ds["z_bed"][:, :, :], 0.0))
        close(ds)
        A = [cellarea(lat[j], abs(lat[2] - lat[1]), abs(lon[2] - lon[1])) for i in eachindex(lon), j in eachindex(lat)]
        reg = [region_of(lon[i], lat[j]) for i in eachindex(lon), j in eachindex(lat)]
        k0 = findfirst(==(0.0), ages)
        ocean0 = (zb[:, :, k0] .< 0) .& (H[:, :, k0] .== 0)
        for (k, a) in enumerate(ages)
            h = H[:, :, k]; b = zb[:, :, k]
            af = max.(0.0, h .- (RHO_W / 910.0) .* max.(0.0, .-b))                         # above flotation, SELEN's 910 kg/m3
            af0 = max.(0.0, H[:, :, k0] .- (RHO_W / 910.0) .* max.(0.0, .-zb[:, :, k0]))
            g = (h .> 0) .& (910.0 .* h .>= RHO_W .* max.(0.0, .-b))
            for r in 1:4
                m = reg .== r
                dv = sum((A .* (h .- H[:, :, k0]))[m]) / 1e9                                # [1e6 km3]
                dvaf = sum((A .* (af .- af0))[m]) / 1e9
                vd0 = sum((A .* max.(0.0, .-zb[:, :, k0]))[m .& g .& ocean0]) / 1e9               # displaced water, present depth
                vdt = sum((A .* max.(0.0, .-b))[m .& g .& ocean0]) / 1e9                         # ... depth of the time
                rows[(grid, a, r)] = (dv = dv, dvaf = dvaf, vd0 = vd0, vdt = vdt)
            end
        end
        rows[grid] = ages
    end
    # Gowan's volume_sorted.txt: time, km3, 1e6 km3, ESL (361e6 km2, 910 kg/m3)
    native = Dict{Tuple{Float64,Int},Float64}()
    for (r, nm) in enumerate(REGIONS)
        for l in eachline(joinpath(bundle, "ice_reconstruction_files", nm, "volume_sorted.txt"))
            v = split(l); length(v) >= 3 || continue
            native[(parse(Float64, v[1]), r)] = parse(Float64, v[3])
        end
    end
    sle(v) = v * 1e6 * 0.91 / 361e6 * 1e3           # [m]: 1e6 km3 -> km of water over 361e6 km2 -> m
    ages = rows["GLOBAL-0.25DEG"]
    open(joinpath(out, "pm_budget.csv"), "w") do c
        println(c, "age,region,native,st1_total,dv_025,dv_05,dvaf_025,vd0_025,vdt_025")
        for a in ages, r in 1:4
            x = rows[("GLOBAL-0.25DEG", a, r)]; y = rows[("GLOBAL-0.5DEG", a, r)]
            println(c, join((a, REGIONS[r], get(native, (a, r), NaN), r == 1 ? ST1_VOL[findfirst(==(a), ST1_AGES)] : NaN,
                             x.dv, y.dv, x.dvaf, x.vd0, x.vdt), ","))
        end
    end
    f(x) = isnan(x) ? "-" : string(round(x, digits = 2))
    open(joinpath(out, "pm_budget.md"), "w") do io
        println(io, "Ice volume change against 0 ka [1e6 km3], by Gowan's domains: native = volume_sorted.txt (5 km ICESHEET grids); ",
                "FesmData 0.25 and 0.5 degree H_ice (total), above flotation (910/1000, the bed of the time), the water displaced by ",
                "grounded ice on cells that are ocean today (present depth / depth of the time). SLE = V 0.91 / 361e6 km2.\n")
        println(io, "| ka | domain | native | FesmData 0.25 | FesmData 0.5 | above flotation | displaced water, depth 0 ka / of the time |")
        println(io, "|---|---|---|---|---|---|---|")
        for a in (80000.0, 60000.0, 22500.0, 20000.0, 15000.0, 10000.0, 2500.0)
            tot = zeros(6)
            for r in 1:4
                x = rows[("GLOBAL-0.25DEG", a, r)]; y = rows[("GLOBAL-0.5DEG", a, r)]; n = get(native, (a, r), NaN)
                tot .+= [n, x.dv, y.dv, x.dvaf, x.vd0, x.vdt]
                println(io, "| ", r == 1 ? a / 1000 : "", " | ", REGIONS[r], " | ", f(n), " | ", f(x.dv), " | ", f(y.dv), " | ", f(x.dvaf),
                        " | ", f(x.vd0), " / ", f(x.vdt), " |")
            end
            i = findfirst(==(a), ST1_AGES)
            println(io, "| | total (Suppl. Table 1: ", ST1_VOL[i], ", ESL ", ST1_ESL[i], " m) | ", f(tot[1]), " | ", f(tot[2]), " (", f(sle(tot[2])),
                    " m) | ", f(tot[3]), " | ", f(tot[4]), " (", f(sle(tot[4])), " m) | ", f(tot[5]), " / ", f(tot[6]), " (", f(sle(tot[5]) / 0.91),
                    " m water) |")
        end
    end
    println("wrote ", out)
end

# --- pm-diag: TABOO's Earth (SELEN of PaleoMIST: calc_sl/data.inc NV = 2, CDE = 0, 120 km lithosphere) ---------------------------
# TABOO's Sbr PREM (selen_evan_mod/taboo/taboo.f90, the ANICE 2.1 version): volume means over [radinf, radsup] of the 200 s PREM table,
# step-wise with the value at the upper node of each interval; rigidity in 1e11 Pa, density in kg/m3
function taboo(premfile)
    l = [parse.(Float64, split(x, ',')) for x in eachline(premfile) if !isempty(strip(x))]
    r = [x[1] for x in l]; mu = [x[3] for x in l]; rho = [x[4] for x in l]
    function avg(top, dep)
        itop = findfirst(x -> x >= top, r)
        i = findfirst(x -> x >= dep, r); idep = r[i] == dep ? i : i - 1
        r[idep] == r[itop] && return (mu[itop] / 1e3, rho[itop] * 1e3)
        w(k) = (r[k+1]^3 - r[k]^3) / (r[itop]^3 - r[idep]^3)
        (sum(w(k) * mu[k+1] for k in idep:itop-1) / 1e3, sum(w(k) * rho[k+1] for k in idep:itop-1) * 1e3)
    end
    for (nm, top, dep) in (("lithosphere", 6371.0, 6251.0), ("upper mantle", 6251.0, 5701.0), ("lower mantle", 5701.0, 3480.0), ("core", 3480.0, 0.0))
        m, d = avg(top, dep)
        println(rpad(nm, 14), " r ", dep, "-", top, " km: mu = ", round(m * 1e11, sigdigits = 6), " Pa, rho = ", round(d, digits = 2), " kg/m3")
    end
end

# --- pm-diag: SELEN sea level of PaleoMIST from its Stokes coefficients ----------------------------------------------------------
# The sea-level SH coefficients S_j(k) of SELEN 2.9 (Fortran unformatted, one record COMPLEX*16 S(JMAX, 0:NN), j = l(l+1)/2 + m + 1,
# time index k = 0..NN from the oldest record to the present), synthesised as Gowan's calc_sl (rslc.f90, rslc-grid.f90):
# S(lon, lat) = sum_j (2 - delta_m0) Re(S_j Y_j), Y_lm the 4-pi normalised complex harmonics with the Condon-Shortley phase
# (HARMONICS.f90: PlmBar / sqrt(2 - delta_m0) times exp(i m lon)), and rsl(age k DELTA) = S(NN - k) - S(NN). Here by latitude (the
# Legendre functions once per latitude), on the 0.5 degree cell centres; the records at multiples of 2.5 kyr. Checked against
# Gowan's program at points (pm_selen.sh, rsl_spreadsheet.dat: max difference printed).
function plm_complex!(p, L, x)          # p[l+1, m+1] = (-1)^m sqrt((2l+1)(l-m)!/(l+m)!) P_lm(x)
    s = sqrt(max(0.0, 1 - x^2)); fill!(p, 0.0); p[1, 1] = 1.0
    for m in 0:L
        m > 0 && (p[m+1, m+1] = -sqrt((2m + 1) / (2m)) * s * p[m, m])
        m < L && (p[m+2, m+1] = sqrt(2m + 3) * x * p[m+1, m+1])
        for l in m+2:L
            a = sqrt((4l^2 - 1) / (l^2 - m^2)); b = sqrt(((l - 1)^2 - m^2) / (4(l - 1)^2 - 1))
            p[l+1, m+1] = a * (x * p[l, m+1] - b * p[l-1, m+1])
        end
    end
    p
end
function selen_read(file, L, nn)
    J = (L + 1) * (L + 2) ÷ 2
    open(file) do io
        n = read(io, Int32); @assert n == J * (nn + 1) * 16 "record length $n != JMAX (NN+1) 16 for L = $L, NN = $nn"
        S = Array{ComplexF64}(undef, J, nn + 1); read!(io, S); S
    end
end
function selen_synth(S, L, nn, lon, lat)
    F = zeros(length(lon), length(lat), nn + 1); p = zeros(L + 1, L + 1)
    E = [cis(m * deg2rad(x)) for x in lon, m in 0:L]
    jx(l, m) = l * (l + 1) ÷ 2 + m + 1
    for (j, la) in enumerate(lat)
        plm_complex!(p, L, sind(la))
        c = zeros(ComplexF64, L + 1, nn + 1)
        for m in 0:L, l in m:L
            @views c[m+1, :] .+= S[jx(l, m), :] .* p[l+1, m+1]
        end
        c[2:end, :] .*= 2
        F[:, j, :] = real.(E * c)
    end
    F
end
function selen_grid(file, nn, delta, out, check = "")
    L = 256; nn = parse(Int, nn); delta = parse(Float64, delta)
    lon = collect(-179.75:0.5:179.75); lat = collect(-89.75:0.5:89.75)
    Sx = selen_read(file, L, nn)
    F = selen_synth(Sx, L, nn, lon, lat)
    rsl(k) = F[:, :, nn - k + 1] .- F[:, :, nn + 1]                    # age k DELTA, relative to present
    keep = [k for k in 0:nn if isapprox(k * delta / 2.5, round(k * delta / 2.5); atol = 1e-6)]
    ages = [k * delta * 1000 for k in keep]
    if !isempty(check)                    # Gowan's calc_sl at points: age[kyr] lon lat rsl
        dmax = 0.0; n = 0; cache = Dict{Tuple{Float64,Float64},Vector{Float64}}()
        for l in eachline(check)
            v = parse.(Float64, split(l)); k = round(Int, v[1] / delta)
            c = get!(() -> selen_synth(Sx, L, nn, [v[2]], [v[3]])[1, 1, :], cache, (v[2], v[3]))
            dmax = max(dmax, abs((c[nn-k+1] - c[nn+1]) - v[4])); n += 1
        end
        println("check against ", check, ": ", n, " values, max |difference| = ", dmax, " m")
    end
    isfile(out) && rm(out)
    NCDataset(out, "c") do nc
        defDim(nc, "lon", length(lon)); defDim(nc, "lat", length(lat)); defDim(nc, "time", length(ages))
        defVar(nc, "lon", lon, ("lon",)); defVar(nc, "lat", lat, ("lat",)); defVar(nc, "time", ages, ("time",); attrib = ["units" => "years ago"])
        defVar(nc, "dz_sl", Float32.(cat([rsl(k) for k in keep]...; dims = 3)), ("lon", "lat", "time");
               attrib = ["units" => "m", "long_name" => "SELEN sea-level change from present, S(t) - S(0)"])
        nc.attrib["source"] = file
    end
    println("wrote ", out, ": ", length(ages), " records")
end

# --- pm-diag: registration of the FesmData ice against Gowan's margins (geographic GMT polygons of the PANGAEA bundle) ------------
#   pm_prep.jl margins <fesmdata> <PaleoMIST bundle> : per domain and record, the ice area of FesmData H_ice (0.25 degree, H > 0), the
#   area inside the margin polygons, their overlap, and the offset of the two centroids (east, north, km)
function read_gmt_polygons(f)
    polys = Vector{Vector{NTuple{2,Float64}}}(); cur = NTuple{2,Float64}[]
    for l in eachline(f)
        if startswith(l, ">")
            isempty(cur) || push!(polys, cur); cur = NTuple{2,Float64}[]
        elseif !startswith(l, "#") && !isempty(strip(l))
            v = parse.(Float64, split(l)[1:2]); push!(cur, (v[1], v[2]))
        end
    end
    isempty(cur) || push!(polys, cur)
    polys
end
function inside(polys, x, y)            # even-odd over all rings (holes included)
    c = false
    for p in polys, i in eachindex(p)
        (x1, y1) = p[i]; (x2, y2) = p[i == 1 ? end : i - 1]
        ((y1 > y) != (y2 > y)) && (x < (x2 - x1) * (y - y1) / (y2 - y1) + x1) && (c = !c)
    end
    c
end
function margins(fesmdata, bundle)
    ds = NCDataset(joinpath(fesmdata, "Global", "GLOBAL-0.25DEG", "GLOBAL-0.25DEG_IceHistory-PaleoMIST.nc"))
    ages = Float64.(ds["time"][:]); lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:])
    lonw = mod.(lon .+ 180, 360) .- 180
    println("| domain | ka | ice area FesmData [1e6 km2] | inside margins | overlap / union | centroid offset east, north [km] |")
    println("|---|---|---|---|---|---|")
    # North America also within the Laurentide core box (140-60 W, 38-68 N: without Greenland and the domain edges)
    core(lo, la) = -140 <= lo <= -60 && 38 <= la <= 68
    for (r, nm, box) in [[(r, nm, nothing) for (r, nm) in enumerate(REGIONS)]; (1, "North_America", core)], a in (60000.0, 20000.0, 10000.0)
        f = joinpath(bundle, "margins", nm, "$(Int(a)).gmt"); isfile(f) || continue
        P = read_gmt_polygons(f)
        H = Float64.(nomissing(ds["H_ice"][:, :, findfirst(==(a), ages)], 0.0))
        xs = [q[1] for p in P for q in p]; ys = [q[2] for p in P for q in p]
        sumi = summ = sumb = sumu = 0.0; ci = zeros(3); cm = zeros(3)
        for j in eachindex(lat), i in eachindex(lon)
            region_of(lonw[i], lat[j]) == r || continue
            box === nothing || box(lonw[i], lat[j]) || continue
            ice = H[i, j] > 0
            mar = (minimum(ys) <= lat[j] <= maximum(ys)) && (minimum(xs) <= lonw[i] <= maximum(xs)) && inside(P, lonw[i], lat[j])
            (ice || mar) || continue
            A = cellarea(lat[j], 0.25, 0.25); v = [cosd(lat[j]) * cosd(lonw[i]), cosd(lat[j]) * sind(lonw[i]), sind(lat[j])]
            ice && (sumi += A; ci .+= A .* v); mar && (summ += A; cm .+= A .* v)
            ice && mar && (sumb += A); sumu += A
        end
        g(c) = (c ./ sqrt(sum(c .^ 2)))
        gi = g(ci); gm = g(cm); lo = atand(gm[2], gm[1]); la = asind(gm[3])
        e = [-sind(lo), cosd(lo), 0.0]; n = [-sind(la) * cosd(lo), -sind(la) * sind(lo), cosd(la)]
        d = 6371.0 .* (gi .- gm)
        println("| ", nm, box === nothing ? "" : " (core box)", " | ", a / 1000, " | ", round(sumi / 1e6, digits = 2), " | ", round(summ / 1e6, digits = 2), " | ", round(sumb / sumu, digits = 3),
                " | ", round(sum(d .* e), digits = 0), ", ", round(sum(d .* n), digits = 0), " |")
    end
    close(ds)
end

function main(args)
    args[1] == "visc1d" && return write_visc_uniform(args[2], args[3])
    args[1] == "budget" && return budget(args[2], args[3], args[4])
    args[1] == "taboo" && return taboo(args[2])
    args[1] == "margins" && return margins(args[2], args[3])
    args[1] == "selen" && return selen_grid(args[2:end]...)
    fesmdata, out = args[1], args[2]; grid = length(args) >= 3 ? args[3] : "GLOBAL-0.25DEG"
    mkpath(out)
    src = joinpath(fesmdata, "Global", grid, "$(grid)_IceHistory-PaleoMIST.nc")
    ds = NCDataset(src)
    ages = Float64.(ds["time"][:]); lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:])
    @assert issorted(ages, rev = true) && ages[end] == 0.0 && all(diff(ages) .== -2500.0)
    Hraw = ds["H_ice"][:, :, :]; zraw = ds["z_bed"][:, :, end]
    zb = [Float32.(nomissing(ds["z_bed"][:, :, k], 0.0f0)) for k in eachindex(ages)]
    close(ds)
    nmH = count(ismissing, Hraw); nmz = count(ismissing, zraw)
    println(src, ": ", length(lon), " x ", length(lat), ", ", length(ages), " records ", ages[1], " to ", ages[end],
            " ka; missing H_ice ", nmH, ", z_bed(0) ", nmz)
    nmz == 0 || error("z_bed of the 0 ka record has missing values")
    H = Float32.(nomissing(Hraw, 0.0f0)); H[H .< 0] .= 0
    ice = [H[:, :, k] for k in eachindex(ages)]

    # pm_ref.nc: the present-day reference; pm_ref80.nc (pm-diag): SELEN's initial state, PaleoMIST's 80 ka bed and ice (the Earth in
    # equilibrium with them, run without spin-up), so that rsl and the ocean load are relative to the 80 ka state as SELEN's S
    for (nm, k, lbl) in (("pm_ref.nc", length(ages), "0 ka"), ("pm_ref80.nc", 1, "80 ka"))
        r = joinpath(out, nm); isfile(r) && rm(r)
        NCDataset(r, "c") do nc
            defDim(nc, "lon", length(lon)); defDim(nc, "lat", length(lat))
            defVar(nc, "lon", Float64, ("lon",); attrib = ["units" => "degrees_east"])[:] = lon
            defVar(nc, "lat", Float64, ("lat",); attrib = ["units" => "degrees_north"])[:] = lat
            defVar(nc, "z_bed", Float64, ("lon", "lat"); attrib = ["units" => "m", "comment" => "PaleoMIST z_bed, $lbl"])[:, :] = Float64.(zb[k])
            defVar(nc, "H_ice", Float64, ("lon", "lat"); attrib = ["units" => "m", "comment" => "PaleoMIST H_ice, $lbl"])[:, :] = Float64.(ice[k])
        end
        println("wrote ", r)
    end

    write_force(joinpath(out, "pm_force_B.nc"), lon, lat, steps(ages, ice)...)
    kA = findfirst(==(22500.0), ages)
    write_force(joinpath(out, "pm_force_A.nc"), lon, lat, [-22500.0, -22500.0 + DT_STEP], [ice[kA], ice[kA]])
    kS = findall(a -> a in (25000.0, 22500.0, 20000.0), ages)
    write_force(joinpath(out, "pm_force_S.nc"), lon, lat, steps(ages[kS], ice[kS])...)
    write_variants(out, lon, lat, ages, ice, zb)
end
main(ARGS)
