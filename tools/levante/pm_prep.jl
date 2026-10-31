# PaleoMIST test of the standalone VILMA2: the inputs of vilma.x from the FesmData PaleoMIST product (minimal MIS 3 scenario).
#   julia --project=tools tools/levante/pm_prep.jl <fesmdata> <out folder> [grid = GLOBAL-0.25DEG]
# pm_ref.nc     : z_bed, H_ice of the 0 ka record (the present-day reference, i_eq = 1).
# pm_force_B.nc : h_ice(lon, lat, time), the records as steps (SELEN's Heaviside history: the load of record t_k applies from t_k
#                 to t_{k+1}). vilma.x ramps the ice linearly between two slices (solid_earth_update), so each record enters with
#                 a ramp of DT_STEP years: slices t_k (ice of t_{k-1}) and t_k + DT_STEP (ice of t_k). Starts at 80 ka (the spin-up
#                 load), ends at 0 ka + DT_STEP. The state at t_k + DT_STEP is the one to compare with record t_k.
# pm_force_A.nc : the 22.5 ka ice held, two slices (-22500, -22500 + DT_STEP): the spin-up under it is the A run.
# pm_force_S.nc : smoke, the records 25, 22.5, 20 ka as for B.
# Model time t = -age (years before 1950).
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

function main(args)
    args[1] == "visc1d" && return write_visc_uniform(args[2], args[3])
    fesmdata, out = args[1], args[2]; grid = length(args) >= 3 ? args[3] : "GLOBAL-0.25DEG"
    mkpath(out)
    src = joinpath(fesmdata, "Global", grid, "$(grid)_IceHistory-PaleoMIST.nc")
    ds = NCDataset(src)
    ages = Float64.(ds["time"][:]); lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:])
    @assert issorted(ages, rev = true) && ages[end] == 0.0 && all(diff(ages) .== -2500.0)
    Hraw = ds["H_ice"][:, :, :]; zraw = ds["z_bed"][:, :, end]
    close(ds)
    nmH = count(ismissing, Hraw); nmz = count(ismissing, zraw)
    println(src, ": ", length(lon), " x ", length(lat), ", ", length(ages), " records ", ages[1], " to ", ages[end],
            " ka; missing H_ice ", nmH, ", z_bed(0) ", nmz)
    nmz == 0 || error("z_bed of the 0 ka record has missing values")
    H = Float32.(nomissing(Hraw, 0.0f0)); H[H .< 0] .= 0
    ice = [H[:, :, k] for k in eachindex(ages)]

    r = joinpath(out, "pm_ref.nc"); isfile(r) && rm(r)
    NCDataset(r, "c") do nc
        defDim(nc, "lon", length(lon)); defDim(nc, "lat", length(lat))
        defVar(nc, "lon", Float64, ("lon",); attrib = ["units" => "degrees_east"])[:] = lon
        defVar(nc, "lat", Float64, ("lat",); attrib = ["units" => "degrees_north"])[:] = lat
        defVar(nc, "z_bed", Float64, ("lon", "lat"); attrib = ["units" => "m", "comment" => "PaleoMIST z_bed, 0 ka"])[:, :] = Float64.(zraw)
        defVar(nc, "H_ice", Float64, ("lon", "lat"); attrib = ["units" => "m", "comment" => "PaleoMIST H_ice, 0 ka"])[:, :] = Float64.(ice[end])
    end
    println("wrote ", r)

    write_force(joinpath(out, "pm_force_B.nc"), lon, lat, steps(ages, ice)...)
    kA = findfirst(==(22500.0), ages)
    write_force(joinpath(out, "pm_force_A.nc"), lon, lat, [-22500.0, -22500.0 + DT_STEP], [ice[kA], ice[kA]])
    kS = findall(a -> a in (25000.0, 22500.0, 20000.0), ages)
    write_force(joinpath(out, "pm_force_S.nc"), lon, lat, steps(ages[kS], ice[kS])...)
end
main(ARGS)
