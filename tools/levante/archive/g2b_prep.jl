# G2b: the inputs of the standalone vilma.x from the global 0.5 degree PaleoMIST file.
#   julia --project=tools tools/levante/g2b_prep.jl <fesmdata> <out folder> <t0> <nstep> [dt = 10]
# g2b_force_<nstep>.nc : h_ice(lon, lat, time), slice k at the model time t0 + k dt holding the ice of t0 + (k - 1) dt (k >= 1;
#                        slice 0 and 1 hold the ice of t0): the load of the Earth in the run sequence is that of the previous
#                        step (lag 1, phase4-plan 4.1), so the standalone run gets the same history. Linear in time between
#                        the records (age = -t), as the datasrc of the globe_ice does (fesm_datasrc, years ago, time_ref 1950).
# g2b_ref.nc           : z_bed, H_ice of the record at 0 ka (the present-day reference, i_eq = 1)
using NCDatasets

function main(args)
    fesmdata, out = args[1], args[2]
    t0 = parse(Float64, args[3]); nstep = parse(Int, args[4]); dt = length(args) >= 5 ? parse(Float64, args[5]) : 10.0
    mkpath(out)
    src = joinpath(fesmdata, "Global", "GLOBAL-0.5DEG", "GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc")
    ds = NCDataset(src)
    ages = Float64.(ds["time"][:]); lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:])
    H = Float64.(nomissing(Array(ds["H_ice"]), 0.0)); zb = Float64.(nomissing(Array(ds["z_bed"]), 0.0))
    close(ds)
    @assert ages[end] == 0.0 && ages[1] > ages[end] && issorted(ages, rev = true)
    function hist_at(t)                         # as g2b_replay hist_at
        age = -t
        age >= ages[1] && return H[:, :, 1]
        age <= ages[end] && return H[:, :, end]
        k = findfirst(k -> age <= ages[k] && age > ages[k+1], 1:length(ages)-1)
        w = (ages[k] - age) / (ages[k] - ages[k+1])
        return (1 - w) .* H[:, :, k] .+ w .* H[:, :, k+1]
    end
    f = joinpath(out, "g2b_force_$(nstep).nc")
    isfile(f) && rm(f)
    nc = NCDataset(f, "c")
    defDim(nc, "lon", length(lon)); defDim(nc, "lat", length(lat)); defDim(nc, "time", nstep + 1)
    vlon = defVar(nc, "lon", Float64, ("lon",); attrib = ["units" => "degrees_east"]); vlon[:] = lon
    vlat = defVar(nc, "lat", Float64, ("lat",); attrib = ["units" => "degrees_north"]); vlat[:] = lat
    vt = defVar(nc, "time", Float64, ("time",); attrib = ["units" => "years"]); vt[:] = [t0 + k * dt for k in 0:nstep]
    vh = defVar(nc, "h_ice", Float64, ("lon", "lat", "time"); chunksizes = (length(lon), length(lat), 1), deflatelevel = 2, attrib = ["units" => "m", "comment" => "slice k = the ice of t0 + (k-1) dt, slice 0 = t0 (lag 1)"])
    for k in 0:nstep
        vh[:, :, k+1] = hist_at(t0 + max(k - 1, 0) * dt)
    end
    close(nc)
    r = joinpath(out, "g2b_ref.nc"); isfile(r) && rm(r)
    nc = NCDataset(r, "c")
    defDim(nc, "lon", length(lon)); defDim(nc, "lat", length(lat))
    vlon = defVar(nc, "lon", Float64, ("lon",); attrib = ["units" => "degrees_east"]); vlon[:] = lon
    vlat = defVar(nc, "lat", Float64, ("lat",); attrib = ["units" => "degrees_north"]); vlat[:] = lat
    defVar(nc, "z_bed", Float64, ("lon", "lat"); attrib = ["units" => "m"])[:, :] = zb[:, :, end]
    defVar(nc, "H_ice", Float64, ("lon", "lat"); attrib = ["units" => "m"])[:, :] = H[:, :, end]
    close(nc)
    println("wrote ", f, " (", nstep + 1, " slices, t = ", t0, " to ", t0 + nstep * dt, ") and ", r)
end
main(ARGS)
