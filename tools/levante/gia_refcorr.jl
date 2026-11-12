# The present-day reference correction r of a glacial-cycle GIA run (pm-ref, docs/dev/vilma2-paleomist.md): the reference bed of a
# run that starts from equilibrium must be today's bed plus the rebound still to come at its end. VILMA2 with i_eq = 3 takes
# z_bed_eq = z_bed_ref + r (r from rsl_restart_file), its bed at 0 ka is z_bed_eq - rsl(0 ka), so the 0 ka bed misfit is r - rsl(0 ka)
# and r is the fixed point r = rsl(0 ka) of the run with reference r. gia_refcorr.sh iterates it; this file does the rest.
#   julia --project=tools/levante tools/levante/gia_refcorr.jl prep  <ice history.nc> <linear|steps> <out folder>
#       force.nc (h_ice(lon, lat, time), model time = -age) and ref.nc (z_bed, H_ice of the history's 0 ka record: VILMA2's
#       present-day reference) from a FesmData IceHistory product. linear: with a slice DT_RATE before present (the ice of the
#       line), so VILMA2's dz_bed_rate at 0 ka is the rate of the last century, as in a kryos run (geo cadence 100 yr), not
#       the mean of the last record interval (2.5 kyr: too high where the rebound decays); the GNSS scores read it.
#   julia --project=tools/levante tools/levante/gia_refcorr.jl next  <vilma_out.nc> <r next.nc> [<r used.nc>]
#       r_{n+1} = rsl(0 ka) of the run, on its own Gauss grid; prints the 0 ka bed misfit r_n - rsl_n(0 ka) (max, location, RMS).
#   julia --project=tools/levante tools/levante/gia_refcorr.jl write <r.nc> <out folder> <X grid>:<history on that grid> <domain>:<hub grid>:<topo file> ...
#       The kryos reference files (bilinear from the Gauss grid; H_ice stays present-day):
#       <out>/<X grid>/<X grid>_REFCORR.nc            z_bed = the history's 0 ka bed + r, H_ice its 0 ka ice, r (globe_ice X)
#       <out>/<domain>/<hub>/<topo file name>-refcorr  the domain topography with z_bed + r and z_srf adjusted to it, r (the other
#                                                      variables, masks and fractions, are copied unchanged: kryos reads z_bed,
#                                                      z_srf, H_ice, z_bed_sd)
#       <out>/GLOBAL-0.5DEG_r.nc                       r on the 0.5 degree grid (the map)
#   Densities of z_srf: RHO_ICE, RHO_SW (environment; kryos phys_const Earth 910, 1028).
using NCDatasets, Printf, Dates

const DT_STEP = 1.0    # [yr] the ramp of a step (steps only)
const DT_RATE = 100.0  # [yr] the last forcing interval (linear only): the interval of dz_bed_rate at 0 ka

function prep(hist, interp, out)
    mkpath(out)
    ds = NCDataset(hist)
    ages = Float64.(ds["time"][:]); lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:])
    issorted(ages, rev = true) && ages[end] == 0.0 || error("$hist: time must be years ago, oldest first, ending at 0")
    H = Float32.(nomissing(ds["H_ice"][:, :, :], 0.0f0)); H[H .< 0] .= 0
    zb = Float64.(nomissing(ds["z_bed"][:, :, end], NaN))
    close(ds)
    any(isnan, zb) && error("$hist: the 0 ka z_bed has missing values")
    if interp == "linear"
        t = -ages; s = [H[:, :, k] for k in eachindex(ages)]
        if ages[end-1] > DT_RATE        # the slice at -DT_RATE, on the line between the last two records
            w = Float32(DT_RATE / ages[end-1])
            insert!(t, length(t), -DT_RATE); insert!(s, length(s), w .* H[:, :, end-1] .+ (1 - w) .* H[:, :, end])
        end
    elseif interp == "steps"            # record t_k from t_k to t_{k+1} (SELEN's Heaviside history; a 1 yr ramp)
        t = [-ages[1]]; s = [H[:, :, 1]]
        for k in 2:length(ages)
            push!(t, -ages[k], -ages[k] + DT_STEP); push!(s, H[:, :, k-1], H[:, :, k])
        end
    else
        error("interp: linear | steps")
    end
    f = joinpath(out, "force.nc"); isfile(f) && rm(f)
    NCDataset(f, "c") do nc
        defDim(nc, "lon", length(lon)); defDim(nc, "lat", length(lat)); defDim(nc, "time", length(t))
        defVar(nc, "lon", lon, ("lon",); attrib = ["units" => "degrees_east"])
        defVar(nc, "lat", lat, ("lat",); attrib = ["units" => "degrees_north"])
        defVar(nc, "time", t, ("time",); attrib = ["units" => "years", "comment" => "model time = -age"])
        v = defVar(nc, "h_ice", Float32, ("lon", "lat", "time"); chunksizes = (length(lon), length(lat), 1), deflatelevel = 2,
                   attrib = ["units" => "m", "comment" => "H_ice of $(hist), $(interp)"])
        for (k, x) in enumerate(s); v[:, :, k] = x; end
    end
    r = joinpath(out, "ref.nc"); isfile(r) && rm(r)
    NCDataset(r, "c") do nc
        defDim(nc, "lon", length(lon)); defDim(nc, "lat", length(lat))
        defVar(nc, "lon", lon, ("lon",); attrib = ["units" => "degrees_east"]); defVar(nc, "lat", lat, ("lat",); attrib = ["units" => "degrees_north"])
        defVar(nc, "z_bed", zb, ("lon", "lat"); attrib = ["units" => "m", "comment" => "0 ka z_bed of $(hist)"])
        defVar(nc, "H_ice", Float64.(H[:, :, end]), ("lon", "lat"); attrib = ["units" => "m", "comment" => "0 ka H_ice of $(hist)"])
    end
    println("wrote ", f, " (", length(t), " slices, ", t[1], " to ", t[end], ") and ", r)
end

function next(outnc, f, prev = "")
    lon, lat, r1 = NCDataset(outnc) do ds
        Float64.(ds["lon"][:]), Float64.(ds["lat"][:]), Float64.(ds["rsl"][:, :, end])
    end
    r0 = isempty(prev) ? zeros(size(r1)) : NCDataset(ds -> Float64.(ds["rsl"][:, :]), prev)
    d = r0 .- r1                                         # modelled 0 ka bed - today's bed
    w = [cosd(lat[j]) for i in eachindex(lon), j in eachindex(lat)]
    k = argmax(abs.(d))
    @printf("bed misfit at 0 ka (%s): max |d| = %.2f m at %.1f E %.1f N, RMS %.3f m; next correction: max %.1f, min %.1f m\n",
            outnc, abs(d[k]), lon[k[1]], lat[k[2]], sqrt(sum(w .* d .^ 2) / sum(w)), maximum(r1), minimum(r1))
    isfile(f) && rm(f)
    NCDataset(f, "c") do nc
        defDim(nc, "lon", length(lon)); defDim(nc, "lat", length(lat))
        defVar(nc, "lon", lon, ("lon",); attrib = ["units" => "degrees_east"]); defVar(nc, "lat", lat, ("lat",); attrib = ["units" => "degrees_north"])
        defVar(nc, "rsl", r1, ("lon", "lat"); attrib = ["units" => "m", "comment" => "rsl at 0 ka of $(outnc): the reference correction (i_eq = 3)"])
    end
end

# bilinear from the Gauss grid (lon periodic, ascending after sorting) to points
function sampler(file)
    lon, lat, r = NCDataset(file) do ds
        Float64.(ds["lon"][:]), Float64.(ds["lat"][:]), Float64.(ds["rsl"][:, :])
    end
    p = sortperm(mod.(lon, 360)); q = sortperm(lat); gl = mod.(lon[p], 360); gt = lat[q]; R = r[p, q]
    n = length(gl)
    function (x, y)
        xm = mod(x, 360); i = searchsortedlast(gl, xm)
        i0 = i == 0 ? n : i; i1 = (i == n || i == 0) ? 1 : i + 1
        x0 = i == 0 ? gl[n] - 360 : gl[i0]; x1 = i == n ? gl[1] + 360 : gl[i1]; wx = (xm - x0) / (x1 - x0)
        j = searchsortedlast(gt, y)
        (j0, j1, wy) = j == 0 ? (1, 1, 0.0) : j == length(gt) ? (j, j, 0.0) : (j, j + 1, (y - gt[j]) / (gt[j+1] - gt[j]))
        (1 - wy) * ((1 - wx) * R[i0, j0] + wx * R[i1, j0]) + wy * ((1 - wx) * R[i0, j1] + wx * R[i1, j1])
    end
end

calc_z_srf(H, zb, ri, rs) = max(zb + H, (1 - ri / rs) * H)    # htopo calc_z_srf, sea level 0

function stamp!(nc, rfile)
    nc.attrib["refcorr"] = "z_bed + r, r = rsl(0 ka) of the converged GIA run ($(rfile)); H_ice present-day; tools/levante/gia_refcorr.sh"
    nc.attrib["history"] = string(now()) * ": gia_refcorr.jl write\n" * get(nc.attrib, "history", "")
end

function write(rfile, out, xspec, doms...)
    ri = parse(Float64, get(ENV, "RHO_ICE", "910.0")); rs = parse(Float64, get(ENV, "RHO_SW", "1028.0"))
    S = sampler(rfile)
    # the 0.5 degree map
    lon = collect(-179.75:0.5:179.75); lat = collect(-89.75:0.5:89.75)
    f = joinpath(out, "GLOBAL-0.5DEG_r.nc"); isfile(f) && rm(f)
    NCDataset(f, "c") do nc
        defDim(nc, "lon", length(lon)); defDim(nc, "lat", length(lat))
        defVar(nc, "lon", lon, ("lon",)); defVar(nc, "lat", lat, ("lat",))
        defVar(nc, "r", [S(x, y) for x in lon, y in lat], ("lon", "lat"); attrib = ["units" => "m", "long_name" => "present-day reference correction"])
    end
    # globe_ice X: the history's 0 ka bed + r and its 0 ka ice on the X grid
    xg, xh = split(xspec, ':'; limit = 2)
    lo, la, zb, H = NCDataset(xh) do ds
        Float64.(ds["lon"][:]), Float64.(ds["lat"][:]), Float64.(nomissing(ds["z_bed"][:, :, end], NaN)), Float64.(nomissing(ds["H_ice"][:, :, end], 0.0))
    end
    any(isnan, zb) && error("$xh: missing 0 ka z_bed")
    r = [S(x, y) for x in lo, y in la]
    mkpath(joinpath(out, xg)); f = joinpath(out, xg, "$(xg)_REFCORR.nc"); isfile(f) && rm(f)
    NCDataset(f, "c") do nc
        defDim(nc, "lon", length(lo)); defDim(nc, "lat", length(la))
        defVar(nc, "lon", lo, ("lon",); attrib = ["units" => "degrees_east"]); defVar(nc, "lat", la, ("lat",); attrib = ["units" => "degrees_north"])
        defVar(nc, "z_bed", zb .+ r, ("lon", "lat"); attrib = ["units" => "m", "long_name" => "reference bed: 0 ka z_bed of the history + r"])
        defVar(nc, "H_ice", H, ("lon", "lat"); attrib = ["units" => "m", "long_name" => "0 ka H_ice of the history"])
        defVar(nc, "r", r, ("lon", "lat"); attrib = ["units" => "m", "long_name" => "present-day reference correction"])
        nc.attrib["grid_name"] = String(xg); nc.attrib["source"] = String(xh); stamp!(nc, rfile)
    end
    println("wrote ", f, ": r in ", extrema(r))
    # the domains: the hub topography + r
    for d in doms
        dom, hub, topo = split(d, ':'; limit = 3)
        gf = joinpath(dirname(topo), "$(hub)_grid.nc")
        lon2, lat2 = NCDataset(ds -> (Float64.(ds["lon2D"][:, :]), Float64.(ds["lat2D"][:, :])), gf)
        rd = S.(lon2, lat2)
        o = joinpath(out, dom, hub); mkpath(o)
        f = joinpath(o, replace(basename(topo), ".nc" => "-refcorr.nc")); cp(topo, f; force = true); chmod(f, 0o644)
        NCDataset(f, "a") do nc
            zb0 = Float64.(nomissing(nc["z_bed"][:, :], NaN)); H0 = Float64.(nomissing(nc["H_ice"][:, :], 0.0))
            zs0 = Float64.(nomissing(nc["z_srf"][:, :], NaN))
            size(zb0) == size(rd) || error("$topo and $gf differ in size")
            zb1 = zb0 .+ rd
            nc["z_bed"][:, :] = zb1
            nc["z_srf"][:, :] = zs0 .+ calc_z_srf.(H0, zb1, ri, rs) .- calc_z_srf.(H0, zb0, ri, rs)
            defVar(nc, "r", rd, dimnames(nc["z_bed"]); attrib = ["units" => "m", "long_name" => "present-day reference correction (added to z_bed)"])
            stamp!(nc, rfile)
        end
        println("wrote ", f, ": r in ", extrema(rd))
    end
end

function main(args)
    c = args[1]
    c == "prep"  && return prep(args[2:4]...)
    c == "next"  && return next(args[2:end]...)
    c == "write" && return write(args[2:end]...)
    error("gia_refcorr.jl prep | next | write")
end
main(ARGS)
