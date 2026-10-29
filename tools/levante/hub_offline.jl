# G1b (4) on a real run: the last record of globe_ice.nc against the merge computed offline (KryosTools.GlobeIce, independent of the
# Fortran code), from the domains' ice of the run, the background history and the coupler's maps.
#   julia --project=tools tools/levante/hub_offline.jl <run folder> <hub_areas.x> <bg history file> [--dt 10] [--tol-abs X]
# The run: gate sci-hub (domains North NH-32KM, Antarctica ANT-32KM; globe_ice on GLOBAL-0.5DEG, ext = file held at its 0 ka
# slice, so the external ice is the background X). The ice of the last gather is the ice exports of the step before the output time:
# H_ice of `<ice>_yelmo_sm.nc` at the output time minus dt (small 2D set, every dtt in the gate); both candidate times (output time
# and minus dt) are tried and reported. The areas are the hub's model areas (hub_areas.x on the run's map folder).
using NCDatasets, Printf
using KryosTools.GlobeIce

function main(args)
    rundir, areasx, bgfile = args[1:3]
    dt = 10.0; tolabs = 0.0
    i = findfirst(==("--dt"), args); i === nothing || (dt = parse(Float64, args[i+1]))
    i = findfirst(==("--tol-abs"), args); i === nothing || (tolabs = parse(Float64, args[i+1]))
    maps = realpath(joinpath(rundir, "maps"))
    doms = [("North", "NH-32KM", "ice_nor"), ("Antarctica", "ANT-32KM", "ice_ant")]
    bg = "GLOBAL-0.5DEG"
    tmp = mktempdir()
    areas = joinpath(tmp, "areas.nc")
    run(`$areasx $maps $areas $bg $([d[2] for d in doms])`)
    tlast = NCDataset(joinpath(rundir, "globe_ice.nc")) do g
        haskey(g, "time") ? Float64(Array(g["time"])[end]) : NaN
    end
    results = []
    for dtoff in (0.0, dt)
        inp = joinpath(tmp, "inputs_$(Int(dtoff)).nc")
        NCDataset(inp, "c") do ds
            NCDatasets.defDim(ds, "lon", 720); NCDatasets.defDim(ds, "lat", 360)
            ab = NCDataset(areas) do a; Float64.(a["area_$bg"][:, :]); end
            v = defVar(ds, "area_bg", Float64, ("lon", "lat")); v[:, :] = ab
            hb = NCDataset(bgfile) do b; Float64.(b["H_ice"][:, :, end]); end     # 0 ka: the last record
            v = defVar(ds, "H_ext", Float64, ("lon", "lat")); v[:, :] = hb
            ds.attrib["domains"] = join(first.(doms), ",")
            for (d, grid, ice) in doms
                f = joinpath(rundir, ice * "_yelmo_sm.nc")
                H = NCDataset(f) do s
                    t = Float64.(Array(s["time"])); k = argmin(abs.(t .- (tlast - dtoff)))
                    Float64.(s["H_ice"][:, :, k])
                end
                nx, ny = size(H)
                defDim(ds, "x_$d", nx); defDim(ds, "y_$d", ny)
                v = defVar(ds, "H_$d", Float64, ("x_$d", "y_$d")); v[:, :] = H
                ad = NCDataset(areas) do a; Float64.(a["area_$grid"][:, :]); end
                v = defVar(ds, "area_$d", Float64, ("x_$d", "y_$d")); v[:, :] = ad
                ds.attrib["map_$d"] = joinpath(maps, "scrip-con_$(grid)_$(bg).nc")
            end
        end
        r = check_globe_ice(joinpath(rundir, "globe_ice.nc"), inp; tol_abs = tolabs)
        push!(results, (dtoff, r))
        @printf("ice of t_out - %g yr: max |H_ice| error %.3e, max |cov| error %.3e, domains %d -> %s\n", dtoff, r.max_err_H, r.max_err_cov, r.n_domains, r.ok ? "equal" : "DIFFERENT")
    end
    return any(x -> x[2].ok, results) ? 0 : 1
end
exit(main(ARGS))
