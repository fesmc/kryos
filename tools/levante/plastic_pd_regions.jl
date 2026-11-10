# The present-day ice thickness of the plastic backend against the observations, per region (docs/dev/plastic-backend.md, "Calibration
# v2"): the topography hubs of FesmData (NH-32KM: TOPO-GEBCO2026, Greenland from BedMachine v6, elsewhere IceBoost v2 or GEBCO;
# ANT-16KM: TOPO-BedMachine-v4; PAT-PAL-8KM: TOPO-GEBCO2026, IceBoost v2 or GEBCO) and the regions of their REGIONS maps. Included
# by plastic_taub_fit_report.jl and plastic_esmda_report.jl; as a script, the table of a kryos run (its ice_<instance>.nc at 0 ka):
#   julia --project=tools/levante tools/levante/plastic_pd_regions.jl <run dir> [<fesmdata>]
# Compared over the cells of the model's ice (the 0 ka margin) where the observation is not floating (the fit cells of Stage 1):
# volume [m SLE] of both, mean and rms of H - H_obs [m]; and the observed ice outside the margin (not representable).
using NCDatasets, Printf, Statistics

const PD_RHO_I, PD_RHO_SW, PD_A_OCN = 910.0, 1028.0, 3.618e14
const PD_GRIDS = Dict(
    "NH-32KM" => (domain = "North", topo = "NH-32KM_TOPO-GEBCO2026.nc", instance = "ice_nor",
                  regions = [("Greenland", "region_2", [103]), ("Canadian Arctic", "region_3", [10101]),
                             ("Cordillera and Alaska", "region_3", [10103]), ("Hudson Bay and Strait", "region_3", [10102]),
                             ("Iceland", "region_2", [105]), ("Svalbard", "region_3", [10203]),
                             ("Russian Arctic (Barents-Kara)", "region_3", [10201]), ("Fennoscandia", "region_3", [10204]),
                             ("British Isles", "region_3", [10202]), ("Asia", "region_2", [104])]),
    "ANT-16KM" => (domain = "Antarctica", topo = "ANT-16KM_TOPO-BedMachine-v4.nc", instance = "ice_ant",
                   regions = [("West Antarctica", "region_3", [20101]), ("East Antarctica", "region_3", [20102]),
                              ("Antarctic Peninsula", "region_3", [20103])]),
    "PAT-PAL-8KM" => (domain = "PatagoniaPaleo", topo = "PAT-PAL-8KM_TOPO-GEBCO2026.nc", instance = "ice_pat",
                      regions = [("Patagonia", "region_2", [202])]))

# the observations and regions of a grid: H_obs, the observed mask (3 = floating), the source of each cell, the region masks
function pd_obs(fes, grid)
    g = PD_GRIDS[grid]; d = joinpath(fes, g.domain, grid)
    H, mo, src, srcnames, x, y = NCDataset(joinpath(d, g.topo)) do ds
        Float64.(nomissing(ds["H_ice"][:, :], 0.0)), Int.(ds["mask"][:, :]), Int.(ds["src_id"][:, :]),
        split(ds["src_id"].attrib["flag_meanings"]), Float64.(ds["xc"][:]), Float64.(ds["yc"][:])
    end
    H[H .< 0] .= 0.0
    R = NCDataset(joinpath(d, "$(grid)_REGIONS.nc")) do ds
        Dict(v => Int.(nomissing(ds[v][:, :], -1)) for v in ("region_2", "region_3"))
    end
    regs = [(n, in(codes).(R[v])) for (n, v, codes) in g.regions]
    rest = .!reduce(.|, last.(regs))
    push!(regs, ("other", rest))
    (H = H, mask = mo, src = src, srcnames = srcnames, regions = regs, dxdy = (x[2] - x[1]) * (y[2] - y[1]) * 1e6)
end

pd_sle(H, m, dxdy) = sum(H[m]; init = 0.0) * dxdy * PD_RHO_I / PD_RHO_SW / PD_A_OCN

# rows per region: the model fields (label => H) against the observations, over `margin` (the model's 0 ka ice) where not floating
function pd_table(io, ob, margin, models; title = "")
    fit = margin .& (ob.mask .!= 3)
    isempty(title) || println(io, title, "\n")
    print(io, "| region | fit cells | source of H_obs | H_obs [m SLE] | outside the margin [m SLE] |")
    for (l, _) in models; print(io, " ", l, ": volume (ratio) | bias, rms [m] |"); end
    println(io, "\n|---|---|---|---|---|", repeat("---|---|", length(models)))
    regs = [ob.regions; ("all", trues(size(margin)))]
    for (n, r) in regs
        k = fit .& r; out = .!margin .& r .& (ob.H .> 0)
        (count(k) == 0 && count(out) == 0) && continue
        srcs = join(["$(ob.srcnames[s + 1]) $(round(Int, 100 * count(k .& (ob.src .== s)) / max(count(k), 1))) %"
                     for s in sort(unique(ob.src[k]))], ", ")
        vo = pd_sle(ob.H, k, ob.dxdy)
        @printf(io, "| %s | %d | %s | %.3f | %.3f |", n, count(k), srcs, vo, pd_sle(ob.H, out, ob.dxdy))
        for (_, H) in models
            dH = (H .- ob.H)[k]; v = pd_sle(H, k, ob.dxdy)
            isempty(dH) ? print(io, " - | - |") : @printf(io, " %.3f (%.2f) | %+.0f, %.0f |", v, v / vo, mean(dH), sqrt(mean(dH .^ 2)))
        end
        println(io)
    end
    println(io)
end

# H_ice of a kryos run's instance at the record nearest 0 ka
function pd_run_H(dir, instance)
    f = joinpath(dir, "$(instance).nc"); isfile(f) || return nothing
    NCDataset(f) do ds
        t = Float64.(ds["time"][:]); Float64.(nomissing(ds["H_ice"][:, :, argmin(abs.(t))], 0.0))
    end
end

function pd_run(io, dir, fes; label = basename(dir))
    for grid in ("NH-32KM", "ANT-16KM")
        H = pd_run_H(dir, PD_GRIDS[grid].instance); H === nothing && continue
        pd_table(io, pd_obs(fes, grid), H .> 0, [label => H]; title = "$(grid), $(label), 0 ka")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    pd_run(stdout, ARGS[1], length(ARGS) >= 2 ? ARGS[2] : get(ENV, "FESMDATA", "/work/ba1442/fesmdata"))
end
