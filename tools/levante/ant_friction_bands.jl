#!/usr/bin/env julia
# ant-friction (docs/dev/ant-friction.md): regional bands of ANT PD states, one table row per run (the bands of docs/dev/setup-diagnosis.md).
#   julia +1.13 --project=$G/julia-plot ant_friction_bands.jl <out.md> <run dir> ...      (sbatch ant_friction_bands.sh)
# Final Yelmo restart of each run. Masks fixed for all runs, from BedMachine v4 on ANT-32KM: grounded reference ice (H > 10 m, grounded by
# flotation), land / marine bed (z_bed >= 0 / < 0), distance to the edge of the grounded reference ice (interior > 300 km, 100-300 km, margin < 100 km);
# WAIS / EAIS from REGIONS region_3 (20101 / 20102). H bias: model H - BedMachine H on the grounded reference ice (model ice-free cells count).
# Velocity: log10(u_s model / Rignot 2011) median on reference cells with model and obs u_s > 1 m/yr. Sliding: u_b/u_s > 0.5 area fraction (model grounded ice).
using NCDatasets, Statistics, Printf

const FA = "/work/ba1442/fesmdata/Antarctica/ANT-32KM/ANT-32KM"
const RI, RSW = 917.0, 1028.0
rd(f, v) = (ds = NCDataset(f); x = ds[v]; a = ndims(x) == 3 ? x[:, :, size(x, 3)] : x[:, :]; a = Float64.(coalesce.(a, NaN)); close(ds); a)

function dist_edge(ice::BitMatrix, x, y)
    nx, ny = size(ice); bx = Float64[]; by = Float64[]
    for j in 1:ny, i in 1:nx
        ice[i, j] && continue
        any(1 <= i + di <= nx && 1 <= j + dj <= ny && ice[i+di, j+dj] for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1))) && (push!(bx, x[i]); push!(by, y[j]))
    end
    d = fill(NaN, nx, ny)
    for j in 1:ny, i in 1:nx; ice[i, j] && (d[i, j] = minimum(hypot.(bx .- x[i], by .- y[j]))); end
    d
end
med(x, m) = (k = m .& isfinite.(x); any(k) ? median(x[k]) : NaN)
afrac(c, w, m) = (s = sum(w[m]); s > 0 ? sum(w[m .& c]) / s : NaN)

function main(out, runs)
    bm = "$(FA)_TOPO-BedMachine-v4.nc"
    Hb, zb = rd(bm, "H_ice"), rd(bm, "z_bed"); zb[isnan.(zb)] .= -5000
    ds = NCDataset(bm); x = Float64.(ds["xc"][:]); y = Float64.(ds["yc"][:]); close(ds)
    maximum(abs, x) > 1e5 && (x ./= 1e3; y ./= 1e3)
    ref = (Hb .> 10) .& (Hb .* RI ./ RSW .> -zb)
    d = dist_edge(ref, x, y)
    ds = NCDataset("$(FA)_REGIONS.nc"); r3 = Int.(coalesce.(ds["region_3"][:, :], -1)); close(ds)
    wa, ea, land = r3 .== 20101, r3 .== 20102, zb .>= 0
    uo = rd("$(FA)_Velocity-Rignot2011.nc", "uxy_srf"); uo[.!(rd("$(FA)_Velocity-Rignot2011.nc", "f_valid_uxy_srf") .> 0.5)] .= NaN
    bands = ["all" => ref, "WAIS" => ref .& wa, "EAIS land" => ref .& ea .& land, "EAIS marine" => ref .& ea .& .!land,
             "int > 300 km" => ref .& (d .> 300), "100-300 km" => ref .& (d .>= 100) .& (d .<= 300), "margin < 100 km" => ref .& (d .< 100)]
    ubands = ["WAIS" => ref .& wa, "EAIS" => ref .& ea, "margin < 100 km" => ref .& (d .< 100)]
    io = IOBuffer()
    println(io, "H bias mean [m] on the grounded BedMachine ice (rms for all); log10(u_s model / Rignot 2011) median; sliding-dominated area (u_b/u_s > 0.5) of WAIS / EAIS.\n")
    println(io, "| run | H bias all (rms) | " * join(first.(bands[2:end]), " | ") * " | log10 u ratio " * join(first.(ubands), " / ") * " | sliding WAIS / EAIS |")
    println(io, "|---|" * repeat("---|", length(bands) + 2))
    for r in runs
        rs = filter(s -> startswith(s, "restart-"), readdir(r))
        f = joinpath(r, rs[argmax(parse.(Int, replace.(rs, "restart-" => "")))], "Antarctica", "ice_ant_restart.nc")
        H, us, ub, fg, a = rd(f, "H_ice"), rd(f, "uxy_s"), rd(f, "uxy_b"), rd(f, "f_grnd"), rd(f, "area")
        dH = H .- Hb; g = (H .> 10) .& (fg .>= 0.99)
        hb = [@sprintf("%+.0f", mean(dH[m])) for (_, m) in bands]
        hb[1] *= @sprintf(" (%.0f)", sqrt(mean(dH[ref] .^ 2)))
        lr = log10.(max.(us, 1e-3) ./ max.(uo, 1e-3)); ok = isfinite.(uo) .& (uo .> 1) .& (us .> 1)
        ur = join((@sprintf("%+.2f", med(lr, m .& ok)) for (_, m) in ubands), " / ")
        sl = ub ./ max.(us, 1e-6) .> 0.5
        println(io, "| $(basename(r)) | " * join(hb, " | ") * " | $ur | " * @sprintf("%.2f / %.2f", afrac(sl, a, g .& wa), afrac(sl, a, g .& ea)) * " |")
    end
    s = String(take!(io)); print(s); open(out, "w") do o; write(o, s); end
    return 0
end
exit(main(ARGS[1], ARGS[2:end]))
