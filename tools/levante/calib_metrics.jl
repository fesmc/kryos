#!/usr/bin/env julia
# calib-cf metrics (docs/dev/calib-cf.md): the final state of a calib run against the references, on Levante (NCDatasets only).
#   julia calib_metrics.jl <run dir> <Antarctica|North> <pd|lgm> <out prefix>
# Writes <prefix>.md (tables), <prefix>_ts.tsv (V_sle, V_ice_g, A_ice_g of the 1-D file every 500 yr), North: <prefix>_reg_ts.tsv (V_af of GRL/NA/EUR
# of every 2-D record), and <prefix>_fields.nc (model and reference
# grounded thickness, the ice masks of the model and (pd) the reference, for the maps of tools/calib_maps.jl). Grounded ice: H_ice where f_grnd >= 0.5 (model, BedMachine, GEBCO); the history files hold
# grounded ice only (PaleoMIST) or ice thickness (GLAC-1D, ICE-6G_C: every H > 0 is taken as grounded; the bed there is z_srf - H_ice).
# V_af: ice above flotation as sea level equivalent, rho_i 917, rho_sw 1028, ocean area 3.625e14 m2, flotation against the sea level of the state
# (model: z_sl of the output; references: 0 for the present day, -120 m for the LGM). Cell areas: Yelmo's `area` of the run (the same for the references).
using NCDatasets, Printf, Statistics
const RI, RSW, AOC = 917.0, 1028.0, 3.625e14
rd(f, v; kw...) = (ds = NCDataset(f); x = nomissing(ds[v][:, :, kw...], NaN); close(ds); x)

function lastrec(f, vars)
    ds = NCDataset(f); t = ds["time"][:]; i = length(t)
    out = Dict(v => Float64.(nomissing(ds[v][:, :, i], NaN)) for v in vars); close(ds); (t[i], out)
end
function nearest(f, v, ta)       # ta [yr ago]; linear interpolation in time between the two bracketing records; returns the field
    ds = NCDataset(f); t = Float64.(ds["time"][:])
    rec(j) = Float64.(nomissing(ds[v][:, :, j], 0.0))
    p = sortperm(t); ts = t[p]
    ta = clamp(ta, ts[1], ts[end]); k = searchsortedlast(ts, ta)
    if ts[k] == ta || k == length(ts); x = rec(p[k])
    else w = (ta - ts[k]) / (ts[k+1] - ts[k]); x = (1 - w) .* rec(p[k]) .+ w .* rec(p[k+1]); end
    close(ds); x
end
clean(x) = replace(x, NaN => 0.0)

function margin(g)                 # grounded cells with a non-grounded 4-neighbour (outside the array counts as non-grounded)
    nx, ny = size(g); m = falses(nx, ny)
    for j in 1:ny, i in 1:nx
        g[i, j] || continue
        for (di, dj) in ((1,0),(-1,0),(0,1),(0,-1))
            ii, jj = i + di, j + dj
            if ii < 1 || ii > nx || jj < 1 || jj > ny || !g[ii, jj]; m[i, j] = true; break; end
        end
    end
    m
end
function gl_dist(ga, gb, xc, yc)   # km, from the margin cells of ga to the nearest margin cell of gb
    ma, mb = findall(margin(ga)), findall(margin(gb))
    (isempty(ma) || isempty(mb)) && return Float64[]
    xb = [xc[c[1]] for c in mb]; yb = [yc[c[2]] for c in mb]
    [minimum(hypot.(xb .- xc[c[1]], yb .- yc[c[2]])) / 1e3 for c in ma]
end

function main(args)
    run, dom, case, pre = args; sfx = dom == "Antarctica" ? "ant" : "nor"; grid = dom == "Antarctica" ? "ANT-32KM" : "NH-32KM"
    fd = "/work/ba1442/fesmdata/$dom/$grid/$grid"
    # model
    fice = "$run/ice_$sfx.nc"; fy = "$run/ice_$(sfx)_yelmo.nc"; fts = "$run/ice_$(sfx)_yelmo_ts.nc"
    t_end, M = lastrec(fice, ["H_ice", "f_grnd", "z_bed", "z_sl"])
    area = Float64.(nomissing(NCDataset(fy)["area"][:, :, 1], 0.0)); area .= area .* (maximum(area) < 1e6 ? 1e6 : 1.0)   # km2 -> m2 if so
    ds = NCDataset(fy); xc = Float64.(ds["xc"][:]) .* (maximum(abs.(ds["xc"][:])) < 1e4 ? 1e3 : 1.0); yc = Float64.(ds["yc"][:]) .* (maximum(abs.(ds["yc"][:])) < 1e4 ? 1e3 : 1.0)
    lat = Float64.(nomissing(ds["lat2D"][:, :], 0.0)); lon = mod.(Float64.(nomissing(ds["lon2D"][:, :], 0.0)) .+ 180, 360) .- 180; close(ds)
    Hm = clean(M["H_ice"]); fg = clean(M["f_grnd"]); zb = clean(M["z_bed"]); zsl = clean(M["z_sl"])
    gm = (fg .>= 0.5) .& (Hm .> 0); Hg_m = Hm .* gm
    sle(H, g, zb, zsl) = sum(@. area * g * max(H - max(zsl - zb, 0.0) * RSW / RI, 0.0)) * RI / RSW / AOC
    stats(H, g, zb, zsl) = (V = sum(area .* H .* g) / 1e15, Vaf = sle(H, g, zb, zsl), A = sum(area .* g) / 1e12, Hmax = maximum(H .* g))
    # references: name => (H grounded, grounded mask, z_bed, z_sl)
    refs = Vector{Tuple{String,Matrix{Float64},BitMatrix,Matrix{Float64},Float64}}()
    topo = "$fd" * (dom == "Antarctica" ? "_TOPO-BedMachine-v4.nc" : "_TOPO-GEBCO2026.nc")
    if case == "pd"
        H = clean(rd(topo, "H_ice")); g = (clean(rd(topo, "f_grnd")) .>= 0.5) .& (H .> 0)
        push!(refs, ((dom == "Antarctica" ? "BedMachine v4" : "GEBCO 2026"), H .* g, g, clean(rd(topo, "z_bed")), 0.0))
    else
        zl = -120.0
        for (nm, ta) in (("PaleoMIST 20 ka", 20e3), ("PaleoMIST 22.5 ka", 22.5e3))
            f = "$(fd)_IceHistory-PaleoMIST.nc"; H = nearest(f, "H_ice", ta); zbr = nearest(f, "z_bed", ta)
            push!(refs, (nm, H, H .> 0, zbr, zl))
        end
        for (nm, f, ta) in (("GLAC-1D 21 ka", "$(fd)_IceHistory-GLAC1D.nc", 21e3), ("ICE-6G_C 21 ka", "$(fd)_IceHistory-ICE6G_C.nc", 21e3))
            H = nearest(f, "H_ice", ta); zs = nearest(f, "z_srf", ta)
            push!(refs, (nm, H, H .> 0, zs .- H, zl))
        end
    end
    io = IOBuffer()
    s0 = stats(Hg_m, gm, zb, zsl)
    ts = NCDataset(fts); tt = Float64.(ts["time"][:]); vsle = Float64.(nomissing(ts["V_sle"][:], NaN)); vg = Float64.(nomissing(ts["V_ice_g"][:], NaN)); ag = Float64.(nomissing(ts["A_ice_g"][:], NaN)); close(ts)
    println(io, "## $dom $case: $run (final t = $t_end yr)\n")
    println(io, "Model (grounded ice, f_grnd >= 0.5): V = $(@sprintf("%.3f", s0.V)) 1e6 km3, V_af = $(@sprintf("%.2f", s0.Vaf)) m SLE (Yelmo V_sle of the last 1-D record: $(@sprintf("%.2f", vsle[end])) m), A = $(@sprintf("%.3f", s0.A)) 1e6 km2, H_max = $(@sprintf("%.0f", s0.Hmax)) m.\n")
    nT = length(tt); sel(dt) = (j = findlast(<=(tt[end] - dt), tt); j === nothing ? 1 : j)
    rate(dt) = (vsle[end] - vsle[sel(dt)]) / ((tt[end] - tt[sel(dt)]) / 1e3)
    println(io, @sprintf("Convergence (Yelmo V_sle): %.2f m/kyr over the last 1 kyr, %.2f over 2 kyr, %.2f over 5 kyr; V_sle at the start %.2f m, at 5-kyr marks: %s.\n", rate(1e3), rate(2e3), rate(5e3), vsle[1],
        join([@sprintf("%.1f", vsle[findmin(abs.(tt .- (tt[1] + k*5e3)))[2]]) for k in 0:6 if tt[1] + k*5e3 <= tt[end]], " ")))
    println(io, "| reference | V ref [1e6 km3] | V_af ref [m SLE] | A ref [1e6 km2] | dV_af (model - ref) [m] | dA [1e6 km2] | model-only / ref-only area [1e6 km2] | thickness bias / rms, union [m] | bias / rms, overlap [m] | GL dist model->ref mean / median / p90 [km] | ref->model mean [km] |")
    println(io, "|---|---|---|---|---|---|---|---|---|---|---|")
    flds = Dict{String,Matrix{Float64}}("H_model" => Hg_m, "z_bed_model" => zb, "ice_model" => Float64.(Hm .> 1))   # ice_*: ice mask (H > 1 m, grounded and floating)
    case == "pd" && (flds["ice_ref"] = Float64.(clean(rd(topo, "H_ice")) .> 1))
    for (nm, H, g, zbr, zlr) in refs
        sr = stats(H, g, zbr, fill(zlr, size(H)))
        un = gm .| g; ov = gm .& g; d = Hg_m .- H
        dmo = gl_dist(gm, g, xc, yc); dro = gl_dist(g, gm, xc, yc)
        println(io, @sprintf("| %s | %.3f | %.2f | %.3f | %+.2f | %+.3f | %.3f / %.3f | %+.1f / %.1f | %+.1f / %.1f | %.0f / %.0f / %.0f | %.0f |", nm, sr.V, sr.Vaf, sr.A, s0.Vaf - sr.Vaf, s0.A - sr.A,
            sum(area .* (gm .& .!g)) / 1e12, sum(area .* (g .& .!gm)) / 1e12, mean(d[un]), sqrt(mean(d[un] .^ 2)), mean(d[ov]), sqrt(mean(d[ov] .^ 2)),
            isempty(dmo) ? NaN : mean(dmo), isempty(dmo) ? NaN : median(dmo), isempty(dmo) ? NaN : quantile(dmo, 0.9), isempty(dro) ? NaN : mean(dro)))
        flds["H_" * replace(nm, r"[ .\-]" => "_")] = H
    end
    # regional split (NH: boxes on lon/lat; ANT: none)
    if dom == "North"
        grl = (lat .>= 58) .& (lon .>= -75) .& (lon .<= -10)
        na = (lon .>= -170) .& (lon .<= -50) .& (lat .>= 35) .& .!grl
        eu = .!grl .& .!na .& (lat .>= 45)
        println(io, "\nRegions (boxes, approximate; GRL: lat >= 58, lon -75..-10; NA: lon -170..-50, lat >= 35 without GRL; EUR: rest north of 45 N): V_af [m SLE] model | refs\n")
        println(io, "| region | model | " * join([r[1] for r in refs], " | ") * " |"); println(io, "|---|---|" * repeat("---|", length(refs)))
        for (rn, rm) in (("GRL", grl), ("NA", na), ("EUR", eu))
            row = [@sprintf("%.2f", sle(Hg_m, gm .& rm, zb, zsl))]
            for (nm, H, g, zbr, zlr) in refs; push!(row, @sprintf("%.2f", sle(H, g .& rm, zbr, fill(zlr, size(H))))); end
            println(io, "| $rn | " * join(row, " | ") * " |")
        end
        # regional V_af of every 2-D record (<prefix>_reg_ts.tsv): where the ice sheet drifts
        ds = NCDataset(fice); t2 = Float64.(ds["time"][:])
        open(pre * "_reg_ts.tsv", "w") do f
            println(f, "time\tGRL\tNA\tEUR")
            for j in eachindex(t2)
                r = Dict(v => clean(Float64.(nomissing(ds[v][:, :, j], NaN))) for v in ("H_ice", "f_grnd", "z_bed", "z_sl"))
                g = (r["f_grnd"] .>= 0.5) .& (r["H_ice"] .> 0); H = r["H_ice"] .* g
                println(f, @sprintf("%g\t%.3f\t%.3f\t%.3f", t2[j], [sle(H, g .& rm, r["z_bed"], r["z_sl"]) for rm in (grl, na, eu)]...))
            end
        end
        close(ds)
    end
    # region split of the ANT: WAIS / EAIS / APIS from the `regions` of the output if present
    open(pre * ".md", "w") do f; write(f, String(take!(io))); end
    open(pre * "_ts.tsv", "w") do f
        println(f, "time\tV_sle\tV_ice_g\tA_ice_g")
        for j in 1:nT; (j == 1 || j == nT || mod(tt[j] - tt[1], 500.0) == 0) && println(f, @sprintf("%g\t%.4f\t%.5g\t%.5g", tt[j], vsle[j], vg[j], ag[j])); end
    end
    nc = NCDataset(pre * "_fields.nc", "c"); nx, ny = size(Hg_m)
    defDim(nc, "xc", nx); defDim(nc, "yc", ny)
    defVar(nc, "xc", xc, ("xc",)); defVar(nc, "yc", yc, ("yc",))
    for (k, v) in flds; defVar(nc, k, Float32.(v), ("xc", "yc")); end
    defVar(nc, "f_grnd_model", Float32.(fg), ("xc", "yc")); defVar(nc, "lat", Float32.(lat), ("xc", "yc")); defVar(nc, "lon", Float32.(lon), ("xc", "yc"))
    close(nc)
    return 0
end
exit(main(ARGS))
