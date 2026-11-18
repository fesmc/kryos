# The calibration surface of the Coulomb law's c_f (taub_method = "coulomb"; docs/backend-plastic.md): c_f*(H_N, lambda_min) per
# region, interpolated from an LGM calibration (tests/plastic_coulomb_calib.f90), the centre of the c_f of the joint ensemble
# (experiments/joint-ens.sh, docs/science/joint-ensemble.md).
#   julia --project=tools/levante tools/levante/plastic_cf_surface.jl at <calib NH-32KM.nc> <calib ANT-16KM.nc> <H_N> <lam_N> <lam_A>
#       c_f* of every region at one point
#   julia --project=tools/levante tools/levante/plastic_cf_surface.jl check <calib NH-32KM.nc> <calib ANT-16KM.nc> <outdir> [ref=<H_N>:<lam_N>:<lam_A>:<c_f x 5>] ...
#       the checks (below): <outdir>/cf_surface.md and cf_surface.png
# The calibration: the north run (pd_region 2, target_region 1, 3, 6: NAIS, Fennoscandia, BarentsKara at time_ka) and the
# Antarctic run (pd_region 4, its present-day fit) of plastic_coulomb_calib.x, each on a grid H_N x lambda_min; the north's
# lambda_min is that of ice_nor and ice_pat (lam_N), the Antarctic one that of ice_ant (lam_A, its zb0 / zb1 those of the run).
# c_f* per region: NAIS (1), Fennoscandia (3), BarentsKara (6) the c_f that gives their target volume at time_ka; Greenland (2)
# the present-day fit of the north run; Antarctica (4) the present-day fit of the Antarctic run (the least chi2 of the
# present-day thickness: not ant-calib's pick of an LGM excess near 12.5 m, which C7 used; the check reports both);
# Patagonia (5) has no calibration and takes Fennoscandia's.
# Interpolation: ln c_f bilinear in (ln H_N, lambda_min) on the grid; a point outside the grid is an error (no extrapolation).
# Checks: (1) the grid nodes are reproduced exactly; (2) every interior node left out along each axis, ln c_f interpolated from
# its two neighbours: the interpolation error at half the grid's spacing bounds it at the nodes' spacing (the table: max |error|
# per region, %); (3) reference points (ref=...: H_N, lam_N, lam_A and the c_f of NAIS, Greenland, Fennoscandia,
# Antarctica, BarentsKara, as the par files of experiments/coulomb-band.sh): the surface against them.
using CairoMakie, NCDatasets, Printf
include("plastic_cf_regions.jl")          # REGION_NAMES, write_cf

const CF_REGIONS = ["NAIS", "Greenland", "Fennoscandia", "BarentsKara", "Antarctica"]   # the calibrated regions, this order
const CF_IDS = Dict("NAIS" => 1, "Greenland" => 2, "Fennoscandia" => 3, "Antarctica" => 4, "BarentsKara" => 6)

# a calibration grid: ln H_N, lambda_min and per region ln c_f (H_N x lambda_min)
struct CfGrid
    lh::Vector{Float64}; lam::Vector{Float64}; lcf::Dict{String,Matrix{Float64}}; file::String
end
function read_calib(f, regions)
    NCDataset(f) do ds
        HN = Float64.(ds["HN"][:]); lam = Float64.(ds["lam_min"][:])
        (issorted(HN) && issorted(lam)) || error("$f: H_N and lambda_min must be increasing")
        rid = haskey(ds, "region_id") ? Int.(ds["region_id"][:]) : collect(1:size(ds["lcf"], 1))
        lcf = Float64.(ds["lcf"][:, :, :]); lpd = Float64.(ds["lcf_pd"][:, :])
        tgt = haskey(ds, "target_region") ? [Int(x) for x in ds["target_region"][:] if x > 0] : Int[]
        out = Dict{String,Matrix{Float64}}()
        for r in regions
            id = CF_IDS[r]
            if id in tgt                                 # a target region: its c_f at time_ka
                out[r] = lcf[findfirst(==(id), rid), :, :]
            else                                         # the present-day fit (pd_region: lcf_pd)
                out[r] = lpd
            end
        end
        CfGrid(log.(HN), lam, out, f)
    end
end
# the surface: the north run for NAIS, Greenland, Fennoscandia, BarentsKara; the Antarctic run for Antarctica
struct CfSurface
    N::CfGrid; A::CfGrid
end
function cf_surface(fN, fA)
    N = read_calib(fN, ["NAIS", "Greenland", "Fennoscandia", "BarentsKara"])
    tN = NCDataset(ds -> [Int(x) for x in ds["target_region"][:] if x > 0], fN)
    sort(tN) == [1, 3, 6] || error("$fN: the north run must have the target regions 1, 3, 6 (NAIS, Fennoscandia, BarentsKara), not $tN")
    NCDataset(ds -> haskey(ds, "target_region") && any(>(0), ds["target_region"][:]), fA) && error("$fA: the Antarctic run has targets")
    CfSurface(N, read_calib(fA, ["Antarctica"]))
end
# bilinear interpolation of M on (x, y) at (a, b), within the grid
function bilin(x, y, M, a, b)
    (x[1] - 1e-12 <= a <= x[end] + 1e-12 && y[1] - 1e-12 <= b <= y[end] + 1e-12) ||
        error(@sprintf("(H_N %.1f, lambda_min %.3f) outside the calibration grid (H_N %.0f-%.0f, lambda_min %.2f-%.2f)",
                       exp(a), b, exp(x[1]), exp(x[end]), y[1], y[end]))
    i = clamp(searchsortedlast(x, a), 1, length(x) - 1); j = clamp(searchsortedlast(y, b), 1, length(y) - 1)
    u = clamp((a - x[i]) / (x[i+1] - x[i]), 0, 1); v = clamp((b - y[j]) / (y[j+1] - y[j]), 0, 1)
    (1 - u) * (1 - v) * M[i, j] + u * (1 - v) * M[i+1, j] + (1 - u) * v * M[i, j+1] + u * v * M[i+1, j+1]
end
# c_f* at (H_N, lam_N, lam_A): region => c_f
function cf_star(S::CfSurface, HN, lamN, lamA)
    d = Dict(r => exp(bilin(S.N.lh, S.N.lam, S.N.lcf[r], log(HN), lamN)) for r in ("NAIS", "Greenland", "Fennoscandia", "BarentsKara"))
    d["Antarctica"] = exp(bilin(S.A.lh, S.A.lam, S.A.lcf["Antarctica"], log(HN), lamA))
    d
end
# the c_f per region id of plastic_cf_regions.jl (1-6) from region => c_f; Patagonia (5) takes Fennoscandia's c_f*
cf_by_id(cf, cfstar) = [cf["NAIS"], cf["Greenland"], cf["Fennoscandia"], cf["Antarctica"], cfstar["Fennoscandia"], cf["BarentsKara"]]

# the leave-one-out along each axis: per region the max |relative error| [%] of the interior nodes, and where
function loo_axes(G::CfGrid, r)
    M = G.lcf[r]; eh = 0.0; el = 0.0; wh = (0, 0); wl = (0, 0)
    for i in 2:length(G.lh)-1, j in axes(M, 2)
        u = (G.lh[i] - G.lh[i-1]) / (G.lh[i+1] - G.lh[i-1]); e = abs(expm1((1 - u) * M[i-1, j] + u * M[i+1, j] - M[i, j]))
        e > eh && (eh = e; wh = (i, j))
    end
    for i in axes(M, 1), j in 2:length(G.lam)-1
        u = (G.lam[j] - G.lam[j-1]) / (G.lam[j+1] - G.lam[j-1]); e = abs(expm1((1 - u) * M[i, j-1] + u * M[i, j+1] - M[i, j]))
        e > el && (el = e; wl = (i, j))
    end
    (h = 100eh, wh = wh, l = 100el, wl = wl)
end

function check(fN, fA, outdir, refs)
    mkpath(outdir); S = cf_surface(fN, fA)
    io = open(joinpath(outdir, "cf_surface.md"), "w")
    println(io, "# c_f* surface\n\nnorth: $fN (H_N ", join(round.(Int, exp.(S.N.lh)), ", "), "; lambda_min ", join(S.N.lam, ", "),
            ")\n\nAntarctica: $fA (H_N ", join(round.(Int, exp.(S.A.lh)), ", "), "; lambda_min ", join(S.A.lam, ", "), ")\n")
    # (1) the nodes
    dn = 0.0
    for (G, rs) in ((S.N, ["NAIS", "Greenland", "Fennoscandia", "BarentsKara"]), (S.A, ["Antarctica"])), r in rs, i in eachindex(G.lh), j in eachindex(G.lam)
        dn = max(dn, abs(bilin(G.lh, G.lam, G.lcf[r], G.lh[i], G.lam[j]) - G.lcf[r][i, j]))
    end
    println(io, @sprintf("Nodes: max |ln c_f(surface) - ln c_f(calibration)| = %.1e\n", dn))
    dn < 1e-12 || error("the surface does not reproduce the nodes")
    # (2) leave-one-out along each axis
    println(io, "Leave-one-out along each axis (an interior node from its two neighbours, twice the grid's spacing): max |c_f error| [%]\n")
    println(io, "| region | along H_N | at (H_N, lambda_min) | along lambda_min | at |\n|---|---|---|---|---|")
    for (G, rs) in ((S.N, ["NAIS", "Greenland", "Fennoscandia", "BarentsKara"]), (S.A, ["Antarctica"])), r in rs
        l = loo_axes(G, r); at(w) = w == (0, 0) ? "-" : @sprintf("(%.0f, %.2f)", exp(G.lh[w[1]]), G.lam[w[2]])
        println(io, @sprintf("| %s | %.2f | %s | %.2f | %s |", r, l.h, at(l.wh), l.l, at(l.wl)))
    end
    # (3) reference points
    if !isempty(refs)
        println(io, "\nReference points: c_f of the reference / c_f* of the surface (difference %)\n")
        println(io, "| H_N | lam_N | lam_A | ", join(CF_REGIONS, " | "), " |\n|---|---|---|", repeat("---|", length(CF_REGIONS)))
        for rf in refs
            v = parse.(Float64, split(rf, ":")); length(v) == 3 + length(CF_REGIONS) || error("ref=$rf: H_N:lam_N:lam_A and $(length(CF_REGIONS)) c_f")
            c = cf_star(S, v[1], v[2], v[3])
            println(io, @sprintf("| %.0f | %.2f | %.2f | ", v[1], v[2], v[3]),
                    join([@sprintf("%.5f / %.5f (%+.2f)", v[3+k], c[r], 100 * (c[r] / v[3+k] - 1)) for (k, r) in enumerate(CF_REGIONS)], " | "), " |")
        end
    end
    close(io)
    print(read(joinpath(outdir, "cf_surface.md"), String))
    # the figure: c_f* against lambda_min per region, a line per H_N node, the surface between the nodes dashed at mid H_N
    fig = Figure(size = (1500, 320))
    for (k, r) in enumerate(CF_REGIONS)
        G = r == "Antarctica" ? S.A : S.N
        ax = Axis(fig[1, k]; title = r, xlabel = r == "Antarctica" ? "λ_min Antarctica" : "λ_min north",
                  ylabel = k == 1 ? "c_f*" : "", yscale = log10)
        cm = cgrad(:viridis, length(G.lh); categorical = true)
        for i in eachindex(G.lh)
            scatterlines!(ax, G.lam, exp.(G.lcf[r][i, :]); color = cm[i], markersize = 5, label = @sprintf("%.0f m", exp(G.lh[i])))
        end
        lf = range(G.lam[1], G.lam[end], length = 60)
        for i in 1:length(G.lh)-1
            hm = (G.lh[i] + G.lh[i+1]) / 2
            lines!(ax, lf, [exp(bilin(G.lh, G.lam, G.lcf[r], hm, b)) for b in lf]; color = cm[i], linestyle = :dash, linewidth = 1)
        end
        k == length(CF_REGIONS) && axislegend(ax; position = :rt, labelsize = 9, framevisible = false, rowgap = 0)
    end
    save(joinpath(outdir, "cf_surface.png"), fig)
end

if abspath(PROGRAM_FILE) == @__FILE__
    cmd = isempty(ARGS) ? "" : ARGS[1]
    if cmd == "at" && length(ARGS) == 6
        c = cf_star(cf_surface(ARGS[2], ARGS[3]), parse.(Float64, ARGS[4:6])...)
        println(join([@sprintf("%s %.6f", r, c[r]) for r in CF_REGIONS], ", "))
    elseif cmd == "check" && length(ARGS) >= 4
        check(ARGS[2], ARGS[3], ARGS[4], [a[5:end] for a in ARGS[5:end] if startswith(a, "ref=")])
    else
        println("usage: plastic_cf_surface.jl at <calib NH.nc> <calib ANT.nc> <H_N> <lam_N> <lam_A> | check <calib NH.nc> <calib ANT.nc> <outdir> [ref=...]")
        exit(1)
    end
end
