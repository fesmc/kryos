# The LGM calibration of the Coulomb law (tests/plastic_coulomb_calib.f90, one run per grid on the same H_N x lambda_min
# lists): the targets, the feasible (H_N, lambda_min), the c_f per ice sheet and a prior box.
#   julia --project=tools/levante tools/levante/plastic_coulomb_calib.jl <out folder> <calib NH-32KM.nc> <calib ANT-16KM.nc>
# The NH run: pd_region 2 (Greenland, its present-day fit), target regions at 20 ka NAIS (1) and the Eurasian ones: with the
# v2 maps Fennoscandia (3) and BarentsKara (6), with v1 Eurasia (3) (plastic_cf_regions.jl); the ANT run: pd_region 4.
# Targets (volume above flotation, m SLE): Antarctic excess over its own (modelled) present day in [5, 20], Greenland's in
# [1, 5], every target region reaching its target_sle (the c_f search succeeded, within 0.5 m); soft: maximum thickness of
# the northern ice sheets <= 4500 m, not applied at H_N >= 2000 m (there the domes are reported only; Antarctica: reported,
# its observed present-day maximum on the grid already exceeds it). Writes report.md, calib_plane.png (chi2/n of the present
# day over the plane, the feasible cells marked, the excess as contours), calib_cf.png (c_f per ice sheet, maximum
# thickness).
using CairoMakie, NCDatasets, Printf

const DOME_MAX = 4500.0      # [m] soft target on the northern domes
const DOME_HN = 2000.0       # [m] from this H_N on the dome criterion is dropped (reported only)
const NAMES = Dict(1 => "NAIS", 2 => "GrIS", 3 => "EIS", 4 => "AIS", 5 => "Patagonia", 6 => "BarentsKara")

out = ARGS[1]; mkpath(out)
rd(f) = NCDataset(f) do ds
    (HN = Float64.(ds["HN"][:]), lam = Float64.(ds["lam_min"][:]), lcf_pd = Float64.(ds["lcf_pd"][:, :]),
     chi2 = Float64.(ds["chi2n_pd"][:, :]), sle_pd = Float64.(ds["sle_pd"][:, :]), obs = ds["sle_pd_obs"][1],
     lcf = Float64.(ds["lcf"][:, :, :]), sle = Float64.(ds["sle"][:, :, :]), hmax = Float64.(ds["hmax"][:, :, :]),
     rid = haskey(ds, "region_id") ? Int.(ds["region_id"][:]) : collect(1:size(ds["lcf"], 1)),
     tgt = haskey(ds, "target_region") ? [Int(x) for x in ds["target_region"][:] if x > 0] : Int[],
     tsle = haskey(ds, "target_sle") ? Float64.(ds["target_sle"][:]) : Float64[])
end
N = rd(ARGS[2]); A = rd(ARGS[3])
N.HN == A.HN && N.lam == A.lam || error("the two runs are not on the same H_N x lambda_min")
HN, lam = N.HN, N.lam; nh, nl = length(HN), length(lam)
6 in N.rid && (NAMES[3] = "Fennoscandia")
ix(R, r) = (k = findfirst(==(r), R.rid); k === nothing && error("region $r is not in the run"); k)
ex_G = N.sle[ix(N, 2), :, :] .- N.sle_pd; ex_A = A.sle[ix(A, 4), :, :] .- A.sle_pd
reached(r) = (j = findfirst(==(r), N.tgt); abs.(N.sle[ix(N, r), :, :] .- N.tsle[j]) .< 0.5)
okT = reduce(.&, [reached(r) for r in N.tgt])
north = [r for r in (1, 2, 3, 6) if r in N.rid]
dome = maximum(N.hmax[[ix(N, r) for r in north], :, :], dims = 1)[1, :, :]; domeA = A.hmax[ix(A, 4), :, :]
hard = (5 .<= ex_A .<= 20) .& (1 .<= ex_G .<= 5) .& okT
nodome = repeat(HN .>= DOME_HN, 1, nl)                     # the dome criterion not applied
soft = hard .& ((dome .<= DOME_MAX) .| nodome)
tg = N.tgt                                                 # the target regions, in the run's order

io = open(joinpath(out, "report.md"), "w")
println(io, "# Coulomb law: LGM calibration (20 ka, no Earth)\n\nObserved present day: Greenland ", @sprintf("%.2f", N.obs),
        " m SLE, Antarctica ", @sprintf("%.2f", A.obs), " m SLE. Targets: Antarctic excess 5-20 m, Greenland 1-5 m, ",
        join(["$(NAMES[r]) $(N.tsle[j]) m" for (j, r) in enumerate(tg)], ", "), "; maximum thickness of the northern ice sheets <= ",
        Int(DOME_MAX), " m (soft; not applied at H_N >= ", Int(DOME_HN), " m).\n")
println(io, "| H_N | λ_min | χ²/n GrIS | χ²/n AIS | c_f GrIS | c_f AIS | " * join(["c_f $(NAMES[r]) | " for r in tg]) *
        "GrIS excess | AIS excess | " * join(["$(NAMES[r]) | " for r in tg]) * "max H NH [m] | max H AIS [m] | feasible |\n|" *
        repeat("---|", 11 + 2length(tg)))
for il in 1:nl, ih in 1:nh
    @printf(io, "| %g | %g | %.2f | %.2f | %.4f | %.4f | ", HN[ih], lam[il], N.chi2[ih, il], A.chi2[ih, il], exp(N.lcf_pd[ih, il]),
            exp(A.lcf_pd[ih, il]))
    for r in tg; @printf(io, "%.4f | ", exp(N.lcf[ix(N, r), ih, il])); end
    @printf(io, "%.1f | %.1f | ", ex_G[ih, il], ex_A[ih, il])
    for r in tg; @printf(io, "%.2f | ", N.sle[ix(N, r), ih, il]); end
    @printf(io, "%.0f | %.0f | %s |\n", dome[ih, il], domeA[ih, il],
            soft[ih, il] ? (nodome[ih, il] && dome[ih, il] > DOME_MAX ? "all (dome not applied)" : "all") :
            hard[ih, il] ? "hard (dome > 4.5 km)" : "")
end
for (nm, m) in (("all targets", soft), ("the hard targets", hard))
    k = findall(m)
    if isempty(k)
        println(io, "\nNo (H_N, λ_min) meets ", nm, ".")
        continue
    end
    rng(v) = @sprintf("%.4g - %.4g", minimum(v[k]), maximum(v[k]))
    println(io, "\n", length(k), " of ", nh * nl, " cells meet ", nm, ": H_N ", rng([HN[c[1]] for c in CartesianIndices(m)]),
            " m, λ_min ", rng([lam[c[2]] for c in CartesianIndices(m)]), "; c_f GrIS ", rng(exp.(N.lcf_pd)), ", AIS ",
            rng(exp.(A.lcf_pd)), join([", $(NAMES[r]) " * rng(exp.(N.lcf[ix(N, r), :, :])) for r in tg]),
            "; χ²/n GrIS ", rng(N.chi2), ", AIS ", rng(A.chi2), ".")
end
println(io, "\nMaximum thickness at 20 ka per northern ice sheet [m], H_N >= ", Int(DOME_HN), " m:\n")
println(io, "| H_N | λ_min | " * join(["$(NAMES[r]) | " for r in north]) * "\n|" * repeat("---|", 2 + length(north)))
for il in 1:nl, ih in 1:nh
    HN[ih] >= DOME_HN || continue
    println(io, "| ", @sprintf("%g | %g | ", HN[ih], lam[il]), join([@sprintf("%.0f | ", N.hmax[ix(N, r), ih, il]) for r in north]))
end
close(io)

# figure: the plane (indices, labelled), chi2/n, feasible cells, excess contours
xt = (1:nh, [@sprintf("%g", h) for h in HN]); yt = (1:nl, [@sprintf("%g", l) for l in lam])
fig = Figure(size = (1250, 520))
for (p, (nm, R, ex, lv)) in enumerate((("Greenland", N, ex_G, [1, 3, 5, 10]), ("Antarctica", A, ex_A, [5, 10, 20, 30])))
    ax = Axis(fig[1, 2p - 1], xticks = xt, yticks = yt, xlabel = "H_N [m]", ylabel = "λ_min",
              title = "$nm: present-day χ²/n; contours: 20 ka excess [m SLE]")
    hm = heatmap!(ax, 1:nh, 1:nl, R.chi2, colormap = :viridis, colorrange = (minimum(R.chi2), min(maximum(R.chi2), 4minimum(R.chi2))))
    contour!(ax, 1:nh, 1:nl, ex, levels = lv, color = :white, labels = true, labelsize = 11)
    for c in CartesianIndices(hard)
        hard[c] && scatter!(ax, [c[1]], [c[2]], marker = soft[c] ? :star5 : :circle, color = :red, markersize = soft[c] ? 18 : 10)
    end
    Colorbar(fig[1, 2p], hm, label = "χ²/n")
end
Label(fig[2, 1:4], "red star: every target met (the dome criterion below H_N $(Int(DOME_HN)) m only); red circle: the hard targets, dome > 4.5 km",
      fontsize = 12)
save(joinpath(out, "calib_plane.png"), fig)

pan = [("c_f Greenland (PD)", exp.(N.lcf_pd)), ("c_f Antarctica (PD)", exp.(A.lcf_pd)),
       [("c_f $(NAMES[r]) (target)", exp.(N.lcf[ix(N, r), :, :])) for r in tg]..., ("max thickness, northern ice sheets, 20 ka [m]", dome)]
fig = Figure(size = (300 * length(pan), 380))
for (p, (nm, v)) in enumerate(pan)
    cf = p < length(pan)
    ax = Axis(fig[1, 2p - 1], xticks = xt, yticks = yt, xlabel = "H_N [m]", ylabel = "λ_min", title = nm, xticklabelrotation = 0.8)
    hm = heatmap!(ax, 1:nh, 1:nl, cf ? log10.(v) : v, colormap = cf ? :plasma : :thermal)
    for c in CartesianIndices(hard)
        hard[c] && scatter!(ax, [c[1]], [c[2]], marker = soft[c] ? :star5 : :circle, color = :cyan, markersize = soft[c] ? 14 : 8)
    end
    Colorbar(fig[1, 2p], hm, label = cf ? "log10" : "m")
end
save(joinpath(out, "calib_cf.png"), fig)
println(read(joinpath(out, "report.md"), String))
