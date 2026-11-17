# The LGM calibration of the Coulomb law (tests/plastic_coulomb_calib.f90, one run per grid on the same H_N x lambda_min
# lists): the targets, the feasible (H_N, lambda_min), the c_f per ice sheet and a prior box.
#   julia --project=tools/levante tools/levante/plastic_coulomb_calib.jl <out folder> <calib NH-32KM.nc> <calib ANT-16KM.nc>
# The NH run: pd_region 2 (Greenland, its present-day fit), targets NAIS (1) and EIS (3) at 20 ka; the ANT run:
# pd_region 4. Targets (volume above flotation, m SLE): Antarctic excess over its own (modelled) present day in
# [5, 20], Greenland's in [1, 5], NAIS and EIS reaching their target (the midpoints of [60, 85], [14, 25]: the c_f
# search succeeded); soft: maximum thickness of NAIS, EIS and Greenland <= 4500 m (Antarctica: reported, its observed
# present-day maximum on the grid already exceeds it). Writes report.md, calib_plane.png (chi2/n of the present day
# over the plane, the feasible cells marked, the excess as contours), calib_cf.png (c_f per ice sheet, maximum
# thickness).
using CairoMakie, NCDatasets, Printf

out = ARGS[1]; mkpath(out)
rd(f) = NCDataset(f) do ds
    (HN = Float64.(ds["HN"][:]), lam = Float64.(ds["lam_min"][:]), lcf_pd = Float64.(ds["lcf_pd"][:, :]),
     chi2 = Float64.(ds["chi2n_pd"][:, :]), sle_pd = Float64.(ds["sle_pd"][:, :]), obs = ds["sle_pd_obs"][1],
     lcf = Float64.(ds["lcf"][:, :, :]), sle = Float64.(ds["sle"][:, :, :]), hmax = Float64.(ds["hmax"][:, :, :]),
     tgt = Int.(ds["target_region"][:]), tsle = Float64.(ds["target_sle"][:]))
end
N = rd(ARGS[2]); A = rd(ARGS[3])
N.HN == A.HN && N.lam == A.lam || error("the two runs are not on the same H_N x lambda_min")
HN, lam = N.HN, N.lam; nh, nl = length(HN), length(lam)
ex_G = N.sle[2, :, :] .- N.sle_pd; ex_A = A.sle[4, :, :] .- A.sle_pd
reached(r) = (j = findfirst(==(r), N.tgt); abs.(N.sle[r, :, :] .- N.tsle[j]) .< 0.5)
okN = reached(1); okE = reached(3)
dome = maximum(N.hmax[[1, 2, 3], :, :], dims = 1)[1, :, :]; domeA = A.hmax[4, :, :]
hard = (5 .<= ex_A .<= 20) .& (1 .<= ex_G .<= 5) .& okN .& okE
soft = hard .& (dome .<= 4500)

io = open(joinpath(out, "report.md"), "w")
println(io, "# Coulomb law: LGM calibration (20 ka, no Earth)\n\nObserved present day: Greenland ", @sprintf("%.2f", N.obs),
        " m SLE, Antarctica ", @sprintf("%.2f", A.obs), " m SLE. Targets: Antarctic excess 5-20 m, Greenland 1-5 m, NAIS ",
        N.tsle[findfirst(==(1), N.tgt)], " m, EIS ", N.tsle[findfirst(==(3), N.tgt)], " m (midpoints), maximum thickness of NAIS, EIS, Greenland <= 4500 m (soft).\n")
println(io, "| H_N | λ_min | χ²/n GrIS | χ²/n AIS | c_f GrIS | c_f AIS | c_f NAIS | c_f EIS | GrIS excess | AIS excess | NAIS | EIS | max H NH [m] | max H AIS [m] | feasible |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for il in 1:nl, ih in 1:nh
    @printf(io, "| %g | %g | %.2f | %.2f | %.4f | %.4f | %.4f | %.4f | %.1f | %.1f | %.1f | %.1f | %.0f | %.0f | %s |\n", HN[ih], lam[il],
            N.chi2[ih, il], A.chi2[ih, il], exp(N.lcf_pd[ih, il]), exp(A.lcf_pd[ih, il]), exp(N.lcf[1, ih, il]),
            exp(N.lcf[3, ih, il]), ex_G[ih, il], ex_A[ih, il], N.sle[1, ih, il], N.sle[3, ih, il], dome[ih, il], domeA[ih, il],
            soft[ih, il] ? "all" : hard[ih, il] ? "hard (dome > 4.5 km)" : "")
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
            rng(exp.(A.lcf_pd)), ", NAIS ", rng(exp.(N.lcf[1, :, :])), ", EIS ", rng(exp.(N.lcf[3, :, :])),
            "; χ²/n GrIS ", rng(N.chi2), ", AIS ", rng(A.chi2), ".")
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
Label(fig[2, 1:4], "red star: every target met; red circle: the hard targets, dome > 4.5 km", fontsize = 12)
save(joinpath(out, "calib_plane.png"), fig)

fig = Figure(size = (1500, 380))
for (p, (nm, v)) in enumerate((("c_f Greenland (PD)", exp.(N.lcf_pd)), ("c_f Antarctica (PD)", exp.(A.lcf_pd)),
                               ("c_f NAIS (target)", exp.(N.lcf[1, :, :])), ("c_f EIS (target)", exp.(N.lcf[3, :, :])),
                               ("max thickness NAIS, EIS, GrIS at 20 ka [m]", dome)))
    ax = Axis(fig[1, 2p - 1], xticks = xt, yticks = yt, xlabel = "H_N [m]", ylabel = "λ_min", title = nm, xticklabelrotation = 0.8)
    hm = heatmap!(ax, 1:nh, 1:nl, p <= 4 ? log10.(v) : v, colormap = p <= 4 ? :plasma : :thermal)
    for c in CartesianIndices(hard)
        hard[c] && scatter!(ax, [c[1]], [c[2]], marker = soft[c] ? :star5 : :circle, color = :cyan, markersize = soft[c] ? 14 : 8)
    end
    Colorbar(fig[1, 2p], hm, label = p <= 4 ? "log10" : "m")
end
save(joinpath(out, "calib_cf.png"), fig)
println(read(joinpath(out, "report.md"), String))
