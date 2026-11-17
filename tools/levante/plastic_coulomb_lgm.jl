# The ice-sheet sizes of a parameter scan of the Coulomb yield law at one time (tests/plastic_coulomb_scan.f90, no Earth:
# the history's bed of the slice, z_sl = 0) against PaleoMIST and GLAC-1D at that time:
#   julia --project=tools/levante tools/levante/plastic_coulomb_lgm.jl <out folder> <time_ka> <scan.nc> ...
# Per case and ice sheet (the regions of the scan's region map: 1 Laurentide = NAIS, 2 Greenland, 3 Eurasia = EAIS,
# 4 Antarctica, 5 Patagonia): the volume, the volume above flotation and its sea-level equivalent SLE = V_af rho_i /
# (rho_sw A_ocean) (910, 1028 kg m-3, A_ocean = 3.618e14 m2), the maximum thickness. The references on the same grid
# (FesmData IceHistory-PaleoMIST, -GLAC1D of the grid, the slice at time_ka, their own bed): their ice cells outside
# every region take the region of the nearest region cell. Writes <out>/lgm_cases.csv, lgm_report.md and, for a scan
# with several cases, lgm_volume_<grid>.png (SLE against c_f, a colour per H_N, a line style per lambda_min, one panel per
# ice sheet) and lgm_maps_<grid>.png (PaleoMIST, GLAC-1D and three Coulomb cases).
using CairoMakie, NCDatasets, Printf, Statistics

const RHO_I = 910.0; const RHO_SW = 1028.0; const A_OC = 3.618e14
const NAMES = ["NAIS (Laurentide)", "Greenland", "EAIS (Eurasia)", "Antarctica", "Patagonia"]
const FES = get(ENV, "FESMDATA", "/albedo/work/projects/p_forclima/fesmdata")
const DOMAIN = Dict("NH-32KM" => "North", "ANT-16KM" => "Antarctica", "PAT-PAL-8KM" => "PatagoniaPaleo")

out = ARGS[1]; mkpath(out); tka = parse(Float64, ARGS[2]); files = ARGS[3:end]

# volume [m3], volume above flotation [m3], SLE [m], max H [m] over the cells of region r
function sizes(H, zb, reg, r, dA)
    k = (reg .== r) .& (H .> 0)
    hf = (RHO_SW / RHO_I) .* max.(0.0, .-zb)
    V = sum(H[k]) * dA; Vaf = sum(max.(H[k] .- hf[k], 0.0)) * dA
    (V = V, Vaf = Vaf, sle = Vaf * RHO_I / (RHO_SW * A_OC), Hmax = any(k) ? maximum(H[k]) : 0.0)
end

# the region map, with the cells of `extra` outside every region set to the nearest region cell's
function fill_regions(reg, extra)
    r = copy(reg); src = findall(>(0), reg)
    for c in findall(extra .& (reg .<= 0))
        _, k = findmin(s -> (s[1] - c[1])^2 + (s[2] - c[2])^2, src)
        r[c] = reg[src[k]]
    end
    r
end

function reference(grid, product, tka)
    f = joinpath(FES, DOMAIN[grid], grid, "$(grid)_IceHistory-$(product).nc")
    NCDataset(f) do ds
        t = Float64.(ds["time"][:]); it = argmin(abs.(t .- 1000tka))
        abs(t[it] - 1000tka) < 1 || error("$f: no slice at $tka ka")
        H = Float64.(coalesce.(ds["H_ice"][:, :, it], 0.0)); zb = Float64.(coalesce.(ds["z_bed"][:, :, it], 0.0))
        max.(H, 0.0), zb
    end
end

csv = open(joinpath(out, "lgm_cases.csv"), "w")
println(csv, "grid,case,cf,HN,lam_min,region,name,V_1e6km3,Vaf_1e6km3,SLE_m,Hmax_m")
rep = open(joinpath(out, "lgm_report.md"), "w")
for f in files
    ds = NCDataset(f)
    x = Float64.(ds["xc"][:]); y = Float64.(ds["yc"][:]); dA = (x[2] - x[1]) * (y[2] - y[1]) * 1e6
    cf = Float64.(ds["case_cf"][:]); HN = Float64.(ds["case_HN"][:]); lam = Float64.(ds["case_lam_min"][:])
    zb = Float64.(ds["z_bed"][:, :]); reg = Int.(ds["region"][:, :]); H = Float64.(ds["H"][:, :, :])
    close(ds)
    grid = Dict((276, 276) => "NH-32KM", (381, 381) => "ANT-16KM", (131, 433) => "PAT-PAL-8KM")[size(reg)]
    regs = [r for r in 1:5 if any(==(r), reg)]
    Hp, zbp = reference(grid, "PaleoMIST", tka); Hg, zbg = reference(grid, "GLAC1D", tka)
    refs = [("PaleoMIST", Hp, zbp, fill_regions(reg, Hp .> 0)), ("GLAC-1D", Hg, zbg, fill_regions(reg, Hg .> 0))]
    println(rep, "## ", basename(f), " (", grid, ", ", tka, " ka, ", length(cf), " cases)\n")
    for (nm, Hr, zr, rr) in refs, r in regs
        s = sizes(Hr, zr, rr, r, dA)
        @printf(rep, "- %s %s: V %.2f 10^6 km3, V_af %.2f 10^6 km3, %.1f m SLE, max H %.0f m\n", nm, NAMES[r], s.V * 1e-15, s.Vaf * 1e-15, s.sle, s.Hmax)
        @printf(csv, "%s,%s,,,,%d,%s,%.4f,%.4f,%.3f,%.1f\n", grid, nm, r, NAMES[r], s.V * 1e-15, s.Vaf * 1e-15, s.sle, s.Hmax)
    end
    S = Dict(r => [sizes(H[:, :, c], zb, reg, r, dA) for c in eachindex(cf)] for r in regs)
    for c in eachindex(cf), r in regs
        s = S[r][c]
        @printf(csv, "%s,%d,%.5f,%g,%g,%d,%s,%.4f,%.4f,%.3f,%.1f\n", grid, c, cf[c], HN[c], lam[c], r, NAMES[r], s.V * 1e-15, s.Vaf * 1e-15, s.sle, s.Hmax)
    end
    if length(cf) == 1
        println(rep, "\nThe case (c_f ", cf[1], ", H_N ", HN[1], " m, lambda_min ", lam[1], "): ",
                join([@sprintf("%s %.1f m SLE (V %.2f 10^6 km3, max H %.0f m)", NAMES[r], S[r][1].sle, S[r][1].V * 1e-15, S[r][1].Hmax) for r in regs], "; "), "\n")
        continue
    end
    ucf = sort(unique(cf)); uHN = sort(unique(HN)); ulam = sort(unique(lam))
    idx(a, b, c) = findfirst(k -> cf[k] == a && HN[k] == b && lam[k] == c, eachindex(cf))
    sheets = [r for r in (1, 3) if r in regs]
    # sensitivity: the SLE range over each parameter with the others at the middle of their lists
    cm = ucf[cld(length(ucf), 2)]; hm = uHN[cld(length(uHN), 2)]; lm = ulam[cld(length(ulam), 2)]
    println(rep, "\nSLE [m] of the scan (min - max over all cases), and its range over one parameter with the others at c_f ", cm,
            ", H_N ", hm, " m, lambda_min ", lm, ":\n\n| ice sheet | all cases | c_f | H_N | lambda_min | PaleoMIST | GLAC-1D |\n|---|---|---|---|---|---|---|")
    for r in sheets
        v = [s.sle for s in S[r]]
        rg(ks) = (w = [S[r][k].sle for k in ks]; @sprintf("%.1f - %.1f", minimum(w), maximum(w)))
        sp = sizes(refs[1][2], refs[1][3], refs[1][4], r, dA).sle; sg = sizes(refs[2][2], refs[2][3], refs[2][4], r, dA).sle
        @printf(rep, "| %s | %.1f - %.1f | %s | %s | %s | %.1f | %.1f |\n", NAMES[r], minimum(v), maximum(v),
                rg([idx(a, hm, lm) for a in ucf]), rg([idx(cm, b, lm) for b in uHN]), rg([idx(cm, hm, c) for c in ulam]), sp, sg)
    end
    # the cases within the published range [min, max of PaleoMIST and GLAC-1D] of both NAIS and EAIS
    lo = Dict(r => minimum(sizes(R[2], R[3], R[4], r, dA).sle for R in refs) for r in sheets)
    hi = Dict(r => maximum(sizes(R[2], R[3], R[4], r, dA).sle for R in refs) for r in sheets)
    inr = [k for k in eachindex(cf) if all(lo[r] <= S[r][k].sle <= hi[r] for r in sheets)]
    println(rep, "\nCases with every ice sheet (", join(NAMES[sheets], ", "), ") within [min, max] of PaleoMIST and GLAC-1D: ",
            length(inr), " of ", length(cf), "\n")
    for k in inr
        @printf(rep, "- c_f %.4f, H_N %g, lambda_min %g: %s\n", cf[k], HN[k], lam[k], join([@sprintf("%s %.1f m", NAMES[r], S[r][k].sle) for r in sheets], ", "))
    end
    println(rep)

    # figure: SLE against c_f
    cols = Makie.wong_colors(); sty = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    fig = Figure(size = (520 * length(sheets), 460))
    for (p, r) in enumerate(sheets)
        ax = Axis(fig[1, p], xscale = log10, xlabel = "c_f", ylabel = "volume above flotation [m SLE]", title = NAMES[r] * " at $(tka) ka")
        for (ih, b) in enumerate(uHN), (il, c) in enumerate(ulam)
            lines!(ax, ucf, [S[r][idx(a, b, c)].sle for a in ucf], color = cols[ih], linestyle = sty[il])
        end
        for (k, R) in enumerate(refs)
            hlines!(ax, [sizes(R[2], R[3], R[4], r, dA).sle], color = :black, linestyle = k == 1 ? :solid : :dash, linewidth = 2)
        end
    end
    Legend(fig[1, length(sheets) + 1],
           vcat([LineElement(color = cols[i]) for i in eachindex(uHN)], [LineElement(color = :gray, linestyle = sty[i]) for i in eachindex(ulam)],
                [LineElement(color = :black, linewidth = 2), LineElement(color = :black, linestyle = :dash, linewidth = 2)]),
           vcat(["H_N $(Int(b)) m" for b in uHN], ["λ_min $c" for c in ulam], ["PaleoMIST", "GLAC-1D"]))
    save(joinpath(out, "lgm_volume_$(grid).png"), fig)

    # maps: PaleoMIST, GLAC-1D and, at the middle lambda_min, the c_f nearest PaleoMIST's NAIS SLE for the smallest,
    # the middle and the largest H_N
    r0 = sheets[1]; sp0 = sizes(refs[1][2], refs[1][3], refs[1][4], r0, dA).sle
    picks = [argmin(k -> HN[k] == b && lam[k] == lm ? abs(S[r0][k].sle - sp0) : Inf, eachindex(cf)) for b in (uHN[1], hm, uHN[end])]
    fig = Figure(size = (1700, 420))
    panels = vcat([(R[1] * @sprintf(" (NAIS %.0f, EAIS %.0f m)", [sizes(R[2], R[3], R[4], r, dA).sle for r in (1, 3)]...), R[2]) for R in refs],
                  [(@sprintf("c_f %.4f, H_N %g, λ_min %g (NAIS %.0f, EAIS %.0f m)", cf[k], HN[k], lam[k], S[1][k].sle, S[3][k].sle), H[:, :, k]) for k in picks])
    local hm1
    for (p, (t, F)) in enumerate(panels)
        ax = Axis(fig[1, p], title = t, aspect = DataAspect(), titlesize = 11)
        hidedecorations!(ax)
        hm1 = heatmap!(ax, x, y, ifelse.(F .> 0, F, NaN), colormap = :ice, colorrange = (0, 4500))
    end
    Colorbar(fig[1, length(panels) + 1], hm1, label = "H [m]")
    save(joinpath(out, "lgm_maps_$(grid).png"), fig)
end
close(csv); close(rep)
println(read(joinpath(out, "lgm_report.md"), String))
