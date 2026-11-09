# Where the sea-level data constrain the ice and where they do not: the spread of the surface elevation across an ensemble of a tau_b
# calibration (plastic_esmda.jl), all members against the members that fit the data about equally well (docs/dev/plastic-backend.md).
#   julia --project=<julia-plot> tools/levante/plastic_ensemble_spread.jl <calibration folder> <iteration> <2-D folder> <out folder>
#         [age=20000] [ngood=10]
# <2-D folder>/mNNN: the members of the iteration rerun with ice_nor.nc, ice_ant.nc (2-D, the 2.5 kyr records); their scores are those
# of the calibration (iter_kk/members.csv: chi2 / n = d' R^-1 d / n over the assimilated index points, R with the correlated model
# error). "Good": the ngood members with the lowest chi2 / n. Surface elevation above the sea level of the age, z_srf - z_sl, on the
# cells with ice in every member (the margin is prescribed, so these are the ice cells). Writes to <out>:
#   spread_<age>ka.png   per grid: mean of the good members; sd of all members; sd of the good members; their ratio (1: the data do not
#                        narrow the spread; small: they do)
#   members.png          chi2 / n against the ice volume per scale group at the age and against the ESL, good members marked
#   spread_<age>ka.md    the numbers: the members by chi2 / n, the effective ensemble size of likelihood weights exp(-n chi2 / 2)
#                        and of tempered ones, the median sd per grid
using CairoMakie, NCDatasets, Printf, Statistics, DelimitedFiles

C, k, D2, out = ARGS[1], parse(Int, ARGS[2]), ARGS[3], ARGS[4]
opt = Dict(split(a, "=", limit = 2) for a in ARGS[5:end])
age = parse(Float64, get(opt, "age", "20000")); ngood = parse(Int, get(opt, "ngood", "10"))
mkpath(out)
const RHO_I, RHO_SW, A_OCN = 910.0, 1028.0, 3.618e14
const FES = "/work/ba1442/fesmdata"
const GRIDS = [("ice_nor", "North", "NH-32KM"), ("ice_ant", "Antarctica", "ANT-16KM")]

# --- members and their scores ---------------------------------------------------------------------------------------------------------
M = readdlm(joinpath(C, @sprintf("iter_%02d", k), "members.csv"), ',', Any; skipstart = 1)
chi = Dict(Int(r[1]) => Float64(r[3]) for r in eachrow(M) if r[2] == 1 && Int(r[1]) != 0)
esl = Dict(Int(r[1]) => Float64(r[4]) for r in eachrow(M) if r[2] == 1)
mems = sort([m for m in keys(chi) if isfile(joinpath(D2, @sprintf("m%03d", m), "status")) &&
             split(read(joinpath(D2, @sprintf("m%03d", m), "status"), String))[1] == "0"])
length(mems) >= ngood + 2 || error("only $(length(mems)) members with 2-D output")
order = sort(mems, by = m -> chi[m]); good = order[1:ngood]
nobs = count(==(1), readdlm(joinpath(C, "obs_kept.csv"), ',', Any; skipstart = 2)[:, end])

# --- surface elevation at the age -----------------------------------------------------------------------------------------------------
function field(m, inst)
    NCDataset(joinpath(D2, @sprintf("m%03d", m), "$(inst).nc")) do ds
        t = Float64.(ds["time"][:]); i = findfirst(x -> abs(x + age) < 1, t)
        g(v) = Float64.(nomissing(ds[v][:, :, i], NaN))
        H = g("H_ice"); zs = g("z_srf") .- g("z_sl"); zs[H .<= 1] .= NaN
        zs, H, Float64.(ds["xc"][:]), Float64.(ds["yc"][:])
    end
end
gname = Dict(Int(r[1]) => String(r[6]) for r in eachrow(readdlm(joinpath(C, "groups.csv"), ',', Any; skipstart = 1)))
groupmap(dom, grid) = NCDataset(joinpath(FES, dom, grid, "$(grid)_IceHistory-PaleoMIST_taub.nc")) do ds
    [ismissing(d) ? "other" : get(gname, Int(d), "other") for d in ds["taub_domain"][:, :]]
end
res = Dict()
vol = Dict{String,Dict{Int,Float64}}()
for (inst, dom, grid) in GRIDS
    Z = Dict(); x = y = nothing; G = groupmap(dom, grid)
    for m in mems
        zs, H, x, y = field(m, inst); Z[m] = zs
        w = (x[2] - x[1]) * (y[2] - y[1]) * 1e6 * RHO_I / RHO_SW / A_OCN
        for g in unique(G)
            g == "other" && continue
            get!(vol, g, Dict{Int,Float64}())[m] = get(get(vol, g, Dict()), m, 0.0) + sum(H[(G .== g) .& (H .> 1)]) * w
        end
    end
    A = cat([Z[m] for m in mems]..., dims = 3); Ag = cat([Z[m] for m in good]..., dims = 3)
    ice = dropdims(all(isfinite, A, dims = 3), dims = 3)
    sdall = dropdims(std(A, dims = 3), dims = 3); sdg = dropdims(std(Ag, dims = 3), dims = 3)
    mg = dropdims(mean(Ag, dims = 3), dims = 3)
    for F in (sdall, sdg, mg); F[.!ice] .= NaN; end
    res[grid] = (x = x, y = y, sdall = sdall, sdg = sdg, mg = mg, ratio = sdg ./ sdall, ice = ice)
end

# --- figures ------------------------------------------------------------------------------------------------------------------------
set_theme!(fontsize = 20)
fig = Figure(size = (2200, 1100))
for (r, (_, _, grid)) in enumerate(GRIDS)
    R = res[grid]; xk, yk = R.x ./ 1000, R.y ./ 1000
    panels = [("$grid: surface elevation, mean of the $ngood best [m]", R.mg, :grays, (0, 4000)),
              ("sd, all $(length(mems)) members [m]", R.sdall, :viridis, (0, quantile(filter(isfinite, R.sdall), 0.98))),
              ("sd, the $ngood best [m]", R.sdg, :viridis, (0, quantile(filter(isfinite, R.sdall), 0.98))),
              ("ratio sd best / sd all", R.ratio, Reverse(:RdYlBu), (0, 1.2))]
    for (c, (t, F, cm, cr)) in enumerate(panels)
        ax = Axis(fig[2r - 1, c], title = t, aspect = DataAspect(), titlesize = 18)
        hidedecorations!(ax)
        hm = heatmap!(ax, xk, yk, F, colormap = cm, colorrange = cr, nan_color = RGBf(0.92, 0.92, 0.92), highclip = cm == :viridis ? :yellow : nothing)
        contour!(ax, xk, yk, Float64.(R.ice), levels = [0.5], color = :black, linewidth = 0.8)
        Colorbar(fig[2r, c], hm, vertical = false, flipaxis = false, width = Relative(0.8))
    end
end
Label(fig[0, :], @sprintf("Surface elevation at %g kyr: spread across the prior ensemble (iteration %d, %d members) and across the %d with the lowest chi2/n",
                          age / 1000, k, length(mems), ngood), fontsize = 24)
save(joinpath(out, @sprintf("spread_%gka.png", age / 1000)), fig)

groups = sort(collect(keys(vol)))
fig = Figure(size = (400 * (length(groups) + 1), 450))
for (c, g) in enumerate([groups; "ESL"])
    ax = Axis(fig[1, c], xlabel = g == "ESL" ? "ESL at 20 kyr [m]" : "$g at $(round(Int, age / 1000)) kyr [m SLE]",
              ylabel = c == 1 ? "chi2 / n" : "")
    xs = g == "ESL" ? [esl[m] for m in mems] : [vol[g][m] for m in mems]
    scatter!(ax, xs, [chi[m] for m in mems], color = [m in good ? :firebrick : :gray50 for m in mems], markersize = 12)
end
save(joinpath(out, "members.png"), fig)

# --- numbers ------------------------------------------------------------------------------------------------------------------------
ess(w) = sum(w)^2 / sum(w .^ 2)
cs = [chi[m] for m in mems]
open(joinpath(out, @sprintf("spread_%gka.md", age / 1000)), "w") do io
    println(io, "Members (iteration $k): $(length(mems)); assimilated index points n = $nobs; chi2/n best ", @sprintf("%.3f", minimum(cs)),
            ", median ", @sprintf("%.3f", median(cs)), ", worst ", @sprintf("%.3f", maximum(cs)), "\n")
    println(io, "Effective ensemble size of weights exp(-n chi2/n / (2 T)) (T = 1: the likelihood with this R; T > 1 tempered):\n")
    println(io, "| T | ESS |\n|---|---|")
    for T in (1, 10, 100, 1000)
        w = exp.(-(nobs / 2T) .* (cs .- minimum(cs))); @printf(io, "| %d | %.2f |\n", T, ess(w))
    end
    println(io, "\nMedian sd of the surface elevation over the ice [m]:\n\n| grid | all | best $ngood | median ratio |\n|---|---|---|---|")
    for (_, _, grid) in GRIDS
        R = res[grid]
        @printf(io, "| %s | %.0f | %.0f | %.2f |\n", grid, median(filter(isfinite, R.sdall)), median(filter(isfinite, R.sdg)),
                median(filter(isfinite, R.ratio)))
    end
    println(io, "\nBest members: ", join([@sprintf("m%03d (%.3f)", m, chi[m]) for m in good], ", "))
end
println("wrote ", out)
