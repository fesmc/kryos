# Radial viscosity profiles of the Earth search in one figure (docs/science/earth-search.md, "The parameters"): the prior range
# per layer, Gowan's 1-D Earth (PaleoMIST), the area-mean profiles of Bagge's 8 structures (earth basis), and sampled Earths:
# wave A 1-D Earths from the scored cases and wave B 3-D Earths (their radial profile P and the area mean of the 3-D field).
#   julia --project=tools/levante tools/levante/plastic_earth_profiles.jl <basis.nc> <waveA cases.csv> <earthsB dir> <out.png> \
#         A=E165,E117,E050,E120 B=E014
# The elastic lithosphere is drawn at log10 eta = 24 ("elastic"); 3-D fields are clamped to [19.5, 24] before the area mean.
using CairoMakie, NCDatasets, Printf, Statistics

basisf, casesf, dirB, outf = ARGS[1:4]
opt = Dict(String(split(a, "=")[1]) => split(split(a, "=", limit = 2)[2], ",") for a in ARGS[5:end])
const LO, HI = 19.5, 24.0
const PRIOR = [(50.0, 410.0, 19.5, 21.5, "UM"), (410.0, 670.0, 20.0, 22.5, "TZ"), (670.0, 1500.0, 21.0, 23.5, "LM1"),
               (1500.0, 2891.0, 21.0, 23.5, "LM2")]
const COLS = [:royalblue, :darkorange, :purple, :teal, :olive, :brown]

dep = range(0.0, 2891.0, length = 2892)
radial(θ) = d -> d < θ[1] ? HI : d < 410 ? θ[2] : d < 670 ? θ[3] : d < 1500 ? θ[4] : θ[5]
gowan(d) = d < 120 ? HI : d < 670 ? log10(4e20) : log10(4e22)

# Bagge's structures: the area-mean profiles of the basis
b = NCDataset(basisf)
bdep = Float64.(b["depth"][:]) ./ 1e3
names = [String(first(split(strip(s)))) for s in split(b.attrib["fields"], ";")]
prof = Float64.(b["profile"][:, :])
close(b)

# wave A: the radial parameters of each member (the same under both ices)
waveA = Dict{String, Vector{Float64}}()
for l in readlines(casesf)[2:end]
    f = split(l, ","); waveA[f[2]] = parse.(Float64, f[3:7])
end

# wave B: parameters and the area-mean profile of the 3-D field
function waveB(m)
    NCDataset(joinpath(dirB, m * ".nc")) do e
        θ = parse.(Float64, split(e.attrib["parameters"]))
        lat = Float64.(e["lat"][:]); w = cosd.(lat) ./ sum(cosd.(lat))
        r = Float64.(e["r"][:]) ./ 1e3; E = e["eta"]
        pm = [sum(mean(clamp.(Float64.(E[:, :, k]), LO, HI), dims = 1)[1, :] .* w) for k in eachindex(r)]
        θ, 6371.0 .- r, pm
    end
end

fig = Figure(size = (1200, 800), fontsize = 15)
axs = [Axis(fig[1, i], xlabel = "log10 η [Pa s]", ylabel = i == 1 ? "depth [km]" : "", yreversed = true,
            title = i == 1 ? "a. upper 1000 km" : "b. whole mantle", limits = ((19.3, 24.3), (0, i == 1 ? 1000 : 2891)),
            xticks = ([20, 21, 22, 23, 24], ["20", "21", "22", "23", "elastic"])) for i in 1:2]
for ax in axs
    for (d0, d1, x0, x1, nm) in PRIOR
        poly!(ax, Rect(x0, d0, x1 - x0, d1 - d0), color = (:gray80, 0.45), strokewidth = 0, label = "prior range")
    end
    poly!(ax, Rect(23.75, 50.0, 0.5, 150.0), color = (:gray60, 0.45), strokewidth = 0, label = "prior lithosphere (50-200 km)")
    for (m, nm) in enumerate(names)
        startswith(nm, "B21") || continue
        lines!(ax, clamp.(prof[m, :], LO, HI), bdep, color = (:gray40, 0.8), linewidth = 1, label = "Bagge's 8 structures (area mean)")
    end
    lines!(ax, gowan.(dep), dep, color = :green4, linewidth = 3, linestyle = :dash, label = "Gowan 1-D (PaleoMIST)")
    for (i, m) in enumerate(get(opt, "A", String[]))
        θ = waveA[m]
        lines!(ax, radial(θ).(dep), dep, color = COLS[i], linewidth = 2, label = "wave A " * m)
    end
    for m in get(opt, "B", String[])
        θ, d3, pm = waveB(m)
        lines!(ax, radial(θ).(dep), dep, color = :crimson, linewidth = 2, linestyle = :dot, label = "wave B $m, radial profile")
        lines!(ax, pm, d3, color = :crimson, linewidth = 3.5, label = "wave B $m, 3-D (area mean)")
    end
end
Legend(fig[1, 3], axs[1], unique = true, merge = true, framevisible = false, labelsize = 13)
save(outf, fig, px_per_unit = 2)
println("wrote ", outf)
