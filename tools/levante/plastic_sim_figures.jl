# The figures of docs/simulations.md, from the scores of its tables.
#   julia --project=tools/levante tools/levante/plastic_sim_figures.jl <scores.csv> <outdir>
# scores.csv: plastic_ensemble.jl scores of the references and the ens-v2 members (refs-test: experiments/scores.sh).
# The figures of the Earth search come from plastic_earth_emulate.jl. Writes in outdir:
#   ladder.png   per simulation of the ladder log L (ens-v2: every member) and I per data group against the bar
#   bagge.png    I per data group of PaleoMIST's ice on Bagge's structures 1-6, against their lateral factor
using CairoMakie, DelimitedFiles, Printf

const GROUPS = ["Laurentide", "Greenland", "Fennoscandia", "Barents", "Antarctica", "far field", "ESL"]
const SHORT  = ["Laur.", "Greenl.", "Fenn.", "Bar.", "Antarc.", "FF", "ESL"]
const CUT    = Dict(g => (g == "Greenland" ? 7.88 : 3.0) for g in GROUPS)     # the bar of simulations.md
const LADDER = [                      # (label, series), top to bottom
    ("P1 pub-PMIST", "pub-PMIST"), ("P2 pub-GLAC1D", "pub-GLAC1D"),
    ("I1 ice-PMIST | E1D-PMIST", "ice-PMIST | E1D-PMIST"), ("I2 ice-PMIST | E1D-B21-6", "ice-PMIST | E1D-B21-6"),
    ("   ice-PMIST | E1D-B21-6-layers", "ice-PMIST | E1D-B21-6-layers"), ("I3 ice-PMIST | E3D-B21-6", "ice-PMIST | E3D-B21-6"),
    [("I4 ice-PMIST | E3D-B21-$k", "ice-PMIST | E3D-B21-$k") for k in (1, 2, 3, 4, 5, 7, 8)]...,
    ("I5 ice-GLAC1D | E1D-PMIST", "ice-GLAC1D | E1D-PMIST"), ("I6 ice-GLAC1D | E1D-B21-6-layers", "ice-GLAC1D | E1D-B21-6-layers"),
    ("I7 ice-GLAC1D | E3D-B21-6", "ice-GLAC1D | E3D-B21-6"), ("I8 ice-GLAC1D-9927 | E1D-PMIST", "ice-GLAC1D-9927 | E1D-PMIST"),
    ("   ice-GLAC1D-9927 | E3D-B21-6", "ice-GLAC1D-9927 | E3D-B21-6"),
    ("T1 tau-PMIST | E1D-PMIST", "tau-PMIST | E1D-PMIST"), ("D1 tau-PMIST-PD | E1D-PMIST", "tau-PMIST-PD | E1D-PMIST"),
    ("D2 tau-PMIST-PD | E3D-B21-6", "m000"), ("D3 tau-PMIST-PD | E3D-B21-6 @NH16", "tau-PMIST-PD | E3D-B21-6 @NH16"),
    ("E2 ens-v2, the 50 members", ""), ("   best log L (m048)", "m048"), ("   most passed (m036)", "m036")]
const BAGGE_F = ["0.29", "0.3", "0.4", "0.5", "0.75", "1.0"]           # lateral factor of structures 1-6 (v_<f>_s16)

# (series, group) => (I, logL) of the window rows
function read_scores(f)
    M, h = readdlm(f, ',', Any; header = true); h = String.(vec(h)); c(n) = findfirst(==(n), h)
    out = Dict{Tuple{String,String},Tuple{Float64,Float64}}()
    for i in axes(M, 1)
        M[i, c("window")] == "window" || continue
        I = M[i, c("I")]; out[(string(M[i, c("series")]), string(M[i, c("group")]))] =
            (I isa Number ? Float64(I) : NaN, Float64(M[i, c("logL")]))
    end
    out
end
logL(S, s) = S[(s, "all")][2]
Ig(S, s, g) = get(S, (s, g), (NaN, NaN))[1]       # NaN: no such row (the E2 members' row, a group without data)

function fig_ladder(S, file)
    n = length(LADDER); y = n:-1:1
    fig = Figure(size = (1250, 26n + 120))
    a1 = Axis(fig[1, 1]; xlabel = "log L (higher is better)", yticks = (collect(y), first.(LADDER)), width = 330,
              title = "log-likelihood")
    members = [s for (s, g) in keys(S) if g == "all" && occursin(r"^m\d{3}$", s)]
    ye = y[findfirst(isempty ∘ last, LADDER)]
    scatter!(a1, logL.(Ref(S), members), fill(ye, length(members)); color = (:gray, 0.5), markersize = 6)
    for (lab, ls) in (("P1", :dash), ("P2", :dot))
        vlines!(a1, logL(S, LADDER[lab == "P1" ? 1 : 2][2]); color = :black, linestyle = ls, linewidth = 1)
    end
    r = .!isempty.(last.(LADDER))
    scatter!(a1, logL.(Ref(S), last.(LADDER)[r]), collect(y)[r]; color = :black, markersize = 9)
    a2 = Axis(fig[1, 2]; xticks = (1:7, SHORT), xaxisposition = :top, yticksvisible = false, yticklabelsvisible = false,
              title = "I per data group (blue: within the bar)", width = 7 * 62)
    Z = [Ig(S, s, g) - CUT[g] for g in GROUPS, (_, s) in reverse(LADDER)]
    heatmap!(a2, 1:7, 1:n, clamp.(Z, -8, 8); colormap = Reverse(:RdBu), colorrange = (-8, 8), nan_color = :white)
    for (j, (_, s)) in enumerate(reverse(LADDER)), (i, g) in enumerate(GROUPS)
        v = Ig(S, s, g); isnan(v) || text!(a2, i, j; text = @sprintf("%.1f", v), align = (:center, :center), fontsize = 11)
    end
    linkyaxes!(a1, a2); ylims!(a1, 0.5, n + 1)
    Label(fig[2, 1:2], "Dashed: P1, dotted: P2. Grey dots: the members of ens-v2. Window ≤ 25 ka; bar: I ≤ 3, Greenland ≤ 7.9.";
          fontsize = 12, halign = :left)
    resize_to_layout!(fig); save(file, fig)
end

function fig_bagge(S, file)
    fig = Figure(size = (760, 440))
    a = Axis(fig[1, 1]; xlabel = "Bagge's structure (its lateral factor)", ylabel = "I",
             title = "PaleoMIST's ice on Bagge's structures 1-6 (1-D: the radial mean of 6)",
             xticks = (0:6, ["1-D"; ["$k\n($f)" for (k, f) in enumerate(BAGGE_F)]]))
    hlines!(a, 3; color = :gray, linestyle = :dash)
    for (g, col) in zip(GROUPS[1:6], Makie.wong_colors())
        v = [Ig(S, "ice-PMIST | E3D-B21-$k", g) for k in 1:6]
        scatterlines!(a, 1:6, v; color = col, label = g)
        scatter!(a, [0.0], [Ig(S, "ice-PMIST | E1D-B21-6", g)]; color = col, marker = :diamond, markersize = 12)
    end
    Legend(fig[1, 2], a; framevisible = false); save(file, fig)
end

function main(scores, outdir)
    mkpath(outdir); S = read_scores(scores)
    fig_ladder(S, joinpath(outdir, "ladder.png"))
    fig_bagge(S, joinpath(outdir, "bagge.png"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS...)
end
