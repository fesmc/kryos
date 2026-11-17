# Jumps of the ice mass of coupled runs (experiments/coulomb-band.sh and other kryos plastic runs): per run and instance (ant,
# nor, pat) the step changes ΔM of ice_<i>_ts.nc (M_ice [Gt], one record per coupling step), their maximum |ΔM| (and its age),
# the median |ΔM|, and the alternating pairs: two consecutive steps, both above the threshold thr, of opposite sign (the solve
# flipping between two states and back).
#   julia --project=tools/levante tools/levante/plastic_mass_jumps.jl <out folder> [thr=<Gt>] <label>=<run folder> ...
# thr defaults to 5e4 Gt. Writes jumps.md (one table per instance) and jumps.png (ΔM per step against age, a row per instance).
using CairoMakie, NCDatasets, Printf, Statistics

const INST = [("ant", "Antarctica"), ("nor", "North"), ("pat", "Patagonia")]

out = ARGS[1]; mkpath(out)
thr = 5.0e4; args = ARGS[2:end]
if !isempty(args) && startswith(args[1], "thr=")
    thr = parse(Float64, split(args[1], "=")[2]); args = args[2:end]
end
runs = [(String(split(a, "="; limit = 2)[1]), String(split(a, "="; limit = 2)[2])) for a in args]

# (age [ka] of the end of each step, ΔM [Gt], step [yr]) of one instance of a run; nothing without its file
function steps(dir, i)
    f = joinpath(dir, "ice_$(i)_ts.nc"); isfile(f) || return nothing
    t, M = NCDataset(ds -> (Float64.(ds["time"][:]), Float64.(ds["M_ice"][:])), f)
    (-t[2:end] ./ 1000, diff(M), median(diff(t)))
end
alternating(dM) = count(k -> abs(dM[k]) > thr && abs(dM[k + 1]) > thr && sign(dM[k]) != sign(dM[k + 1]), 1:length(dM) - 1)

io = open(joinpath(out, "jumps.md"), "w")
println(io, @sprintf("Step changes of the ice mass M_ice [Gt]; alternating pairs: consecutive steps both above %.0e Gt, opposite sign.", thr))
fig = Figure(size = (1400, 330 * length(INST)))
for (p, (i, nm)) in enumerate(INST)
    println(io, "\n", nm, ":\n\n| run | step [yr] | max \\|ΔM\\| | at [ka] | median \\|ΔM\\| | steps > thr | alternating pairs |\n|---|---|---|---|---|---|---|")
    ax = Axis(fig[p, 1], xreversed = true, xlabel = "ka", ylabel = "ΔM per step [Gt]", title = nm)
    for (k, (lbl, dir)) in enumerate(runs)
        s = steps(dir, i); s === nothing && continue
        a, dM, dt = s; m = argmax(abs.(dM))
        println(io, @sprintf("| %s | %.0f | %.2e | %.1f | %.1e | %d | %d |", lbl, dt, abs(dM[m]), a[m], median(abs.(dM)),
                             count(abs.(dM) .> thr), alternating(dM)))
        lines!(ax, a, dM, color = Makie.wong_colors()[mod1(k, 7)], linewidth = 1.2, label = lbl)
    end
    hlines!(ax, [-thr, thr], color = :gray, linestyle = :dot)
    p == 1 && Legend(fig[1:length(INST), 2], ax, framevisible = false, labelsize = 11)
end
close(io)
save(joinpath(out, "jumps.png"), fig)
println(read(joinpath(out, "jumps.md"), String))
