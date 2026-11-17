# The ice sheets of coupled runs (experiments/coulomb-band.sh, and other kryos plastic runs) through time against PaleoMIST
# and GLAC-1D: the volume above flotation per ice sheet [m SLE] relative to the present day, every out_2D record.
#   julia --project=tools/levante tools/levante/plastic_coulomb_coupled.jl <out folder> <label>=<run folder> ...
# Ice sheets: the regions of the maps v2 $KRYOS_DATA/plastic/{grid}_taub_regions_v2.nc for every run (geography only, whatever
# map a run used; every cell has a region): NAIS = Laurentide, Fennoscandia (with the British-Irish ice), BarentsKara, Greenland
# on NH-32KM (EIS = Fennoscandia + BarentsKara in the table); Antarctica on ANT-16KM; Patagonia on PAT-PAL-8KM, in the total
# only; a run's ice_nor.nc, ice_ant.nc, ice_pat.nc (H_ice, z_bed, z_sl; any record times); the references the FesmData
# IceHistory products of the grids. SLE = V_af rho_i / (rho_sw A_ocean) (910, 1028 kg m-3, 3.618e14 m2). A Coulomb run labelled
# "tau-<set>-H<H_N>..." (set CB, CBS, CBD, CBF, ...) is coloured by its H_N and styled by its set, thick (the band set CB thin,
# dashed on E1D-ES1-E117); other labels grey dotted. With ice_*_plastic.nc: the mean tau_b over each ice sheet. coupled.md: per series the 20 ka and MIS 3 (30-60 ka, mean of
# the records) volume per ice sheet and the maximum thickness at 20 ka. Writes coupled_esl.png, coupled_taub.png, coupled.md.
using CairoMakie, NCDatasets, Printf, Statistics

const RHO_I = 910.0; const RHO_SW = 1028.0; const A_OC = 3.618e14
const FES = get(ENV, "FESMDATA", "/albedo/work/projects/p_forclima/fesmdata")
const KD = get(ENV, "KRYOS_DATA", "/albedo/work/projects/p_forclima/kryos_data")
const INST = [("nor", "North", "NH-32KM", [(1, "NAIS"), (3, "Fennoscandia"), (6, "BarentsKara"), (2, "Greenland")]),
              ("ant", "Antarctica", "ANT-16KM", [(4, "Antarctica")]), ("pat", "PatagoniaPaleo", "PAT-PAL-8KM", [(5, "Patagonia")])]
const SHEETS = ["NAIS", "Fennoscandia", "BarentsKara", "Greenland", "Antarctica", "Patagonia"]

out = ARGS[1]; mkpath(out)
runs = [(String(split(a, "=")[1]), String(split(a, "="; limit = 2)[2])) for a in ARGS[2:end]]

vaf(H, zb, zsl, m) = sum(max.(H .- (RHO_SW / RHO_I) .* max.(zsl .- zb, 0.0), 0.0)[m .& (H .> 0)]) * RHO_I / (RHO_SW * A_OC)
regions(g) = NCDataset(ds -> Int.(ds["region"][:, :]), joinpath(KD, "plastic", "$(g)_taub_regions_v2.nc"))
dA(ds) = (x = Float64.(ds["xc"][:]); y = Float64.(ds["yc"][:]); (x[2] - x[1]) * (y[2] - y[1]) * 1e6)

# SLE per sheet: Dict(sheet => (t_ka, sle)) of a run, the mean tau_b per sheet, the maximum thickness per sheet at 20 ka
function run_series(dir)
    S = Dict{String,Any}(); TB = Dict{String,Any}(); D = Dict{String,Float64}()
    for (i, dom, g, sh) in INST
        f = joinpath(dir, "ice_$(i).nc"); isfile(f) || continue
        reg = regions(g)
        NCDataset(f) do ds
            t = -Float64.(ds["time"][:]) ./ 1000; a = dA(ds)
            H = Float64.(coalesce.(ds["H_ice"][:, :, :], 0.0)); zb = Float64.(ds["z_bed"][:, :, :]); zs = Float64.(ds["z_sl"][:, :, :])
            k20 = argmin(abs.(t .- 20))
            for (r, s) in sh
                S[s] = (t, [vaf(H[:, :, k], zb[:, :, k], zs[:, :, k], reg .== r) * a for k in eachindex(t)])
                D[s] = maximum(H[:, :, k20][reg .== r]; init = 0.0)
            end
        end
        fp = joinpath(dir, "ice_$(i)_plastic.nc"); isfile(fp) || continue
        NCDataset(fp) do ds
            t = -Float64.(ds["time"][:]) ./ 1000; tb = Float64.(ds["tau_b"][:, :, :]); m = Float64.(ds["mask"][:, :, :]) .> 0.5
            for (r, s) in sh
                TB[s] = (t, [(k = m[:, :, j] .& (reg .== r); any(k) ? mean(tb[:, :, j][k]) / 1e3 : NaN) for j in eachindex(t)])
            end
        end
    end
    S, TB, D
end
function ref_series(product)
    S = Dict{String,Any}()
    for (i, dom, g, sh) in INST
        f = joinpath(FES, dom, g, "$(g)_IceHistory-$(product).nc"); isfile(f) || continue
        reg = regions(g)
        NCDataset(f) do ds
            t = Float64.(ds["time"][:]) ./ 1000; k = findall(x -> 0 <= x <= 80, t); a = dA(ds)
            H = max.(Float64.(coalesce.(ds["H_ice"][:, :, k], 0.0)), 0.0); zb = Float64.(coalesce.(ds["z_bed"][:, :, k], 0.0))
            for (r, s) in sh
                S[s] = (t[k], [vaf(H[:, :, j], zb[:, :, j], zeros(size(reg)), reg .== r) * a for j in eachindex(k)])
            end
        end
    end
    S
end
rel(ts) = (ts[1], ts[2] .- ts[2][argmin(ts[1])])          # relative to the record nearest the present day
function total(S)
    ks = [s for s in SHEETS if haskey(S, s)]; t = S[ks[1]][1]
    (t, sum(rel(S[s])[2] for s in ks))
end

refs = [("PaleoMIST", ref_series("PaleoMIST"), :black, :solid), ("GLAC-1D", ref_series("GLAC1D"), :black, :dash)]
R = [(lbl, run_series(dir)...) for (lbl, dir) in runs]
# Coulomb runs "tau-<set>-H<H_N>[...]": colour by H_N, line style by set (in order of first appearance; the band set "CB"
# thin, dashed on E1D-ES1-E117); any other label grey dotted.
cbm(l) = match(r"^tau-(CB[A-Z0-9]*)-H(\d+)", l)
hns = sort(unique([parse(Int, cbm(l)[2]) for (l, _, _, _) in R if cbm(l) !== nothing]))
cols = Dict(h => Makie.wong_colors()[mod1(k, 7)] for (k, h) in enumerate(hns))
sets = unique([cbm(l)[1] for (l, _, _, _) in R if cbm(l) !== nothing && cbm(l)[1] != "CB"])
const SET_LS = [:solid, :dashdot, :dot, :dashdotdot, (:dash, :dense), (:dot, :loose)]
function style(l)
    m = cbm(l); m === nothing && return (:gray50, :dot, 1.6)
    c = cols[parse(Int, m[2])]
    m[1] == "CB" && return (c, occursin("E117", l) ? :dash : :solid, 1.6)
    (c, SET_LS[mod1(findfirst(==(m[1]), sets), length(SET_LS))], 3.0)
end

fig = Figure(size = (1900, 820))
panels = ["NAIS", "Fennoscandia", "BarentsKara", "Greenland", "Antarctica", "total"]
for (p, s) in enumerate(panels)
    ax = Axis(fig[(p - 1) ÷ 3 + 1, (p - 1) % 3 + 1], xreversed = true, xlabel = "ka", ylabel = "m SLE above present",
              title = s == "total" ? "total (NAIS + EIS + Greenland + Antarctica + Patagonia)" : s)
    vspan!(ax, 30, 60, color = (:gray, 0.12)); vlines!(ax, [20], color = :gray, linestyle = :dot)
    for (nm, S, c, ls) in refs
        (s == "total" || haskey(S, s)) || continue
        t, v = s == "total" ? total(S) : rel(S[s]); lines!(ax, t, v, color = c, linestyle = ls, linewidth = 2.5, label = nm)
    end
    for (l, S, _, _) in R
        (s == "total" || haskey(S, s)) || continue
        t, v = s == "total" ? total(S) : rel(S[s]); c, ls, lw = style(l)
        lines!(ax, t, v, color = c, linestyle = ls, linewidth = lw, label = l)
    end
    p == 1 && Legend(fig[1:2, 4], ax, framevisible = false, labelsize = 11, merge = true, nbanks = 1)
end
save(joinpath(out, "coupled_esl.png"), fig)

if any(!isempty(TB) for (_, _, TB, _) in R)
    fig = Figure(size = (1750, 420))
    for (p, s) in enumerate(SHEETS[1:5])
        ax = Axis(fig[1, p], xreversed = true, xlabel = "ka", ylabel = "mean τ_b on the ice [kPa]", title = s)
        for (l, _, TB, _) in R
            haskey(TB, s) || continue
            c, ls, lw = style(l); lines!(ax, TB[s][1], TB[s][2], color = c, linestyle = ls, linewidth = lw, label = l)
        end
    end
    save(joinpath(out, "coupled_taub.png"), fig)
end

io = open(joinpath(out, "coupled.md"), "w")
at(ts, a) = ts[2][argmin(abs.(ts[1] .- a))]
mis3(ts) = mean(ts[2][(ts[1] .>= 30) .& (ts[1] .<= 60)])
eis(S) = (t = rel(S["Fennoscandia"])[1]; (t, rel(S["Fennoscandia"])[2] .+ rel(S["BarentsKara"])[2]))
cols_md = ["NAIS", "Fennoscandia", "BarentsKara", "EIS", "Greenland", "Antarctica", "total"]
println(io, "Volume above flotation above present [m SLE] at 20 ka and its MIS 3 mean (30-60 ka):\n")
println(io, "| series | " * join(["$s 20 ka | $s MIS 3" for s in cols_md], " | ") * " |\n|---|" * repeat("---|", 2 * length(cols_md)))
for (nm, S) in [[(n, S) for (n, S, _, _) in refs]; [(l, S) for (l, S, _, _) in R]]
    vals = String[]
    for s in cols_md
        (s in ("total", "EIS") || haskey(S, s)) || (push!(vals, "–", "–"); continue)
        ts = s == "total" ? total(S) : s == "EIS" ? eis(S) : rel(S[s]); push!(vals, @sprintf("%.1f", at(ts, 20.0)), @sprintf("%.1f", mis3(ts)))
    end
    println(io, "| ", nm, " | ", join(vals, " | "), " |")
end
println(io, "\nMaximum thickness at 20 ka [m]:\n")
println(io, "| run | " * join(SHEETS[1:5], " | ") * " |\n|---|" * repeat("---|", 5))
for (l, _, _, D) in R
    println(io, "| ", l, " | ", join([haskey(D, s) ? @sprintf("%.0f", D[s]) : "–" for s in SHEETS[1:5]], " | "), " |")
end
close(io)
println(read(joinpath(out, "coupled.md"), String))
