# Global mean temperature of runs on imposed geometry (par/kryos_NH-ANT_data_paleoemu*.nml) against their geometry.
# Per run, from clim_glob_emu.nc (the emulator on its CESM grid, area weighted, annual means of the months):
#   dGMT*   = the emulator's dT* (at the reference surface)
#           + bias_scale w B, the bias term (B the global annual mean of the bias file's dtas_b,
#             w = max(0, (dV North America + dV Eurasia) / dV_ref), as paleoemu_bias_weight)
#   dGMT    = dGMT* - Gamma dz, at the surface of the time (dz the emulator's elevation change, max(z_srf, 0) minus
#             the reference's; Gamma = 6.5e-3 K/m)
# all against the pre-industrial reference. Check: dGMT* from clim_glob.nc (tas on GLOBAL-0.5DEG, after the bias)
# minus the reference file's tas_ref. Geometry: the emulator's dV per region (m SLE) and global-mean dz, and from ice_nor.nc and
# ice_ant.nc the ice area and the mean surface elevation of the ice.
#   kryos_jl tools/levante/geom_gmt.jl <out prefix> <label>=<run dir> ...        (tools/levante/env.sh)
# Writes <out prefix>.csv (every output time of every run), <out prefix>.png and prints the table at 20 ka.
using NCDatasets, CairoMakie, Printf, Statistics

const GAMMA = 6.5e-3
out = ARGS[1]
runs = [(String(a[1:findfirst('=', a)-1]), a[findfirst('=', a)+1:end]) for a in ARGS[2:end]]

rd(f, v) = NCDataset(ds -> Float64.(coalesce.(Array(ds[v]), NaN)), f)
function nml_value(R, key)              # the value of key in the run's par file (first match), as a string
    f = only(filter(f -> endswith(f, ".nml") && startswith(f, "kryos_"), readdir(R)))
    m = match(Regex("^\\s*" * key * "\\s*=\\s*[\"']?([^\"'!\\s]+)", "m"), read(joinpath(R, f), String))
    m === nothing ? "" : m[1]
end
amean(X, w) = sum(X .* w) / sum(w)       # area mean of a lon x lat field with weights w

function analyse(R)
    e = joinpath(R, "clim_glob_emu.nc")
    t = rd(e, "time"); lat = rd(e, "lat"); w = repeat(cosd.(lat)', length(rd(e, "lon")), 1)
    dV = rd(e, "dV")                                      # [region, time]
    dz = rd(e, "dz"); dtas = rd(e, "dtas")                # [lon, lat, time], [lon, lat, month, time]
    nt = length(t)
    dTs = [amean(dropdims(mean(dtas[:, :, :, k]; dims=3); dims=3), w) for k in 1:nt]
    dzm = [amean(dz[:, :, k], w) for k in 1:nt]
    bfile = nml_value(R, "bias"); bs = parse(Float64, nml_value(R, "bias_scale"))
    B, dVref = NCDataset(joinpath(R, bfile)) do ds
        b = Float64.(coalesce.(Array(ds["dtas_b"]), NaN))
        blat = Float64.(ds["lat"][:]); bw = repeat(cosd.(blat)', size(b, 1), 1)
        amean(dropdims(mean(b; dims=3); dims=3), bw), Float64(ds.attrib["dV_ref"])
    end
    wgt = max.(0.0, (dV[1, :] .+ dV[2, :]) ./ dVref)
    dTb = bs .* wgt .* B
    dTl = -GAMMA .* dzm
    # check: tas of clim_glob.nc (after the bias) minus the reference's tas_ref, area means of the annual means
    ann(X) = ndims(X) == 3 ? dropdims(mean(X; dims=3); dims=3) : X
    c = joinpath(R, "clim_glob.nc")
    chk = NCDataset(c) do ds
        clat = Float64.(ds["lat"][:]); cw = repeat(cosd.(clat)', length(ds["lon"]), 1)
        Tref = amean(ann(rd(joinpath(R, nml_value(R, "reference")), "tas_ref")), cw)
        [amean(ann(Float64.(ds["tas"][:, :, :, k])), cw) - Tref for k in 1:nt]
    end
    # ice of the domains (the run's data ice): area [1e6 km2] and mean surface elevation of the ice [m]
    dom = Dict{String,Any}()
    for d in ("nor", "ant")
        f = joinpath(R, "ice_$d.nc"); ti = rd(f, "time"); H = rd(f, "H_ice"); zs = rd(f, "z_srf"); x = rd(f, "xc")
        A = ((x[2] - x[1]) * 1e3)^2
        area = [count(H[:, :, k] .> 0) * A / 1e12 for k in eachindex(ti)]
        zice = [(m = H[:, :, k] .> 0; any(m) ? mean(zs[:, :, k][m]) : NaN) for k in eachindex(ti)]
        idx = [argmin(abs.(ti .- tv)) for tv in t]        # the ice output nearest each climate time
        dom[d] = (area=area[idx], zice=zice[idx])
    end
    (t=t, dTs_emu=dTs, dTb=dTb, dTl=dTl, dGs=dTs .+ dTb, dG=dTs .+ dTb .+ dTl, chk=chk, w=wgt,
     dV_NA=dV[1, :], dV_EU=dV[2, :], dV_ANT=dV[3, :], dz=dzm,
     A_nor=dom["nor"].area, z_nor=dom["nor"].zice, A_ant=dom["ant"].area, z_ant=dom["ant"].zice)
end

res = [(lab, analyse(R)) for (lab, R) in runs]

cols = [:t, :dG, :dGs, :dTs_emu, :dTb, :dTl, :chk, :w, :dV_NA, :dV_EU, :dV_ANT, :dz, :A_nor, :z_nor, :A_ant, :z_ant]
open(out * ".csv", "w") do io
    println(io, "run,", join(string.(cols), ","))
    for (lab, r) in res, k in eachindex(r.t)
        println(io, lab, ",", join([@sprintf("%.4g", getfield(r, c)[k]) for c in cols], ","))
    end
end

println("At 20 ka (K against pre-industrial; dV in m SLE, A in 1e6 km2, z in m):")
@printf("%-10s %7s %7s %7s %7s %7s %7s %6s %7s %7s %7s %6s %6s %6s %6s\n", "run", "dGMT", "dGMT*", "emu", "bias",
        "lapse", "check", "w", "dV_NA", "dV_EU", "dV_ANT", "dz", "A_NH", "zi_NH", "zi_ANT")
for (lab, r) in res
    k = argmin(abs.(r.t .+ 20e3))
    @printf("%-10s %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %6.2f %7.1f %7.1f %7.1f %6.1f %6.2f %6.0f %6.0f\n", lab,
            r.dG[k], r.dGs[k], r.dTs_emu[k], r.dTb[k], r.dTl[k], r.chk[k], r.w[k], r.dV_NA[k], r.dV_EU[k],
            r.dV_ANT[k], r.dz[k], r.A_nor[k], r.z_nor[k], r.z_ant[k])
end

palette = Makie.wong_colors()
fig = Figure(size=(1300, 1150), fontsize=14)
ka(r) = -r.t ./ 1e3
function panel(pos, ylab, series; ttl="")
    ax = Axis(fig[pos...], xlabel="age [ka]", ylabel=ylab, title=ttl, xreversed=true)
    for (i, (lab, r)) in enumerate(res), (f, ls) in series
        lines!(ax, ka(r), getfield(r, f), color=palette[i], linestyle=ls, label=lab)
    end
    ax
end
ax1 = panel((1, 1), "ΔGMT [K]", [(:dG, :solid), (:dGs, :dash)]; ttl="ΔGMT at the surface (solid) and at the reference surface (dashed)")
axislegend(ax1, unique=true, position=:lb)
panel((1, 2), "contribution [K]", [(:dTs_emu, :solid), (:dTb, :dash), (:dTl, :dot)]; ttl="emulator dT* (solid), bias (dashed), lapse −Γ dz (dotted)")
panel((2, 1), "ΔV [m SLE]", [(:dV_NA, :solid), (:dV_EU, :dash), (:dV_ANT, :dot)]; ttl="emulator ice volume: North America (solid), Eurasia (dashed), Antarctica (dotted)")
panel((2, 2), "z [m]", [(:z_nor, :solid), (:z_ant, :dash)]; ttl="mean surface elevation of the ice: North (solid), Antarctica (dashed)")
panel((3, 1), "area [10⁶ km²]", [(:A_nor, :solid), (:A_ant, :dash)]; ttl="ice area: North (solid), Antarctica (dashed)")
ax6 = Axis(fig[3, 2], xlabel="ΔV North America + Eurasia [m SLE]", ylabel="ΔGMT [K]", title="ΔGMT against NH ice volume")
for (i, (lab, r)) in enumerate(res)
    scatter!(ax6, r.dV_NA .+ r.dV_EU, r.dG, color=palette[i], markersize=5, label=lab)
end
save(out * ".png", fig)
println("wrote $(out).csv, $(out).png")
