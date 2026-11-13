# NH time slices of a run on imposed ice (par/kryos_NH-ANT_data_paleoemu.nml, North on NH-32KM), MIS 3 (30 ka), LGM
# (21 ka), Early Holocene (9 ka) and pre-industrial (0 ka), the nearest output time to each:
#   mode abs:  JJA near-surface air temperature at the surface (the climate's tas taken to z_srf by its lapse rate),
#              annual precipitation (the climate's) and the SMB on the ice;
#   mode anom: pre-industrial as in abs (first column, its colour bars on the left), then for each slice the JJA
#              temperature change and the precipitation ratio against pre-industrial, and the SMB (not its change).
# Coastline of the time (z_bed = z_sl) and the ice margin drawn over a white halo.
#   kryos_jl tools/levante/nh_slices.jl <run dir> <out.png> <abs|anom>        (tools/levante/env.sh)
using NCDatasets, CairoMakie, Printf, Statistics
R, out, mode = ARGS
rd(f, v) = NCDataset(ds -> Float64.(coalesce.(Array(ds[v]), NaN)), joinpath(R, f))
slices = [(-30e3, "MIS 3"), (-21e3, "LGM"), (-9e3, "Early Holocene"), (0.0, "Pre-industrial")]
ti, ts, tc = rd("ice_nor.nc", "time"), rd("surface_nor.nc", "time"), rd("climate_nor.nc", "time")
H, zb, zs, zsl = (rd("ice_nor.nc", v) for v in ("H_ice", "z_bed", "z_srf", "z_sl"))
smb = rd("surface_nor.nc", "smb")
tas, lap, zc, prc = rd("climate_nor.nc", "tas"), rd("climate_nor.nc", "lapse"), rd("climate_nor.nc", "z_srf_clim"), rd("climate_nor.nc", "pr")
x = rd("ice_nor.nc", "xc"); A = ((x[2] - x[1]) * 1e3)^2
near(t, v) = argmin(abs.(t .- v))
function state(tv)
    ki, ks, kc = near(ti, tv), near(ts, tv), near(tc, tv)
    T = dropdims(mean(tas[:, :, 6:8, kc]; dims=3); dims=3) .- dropdims(mean(lap[:, :, 6:8, kc]; dims=3); dims=3) .*
        (zs[:, :, ki] .- zc[:, :, kc]) .- 273.15
    P = dropdims(mean(prc[:, :, :, kc]; dims=3); dims=3)
    ice = H[:, :, ki] .> 0; land = zb[:, :, ki] .> zsl[:, :, ki]; S = smb[:, :, ks]
    (t=ti[ki], T=T, P=P, S=S, ice=ice, land=land, tot=sum(S[ice .& isfinite.(S)]) * A * 1e-12,
     vol=sum(H[:, :, ki][ice]) * A * 917 / 1000 / 3.618e14)
end
st = [state(tv) for (tv, _) in slices]; pi = st[end]
cols = mode == "abs" ? collect(1:4) : [4, 1, 2, 3]                   # anom: pre-industrial first
lab = Dict("T" => "JJA air temperature [°C]", "P" => "precipitation [mm w.e./yr]", "S" => "SMB [kg m⁻² yr⁻¹]",
           "dT" => "JJA temperature change [K]", "rP" => "precipitation / pre-industrial")
spec(k, s) = k == "T"  ? (s.T, (-30, 25), Reverse(:RdYlBu), identity) :
             k == "P"  ? (max.(s.P, 20), (50, 3000), :YlGnBu, log10) :
             k == "S"  ? (ifelse.(s.ice, s.S, NaN), (-2000, 2000), :RdBu, identity) :
             k == "dT" ? (s.T .- pi.T, (-25, 25), Reverse(:RdBu), identity) :
                         (max.(s.P ./ pi.P, 0.1), (0.25, 4.0), :BrBG, log10)
fig = Figure(size=(330 * 4 + 260, 1000), fontsize=13, backgroundcolor=:white)
off = mode == "abs" ? 0 : 1                                           # anom: a colour-bar column on the left
hms = Dict{String,Any}()
for (j, c) in enumerate(cols)
    s = st[c]; keys_ = mode == "abs" || c == 4 ? ["T", "P", "S"] : ["dT", "rP", "S"]
    for (r, k) in enumerate(keys_)
        a, cr, cm, sc = spec(k, s)
        ax = Axis(fig[r, j + off]; aspect=DataAspect(), backgroundcolor=RGBf(0.86, 0.90, 0.95),
                  title=r == 1 ? @sprintf("%s, %.1f ka", slices[c][2], -s.t / 1e3) : (r == 3 ? @sprintf("∫SMB %.0f Gt/yr", s.tot) : ""),
                  titlesize=r == 1 ? 15 : 12, titlefont=r == 1 ? :bold : :regular)
        hidedecorations!(ax); hidespines!(ax)
        hm = heatmap!(ax, a; colorrange=cr, colormap=cm, colorscale=sc, nan_color=:transparent,
                      lowclip=first(to_colormap(cm)), highclip=last(to_colormap(cm)))
        k == "S" && heatmap!(ax, ifelse.(s.land .& .!s.ice, 1.0, NaN); colormap=[RGBf(0.80, 0.78, 0.74)], nan_color=:transparent)
        coast, margin = Float64.(s.land .| s.ice), Float64.(s.ice)
        contour!(ax, coast; levels=[0.5], color=(:white, 0.9), linewidth=2.0)
        contour!(ax, coast; levels=[0.5], color=:grey15, linewidth=0.7)
        contour!(ax, margin; levels=[0.5], color=(:white, 0.95), linewidth=3.2)
        contour!(ax, margin; levels=[0.5], color=:black, linewidth=1.5)
        r == 1 && text!(ax, 0.03, 0.03; text=@sprintf("%.0f m SLE", s.vol), space=:relative, fontsize=11, color=:black)
        hms[k] = hm
    end
end
cb(pos, k; kw...) = Colorbar(fig[pos...], hms[k]; label=lab[k], width=12,
                             ticks=k == "P" ? [50, 100, 300, 1000, 3000] : k == "rP" ? [0.25, 0.5, 1, 2, 4] : Makie.automatic, kw...)
if mode == "abs"
    for (r, k) in enumerate(["T", "P", "S"]); cb((r, 5), k); end
else
    for (r, k) in enumerate(["T", "P", "S"]); cb((r, 1), k; flipaxis=false); end
    for (r, k) in enumerate(["dT", "rP"]); cb((r, 6), k); end
    colgap!(fig.layout, 1, 10); colgap!(fig.layout, 2, 18)
end
rowgap!(fig.layout, 4)
Label(fig[0, (1+off):(4+off)], mode == "abs" ? "Northern Hemisphere on imposed PaleoMIST ice: PaleoEmulator climate, chion-BESSI SMB" :
      "Change against pre-industrial: JJA temperature, precipitation ratio; SMB as is"; fontsize=18, font=:bold)
save(out, fig; px_per_unit=2)
