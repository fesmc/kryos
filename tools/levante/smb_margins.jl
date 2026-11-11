# Where does the NH ablation of a run on imposed ice come from? Per region of NH-32KM (FesmData REGIONS, level 2 and 3) at
# 20 ka and 0 ka (the output times present): ∫smb [Gt/yr], the ablation band (grounded ice with smb < -1000 kg m-2 yr-1)
# as a share of the region's ice, its height and latitude, the JJA temperature at the ice surface split into the climate
# at z_srf_clim (the emulator on the ERA5 pre-industrial reference) and the lapse term, the anomaly of that climate
# against 0 ka (NaN without a 0 ka output, e.g. an LGM snapshot), and chion's insolation, albedo, PDDs and snowfall.
#
#   julia --project=tools/levante tools/levante/smb_margins.jl <run dir>       (the North domain: ice_nor, surface_nor,
#                                                                               climate_nor on NH-32KM; $FESMDATA)
using NCDatasets, Statistics, Printf
R = ARGS[1]
F = joinpath(get(ENV, "FESMDATA", "/work/ba1442/fesmdata"), "North/NH-32KM")
rd(f, v) = NCDataset(ds -> Float64.(coalesce.(Array(ds[v]), NaN)), joinpath(R, f))
lat = NCDataset(ds -> Float64.(Array(ds["lat2D"])), "$F/NH-32KM_grid.nc")
reg2 = NCDataset(ds -> Array(ds["region_2"]), "$F/NH-32KM_REGIONS.nc")
reg3 = NCDataset(ds -> Array(ds["region_3"]), "$F/NH-32KM_REGIONS.nc")
regions = [("Laurentide (N America excl. CAA, Cordillera)", (reg2 .== 101) .& .!in.(reg3, Ref((10101, 10103)))),
           ("Cordillera and Alaska", reg3 .== 10103), ("Canadian Arctic Archipelago", reg3 .== 10101),
           ("Greenland", reg2 .== 103), ("Fennoscandia", reg3 .== 10204), ("Barents-Kara", reg3 .== 10201),
           ("British Isles", reg3 .== 10202), ("Svalbard", reg3 .== 10203), ("Asia", reg2 .== 104)]
ti = rd("ice_nor.nc", "time"); ts = rd("surface_nor.nc", "time"); tc = rd("climate_nor.nc", "time")
H, fg, zs = rd("ice_nor.nc", "H_ice"), rd("ice_nor.nc", "f_grnd"), rd("ice_nor.nc", "z_srf")
smb, melt, sf, S, alb, pdd = (rd("surface_nor.nc", v) for v in ("smb", "melt", "sf", "S", "alb_s", "PDDs"))
tas, lap, zc = rd("climate_nor.nc", "tas"), rd("climate_nor.nc", "lapse"), rd("climate_nor.nc", "z_srf_clim")
x = rd("ice_nor.nc", "xc"); A = ((x[2] - x[1]) * 1e3)^2
jja(a, k) = dropdims(mean(a[:, :, 6:8, k]; dims=3); dims=3)       # tas(xc, yc, month, time)
k0c = findfirst(==(0.0), tc)                                     # nothing in a snapshot run: dT0ka NaN
for tk in (-20000.0, 0.0)
    ki, ks, kc = findfirst(==(tk), ti), findfirst(==(tk), ts), findfirst(==(tk), tc)
    any(isnothing, (ki, ks, kc)) && continue
    G = (H[:, :, ki] .> 0) .& (fg[:, :, ki] .> 0.5)
    Tc = jja(tas, kc) .- 273.15; dz = zs[:, :, ki] .- zc[:, :, kc]; Tl = .-jja(lap, kc) .* dz; Ts = Tc .+ Tl
    dT = isnothing(k0c) ? fill(NaN, size(Tc)) : jja(tas, kc) .- jja(tas, k0c)
    @printf("\n== %.0f ka   (band: grounded ice with smb < -1000 kg m-2 yr-1; JJA T at the ice surface = T(z_clim) + lapse term)\n", -tk / 1e3)
    @printf("%-46s %9s %9s %8s %7s %6s %6s %7s %7s %7s %6s %6s %6s %6s\n", "region", "area Mkm2", "smb Gt/yr", "band %", "band z",
            "lat", "T_JJA", "T(zc)", "lapse", "dT0ka", "S", "alb", "PDD", "sf")
    for (name, m) in regions
        g = G .& m; b = g .& (smb[:, :, ks] .< -1000)
        count(g) == 0 && continue
        mb(a) = count(b) > 0 ? mean(a[b]) : NaN
        @printf("%-46s %9.2f %9.0f %8.1f %7.0f %6.1f %6.1f %7.1f %7.1f %7.1f %6.0f %6.2f %6.0f %6.0f\n", name, count(g) * A / 1e12,
                sum(smb[:, :, ks][g]) * A * 1e-12, 100 * count(b) / count(g), mb(zs[:, :, ki]), mb(lat), mb(Ts), mb(Tc), mb(Tl),
                mb(dT), mb(S[:, :, ks]), mb(alb[:, :, ks]), mb(pdd[:, :, ks]), mb(sf[:, :, ks]))
    end
end
