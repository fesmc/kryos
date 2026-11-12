# The SMB target of the Greenland PD run on GRL-16KM: MAR v3.14 (ERA5, 1 km monthly) 1981-2010, annual means over the
# ICE of each cell only (MSK >= 2: the ice sheet, 4, and the peripheral ice), binned onto the GRL-16KM cells (the same
# polar stereographic projection; centres xc = -720 + 16(i-1), yc = -3450 + 16(j-1) km): smb, melt, runoff, sf, rf
# [kg m-2 yr-1], z_srf [m]; f_ice, the ice share of the cell (n_ice/256), and f_gris, the ice-sheet share (MSK 4).
# FesmData's remap averages MAR over the whole cell, the tundra too, whose snowfall and seasonal melt are not the ice's.
# Prints the mean of each MSK class. The input: MAR v3.14 1 km monthly on Levante,
# $FESMDATA_SOURCES/mar314_era5_1km_monthly; the output goes to $KRYOS_DATA/smb.
#   julia --project=tools/levante tools/levante/mar_ice16.jl <mar dir> <out.nc>
using NCDatasets, Statistics, Printf
D, out = ARGS
vars = ("SMB" => "smb", "ME" => "melt", "RU" => "runoff", "SF" => "sf", "RF" => "rf")
x, y, msk, srf = NCDataset(ds -> (Float64.(ds["x"][:]), Float64.(ds["y"][:]), Array(ds["MSK"]), Float64.(coalesce.(Array(ds["SRF"]), NaN))),
                           joinpath(D, "MARv3.14-monthly-ERA5-1981.nc"))
x = maximum(abs, x) > 1e4 ? x ./ 1e3 : x; y = maximum(abs, y) > 1e4 ? y ./ 1e3 : y
nx, ny = 106, 181
I = [round(Int, (v + 720) / 16) + 1 for v in x]; J = [round(Int, (v + 3450) / 16) + 1 for v in y]
ice = coalesce.(msk, 0) .>= 2
inb(i, j) = 1 <= I[i] <= nx && 1 <= J[j] <= ny
acc = Dict(n => zeros(nx, ny) for (_, n) in vars); cnt = zeros(nx, ny); cg = zeros(nx, ny); zs = zeros(nx, ny)
for j in eachindex(y), i in eachindex(x)
    (ice[i, j] && inb(i, j)) || continue
    cnt[I[i], J[j]] += 1; zs[I[i], J[j]] += srf[i, j]; coalesce(msk[i, j], 0) == 4 && (cg[I[i], J[j]] += 1)
end
cls = Dict(k => Float64[] for k in 0:4); first_year = true
for yr in 1981:2010
    NCDataset(joinpath(D, "MARv3.14-monthly-ERA5-$yr.nc")) do ds
        for (v, n) in vars
            a = dropdims(sum(Float64.(coalesce.(ds[v][:, :, 1:12], NaN)); dims=3); dims=3)    # [mm w.e./yr]
            A = acc[n]
            for j in eachindex(y), i in eachindex(x)
                (ice[i, j] && inb(i, j)) || continue
                A[I[i], J[j]] += a[i, j]
            end
            if yr == 1981
                for k in 0:4
                    b = coalesce.(msk, -1) .== k
                    any(b) && @printf("1981 %-3s MSK %d: n %8d  mean %8.1f  SRF %6.0f\n", v, k, count(b), mean(filter(isfinite, a[b])), mean(filter(isfinite, srf[b])))
                end
            end
        end
    end
end
NCDataset(out, "c") do ds
    defDim(ds, "xc", nx); defDim(ds, "yc", ny)
    defVar(ds, "xc", collect(-720.0:16:-720+16(nx-1)), ("xc",); attrib=["units" => "km"])
    defVar(ds, "yc", collect(-3450.0:16:-3450+16(ny-1)), ("yc",); attrib=["units" => "km"])
    for (_, n) in vars
        defVar(ds, n, Float32.(ifelse.(cnt .> 0, acc[n] ./ (30 .* max.(cnt, 1)), NaN)), ("xc", "yc");
               attrib=["units" => "kg m-2 yr-1", "long_name" => "$n, 1981-2010 mean over the ice of the cell"])
    end
    defVar(ds, "z_srf", Float32.(ifelse.(cnt .> 0, zs ./ max.(cnt, 1), NaN)), ("xc", "yc"); attrib=["units" => "m"])
    defVar(ds, "f_ice", Float32.(cnt ./ 256), ("xc", "yc"); attrib=["units" => "1", "long_name" => "ice share of the cell (MAR MSK >= 2)"])
    defVar(ds, "f_gris", Float32.(cg ./ 256), ("xc", "yc"); attrib=["units" => "1", "long_name" => "ice-sheet share of the cell (MAR MSK = 4)"])
    ds.attrib["source"] = "MARv3.14 ERA5 1 km monthly (Fettweis), 1981-2010; ice-only means binned to GRL-16KM (kryos tools/levante/mar_ice16.jl)"
    ds.attrib["grid_name"] = "GRL-16KM"
end
