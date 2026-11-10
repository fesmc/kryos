#!/usr/bin/env julia
# setup-diag (docs/dev/setup-diagnosis.md): basal and flow state of existing calib-cf final states, analysis only (no runs).
#   julia +1.13 --project=$G/julia-plot setup_diag.jl <out dir>      (sbatch setup_diag.sh)
# Reads the Yelmo restart of each state (all diagnostic fields: T_prime_b, W_til, N_eff, c_bed, f_slide, taud, taub, u, Q_geo, smb, ...)
# and FesmData observations (velocity, SMB, GHF; ANT basal state). Greenland obs (GRL-32KM) are mapped to NH-32KM by nearest neighbour.
# Writes <out>/setup_diag.md (tables) and <out>/setup_diag_<case>.png (maps).
using NCDatasets, Statistics, Printf, CairoMakie

const S9 = "/work/ba1442/robinson/kryos-gate/runs/calib/kx-624f686/S9"
const FD = "/work/ba1442/fesmdata"
const GRLSA = "/work/ba1442/robinson/kryos-gate/runs/grl/kx-af8f1b9/S3/base/restart-15000/Greenland/ice_grl_restart.nc"
const RI, RSW, GRAV = 917.0, 1028.0, 9.81

const VARS = ["H_ice", "H_ice_ref", "f_grnd", "z_bed", "z_sl", "T_prime_b", "hyd_W_til", "hyd_W_til_max", "N_eff", "cb_ref",
              "c_bed", "f_slide", "beta", "taub", "taud", "uxy_s", "uxy_b", "bmb_grnd", "Q_geo", "Q_b", "smb", "T_srf", "enh_bar",
              "ATT_bar", "dHidt", "x2D", "y2D", "lon2D", "lat2D", "area"]

function rdnc(f, vars)
    ds = NCDataset(f); d = Dict{String,Matrix{Float64}}()
    for v in vars
        haskey(ds, v) || (println("  missing $v in $f"); continue)
        x = ds[v]; a = ndims(x) == 3 ? x[:, :, size(x, 3)] : x[:, :]
        d[v] = Float64.(coalesce.(a, NaN))
    end
    close(ds); d
end
function rdvar(f, v)
    ds = NCDataset(f); x = ds[v]; a = ndims(x) == 3 ? x[:, :, size(x, 3)] : x[:, :]
    u = get(x.attrib, "units", ""); a = Float64.(coalesce.(a, NaN)); close(ds); (a, u)
end
rdint(f, v) = (ds = NCDataset(f); r = Int.(coalesce.(ds[v][:, :], -1)); close(ds); r)

function state(f)
    d = rdnc(f, VARS); sz = size(d["H_ice"])
    for v in VARS; haskey(d, v) || (d[v] = fill(NaN, sz)); end
    maximum(abs.(filter(isfinite, d["x2D"]))) > 1e5 && (d["x2D"] ./= 1e3; d["y2D"] ./= 1e3)    # km
    maximum(filter(isfinite, d["area"])) < 1e6 && (d["area"] .*= 1e6)                          # m2
    ds = NCDataset(f); d["xc1"] = reshape(Float64.(ds["xc"][:]), :, 1); d["yc1"] = reshape(Float64.(ds["yc"][:]), :, 1); close(ds)
    maximum(abs.(d["xc1"])) > 1e5 && (d["xc1"] ./= 1e3; d["yc1"] ./= 1e3)
    d
end

# distance [km] of the cells of `ice` to the nearest non-ice cell touching ice
function dist_edge(ice::BitMatrix, x, y)
    nx, ny = size(ice); bx = Float64[]; by = Float64[]
    for j in 1:ny, i in 1:nx
        ice[i, j] && continue
        if any(1 <= i + di <= nx && 1 <= j + dj <= ny && ice[i+di, j+dj] for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1)))
            push!(bx, x[i, j]); push!(by, y[i, j])
        end
    end
    d = fill(NaN, nx, ny)
    for j in 1:ny, i in 1:nx
        ice[i, j] && (d[i, j] = minimum(hypot.(bx .- x[i, j], by .- y[i, j])))
    end
    d
end

# nearest-neighbour index map source -> target (unit vectors), within dmax km
function nnindex(lons, lats, lont, latt, mt; dmax = 25.0)
    xs = vec(cosd.(lats) .* cosd.(lons)); ys = vec(cosd.(lats) .* sind.(lons)); zs = vec(sind.(lats))
    idx = zeros(Int, size(lont))
    for I in findall(mt)
        xt, yt, zt = cosd(latt[I]) * cosd(lont[I]), cosd(latt[I]) * sind(lont[I]), sind(latt[I])
        d2 = (xs .- xt) .^ 2 .+ (ys .- yt) .^ 2 .+ (zs .- zt) .^ 2; k = argmin(d2)
        sqrt(d2[k]) * 6371 <= dmax && (idx[I] = k)
    end
    idx
end
applyidx(idx, f) = [k == 0 ? NaN : f[k] for k in idx]

smb_mi(a, u) = occursin("kg", u) || occursin("mm", u) ? a ./ RI : a          # -> m ice / yr

fin(x, m) = (k = m .& isfinite.(x); x[k])
wmean(x, w, m) = (k = m .& isfinite.(x); any(k) ? sum(x[k] .* w[k]) / sum(w[k]) : NaN)
med(x, m) = (v = fin(x, m); isempty(v) ? NaN : median(v))
q(x, m, p) = (v = fin(x, m); isempty(v) ? NaN : quantile(v, p))
afrac(c, w, m) = (s = sum(w[m]); s > 0 ? sum(w[m .& c]) / s : NaN)

# the rows of the per-state table; s: state, obs: Dict (uobs, smbobs, tsobs on the model grid), msk: masks (name => BitMatrix),
# till_rate [m/yr]
function metrics(s, obs, masks, till_rate)
    H, Hr, w = s["H_ice"], s["H_ice_ref"], s["area"]
    gm = (H .> 10) .& (s["f_grnd"] .>= 0.99)
    P0 = RI * GRAV .* H; nP = s["N_eff"] ./ P0
    temp = s["T_prime_b"] .> -0.5
    melt = -s["bmb_grnd"]                       # Yelmo bmb < 0: melt [m ice/yr]
    sat = s["hyd_W_til"] ./ s["hyd_W_til_max"]
    rcd = s["c_bed"] ./ s["taud"]
    cfn = s["taud"] ./ s["N_eff"]
    fb = s["uxy_b"] ./ s["uxy_s"]
    uo = obs["uobs"]; ok_u = isfinite.(uo) .& (uo .> 1) .& (s["uxy_s"] .> 0)
    lr = log10.(max.(s["uxy_s"], 1e-3) ./ max.(uo, 1e-3))
    rows = Vector{Pair{String,Vector{String}}}()
    add(n, f) = push!(rows, n => [f(m) for (_, m) in masks])
    add("grounded ice area, model / ref [1e3 km2]", m -> @sprintf("%.0f / %.0f", sum(w[m .& gm]) / 1e9, sum(w[m .& (Hr .> 10)]) / 1e9))
    add("H bias mean / rms on ref ice [m]", m -> (k = m .& (Hr .> 10); @sprintf("%+.0f / %.0f", mean(H[k] .- Hr[k]), sqrt(mean((H[k] .- Hr[k]) .^ 2)))))
    add("H mean model / ref on ref ice [m]", m -> (k = m .& (Hr .> 10); @sprintf("%.0f / %.0f", mean(H[k]), mean(Hr[k]))))
    add("dH/dt rms [m/yr]", m -> @sprintf("%.3f", sqrt(mean(fin(s["dHidt"], m .& gm) .^ 2))))
    add("temperate bed (T'_b > -0.5 K), area fraction", m -> @sprintf("%.2f", afrac(temp, w, m .& gm)))
    add("wet till (W_til > 1 cm), area fraction", m -> @sprintf("%.2f", afrac(s["hyd_W_til"] .> 0.01, w, m .& gm)))
    add("till saturation W_til/W_til_max, temperate: mean", m -> @sprintf("%.2f", wmean(sat, w, m .& gm .& temp)))
    add("basal melt, temperate: median / p90 [mm/yr]", m -> @sprintf("%.1f / %.1f", 1e3 * med(melt, m .& gm .& temp), 1e3 * q(melt, m .& gm .& temp, 0.9)))
    add("temperate cells with melt > till_rate ($(1e3*till_rate) mm/yr)", m -> @sprintf("%.2f", afrac(melt .> till_rate, w, m .& gm .& temp)))
    add("N_eff/P0 median: all / temperate / cold", m -> @sprintf("%.2f / %.2f / %.2f", med(nP, m .& gm), med(nP, m .& gm .& temp), med(nP, m .& gm .& .!temp)))
    add("N_eff/P0 > 0.9, area fraction", m -> @sprintf("%.2f", afrac(nP .> 0.9, w, m .& gm)))
    add("cb_ref median", m -> @sprintf("%.3f", med(s["cb_ref"], m .& gm)))
    add("c_bed median / p10 [kPa]", m -> @sprintf("%.0f / %.0f", med(s["c_bed"], m .& gm) / 1e3, q(s["c_bed"], m .& gm, 0.1) / 1e3))
    add("taud median [kPa]", m -> @sprintf("%.0f", med(s["taud"], m .& gm) / 1e3))
    add("taub median [kPa]", m -> @sprintf("%.0f", med(s["taub"], m .& gm) / 1e3))
    add("c_bed/taud median", m -> @sprintf("%.1f", med(rcd, m .& gm)))
    add("c_bed < 2 taud (friction-limited), area fraction", m -> @sprintf("%.2f", afrac(rcd .< 2, w, m .& gm)))
    add("cf for c_bed = taud (taud/N_eff): median / temperate", m -> @sprintf("%.4f / %.4f", med(cfn, m .& gm), med(cfn, m .& gm .& temp)))
    add("f_slide median; < 0.1 area fraction", m -> @sprintf("%.2f; %.2f", med(s["f_slide"], m .& gm), afrac(s["f_slide"] .< 0.1, w, m .& gm)))
    add("beta median [Pa yr/m]", m -> @sprintf("%.2g", med(s["beta"], m .& gm)))
    add("u_b/u_s median (u_s > 1 m/yr)", m -> @sprintf("%.2f", med(fb, m .& gm .& (s["uxy_s"] .> 1))))
    add("sliding-dominated (u_b/u_s > 0.5), area fraction", m -> @sprintf("%.2f", afrac(fb .> 0.5, w, m .& gm)))
    add("u_s median model / obs, obs cells [m/yr]", m -> @sprintf("%.1f / %.1f", med(s["uxy_s"], m .& gm .& ok_u), med(uo, m .& gm .& ok_u)))
    add("u_s p90 model / obs [m/yr]", m -> @sprintf("%.0f / %.0f", q(s["uxy_s"], m .& gm .& ok_u, 0.9), q(uo, m .& gm .& ok_u, 0.9)))
    add("log10(u_s model / obs) median", m -> @sprintf("%+.2f", med(lr, m .& gm .& ok_u)))
    add("Q_geo mean [mW/m2]", m -> @sprintf("%.1f", wmean(s["Q_geo"], w, m .& gm)))
    add("Q_b (friction heat) mean [mW/m2]", m -> @sprintf("%.1f", wmean(s["Q_b"], w, m .& gm)))
    add("T_srf mean model / obs on ref ice [degC]", m -> (k = m .& (Hr .> 10); @sprintf("%.1f / %.1f", wmean(s["T_srf"], w, k) - 273.15, wmean(obs["tsobs"], w, k) - 273.15)))
    add("SMB mean model / obs on ref ice [m ice/yr]", m -> (k = m .& (Hr .> 10); @sprintf("%.3f / %.3f", wmean(s["smb"], w, k), wmean(obs["smbobs"], w, k))))
    add("SMB integral model / obs on ref ice [Gt/yr]", m -> (k = m .& (Hr .> 10); kk = k .& isfinite.(s["smb"]) .& isfinite.(obs["smbobs"]);
        @sprintf("%.0f / %.0f", sum(s["smb"][kk] .* w[kk]) * RI / 1e12, sum(obs["smbobs"][kk] .* w[kk]) * RI / 1e12)))
    add("ablation area (smb < 0) model / obs, obs cells, area fraction", m -> (k = m .& (Hr .> 10) .& isfinite.(obs["smbobs"]); @sprintf("%.2f / %.2f", afrac(s["smb"] .< 0, w, k), afrac(obs["smbobs"] .< 0, w, k))))
    add("SMB mean model / obs, obs cells [m ice/yr]", m -> (k = m .& (Hr .> 10) .& isfinite.(obs["smbobs"]) .& isfinite.(s["smb"]); @sprintf("%.3f / %.3f", wmean(s["smb"], w, k), wmean(obs["smbobs"], w, k))))
    add("enh_bar mean", m -> @sprintf("%.2f", wmean(s["enh_bar"], w, m .& gm)))
    add("ATT_bar mean [1e-17 Pa-3 yr-1]", m -> @sprintf("%.2f", wmean(s["ATT_bar"], w, m .& gm) * 1e17))
    rows
end

function table(io, title, masks, rows)
    println(io, "### $title\n")
    println(io, "| metric | " * join(first.(masks), " | ") * " |"); println(io, "|---|" * repeat("---|", length(masks)))
    for (n, v) in rows; println(io, "| $n | " * join(v, " | ") * " |"); end
    println(io)
end

function sens(io, title, a, b, masks)          # a: higher cf_ref, b: lower
    println(io, "### $title\n")
    println(io, "| mask | dH mean, b - a [m] | du_s median ratio b/a | du_b median ratio b/a | dV ice [1e3 km3] |"); println(io, "|---|---|---|---|---|")
    for (n, m) in masks
        g = m .& (a["H_ice"] .> 10) .& (b["H_ice"] .> 10)
        rs = b["uxy_s"] ./ a["uxy_s"]; rb = b["uxy_b"] ./ max.(a["uxy_b"], 1e-6)
        println(io, @sprintf("| %s | %+.0f | %.2f | %.2f | %+.1f |", n, mean(b["H_ice"][g] .- a["H_ice"][g]), med(rs, g), med(rb, g),
            sum((b["H_ice"] .- a["H_ice"])[m] .* a["area"][m]) / 1e12))
    end
    println(io)
end

function maps(png, s, obs, box, title)
    i1, i2, j1, j2 = box
    x = vec(s["xc1"])[i1:i2]; y = vec(s["yc1"])[j1:j2]
    sub(a) = a[i1:i2, j1:j2]
    gm = (s["H_ice"] .> 10) .& (s["f_grnd"] .>= 0.99); msk(a) = sub(ifelse.(gm, a, NaN))
    ref = s["H_ice_ref"] .> 10
    P0 = RI * GRAV .* s["H_ice"]
    pan = [
        ("H - H_ref [m]", sub(ifelse.(ref .| (s["H_ice"] .> 10), s["H_ice"] .- s["H_ice_ref"], NaN)), (-800, 800), :RdBu, true),
        ("T'_b [K]", msk(s["T_prime_b"]), (-15, 0), :viridis, false),
        ("W_til/W_til_max", msk(s["hyd_W_til"] ./ s["hyd_W_til_max"]), (0, 1), :Blues, false),
        ("log10 N_eff/P0", msk(log10.(max.(s["N_eff"] ./ P0, 1e-4))), (-2, 0), :viridis, false),
        ("log10 c_bed/taud", msk(log10.(max.(s["c_bed"] ./ s["taud"], 1e-3))), (-1, 2), :RdBu, true),
        ("u_b/u_s", msk(s["uxy_b"] ./ max.(s["uxy_s"], 1e-6)), (0, 1), :magma, false),
        ("log10 u_s model/obs", msk(log10.(max.(s["uxy_s"], 1e-2) ./ max.(obs["uobs"], 1e-2))), (-1.5, 1.5), :RdBu, true),
        ("SMB model - obs [m ice/yr]", sub(ifelse.(ref, s["smb"] .- obs["smbobs"], NaN)), (-1, 1), :RdBu, true),
        ("Q_geo [mW/m2]", msk(s["Q_geo"]), (30, 90), :inferno, false),
    ]
    fig = Figure(size = (1500, 1450)); Label(fig[0, 1:3], title, fontsize = 18)
    for (k, (t, a, cr, cm, rev)) in enumerate(pan)
        r, c = divrem(k - 1, 3)
        ax = Axis(fig[r+1, c+1][1, 1], title = t, aspect = DataAspect()); hidedecorations!(ax)
        hm = heatmap!(ax, x, y, a, colorrange = cr, colormap = rev ? Reverse(cm) : cm, nan_color = :grey85)
        contour!(ax, x, y, Float64.(sub(ref)), levels = [0.5], color = :black, linewidth = 0.7)
        Colorbar(fig[r+1, c+1][1, 2], hm)
    end
    save(png, fig)
end

function bbox(m; pad = 2)
    I = findall(m); is = [c[1] for c in I]; js = [c[2] for c in I]
    (max(minimum(is) - pad, 1), min(maximum(is) + pad, size(m, 1)), max(minimum(js) - pad, 1), min(maximum(js) + pad, size(m, 2)))
end

function main(out)
    mkpath(out); io = IOBuffer()
    println(io, "# setup-diag: field diagnosis of the calib-cf final states\n")
    println(io, "Generated by `tools/levante/setup_diag.jl` from the Yelmo restarts (final state). Grounded ice: H > 10 m, f_grnd >= 0.99. ",
        "Temperate: T'_b > -0.5 K. Area-weighted fractions. Reference thickness: H_ice_ref of the run (the hub: GEBCO 2026 North, BedMachine ANT; ",
        "BedMachine v6 standalone Greenland). Interior / margin: distance to the edge of the reference ice > 300 km / < 100 km.\n")

    # ---------------- North PD: Greenland ----------------
    fn = "$FD/North/NH-32KM/NH-32KM"
    r2 = rdint("$(fn)_REGIONS.nc", "region_2")
    nor10 = state("$S9/North-pd-cf10/restart-30000/North/ice_nor_restart.nc")
    nor20 = state("$S9/North-pd/restart-30000/North/ice_nor_restart.nc")
    grl = r2 .== 103
    gref = grl .& (nor10["H_ice_ref"] .> 10)
    dg = dist_edge(gref, nor10["x2D"], nor10["y2D"])
    # Greenland obs (GRL-32KM) -> NH-32KM
    fg = "$FD/Greenland/GRL-32KM/GRL-32KM"
    lon_g, _ = rdvar("$(fg)_grid.nc", "lon2D"); lat_g, _ = rdvar("$(fg)_grid.nc", "lat2D")
    idx = nnindex(lon_g, lat_g, nor10["lon2D"], nor10["lat2D"], grl)
    uo, _ = rdvar("$(fg)_Velocity-Joughin2018.nc", "uxy_srf"); fv, _ = rdvar("$(fg)_Velocity-Joughin2018.nc", "f_valid_uxy_srf"); uo[.!(fv .> 0.5)] .= NaN
    sm, su = rdvar("$(fg)_SMB-MARv3.14-1981-2010-annual-mean.nc", "smb"); ts, _ = rdvar("$(fg)_SMB-MARv3.14-1981-2010-annual-mean.nc", "T_srf")
    mk, _ = rdvar("$(fg)_SMB-MARv3.14-1981-2010-annual-mean.nc", "mask"); sm[mk .!= 4] .= NaN; ts[mk .!= 4] .= NaN     # MAR ice sheet only
    gA, _ = rdvar("$(fg)_grid.nc", "area"); kk = isfinite.(sm); println(io, @sprintf("MAR v3.14 1981-2010 SMB on its GRL-32KM ice-sheet mask (native): %.0f Gt/yr over %.0f 1e3 km2.\n", sum(sm[kk] .* gA[kk]) / 1e12, sum(gA[kk]) / 1e9))
    obsn = Dict("uobs" => applyidx(idx, uo), "smbobs" => applyidx(idx, smb_mi(sm, su)), "tsobs" => applyidx(idx, ts))
    mg = ["Greenland (103)" => grl, "interior > 300 km" => grl .& (dg .> 300), "100-300 km" => grl .& (dg .>= 100) .& (dg .<= 300), "margin < 100 km" => grl .& (dg .< 100)]
    println(io, "## North PD, Greenland (region 103): North-pd-cf10, cf_ref 0.1, t = 30 kyr\n")
    table(io, "North-pd-cf10 (TILL closure, W_til_max 0.4, till_rate 8 mm/yr; obs: Joughin 2018, MAR v3.14 1981-2010)", mg, metrics(nor10, obsn, mg, 0.008))
    table(io, "North-pd (cf_ref 0.2)", mg, metrics(nor20, obsn, mg, 0.008))
    sens(io, "cf_ref 0.2 (a) -> 0.1 (b), Greenland", nor20, nor10, mg)
    # GHF datasets over Greenland (ref ice)
    println(io, "### GHF over the Greenland reference ice [mW/m2]: mean (interior)\n")
    println(io, "| dataset | mean | interior > 300 km |"); println(io, "|---|---|---|")
    w = nor10["area"]
    println(io, @sprintf("| Shapiro2004 (model Q_geo) | %.1f | %.1f |", wmean(nor10["Q_geo"], w, gref), wmean(nor10["Q_geo"], w, gref .& (dg .> 300))))
    for ds in ("Colgan2021", "Martos2018", "FoxMaule2005", "Lucazeau2019", "Shapiro2004", "Davies2013")
        f = "$(fg)_GHF-$ds.nc"; isfile(f) || continue
        a, _ = rdvar(f, "ghf"); b = applyidx(idx, a)
        println(io, @sprintf("| %s (GRL-32KM) | %.1f | %.1f |", ds, wmean(b, w, gref), wmean(b, w, gref .& (dg .> 300))))
    end
    println(io)
    maps("$out/setup_diag_nor_pd_cf10.png", nor10, obsn, bbox(gref; pad = 3), "North-pd-cf10, Greenland (cf_ref 0.1, t = 30 kyr)")

    # ---------------- standalone Greenland (yelmox_Greenland settings, GRL-16KM) ----------------
    if isfile(GRLSA)
        f16 = "$FD/Greenland/GRL-16KM/GRL-16KM"
        sa = state(GRLSA)
        g16 = rdint("$(f16)_REGIONS.nc", "region_2") .== 103
        u16, _ = rdvar("$(f16)_Velocity-Joughin2018.nc", "uxy_srf"); v16, _ = rdvar("$(f16)_Velocity-Joughin2018.nc", "f_valid_uxy_srf"); u16[.!(v16 .> 0.5)] .= NaN
        s16, su16 = rdvar("$(f16)_SMB-MARv3.14-1981-2010-annual-mean.nc", "smb"); t16, _ = rdvar("$(f16)_SMB-MARv3.14-1981-2010-annual-mean.nc", "T_srf")
        mk16, _ = rdvar("$(f16)_SMB-MARv3.14-1981-2010-annual-mean.nc", "mask"); s16[mk16 .!= 4] .= NaN; t16[mk16 .!= 4] .= NaN
        obs16 = Dict("uobs" => u16, "smbobs" => smb_mi(s16, su16), "tsobs" => t16)
        gr16 = g16 .& (sa["H_ice_ref"] .> 10); d16 = dist_edge(gr16, sa["x2D"], sa["y2D"])
        m16 = ["Greenland (103)" => g16, "interior > 300 km" => g16 .& (d16 .> 300), "100-300 km" => g16 .& (d16 .>= 100) .& (d16 .<= 300), "margin < 100 km" => g16 .& (d16 .< 100)]
        println(io, "## Standalone Greenland, yelmox_Greenland settings (GRL-16KM, cf_ref 0.8, W_til_max 2, till_rate 1 mm/yr, GHF 50, MAR+ITM), gate run S3 base, t = 15 kyr\n")
        table(io, "grl S3 base (obs: Joughin 2018, MAR v3.14 1981-2010 on GRL-16KM)", m16, metrics(sa, obs16, m16, 0.001))
        maps("$out/setup_diag_grl_standalone.png", sa, obs16, bbox(gr16; pad = 3), "standalone Greenland (yelmox settings, GRL-16KM, t = 15 kyr)")
    end

    # ---------------- North lgmrecon (frozen 21 ka) ----------------
    lr10 = state("$S9/North-lgmrecon-cf10/restart-20000/North/ice_nor_restart.nc")
    lr20f = "$S9/North-lgmrecon/restart-20000/North/ice_nor_restart.nc"
    nanobs = Dict(k => fill(NaN, size(lr10["H_ice"])) for k in ("uobs", "smbobs", "tsobs"))
    ml = ["North America (101)" => r2 .== 101, "W Eurasia (102)" => r2 .== 102, "Greenland (103)" => grl]
    println(io, "## North lgmrecon (frozen 21 ka forcing, from PaleoMIST 20 ka): North-lgmrecon-cf10, t = 20 kyr\n")
    println(io, "H bias here is against the PD reference (not meaningful); no obs.\n")
    table(io, "North-lgmrecon-cf10", ml, metrics(lr10, nanobs, ml, 0.008))
    isfile(lr20f) && sens(io, "lgmrecon cf_ref 0.2 (a) -> 0.1 (b)", state(lr20f), lr10, ml)
    maps("$out/setup_diag_nor_lgmrecon_cf10.png", lr10, nanobs, (1, size(r2, 1), 1, size(r2, 2)), "North-lgmrecon-cf10 (t = 20 kyr; H bias vs PD)")

    # ---------------- Antarctica PD ----------------
    fa = "$FD/Antarctica/ANT-32KM/ANT-32KM"
    r3 = rdint("$(fa)_REGIONS.nc", "region_3")
    a05 = state("$S9/Antarctica-pd-cf05/restart-30000/Antarctica/ice_ant_restart.nc")
    a10 = state("$S9/Antarctica-pd-cf10/restart-30000/Antarctica/ice_ant_restart.nc")
    ua, _ = rdvar("$(fa)_Velocity-Rignot2011.nc", "uxy_srf"); va, _ = rdvar("$(fa)_Velocity-Rignot2011.nc", "f_valid_uxy_srf"); ua[.!(va .> 0.5)] .= NaN
    sa_, sau = rdvar("$(fa)_SMB-RACMO2.3p2-1981-2010-annual-mean.nc", "smb"); ta, _ = rdvar("$(fa)_SMB-RACMO2.3p2-1981-2010-annual-mean.nc", "T_srf")
    obsa = Dict("uobs" => ua, "smbobs" => smb_mi(sa_, sau), "tsobs" => ta)
    aref = (a05["H_ice_ref"] .> 10) .& (a05["H_ice_ref"] .* RI ./ RSW .> -a05["z_bed"])      # grounded reference
    da = dist_edge(aref, a05["x2D"], a05["y2D"])
    wa, ea = r3 .== 20101, r3 .== 20102
    land = a05["z_bed"] .>= 0
    ma = ["WAIS" => wa, "WAIS interior > 300 km" => wa .& (da .> 300), "EAIS" => ea, "EAIS interior > 300 km" => ea .& (da .> 300),
          "EAIS bed >= 0" => ea .& land, "EAIS bed < 0" => ea .& .!land, "EAIS margin < 100 km" => ea .& (da .< 100)]
    println(io, "## Antarctica PD: Antarctica-pd-cf05 (cf_ref 0.05, scale_zb exp -400/0, MARINE closure p = 1), t = 30 kyr\n")
    table(io, "Antarctica-pd-cf05 (obs: Rignot 2011, RACMO2.3p2 1981-2010)", ma, metrics(a05, obsa, ma, 0.001))
    sens(io, "cf_ref 0.1 (a) -> 0.05 (b), Antarctica", a10, a05, ma)
    # Dawson 2022 basal state and GHF datasets
    fdw = "$(fa)_BasalState-Dawson2022.nc"
    if isfile(fdw)
        bs, _ = rdvar(fdw, "basal_state"); bd, _ = rdvar(fdw, "beta"); bm, _ = rdvar(fdw, "bmelt_grnd")
        gm = (a05["H_ice"] .> 10) .& (a05["f_grnd"] .>= 0.99)
        println(io, "### Antarctica: Dawson et al. (2022) inversion vs model (cf05)\n")
        println(io, "| mask | thawed / thawable / frozen (Dawson) | temperate (model) | beta median Dawson (linear law) / model [Pa yr/m] | bmelt median Dawson / model [mm/yr] |"); println(io, "|---|---|---|---|---|")
        for (n, m) in ma
            k = m .& aref .& isfinite.(bs)
            println(io, @sprintf("| %s | %.2f / %.2f / %.2f | %.2f | %.2g / %.2g | %.1f / %.1f |", n, afrac(bs .== 3, a05["area"], k), afrac(bs .== 2, a05["area"], k), afrac(bs .== 1, a05["area"], k),
                afrac(a05["T_prime_b"] .> -0.5, a05["area"], k .& gm), med(bd, k), med(a05["beta"], k .& gm), 1e3 * med(bm, k), -1e3 * med(a05["bmb_grnd"], k .& gm)))
        end
        println(io)
    end
    println(io, "### GHF over the Antarctic grounded reference ice [mW/m2]\n")
    println(io, "| dataset | WAIS | EAIS | EAIS interior > 300 km |"); println(io, "|---|---|---|---|")
    w = a05["area"]
    println(io, @sprintf("| model Q_geo (Shapiro2004) | %.1f | %.1f | %.1f |", wmean(a05["Q_geo"], w, aref .& wa), wmean(a05["Q_geo"], w, aref .& ea), wmean(a05["Q_geo"], w, aref .& ea .& (da .> 300))))
    for ds in ("Martos2017", "Stal2020", "Loesing2021", "HazzardRichards2024", "Lucazeau2019", "FoxMaule2005", "Shapiro2004")
        f = "$(fa)_GHF-$ds.nc"; isfile(f) || continue
        ds_ = NCDataset(f); cand = [k for k in keys(ds_) if occursin("ghf", lowercase(k)) && !occursin("valid", k) && !occursin("_sd", k)]; close(ds_)
        vn = "ghf" in cand ? "ghf" : (isempty(cand) ? "none" : cand[1])
        vn == "none" && continue
        a, _ = rdvar(f, vn)
        println(io, @sprintf("| %s | %.1f | %.1f | %.1f |", ds, wmean(a, w, aref .& wa), wmean(a, w, aref .& ea), wmean(a, w, aref .& ea .& (da .> 300))))
    end
    println(io)
    maps("$out/setup_diag_ant_pd_cf05.png", a05, obsa, (1, size(r3, 1), 1, size(r3, 2)), "Antarctica-pd-cf05 (t = 30 kyr)")

    open("$out/setup_diag.md", "w") do f; write(f, String(take!(io))); end
    println("wrote $out/setup_diag.md")
    return 0
end
exit(main(ARGS[1]))
