# G2b diagnosis: where do the load and ocean function differ between kryos (host-resolved flotation) and the standalone vilma.x,
# and does the load difference explain the dz_bed difference.
#   julia --project=tools tools/levante/g2b_diag.jl <run base: .../kx-<sha>/S9> <tag> <maps folder> <out folder> [--ref <tag>]
# Reads <base>/<tag>/{kryos,vilma}. Writes <out>/g2b_diag_<tag>.nc (Gauss grid, selected times: dz_bed, C, the load and its parts
# for both runs, host_w) and <out>/g2b_diag_<tag>.md (series of the load difference by region, the local isostatic estimate of the
# dz_bed difference it implies, point series). --ref <tag>: also the dz_bed difference of <tag> (e.g. a run pair without spin-up).
#
# The load, in metres of sea water (VILMA2 vilma_sle, host mode): L = (rho_i/rho_w) dI_g + C rsl - fill, with
#   standalone: dI_g = h (1 - C) - h_eq (1 - C0), fill = z_bed_eq (C - C0), C0 = ocean(z_bed_eq, h_eq)   (h = all ice, VILMA2's flotation)
#   kryos:      dI_g = (host_Ig - host_Ig_ref) + a h_o (1 - C_own) - a0 h_eq (1 - C0), fill = host_F + a z_bed_eq (C_own - C0),
#               a = 1 - w, C_own = (C - host_C) / a; host_Ig, host_C, host_F recomputed here from the ice-grid records
#               (H_ice f_grnd, 1 - f_grnd, z_bed_ref (f_grnd,ref - f_grnd), con maps ice -> Gauss, per unit Gauss area),
#               f_grnd,ref the binary rule on H_ice(0 ka), z_bed_ref, z_sl = 0; the own part's ice and reference from the
#               standalone's fields (the same PaleoMIST file, VILMA2's remap instead of kryos').
# The record at t holds the ice of t (published after the scatter of step t); the Earth at t solved with the ice of t - 10 yr: the
# 10-yr offset is ignored here.
using NCDatasets, Printf, Statistics, SparseArrays

const RI = 910.0; const RW = 1028.0; const RM = 3300.0
ocean(topo, ice) = Float64.((topo .< 0) .& (RI .* ice .< -RW .* topo))

function conmap(f)
    ds = NCDataset(f)
    s = Int.(ds["src_address"][:]); d = Int.(ds["dst_address"][:]); wt = Float64.(ds["remap_matrix"][1, :])
    ns = ds.dim["src_grid_size"]; nd = ds.dim["dst_grid_size"]; close(ds)
    sparse(d, s, wt, nd, ns)                         # weights are overlap areas [m2]
end

wsum(x, a, m) = sum(x[m] .* a[m])
wrms(x, a, m) = sqrt(wsum(x .^ 2, a, m) / wsum(ones(size(x)), a, m))

function main(args)
    ref = nothing
    i = findfirst(==("--ref"), args); i === nothing || (ref = args[i+1]; deleteat!(args, i:i+1))
    base, tag, maps, out = args[1:4]; mkpath(out)
    kd = joinpath(base, tag, "kryos"); vd = joinpath(base, tag, "vilma")
    k = NCDataset(joinpath(kd, "kryos_out.nc")); v = NCDataset(joinpath(vd, "vilma_out.nc"))
    lon = Float64.(k["lon"][:]); lat = Float64.(k["lat"][:]); nphi = length(lon); nlat = length(lat); lmax = nphi ÷ 4
    flip = Float64(v["lat"][1]) > Float64(v["lat"][end])
    fix(a) = flip ? a[:, end:-1:1] : a
    rs(name, n) = fix(Float64.(nomissing(v[name][:, :, n], 0.0)))
    rk(name, r) = Float64.(nomissing(k[name][:, :, r], 0.0))
    tr = Float64.(k["rec"][:]); nrec = length(tr)
    stepof(r) = round(Int, tr[r] - tr[1]) ÷ 10
    w = Float64.(k["host_w"][:, :]); area = Float64.(k["area"][:, :]); LAT = repeat(reshape(lat, 1, :), nphi, 1)
    zeq = fix(Float64.(v["z_bed_eq"][:, :])); heq = fix(Float64.(v["h_ice_eq"][:, :])); C0 = ocean(zeq, heq)
    a = max.(0.0, 1.0 .- w)

    # the domains: ice-grid records and the con maps to the Gauss grid
    doms = []
    for (dn, gn) in (("North", "NH-32KM"), ("Antarctica", "ANT-32KM"))
        M = conmap(joinpath(maps, "scrip-con_$(gn)_GAUSS-L$(lmax).nc"))
        ds = NCDataset(joinpath(kd, "g2b_$(dn).nc"))
        zbr = Float64.(ds["z_bed_ref"][:, :]); Hpd = Float64.(ds["H_ice"][:, :, nrec])
        fgr = Float64.((zbr .>= 0) .| (RI .* Hpd .> -RW .* zbr))
        push!(doms, (name = dn, M = M, ds = ds, zbr = zbr, fgr = fgr))
    end
    tog(M, x) = reshape(M * vec(x), nphi, nlat) ./ area
    Igref = sum(tog(d.M, d.ds["H_ice"][:, :, nrec] .* d.fgr) for d in doms)

    function loads(r)
        n = stepof(r) + 1
        Cs = rs("C_ocean", n); rsl_s = rs("rsl", n); h = rs("h_ice", n)
        dIs = h .* (1 .- Cs) .- heq .* (1 .- C0); Ws = Cs .* rsl_s .- zeq .* (Cs .- C0)
        Ck = rk("C_ocn", r); rsl_k = rk("rsl", r)
        hIg = zeros(nphi, nlat); hC = zeros(nphi, nlat); hF = zeros(nphi, nlat)
        for d in doms
            H = Float64.(d.ds["H_ice"][:, :, r]); fg = Float64.(d.ds["f_grnd"][:, :, r])
            hIg .+= tog(d.M, H .* fg); hC .+= tog(d.M, 1 .- fg); hF .+= tog(d.M, d.zbr .* (d.fgr .- fg))
        end
        Cown = ifelse.(a .> 1e-6, clamp.((Ck .- hC) ./ max.(a, 1e-6), 0.0, 1.0), 0.0)
        dIk = (hIg .- Igref) .+ ifelse.(a .> 1e-6, a .* h .* (1 .- Cown) .- a .* heq .* (1 .- C0), 0.0)
        Wk = Ck .* rsl_k .- (hF .+ ifelse.(a .> 1e-6, a .* zeq .* (Cown .- C0), 0.0))
        (Cs = Cs, Ck = Ck, dIs = dIs, dIk = dIk, Ws = Ws, Wk = Wk, Ls = (RI / RW) .* dIs .+ Ws, Lk = (RI / RW) .* dIk .+ Wk,
         bs = rs("dz_bed", n), bk = rk("dz_bed", r), hC = hC)
    end

    regs = [("host NH", (w .> 0.5) .& (LAT .> 0)), ("host ANT", (w .> 0.5) .& (LAT .< 0)), ("rim", (w .> 0) .& (w .<= 0.5)), ("outside", w .== 0)]
    io = IOBuffer()
    println(io, "# G2b diagnosis, ", tag, " (lmax ", lmax, ")\n\nbase: ", base, "\n")
    println(io, "Load L in m of sea water (kryos reconstructed from the ice-grid records, see the script header); differences kryos - standalone.\n")
    println(io, "| t | region | int dL [1e12 m3] | int abs dL [1e12 m3] | of which ice (rho_i/rho_w dI_g): int / int abs | water: int / int abs | int dC [1e9 m2] | int abs dC | rms d dz_bed [m] | rms local isostatic -(rho_w/rho_m) dL [m] |")
    println(io, "|---|---|---|---|---|---|---|---|---|---|")
    sel = unique([argmin(abs.(tr .- t)) for t in (-21000, -20000, -18000, -15000, -12000, -10000, -9000, -8000, -7000, -5000, -2500, 0)])
    # the time integral of the load difference (m yr) and the series of the regional integrals
    cum = zeros(nphi, nlat); ser = Dict(nm => zeros(nrec, 4) for (nm, _) in regs)
    snaps = Dict{Int,Any}()
    tmap = [argmin(abs.(tr .- t)) for t in (-21000, -15000, -12000, -10000, -8000, -5000, 0)]
    for r in 1:nrec
        L = loads(r); dL = L.Lk .- L.Ls
        r > 1 && (cum .+= dL .* (tr[r] - tr[r-1]))
        for (nm, m) in regs
            ser[nm][r, :] .= (wsum(dL, area, m), wsum(abs.(dL), area, m), wsum(L.Ck .- L.Cs, area, m), wrms(L.bk .- L.bs, area, m))
        end
        r in tmap && (snaps[r] = (L, copy(cum)))
        r in sel || continue
        for (nm, m) in regs[1:2]
            dI = (RI / RW) .* (L.dIk .- L.dIs); dW = L.Wk .- L.Ws; dC = L.Ck .- L.Cs
            println(io, @sprintf("| %.0f | %s | %.1f | %.1f | %.1f / %.1f | %.1f / %.1f | %.0f | %.0f | %.2f | %.2f |", tr[r], nm,
                wsum(dL, area, m) / 1e12, wsum(abs.(dL), area, m) / 1e12, wsum(dI, area, m) / 1e12, wsum(abs.(dI), area, m) / 1e12,
                wsum(dW, area, m) / 1e12, wsum(abs.(dW), area, m) / 1e12, wsum(dC, area, m) / 1e9, wsum(abs.(dC), area, m) / 1e9,
                wrms(L.bk .- L.bs, area, m), wrms((RW / RM) .* dL, area, m)))
        end
    end

    # the ocean function difference at the seed: where (present-day shelf, ice, other)
    L1 = loads(1); dC = L1.Ck .- L1.Cs
    shelf = (zeq .< 0) .& (zeq .> -200) .& (L1.Cs .< 0.5); ice = (L1.dIs .> 1) .| (L1.dIk .> 1)
    println(io, "\nOcean function at the seed, kryos - standalone [1e9 m2] (int / int abs), by cell class:\n")
    println(io, "| region | all | standalone land on the present-day shelf (z_bed_eq in [-200, 0]) | ice in either run | other |")
    println(io, "|---|---|---|---|---|")
    for (nm, m) in regs[1:2]
        f(mm) = @sprintf("%.0f / %.0f", wsum(dC, area, mm) / 1e9, wsum(abs.(dC), area, mm) / 1e9)
        println(io, "| ", nm, " | ", f(m), " | ", f(m .& shelf .& .!ice), " | ", f(m .& ice), " | ", f(m .& .!shelf .& .!ice), " |")
    end

    # correlation of the time-integrated load difference with the dz_bed difference at the end
    rN = nrec; LN, cN = snaps[tmap[end]]
    dB = LN.bk .- LN.bs
    for (nm, m) in regs[1:2]
        c = cor(vec(cN[m]), vec(dB[m]))
        println(io, @sprintf("\n%s: spatial correlation of int_t dL dt (21 ka to 0) with d dz_bed at 0 ka: %.2f", nm, c))
    end
    if ref !== nothing
        k2 = NCDataset(joinpath(base, ref, "kryos", "kryos_out.nc")); v2 = NCDataset(joinpath(base, ref, "vilma", "vilma_out.nc"))
        println(io, "\n## dz_bed difference, ", tag, " against ", ref, "\n\n| t | region | rms ", tag, " | rms ", ref, " |\n|---|---|---|---|")
        for r in sel
            n = stepof(r) + 1
            d1 = rk("dz_bed", r) .- rs("dz_bed", n)
            d2 = Float64.(nomissing(k2["dz_bed"][:, :, r], 0.0)) .- fix(Float64.(nomissing(v2["dz_bed"][:, :, n], 0.0)))
            for (nm, m) in regs[1:2]
                println(io, @sprintf("| %.0f | %s | %.2f | %.2f |", tr[r], nm, wrms(d1, area, m), wrms(d2, area, m)))
            end
        end
        close(k2); close(v2)
    end
    md = joinpath(out, "g2b_diag_$(tag).md"); write(md, String(take!(io)))

    # the maps
    f = joinpath(out, "g2b_diag_$(tag).nc"); isfile(f) && rm(f)
    nc = NCDataset(f, "c")
    defDim(nc, "lon", nphi); defDim(nc, "lat", nlat); defDim(nc, "time", length(tmap))
    defVar(nc, "lon", lon, ("lon",)); defVar(nc, "lat", lat, ("lat",)); defVar(nc, "time", tr[tmap], ("time",))
    defVar(nc, "host_w", w, ("lon", "lat")); defVar(nc, "z_bed_eq", zeq, ("lon", "lat"))
    for (nm, fn) in (("dz_bed_k", L -> L.bk), ("dz_bed_s", L -> L.bs), ("C_k", L -> L.Ck), ("C_s", L -> L.Cs), ("host_C", L -> L.hC),
                     ("L_k", L -> L.Lk), ("L_s", L -> L.Ls), ("dIg_k", L -> L.dIk), ("dIg_s", L -> L.dIs), ("W_k", L -> L.Wk), ("W_s", L -> L.Ws))
        x = zeros(Float32, nphi, nlat, length(tmap))
        for (q, r) in enumerate(tmap); x[:, :, q] .= fn(snaps[r][1]); end
        defVar(nc, nm, x, ("lon", "lat", "time"); deflatelevel = 2)
    end
    x = zeros(Float32, nphi, nlat, length(tmap))
    for (q, r) in enumerate(tmap); x[:, :, q] .= snaps[r][2]; end
    defVar(nc, "int_dL", x, ("lon", "lat", "time"); deflatelevel = 2, attrib = ["units" => "m yr", "long_name" => "time integral of L_k - L_s since the seed"])
    s = zeros(nrec, length(regs), 4)
    for (q, (nm, _)) in enumerate(regs); s[:, q, :] .= ser[nm]; end
    defDim(nc, "rec", nrec); defDim(nc, "region", length(regs)); defDim(nc, "stat", 4)
    defVar(nc, "rec", tr, ("rec",))
    defVar(nc, "series", s, ("rec", "region", "stat"); attrib = ["comment" => "regions: " * join(first.(regs), ", ") * "; stats: int dL [m3], int abs dL [m3], int dC [m2], rms d dz_bed [m]"])
    close(nc)
    println("wrote ", md, " and ", f)
end
main(copy(ARGS))
