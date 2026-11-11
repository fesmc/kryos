# The mass budget of an imposed ice sheet per drainage basin: is a surface mass balance consistent with the ice it is
# computed on? A run with the ice backend data (par/kryos_NH-ANT_data_paleoemu.nml, par/kryos_GRL-ANT_data_paleoemu_pd.nml)
# gives smb on a prescribed geometry; without velocities the test is the budget of each basin:
#   B = ∫ (smb - ρ_i ∂H/∂t) dA      [Gt/yr], the mass that must leave the basin across its margin;
#   B ≈ 0 for a land-terminating basin (no calving), B ≥ 0 for a marine one (calving and the shelves take it).
# Basins: steepest descent on z_srf over the grounded ice, after filling its pits (priority flood), each cell draining to
# the margin cell its path leaves the ice from; basins below the minimum area are merged into the neighbour they share the
# longest border with. A basin is marine when most of its balance flux leaves into the ocean (bed below sea level, or
# floating ice), edge when it leaves across the grid edge, land otherwise.
# Also per basin: the accumulation-area ratio AAR, and the balance flux q accumulated downslope with its speed
# u = q/(ρ_i H dx); q < 0 inside the ice (more net ablation upstream than supply) marks an inconsistent smb or geometry.
#
#   julia --project=tools/levante tools/levante/smb_basin_budget.jl <run dir> <out dir> [amin_km2=20000] [rho_i=910]
#
# Reads every pair ice_<s>.nc + surface_<s>.nc of the run (the same grid: the ice instance on the hub grid of chion) at
# their common output times. Writes to <out dir>:
#   basins.csv   per domain, time and basin: type, marine flux fraction, area, ∫smb, ∫ρ_i dH/dt, B, B/area, AAR, the area
#                fraction with q < 0
#   domains.csv  per domain and time: totals, and the misfits rms_land = sqrt(Σ A (B/A)² / Σ A) over land basins and
#                rms_marine (the same of min(B, 0)) [kg m-2 yr-1], the terms of the fit (smb_itm_fit.jl)
#   basins_<s>.nc  basin id, its type, the balance flux q [Gt/yr] and speed u_bal [m/yr] per time
#   budget_<s>.png the budget through time and the maps at the time nearest 20 ka (or the last)
using NCDatasets, CairoMakie, Printf, Statistics

const NB = ((-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1))
const TYPES = ("land", "marine", "edge")

# ---- a minimal binary min-heap of (key, index) -------------------------------------------------------------------
struct Heap
    k::Vector{Float64}
    v::Vector{Int}
end
Heap() = Heap(Float64[], Int[])
Base.isempty(h::Heap) = isempty(h.k)
function Base.push!(h::Heap, key::Float64, val::Int)
    push!(h.k, key); push!(h.v, val); i = length(h.k)
    while i > 1 && h.k[i ÷ 2] > h.k[i]
        p = i ÷ 2; h.k[i], h.k[p] = h.k[p], h.k[i]; h.v[i], h.v[p] = h.v[p], h.v[i]; i = p
    end
end
function Base.pop!(h::Heap)
    key, val = h.k[1], h.v[1]
    h.k[1] = h.k[end]; h.v[1] = h.v[end]; pop!(h.k); pop!(h.v)
    n = length(h.k); i = 1
    while true
        l, r, m = 2i, 2i + 1, i
        l <= n && h.k[l] < h.k[m] && (m = l)
        r <= n && h.k[r] < h.k[m] && (m = r)
        m == i && break
        h.k[i], h.k[m] = h.k[m], h.k[i]; h.v[i], h.v[m] = h.v[m], h.v[i]; i = m
    end
    return key, val
end

# ---- the drainage of one slice ------------------------------------------------------------------------------------
"""
    drainage(G, zs, ocean) -> (rcv, exit, zf)

On the grounded ice G (nx×ny): the filled surface zf (priority flood from the margin, a strictly increasing ε inland),
each cell's receiver rcv (linear index, 0 = it leaves the ice) by steepest descent on zf, and for a cell that leaves
the ice its exit type (1 land, 2 ocean, 3 grid edge).
"""
function drainage(G, zs, ocean)
    nx, ny = size(G); L = LinearIndices(G)
    zf = fill(NaN, nx, ny); seen = falses(nx, ny); h = Heap()
    margin(i, j) = any(((di, dj),) -> !(1 <= i + di <= nx && 1 <= j + dj <= ny) || !G[i+di, j+dj], NB)
    for j in 1:ny, i in 1:nx
        if G[i, j] && margin(i, j)
            seen[i, j] = true; zf[i, j] = zs[i, j]; push!(h, zf[i, j], L[i, j])
        end
    end
    C = CartesianIndices(G)
    while !isempty(h)
        z, c = pop!(h); i, j = Tuple(C[c])
        for (di, dj) in NB
            a, b = i + di, j + dj
            (1 <= a <= nx && 1 <= b <= ny && G[a, b] && !seen[a, b]) || continue
            seen[a, b] = true; zf[a, b] = max(zs[a, b], z + 1e-3); push!(h, zf[a, b], L[a, b])
        end
    end
    rcv = zeros(Int, nx, ny); exit = zeros(Int8, nx, ny)
    for j in 1:ny, i in 1:nx
        G[i, j] || continue
        best, r, e = 0.0, 0, Int8(0)
        for (di, dj) in NB
            a, b = i + di, j + dj; d = hypot(di, dj)
            if !(1 <= a <= nx && 1 <= b <= ny)
                s, rr, ee = Inf, 0, Int8(3)                         # the grid edge: always downhill
            elseif G[a, b]
                s, rr, ee = (zf[i, j] - zf[a, b]) / d, L[a, b], Int8(0)
            else
                s, rr, ee = (zf[i, j] - zs[a, b]) / d, 0, ocean[a, b] ? Int8(2) : Int8(1)
            end
            if s > best
                best, r, e = s, rr, ee
            end
        end
        if best == 0.0                                              # a margin pit: it leaves to its lowest outside neighbour
            zmin = Inf
            for (di, dj) in NB
                a, b = i + di, j + dj
                (1 <= a <= nx && 1 <= b <= ny && !G[a, b]) || continue
                zs[a, b] < zmin && (zmin = zs[a, b]; e = ocean[a, b] ? Int8(2) : Int8(1))
            end
            r = 0
        end
        rcv[i, j] = r; exit[i, j] = r == 0 ? e : Int8(0)
    end
    return rcv, exit, zf
end

"""Basin labels (the outlet cell's index), the cells ordered from the outlets inland."""
function basins_raw(G, rcv, zf)
    order = sort(findall(vec(G)); by=c -> zf[c])                  # receivers first (strictly lower zf)
    lab = zeros(Int, size(G))
    for c in order
        lab[c] = rcv[c] == 0 ? c : lab[rcv[c]]
    end
    return lab, order
end

"""Merge basins below amin cells into the neighbour with the longest shared border (smallest first)."""
function merge_small!(lab, G, amin_cells)
    nx, ny = size(G)
    while true
        ids = unique(lab[G]); n = Dict(b => count(==(b), lab) for b in ids)
        small = sort([b for b in ids if n[b] < amin_cells]; by=b -> n[b])
        merged = false
        for b in small
            border = Dict{Int,Int}()
            for c in findall(==(b), lab)
                i, j = Tuple(c)
                for (di, dj) in NB
                    a, bb = i + di, j + dj
                    (1 <= a <= nx && 1 <= bb <= ny && G[a, bb] && lab[a, bb] != b) || continue
                    border[lab[a, bb]] = get(border, lab[a, bb], 0) + 1
                end
            end
            isempty(border) && continue                              # an isolated ice cap: stays its own basin
            to = argmax(border); lab[lab .== b] .= to; merged = true
            break                                                    # sizes changed: recount
        end
        merged || return lab
    end
end

# ---- one domain ---------------------------------------------------------------------------------------------------
rd(ds, v) = Float64.(coalesce.(Array(ds[v]), NaN))

function domain_budget(run, s, out, amin_km2, rho_i, io_b, io_d)
    di = NCDataset(joinpath(run, "ice_$s.nc")); ds = NCDataset(joinpath(run, "surface_$s.nc"))
    x = Float64.(ds["xc"][:]); dxm = (x[2] - x[1]) * 1e3; A = dxm^2
    ti, ts = Float64.(di["time"][:]), Float64.(ds["time"][:])
    size(di["H_ice"])[1:2] == size(ds["smb"])[1:2] ||
        error("ice_$s.nc and surface_$s.nc are on different grids: run the ice on the grid of the surface")
    times = intersect(ti, ts); nt = length(times)
    nx, ny = size(di["H_ice"])[1:2]
    lab_t = zeros(Int32, nx, ny, nt); typ_t = zeros(Int8, nx, ny, nt); q_t = fill(NaN32, nx, ny, nt); u_t = fill(NaN32, nx, ny, nt)
    series = NamedTuple[]
    amin_cells = max(1, round(Int, amin_km2 * 1e6 / A))
    for (k, t) in enumerate(times)
        ki, ks = findfirst(==(t), ti), findfirst(==(t), ts)
        H, zs, zb, zl = (di[v][:, :, ki] |> a -> Float64.(coalesce.(a, NaN)) for v in ("H_ice", "z_srf", "z_bed", "z_sl"))
        fg = Float64.(coalesce.(di["f_grnd"][:, :, ki], 0.0)); dH = Float64.(coalesce.(di["dHidt"][:, :, ki], 0.0))
        smb = Float64.(coalesce.(ds["smb"][:, :, ks], NaN))          # mm w.e./yr = kg m-2 yr-1
        G = (H .> 0) .& (fg .> 0.5) .& isfinite.(smb)
        ocean = (zb .< zl) .| ((H .> 0) .& (fg .<= 0.5))
        rcv, exit, zf = drainage(G, zs, ocean)
        lab, order = basins_raw(G, rcv, zf)
        merge_small!(lab, G, amin_cells)
        net = smb .- rho_i .* dH                                     # kg m-2 yr-1
        q = zeros(nx, ny)                                            # kg/yr through each cell, donors first
        for c in Iterators.reverse(order)
            q[c] += net[c] * A
            rcv[c] > 0 && (q[rcv[c]] += q[c])
        end
        u = ifelse.(G, q ./ (rho_i .* max.(H, 1.0) .* dxm), NaN)
        ids = sort(unique(lab[G]))
        tot = (A=0.0, smb=0.0, dh=0.0, Bl=0.0, Bm=0.0, Al=0.0, Am=0.0, wl=0.0, wm=0.0, acc=0.0)
        for b in ids
            m = lab .== b
            outl = m .& (rcv .== 0)
            qo = [sum(max.(q[outl .& (exit .== e)], 0.0)) for e in 1:3]
            sq = sum(qo); ty = sq > 0 ? argmax(qo) : (any(exit[outl] .== 2) ? 2 : 1)
            mf = sq > 0 ? qo[2] / sq : 0.0
            ar = count(m) * A; Ssmb = sum(smb[m]) * A; Sdh = rho_i * sum(dH[m]) * A; B = Ssmb - Sdh
            aar = count(m .& (smb .> 0)) / count(m); qneg = count(m .& (q .< 0)) / count(m)
            @printf(io_b, "%s,%.1f,%d,%s,%.3f,%.1f,%.3f,%.3f,%.3f,%.2f,%.3f,%.3f\n", s, t / 1e3, b, TYPES[ty], mf,
                    ar / 1e6, Ssmb * 1e-12, Sdh * 1e-12, B * 1e-12, B / ar, aar, qneg)
            view(typ_t, :, :, k)[m] .= ty
            tot = merge(tot, (A=tot.A + ar, smb=tot.smb + Ssmb, dh=tot.dh + Sdh, acc=tot.acc + count(m .& (smb .> 0)) * A))
            if ty == 2
                tot = merge(tot, (Bm=tot.Bm + B, Am=tot.Am + ar, wm=tot.wm + ar * (min(B, 0.0) / ar)^2))
            else
                tot = merge(tot, (Bl=tot.Bl + B, Al=tot.Al + ar, wl=tot.wl + ar * (B / ar)^2))
            end
        end
        rl = tot.Al > 0 ? sqrt(tot.wl / tot.Al) : 0.0; rm = tot.Am > 0 ? sqrt(tot.wm / tot.Am) : 0.0
        @printf(io_d, "%s,%.1f,%d,%.1f,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f,%.2f\n", s, t / 1e3, length(ids), tot.A / 1e6,
                tot.smb * 1e-12, tot.dh * 1e-12, tot.Bl * 1e-12, tot.Bm * 1e-12, tot.acc / max(tot.A, 1.0), rl, rm)
        push!(series, (t=t / 1e3, smb=tot.smb * 1e-12, dh=tot.dh * 1e-12, Bl=tot.Bl * 1e-12, Bm=tot.Bm * 1e-12, rl=rl, rm=rm))
        lab_t[:, :, k] .= ifelse.(G, lab, 0); q_t[:, :, k] .= ifelse.(G, q .* 1e-12, NaN); u_t[:, :, k] .= u
    end
    # the maps
    NCDataset(joinpath(out, "basins_$s.nc"), "c") do o
        defDim(o, "xc", nx); defDim(o, "yc", ny); defDim(o, "time", nt)
        defVar(o, "xc", x, ("xc",); attrib=["units" => "km"]); defVar(o, "yc", Float64.(ds["yc"][:]), ("yc",); attrib=["units" => "km"])
        defVar(o, "time", times, ("time",); attrib=["units" => "years"])
        defVar(o, "basin", lab_t, ("xc", "yc", "time"); attrib=["long_name" => "basin (its outlet cell; 0 no grounded ice)"])
        defVar(o, "type", typ_t, ("xc", "yc", "time"); attrib=["long_name" => "1 land, 2 marine, 3 edge"])
        defVar(o, "q", q_t, ("xc", "yc", "time"); attrib=["units" => "Gt yr-1", "long_name" => "balance flux (smb - rho_i dH/dt accumulated downslope)"])
        defVar(o, "u_bal", u_t, ("xc", "yc", "time"); attrib=["units" => "m yr-1", "long_name" => "balance speed q/(rho_i H dx)"])
    end
    # the figure
    k = argmin(abs.(times ./ 1e3 .+ 20.0)); t = times[k]
    smb = Float64.(coalesce.(ds["smb"][:, :, findfirst(==(t), ts)], NaN)); lab = lab_t[:, :, k]
    G = lab .> 0
    fig = Figure(size=(1600, 900))
    ax = Axis(fig[1, 1:4]; xlabel="time [ka]", ylabel="Gt/yr", title="$s: ∫smb, ∫ρ_i dH/dt, B over land and marine basins")
    for (f, l) in ((:smb, "∫smb"), (:dh, "∫ρ_i dH/dt"), (:Bl, "B land (→ 0)"), (:Bm, "B marine (≥ 0)"))
        scatterlines!(ax, [r.t for r in series], [getfield(r, f) for r in series]; label=l)
    end
    hlines!(ax, 0; color=:black, linewidth=0.5); axislegend(ax; position=:lb)
    Bmap = fill(NaN, nx, ny)
    for b in unique(lab[G])
        m = lab .== b; Bmap[m] .= (sum(smb[m]) - rho_i * sum(Float64.(coalesce.(di["dHidt"][:, :, findfirst(==(t), ti)], 0.0))[m])) / count(m)
    end
    for (c, (a, ti_, cr, cm)) in enumerate(((ifelse.(G, smb, NaN), "smb [kg m-2 yr-1], $(t / 1e3) ka", (-1500, 1500), :RdBu),
                                            (ifelse.(G, Float64.(lab), NaN), "basins", nothing, :tab20),
                                            (Bmap, "B/area per basin [kg m-2 yr-1]", (-500, 500), :RdBu),
                                            (log10.(abs.(u_t[:, :, k]) .+ 1), "log10 |u_bal| [m/yr] (q < 0 hatched grey)", (0, 3.5), :viridis)))
        ax = Axis(fig[2, c]; title=ti_, aspect=DataAspect()); hidedecorations!(ax)
        h = cr === nothing ? heatmap!(ax, a; colormap=cm) : heatmap!(ax, a; colorrange=cr, colormap=cm)
        if c == 4
            neg = ifelse.(G .& (q_t[:, :, k] .< 0), 1.0, NaN); heatmap!(ax, neg; colormap=[:grey60, :grey60])
        end
        cr === nothing || Colorbar(fig[3, c], h; vertical=false)
    end
    save(joinpath(out, "budget_$s.png"), fig)
    close(di); close(ds)
end

function main(args)
    length(args) >= 2 || (println("usage: smb_basin_budget.jl <run dir> <out dir> [amin_km2=20000] [rho_i=910]"); exit(1))
    run, out = args[1], args[2]
    amin = length(args) >= 3 ? parse(Float64, args[3]) : 2.0e4
    rho_i = length(args) >= 4 ? parse(Float64, args[4]) : 910.0
    mkpath(out)
    sfx = [m[1] for m in (match(r"^ice_(\w+)\.nc$", f) for f in readdir(run)) if m !== nothing &&
           !endswith(m[1], "_ts") && !endswith(m[1], "_plastic") && isfile(joinpath(run, "surface_$(m[1]).nc"))]
    isempty(sfx) && error("no pair ice_<s>.nc + surface_<s>.nc in $run")
    open(joinpath(out, "basins.csv"), "w") do io_b
        open(joinpath(out, "domains.csv"), "w") do io_d
            println(io_b, "domain,time_ka,basin,type,marine_frac,area_km2,smb_Gt,rho_dHdt_Gt,B_Gt,B_per_area,AAR,frac_q_neg")
            println(io_d, "domain,time_ka,n_basins,area_km2,smb_Gt,rho_dHdt_Gt,B_land_Gt,B_marine_Gt,AAR,rms_land,rms_marine")
            for s in sfx
                println("domain $s"); domain_budget(run, s, out, amin, rho_i, io_b, io_d)
            end
        end
    end
    println("wrote $(out)/basins.csv, domains.csv, basins_<s>.nc, budget_<s>.png")
end

main(ARGS)
