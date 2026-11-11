# The fit of chion's ITM (itm_c, itm_t) per domain from an ensemble of runs on imposed ice (ice backend data):
#   <ens>/params.txt   runme's parameter file (header: group.key ...; one row per member, member i = row i, 0-based)
#   <ens>/pm/<i>/      par/kryos_NH-ANT_data_paleoemu.nml (PaleoMIST), with <ens>/pm/<i>/budget/basins.csv from
#                      smb_basin_budget.jl (tmin, tmax around the LGM window are enough)
#   <ens>/pd/<i>/      par/kryos_GRL-ANT_data_paleoemu_pd.nml (the present day on the reference topography)
# The cost of a domain (each domain has its own ITM group, so the domains are fitted independently):
#   J = (rms_land/sB)² + (rms_marine/sB)² + wPD (rmse_PD/sPD)²
#   rms_land   = sqrt(Σ A (B/A)² / Σ A) over the land-terminating basins of the LGM slices (B = ∫(smb - ρ_i ∂H/∂t): ≈ 0)
#   rms_marine = the same of min(B, 0) over the marine basins (B ≥ 0)
#   rmse_PD    = the rms difference of smb from the regional climate model on the ice of the PD run [kg m-2 yr-1]:
#                Greenland MAR v3.14 1961-1990 (the emulator's climate is pre-industrial), Antarctica RACMO2.3p2 1981-2010
# Basins in the regions of `exclude` (smb_basin_budget.jl's region names) are left out of the LGM terms.
#
#   julia --project=tools/levante tools/levante/smb_itm_fit.jl <ens dir> <out dir> [lgm=-23:-19] [sB=200] [sPD=200]
#         [wPD=1] [exclude=Region,Region]
# The budgets first, one job per member (tools/levante/jl.sh):
#   for i in $(seq 0 14); do tools/levante/jl.sh <ens>/pm/$i/budget.log tools/levante/smb_basin_budget.jl \
#       <ens>/pm/$i <ens>/pm/$i/budget tmin=-25 tmax=-17.5; done
# Writes <out>/fit.csv (per member and domain: the terms and J), <out>/fit.png (J and its terms over itm_c, itm_t), and
# prints the best member per domain and the minimum of a quadratic fit of J(itm_c, itm_t) when it lies inside the grid.
using NCDatasets, CairoMakie, Printf, Statistics, DelimitedFiles, LinearAlgebra

const FESMDATA = get(ENV, "FESMDATA", "/work/ba1442/fesmdata")
# the domains: suffix of the PaleoMIST run, suffix of the PD run, the ITM group, the RCM of the PD run
const DOMAINS = [(pm="nor", pd="grl", group="itm_chion",
                  rcm="$FESMDATA/Greenland/GRL-16KM/GRL-16KM_SMB-MARv3.14-1961-1990-annual-mean.nc", rcm_name="MAR 1961-1990"),
                 (pm="ant", pd="ant", group="itm_chion_ant",
                  rcm="$FESMDATA/Antarctica/ANT-32KM/ANT-32KM_SMB-RACMO2.3p2-1981-2010-annual-mean.nc", rcm_name="RACMO2.3p2 1981-2010")]

last2(ds, v) = (a = Array(ds[v]); Float64.(coalesce.(ndims(a) == 3 ? a[:, :, end] : a, NaN)))

"""The LGM terms of one member and domain from basins.csv."""
function lgm_terms(file, s, lgm, exclude)
    m = readdlm(file, ','; header=true); d, h = m[1], vec(m[2]); col(n) = findfirst(==(n), h)
    rows = [r for r in axes(d, 1) if d[r, col("domain")] == s && lgm[1] <= d[r, col("time_ka")] <= lgm[2] &&
            !(String(d[r, col("region")]) in exclude)]
    isempty(rows) && error("$file: no basin of $s in $(lgm) ka")
    A = Float64.(d[rows, col("area_km2")]); B = Float64.(d[rows, col("B_per_area")]); land = d[rows, col("type")] .!= "marine"
    rl = any(land) ? sqrt(sum(A[land] .* B[land] .^ 2) / sum(A[land])) : 0.0
    rm = any(.!land) ? sqrt(sum(A[.!land] .* min.(B[.!land], 0.0) .^ 2) / sum(A[.!land])) : 0.0
    nt = length(unique(d[rows, col("time_ka")]))
    smb = sum(Float64.(d[rows, col("smb_Gt")]); init=0.0) / nt; Bl = sum(Float64.(d[rows[land], col("B_Gt")]); init=0.0) / nt
    return (rms_land=rl, rms_marine=rm, smb_lgm=smb, B_land_lgm=Bl)
end

"""The PD terms of one member and domain: smb against the RCM on the ice of both."""
function pd_terms(dir, dom)
    m = NCDataset(joinpath(dir, "surface_$(dom.pd).nc")); i = NCDataset(joinpath(dir, "ice_$(dom.pd).nc")); o = NCDataset(dom.rcm)
    x = Float64.(m["xc"][:]); A = ((x[2] - x[1]) * 1e3)^2
    smb, so, H = last2(m, "smb"), last2(o, "smb"), last2(i, "H_ice")
    ok = (H .> 0) .& isfinite.(smb) .& isfinite.(so)
    r = (rmse_pd=sqrt(mean((smb[ok] .- so[ok]) .^ 2)), smb_pd=sum(smb[ok]) * A * 1e-12, smb_rcm=sum(so[ok]) * A * 1e-12,
         abl_pd=count(smb[ok] .< 0) / count(ok), abl_rcm=count(so[ok] .< 0) / count(ok))
    close(m); close(i); close(o)
    return r
end

"""Least-squares quadratic J(c, t); its stationary point if it is a minimum inside the box, else nothing."""
function quad_min(c, t, J)
    X = hcat(ones(length(c)), c, t, c .^ 2, c .* t, t .^ 2); length(unique(c)) >= 3 && length(unique(t)) >= 3 || return nothing
    a = X \ J
    Hs = [2a[4] a[5]; a[5] 2a[6]]; isposdef(Symmetric(Hs)) || return nothing
    p = -(Hs \ [a[2], a[3]])
    (minimum(c) <= p[1] <= maximum(c) && minimum(t) <= p[2] <= maximum(t)) || return nothing
    return (itm_c=p[1], itm_t=p[2], J=dot(a, [1, p[1], p[2], p[1]^2, p[1] * p[2], p[2]^2]))
end

function main(args)
    length(args) >= 2 || (println("usage: smb_itm_fit.jl <ens dir> <out dir> [lgm=-23:-19] [sB=200] [sPD=200] [wPD=1] [exclude=R,R]"); exit(1))
    ens, out = args[1], args[2]; mkpath(out)
    kw = Dict(String(k) => String(v) for (k, v) in (split(a, "="; limit=2) for a in args[3:end]))
    issubset(keys(kw), ("lgm", "sB", "sPD", "wPD", "exclude")) || error("unknown option(s): ", setdiff(keys(kw), ("lgm", "sB", "sPD", "wPD", "exclude")))
    lgm = Tuple(parse.(Float64, split(get(kw, "lgm", "-23:-19"), ":")))
    sB, sPD, wPD = (parse(Float64, get(kw, k, d)) for (k, d) in (("sB", "200"), ("sPD", "200"), ("wPD", "1")))
    exclude = filter(!isempty, split(get(kw, "exclude", ""), ","))
    P = readdlm(joinpath(ens, "params.txt"); header=true); pv, ph = Float64.(P[1]), String.(vec(P[2]))
    io = open(joinpath(out, "fit.csv"), "w")
    println(io, "member,domain,itm_c,itm_t,rms_land,rms_marine,smb_lgm_Gt,B_land_lgm_Gt,rmse_pd,smb_pd_Gt,smb_rcm_Gt,abl_pd,abl_rcm,J")
    fig = Figure(size=(1500, 420 * length(DOMAINS)))
    for (r, dom) in enumerate(DOMAINS)
        ic, it = findfirst(==("$(dom.group).itm_c"), ph), findfirst(==("$(dom.group).itm_t"), ph)
        (ic === nothing || it === nothing) && error("params.txt has no $(dom.group).itm_c / itm_t")
        res = NamedTuple[]
        for i in axes(pv, 1)
            mb = joinpath(ens, "pm", string(i - 1), "budget", "basins.csv")
            isfile(mb) || error("$mb missing: run smb_basin_budget.jl on member $(i - 1) first")
            l = lgm_terms(mb, dom.pm, lgm, exclude); p = pd_terms(joinpath(ens, "pd", string(i - 1)), dom)
            J = (l.rms_land / sB)^2 + (l.rms_marine / sB)^2 + wPD * (p.rmse_pd / sPD)^2
            push!(res, merge((member=i - 1, itm_c=pv[i, ic], itm_t=pv[i, it]), l, p, (J=J,)))
            q = res[end]
            @printf(io, "%d,%s,%.2f,%.2f,%.1f,%.1f,%.0f,%.0f,%.1f,%.0f,%.0f,%.3f,%.3f,%.4f\n", q.member, dom.pm, q.itm_c, q.itm_t,
                    q.rms_land, q.rms_marine, q.smb_lgm, q.B_land_lgm, q.rmse_pd, q.smb_pd, q.smb_rcm, q.abl_pd, q.abl_rcm, q.J)
        end
        b = res[argmin([q.J for q in res])]
        @printf("%s (%s): best member %d  itm_c %.1f itm_t %.1f  J %.3f  [LGM rms land %.0f marine %.0f, ∫smb %.0f Gt/yr; PD rmse %.0f, ∫smb %.0f vs %s %.0f Gt/yr]\n",
                dom.pm, dom.group, b.member, b.itm_c, b.itm_t, b.J, b.rms_land, b.rms_marine, b.smb_lgm, b.rmse_pd, b.smb_pd, dom.rcm_name, b.smb_rcm)
        c, t = [q.itm_c for q in res], [q.itm_t for q in res]
        qm = quad_min(c, t, [q.J for q in res])
        qm === nothing ? println("  quadratic fit: no minimum inside the grid (widen the ranges)") :
                         @printf("  quadratic fit: itm_c %.1f itm_t %.2f (J %.3f)\n", qm.itm_c, qm.itm_t, qm.J)
        for (k, (f, ti)) in enumerate(((:J, "$(dom.pm): J"), (:rms_land, "LGM rms land [kg m-2 yr-1]"),
                                        (:rms_marine, "LGM rms marine"), (:rmse_pd, "PD rmse vs $(dom.rcm_name)")))
            ax = Axis(fig[r, 2k - 1]; xlabel="itm_c [W m-2]", ylabel="itm_t [W m-2 K-1]", title=ti)
            v = [getfield(q, f) for q in res]
            sc = scatter!(ax, c, t; color=v, markersize=28, marker=:rect, colormap=:viridis)
            k == 1 && scatter!(ax, [b.itm_c], [b.itm_t]; marker=:star5, color=:red, markersize=20)
            k == 1 && qm !== nothing && scatter!(ax, [qm.itm_c], [qm.itm_t]; marker=:xcross, color=:white, markersize=18)
            Colorbar(fig[r, 2k], sc)
        end
    end
    close(io); save(joinpath(out, "fit.png"), fig)
    println("wrote $(out)/fit.csv, fit.png")
end

main(ARGS)
