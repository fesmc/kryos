# The fit of chion's parameters per domain from an ensemble of runs on imposed ice (ice backend data):
#   <ens>/params.txt   runme's parameter file (header: group.key ...; one row per member, member i = row i, 0-based)
#   <ens>/pm/<i>/      par/kryos_NH-ANT_data_paleoemu.nml (PaleoMIST: an LGM snapshot or the transient), with
#                      <ens>/pm/<i>/budget/basins.csv from smb_basin_budget.jl (smb_ensemble.sh post)
#   <ens>/pd/<i>/      par/kryos_GRL-ANT_data_paleoemu_pd.nml (the present day on the reference topography)
# The parameters of a domain are the columns of params.txt: a group ending in _ant is Antarctica's, its counterpart
# without _ant North's (with Greenland's PD anchor), any other group (e.g. clim_glob) both domains'. The domains are fitted
# independently. The cost of a domain:
#   J = (rms_land/sB)² + (rms_marine/sB)² + wPD (rmse_PD/sPD)² + wI ((∫smb_PD - ∫smb_RCM)/sI)² + wA ((A_abl - A_abl,RCM)/sA)²
#   rms_land   = sqrt(Σ A (B/A)² / Σ A) over the land-terminating basins of the LGM slices (B = ∫(smb - ρ_i ∂H/∂t): ≈ 0)
#   rms_marine = the same of min(B, 0) over the marine basins (B ≥ 0)
#   rmse_PD    = the rms difference of smb from the regional climate model on the ice of the PD run [kg m-2 yr-1]:
#                Greenland MAR v3.14 and Antarctica RACMO2.3p2, 1981-2010 like the PD run's reference climate
#   ∫smb_PD    = the integrated smb of the PD run on the ice of both [Gt/yr]; A_abl its ablation area (smb < 0) as a share
#                of that ice: the gridpoint rmse alone can prefer a too narrow ablation zone (wI, wA default 0)
#   With fice=true (default) the PD run's smb is weighted by the RCM's ice-covered fraction f_ice of the cell where the
#   RCM has values off the ice (MAR: the tundra): its value on the grid is then the mean over the cell (FesmData's
#   remap), ice and ice-free, and the ice-free part's smb is about 0 (seasonal snow), so f_ice smb is its counterpart.
#   Unweighted, Greenland's margin cells, partly ice-free, read as too much melt (up to 25 % at 600 m). RACMO2.3p2 has
#   values on the ice only (grounded and floating; its f_ice is the grounded share), so Antarctica is not weighted.
#   fice=false compares the PD run's smb as it is everywhere.
# Basins in the regions of `exclude` (smb_basin_budget.jl's region names, e.g. Cordillera_and_Alaska,British_Isles: small
# maritime ice sheets whose precipitation the emulator's 3.75 deg grid does not resolve) are left out of the LGM terms.
#
#   julia --project=tools/levante tools/levante/smb_fit.jl <ens dir> <out dir> [lgm=-23:-19] [sB=200] [sPD=200]
#         [wPD=1] [wI=0] [sI=100] [wA=0] [sA=0.05] [exclude=Region,Region] [fice=true]
# Writes <out>/fit.csv (per member and domain: its parameters, the terms and J), <out>/fit.png (J and its terms against
# each parameter), and prints the best member per domain and the minimum of a quadratic fit of J over the domain's
# parameters (scaled to [0, 1]) when it lies inside the box.
using NCDatasets, CairoMakie, Printf, Statistics, DelimitedFiles, LinearAlgebra

const FESMDATA = get(ENV, "FESMDATA", "/work/ba1442/fesmdata")
# the domains: suffix of the PaleoMIST run, suffix of the PD run, whether its groups end in _ant, the RCM of the PD run
# and whether the RCM has values off the ice, so that its cell means are diluted by the ice-free part (fice)
const DOMAINS = [(pm="nor", pd="grl", ant=false,
                  rcm="$FESMDATA/Greenland/GRL-16KM/GRL-16KM_SMB-MARv3.14-1981-2010-annual-mean.nc", rcm_name="MAR 1981-2010",
                  fice=true),
                 (pm="ant", pd="ant", ant=true,
                  rcm="$FESMDATA/Antarctica/ANT-32KM/ANT-32KM_SMB-RACMO2.3p2-1981-2010-annual-mean.nc", rcm_name="RACMO2.3p2 1981-2010",
                  fice=false)]

"""The columns of params.txt that are the domain's parameters (see the header)."""
function domain_columns(ph, dom)
    grp(c) = String(split(c, ".")[1])
    antgroups = Set(grp(c) for c in ph if endswith(grp(c), "_ant"))
    northgroups = Set(g[1:end-4] for g in antgroups)
    keep(c) = dom.ant ? !(grp(c) in northgroups) : !(grp(c) in antgroups)
    return [k for (k, c) in enumerate(ph) if keep(c)]
end

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

"""The PD terms of one member and domain: smb (times the RCM's f_ice with `fice`) against the RCM on the PD run's ice."""
function pd_terms(dir, dom, fice)
    m = NCDataset(joinpath(dir, "surface_$(dom.pd).nc")); i = NCDataset(joinpath(dir, "ice_$(dom.pd).nc")); o = NCDataset(dom.rcm)
    x = Float64.(m["xc"][:]); A = ((x[2] - x[1]) * 1e3)^2
    smb, so, H = last2(m, "smb"), last2(o, "smb"), last2(i, "H_ice")
    fice && dom.fice && (smb = smb .* replace(last2(o, "f_ice"), NaN => 0.0))
    ok = (H .> 0) .& isfinite.(smb) .& isfinite.(so)
    r = (rmse_pd=sqrt(mean((smb[ok] .- so[ok]) .^ 2)), smb_pd=sum(smb[ok]) * A * 1e-12, smb_rcm=sum(so[ok]) * A * 1e-12,
         abl_pd=count(smb[ok] .< 0) / count(ok), abl_rcm=count(so[ok] .< 0) / count(ok))
    close(m); close(i); close(o)
    return r
end

"""
Least-squares quadratic J(x) over the parameters `X` [member, parameter] (each scaled to [0, 1]); its stationary point,
in the parameters' units, if it is a minimum inside the box, else nothing. Parameters with fewer than 3 values are left out.
"""
function quad_min(X, J)
    v = [k for k in axes(X, 2) if length(unique(X[:, k])) >= 3]
    isempty(v) && return nothing
    lo, hi = minimum(X[:, v]; dims=1), maximum(X[:, v]; dims=1); U = (X[:, v] .- lo) ./ (hi .- lo); n = length(v)
    pairs = [(i, j) for i in 1:n for j in i:n]
    basis(u) = vcat(1.0, u, [u[i] * u[j] for (i, j) in pairs])
    M = reduce(vcat, (basis(U[m, :])' for m in axes(U, 1)))
    size(M, 1) > size(M, 2) || return nothing
    a = M \ J
    g = a[2:n+1]; H = zeros(n, n)
    for (k, (i, j)) in enumerate(pairs)
        c = a[n+1+k]; i == j ? (H[i, i] = 2c) : (H[i, j] = c; H[j, i] = c)
    end
    isposdef(Symmetric(H)) || return nothing
    u = -(H \ g); all(0 .<= u .<= 1) || return nothing
    x = copy(X[1, :]); x[v] = vec(lo) .+ u .* vec(hi .- lo)
    return (x=x, vars=v, J=dot(a, basis(u)))
end

function main(args)
    length(args) >= 2 || (println("usage: smb_fit.jl <ens dir> <out dir> [lgm=-23:-19] [sB=200] [sPD=200] [wPD=1] [wI=0] [sI=100] [wA=0] [sA=0.05] [exclude=R,R] [fice=true]"); exit(1))
    ens, out = args[1], args[2]; mkpath(out)
    kw = Dict(String(k) => String(v) for (k, v) in (split(a, "="; limit=2) for a in args[3:end]))
    known = ("lgm", "sB", "sPD", "wPD", "wI", "sI", "wA", "sA", "exclude", "fice")
    issubset(keys(kw), known) || error("unknown option(s): ", setdiff(keys(kw), known))
    lgm = Tuple(parse.(Float64, split(get(kw, "lgm", "-23:-19"), ":")))
    sB, sPD, wPD, wI, sI, wA, sA = (parse(Float64, get(kw, k, d)) for (k, d) in
        (("sB", "200"), ("sPD", "200"), ("wPD", "1"), ("wI", "0"), ("sI", "100"), ("wA", "0"), ("sA", "0.05")))
    exclude = filter(!isempty, split(get(kw, "exclude", ""), ","))
    fice = parse(Bool, get(kw, "fice", "true"))
    P = readdlm(joinpath(ens, "params.txt"); header=true); pv, ph = Float64.(P[1]), String.(vec(P[2]))
    io = open(joinpath(out, "fit.csv"), "w")
    println(io, "member,domain,", join(ph, ","), ",rms_land,rms_marine,smb_lgm_Gt,B_land_lgm_Gt,rmse_pd,smb_pd_Gt,smb_rcm_Gt,abl_pd,abl_rcm,J")
    cols = [domain_columns(ph, dom) for dom in DOMAINS]
    fig = Figure(size=(380 * (maximum(length, cols) + 1), 380 * length(DOMAINS)))
    for (r, dom) in enumerate(DOMAINS)
        cs = cols[r]; isempty(cs) && error("params.txt has no parameter of $(dom.pm)")
        res = NamedTuple[]
        for i in axes(pv, 1)
            mb = joinpath(ens, "pm", string(i - 1), "budget", "basins.csv")
            isfile(mb) || error("$mb missing: run smb_basin_budget.jl on member $(i - 1) first (smb_ensemble.sh post)")
            l = lgm_terms(mb, dom.pm, lgm, exclude); p = pd_terms(joinpath(ens, "pd", string(i - 1)), dom, fice)
            J = (l.rms_land / sB)^2 + (l.rms_marine / sB)^2 + wPD * (p.rmse_pd / sPD)^2 +
                wI * ((p.smb_pd - p.smb_rcm) / sI)^2 + wA * ((p.abl_pd - p.abl_rcm) / sA)^2
            push!(res, merge((member=i - 1,), l, p, (J=J,)))
            q = res[end]
            @printf(io, "%d,%s,%s,%.1f,%.1f,%.0f,%.0f,%.1f,%.0f,%.0f,%.3f,%.3f,%.4f\n", q.member, dom.pm,
                    join([k in cs ? @sprintf("%g", pv[i, k]) : "" for k in eachindex(ph)], ","),
                    q.rms_land, q.rms_marine, q.smb_lgm, q.B_land_lgm, q.rmse_pd, q.smb_pd, q.smb_rcm, q.abl_pd, q.abl_rcm, q.J)
        end
        Jv = [q.J for q in res]; kb = argmin(Jv); b = res[kb]
        @printf("%s: best member %d  %s  J %.3f  [LGM rms land %.0f marine %.0f, ∫smb %.0f Gt/yr; PD rmse %.0f, ∫smb %.0f vs %s %.0f Gt/yr]\n",
                dom.pm, b.member, join(["$(ph[k]) $(pv[kb, k])" for k in cs], ", "), b.J, b.rms_land, b.rms_marine, b.smb_lgm,
                b.rmse_pd, b.smb_pd, dom.rcm_name, b.smb_rcm)
        qm = quad_min(pv[:, cs], Jv)
        qm === nothing ? println("  quadratic fit: no minimum inside the box (widen the ranges)") :
                         @printf("  quadratic fit: %s (J %.3f)\n", join(["$(ph[cs[k]]) $(round(qm.x[k]; sigdigits=4))" for k in qm.vars], ", "), qm.J)
        for (c, k) in enumerate(cs)
            ax = Axis(fig[r, c]; xlabel=ph[k], ylabel="J", title=c == 1 ? "$(dom.pm): J (colour: LGM rms land)" : "")
            sc = scatter!(ax, pv[:, k], Jv; color=[q.rms_land for q in res], colormap=:viridis)
            scatter!(ax, [pv[kb, k]], [b.J]; marker=:star5, color=:red, markersize=18)
            c == length(cs) && Colorbar(fig[r, c+1], sc)
        end
    end
    close(io); save(joinpath(out, "fit.png"), fig)
    println("wrote $(out)/fit.csv, fit.png")
end

main(ARGS)
