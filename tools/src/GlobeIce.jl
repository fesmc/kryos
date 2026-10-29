"""
    KryosTools.GlobeIce

The offline merge of the global geometry hub (phase4-plan 5.4, gate G1b (4)): what
`kryos_globe_ice` writes to `globe_ice.nc`, recomputed from its inputs and the
coupler's conservative map files, independently of the Fortran code.

- `read_scrip_con(file)`: a `con` map of the coupler's cache (`scrip-con_<src>_<dst>.nc`):
  links, weights and the coverage payload (`mp_wfull_*`, `mp_wcov_*`).
- `merge_offline(maps, bg, amount_kind)`: sources added as `OfflineSource` (a map, the
  model areas of the source and target grids, the field, owner weights) and combined
  with the background as `coupler_merge_end` does:

      extensive:  var = (1 - min(F,1)) bg + sum_d A_d,   A_dk = sum_j phi_jk o_j area_j v_j / area_k
      intensive:  var = (1 - min(F,1)) bg + sum_d A_d / max(F,1),  A_dk = r_k sum_j w_jk o_j v_j / wfull_k
      coverage:   F_dk = r_k sum_j w_jk o_j / wfull_k,   F = sum_d F_dk

  with `phi_jk = q_j w_jk / wfull_j`, `q_j = wcov_j / (sum of the links of j)` and
  `r_k = wcov_k / (sum of the links of k)`.
- `check_globe_ice(globe_nc, inputs_nc; tol)`: the last record of `H_ice` and the
  coverage `cov` of `globe_ice.nc` against `merge_offline` of the inputs file written
  next to it (variables `area_bg`, `H_ext`, and per domain `<d>` the fields
  `H_<d>`, `area_<d>`, `owner_<d>` (optional) and the global attribute `map_<d>`).
"""
module GlobeIce

using NCDatasets

export ScripMap, read_scrip_con, OfflineSource, merge_offline, check_globe_ice

struct ScripMap
    src::Vector{Int}          # 1-based source index of each link
    dst::Vector{Int}
    w::Vector{Float64}
    n_src::Int
    n_dst::Int
    wfull_src::Vector{Float64}
    wcov_src::Vector{Float64}
    wfull_dst::Vector{Float64}
    wcov_dst::Vector{Float64}
end

function read_scrip_con(file::AbstractString)
    NCDataset(file) do ds
        src = Int.(ds["src_address"][:])
        dst = Int.(ds["dst_address"][:])
        w = Float64.(vec(ds["remap_matrix"][:, :]))
        ns = length(ds["mp_wfull_src"][:])
        nd = length(ds["mp_wfull_dst"][:])
        (minimum(src) >= 1 && maximum(src) <= ns && minimum(dst) >= 1 && maximum(dst) <= nd) ||
            error("read_scrip_con: addresses of $file are not 1-based indices of the grids")
        ScripMap(src, dst, w, ns, nd, Float64.(ds["mp_wfull_src"][:]), Float64.(ds["mp_wcov_src"][:]),
                 Float64.(ds["mp_wfull_dst"][:]), Float64.(ds["mp_wcov_dst"][:]))
    end
end

"One domain of a merge: its map to the target grid, the model areas, the field and the owner weights."
struct OfflineSource
    map::ScripMap
    area_src::Vector{Float64}
    area_dst::Vector{Float64}
    v::Vector{Float64}
    owner::Vector{Float64}
end

function amounts(s::OfflineSource, kind::Symbol)
    m = s.map
    s_src = zeros(m.n_src); s_dst = zeros(m.n_dst)
    for l in eachindex(m.w)
        m.w[l] > 0 || continue
        s_src[m.src[l]] += m.w[l]
        s_dst[m.dst[l]] += m.w[l]
    end
    q = ones(m.n_src); r = ones(m.n_dst)
    for j in 1:m.n_src
        s_src[j] > 0 && (q[j] = m.wcov_src[j] / s_src[j])
    end
    for k in 1:m.n_dst
        s_dst[k] > 0 && (r[k] = m.wcov_dst[k] / s_dst[k])
    end
    amt = zeros(m.n_dst); cov = zeros(m.n_dst); wsum_v = zeros(m.n_dst)
    for l in eachindex(m.w)
        w = m.w[l]
        w > 0 || continue
        j = m.src[l]; k = m.dst[l]
        o = s.owner[j]
        cov[k] += w * o
        if kind == :extensive
            phi = m.wfull_src[j] > 0 ? q[j] * w / m.wfull_src[j] : 0.0
            amt[k] += phi * o * s.area_src[j] * s.v[j]
        else
            wsum_v[k] += w * o * s.v[j]
        end
    end
    for k in 1:m.n_dst
        if m.wfull_dst[k] > 0
            cov[k] *= r[k] / m.wfull_dst[k]
            kind == :intensive && (amt[k] = r[k] * wsum_v[k] / m.wfull_dst[k])
        end
        kind == :extensive && (amt[k] /= s.area_dst[k])
    end
    amt, cov
end

"""
    merge_offline(sources, background; kind = :extensive) -> (var, F)

The merge of `sources` (a vector of `OfflineSource`, any order: summed in the order given, the
caller sorts by name) onto the target grid with `background` (a vector of the target size).
"""
function merge_offline(sources::AbstractVector{OfflineSource}, background::AbstractVector{<:Real};
                       kind::Symbol = :extensive)
    kind in (:extensive, :intensive) || error("merge_offline: kind must be :extensive or :intensive")
    n = length(background)
    atot = zeros(n); ftot = zeros(n)
    for s in sources
        a, f = amounts(s, kind)
        atot .+= a
        ftot .+= f
    end
    var = if kind == :extensive
        (1 .- min.(ftot, 1)) .* background .+ atot
    else
        (1 .- min.(ftot, 1)) .* background .+ atot ./ max.(ftot, 1)
    end
    var, ftot
end

"""
    check_globe_ice(globe_nc, inputs_nc; tol_abs = 0.0, rtol = 1e-12) -> NamedTuple

The last record of `H_ice` and the static `cov` of `globe_ice.nc` against the offline merge of
the inputs file (see the module). `tol_abs`: the rounding of the file (0 for a lossless run).
Returns `(ok, max_err_H, max_err_cov, n_domains)`.
"""
function check_globe_ice(globe_nc::AbstractString, inputs_nc::AbstractString; tol_abs::Real = 0.0, rtol::Real = 1e-12)
    H_file, cov_file, nlon, nlat = NCDataset(globe_nc) do ds
        h = ds["H_ice"]
        nt = size(h, 3)
        Float64.(h[:, :, nt]), Float64.(ds["cov"][:, :]), size(h, 1), size(h, 2)
    end
    NCDataset(inputs_nc) do ds
        names = split(String(ds.attrib["domains"]), ",")
        abg = vec(Float64.(ds["area_bg"][:, :]))
        Hext = vec(Float64.(ds["H_ext"][:, :]))
        srcs = OfflineSource[]
        for d in sort(String.(names))
            m = read_scrip_con(String(ds.attrib["map_" * d]))
            own = haskey(ds, "owner_" * d) ? vec(Float64.(ds["owner_" * d][:, :])) :
                  ones(length(ds["H_" * d][:, :]))
            push!(srcs, OfflineSource(m, vec(Float64.(ds["area_" * d][:, :])), abg, vec(Float64.(ds["H_" * d][:, :])), own))
        end
        # the coupler maps each domain to the background grid; the merge target areas are the background's
        H_off, F = merge_offline(srcs, Hext; kind = :extensive)
        eh = maximum(abs.(vec(H_file) .- H_off))
        ec = maximum(abs.(vec(cov_file) .- min.(F, 1)))
        scale = max(maximum(abs.(H_off)), 1.0)
        ok = eh <= tol_abs + rtol * scale && ec <= max(tol_abs, 0.0) + 1e-9
        (ok = ok, max_err_H = eh, max_err_cov = ec, n_domains = length(srcs))
    end
end

end # module
