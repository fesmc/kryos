"""
    KryosTools.Compare

Compare two run directories bit for bit (regression-design.md §3): a yelmox
and a kryos run through a mapping table (`tools/gate/map_phase1.toml`), or two
runs of the same model (`map = "identity"`).

- `compare(ref, run; map, gate, deviations, restarts)` -> `Report`
- `compare_chunks(cont, chunk; gate, restarts)` -> `Report`: the restart round
  trip of one model (regression-design.md §1.2)
- `report_markdown(io, rep)`, `report_json(rep)`, `exit_code(rep)`: 0 pass,
  1 fail, 2 error (missing file or record, unreadable map or registry)
- `plot_variable(rep, var; dir)`: figures of one variable, defined by the
  CairoMakie extension (`using CairoMakie` first)

Classes (§3.3): E exact (type, shape and bits; `-0.0` vs `+0.0` is a
"signed zero" failure), C conversion-exact (`run` equals bit for bit the
declared expression of `ref`, with the declared `units`), U (`n` ulp) and T
(`atol + rtol*|ref|`) only through a registered deviation, A absent (listed),
N new (kryos-only, listed). The time axes are compared first: every record
time of a ref file must be in the run file, equal bit for bit.
"""
module Compare

using NCDatasets
using TOML
using JSON
using Printf
using Dates

export compare, compare_chunks, Report, VarResult, report_markdown, report_json, exit_code, plot_variable,
       load_map, load_gate, load_deviations

const MISSING_RECORD = "missing record"

# ---------------------------------------------------------------------------
# Inputs: the map, the gate, the deviations
# ---------------------------------------------------------------------------

"One row of the map (see the header of map_phase1.toml)."
Base.@kwdef struct Row
    ref::String
    run::String = ""
    vars::Vector{String} = ["*"]
    exclude::Vector{String} = String[]
    rename::Dict{String,String} = Dict{String,String}()
    run_var::String = ""
    class::String = "E"
    expr::String = ""
    units::String = ""
    comp::String = ""
end

struct MapSpec
    path::String
    rows::Vector{Row}
    restart_rows::Vector{Row}
    ref_dir::String
    run_dir::String
    constants::Dict{String,Float64}
    order::Vector{String}
    drift::Vector{String}
end

"""
    load_map(path) -> MapSpec

Read a mapping table; errors (exit 2) on unknown classes, missing fields or
an expression that does not parse.
"""
function load_map(path::AbstractString)
    d = TOML.parsefile(path)
    function row(r)
        haskey(r, "ref") || error("$path: a row without ref")
        cls = get(r, "class", "E")
        cls in ("E", "C", "A") || error("$path: row $(r["ref"]): class $cls (E, C or A)")
        cls != "A" && !haskey(r, "run") && error("$path: row $(r["ref"]): no run file")
        cls == "C" && !haskey(r, "expr") && error("$path: row $(r["ref"]): class C without expr")
        v = get(r, "vars", "*")
        vars = v isa AbstractString ? [String(v)] : String.(v)
        ex = get(r, "expr", "")
        isempty(ex) || parse_expr(ex)
        Row(ref=r["ref"], run=get(r, "run", ""), vars=vars, exclude=String.(get(r, "exclude", String[])),
            rename=Dict{String,String}(String(k) => String(v) for (k, v) in get(r, "rename", Dict())),
            run_var=get(r, "run_var", ""), class=cls, expr=ex, units=get(r, "units", ""),
            comp=get(r, "comp", ""))
    end
    rows = [row(r) for r in get(d, "row", [])]
    rst = get(d, "restart", Dict())
    rrows = [row(r) for r in get(rst, "row", [])]
    consts = Dict{String,Float64}(k => Float64(v) for (k, v) in get(d, "constants", Dict()))
    return MapSpec(String(path), rows, rrows, get(rst, "ref_dir", "restart-{kyr}-kyr"),
                   get(rst, "run_dir", "restart-{time}"), consts,
                   String.(get(d, "order", String[])), String.(get(get(d, "drift", Dict()), "vars", String[])))
end

"""
    load_gate(name_or_path; dir) -> Dict

The gate file (`tools/gate/<gate>.toml`): `gate`, `base`, `p`, `instances`
(placeholders), `chunks`, `restart_times`. A bare name is looked up in `dir`.
"""
function load_gate(g::AbstractString; dir::AbstractString=joinpath(@__DIR__, "..", "gate"))
    path = isfile(g) ? g : joinpath(dir, g * ".toml")
    isfile(path) || error("gate $g: no file $path")
    d = TOML.parsefile(path)
    haskey(d, "gate") || error("$path: no gate name")
    return d
end

Base.@kwdef struct Deviation
    id::String
    rows::Vector{String}
    class::String
    n::Int = 1
    atol::Float64 = 0.0
    rtol::Float64 = 0.0
    scope::String
    since::String
    changelog::String
    approved::String
    gates::Vector{String} = String[]
end

"""
    load_deviations(path; repo) -> Vector{Deviation}

Read the registry (§3.5) and check every entry: the fields, class U/T/skip,
scope local/rebaseline, and that the changelog file exists (relative to the
kryos repository `repo`, by default two levels above the registry).
"""
function load_deviations(path::AbstractString; repo::AbstractString=normpath(joinpath(dirname(abspath(path)), "..", "..")))
    d = TOML.parsefile(path)
    out = Deviation[]
    for e in get(d, "deviation", [])
        for k in ("id", "rows", "class", "scope", "since", "changelog", "approved")
            haskey(e, k) || error("$path: deviation $(get(e, "id", "?")): no $k")
        end
        e["class"] in ("U", "T", "skip") || error("$path: deviation $(e["id"]): class $(e["class"]) (U, T or skip)")
        e["scope"] in ("local", "rebaseline") || error("$path: deviation $(e["id"]): scope $(e["scope"])")
        e["class"] == "T" && !(haskey(e, "atol") || haskey(e, "rtol")) &&
            error("$path: deviation $(e["id"]): class T without atol, rtol")
        cl = split(e["changelog"], '#')[1]
        isfile(joinpath(repo, cl)) || error("$path: deviation $(e["id"]): changelog $cl does not exist in $repo")
        any(x -> x.id == e["id"], out) && error("$path: deviation $(e["id"]) twice")
        push!(out, Deviation(id=e["id"], rows=String.(e["rows"]), class=e["class"], n=get(e, "n", 1),
                             atol=Float64(get(e, "atol", 0.0)), rtol=Float64(get(e, "rtol", 0.0)),
                             scope=e["scope"], since=e["since"], changelog=e["changelog"],
                             approved=e["approved"], gates=String.(get(e, "gates", String[]))))
    end
    return out
end

# ---------------------------------------------------------------------------
# Expressions of class C (a small, safe evaluator over Julia syntax)
# ---------------------------------------------------------------------------

const OPS = Dict{Symbol,Function}(:+ => +, :- => -, :* => *, :/ => /)

function parse_expr(s::AbstractString)
    ex = Meta.parse(s)
    check(e) = e isa Symbol || e isa Number ||
               (e isa Expr && e.head == :call && (e.args[1] in keys(OPS) || e.args[1] in (:Float64, :Float32)) &&
                all(check, e.args[2:end])) ||
               error("expression '$s': only names, numbers, + - * / and Float64(), Float32()")
    check(ex)
    return ex
end

"Evaluate an expression elementwise; names are arrays (the ref variables) or constants."
function eval_expr(ex, env::Dict{Symbol,Any})
    ex isa Symbol && (haskey(env, ex) ? (return env[ex]) : error("expression: unknown name $ex"))
    ex isa Number && return ex
    f = ex.args[1]
    args = [eval_expr(a, env) for a in ex.args[2:end]]
    f == :Float64 && return Float64.(args[1])
    f == :Float32 && return Float32.(args[1])
    length(args) == 1 && f == :- && return (-).(args[1])
    return broadcast(OPS[f], args...)
end

# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

"The comparison of one run variable with its ref variable(s)."
Base.@kwdef mutable struct VarResult
    key::String                       # "<ref file>:<var>" (deviations match it)
    ref_file::String                  # relative to REF
    ref_vars::Vector{String}
    run_file::String = ""             # relative to RUN
    run_var::String = ""
    class::String                     # declared: E C A N
    applied::String = ""              # E C U T skip (after deviations)
    expr::String = ""
    units::String = ""
    comp::String = ""
    status::String = "pass"           # pass fail error absent new skipped drift
    note::String = ""
    deviation::String = ""
    ncells::Int = 0
    ndiff::Int = 0
    nsigned_zero::Int = 0
    nnan::Int = 0
    max_abs::Float64 = 0.0
    max_rel::Float64 = 0.0
    max_ulp::Int64 = 0
    first_rec::Int = 0                # aligned record (1-based), 0 = none / static
    first_time::Float64 = NaN
    loc_max::Vector{Int} = Int[]      # index of the maximum |difference| within its record
    time_max::Float64 = NaN
    diff_times::Vector{Float64} = Float64[]   # record times with a difference (-Inf: static)
    times::Vector{Float64} = Float64[]        # aligned record times
    maxdiff::Vector{Float64} = Float64[]      # max |difference| per aligned record
    # to reload the data (plots)
    ref_path::String = ""
    run_path::String = ""
    ref_recs::Vector{Int} = Int[]
    run_recs::Vector{Int} = Int[]
end

mutable struct Report
    ref::String
    run::String
    map::String
    gate::String
    mode::String                      # mapped | identity | drift
    results::Vector{VarResult}
    errors::Vector{String}            # exit 2
    unmapped::Vector{String}          # ref files and variables without a row: failures
    new_files::Vector{String}         # run files no row reads (N)
    deviations_used::Vector{String}
    deviations_unmatched::Vector{String}
    deviation_notes::Vector{String}
    first_divergence::Union{Nothing,Dict{String,Any}}
    order::Vector{String}
    constants::Dict{String,Float64}
end

Report(ref, run, map, mode) = Report(ref, run, map, "", mode, VarResult[], String[], String[], String[],
                                     String[], String[], String[], nothing, String[], Dict{String,Float64}())

"""
    exit_code(rep) -> Int

2 on an error, 1 on a failure (a variable, an unmapped yelmox file or variable,
a deviation that is not local), else 0. A drift report is 0 unless an error.
"""
function exit_code(rep::Report)
    isempty(rep.errors) || return 2
    rep.mode == "drift" && return 0
    (any(r -> r.status == "fail", rep.results) || !isempty(rep.unmapped)) && return 1
    return 0
end

# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

"Raw values of a variable (no scaling, no fill-value conversion), and its dimension names."
function read_raw(ds, name::AbstractString)
    v = ds[name]
    raw = v.var
    nd = ndims(raw)
    a = nd == 0 ? fill(raw[]) : raw[ntuple(_ -> Colon(), nd)...]
    return a, collect(String, dimnames(v))
end

"The `_FillValue` (else `missing_value`) attribute of a variable, or `nothing`."
function fill_value(v)
    at = v.attrib
    haskey(at, "_FillValue") && return at["_FillValue"]
    haskey(at, "missing_value") && return at["missing_value"]
    return nothing
end

"The raw time values of a file (empty without a time variable)."
function time_values(ds)
    haskey(ds, "time") || return Float64[]
    a, _ = read_raw(ds, "time")
    return Float64.(vec(a))
end

"Dimension coordinates of a file (1-D variables named as their dimension), time excluded."
coord_vars(ds) = [n for n in keys(ds) if n != "time" && haskey(ds.dim, n) &&
                  collect(dimnames(ds[n])) == [n]]

glob_regex(g::AbstractString) = Regex("^" * replace(replace(g, r"[.+^$()|\[\]{}\\]" => s"\\\0"), "*" => ".*", "?" => ".") * "\$")
globmatch(g::AbstractString, s::AbstractString) = occursin(glob_regex(g), s)

# ---------------------------------------------------------------------------
# Comparing arrays
# ---------------------------------------------------------------------------

ordered(x::Float32) = (i = Int64(reinterpret(Int32, x)); i >= 0 ? i : Int64(typemin(Int32)) - i)
ordered(x::Float64) = (i = Int128(reinterpret(Int64, x)); i >= 0 ? i : Int128(typemin(Int64)) - i)
ulps(a::T, b::T) where {T<:Union{Float32,Float64}} = abs(ordered(a) - ordered(b))

bitsequal(a::T, b::T) where {T<:Union{Float32,Float64}} = reinterpret(Unsigned, a) == reinterpret(Unsigned, b)
bitsequal(a, b) = isequal(a, b)

"""
Accumulate the differences of one record: `a` the expected (ref or its
conversion), `b` the run values, same type and shape. `tol` decides whether a
differing cell is accepted (U, T): returns (ndiff, nfail, nsz, nnan, maxabs, imax, maxrel, maxulp).
"""
function diff_record(a::AbstractArray, b::AbstractArray, tol)
    ndiff = 0; nfail = 0; nsz = 0; nnan = 0
    maxabs = 0.0; imax = 0; maxrel = 0.0; maxulp = Int128(0)
    isfl = eltype(a) <: Union{Float32,Float64}
    @inbounds for i in eachindex(a, b)
        x = a[i]; y = b[i]
        bitsequal(x, y) && continue
        ndiff += 1
        if isfl
            if isnan(x) || isnan(y)
                nnan += 1
                nfail += 1
                continue
            end
            x == y && (nsz += 1)                 # +0.0 vs -0.0
            d = abs(Float64(y) - Float64(x))
            if d > maxabs || imax == 0
                maxabs = d; imax = i
            end
            x != 0 && (maxrel = max(maxrel, d / abs(Float64(x))))
            u = ulps(x, y)
            maxulp = max(maxulp, u)
            tol(x, y, u) || (nfail += 1)
        else
            imax == 0 && (imax = i)
            nfail += 1
        end
    end
    return ndiff, nfail, nsz, nnan, maxabs, imax, maxrel, maxulp
end

"The acceptance of a differing cell under the applied class."
function tolerance(class::AbstractString, dev)
    class == "U" && return (x, y, u) -> u <= dev.n
    class == "T" && return (x, y, u) -> abs(Float64(y) - Float64(x)) <= dev.atol + dev.rtol * abs(Float64(x))
    return (x, y, u) -> false
end

# ---------------------------------------------------------------------------
# One comparison: ref variable(s) -> run variable
# ---------------------------------------------------------------------------

"""
Compare `r.ref_vars` of the ref dataset with `r.run_var` of the run dataset
through the aligned records (`rmap`: ref record => run record). Fills `r`.
"""
function compare_var!(r::VarResult, dsr, dsu, rmap, consts::Dict{String,Float64}, dev)
    tol = tolerance(r.applied, dev)
    # the expected values
    arrays = Dict{Symbol,Any}()
    dims_ref = String[]
    for v in r.ref_vars
        haskey(dsr, v) || (r.status = "error"; r.note = "no variable $v in the ref file"; return r)
        a, dn = read_raw(dsr, v)
        arrays[Symbol(v)] = a
        isempty(dims_ref) && (dims_ref = dn)
        dn == dims_ref || (r.status = "fail"; r.note = "ref variables $(join(r.ref_vars, ", ")) on different dimensions"; return r)
    end
    if !haskey(dsu, r.run_var)
        r.status = "fail"; r.note = "no variable $(r.run_var) in $(r.run_file)"
        return r
    end
    b, dims_run = read_raw(dsu, r.run_var)
    if r.class == "C"
        env = Dict{Symbol,Any}(Symbol(k) => v for (k, v) in consts)
        merge!(env, arrays)
        length(r.ref_vars) == 1 && (env[:x] = arrays[Symbol(r.ref_vars[1])])
        a = eval_expr(parse_expr(r.expr), env)
        a isa AbstractArray || (a = fill(a, size(arrays[Symbol(r.ref_vars[1])])))
        # missing values: a ref fill value is the run's fill value, not converted
        for v in r.ref_vars
            fv = fill_value(dsr[v])
            fv === nothing && continue
            x = arrays[Symbol(v)]
            mask = [bitsequal(e, convert(eltype(x), fv)) for e in x]
            any(mask) || continue
            ufv = fill_value(dsu[r.run_var])
            if ufv === nothing
                r.note = "the ref has fill values, the run variable no _FillValue"
            else
                a[mask] .= convert(eltype(a), ufv)
            end
        end
        if !isempty(r.units)
            u = get(dsu[r.run_var].attrib, "units", "")
            u == r.units || (r.status = "fail"; r.note = "units \"$u\", expected \"$(r.units)\"")
        end
    else
        a = arrays[Symbol(r.ref_vars[1])]
    end
    if eltype(a) != eltype(b)
        r.status = "fail"; r.note = "type $(eltype(b)), expected $(eltype(a))"
        return r
    end
    # the record dimension and the other dimensions
    ta = findfirst(==("time"), dims_ref)
    tb = findfirst(==("time"), dims_run)
    (ta === nothing) != (tb === nothing) &&
        (r.status = "fail"; r.note = "time dimension in one file only"; return r)
    other(dn, t) = t === nothing ? dn : dn[setdiff(1:length(dn), t)]
    sa = t -> t === nothing ? size(a) : size(a)[setdiff(1:ndims(a), t)]
    if other(dims_ref, ta) != other(dims_run, tb) || sa(ta) != (tb === nothing ? size(b) : size(b)[setdiff(1:ndims(b), tb)])
        r.status = "fail"
        r.note = "shape $(size(b)) $(dims_run), expected $(size(a)) $(dims_ref)"
        return r
    end
    recs = ta === nothing ? [(0, 0, -Inf)] : rmap
    nfail_total = 0
    for (ia, ib, t) in recs
        x = ia == 0 ? a : selectdim(a, ta, ia)
        y = ib == 0 ? b : selectdim(b, tb, ib)
        ndiff, nfail, nsz, nnan, maxabs, imax, maxrel, maxulp = diff_record(x, y, tol)
        r.ncells += length(x)
        push!(r.times, t)
        push!(r.maxdiff, maxabs)
        ia == 0 || (push!(r.ref_recs, ia); push!(r.run_recs, ib))
        ndiff == 0 && continue
        r.ndiff += ndiff
        r.nsigned_zero += nsz
        r.nnan += nnan
        nfail_total += nfail
        push!(r.diff_times, t)
        if r.first_rec == 0 && ia != 0 || (ia == 0 && isnan(r.first_time))
            r.first_rec = ia; r.first_time = t
        end
        if maxabs > r.max_abs || isempty(r.loc_max)
            r.max_abs = maxabs; r.time_max = t
            r.loc_max = imax == 0 ? Int[] : collect(Tuple(CartesianIndices(x)[imax]))
        end
        r.max_rel = max(r.max_rel, maxrel)
        r.max_ulp = max(r.max_ulp, Int64(min(maxulp, Int128(typemax(Int64)))))
    end
    if r.status == "pass"
        if r.ndiff > 0 && (nfail_total > 0 || r.applied in ("E", "C"))
            r.status = "fail"
            r.nsigned_zero > 0 && (r.note = "signed zero: $(r.nsigned_zero) cells")
            r.nnan > 0 && (r.note = strip(r.note * " NaN mismatch: $(r.nnan) cells"))
        end
    end
    return r
end

# ---------------------------------------------------------------------------
# File pairs and rows
# ---------------------------------------------------------------------------

"Placeholders `{name}` of a pattern: instances, else captures."
function pattern_regex(pat::AbstractString, inst::Dict{String,String})
    out = IOBuffer()
    caps = String[]
    last = 1
    for m in eachmatch(r"\{(\w+)\}", pat)
        print(out, escape_re(pat[last:prevind(pat, m.offset)]))
        n = m.captures[1]
        if haskey(inst, n)
            print(out, escape_re(inst[n]))
        else
            n in caps ? print(out, "\\g{$n}") : print(out, "(?<$n>[^/]+?)")
            push!(caps, n)
        end
        last = m.offset + length(m.match)
    end
    print(out, escape_re(pat[last:end]))
    return Regex("^" * String(take!(out)) * "\$"), caps
end
escape_re(s::AbstractString) = replace(s, r"[.+*?^$()|\[\]{}\\]" => s"\\\0")

function substitute(pat::AbstractString, vals::Dict{String,String})
    s = pat
    for m in eachmatch(r"\{(\w+)\}", pat)
        haskey(vals, m.captures[1]) || error("pattern $pat: no value for {$(m.captures[1])} (gate instances?)")
        s = replace(s, m.match => vals[m.captures[1]])
    end
    return s
end

"""
The time axis check of a file pair: (ref record, run record, time) or the
error. With `window = (lo, hi)` only the ref records with `lo < t <= hi` are
aligned.
"""
function align_records(dsr, dsu, ref_file, run_file; window=nothing)
    tr = time_values(dsr)
    isempty(tr) && return Tuple{Int,Int,Float64}[], nothing
    tu = time_values(dsu)
    out = Tuple{Int,Int,Float64}[]
    for (i, t) in enumerate(tr)
        window === nothing || window[1] < t <= window[2] || continue
        j = findfirst(u -> bitsequal(u, t), tu)
        j === nothing && return out, "$MISSING_RECORD: time $t of $ref_file (record $i) is not in $run_file"
        push!(out, (i, j, t))
    end
    return out, nothing
end

"Open both files and run the comparisons of one file pair."
function compare_pair!(rep::Report, results::Vector{VarResult}, refdir, rundir, ref_rel, run_rel, devs, mode;
                       window=nothing)
    rp = joinpath(refdir, ref_rel); up = joinpath(rundir, run_rel)
    if !isfile(up)
        push!(rep.errors, "missing file: $run_rel (for $ref_rel)")
        foreach(r -> (r.status = "error"; r.note = "missing run file"), results)
        append!(rep.results, results)
        return
    end
    NCDataset(rp) do dsr
        NCDataset(up) do dsu
            rmap, err = align_records(dsr, dsu, ref_rel, run_rel; window=window)
            if err !== nothing
                push!(rep.errors, err)
                foreach(r -> (r.status = "error"; r.note = "time axis"), results)
                append!(rep.results, results)
                return
            end
            for r in results
                r.ref_path = rp; r.run_path = up
                dev = apply_deviation!(r, rep, devs, mode)
                if r.class == "A"
                    r.status = "absent"
                elseif r.applied == "skip"
                    r.status = "skipped"
                else
                    try
                        compare_var!(r, dsr, dsu, rmap, rep.constants, dev)
                    catch e
                        r.status = "error"; r.note = sprint(showerror, e)
                        push!(rep.errors, "$(r.key): $(r.note)")
                    end
                end
                mode == "drift" && r.status in ("pass", "fail") && (r.status = "drift")
                push!(rep.results, r)
            end
        end
    end
end

"The deviation of a row (the first whose rows match), recorded in the result."
function apply_deviation!(r::VarResult, rep::Report, devs, mode)
    r.applied = r.class
    mode == "mapped" || return nothing
    for d in devs
        d.scope == "local" || continue
        (isempty(d.gates) || rep.gate in d.gates) || continue
        if any(p -> globmatch(p, r.key) || any(v -> globmatch(p, "$(r.ref_file):$v"), r.ref_vars), d.rows)
            r.deviation = d.id
            r.applied = d.class
            d.id in rep.deviations_used || push!(rep.deviations_used, d.id)
            return d
        end
    end
    return nothing
end

"Variables of a ref file a row compares, as (ref vars, run var) pairs."
function row_vars(row::Row, refvars::Vector{String}, coords::Vector{String}, named_same::Set{String})
    data = [v for v in refvars if v != "time" && !(v in coords)]
    if row.vars == ["*"]
        sel = [v for v in data if !(v in named_same) && !any(g -> globmatch(g, v), row.exclude)]
    else
        sel = [v for v in data if any(g -> globmatch(g, v), row.vars)]
    end
    isempty(row.run_var) || return isempty(sel) ? Tuple{Vector{String},String}[] : [(sel, row.run_var)]
    return [([v], get(row.rename, v, v)) for v in sel]
end

function list_nc(dir::AbstractString)
    isdir(dir) || return String[]
    return sort([f for f in readdir(dir) if endswith(f, ".nc") && isfile(joinpath(dir, f))])
end

function list_nc_rec(dir::AbstractString)
    out = String[]
    for (root, _, files) in walkdir(dir), f in files
        endswith(f, ".nc") && push!(out, relpath(joinpath(root, f), dir))
    end
    return sort(out)
end

"""
Run the rows of a map over the files of `refdir` (relative names), the run
files in `rundir`. `prefix` is prepended to the reported file names (restart
bundles).
"""
function compare_rows!(rep::Report, rows::Vector{Row}, refdir, rundir, inst, devs, mode; prefix="")
    reffiles = list_nc(refdir)
    used_run = Set{String}()
    covered = Dict{String,Set{String}}()           # ref file => variables matched by a row
    for f in reffiles
        matched = false
        NCDataset(joinpath(refdir, f)) do ds
            refvars = collect(String, keys(ds))
            coords = coord_vars(ds)
            cov = get!(covered, f, Set{String}(vcat(coords, ["time"])))
            # rows of this file, grouped by run file
            pairs = Dict{String,Vector{VarResult}}()
            porder = String[]
            for row in rows
                re, caps = pattern_regex(row.ref, inst)
                m = match(re, f)
                m === nothing && continue
                matched = true
                vals = copy(inst)
                for c in caps
                    vals[c] = String(m[Symbol(c)])
                end
                run_rel = row.class == "A" ? "" : substitute(row.run, vals)
                # "*": the variables no other row of this ref and run file names
                named_same = Set{String}()
                for o in rows
                    (o === row || o.vars == ["*"]) && continue
                    ore, ocaps = pattern_regex(o.ref, inst)
                    om = match(ore, f)
                    om === nothing && continue
                    if o.class != "A"
                        ovals = copy(inst)
                        foreach(c -> ovals[c] = String(om[Symbol(c)]), ocaps)
                        substitute(o.run, ovals) == run_rel || continue
                    end
                    for v in refvars
                        any(g -> globmatch(g, v), o.vars) && push!(named_same, v)
                    end
                end
                for (rv, uv) in row_vars(row, refvars, coords, named_same)
                    union!(cov, rv)
                    r = VarResult(key="$(prefix)$(f):$(join(rv, ","))", ref_file=prefix * f, ref_vars=rv,
                                  run_file=prefix * run_rel, run_var=uv, class=row.class, expr=row.expr,
                                  units=row.units, comp=row.comp)
                    push!(get!(() -> (push!(porder, run_rel); VarResult[]), pairs, run_rel), r)
                end
                if row.class != "A" && !haskey(pairs, run_rel)
                    pairs[run_rel] = VarResult[]; push!(porder, run_rel)
                end
            end
            matched || (push!(rep.unmapped, "$(prefix)$(f): no row"); return)
            for v in refvars
                v in cov || push!(rep.unmapped, "$(prefix)$(f):$v: no row")
            end
            for run_rel in porder
                rs = pairs[run_rel]
                if isempty(run_rel)               # A rows
                    for r in rs
                        r.status = "absent"; r.applied = "A"; push!(rep.results, r)
                    end
                    continue
                end
                push!(used_run, run_rel)
                # coordinates with every file pair
                if isfile(joinpath(rundir, run_rel))
                    cs = NCDataset(ds2 -> coord_vars(ds2), joinpath(rundir, run_rel))
                    for c in coords
                        c in cs || continue
                        any(r -> r.ref_vars == [c], rs) && continue
                        push!(rs, VarResult(key="$(prefix)$(f):$c", ref_file=prefix * f, ref_vars=[c],
                                            run_file=prefix * run_rel, run_var=c, class="E",
                                            comp=isempty(rs) ? "" : rs[1].comp))
                    end
                end
                compare_pair!(rep, rs, refdir, rundir, f, run_rel, devs, mode)
                foreach(r -> (r.ref_file = prefix * f; r.run_file = prefix * run_rel), rs)
            end
        end
    end
    return used_run
end

"N: run variables of compared run files that no comparison reads, and run files no row reads."
function list_new!(rep::Report, rundir, used_run; prefix="")
    for f in list_nc(rundir)
        if !(f in used_run)
            push!(rep.new_files, prefix * f)
            continue
        end
        read_vars = Set(r.run_var for r in rep.results if r.run_file == prefix * f)
        NCDataset(joinpath(rundir, f)) do ds
            for v in keys(ds)
                (v in read_vars || v == "time") && continue
                push!(rep.results, VarResult(key="$(prefix)$(f):$v", ref_file="", ref_vars=String[],
                                             run_file=prefix * f, run_var=v, class="N", applied="N",
                                             status="new"))
            end
        end
    end
end

# ---------------------------------------------------------------------------
# Restart bundles
# ---------------------------------------------------------------------------

"kryos bundle time: shortest exact decimal of the time in years."
time_label(t::Real) = isinteger(t) ? string(Int(t)) : string(Float64(t))
"yelmox bundle time: time/1000 with three decimals."
kyr_label(t::Real) = @sprintf("%.3f", t / 1000)

bundle_regex(pat) = Regex("^" * replace(escape_re(pat), "\\{kyr\\}" => "(?<kyr>[-0-9.]+)", "\\{time\\}" => "(?<time>[-0-9.e]+)") * "\$")

"Time of the bundle folder `d` named by `pat`, or `nothing`."
function bundle_time(pat, d)
    m = match(bundle_regex(pat), d)
    m === nothing && return nothing
    return m[:kyr] !== nothing ? parse(Float64, m[:kyr]) * 1000 : parse(Float64, m[:time])
end

bundle_times(refdir, pat) =
    sort([t for d in readdir(refdir) if isdir(joinpath(refdir, d)) for t in (bundle_time(pat, d),) if t !== nothing])

"Restart bundle folders of a run directory (either model)."
restart_dirs(dir) = sort([d for d in readdir(dir) if startswith(d, "restart-") && isdir(joinpath(dir, d))])

"Bundle names of yelmox and kryos (the defaults of the map's `[restart]`)."
const BUNDLE_PATTERNS = ("restart-{kyr}-kyr", "restart-{time}")

bundle_name(pat, t) = replace(pat, "{kyr}" => kyr_label(t), "{time}" => time_label(t))

# ---------------------------------------------------------------------------
# The comparison
# ---------------------------------------------------------------------------

const YELMOX_COMPONENTS = [r"^yelmo" => "ice", r"^isos" => "geo", r"^bsl" => "geo", r"^htopo" => "hub",
                           r"^snap" => "climate", r"^esm" => "climate", r"^rembo" => "climate",
                           r"^chion" => "surface", r"^smbpal" => "surface", r"^mshlf" => "shelf",
                           r"^marine_shelf" => "shelf", r"^coef_" => "forcing"]

"Component of a file in identity mode: an instance prefix of the gate, else the yelmox name."
function component_of(file::AbstractString, inst::Dict{String,String})
    base = basename(file)
    best = ""; comp = ""
    for (k, v) in inst
        k == "domain" && continue
        if startswith(base, v) && length(v) > length(best)
            best = v; comp = k == "geo_glob" ? "geo" : k
        end
    end
    isempty(comp) || return comp
    startswith(base, get(inst, "domain", "\0")) && return "hub"
    for (re, c) in YELMOX_COMPONENTS
        occursin(re, base) && return c
    end
    return ""
end

"""
    compare(ref, run; map="identity", gate="", deviations="", restarts=false,
            restart_times=nothing) -> Report

Compare the run directory `run` against `ref`. `map` is a mapping table path
or `"identity"` (every `*.nc` of `ref` against the same file of `run`, class
E). `gate` (a name of `tools/gate` or a path) gives the instance placeholders
and the restart times; `deviations` the registry (mapped mode only).
`restarts` also compares the restart bundles.
"""
function compare(refdir::AbstractString, rundir::AbstractString; map::AbstractString="identity",
                 gate::AbstractString="", deviations::AbstractString="", restarts::Bool=false,
                 restart_times=nothing, repo::Union{Nothing,AbstractString}=nothing)
    rep = Report(refdir, rundir, map, map == "identity" ? "identity" : "mapped")
    rep.gate = gate
    inst = Dict{String,String}()
    gd = Dict{String,Any}()
    devs = Deviation[]
    spec = nothing
    try
        isdir(refdir) || error("no directory $refdir")
        isdir(rundir) || error("no directory $rundir")
        if !isempty(gate)
            gd = load_gate(gate)
            rep.gate = gd["gate"]
            inst = Dict{String,String}(k => String(v) for (k, v) in get(gd, "instances", Dict()))
        end
        restart_times === nothing && (restart_times = Float64.(get(gd, "restart_times", Float64[])))
        if map != "identity"
            spec = load_map(map)
            rep.order = spec.order
            rep.constants = spec.constants
            if !isempty(deviations)
                devs = repo === nothing ? load_deviations(deviations) : load_deviations(deviations; repo=repo)
            end
            if any(d -> d.scope == "rebaseline" && (isempty(d.gates) || rep.gate in d.gates), devs)
                rep.mode = "drift"
                ids = [d.id for d in devs if d.scope == "rebaseline" && (isempty(d.gates) || rep.gate in d.gates)]
                append!(rep.deviations_used, ids)
                push!(rep.deviation_notes, "rebaselined ($(join(ids, ", "))): drift report of $(join(spec.drift, ", ")), no pass/fail")
            end
        else
            rep.order = copy(IDENTITY_ORDER)
            isempty(deviations) || push!(rep.deviation_notes, "identity mode: deviations not applied")
        end
    catch e
        push!(rep.errors, sprint(showerror, e))
        return rep
    end

    try
        if rep.mode == "identity"
            compare_identity!(rep, refdir, rundir, inst)
            restarts && compare_identity_restarts!(rep, refdir, rundir, inst, restart_dirs(refdir))
        else
            rows = spec.rows
            if rep.mode == "drift"
                rows = drift_rows(spec)
            end
            used = compare_rows!(rep, rows, refdir, rundir, inst, devs, rep.mode)
            rep.mode == "drift" && empty!(rep.unmapped)
            rep.mode == "mapped" && list_new!(rep, rundir, used)
            if restarts && rep.mode == "mapped"
                compare_restarts!(rep, spec, refdir, rundir, inst, devs, restart_times)
            end
        end
    catch e
        push!(rep.errors, sprint(showerror, e))
    end

    # local deviations: the other rows must still pass; unmatched ones are listed
    for d in devs
        d.scope == "local" || continue
        (isempty(d.gates) || rep.gate in d.gates) || continue
        if !(d.id in rep.deviations_used)
            push!(rep.deviations_unmatched, d.id)
        elseif any(r -> r.status == "fail" && r.deviation != d.id, rep.results)
            push!(rep.deviation_notes, "deviation $(d.id) is not local: rows outside it fail")
        end
    end
    rep.first_divergence = first_divergence(rep)
    return rep
end

"Chunk folders of a chunk run (yelmox `scripts/run_chunks.sh`)."
chunk_dir(k::Integer) = @sprintf("c%02d", k)

"""
    compare_chunks(cont, chunk; gate, restarts=true) -> Report

The restart round trip of one model (regression-design.md §1.2): the run
`chunk` is split at the gate's `chunks = [t0, t1, ..., tn]` into the folders
`c00`, `c01`, ...; chunk `k` (0-based) covers `(t_k, t_k+1]`, the first also
`t0`. Every record of `cont` in that interval must be in chunk `k` and equal
(class E), and with `restarts` the bundle of `cont` at `t_k+1` must equal the
one chunk `k` writes at `t_k+1`. The record at `t_k` and the bundle a chunk
writes at its own start are diagnostics of a state not yet stepped (not every
field is restored) and are not compared.
"""
function compare_chunks(contdir::AbstractString, chunkdir::AbstractString; gate::AbstractString,
                        restarts::Bool=true)
    rep = Report(contdir, chunkdir, "identity", "identity")
    rep.order = copy(IDENTITY_ORDER)
    inst = Dict{String,String}()
    ts = Float64[]
    try
        isdir(contdir) || error("no directory $contdir")
        isdir(chunkdir) || error("no directory $chunkdir")
        gd = load_gate(gate)
        rep.gate = gd["gate"]
        inst = Dict{String,String}(k => String(v) for (k, v) in get(gd, "instances", Dict()))
        ts = Float64.(get(gd, "chunks", Float64[]))
        length(ts) >= 3 || error("gate $(rep.gate): chunks $(ts) (at least two chunks)")
    catch e
        push!(rep.errors, sprint(showerror, e))
        return rep
    end
    try
        for k in 0:length(ts)-2
            cdir = joinpath(chunkdir, chunk_dir(k))
            isdir(cdir) || (push!(rep.errors, "missing chunk folder $(chunk_dir(k))"); continue)
            lo = k == 0 ? -Inf : ts[k+1]
            compare_identity!(rep, contdir, cdir, inst; prefix="$(chunk_dir(k))/", window=(lo, ts[k+2]))
            restarts || continue
            ds = [d for d in restart_dirs(contdir) if any(p -> bundle_time(p, d) == ts[k+2], BUNDLE_PATTERNS)]
            isempty(ds) && (push!(rep.errors, "no restart bundle at $(ts[k+2]) in $contdir"); continue)
            compare_identity_restarts!(rep, contdir, cdir, inst, ds; prefix="$(chunk_dir(k))/")
        end
    catch e
        push!(rep.errors, sprint(showerror, e))
    end
    rep.first_divergence = first_divergence(rep)
    return rep
end

function drift_rows(spec::MapSpec)
    rows = Row[]
    for dv in spec.drift
        f, v = split(dv, ':', limit=2)
        for row in spec.rows
            row.class == "A" && continue
            row.ref == f || continue
            (row.vars == ["*"] || any(g -> globmatch(g, v), row.vars)) || continue
            push!(rows, Row(ref=row.ref, run=row.run, vars=[String(v)], rename=row.rename, class=row.class,
                            expr=row.expr, units=row.units, comp=row.comp))
            break
        end
    end
    return rows
end

const IDENTITY_ORDER = ["geo", "ice", "hub", "climate", "surface", "shelf", "forcing"]

function compare_identity!(rep::Report, refdir, rundir, inst; prefix="", window=nothing)
    files = list_nc(refdir)
    for f in files
        NCDataset(joinpath(refdir, f)) do ds
            vs = [v for v in keys(ds) if v != "time"]
            comp = component_of(f, inst)
            rs = [VarResult(key="$(prefix)$(f):$v", ref_file=prefix * f, ref_vars=[v], run_file=prefix * f,
                            run_var=v, class="E", comp=comp) for v in vs]
            compare_pair!(rep, rs, refdir, rundir, f, f, Deviation[], "identity"; window=window)
        end
    end
    for f in list_nc(rundir)
        f in files || push!(rep.new_files, prefix * f)
    end
    # run variables the ref file does not have
    for f in files
        isfile(joinpath(rundir, f)) || continue
        refv = NCDataset(ds -> Set(keys(ds)), joinpath(refdir, f))
        NCDataset(joinpath(rundir, f)) do ds
            for v in keys(ds)
                v in refv && continue
                push!(rep.results, VarResult(key="$(prefix)$(f):$v", ref_file="", ref_vars=String[],
                                             run_file=prefix * f, run_var=v, class="N", applied="N", status="new"))
            end
        end
    end
end

"The bundle folders `dirs` of `refdir` against the same folders of `rundir`."
function compare_identity_restarts!(rep::Report, refdir, rundir, inst, dirs; prefix="")
    for d in dirs
        isdir(joinpath(rundir, d)) || (push!(rep.errors, "missing restart bundle: $(prefix)$d"); continue)
        for f in list_nc_rec(joinpath(refdir, d))
            compare_identity_file!(rep, joinpath(refdir, d), joinpath(rundir, d), f, inst, "$(prefix)$d/")
        end
    end
end

function compare_identity_file!(rep::Report, refdir, rundir, f, inst, prefix)
    NCDataset(joinpath(refdir, f)) do ds
        vs = [v for v in keys(ds) if v != "time"]
        rs = [VarResult(key="$(prefix)$(f):$v", ref_file=prefix * f, ref_vars=[v], run_file=prefix * f,
                        run_var=v, class="E", comp=component_of(f, inst)) for v in vs]
        compare_pair!(rep, rs, refdir, rundir, f, f, Deviation[], "identity")
        foreach(r -> (r.ref_file = prefix * f; r.run_file = prefix * f), rs)
    end
end

function compare_restarts!(rep::Report, spec::MapSpec, refdir, rundir, inst, devs, times)
    ts = isempty(times) ? bundle_times(refdir, spec.ref_dir) : times
    isempty(ts) && push!(rep.errors, "--restarts: no restart bundle in $refdir")
    for t in ts
        rd = bundle_name(spec.ref_dir, t)
        ud = bundle_name(spec.run_dir, t)
        isdir(joinpath(refdir, rd)) || (push!(rep.errors, "missing restart bundle: $rd (ref)"); continue)
        isdir(joinpath(rundir, ud)) || (push!(rep.errors, "missing restart bundle: $ud (run)"); continue)
        n0 = length(rep.results)
        compare_rows!(rep, spec.restart_rows, joinpath(refdir, rd), joinpath(rundir, ud), inst, devs,
                      rep.mode; prefix="$rd/")
        # run file names relative to the run bundle
        for r in rep.results[n0+1:end]
            startswith(r.run_file, "$rd/") && (r.run_file = "$ud/" * r.run_file[length(rd)+2:end])
            isfinite(r.first_time) || (r.ndiff > 0 && (r.first_time = t))
            replace!(r.diff_times, -Inf => t)
        end
    end
end

"""
The earliest time at which a compared variable fails, and every variable that
fails at that time, in the coupling order (`order` of the map).
"""
function first_divergence(rep::Report)
    fails = [r for r in rep.results if r.status in ("fail", "drift") && !isempty(r.diff_times)]
    isempty(fails) && return nothing
    t = minimum(minimum(r.diff_times) for r in fails)
    at = [r for r in fails if t in r.diff_times]
    rank(r) = (c = findfirst(==(r.comp), rep.order); c === nothing ? length(rep.order) + 1 : c)
    sort!(at, by=r -> (rank(r), r.key))
    return Dict{String,Any}("time" => t, "static" => t == -Inf,
                            "variables" => [Dict("key" => r.key, "comp" => r.comp, "run" => "$(r.run_file):$(r.run_var)",
                                                 "max_abs" => r.max_abs) for r in at])
end

# ---------------------------------------------------------------------------
# Reports
# ---------------------------------------------------------------------------

outcome(rep::Report) = (c = exit_code(rep); c == 2 ? "ERROR" : rep.mode == "drift" ? "DRIFT REPORT" : c == 1 ? "FAIL" : "PASS")

fmt(x::Real) = isnan(x) ? "" : x == -Inf ? "static" : x == 0 ? "0" : @sprintf("%.4g", x)

"""
    report_markdown(io, rep)

The report as Markdown: outcome, first divergence, errors, unmapped yelmox
files and variables, per variable (class, status, cells differing, max abs,
max rel, max ulp, first differing record, location of the maximum), A and N
lists, deviations.
"""
function report_markdown(io::IO, rep::Report)
    println(io, "# kryos-compare: ", outcome(rep), " (exit ", exit_code(rep), ")\n")
    println(io, "- ref: `", rep.ref, "`")
    println(io, "- run: `", rep.run, "`")
    println(io, "- map: `", rep.map, "` (", rep.mode, ")", isempty(rep.gate) ? "" : ", gate `$(rep.gate)`")
    println(io, "- date: ", Dates.format(now(), "yyyy-mm-dd HH:MM"))
    cmp = [r for r in rep.results if !(r.status in ("absent", "new"))]
    nf = count(r -> r.status == "fail", cmp)
    println(io, "- compared: ", length(cmp), " variables, ", nf, " failing, ",
            count(r -> r.status == "absent", rep.results), " absent (A), ",
            count(r -> r.status == "new", rep.results), " new (N)\n")
    fd = rep.first_divergence
    println(io, "## First divergence\n")
    if fd === nothing
        println(io, "None.\n")
    else
        println(io, fd["static"] ? "Static fields" : "Time $(fmt(fd["time"]))", ", in coupling order:\n")
        for v in fd["variables"]
            println(io, "- `", v["key"], "` -> `", v["run"], "` (", v["comp"], "), max |diff| ", fmt(v["max_abs"]))
        end
        println(io)
    end
    if !isempty(rep.errors)
        println(io, "## Errors\n")
        foreach(e -> println(io, "- ", e), rep.errors)
        println(io)
    end
    if !isempty(rep.unmapped)
        println(io, "## Without a row (failures)\n")
        foreach(e -> println(io, "- `", e, "`"), rep.unmapped)
        println(io)
    end
    println(io, "## Variables\n")
    println(io, "| ref | run | class | status | cells | differing | max abs | max rel | max ulp | first diff (t, rec) | loc max | note |")
    println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in sort(cmp, by=r -> (r.status == "pass", r.key))
        cls = r.applied == r.class || isempty(r.applied) ? r.class : "$(r.class)->$(r.applied) ($(r.deviation))"
        fr = r.ndiff == 0 ? "" : r.first_rec == 0 ? "static" : "$(fmt(r.first_time)), $(r.first_rec)"
        println(io, "| `", r.key, "` | `", r.run_file, ":", r.run_var, "` | ", cls, " | ", r.status, " | ",
                r.ncells, " | ", r.ndiff, " | ", fmt(r.max_abs), " | ", fmt(r.max_rel), " | ",
                r.ndiff == 0 ? "" : r.max_ulp, " | ", fr, " | ", isempty(r.loc_max) ? "" : join(r.loc_max, ","),
                " | ", r.note, " |")
    end
    println(io)
    ab = [r for r in rep.results if r.status == "absent"]
    if !isempty(ab)
        println(io, "## Absent (A)\n")
        println(io, join(["`$(r.key)`" for r in ab], ", "), "\n")
    end
    nw = [r for r in rep.results if r.status == "new"]
    if !isempty(nw) || !isempty(rep.new_files)
        println(io, "## New (N)\n")
        isempty(rep.new_files) || println(io, "Files: ", join(["`$f`" for f in rep.new_files], ", "), "\n")
        isempty(nw) || println(io, "Variables: ", join(["`$(r.run_file):$(r.run_var)`" for r in nw], ", "), "\n")
    end
    println(io, "## Deviations\n")
    println(io, "- used: ", isempty(rep.deviations_used) ? "none" : join(rep.deviations_used, ", "))
    println(io, "- matched nothing (remove them): ", isempty(rep.deviations_unmatched) ? "none" : join(rep.deviations_unmatched, ", "))
    foreach(n -> println(io, "- ", n), rep.deviation_notes)
    return nothing
end
report_markdown(rep::Report) = sprint(report_markdown, rep)

jnum(x::Real) = isfinite(x) ? x : x == -Inf ? "static" : nothing

"""
    report_json(rep) -> String

The report as JSON (the same content as the Markdown, every result with its
per-record maximum difference).
"""
function report_json(rep::Report)
    res = [Dict("key" => r.key, "ref_file" => r.ref_file, "ref_vars" => r.ref_vars, "run_file" => r.run_file,
                "run_var" => r.run_var, "class" => r.class, "applied" => r.applied, "deviation" => r.deviation,
                "expr" => r.expr, "units" => r.units, "comp" => r.comp, "status" => r.status, "note" => r.note,
                "cells" => r.ncells, "differing" => r.ndiff, "signed_zero" => r.nsigned_zero, "nan" => r.nnan,
                "max_abs" => jnum(r.max_abs), "max_rel" => jnum(r.max_rel), "max_ulp" => r.max_ulp,
                "first_record" => r.first_rec, "first_time" => jnum(r.first_time), "loc_max" => r.loc_max,
                "time_max" => jnum(r.time_max), "times" => jnum.(r.times), "maxdiff" => r.maxdiff)
           for r in rep.results]
    fd = rep.first_divergence
    d = Dict("outcome" => outcome(rep), "exit" => exit_code(rep), "ref" => rep.ref, "run" => rep.run,
             "map" => rep.map, "mode" => rep.mode, "gate" => rep.gate, "errors" => rep.errors,
             "unmapped" => rep.unmapped, "new_files" => rep.new_files,
             "first_divergence" => fd === nothing ? nothing :
                 Dict("time" => jnum(fd["time"]), "variables" => fd["variables"]),
             "deviations" => Dict("used" => rep.deviations_used, "unmatched" => rep.deviations_unmatched,
                                  "notes" => rep.deviation_notes),
             "results" => res)
    return JSON.json(d)
end

# ---------------------------------------------------------------------------
# Plots (CairoMakie extension)
# ---------------------------------------------------------------------------

"""
    plot_variable(rep, var; dir=".") -> Vector{String}

Figures of a compared variable (`var` a key `file:var`, or a variable name):
the ref (converted), the run and their difference at the first differing
record and at the last, and the maximum |difference| per record. Defined by
the CairoMakie extension: `using CairoMakie` first.
"""
function plot_variable end

"The result of a variable for plots: by key, else by ref or run variable name."
function find_result(rep::Report, var::AbstractString)
    for r in rep.results
        r.key == var && return r
    end
    for r in rep.results
        (var in r.ref_vars || r.run_var == var) && !(r.status in ("absent", "new")) && return r
    end
    error("plot: no compared variable $var")
end

"""
    load_record(rep, r, k) -> (expected, run)

The expected values (ref, converted) and the run values of result `r` at its
aligned record `k` (1-based; static variables: 1).
"""
function load_record(rep::Report, r::VarResult, k::Integer)
    a = NCDataset(r.ref_path) do ds
        arrs = Dict{Symbol,Any}()
        dn = String[]
        for v in r.ref_vars
            x, dn = read_raw(ds, v)
            arrs[Symbol(v)] = x
        end
        if r.class == "C"
            env = Dict{Symbol,Any}(Symbol(k) => v for (k, v) in rep.constants)
            merge!(env, arrs)
            length(r.ref_vars) == 1 && (env[:x] = arrs[Symbol(r.ref_vars[1])])
            x = eval_expr(parse_expr(r.expr), env)
        else
            x = arrs[Symbol(r.ref_vars[1])]
        end
        t = findfirst(==("time"), dn)
        t === nothing ? x : collect(selectdim(x, t, r.ref_recs[k]))
    end
    b = NCDataset(r.run_path) do ds
        y, dn = read_raw(ds, r.run_var)
        t = findfirst(==("time"), dn)
        t === nothing ? y : collect(selectdim(y, t, r.run_recs[k]))
    end
    return a, b
end

end # module
