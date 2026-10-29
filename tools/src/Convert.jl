"""
    KryosTools.Convert

A yelmox par file converted into a kryos par file (regression-design.md §2,
core-design.md §5), fully expanded, with its accounting.

- `convert_yelmox(; yelmox, par, overrides, kryos, keep_library)` -> `Conversion`:
  resolves the yelmox configuration (`FesmUtils.Namelist.resolve`: the par file
  with runme's `-p` overrides applied, the groups it names from the library and
  the defaults), refuses what phase 1 does not do, builds the kryos groups,
  accounts for every yelmox key, validates every group against the kryos
  templates and checks the written file by reading it back.
- `write_par(path, conv)`, `write_report(path, conv)`: the par file and the
  accounting report (markdown).
- `account_summary(conv)`: counts of mapped, dropped and default-equal keys.
- `validate_conversion(conv)`, `check_roundtrip(conv, path)`: the checks of §2.6
  (both also run by `convert_yelmox` / `write_par`).

The rules are data: `tools/data/yelmox2kryos.toml` (refusals, group and key
mapping, output switches, names, templates); the rules it names (`rule = ...`)
are the functions `rule_*` below. The schema of the kryos side is read from the
kryos checkout: `input/kryos_defaults.nml` (core and backend templates),
`input/snapshots_defaults.nml`, `input/forcing_defaults.nml`,
`input/datasrc_defaults.nml` (else `fesm-utils/par/datasrc_defaults.nml`).

Every key of every resolved yelmox group is accounted as `:mapped`
(carried to the kryos file with the resolved value), `:default` (carried, and
equal to the default: the yelmox value came from its defaults file, or equals
the kryos template default) or `:dropped` (with a reason). Keys of a dropped
group that yelmox filled from its defaults are not counted.
"""
module Convert

using TOML
using SHA
using Printf
import JSON
import FesmUtils
import FesmUtils.Namelist as NL

export convert_yelmox, convert_file, Conversion, ConversionError, write_par, write_report,
       account_summary, validate_conversion, check_roundtrip, par_header, CONVERTER_VERSION, DEFAULT_SPEC
export merge_par, write_merged, Merged, merge_report, split_merged

const CONVERTER_VERSION = "0.1"
const DEFAULT_SPEC = normpath(joinpath(@__DIR__, "..", "data", "yelmox2kryos.toml"))

struct ConversionError <: Exception
    msg::String
end
Base.showerror(io::IO, e::ConversionError) = print(io, "ConversionError: ", e.msg)
fail(msg) = throw(ConversionError(msg))

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

"One key of a kryos group: `src` is the resolved yelmox value it carries (`nothing` if computed)."
mutable struct OutKey
    name::String
    raw::String
    value::Any
    origin::String            # "yelmox g.k", "fixed", "template", "rule: ..."
    comment::String
    src::Any
    set::Bool                 # set by a mapping or a fixed value (not only the template default)
    default::Any              # the template default (nothing: none)
end

"One kryos group, in file order. `templates`: the (file id, group) pairs that declare its keys."
mutable struct OutGroup
    name::String
    role::String              # instance role, "" for core and physics groups
    kind::String              # kryos | clock | schedule | globe | domain | instance | datasrc | physics
    component::String
    templates::Vector{Tuple{String,String}}
    keys::Vector{OutKey}
    yelmox::String            # the yelmox group a physics group is copied from
end

"The accounting of one yelmox key."
struct Acc
    group::String
    key::String
    status::Symbol            # :mapped | :default | :dropped
    dst::String
    reason::String
    source::String            # base | library | default | override
end

struct Conversion
    yelmox::String            # root
    par::String               # the source par file
    overrides::Vector{String}
    commit::String
    kryos::String             # kryos checkout the templates were read from
    keep_library::Bool
    domain::String
    sfx::String
    names::Dict{String,String}
    groups::Vector{OutGroup}
    acc::Vector{Acc}
    unfilled::Dict{String,Int}   # dropped group => keys yelmox filled from its defaults (not counted)
    warnings::Vector{String}
    notes::Vector{String}
    files::Vector{Pair{String,String}}   # relative path => sha256 (the files the par file names)
    checks::Vector{String}               # validation lines
    refs::Vector{Tuple{String,String}}   # (naming group, named group) pairs followed by the resolution
    library::Vector{String}
    spec::Dict{String,Any}
end

# ---------------------------------------------------------------------------
# Spec and templates
# ---------------------------------------------------------------------------

load_spec(path::AbstractString=DEFAULT_SPEC) = TOML.parsefile(path)

"The groups of the kryos template files, by id (`kryos_defaults`, ...)."
struct Templates
    files::Dict{String,Vector{NL.Group}}
    paths::Dict{String,String}
end

function template_paths(kryos::AbstractString)
    inp = joinpath(kryos, "input")
    ds = joinpath(inp, "datasrc_defaults.nml")
    isfile(ds) || (ds = joinpath(kryos, "fesm-utils", "par", "datasrc_defaults.nml"))
    return Dict("kryos_defaults" => joinpath(inp, "kryos_defaults.nml"),
                "snapshots_defaults" => joinpath(inp, "snapshots_defaults.nml"),
                "forcing_defaults" => joinpath(inp, "forcing_defaults.nml"),
                "datasrc_defaults" => ds)
end

function load_templates(kryos::AbstractString)
    paths = template_paths(kryos)
    isfile(paths["kryos_defaults"]) || fail("no $(paths["kryos_defaults"]): is $kryos a kryos checkout?")
    files = Dict{String,Vector{NL.Group}}()
    for (id, p) in paths
        isfile(p) && (files[id] = NL.read_namelist(p))
    end
    return Templates(files, paths)
end

function tgroup(T::Templates, id::AbstractString, group::AbstractString)
    haskey(T.files, id) || return nothing
    for g in T.files[id]
        g.name == group && return g
    end
    return nothing
end

# ---------------------------------------------------------------------------
# yelmox resolution
# ---------------------------------------------------------------------------

struct Resolved
    cfg::NL.Config
    groups::Dict{String,NL.RGroup}
    order::Vector{String}
    refs::Vector{NL.GroupRef}
end

function build_refs(spec)
    refs = NL.GroupRef[]
    for r in get(get(spec, "yelmox", Dict()), "ref", [])
        from = Regex(r["from"])
        when = haskey(r, "when") ? (String(r["when"][1]) => r["when"][2]) : nothing
        prefix = get(r, "prefix", "")
        tmpl = haskey(r, "template") ? (String(r["template"][1]), String(r["template"][2])) : nothing
        if haskey(r, "keys")
            for k in r["keys"]
                push!(refs, NL.GroupRef(from, k; prefix=prefix, when=when, template=tmpl))
            end
        elseif haskey(r, "map")
            for (k, g) in r["map"]
                push!(refs, NL.GroupRef(from, k; prefix=prefix, template=(String(r["file"]), String(g))))
            end
        else
            push!(refs, NL.GroupRef(from, r["key"]; prefix=prefix, when=when, template=tmpl))
        end
    end
    return refs
end

"The names a group lists under the rule `r` (`when` honoured), or nothing."
function ref_names(rg::NL.RGroup, r::NL.GroupRef)
    if r.when !== nothing
        v = NL.value_of(rg, r.when.first, r.when.second)
        v == r.when.second || return String[]
    end
    key = NL.ref_key(rg, r.key)
    key === nothing && return String[]
    return [r.prefix * n for n in NL.names_of(rg, key)]
end

"""
    resolve_yelmox(root, par, overrides, spec) -> Resolved

The configuration yelmox sees: `par` with the runme `-p` `overrides` applied
(`["g.k=v", "g.k==a,b"]`), the groups it names from the library (the files of
`[forcing] var_defs`) and from the defaults files of `.runme/info.json`.
"""
function resolve_yelmox(root::AbstractString, par::AbstractString, overrides::AbstractVector{<:AbstractString}, spec)
    info = joinpath(root, ".runme", "info.json")
    isfile(info) || fail("no $info: is $root a yelmox checkout?")
    cfg0 = NL.load_config(root)
    schema = Pair{Any,Tuple{String,String}}[p for p in cfg0.schema]
    for s in get(get(spec, "yelmox", Dict()), "schema", [])
        push!(schema, String(s["name"]) => (String(s["file"]), String(s["group"])))
    end
    # runme's group aliases (-p bsl.k=v)
    aliases = Dict{String,String}()
    try
        j = JSON.parsefile(info)
        for (a, g) in get(j, "grp_aliases", Dict())
            aliases[String(a)] = String(g)
        end
    catch
    end
    ovr = String[]
    for o in overrides
        k = split(o, '=', limit=2)[1]
        g = split(k, '.', limit=2)[1]
        push!(ovr, haskey(aliases, g) ? replace(o, g * "." => aliases[g] * "."; count=1) : String(o))
    end
    parsed = NL.parse_overrides(ovr)
    # the library: [forcing] var_defs of the par file, an override, else the defaults
    lib = var_defs(root, par, parsed, schema, spec)
    cfg = NL.Config(root; schema=schema, library=lib)
    refs = build_refs(spec)
    groups, order = try
        NL.resolve(cfg, par, parsed; follow=refs)
    catch e
        fail("cannot resolve $par: " * sprint(showerror, e))
    end
    return Resolved(cfg, groups, order, refs)
end

function var_defs(root, par, parsed, schema, spec)
    key = get(get(get(spec, "yelmox", Dict()), "library", Dict()), "key", "forcing.var_defs")
    g, k = split(key, '.', limit=2)
    for (ok, toks, _) in parsed
        ok == key && return String.(strip.(toks, '"'))
    end
    for grp in NL.read_namelist(par)
        grp.name == g || continue
        for kk in grp.keys
            kk.name == k && return string.(_as_vec(NL.parse_value(kk.raw)))
        end
    end
    for (pat, t) in schema
        (pat isa Regex ? occursin(pat, g) : pat == g) || continue
        f = joinpath(root, t[1])
        isfile(f) || continue
        for grp in NL.read_namelist(f)
            grp.name == t[2] || continue
            for kk in grp.keys
                kk.name == k && return string.(_as_vec(NL.parse_value(kk.raw)))
            end
        end
    end
    return String[]
end

_as_vec(v::AbstractVector) = v
_as_vec(v) = [v]

"The resolved value of `group.key`, else its yelmox default, else `default`."
function getval(R::Resolved, gk::AbstractString, default=nothing)
    g, k = split(gk, '.', limit=2)
    if haskey(R.groups, g)
        haskey(R.groups[g].entries, k) && return R.groups[g].entries[k].value
    end
    t = NL.template(R.cfg, g)
    if t !== nothing
        dg = NL.cfg_find(R.cfg, t[1], t[2])
        if dg !== nothing
            for kk in dg.keys
                kk.name == k && return NL.parse_value(kk.raw)
            end
        end
    end
    return default
end

function referenced_groups(R::Resolved)
    pairs = Tuple{String,String}[]
    for name in R.order, r in R.refs
        occursin(r.from, name) || continue
        rg = R.groups[name]
        for n in ref_names(rg, r)
            push!(pairs, (name, n))
        end
    end
    return pairs
end

# ---------------------------------------------------------------------------
# Refusals (§2.5)
# ---------------------------------------------------------------------------

function same(a, b)
    (a isa Bool) != (b isa Bool) && return false
    return a == b
end

function check_refusals(R::Resolved, spec)
    bad = String[]
    for r in get(spec, "refuse_group", [])
        re = Regex(r["pattern"])
        for name in R.order
            occursin(re, name) &&
                push!(bad, "group &$name: $(r["why"]) (phase $(r["phase"]))")
        end
    end
    for r in get(spec, "refuse", [])
        if haskey(r, "when")                      # only when that key has that value
            wv = getval(R, r["when"][1])
            (wv !== nothing && same(wv, r["when"][2])) || continue
        end
        ks = haskey(r, "keys") ? r["keys"] : [r["key"]]
        for k in ks
            v = getval(R, k)
            v === nothing && continue
            vals = v isa AbstractVector ? v : [v]
            for x in vals
                if haskey(r, "allow") && !any(a -> same(x, a), r["allow"])
                    push!(bad, "$k = $(NL.format_value(x)): $(r["why"]) (phase $(r["phase"]))")
                end
                if haskey(r, "deny")
                    for d in r["deny"]
                        hit = startswith(d, "re:") ? (x isa AbstractString && occursin(Regex(d[4:end]), x)) : same(x, d)
                        hit && push!(bad, "$k = $(NL.format_value(x)): $(r["why"]) (phase $(r["phase"]))")
                    end
                end
            end
        end
    end
    isempty(bad) || fail("yelmox configuration not convertible in phase 1:\n  " * join(unique(bad), "\n  "))
    return nothing
end

# ---------------------------------------------------------------------------
# Building the kryos groups
# ---------------------------------------------------------------------------

mutable struct State
    R::Resolved
    spec::Dict{String,Any}
    T::Templates
    ctx::Dict{String,String}
    groups::Vector{OutGroup}
    acc::Vector{Acc}
    unfilled::Dict{String,Int}
    warnings::Vector{String}
    notes::Vector{String}
    keep_library::Bool
    referenced::Set{String}
    done::Set{String}                 # rule groups already built (ghf, sed)
end

subst(S::State, s::AbstractString) = begin
    out = String(s)
    for _ in 1:3
        for (k, v) in S.ctx
            out = replace(out, "{" * k * "}" => v)
        end
    end
    out
end
subst(S::State, v) = v

function find_group(S::State, name::AbstractString)
    for g in S.groups
        g.name == name && return g
    end
    return nothing
end
function find_role(S::State, role::AbstractString)
    for g in S.groups
        g.role == role && return g
    end
    return nothing
end
function need_group(S::State, name::AbstractString)
    g = find_group(S, name)
    g === nothing && fail("kryos group &$name is not built (a mapping names it)")
    return g
end

function template_keys(S::State, id::AbstractString, group::AbstractString)
    g = tgroup(S.T, id, group)
    g === nothing && return nothing
    return g.keys
end

"Add the keys of template `(id, group)` (with their defaults) to `og` unless it has them."
function add_template!(S::State, og::OutGroup, id::AbstractString, group::AbstractString)
    ks = template_keys(S, id, group)
    ks === nothing && fail("template group &$group not in $(S.T.paths[id])")
    push!(og.templates, (String(id), String(group)))
    for k in ks
        any(x -> x.name == k.name, og.keys) && continue
        push!(og.keys, OutKey(k.name, k.raw, NL.parse_value(k.raw), "template", "", nothing, false, NL.parse_value(k.raw)))
    end
    return og
end

function new_outgroup(S::State, entry)
    name = subst(S, entry["group"])
    role = get(entry, "role", "")
    kind = if role != ""
        "instance"
    elseif startswith(name, "out_")
        "schedule"
    elseif name == S.ctx["clock"]
        "clock"
    elseif name == S.ctx["globe"]
        "globe"
    elseif name == "kryos"
        "kryos"
    elseif name == S.ctx["domain"]
        "domain"
    else
        "datasrc"
    end
    og = OutGroup(name, role, kind, get(entry, "component", ""), Tuple{String,String}[], OutKey[], "")
    for t in entry["templates"]
        id, g = String.(S.spec["templates"][t])
        add_template!(S, og, id, g)
    end
    return og
end

function setkey!(S::State, og::OutGroup, key::AbstractString, raw::AbstractString; origin::AbstractString="",
                 comment::AbstractString="", src=nothing, soft::Bool=false)
    i = findfirst(k -> k.name == key, og.keys)
    if i === nothing
        soft && return false
        decl = join(("$(t[1]):$(t[2])" for t in og.templates), ", ")
        fail("key '$key' of kryos group &$(og.name) is declared by no template ($decl)")
    end
    k = og.keys[i]
    k.raw = String(raw)
    k.value = NL.parse_value(raw)
    k.origin = String(origin)
    k.comment = String(comment)
    k.src = src
    k.set = true
    return true
end

function setval!(S::State, og::OutGroup, key::AbstractString, v; kw...)
    return setkey!(S, og, key, NL.format_value(v); kw...)
end

function split_dst(S::State, d::AbstractString)
    s = subst(S, d)
    i = findfirst('.', s)
    i === nothing && fail("destination '$d' is not group.key")
    return s[1:i-1], s[i+1:end]
end

function entry_of(S::State, g::AbstractString, k::AbstractString)
    return S.R.groups[g].entries[k]
end

"Status of a carried key: `:default` when yelmox took it from its defaults or it equals the kryos default."
function classify(src_source::AbstractString, value, default)
    src_source == "default" && return :default
    default !== nothing && same_value(value, default) && return :default
    return :mapped
end

same_value(a, b) = typeof(a) == typeof(b) ? a == b : (a isa Number && b isa Number && !(a isa Bool) && !(b isa Bool) && a == b)

function account!(S::State, g, k, status, dst, reason, source)
    push!(S.acc, Acc(String(g), String(k), status, String(dst), String(reason), String(source)))
end

"Carry the yelmox key `g.k` to every destination in `dsts` (raw text copied)."
function carry!(S::State, g, k, dsts)
    e = entry_of(S, g, k)
    dsts = dsts isa AbstractVector ? dsts : [dsts]
    status = :default
    names = String[]
    for d in dsts
        gn, kn = split_dst(S, d)
        og = need_group(S, gn)
        i = findfirst(x -> x.name == kn, og.keys)
        dflt = i === nothing ? nothing : og.keys[i].default
        setkey!(S, og, kn, e.raw; origin="yelmox $g.$k", comment="yelmox $g.$k", src=e.value)
        s = classify(e.source, e.value, dflt)
        s == :mapped && (status = :mapped)
        push!(names, "$gn.$kn")
    end
    account!(S, g, k, status, join(names, ", "), "", e.source)
end

# ---------- rules ----------

function rule_calendar_const(S::State, g, k, row)
    e = entry_of(S, g, k)
    method = getval(S.R, "ctrl.tstep_method", "const")
    cal = getval(S.R, "ctrl.calendar", false)
    og = need_group(S, S.ctx["clock"])
    if method == "const"
        key = cal === true ? "time_const_cal" : "time_const"
        # the clock template holds both as optional: only one is a key of a const clock
        push!(og.keys, OutKey(key, "0.0", 0.0, "template", "", nothing, false, 0.0))
        setkey!(S, og, key, e.raw; origin="yelmox $g.$k", comment="yelmox $g.$k", src=e.value)
        account!(S, g, k, classify(e.source, e.value, 0.0), "$(og.name).$key", "", e.source)
    else
        account!(S, g, k, :dropped, "", "unused: tstep_method = $method", e.source)
    end
end

function rule_equil_time(S::State, g, k, row)
    method = getval(S.R, "sim.init_method", "none")
    og = need_group(S, S.ctx["ice"])
    e = entry_of(S, g, k)
    if k == "init_method"
        if method == "equil"
            haskey(S.R.groups["sim"].entries, "init_equil_time") || fail("sim.init_method = equil without sim.init_equil_time")
            t = S.R.groups["sim"].entries["init_equil_time"]
            setkey!(S, og, "init_equil_time", t.raw; origin="yelmox sim.init_equil_time", comment="yelmox sim.init_method = equil, sim.init_equil_time", src=t.value)
            account!(S, g, k, :mapped, "$(og.name).init_equil_time", "init_method = equil: the key carries init_equil_time", e.source)
        else
            setkey!(S, og, "init_equil_time", "0.0"; origin="yelmox sim.init_method = none", comment="yelmox sim.init_method = none", src=nothing)
            account!(S, g, k, :default, "$(og.name).init_equil_time", "init_method = $method: init_equil_time = 0 (none)", e.source)
        end
    else  # init_equil_time
        if method == "equil"
            account!(S, g, k, e.source == "default" ? :default : :mapped, "$(og.name).init_equil_time", "", e.source)
        else
            account!(S, g, k, :dropped, "", "unused: init_method = $method", e.source)
        end
    end
end

function rule_switch(S::State, g, k, row)
    e = entry_of(S, g, k)
    on = e.value === true
    streams = [s for s in get(S.spec, "stream", []) if s["switch"] == "$g.$k"]
    isempty(streams) && fail("no [[stream]] for the output switch $g.$k")
    set = String[]
    for s in streams, (role, key, sched) in s["set"]
        on || continue
        og = find_role(S, role)
        if og === nothing
            role == "domain" || fail("stream of role $role: no such instance")
            og = need_group(S, S.ctx["domain"])
        end
        ok = setkey!(S, og, key, NL.format_value(sched); origin="yelmox $g.$k", comment="yelmox $g.$k", soft=true)
        if ok
            push!(set, "$(og.name).$key")
        else
            push!(S.warnings, "output $g.$k: key `$key` of &$(og.name) is not declared by the templates ($(join(("$(t[1]):$(t[2])" for t in og.templates), ", "))): the stream (schedule `$sched`) is not written")
        end
    end
    if on
        account!(S, g, k, e.source == "default" ? :default : :mapped, join(set, ", "), "", e.source)
    else
        account!(S, g, k, :default, "", "switch off: the streams stay none", e.source)
    end
end

function rule_library(S::State, g, k, row)
    e = entry_of(S, g, k)
    if S.keep_library
        og = need_group(S, "kryos")
        setkey!(S, og, "library", e.raw; origin="yelmox $g.$k", comment="yelmox $g.$k", src=e.value)
        account!(S, g, k, :mapped, "kryos.library", "", e.source)
    else
        account!(S, g, k, :dropped, "", "the library groups are copied into the file (fully expanded; --keep-library keeps the references)", e.source)
    end
end

function datasrc_fields!(S::State, name::AbstractString, kv::Vector{<:Pair})
    og = need_group(S, name)
    for (k, v) in kv
        setval!(S, og, k, v; origin="rule", comment="")
    end
end

function rule_ghf_field(S::State, g, k, row)
    e = entry_of(S, g, k)
    use_obs = getval(S.R, "ghf.use_obs", false) === true
    field = S.ctx["Q_geo"]
    inst = find_role(S, "geothermal")
    rules = get(get(S.spec, "rules", Dict()), "ghf_field", Dict())
    if !("ghf" in S.done)
        push!(S.done, "ghf")
        if use_obs
            datasrc_fields!(S, field, ["filename" => getval(S.R, "ghf.obs_path"), "name" => getval(S.R, "ghf.obs_name"),
                                       "mode" => "static", "remap" => rules["remap"], "fill" => rules["fill"]])
            setval!(S, inst, "Q_geo_min", rules["ghf_min"]; origin="rule: yelmox geothermal.f90 floor ghf_min = 0.1 mW m-2")
        else
            c = getval(S.R, "ghf.ghf_const")
            datasrc_fields!(S, field, ["mode" => "const", "value" => c / 1000.0])
        end
    end
    mk(d, why="") = account!(S, g, k, e.source == "default" ? :default : :mapped, d, why, e.source)
    if k == "use_obs"
        mk("$field.mode", use_obs ? "static field" : "const")
    elseif k == "obs_path"
        use_obs ? mk("$field.filename") : account!(S, g, k, :dropped, "", "unused: use_obs = False", e.source)
    elseif k == "obs_name"
        use_obs ? mk("$field.name") : account!(S, g, k, :dropped, "", "unused: use_obs = False", e.source)
    elseif k == "ghf_const"
        use_obs ? account!(S, g, k, :dropped, "", "unused: use_obs = True", e.source) :
                  mk("$field.value", "mW m-2 -> W m-2 (/1000)")
    elseif k == "obs_err_name"
        account!(S, g, k, :dropped, "", "ensemble perturbation of the heat flux (f_stdev = 0)", e.source)
    elseif k == "f_stdev"
        account!(S, g, k, :dropped, "", "ensemble perturbation of the heat flux; = 0 (any other value refused)", e.source)
    else
        fail("ghf.$k: no rule")
    end
end

function rule_sed_field(S::State, g, k, row)
    e = entry_of(S, g, k)
    method = getval(S.R, "sed.method", 0)
    field = S.ctx["H_sed"]
    rules = get(get(S.spec, "rules", Dict()), "sed_field", Dict())
    if !("sed" in S.done)
        push!(S.done, "sed")
        if method == 1
            datasrc_fields!(S, field, ["filename" => getval(S.R, "sed.obs_path"), "name" => getval(S.R, "sed.obs_name"),
                                       "mode" => "static", "remap" => rules["remap"], "fill" => rules["fill"]])
        else
            datasrc_fields!(S, field, ["mode" => "const", "value" => 0.0])
        end
    end
    mk(d, why="") = account!(S, g, k, e.source == "default" ? :default : :mapped, d, why, e.source)
    if k == "method"
        mk("$field.mode", method == 1 ? "static field" : "const 0 m (no sediments)")
    elseif k == "obs_path"
        method == 1 ? mk("$field.filename") : account!(S, g, k, :dropped, "", "unused: method = 0", e.source)
    elseif k == "obs_name"
        method == 1 ? mk("$field.name") : account!(S, g, k, :dropped, "", "unused: method = 0", e.source)
    else
        fail("sed.$k: no rule")
    end
end

const RULES = Dict{String,Function}("calendar_const" => rule_calendar_const, "equil_time" => rule_equil_time,
                                    "switch" => rule_switch, "library" => rule_library,
                                    "ghf_field" => rule_ghf_field, "sed_field" => rule_sed_field)

# ---------- groups ----------

function match_src(src::AbstractString, name::AbstractString)
    startswith(src, "re:") ? occursin(Regex(src[4:end]), name) : src == name
end

function find_key_row(S::State, g, k)
    gk = "$g.$k"
    for row in get(S.spec, "key", [])
        match_src(row["src"], gk) && return row
    end
    return nothing
end

function apply_keys_group!(S::State, name::AbstractString)
    rg = S.R.groups[name]
    for k in rg.order
        row = find_key_row(S, name, k)
        row === nothing && fail("yelmox key $name.$k is not accounted for (no [[key]] row in the spec)")
        if haskey(row, "defer")
            continue                          # carried or dropped by apply_experiments!
        elseif haskey(row, "drop")
            account!(S, name, k, :dropped, "", row["drop"], rg.entries[k].source)
        elseif haskey(row, "rule")
            haskey(RULES, row["rule"]) || fail("[[key]] $name.$k: unknown rule $(row["rule"])")
            RULES[row["rule"]](S, name, k, row)
        elseif haskey(row, "dst")
            carry!(S, name, k, row["dst"])
        else
            fail("[[key]] row for $name.$k has no dst, rule or drop")
        end
    end
end

function dtt_value(S::State)
    v = getval(S.R, "ctrl.dtt")
    v isa Number || fail("ctrl.dtt missing")
    return Float64(v)
end

function apply_into_group!(S::State, name::AbstractString, dst::AbstractString)
    rg = S.R.groups[name]
    og = need_group(S, dst)
    dtt = dtt_value(S)
    for k in rg.order
        e = rg.entries[k]
        i = findfirst(x -> x.name == k, og.keys)
        i === nothing && fail("yelmox key $name.$k has no key in the kryos schedule template (&$dst)")
        dflt = og.keys[i].default
        if k == "dt" && e.value isa Number && getval(S.R, "$name.method", "const") == "const" && e.value > 0 && e.value < dtt
            setval!(S, og, "dt", dtt; origin="yelmox $name.dt", comment="yelmox $name.dt = $(NL.format_value(e.value)) raised to dtt", src=nothing)
            push!(S.notes, "$name.dt = $(NL.format_value(e.value)) is below dtt = $(NL.format_value(dtt)): [$dst] dt = $(NL.format_value(dtt)) gives the same output times (fesm_clock rejects dt < dtt)")
            account!(S, name, k, :mapped, "$dst.dt", "raised to dtt: same output times", e.source)
            continue
        end
        setkey!(S, og, k, e.raw; origin="yelmox $name.$k", comment="yelmox $name.$k", src=e.value)
        account!(S, name, k, classify(e.source, e.value, dflt), "$dst.$k", "", e.source)
    end
end

function physics_group!(S::State, name::AbstractString, rule)
    rg = S.R.groups[name]
    og = OutGroup(name, "", "physics", "", Tuple{String,String}[], OutKey[], name)
    drops = get(S.spec, "drop_keys", Dict())
    rst = Set(get(get(S.spec, "restart_none", Dict()), "groups", String[]))
    moves = [m for m in get(S.spec, "move", []) if m["from"] == name]   # keys that go to another physics group
    targets = Dict{String,OutGroup}()
    for k in rg.order
        e = rg.entries[k]
        mi = findfirst(m -> occursin(Regex(m["keys"]), k), moves)
        if mi !== nothing && !haskey(drops, "$name.$k")
            tn = subst(S, moves[mi]["to"])
            tog = get!(targets, tn) do
                OutGroup(tn, "", "physics", "", Tuple{String,String}[], OutKey[], name)
            end
            cm = String(first(split(e.comment, '\n')))
            push!(tog.keys, OutKey(k, e.raw, parse_exact(e.raw), "yelmox $name.$k", cm, e.value, true, nothing))
            account!(S, name, k, e.source == "default" ? :default : :mapped, "$tn.$k", moves[mi]["why"], e.source)
            continue
        end
        if haskey(drops, "$name.$k")
            account!(S, name, k, :dropped, "", drops["$name.$k"], e.source)
            continue
        end
        raw, val, srcv = e.raw, e.value, e.value
        why = ""
        if k == "restart" && name in rst && val != "None"
            raw, val, srcv = "\"None\"", "None", nothing
            why = "restart forced to \"None\" (restarts come from the bundle); was $(NL.format_value(e.value))"
            push!(S.notes, "[$name] restart $(NL.format_value(e.value)) -> \"None\"")
        end
        cm = String(first(split(e.comment, '\n')))
        push!(og.keys, OutKey(k, raw, parse_exact(raw), "yelmox $name.$k", cm, srcv, true, nothing))
        status = e.source == "default" ? :default : :mapped
        account!(S, name, k, status, "$name.$k", why, e.source)
    end
    push!(S.groups, og)
    for (_, tog) in sort!(collect(targets); by=first)
        push!(S.groups, tog)
        push!(S.notes, "[$name] keys moved to [$(tog.name)] ($(join((k.name for k in tog.keys), ", "))): the shelf component owns them")
    end
end

parse_exact(raw) = NL.parse_value(raw)

function drop_group!(S::State, name::AbstractString, why::AbstractString)
    rg = S.R.groups[name]
    nfill = 0
    for k in rg.order
        e = rg.entries[k]
        if e.source == "default"
            nfill += 1
        else
            account!(S, name, k, :dropped, "", why, e.source)
        end
    end
    S.unfilled[name] = nfill
end

function apply_groups!(S::State)
    rules = get(S.spec, "group", [])
    for name in S.R.order
        rule = nothing
        for r in rules
            if haskey(r, "when")                  # the rule applies only when that key has that value
                wv = getval(S.R, r["when"][1])
                (wv !== nothing && same(wv, r["when"][2])) || continue
            end
            match_src(r["src"], name) && (rule = r; break)
        end
        rule === nothing && fail("yelmox group &$name has no [[group]] rule in the spec")
        act = rule["action"]
        if act == "keys"
            apply_keys_group!(S, name)
        elseif act == "into"
            apply_into_group!(S, name, subst(S, rule["dst"]))
        elseif act == "drop"
            drop_group!(S, name, rule["why"])
        elseif act == "experiment"
            # built by apply_experiments! (the keys of a group may go to several experiments)
        elseif act == "copy"
            rg = S.R.groups[name]
            if get(rule, "referenced", false) && !(name in S.referenced)
                drop_group!(S, name, "not read by yelmox: no key of the resolved configuration names it")
            elseif get(rule, "library", false) && S.keep_library && rg.source == "library"
                for k in rg.order
                    account!(S, name, k, :mapped, "(library) $name.$k", "stays in the library group", rg.entries[k].source)
                end
            else
                physics_group!(S, name, rule)
            end
        else
            fail("[[group]] $(rule["src"]): unknown action $act")
        end
    end
end

# ---------- experiments ----------

"Is the yelmox group `name` one of the `keys` groups (its keys go through [[key]] rows)?"
function is_keys_group(S::State, name::AbstractString)
    for r in get(S.spec, "group", [])
        match_src(r["src"], name) && return get(r, "action", "") == "keys"
    end
    return false
end

"Is the [[experiment]] row active? Returns (active, why-not)."
function experiment_active(S::State, row)
    if haskey(row, "when")
        v = getval(S.R, row["when"][1], false)
        same(v, row["when"][2]) || return false, "$(row["when"][1]) = $(NL.format_value(v))"
    end
    if haskey(row, "skip_if")
        v = getval(S.R, row["skip_if"][1])
        v !== nothing && same(v, row["skip_if"][2]) && return false, "$(row["skip_if"][1]) = $(NL.format_value(v))"
    end
    return true, ""
end

"""
The experiment groups (core-design 3.9) of the [[experiment]] rows of the spec: one kryos group per
active row (templates `[experiment]` + `[experiment_<kind>]`, `kind` set, the keys of the yelmox
group carried by name or by the row's `map`), `{domain}.experiments` listing them in row order.
A key of a yelmox group named by a row that no active row carries is dropped, with the reason of
the inactive row that would have taken it.
"""
function apply_experiments!(S::State)
    rows = get(S.spec, "experiment", [])
    names = String[]
    carried = Dict{String,Set{String}}()
    why = Dict{String,String}()
    for row in rows
        from = row["from"]
        haskey(S.R.groups, from) || begin
            active, _ = experiment_active(S, row)
            active && fail("experiment `$(row["kind"])` is active ($(row["when"][1]) = $(NL.format_value(row["when"][2]))) but the par file has no [$from] group")
            continue
        end
        active, reason = experiment_active(S, row)
        if !active
            why[from] = get(why, from, "") * (isempty(get(why, from, "")) ? "" : "; ") * "no `$(row["kind"])` experiment ($reason)"
            continue
        end
        name = subst(S, row["name"])
        og = OutGroup(name, "", "experiment", "", Tuple{String,String}[], OutKey[], from)
        add_template!(S, og, "kryos_defaults", "experiment")
        add_template!(S, og, "kryos_defaults", "experiment_" * row["kind"])
        push!(S.groups, og)
        setval!(S, og, "kind", String(row["kind"]); origin="fixed")
        rg = S.R.groups[from]
        map = get(row, "map", nothing)
        got = get!(carried, from, Set{String}())
        for k in rg.order
            kk = map === nothing ? k : get(map, k, nothing)
            kk === nothing && continue
            carry!(S, from, k, "$name.$kk")
            push!(got, k)
        end
        push!(names, name)
    end
    for from in unique(String(r["from"]) for r in rows)
        haskey(S.R.groups, from) || continue
        rg = S.R.groups[from]
        got = get(carried, from, Set{String}())
        nfill = 0
        for k in rg.order
            k in got && continue
            # a `keys` group (sim) has its other keys accounted by their [[key]] rows; only the deferred ones are ours
            if is_keys_group(S, from)
                row = find_key_row(S, from, k)
                (row !== nothing && haskey(row, "defer")) || continue
            end
            if rg.entries[k].source == "default"       # yelmox's own default, as in drop_group!
                nfill += 1
            else
                account!(S, from, k, :dropped, "", get(why, from, "no experiment takes it"), rg.entries[k].source)
            end
        end
        S.unfilled[from] = nfill
    end
    isempty(names) && return
    dom = need_group(S, S.ctx["domain"])
    setval!(S, dom, "experiments", names; origin="rule: experiments")
    push!(S.notes, "experiments of the domain, in call order: " * join(names, ", "))
end

# order of the physics groups in the file
const PHYSICS_ORDER = [r"^topo$", r"^regions$", r"^isos$", r"^barysealevel$", r"^snapshots$", r"^orbit$",
                       r"^surface_chion$", r"^chion$", r"^chion_const$", r"^\w+_chion$", r"^marine_shelf$", r"^shelf_corr_\w+$",
                       r"^yelmo$", r"^y\w+$", r"^yelmo_\w+$", r"^coef_\w+$", r"^term_\w+$", r"^snapshot_\w+$"]
function physics_rank(name)
    for (i, r) in enumerate(PHYSICS_ORDER)
        occursin(r, name) && return i
    end
    return length(PHYSICS_ORDER) + 1
end

# ---------------------------------------------------------------------------
# The conversion
# ---------------------------------------------------------------------------

function yelmox_commit(root::AbstractString)
    try
        sha = strip(read(pipeline(`git -C $root rev-parse HEAD`; stderr=devnull), String))
        br = strip(read(pipeline(`git -C $root rev-parse --abbrev-ref HEAD`; stderr=devnull), String))
        dirty = !isempty(strip(read(pipeline(`git -C $root status --porcelain --untracked-files=no`; stderr=devnull), String)))
        return string(sha, " (", br, dirty ? ", modified" : "", ")")
    catch
        return "unknown"
    end
end

"""
    convert_yelmox(; yelmox, par, overrides=[], kryos=<this checkout>, spec=DEFAULT_SPEC,
                   keep_library=false) -> Conversion

Convert the yelmox par file `par` (a path, absolute or relative to the
yelmox root `yelmox`) with the runme `-p` `overrides` applied. Stops with a
`ConversionError` naming the key and the phase when the configuration is
outside phase 1 (§2.5), when a yelmox key is not accounted for, or when a
group does not validate against the kryos templates.
"""
function convert_yelmox(; yelmox::AbstractString, par::AbstractString, overrides::AbstractVector{<:AbstractString}=String[],
                        kryos::AbstractString=normpath(joinpath(@__DIR__, "..", "..")),
                        spec::AbstractString=DEFAULT_SPEC, keep_library::Bool=false)
    root = abspath(yelmox)
    parpath = isabspath(par) ? par : (isfile(par) ? abspath(par) : joinpath(root, par))
    isfile(parpath) || fail("no par file $par")
    sp = load_spec(spec)
    T = load_templates(kryos)
    R = resolve_yelmox(root, parpath, overrides, sp)
    check_refusals(R, sp)

    domain = getval(R, "domain.name")
    domain isa AbstractString && !isempty(domain) || fail("[domain] name is not set")
    sfx = get(get(sp["names"], "suffix", Dict()), domain, lowercase(first(domain, 3)))
    ctx = Dict{String,String}("domain" => String(domain), "sfx" => String(sfx))
    for (k, v) in sp["names"]
        v isa AbstractString && (ctx[k] = replace(v, "{sfx}" => sfx))
    end

    S = State(R, sp, T, ctx, OutGroup[], Acc[], Dict{String,Int}(), String[], String[], keep_library,
              Set(n for (_, n) in referenced_groups(R)), Set{String}())

    # the kryos groups with their template keys
    for entry in sp["out"]
        push!(S.groups, new_outgroup(S, entry))
    end
    # instance backends (a template per backend): fixed ones, else from yelmox
    be = Dict{String,String}(k => String(v) for (k, v) in get(sp, "backend", Dict()))
    be["climate"] = string(getval(R, "comps.climate", "snapshots"))
    be["surface"] = string(getval(R, "comps.surface_method", "chion"))
    for og in S.groups
        og.kind == "instance" || continue
        b = get(be, og.role, "")
        isempty(b) && fail("no backend for the instance role $(og.role)")
        i = findfirst(k -> k.name == "backend", og.keys)
        og.keys[i].raw = NL.format_value(b); og.keys[i].value = b; og.keys[i].origin = "fixed"; og.keys[i].set = true
        add_template!(S, og, "kryos_defaults", "$(og.component)_$b")
    end
    # scope: keys that have no meaning in an instance of this role
    for (entry, og) in zip(sp["out"], S.groups)
        for k in get(entry, "skip", String[])
            filter!(x -> x.name != k, og.keys)
        end
    end
    # fixed values, and those of the backend each role has
    fixed = Pair{String,Any}[collect(sp["fixed"])...]
    for (role, byb) in get(sp, "fixed_backend", Dict())
        append!(fixed, collect(get(byb, get(be, role, ""), Dict{String,Any}())))
    end
    for (dst, v) in fixed
        gn, kn = split_dst(S, dst)
        og = need_group(S, gn)
        v2 = v isa AbstractString ? subst(S, v) : v
        setval!(S, og, kn, v2; origin="fixed")
    end

    apply_groups!(S)
    apply_experiments!(S)

    # physics groups after the core ones, in a fixed order
    phys = [g for g in S.groups if g.kind == "physics"]
    core = [g for g in S.groups if g.kind != "physics"]
    sort!(phys, by=g -> physics_rank(g.name), alg=MergeSort)
    S.groups = vcat(core, phys)
    names = [g.name for g in S.groups]
    length(unique(names)) == length(names) || fail("duplicate kryos group names: " * join(unique(n for n in names if count(==(n), names) > 1), ", "))

    # accounting: every key of every resolved group
    for name in R.order
        rg = R.groups[name]
        for k in rg.order
            any(a -> a.group == name && a.key == k, S.acc) ||
                (rg.entries[k].source == "default" && haskey(S.unfilled, name)) ||
                fail("yelmox key $name.$k is not accounted for")
        end
    end

    files = provenance_files(root, S, sp)
    conv = Conversion(root, parpath, collect(String, overrides), yelmox_commit(root), abspath(kryos), keep_library,
                      String(domain), String(sfx), ctx, S.groups, S.acc, S.unfilled, S.warnings, S.notes, files,
                      String[], referenced_groups(R), R.cfg.library, sp)
    validate_conversion(conv)
    return conv
end

function sha_hex(path)
    return bytes2hex(open(sha256, path))
end

function provenance_files(root, S::State, sp)
    rels = String[]
    for r in get(get(sp, "provenance", Dict()), "always", String[])
        push!(rels, r)
    end
    fldr = nothing
    for og in S.groups
        for k in og.keys
            vals = k.value isa AbstractVector ? k.value : [k.value]
            for v in vals
                v isa AbstractString || continue
                (occursin('{', v) || isempty(v) || isabspath(v)) && continue
                isfile(joinpath(root, v)) && push!(rels, v)
            end
        end
        if og.name == "orbit"
            for k in og.keys
                k.name == "fldr" && (fldr = k.value)
            end
        end
    end
    if fldr !== nothing
        for f in get(get(sp, "provenance", Dict()), "orbit", String[])
            push!(rels, joinpath(fldr, f))
        end
    end
    out = Pair{String,String}[]
    for r in unique(rels)
        p = joinpath(root, r)
        push!(out, r => (isfile(p) ? sha_hex(p) : "absent"))
    end
    return out
end

# ---------------------------------------------------------------------------
# Counts
# ---------------------------------------------------------------------------

"""
    summary(conv) -> NamedTuple

`mapped`, `default` and `dropped` key counts over every yelmox key, `total`,
and `unfilled` (keys yelmox filled from its defaults in groups that are
dropped; not counted).
"""
function account_summary(c::Conversion)
    n(s) = count(a -> a.status == s, c.acc)
    return (mapped=n(:mapped), default=n(:default), dropped=n(:dropped), total=length(c.acc),
            unfilled=sum(values(c.unfilled); init=0))
end

# ---------------------------------------------------------------------------
# Validation (§2.6 check 3, against the kryos templates)
# ---------------------------------------------------------------------------

function validation_templates(conv::Conversion, og::OutGroup, T::Templates, yroot, sp)
    og.kind != "physics" && return og.templates, ""
    for v in get(sp, "validate", [])
        hit = (haskey(v, "groups") && og.name in v["groups"]) || (haskey(v, "pattern") && occursin(Regex(v["pattern"]), og.name))
        hit || continue
        for (who, file, grp) in v["files"]
            g = grp == "=" ? og.name : grp
            if who == "kryos"
                for (id, p) in T.paths
                    endswith(p, file) || continue
                    tgroup(T, id, g) === nothing || return [(id, g)], ""
                end
            else
                p = joinpath(yroot, file)
                isfile(p) || continue
                any(x -> x.name == g, NL.read_namelist(p)) && return [("yelmox:" * file, g)], p
            end
        end
    end
    return Tuple{String,String}[], ""
end

"""
    validate_conversion(conv) -> Vector{String}

Every key of every group of the kryos file is declared by the template(s) of
its group (the Julia counterpart of `nml_validate`), the required keys are set
(`[clock]` time_init, time_end, dtt; `[domain]` grid_hub; `[instance]` backend;
schedules method), a const clock has exactly one of time_const, time_const_cal,
every schedule `dt` is a multiple of dtt and a whole number of ticks, and every
group a key names exists (in the file, or the library when kept). Throws a
`ConversionError` listing every failure; the lines it checked are in
`conv.checks`.
"""
function validate_conversion(conv::Conversion)
    sp = conv.spec
    T = load_templates(conv.kryos)
    errs = String[]
    lines = conv.checks
    empty!(lines)
    byname = Dict(g.name => g for g in conv.groups)
    cache = Dict{String,Vector{NL.Group}}()
    function declared(og)
        tpl, p = validation_templates(conv, og, T, conv.yelmox, sp)
        keys = Set{String}()
        for (id, g) in tpl
            if startswith(id, "yelmox:")
                f = id[8:end]
                gs = get!(() -> NL.read_namelist(joinpath(conv.yelmox, f)), cache, f)
                for x in gs
                    x.name == g && union!(keys, lowercase(k.name) for k in x.keys)
                end
            else
                x = tgroup(T, id, g)
                x === nothing || union!(keys, lowercase(k.name) for k in x.keys)
            end
        end
        return keys, tpl
    end
    nval = 0
    for og in conv.groups
        keys, tpl = declared(og)
        if isempty(tpl)
            push!(lines, "&$(og.name): no schema file (not validated)")
            continue
        end
        bad = [k.name for k in og.keys if !(lowercase(k.name) in keys)]
        isempty(bad) || push!(errs, "&$(og.name): keys declared by no template ($(join(("$(t[1]):$(t[2])" for t in tpl), ", "))): " * join(bad, ", "))
        nval += 1
        push!(lines, "&$(og.name): $(length(og.keys)) keys against " * join(("$(replace(t[1], "yelmox:" => "yelmox "))&$(t[2])" for t in tpl), " + ") *
                     (og.name == "orbit" && startswith(tpl[1][1], "yelmox:") ? " (kryos_defaults.nml has no &orbit template yet: yelmox's schema)" : ""))
    end
    # required keys and semantic checks
    function val(g, k)
        og = get(byname, g, nothing)
        og === nothing && return nothing
        i = findfirst(x -> x.name == k, og.keys)
        return i === nothing ? nothing : og.keys[i].value
    end
    function isset(g, k)
        og = get(byname, g, nothing)
        og === nothing && return false
        i = findfirst(x -> x.name == k, og.keys)
        return i !== nothing && og.keys[i].set
    end
    clockname = conv.names["clock"]
    for k in ("time_init", "time_end", "dtt")
        isset(clockname, k) || push!(errs, "&$clockname: required key $k is not set")
    end
    if val(clockname, "method") == "const"
        n = count(k -> haskey(Dict(x.name => x for x in byname[clockname].keys), k), ("time_const", "time_const_cal"))
        n == 1 || push!(errs, "&$clockname: method = const needs exactly one of time_const, time_const_cal (has $n)")
    end
    dtt = val(clockname, "dtt")
    tpy = val(clockname, "ticks_per_year")
    for og in conv.groups
        og.kind == "schedule" || continue
        m = val(og.name, "method")
        isset(og.name, "method") || push!(errs, "&$(og.name): required key method is not set")
        if m == "const"
            dt = val(og.name, "dt")
            if dt isa Number && dtt isa Number && tpy isa Number
                r = dt / dtt
                abs(r - round(r)) < 1e-9 && dt > 0 ||
                    push!(errs, "&$(og.name): dt = $dt is not a positive multiple of dtt = $dtt")
                t = dt * tpy
                abs(t - round(t)) < 1e-6 || push!(errs, "&$(og.name): dt = $dt is not a whole number of ticks ($tpy per year)")
            end
        end
    end
    isset(conv.domain, "grid_hub") || push!(errs, "&$(conv.domain): required key grid_hub is not set")
    for og in conv.groups
        og.kind == "instance" || continue
        isempty(val(og.name, "backend")) && push!(errs, "&$(og.name): empty backend")
    end
    # named groups exist
    present = Set(g.name for g in conv.groups)
    libnames = Set{String}()
    if conv.keep_library
        for f in conv.library
            p = isabspath(f) ? f : joinpath(conv.yelmox, f)
            isfile(p) && union!(libnames, g.name for g in NL.read_namelist(p))
        end
    end
    function need(owner, key, names)
        for n in (names isa AbstractVector ? names : [names])
            (n isa AbstractString && !isempty(n) && lowercase(n) != "none") || continue
            (n in present || n in libnames) ||
                push!(errs, "&$owner $key = \"$n\": no such group in the file")
        end
    end
    for og in conv.groups, k in og.keys
        if og.kind == "kryos" && k.name in ("globe", "clock", "restart_out")
            need(og.name, k.name, k.value)
        elseif og.kind == "globe" && og.name == conv.names["globe"]
            need(og.name, k.name, k.value)
        elseif og.kind == "domain" && k.name in ("topo", "regions", "ice", "geo", "climate", "surface", "shelf", "geothermal", "sediments", "experiments")
            need(og.name, k.name, k.value)
        elseif og.kind == "domain" && k.name == "out_2D"
            need(og.name, k.name, k.value)
        elseif og.kind == "instance"
            if k.name in ("yelmo", "isos", "bsl", "snapshots", "surface_chion", "chion", "orbit", "marine_shelf", "corr", "Q_geo", "H_sed") ||
               startswith(k.name, "out_") && k.name != "out_vars"
                need(og.name, k.name, k.value)
            elseif k.name == "coefs"
                need(og.name, k.name, ["coef_" * n for n in _names_vec(k.value)])
            elseif k.name == "reference"
                need(og.name, k.name, ["snapshot_" * n for n in _names_vec(k.value)])
            elseif k.name == "terms"
                need(og.name, k.name, ["term_" * n for n in _names_vec(k.value)])
            end
        end
    end
    # the groups the resolution followed: every one written (or kept in the library)
    for (owner, n) in conv.refs
        owner in present || continue
        (n in present || n in libnames) || push!(errs, "&$owner names &$n, which is not in the file")
    end
    # restart keys
    for g in get(get(sp, "restart_none", Dict()), "groups", String[])
        v = val(g, "restart")
        v === nothing || v == "None" || push!(errs, "&$g restart = $(NL.format_value(v)): must be \"None\"")
    end
    isempty(errs) || fail("the kryos par file does not validate:\n  " * join(errs, "\n  "))
    push!(lines, "all $nval validated groups: every key declared; required keys set; schedules on the clock ticks; named groups present")
    return lines
end

_names_vec(v::AbstractVector) = [string(x) for x in v if !isempty(string(x))]
_names_vec(v) = isempty(string(v)) ? String[] : [string(v)]

# ---------------------------------------------------------------------------
# Writing, reading back
# ---------------------------------------------------------------------------

function par_header(c::Conversion)
    io = IOBuffer()
    println(io, "kryos par file, converted from yelmox by KryosTools.Convert v", CONVERTER_VERSION,
            " (tools/bin/yelmox2kryos.jl). Generated: do not edit; change the source and convert again.")
    println(io, "Fully expanded: every key of every group is set; the file needs no defaults and no library.")
    println(io, "yelmox commit: ", c.commit)
    println(io, "source: ", relpath(c.par, c.yelmox))
    println(io, "-p: ", isempty(c.overrides) ? "(none)" : join(c.overrides, " "))
    println(io, "sha256 of the input files the file names:")
    for (f, h) in c.files
        println(io, "  ", rpad(f, 40), " ", h)
    end
    return rstrip(String(take!(io)))
end

function namelist_groups(c::Conversion)
    out = NL.Group[]
    for og in c.groups
        keys = [NL.Key(k.name, k.raw; comment=k.comment) for k in og.keys]
        push!(out, NL.Group(og.name, keys))
    end
    return out
end

"""
    write_par(path, conv)

Write the kryos par file, read it back and check that every key reads back to
the value it was given (and, for a mapped key, to the resolved yelmox value:
type and value). Throws a `ConversionError` on a difference.
"""
function write_par(path::AbstractString, c::Conversion)
    mkpath(dirname(abspath(path)))
    NL.write_namelist(path, namelist_groups(c); header=par_header(c))
    check_roundtrip(c, path)
    return path
end

function check_roundtrip(c::Conversion, path::AbstractString)
    back = NL.group_dict(NL.read_namelist(path))
    errs = String[]
    names = [g.name for g in c.groups]
    length(back) == length(names) || push!(errs, "file has $(length(back)) groups, expected $(length(names))")
    for og in c.groups
        g = get(back, og.name, nothing)
        if g === nothing
            push!(errs, "group &$(og.name) not in the written file")
            continue
        end
        kd = Dict(k.name => k for k in g.keys)
        length(g.keys) == length(og.keys) || push!(errs, "&$(og.name): $(length(g.keys)) keys written, expected $(length(og.keys))")
        for k in og.keys
            r = get(kd, k.name, nothing)
            if r === nothing
                push!(errs, "&$(og.name) $(k.name): not in the written file")
                continue
            end
            v = NL.parse_value(r.raw)
            same_exact(v, k.value) || push!(errs, "&$(og.name) $(k.name): reads back $(repr(v)), expected $(repr(k.value))")
            if k.src !== nothing
                same_exact(v, k.src) || push!(errs, "&$(og.name) $(k.name): reads back $(repr(v)) but yelmox resolved $(repr(k.src))")
            end
        end
    end
    isempty(errs) || fail("round trip of $path failed:\n  " * join(errs, "\n  "))
    push!(c.checks, "round trip: $(sum(length(g.keys) for g in c.groups)) keys of $(length(c.groups)) groups read back equal (type and value), mapped keys equal the resolved yelmox value")
    return true
end

same_exact(a, b) = typeof(a) == typeof(b) && a == b

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

"""
    write_report(path, conv)

The accounting of the conversion in markdown: counts, every dropped key with
its reason, the keys carried, the keys that exist only on the kryos side,
notes, warnings, validation lines, and the files a run directory needs.
"""
function write_report(path::AbstractString, c::Conversion)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        report_markdown(io, c)
    end
    return path
end

mdcell(s) = replace(string(s), "|" => "\\|")

function inputs_needed(c::Conversion)
    # (what, where, why)
    rows = Tuple{String,String}[]
    files = Set{String}(f for (f, _) in c.files)
    for (f, h) in c.files
        push!(rows, ("yelmox `$f`", h == "absent" ? "absent in this yelmox checkout" : "named by the par file or read by a library"))
    end
    comps = unique(og.component for og in c.groups if og.kind == "instance")
    for f in ("kryos_defaults.nml", "forcing_defaults.nml", "snapshots_defaults.nml", "datasrc_defaults.nml")
        push!(rows, ("kryos `input/$f`", "templates (kryos_config, tsgen, snapshots, datasrc)"))
    end
    for k in sort(comps)
        push!(rows, ("kryos `input/kryos-variables-$k.md`", "variable table of the $k component (default writer)"))
    end
    push!(rows, ("kryos `input/phys_const.nml`", "`[kryos] phys_const_file`: the run's constants (phys_const_load, no defaults file)"))
    push!(rows, ("kryos `input/yelmo_phys_const.nml`", "synchronised copy of yelmo's; Yelmo reads `sec_year` from it"))
    return rows
end

function data_paths(c::Conversion)
    pats = Set{String}()
    for og in c.groups, k in og.keys
        vals = k.value isa AbstractVector ? k.value : [k.value]
        for v in vals
            v isa AbstractString || continue
            occursin(r"^(fesmdata|ice_data|isostasy_data)/", v) && push!(pats, v)
        end
    end
    return sort!(collect(pats))
end

function grids(c::Conversion)
    gs = String[]
    for og in c.groups, k in og.keys
        k.name == "grid_hub" || k.name == "grid" || continue
        v = k.value
        v isa AbstractString && !isempty(v) && push!(gs, v)
    end
    return unique(gs)
end

function report_markdown(io::IO, c::Conversion)
    s = account_summary(c)
    println(io, "# Conversion report: ", relpath(c.par, c.yelmox), " -> kryos")
    println(io)
    println(io, "Generated by `KryosTools.Convert` v", CONVERTER_VERSION, ". Do not edit.")
    println(io)
    println(io, "- yelmox: `", c.commit, "`")
    println(io, "- source: `", relpath(c.par, c.yelmox), "`")
    println(io, "- `-p`: ", isempty(c.overrides) ? "(none)" : "`" * join(c.overrides, " ") * "`")
    println(io, "- kryos schema: `input/kryos_defaults.nml`, `snapshots_defaults.nml`, `forcing_defaults.nml`, `datasrc_defaults.nml` of the checkout the converter ran in")
    println(io, "- mode: ", c.keep_library ? "library groups kept as references (`--keep-library`)" : "fully expanded, library groups copied in")
    println(io)
    println(io, "## Accounting")
    println(io)
    println(io, "Every key of every group of the resolved yelmox configuration (the par file with the `-p` overrides, the groups it names, the defaults of the keys it leaves out).")
    println(io)
    println(io, "| status | keys | meaning |")
    println(io, "|---|---:|---|")
    println(io, "| mapped | ", s.mapped, " | carried to the kryos file with the resolved value (or its conversion) |")
    println(io, "| default-equal | ", s.default, " | carried; the yelmox value is the default of its defaults file, or equals the kryos template default |")
    println(io, "| dropped | ", s.dropped, " | not carried, with a reason (below) |")
    println(io, "| total | ", s.total, " | |")
    println(io)
    s.unfilled > 0 && println(io, "Not counted: ", s.unfilled, " keys that yelmox fills from its defaults in dropped groups (", join(sort([g for (g, n) in c.unfilled if n > 0]), ", "), ").\n")

    dropped = [a for a in c.acc if a.status == :dropped]
    println(io, "## Dropped (", length(dropped), ")")
    println(io)
    println(io, "| yelmox key | reason |")
    println(io, "|---|---|")
    for a in dropped
        println(io, "| `", a.group, ".", a.key, "` | ", mdcell(a.reason), " |")
    end
    println(io)
    gone = sort(collect(keys(c.unfilled)))
    println(io, "Groups dropped (all their keys, listed above): ", join(unique(a.group for a in dropped if a.group in gone), ", "), ".")
    println(io)

    println(io, "## Carried")
    println(io)
    println(io, "| yelmox group | keys | mapped | default-equal | kryos destination |")
    println(io, "|---|---:|---:|---:|---|")
    for g in unique(a.group for a in c.acc)
        as = [a for a in c.acc if a.group == g && a.status != :dropped]
        isempty(as) && continue
        dsts = unique([first(split(a.dst, '.')) for a in as if !isempty(a.dst)])
        println(io, "| `", g, "` | ", length(as), " | ", count(a -> a.status == :mapped, as), " | ", count(a -> a.status == :default, as), " | ",
                join(("`&" * replace(d, r"^\(library\) " => "") * "`" for d in dsts), ", "), " |")
    end
    println(io)
    println(io, "Keys with a conversion or a remark:")
    println(io)
    println(io, "| yelmox key | kryos key | remark |")
    println(io, "|---|---|---|")
    for a in c.acc
        a.status == :dropped && continue
        isempty(a.reason) && continue
        println(io, "| `", a.group, ".", a.key, "` | `", a.dst, "` | ", mdcell(a.reason), " |")
    end
    println(io)
    println(io, "Core keys (every `[kryos]`, `[clock]`, schedule, `[globe]`, domain, instance and field group), by origin:")
    println(io)
    println(io, "| group | from yelmox | fixed or rule | template default |")
    println(io, "|---|---:|---:|---:|")
    for og in c.groups
        og.kind == "physics" && continue
        fy = count(k -> startswith(k.origin, "yelmox"), og.keys)
        fx = count(k -> k.origin == "fixed" || startswith(k.origin, "rule"), og.keys)
        println(io, "| `&", og.name, "` | ", fy, " | ", fx, " | ", length(og.keys) - fy - fx, " |")
    end
    println(io)
    only_ = String[]
    for og in c.groups
        og.kind == "physics" && continue
        for k in og.keys
            (k.origin == "fixed" || startswith(k.origin, "rule")) && push!(only_, "`$(og.name).$(k.name)` = $(k.raw)")
        end
    end
    println(io, "Keys set by the converter, not from a yelmox key (fixed by the schema or by a rule): ")
    println(io)
    for x in only_
        println(io, "- ", x)
    end
    println(io)
    println(io, "Physics groups copied unchanged (fully expanded from the yelmox defaults files): ",
            join(("`&" * g.name * "`" for g in c.groups if g.kind == "physics"), ", "), ".")
    println(io)

    if !isempty(c.notes)
        println(io, "## Notes")
        println(io)
        for n in c.notes
            println(io, "- ", n)
        end
        println(io)
    end
    println(io, "## Warnings")
    println(io)
    if isempty(c.warnings)
        println(io, "None.")
    else
        for n in unique(c.warnings)
            println(io, "- ", n)
        end
    end
    println(io)
    println(io, "## Checks")
    println(io)
    println(io, "- refusals (§2.5): none of the phase 1 refusals applies")
    println(io, "- accounting: every one of the ", s.total, " keys is mapped, default-equal or dropped with a reason")
    for l in c.checks
        println(io, "- ", l)
    end
    println(io)
    println(io, "## Files a run directory needs")
    println(io)
    println(io, "Run directory layout: the par file, and links `input`, `fesmdata`, `ice_data`, `isostasy_data`, `maps`, `fesm-utils`. `input` must hold the kryos and the yelmox input files together (`Gate.input_overlay` builds it).")
    println(io)
    println(io, "| file | why |")
    println(io, "|---|---|")
    for (a, b) in inputs_needed(c)
        println(io, "| ", a, " | ", b, " |")
    end
    println(io)
    gs = grids(c)
    println(io, "Grid descriptions in `maps/`: ", join(("`grid_$g.txt`" for g in gs), ", "), "; the cached maps between them (frozen map cache of the gate).")
    println(io)
    dp = data_paths(c)
    if !isempty(dp)
        println(io, "Data paths read (`{domain}`, `{grid_name}` resolved by kryos):")
        println(io)
        for p in dp
            println(io, "- `", p, "`")
        end
        println(io)
    end
    println(io, "## Input file hashes")
    println(io)
    for (f, h) in c.files
        println(io, "- `", f, "` ", h)
    end
    return nothing
end

"""
    convert_file(; yelmox, par, out, report="", overrides=[], kryos, keep_library=false) -> Conversion

`convert_yelmox`, then `write_par(out, ...)` and `write_report(report, ...)`.
"""
function convert_file(; out::AbstractString, report::AbstractString="", kw...)
    conv = convert_yelmox(; kw...)
    write_par(out, conv)
    isempty(report) || write_report(report, conv)
    return conv
end

# ---------------------------------------------------------------------------
# Merging converted single-domain files into one N-domain file (phase 3 WP1)
# ---------------------------------------------------------------------------

"""
A merged par file: the groups, the domains in file order, and the accounting of
the merge (`renamed`: domain, old, new; `repointed`: group, key, old, new;
`resolved`: `group.key`, the value of each file, the value taken; `runwide`:
the run-wide groups, `sources`: the files).
"""
struct Merged
    groups::Vector{NL.Group}
    domains::Vector{String}
    suffixes::Vector{String}
    sources::Vector{String}
    renamed::Vector{NTuple{3,String}}
    repointed::Vector{NTuple{4,String}}
    resolved::Vector{Tuple{String,Vector{String},String}}
    runwide::Vector{String}
    header::String
end

"""
The groups every domain of a run shares (taken once, must agree between the
files): the run, clock, globe and schedule groups, the global instances, the
bsl and the forcing coefficients.
"""
is_runwide(name::AbstractString, wide::AbstractVector{<:AbstractString}) =
    name in wide || startswith(name, "out_") || startswith(name, "coef_")

# Keys that name a group, as (prefix of the group name, keys). The value `v` of a key with prefix
# `p` names the group `p*v`. Plain: the domain's component keys, the instance keys that name
# backend groups, the `nml_*` keys of the physics groups and the field groups of a snapshot.
const MERGE_REF_KEYS = [
    "" => ["topo", "regions", "ice", "geo", "climate", "surface", "shelf", "geothermal", "sediments", "experiments",
           "yelmo", "isos", "snapshots", "taper", "bias", "surface_chion", "chion", "orbit", "marine_shelf", "corr",
           "Q_geo", "H_sed"],
    "term_" => ["terms"],
    "snapshot_" => ["reference", "snap0", "snap1"],
]

"The prefix with which the value of key `key` of group `group` names a group, or `nothing` (not a reference)."
function ref_prefix(group::AbstractString, key::AbstractString)
    for (p, ks) in MERGE_REF_KEYS
        key in ks && return p
    end
    startswith(key, "nml_") && return ""
    startswith(group, "snapshot_") && key in ("tas", "pr", "zs", "to", "so") && return ""
    return nothing
end

"`raw` with the quoted element `from` replaced by `to` (either quote character)."
requote(raw::AbstractString, from::AbstractString, to::AbstractString) =
    replace(raw, "\"$from\"" => "\"$to\"", "'$from'" => "'$to'")

"The quoted string elements of a key's raw value."
function str_elements(raw::AbstractString)
    els = NL.elements(raw)
    return [t for (t, q) in els if q]
end

"""
The suffix of a converted single-domain file: what the instance names of its domain group share after the
component (`ice_nor`, `geo_nor`, ... give `nor`).
"""
function domain_suffix(dg::NL.Group, file::AbstractString)
    sfx = Set{String}()
    for k in dg.keys
        k.name in ("ice", "geo", "climate", "surface", "shelf", "geothermal", "sediments") || continue
        v = NL.parse_value(k.raw)
        (v isa AbstractString && !isempty(v)) || continue
        i = findfirst('_', v)
        i === nothing && fail("$file: instance $(k.name) = \"$v\" of &$(dg.name) has no `_<suffix>`: cannot name the domain's groups")
        push!(sfx, v[i+1:end])
    end
    length(sfx) == 1 || fail("$file: the instances of &$(dg.name) do not share one suffix ($(join(sort(collect(sfx)), ", "))): " *
                              "cannot name the domain's groups")
    return only(sfx)
end

kv_of(g::NL.Group) = Dict(k.name => k for k in g.keys)

"""
    merge_par(files; resolve=Dict(), wide=[...]) -> Merged

Merge converted single-domain par files into one file for a run with all their domains
(`[kryos] domains` in the order of `files`).

- The run-wide groups (`[kryos]`, `[clock]`, `[globe]`, `[forcing]`, `[geo_glob]`, `[barysealevel]`,
  `[out_*]`, `[coef_*]`) are taken once and must agree: a key with different values in two files
  stops the merge, naming both, unless `resolve` gives its value (`"group.key" => raw text`, recorded in the
  header). `[kryos] domains` lists the domains and `[forcing] coefs` the union of the coefficients.
- Every other group belongs to its file's domain. A name that two files both have is renamed with
  the domain's suffix (`ytopo` -> `ytopo_nor`; a snapshot `pd` -> `pd_nor`, i.e. `snapshot_pd_nor`) and the keys that name
  the group are repointed (`MERGE_REF_KEYS`); a name only one file has is kept.
- Every file must be a single-domain file with the groups of its domain (the instance groups, the groups they name).
"""
function merge_par(files::AbstractVector{<:AbstractString}; resolve::AbstractDict=Dict{String,String}(),
                   wide::AbstractVector{<:AbstractString}=["kryos", "clock", "globe", "forcing", "geo_glob", "barysealevel"])
    length(files) >= 2 || fail("merge_par: need at least two files")
    res = Dict{String,String}(String(k) => String(v) for (k, v) in resolve)
    used_res = Set{String}()
    parsed = [NL.read_namelist(f) for f in files]
    gds = [NL.group_dict(g) for g in parsed]
    doms = String[]; sfxs = String[]
    for (i, f) in enumerate(files)
        haskey(gds[i], "kryos") || fail("$f: no &kryos group")
        d = NL.parse_value(kv_of(gds[i]["kryos"])["domains"].raw)
        (d isa AbstractString && !isempty(d)) || fail("$f: [kryos] domains = $d: a single-domain file is expected")
        haskey(gds[i], d) || fail("$f: the domain group &$d is not in the file")
        push!(doms, d)
        push!(sfxs, domain_suffix(gds[i][d], f))
    end
    length(unique(doms)) == length(doms) || fail("the files have the same domain ($(join(doms, ", "))): $(join(files, ", "))")
    length(unique(sfxs)) == length(sfxs) || fail("the domains $(join(doms, ", ")) have the same instance suffix ($(join(sfxs, ", ")))")

    # run-wide groups: in the order of the first file, then the new groups of the others
    wnames = String[]
    for g in parsed, grp in g
        is_runwide(grp.name, wide) && !(grp.name in wnames) && push!(wnames, grp.name)
    end
    out = NL.Group[]
    resolved = Tuple{String,Vector{String},String}[]
    conflicts = String[]
    for name in wnames
        have = [(i, gds[i][name]) for i in eachindex(files) if haskey(gds[i], name)]
        keys = NL.Key[]
        names = String[]
        for (_, g) in have, k in g.keys
            k.name in names || push!(names, k.name)
        end
        for kn in names
            vals = [(i, get(kv_of(g), kn, nothing)) for (i, g) in have]
            if name == "kryos" && kn == "domains"
                push!(keys, NL.Key(kn, join(["\"$d\"" for d in doms], " "); comment=join(["domain of " * basename(f) for f in files], ", ")))
                continue
            end
            if name == "forcing" && kn == "coefs"
                cs = String[]
                for (_, k) in vals
                    k === nothing && continue
                    for c in str_elements(k.raw)
                        isempty(c) || c in cs || push!(cs, c)
                    end
                end
                raw = isempty(cs) ? "\"\"" : join(["\"$c\"" for c in cs], " ")
                push!(keys, NL.Key(kn, raw; comment="the union of the coefficients of the files"))
                continue
            end
            present = [(i, k) for (i, k) in vals if k !== nothing]
            first_k = present[1][2]
            agree = all(NL.parse_value(k.raw) == NL.parse_value(first_k.raw) for (_, k) in present) && length(present) == length(have)
            if agree
                push!(keys, first_k)
                continue
            end
            id = "$name.$kn"
            if haskey(res, id)
                push!(used_res, id)
                push!(resolved, (id, String[strip(string(k === nothing ? "(absent)" : k.raw)) for (_, k) in vals], res[id]))
                push!(keys, NL.Key(kn, res[id]; comment="resolved: " * join([basename(files[i]) * " " * (k === nothing ? "(absent)" : strip(k.raw)) for (i, k) in vals], ", ")))
            else
                push!(conflicts, "&$name $kn (resolve[\"$id\"]):\n" *
                      join(["    $(files[i]): " * (k === nothing ? "(absent)" : strip(k.raw)) for (i, k) in vals], "\n"))
                push!(keys, first_k)
            end
        end
        push!(out, NL.Group(name, keys; comment=first(have)[2].comment))
    end
    isempty(conflicts) || fail("run-wide keys that differ between the files and are not resolved:\n" * join(conflicts, "\n"))
    for id in keys(res)
        id in used_res || fail("resolve $id = $(res[id]): the files agree on it or do not have it (remove it)")
    end

    # domain-owned groups: the names two files share are renamed
    owned = [[g for g in parsed[i] if !is_runwide(g.name, wide)] for i in eachindex(files)]
    ownedn = [Set(g.name for g in o) for o in owned]
    renamed = NTuple{3,String}[]
    maps = [Dict{String,String}() for _ in files]
    for i in eachindex(files), g in owned[i]
        if any(j -> j != i && g.name in ownedn[j], eachindex(files))
            new = g.name * "_" * sfxs[i]
            new in ownedn[i] && fail("$(files[i]): cannot rename &$(g.name) to &$new: the file has it")
            maps[i][g.name] = new
            push!(renamed, (doms[i], g.name, new))
        end
    end
    repointed = NTuple{4,String}[]
    refd = [Set{String}() for _ in files]
    for i in eachindex(files)
        for g in owned[i]
            keys = NL.Key[]
            for k in g.keys
                p = ref_prefix(g.name, k.name)
                raw = k.raw
                if p !== nothing
                    for e in unique(str_elements(raw))
                        isempty(e) && continue
                        t = p * e
                        t in ownedn[i] && push!(refd[i], t)
                        if haskey(maps[i], t)
                            new = e * "_" * sfxs[i]
                            raw = requote(raw, e, new)
                            push!(repointed, (get(maps[i], g.name, g.name), k.name, e, new))
                        end
                    end
                end
                push!(keys, NL.Key(k.name, raw; comment=k.comment, pre=k.pre))
            end
            name = get(maps[i], g.name, g.name)
            push!(out, NL.Group(name, keys; comment=g.comment, inline=g.inline, post=g.post))
        end
    end
    # a renamed group no key names would be a group the other domain's keys silently cannot reach
    for i in eachindex(files), (old, _) in maps[i]
        old in refd[i] || fail("$(files[i]): the renamed group &$old is named by no key (MERGE_REF_KEYS has no entry for the key that names it)")
    end
    names = [g.name for g in out]
    length(unique(names)) == length(names) || fail("merged file: duplicate groups " * join(unique([n for n in names if count(==(n), names) > 1]), ", "))

    hdr = IOBuffer()
    println(hdr, "kryos par file, merged from ", length(files), " single-domain files by KryosTools.Convert.merge_par (tools/bin/kryos-merge.jl). Generated: do not edit.")
    println(hdr, "domains: ", join(doms, ", "), " (instance suffixes ", join(sfxs, ", "), ")")
    for (i, f) in enumerate(files)
        println(hdr, "from ", doms[i], ": ", basename(f))
        for g in parsed[i][1:1], l in split(g.pre, '\n')
            occursin(r"^(yelmox commit|source|-p):", strip(l)) && println(hdr, "  ", strip(l))
        end
    end
    for (id, vals, v) in resolved
        println(hdr, "resolved ", id, " = ", v, "  (", join(vals, " | "), ")")
    end
    println(hdr, "domain-specific groups the files share are renamed with the instance suffix (see the merge report).")
    return Merged(out, doms, sfxs, String.(files), renamed, repointed, resolved, wnames, rstrip(String(take!(hdr))))
end

"Write a merged file (the header names the sources and the resolutions)."
function write_merged(path::AbstractString, m::Merged)
    mkpath(dirname(abspath(path)))
    gs = copy(m.groups)
    NL.write_namelist(path, gs; header=m.header)
    back = NL.read_namelist(path)
    [g.name for g in back] == [g.name for g in gs] || fail("merged file $path does not read back to its groups")
    return path
end

"""
    split_merged(m::Merged, i) -> Vector{NL.Group}

The single-domain file of domain `i` inside a merged one (the inverse of `merge_par`, for the checks): the run-wide
groups and the groups reachable from the domain group through the keys of `MERGE_REF_KEYS`, with the renamed
groups and the keys that name them restored. Equals the original file's groups except for the run-wide values the other
domain contributed (`[kryos] domains`, the union of `[forcing] coefs`, resolved keys).
"""
function split_merged(m::Merged, i::Integer)
    gd = NL.group_dict(m.groups)
    sfx = m.suffixes[i]
    undo = Dict(r[3] => r[2] for r in m.renamed if r[1] == m.domains[i])
    domain = m.domains[i]
    reach = String[domain]
    queue = [domain]
    while !isempty(queue)
        g = gd[popfirst!(queue)]
        for k in g.keys
            p = ref_prefix(g.name, k.name)
            p === nothing && continue
            for e in str_elements(k.raw)
                isempty(e) && continue
                t = p * e
                (haskey(gd, t) && !is_runwide(t, m.runwide) && !(t in reach)) || continue
                push!(reach, t); push!(queue, t)
            end
        end
    end
    out = NL.Group[g for g in m.groups if is_runwide(g.name, m.runwide)]
    for n in reach
        g = gd[n]
        keys = NL.Key[]
        for k in g.keys
            p = ref_prefix(g.name, k.name)
            raw = k.raw
            if p !== nothing
                for e in unique(str_elements(raw))
                    haskey(undo, p * e) || continue
                    old = undo[p * e][length(p)+1:end]
                    raw = requote(raw, e, old)
                end
            end
            push!(keys, NL.Key(k.name, raw; comment=k.comment, pre=k.pre))
        end
        push!(out, NL.Group(get(undo, n, n), keys; comment=g.comment, inline=g.inline, post=g.post))
    end
    return out
end

"The merge as a markdown report: sources, run-wide groups, renamed groups, repointed keys, resolutions."
function merge_report(io::IO, m::Merged, out::AbstractString="")
    println(io, "# Merged par file", isempty(out) ? "" : " `$(basename(out))`")
    println(io)
    println(io, "Sources (domain, instance suffix):")
    for (f, d, s) in zip(m.sources, m.domains, m.suffixes)
        println(io, "- `", f, "`: ", d, " (`", s, "`)")
    end
    println(io)
    println(io, "Run-wide groups (taken once, the files agree except where resolved): ", join(["`&$n`" for n in m.runwide], ", "), ".")
    println(io)
    if !isempty(m.resolved)
        println(io, "## Resolved keys")
        println(io)
        println(io, "| key | values of the files | taken |\n|---|---|---|")
        for (id, vals, v) in m.resolved
            println(io, "| `", id, "` | ", join(["`$x`" for x in vals], ", "), " | `", v, "` |")
        end
        println(io)
    end
    println(io, "## Renamed groups (", length(m.renamed), ")")
    println(io)
    for d in m.domains
        rs = [r for r in m.renamed if r[1] == d]
        println(io, "- ", d, ": ", isempty(rs) ? "none" : join(["`&$(r[2])` -> `&$(r[3])`" for r in rs], ", "))
    end
    println(io)
    println(io, "## Repointed keys (", length(m.repointed), ")")
    println(io)
    println(io, "| group | key | was | now |\n|---|---|---|---|")
    for r in m.repointed
        println(io, "| `&", r[1], "` | `", r[2], "` | `", r[3], "` | `", r[4], "` |")
    end
    return nothing
end

end # module
