"""
    KryosTools.Gate

What the gate automation (regression-design.md §4, `tools/gate/gate.sh`) needs
from Julia, as functions:

- the stack file (§4.1): `load_stack(path)`, `stack_problems`, `stack_id`
  (the first 8 hex digits of the sha256 of the stack's canonical text, without
  the `kryos` table: the kryos commit is in the run directory name), `stack_text`,
  `stack_commits`;
- the gate's `-p` lists (§1.1) and chunk overrides (§1.2): `plist`,
  `plist_string`, `chunk_overrides(gate, :yelmox | :kryos, k; prev_dir)`,
  `restart_dir(model, t)`;
- the gate's directories (§1.4): `gate_root`, `ref_dir`, `run_dir`,
  `maps_dir`, `stack_dir`, `report_paths`;
- the run directory's `input` (kryos and yelmox input files together):
  `input_overlay(dest; kryos, yelmox)`;
- consistency of a gate file with a converted par file: `check_instances`;
- the `docs/gate.md` row (§4.4): `gate_row`, `append_gate_row`, `GATE_ROW_HEADER`;
- the kryos-against-kryos gates (phase 3, regression-design.md §1.6, `kind = "kryos"` in the gate file): the runs
  to stage (`kk_runs`, `kk_stage`), the checks (`kk_checks`, `run_kk_checks`: per-domain equality of a two-domain
  run with single-domain runs, identity, chunk == cont, the `bsl` budget `bsl_budget`), `kk_markdown`.

The gate files are `tools/gate/<gate>.toml` (read by `Compare.load_gate`) and
`tools/gate/stack.toml`.
"""
module Gate

using TOML
using SHA
using Printf
using Dates
using NCDatasets
import FesmUtils.Namelist as NL
using ..Compare

export Stack, load_stack, stack_problems, stack_text, stack_id, stack_commits, DEFAULT_STACK
export plist, plist_string, chunk_overrides, restart_dir
export gate_root, stack_dir, ref_dir, run_dir, maps_dir, report_paths
export input_overlay, check_instances
export gate_row, append_gate_row, GATE_ROW_HEADER
export KKResult, Budget, bsl_budget, kk_runs, kk_stage, kk_checks, run_kk_checks, kk_markdown, kk_outcome

const GATE_DIR = normpath(joinpath(@__DIR__, "..", "gate"))
const DEFAULT_STACK = joinpath(GATE_DIR, "stack.toml")

# ---------------------------------------------------------------------------
# The stack (§4.1)
# ---------------------------------------------------------------------------

"A parsed `stack.toml` and the file it came from."
struct Stack
    data::Dict{String,Any}
    path::String
end

const STACK_TABLES = ("yelmox", "libs", "tools", "build", "run")
const STACK_LIBS = ("fesm-utils", "yelmo", "FastHydrology", "elsa", "tracer", "FastIsostasy", "chion", "insol", "rembo1")

"""
    load_stack(path=DEFAULT_STACK) -> Stack

Read the stack file and stop on a missing table or field (`stack_problems`
without `full`).
"""
function load_stack(path::AbstractString=DEFAULT_STACK)
    isfile(path) || error("no stack file $path")
    s = Stack(TOML.parsefile(path), String(path))
    bad = stack_problems(s)
    isempty(bad) || error("stack file $path:\n  " * join(bad, "\n  "))
    return s
end

"""
    stack_problems(stack; full=false) -> Vector{String}

Missing tables and fields (yelmox ref and commit, every library's ref and
commit, configme and runme versions, compiler, modules, flags,
`omp_threads`, partition, account). With `full`, also commits that are not a
full 40-digit sha (a pinned stack is resolved with `gate.sh stack`, which writes
`stack.lock.json`).
"""
function stack_problems(s::Stack; full::Bool=false)
    d = s.data
    bad = String[]
    for t in STACK_TABLES
        haskey(d, t) || push!(bad, "no [$t] table")
    end
    isempty(bad) || return bad
    for k in ("ref", "commit")
        haskey(d["yelmox"], k) || push!(bad, "[yelmox] has no $k")
    end
    for lib in STACK_LIBS
        l = get(d["libs"], lib, nothing)
        if l === nothing
            push!(bad, "[libs.$lib] missing")
            continue
        end
        for k in ("ref", "commit")
            haskey(l, k) || push!(bad, "[libs.$lib] has no $k")
        end
    end
    for k in ("configme", "runme")
        haskey(d["tools"], k) || push!(bad, "[tools] has no $k")
    end
    for k in ("machine", "compiler", "modules", "flags")
        haskey(d["build"], k) || push!(bad, "[build] has no $k")
    end
    for k in ("omp_threads", "partition", "account", "walltime")
        haskey(d["run"], k) || push!(bad, "[run] has no $k")
    end
    if full && isempty(bad)
        for (name, c) in stack_commits(s)
            occursin(r"^[0-9a-f]{40}$", c) || push!(bad, "$name: commit $c is not a full sha")
        end
    end
    return bad
end

"""
    stack_commits(stack) -> Vector{Pair{String,String}}

`name => commit` for yelmox and every library, in the order of the file.
"""
function stack_commits(s::Stack)
    out = Pair{String,String}["yelmox" => String(s.data["yelmox"]["commit"])]
    for lib in STACK_LIBS
        push!(out, lib => String(s.data["libs"][lib]["commit"]))
    end
    return out
end

"""
    stack_text(stack) -> String

The canonical text the stack id is taken from: the tables of the stack with
sorted keys, without the `kryos` table (its commit is in the run directory
name, `kx-<sha7>`) and without `note`.
"""
function stack_text(s::Stack)
    d = Dict{String,Any}(k => v for (k, v) in s.data if k in STACK_TABLES)
    io = IOBuffer()
    TOML.print(io, d; sorted=true)
    return String(take!(io))
end

"The stack id: the first 8 hex digits of the sha256 of `stack_text` (`v24` is the id of the R0 stack)."
stack_id(s::Stack) = bytes2hex(sha256(stack_text(s)))[1:8]

# ---------------------------------------------------------------------------
# -p lists and chunks (§1.1, §1.2)
# ---------------------------------------------------------------------------

"The `-p` overrides of a gate (a gate file as `Compare.load_gate` returns it)."
plist(g::AbstractDict) = String.(get(g, "p", String[]))

"The `-p` list as one string, as it is passed to `runme -p`."
plist_string(g::AbstractDict) = join(plist(g), " ")

"""
    restart_dir(model, t) -> String

The restart bundle folder a run writes at time `t` [yr]: yelmox
`restart-<t/1000, 3 decimals>-kyr`, kryos `restart-<t>` (the shortest exact
decimal).
"""
function restart_dir(model::Symbol, t::Real)
    if model == :yelmox
        return @sprintf("restart-%.3f-kyr", t / 1000)
    elseif model == :kryos
        return "restart-" * (isinteger(t) ? string(Int(round(t))) : string(Float64(t)))
    end
    error("restart_dir: model $model (yelmox or kryos)")
end

fmt_t(t::Real) = isinteger(t) ? string(Int(round(t))) : string(Float64(t))

"""
    chunk_overrides(gate, model, k; prev_dir="") -> Vector{String}

The `-p` overrides of chunk `k` (0-based) of the gate's `chunks = [t0, t1, ...]`
on top of the base: `ctrl.time_init`, `ctrl.time_end` (yelmox) or
`clock.time_init`, `clock.time_end` (kryos), and for `k > 0` the restart bundle
of the previous chunk, `ctrl.restart` / `kryos.restart` = `prev_dir/<bundle>`
(`prev_dir`: the absolute run folder of chunk `k-1`). The first chunk of a gate
whose base has `time_init = t0` still sets both times, as `run_chunks.sh` does.
"""
function chunk_overrides(g::AbstractDict, model::Symbol, k::Integer; prev_dir::AbstractString="")
    ts = Float64.(get(g, "chunks", Float64[]))
    length(ts) >= 3 || error("gate $(get(g, "gate", "?")): chunks $ts (at least two chunks)")
    0 <= k <= length(ts) - 2 || error("chunk $k: the gate has $(length(ts) - 1) chunks")
    grp = model == :yelmox ? "ctrl" : model == :kryos ? "clock" : error("chunk_overrides: model $model")
    p = ["$grp.time_init=$(fmt_t(ts[k+1]))", "$grp.time_end=$(fmt_t(ts[k+2]))"]
    if k > 0
        isempty(prev_dir) && error("chunk $k needs prev_dir (the run folder of chunk $(k - 1))")
        key = model == :yelmox ? "ctrl.restart" : "kryos.restart"
        push!(p, "$key=$(joinpath(prev_dir, restart_dir(model, ts[k+1])))")
    end
    return p
end

# ---------------------------------------------------------------------------
# Directories (§1.4)
# ---------------------------------------------------------------------------

"The gate root `\$G`: `ENV[\"GATE_ROOT\"]`, else `/work/ba1442/robinson/kryos-gate`."
gate_root() = get(ENV, "GATE_ROOT", "/work/ba1442/robinson/kryos-gate")

sha7(s::AbstractString) = (length(s) >= 7 || error("sha $s: need at least 7 characters"); String(first(s, 7)))

"`\$G/stacks/<stack>`: the yelmox and kryos installs and the shared dependency checkouts."
stack_dir(stack::AbstractString; root::AbstractString=gate_root()) = joinpath(root, "stacks", stack)

"`\$G/maps/<stack>`: the frozen map cache."
maps_dir(stack::AbstractString; root::AbstractString=gate_root()) = joinpath(root, "maps", stack)

"""
    ref_dir(gate, stack, yelmox_sha, run; root) -> String

`\$G/ref/<gate>/yx-<sha7>/<stack>/<run>`, a yelmox run (`cont`, `cont2` or `chunk`).
"""
function ref_dir(gate::AbstractString, stack::AbstractString, yelmox_sha::AbstractString, run::AbstractString;
                 root::AbstractString=gate_root())
    run in ("cont", "cont2", "chunk") || error("ref_dir: run $run (cont, cont2 or chunk)")
    return joinpath(root, "ref", gate, "yx-" * sha7(yelmox_sha), stack, run)
end

"""
    run_dir(gate, stack, kryos_sha, run; root) -> String

`\$G/runs/<gate>/kx-<sha7>/<stack>/<run>`, a kryos run (`cont` or `chunk`).
"""
function run_dir(gate::AbstractString, stack::AbstractString, kryos_sha::AbstractString, run::AbstractString;
                 root::AbstractString=gate_root())
    run in ("cont", "chunk") || error("run_dir: run $run (cont or chunk)")
    return joinpath(root, "runs", gate, "kx-" * sha7(kryos_sha), stack, run)
end

"""
    report_paths(kryos_sha, stack; date=today(), root) -> (md, json)

`\$G/reports/<date>_kx-<sha7>_<stack>.{md,json}`.
"""
function report_paths(kryos_sha::AbstractString, stack::AbstractString; date::Union{Date,AbstractString}=today(),
                      root::AbstractString=gate_root())
    base = joinpath(root, "reports", "$(date)_kx-$(sha7(kryos_sha))_$stack")
    return (md=base * ".md", json=base * ".json")
end

# ---------------------------------------------------------------------------
# The run directory's input
# ---------------------------------------------------------------------------

"""
    input_overlay(dest; kryos, yelmox) -> Dict{String,String}

Make the folder `dest` (the `input` link target of a kryos run directory) of
symbolic links to every entry of `yelmox/input` and then of `kryos/input`
(kryos wins a name that both have). kryos reads its templates
(`kryos_defaults.nml`, `forcing_defaults.nml`, `snapshots_defaults.nml`,
`kryos-variables-*.md`), the synchronised library copies (`datasrc_defaults.nml`,
`yelmo_defaults.nml`, `chion_defaults.nml`, ...; `input/sync.txt`), the
phys_const file and the series from `input/`.
Returns name => link target. An existing `dest` is emptied of links first.
"""
function input_overlay(dest::AbstractString; kryos::AbstractString, yelmox::AbstractString)
    mkpath(dest)
    for f in readdir(dest; join=true)
        islink(f) && rm(f)
    end
    out = Dict{String,String}()
    link(src) = begin
        name = basename(src)
        t = joinpath(dest, name)
        islink(t) && rm(t)
        ispath(t) && error("input_overlay: $t exists and is not a link")
        symlink(abspath(src), t)
        out[name] = abspath(src)
    end
    for root in (joinpath(yelmox, "input"), joinpath(kryos, "input"))
        isdir(root) || error("input_overlay: no folder $root")
        for n in sort(readdir(root))
            link(joinpath(root, n))
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# A gate file against a converted par file
# ---------------------------------------------------------------------------

"""
    check_instances(gate, par) -> Vector{String}

The instance names of the gate file (`[instances]`, the placeholders of the
comparison map) against the kryos par file: the domain group, the instance each
component key of it names, and the forcing and global geo groups of `[globe]`.
Returns the differences (empty when they agree).
"""
function check_instances(g::AbstractDict, par::AbstractString)
    inst = Dict{String,String}(k => String(v) for (k, v) in get(g, "instances", Dict()))
    groups = NL.group_dict(NL.read_namelist(par))
    bad = String[]
    dom = get(inst, "domain", "")
    if !haskey(groups, dom)
        return ["domain group &$dom is not in $par"]
    end
    dk = Dict(k.name => NL.parse_value(k.raw) for k in groups[dom].keys)
    for c in ("ice", "geo", "climate", "surface", "shelf", "geothermal", "sediments")
        want = get(inst, c, "")
        have = get(dk, c, "")
        want == have || push!(bad, "[instances] $c = \"$want\" but &$dom $c = \"$have\"")
        haskey(groups, want) || push!(bad, "instance group &$want is not in the par file")
    end
    gl = get(groups, "globe", nothing)
    if gl === nothing
        push!(bad, "no &globe group")
    else
        gk = Dict(k.name => NL.parse_value(k.raw) for k in gl.keys)
        for (c, key) in (("forcing", "forcing"), ("geo_glob", "geo"))
            want = get(inst, c, "")
            have = get(gk, key, "")
            want == have || push!(bad, "[instances] $c = \"$want\" but &globe $key = \"$have\"")
            haskey(groups, want) || push!(bad, "instance group &$want is not in the par file")
        end
    end
    return bad
end

# ---------------------------------------------------------------------------
# The docs/gate.md row (§4.4)
# ---------------------------------------------------------------------------

const GATE_ROW_COLUMNS = ["date", "gate", "kryos", "yelmox", "stack", "flags", "length", "outcome",
                          "first divergence", "deviations", "report"]

"The header of the table of `docs/gate.md` (two lines)."
const GATE_ROW_HEADER = "| " * join(GATE_ROW_COLUMNS, " | ") * " |\n|" * join(fill("---", length(GATE_ROW_COLUMNS)), "|") * "|"

cell(s) = replace(replace(string(s), "|" => "\\|"), '\n' => ' ')

"""
    gate_row(; date, gate, kryos, yelmox, stack, flags, length, outcome,
             first_divergence="", deviations=[], report) -> String

One table row of `docs/gate.md`: date, gate, kryos sha, yelmox sha, stack id
(`stack` may be `id => commits` with the commits it resolves to), flags, length
[yr], outcome (`PASS`, ...), first divergence (if any), deviations used and
the report path on Levante.
"""
function gate_row(; date=today(), gate, kryos, yelmox, stack, flags, length, outcome, first_divergence="",
                  deviations=String[], report)
    st = stack isa Pair ? "$(first(stack)) ($(last(stack)))" : stack
    dev = deviations isa AbstractString ? deviations : (isempty(deviations) ? "none" : join(deviations, ", "))
    fd = isempty(string(first_divergence)) ? "none" : first_divergence
    vals = [date, gate, sha7(kryos), sha7(yelmox), st, flags, length isa Real ? "$(fmt_t(length)) yr" : length, outcome, fd, dev, report]
    return "| " * join(cell.(vals), " | ") * " |"
end

"""
    gate_row(reports; kw...) -> String

The row of a gate from its comparison reports (`Compare.Report`s: `K` vs `R`,
`chunk` vs `cont`, ...): outcome PASS when every report passes (else the worst
of FAIL and ERROR), the earliest first divergence, the union of the deviations
used. `kw` are the remaining fields of the keyword method.
"""
function gate_row(reps::AbstractVector{<:Compare.Report}; kw...)
    outs = [Compare.outcome(r) for r in reps]
    outcome = all(==("PASS"), outs) ? "PASS" : "ERROR" in outs ? "ERROR" : "FAIL"
    fds = [(r.first_divergence["time"], r.first_divergence) for r in reps if r.first_divergence !== nothing]
    fd = ""
    if !isempty(fds)
        t, d = fds[argmin([x[1] for x in fds])]
        vars = [v["key"] for v in d["variables"]]
        fd = (t == -Inf ? "static" : "t = $(fmt_t(t))") * ": " * join(first(vars, 3), ", ") * (Base.length(vars) > 3 ? ", ..." : "")
    end
    devs = unique(reduce(vcat, [r.deviations_used for r in reps]; init=String[]))
    return gate_row(; outcome=outcome, first_divergence=fd, deviations=devs, kw...)
end

"""
    append_gate_row(path, row)

Append `row` to the table of `docs/gate.md` (the file is created with a title
and `GATE_ROW_HEADER` if it does not exist).
"""
function append_gate_row(path::AbstractString, row::AbstractString)
    if !isfile(path)
        mkpath(dirname(abspath(path)))
        open(path, "w") do io
            println(io, "# Gate log")
            println(io)
            println(io, "One row per official gate run (regression-design.md §4.4); a phase is signed off by a row with both gates green.")
            println(io)
            println(io, GATE_ROW_HEADER)
        end
    end
    open(path, "a") do io
        println(io, row)
    end
    return path
end

# ---------------------------------------------------------------------------
# Kryos-against-kryos gates (phase 3)
# ---------------------------------------------------------------------------

"""
A kryos-against-kryos gate file (`kind = "kryos"`) names its runs and its checks:

    [domains.<Domain>.instances]      the instance names of each domain (as [instances] of a one-domain gate)
    skip = ["<file>:<var>", ...]      globs of variables not compared (aggregates over the domains)

    [[run]]                           what to stage (`kk_stage`)
    name = "cont"                     the run's folder
    par  = "sci"                      par/gate/<par>.nml
    kp   = ['barysealevel.method="const"']   overrides of the kryos par file (merged group names)
    swap = true                       [kryos] domains in the reverse order

    [[check]]                         what to compare (`run_kk_checks`; the runs are folders of one directory)
    type = "domain"    domain = "North" ref = "nor" run = "cont"     a two-domain run against a single-domain run
    type = "identical" ref = "cont" run = "cont2"                    two runs of the same configuration
    type = "chunks"    ref = "cont" run = "chunk"                    chunk == cont (c00, c01, ...)
    type = "budget"    run = "cont"                                  the `bsl` budget closes

`restarts = false` leaves out the restart bundles of a domain or identical check.
"""
function kk_runs(g::AbstractDict)
    get(g, "kind", "") == "kryos" || error("gate $(get(g, "gate", "?")): not a kryos-against-kryos gate (kind = \"kryos\")")
    return [Dict{String,Any}(r) for r in get(g, "run", [])]
end

"""
    kk_stage(g, name; domains=nothing) -> (par, overrides)

The par file (`par/gate/<par>.nml`, without the folder) and the kryos `-p` overrides of run `name`: the run's `kp`,
and for `swap = true` `kryos.domains` with the domains of `domains` (default: the gate's `[domains]`, in file order)
reversed. Chunk overrides (`chunk_overrides`) come on top.
"""
function kk_stage(g::AbstractDict, name::AbstractString; domains::Union{Nothing,AbstractVector{<:AbstractString}}=nothing)
    runs = kk_runs(g)
    i = findfirst(r -> r["name"] == name, runs)
    i === nothing && error("gate $(g["gate"]): no run $name (", join([r["name"] for r in runs], ", "), ")")
    r = runs[i]
    p = String.(get(r, "kp", String[]))
    if get(r, "swap", false)
        ds = domains === nothing ? collect(String, get(g, "domain_order", String[])) : String.(domains)
        isempty(ds) && error("gate $(g["gate"]): run $name swaps the domains but the gate has no `domain_order`")
        push!(p, "kryos.domains=" * join(["\"$d\"" for d in reverse(ds)], " "))
    end
    return (String(r["par"]), p)
end

kk_checks(g::AbstractDict) = [Dict{String,Any}(c) for c in get(g, "check", [])]

"The result of a kryos-against-kryos check: `outcome` PASS, FAIL or ERROR; `report` for the comparisons, `detail` for the others."
struct KKResult
    name::String
    kind::String
    outcome::String
    report::Union{Compare.Report,Nothing}
    detail::String
end

"""
`bsl` budget of a run: per record, `bsl(t) - bsl(t_prev)` against `-sum_d dV_bsl_d(t) / A_ocn(t)` recomputed
from the exports (`dV_bsl` of every domain's geo instance, `A_ocn` and `bsl` of the global geo), with `ratio` the
residual over the scale `|bsl(t)| + |bsl(t_prev)| + sum_d |dV_bsl_d(t)| / A_ocn(t)`; closes when every ratio is
within `tol` (a few `eps(Float32)`: the values and the reduction are in `sp`).
"""
struct Budget
    times::Vector{Float64}
    resid::Vector{Float64}
    ratio::Vector{Float64}
    max_ratio::Float64
    tol::Float64
    ok::Bool
    errors::Vector{String}
end

"""
    bsl_budget(rundir, geo_instances; geo_glob="geo_glob", tol=4*eps(Float32), dt=Dict()) -> Budget

`geo_instances`: the geo instance of each domain (`geo_nor`, ...; their `<instance>_ts.nc` hold `dV_bsl`).
`dt[instance]`: the cadence of that instance, its `dV_bsl` counts at the times that are multiples of it
(default: every record). The first record has no previous one and is not checked.
"""
function bsl_budget(rundir::AbstractString, geo_instances::AbstractVector{<:AbstractString}; geo_glob::AbstractString="geo_glob",
                    tol::Real=4 * eps(Float32), dt::AbstractDict=Dict{String,Float64}())
    errs = String[]
    gfile = joinpath(rundir, geo_glob * "_ts.nc")
    isfile(gfile) || return Budget(Float64[], Float64[], Float64[], NaN, tol, false, ["no file $gfile"])
    t, bsl, A = NCDataset(gfile) do ds
        for v in ("bsl", "A_ocn")
            haskey(ds, v) || push!(errs, "$gfile has no variable $v")
        end
        isempty(errs) ? (Compare.time_values(ds), Float64.(vec(Compare.read_raw(ds, "bsl")[1])),
                         Float64.(vec(Compare.read_raw(ds, "A_ocn")[1]))) : (Float64[], Float64[], Float64[])
    end
    isempty(errs) || return Budget(Float64[], Float64[], Float64[], NaN, tol, false, errs)
    n = length(t)
    dV = zeros(n)
    adV = zeros(n)
    for inst in geo_instances
        f = joinpath(rundir, inst * "_ts.nc")
        if !isfile(f)
            push!(errs, "no file $f (the geo instance needs out_1D, dV_bsl in its export list)")
            continue
        end
        NCDataset(f) do ds
            haskey(ds, "dV_bsl") || (push!(errs, "$f has no variable dV_bsl"); return)
            tt = Compare.time_values(ds)
            if length(tt) != n || !all(Compare.bitsequal.(tt, t))
                push!(errs, "$f: the records are not those of $gfile")
                return
            end
            v = Float64.(vec(Compare.read_raw(ds, "dV_bsl")[1]))
            d = Float64(get(dt, inst, 0.0))
            for k in 1:n
                (d > 0 && !isinteger(round(t[k] / d; digits=9))) && continue
                dV[k] += v[k]
                adV[k] += abs(v[k])
            end
        end
    end
    isempty(errs) || return Budget(Float64[], Float64[], Float64[], NaN, tol, false, errs)
    tt = t[2:end]
    resid = [(bsl[k] - bsl[k-1]) + dV[k] / A[k] for k in 2:n]
    scale = [abs(bsl[k]) + abs(bsl[k-1]) + adV[k] / A[k] for k in 2:n]
    ratio = [scale[i] == 0 ? abs(resid[i]) : abs(resid[i]) / scale[i] for i in eachindex(resid)]
    mr = isempty(ratio) ? 0.0 : maximum(ratio)
    return Budget(tt, resid, ratio, mr, tol, mr <= tol, errs)
end

function kk_instances(g::AbstractDict, domain::AbstractString)
    d = get(get(g, "domains", Dict()), domain, nothing)
    d === nothing && error("gate $(g["gate"]): no [domains.$domain]")
    return Dict{String,String}(String(k) => String(v) for (k, v) in get(d, "instances", Dict()))
end

"""
    run_kk_checks(g, runs; only=nothing) -> Vector{KKResult}

The `[[check]]`s of the gate file `g` over the run folders in `runs` (a directory with one folder per run name; the
chunk run has c00, c01, ...). `only`: the names of the checks to run.
"""
function run_kk_checks(g::AbstractDict, runs::AbstractString; only=nothing)
    out = KKResult[]
    for (i, c) in enumerate(kk_checks(g))
        ty = String(get(c, "type", ""))
        name = String(get(c, "name", ""))
        if isempty(name)
            name = ty == "domain" ? "$(c["domain"]): $(c["run"]) vs $(c["ref"])" :
                   ty in ("identical", "chunks") ? "$(c["run"]) vs $(c["ref"]) ($ty)" : "$ty $(get(c, "run", ""))"
        end
        only === nothing || name in only || continue
        try
            if ty in ("domain", "identical")
                refdir = joinpath(runs, c["ref"]); rundir = joinpath(runs, c["run"])
                inst = kk_instances(g, ty == "domain" ? c["domain"] : first(sort(collect(keys(g["domains"])))))
                rep = Compare.compare(refdir, rundir; map="identity", gate=g, instances=inst,
                                      restarts=get(c, "restarts", true))
                push!(out, KKResult(name, ty, Compare.outcome(rep), rep, ""))
            elseif ty == "chunks"
                inst = kk_instances(g, first(sort(collect(keys(g["domains"])))))
                rep = Compare.compare_chunks(joinpath(runs, c["ref"]), joinpath(runs, c["run"]); gate=g, instances=inst,
                                             restarts=get(c, "restarts", true))
                push!(out, KKResult(name, ty, Compare.outcome(rep), rep, ""))
            elseif ty == "budget"
                geos = String[]
                for d in get(c, "domains", sort(collect(keys(g["domains"]))))
                    push!(geos, kk_instances(g, d)["geo"])
                end
                dts = Dict{String,Float64}(String(k) => Float64(v) for (k, v) in get(c, "dt", Dict()))
                b = bsl_budget(joinpath(runs, c["run"]), geos; geo_glob=String(get(c, "geo_glob", "geo_glob")),
                               tol=Float64(get(c, "tol", 4 * eps(Float32))), dt=dts)
                detail = isempty(b.errors) ? @sprintf("%d records, max |residual|/scale = %.3g (bound %.3g), max |residual| = %.3g m",
                                                      length(b.times), b.max_ratio, b.tol, isempty(b.resid) ? 0.0 : maximum(abs, b.resid)) :
                         join(b.errors, "; ")
                push!(out, KKResult(name, ty, !isempty(b.errors) ? "ERROR" : b.ok ? "PASS" : "FAIL", nothing, detail))
            else
                push!(out, KKResult(name, ty, "ERROR", nothing, "unknown check type \"$ty\""))
            end
        catch e
            push!(out, KKResult(name, ty, "ERROR", nothing, sprint(showerror, e)))
        end
    end
    return out
end

"PASS when every check passes, else the worst of FAIL and ERROR."
kk_outcome(rs::AbstractVector{KKResult}) =
    isempty(rs) ? "ERROR" : all(r -> r.outcome == "PASS", rs) ? "PASS" : any(r -> r.outcome == "ERROR", rs) ? "ERROR" : "FAIL"

"The checks as a markdown table, then the per-check comparison summaries (failing variables)."
function kk_markdown(io::IO, g::AbstractDict, rs::AbstractVector{KKResult}; runs::AbstractString="")
    println(io, "# ", g["gate"], " (kryos against kryos): ", kk_outcome(rs), "\n")
    isempty(runs) || println(io, "- runs: `", runs, "`")
    println(io, "- date: ", Dates.format(now(), "yyyy-mm-dd HH:MM"), "\n")
    println(io, "| check | outcome | variables | failing | detail |\n|---|---|---|---|---|")
    for r in rs
        if r.report !== nothing
            rep = r.report
            cmp = [x for x in rep.results if !(x.status in ("absent", "new"))]
            nskip = count(x -> x.status == "skipped", cmp)
            println(io, "| ", r.name, " | ", r.outcome, " | ", length(cmp) - nskip, " (", nskip, " skipped) | ",
                    count(x -> x.status == "fail", cmp), " | ",
                    isempty(rep.new_files) ? "" : "files only in the run: " * join(first(rep.new_files, 4), ", ") *
                    (length(rep.new_files) > 4 ? ", ... ($(length(rep.new_files)))" : ""), " |")
        else
            println(io, "| ", r.name, " | ", r.outcome, " | | | ", r.detail, " |")
        end
    end
    for r in rs
        r.report === nothing && continue
        r.outcome == "PASS" && continue
        println(io, "\n## ", r.name, "\n")
        Compare.report_markdown(io, r.report)
    end
    return nothing
end

end # module
