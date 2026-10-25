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
- the `docs/gate.md` row (§4.4): `gate_row`, `append_gate_row`, `GATE_ROW_HEADER`.

The gate files are `tools/gate/<gate>.toml` (read by `Compare.load_gate`) and
`tools/gate/stack.toml`.
"""
module Gate

using TOML
using SHA
using Printf
using Dates
import FesmUtils.Namelist as NL
using ..Compare

export Stack, load_stack, stack_problems, stack_text, stack_id, stack_commits, DEFAULT_STACK
export plist, plist_string, chunk_overrides, restart_dir
export gate_root, stack_dir, ref_dir, run_dir, maps_dir, report_paths
export input_overlay, check_instances
export gate_row, append_gate_row, GATE_ROW_HEADER

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
    input_overlay(dest; kryos, yelmox, extra=[]) -> Dict{String,String}

Make the folder `dest` (the `input` link target of a kryos run directory) of
symbolic links to every entry of `yelmox/input` and then of `kryos/input`
(kryos wins a name that both have), and to the `extra` files (by default
`kryos/fesm-utils/par/datasrc_defaults.nml`, until the kryos `input/` has its
synchronised copy). kryos reads its templates (`kryos_defaults.nml`,
`forcing_defaults.nml`, `snapshots_defaults.nml`, `datasrc_defaults.nml`,
`kryos-variables-*.md`) and Yelmo and chion their defaults (`yelmo_defaults.nml`,
`chion_defaults.nml`), the phys_const file and the series from `input/`.
Returns name => link target. An existing `dest` is emptied of links first.
"""
function input_overlay(dest::AbstractString; kryos::AbstractString, yelmox::AbstractString, extra=nothing)
    ex = extra === nothing ? [joinpath(kryos, "fesm-utils", "par", "datasrc_defaults.nml")] : String.(extra)
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
    for f in ex
        haskey(out, basename(f)) && continue
        isfile(f) && link(f)
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

end # module
