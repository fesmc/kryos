using KryosTools
using KryosTools.Gate
using KryosTools.Compare
import FesmUtils.Namelist as NL
using TOML

@testset "Gate: stack" begin
    s = load_stack()
    @test isempty(stack_problems(s))
    id = stack_id(s)
    @test length(id) == 8 && all(c -> c in "0123456789abcdef", id)
    @test stack_id(load_stack()) == id
    # the id ignores the kryos commit and the notes, not a library or a flag
    d = deepcopy(s.data)
    d["kryos"]["commit"] = "deadbeef"; d["note"] = "other"; d["name"] = "other"
    @test stack_id(Stack(d, "")) == id
    d["libs"]["yelmo"]["commit"] = "0000000"
    @test stack_id(Stack(d, "")) != id
    d2 = deepcopy(s.data); d2["build"]["flags"] = "-Ofast"
    @test stack_id(Stack(d2, "")) != id
    @test occursin("[libs.yelmo]", stack_text(s)) && !occursin("[kryos]", stack_text(s))
    # problems
    d3 = deepcopy(s.data); delete!(d3["libs"], "insol"); delete!(d3["run"], "omp_threads")
    bad = stack_problems(Stack(d3, ""))
    @test any(occursin("libs.insol", b) for b in bad) && any(occursin("omp_threads", b) for b in bad)
    @test_throws ErrorException load_stack(joinpath(mktempdir(), "none.toml"))
    @test any(occursin("full sha", b) for b in stack_problems(s; full=true))
    @test first(first(stack_commits(s))) == "yelmox"
    @test length(stack_commits(s)) == 10
    @test s.data["run"]["omp_threads"] == 16
end

@testset "Gate: -p lists and chunks" begin
    for (name, extra) in (("grl", String[]), ("ant", ["sim.relax=False", "sim.opt=False"]))
        g = load_gate(name)
        p = plist(g)
        @test "ctrl.time_end=1000" in p && "ctrl.write_lossless=True" in p && "tm_2Dsm.dt=10" in p
        @test all(e -> e in p, extra)
        @test (name == "grl") == !any(startswith("sim."), p)
        @test plist_string(g) == join(p, " ")
        @test g["instances"]["ice"] == "ice_$name"
    end
    g = load_gate("grl")
    @test chunk_overrides(g, :yelmox, 0) == ["ctrl.time_init=0", "ctrl.time_end=500"]
    @test chunk_overrides(g, :yelmox, 1; prev_dir="/r/c00") ==
          ["ctrl.time_init=500", "ctrl.time_end=1000", "ctrl.restart=/r/c00/restart-0.500-kyr"]
    @test chunk_overrides(g, :kryos, 1; prev_dir="/r/c00") ==
          ["clock.time_init=500", "clock.time_end=1000", "kryos.restart=/r/c00/restart-500"]
    @test_throws ErrorException chunk_overrides(g, :kryos, 1)
    @test_throws ErrorException chunk_overrides(g, :kryos, 2; prev_dir="/x")
    @test restart_dir(:yelmox, 1000) == "restart-1.000-kyr"
    @test restart_dir(:kryos, 1000.0) == "restart-1000" && restart_dir(:kryos, 0.5) == "restart-0.5"
    @test_throws ErrorException restart_dir(:other, 1)
end

@testset "Gate: directories" begin
    @test gate_root() == get(ENV, "GATE_ROOT", "/work/ba1442/robinson/kryos-gate")
    r = "/g"
    @test ref_dir("grl", "ab12cd34", "74f3a2bac47b", "cont"; root=r) == "/g/ref/grl/yx-74f3a2b/ab12cd34/cont"
    @test run_dir("ant", "ab12cd34", "5e46158abc", "chunk"; root=r) == "/g/runs/ant/kx-5e46158/ab12cd34/chunk"
    @test maps_dir("S2"; root=r) == "/g/maps/S2" && stack_dir("S2"; root=r) == "/g/stacks/S2"
    @test report_paths("5e46158abc", "ab12cd34"; date="2026-10-25", root=r).md == "/g/reports/2026-10-25_kx-5e46158_ab12cd34.md"
    @test_throws ErrorException ref_dir("grl", "s", "74f3a2b", "chunks"; root=r)
    @test_throws ErrorException run_dir("grl", "s", "74f3a2b", "cont2"; root=r)
    @test_throws ErrorException ref_dir("grl", "s", "74f", "cont"; root=r)
end

@testset "Gate: input overlay" begin
    d = mktempdir()
    yx = joinpath(d, "yx"); kr = joinpath(d, "kr"); fu = joinpath(kr, "fesm-utils", "par")
    mkpath(joinpath(yx, "input", "esm")); mkpath(joinpath(kr, "input")); mkpath(fu)
    write(joinpath(yx, "input", "yelmo_defaults.nml"), "y")
    write(joinpath(yx, "input", "shared.nml"), "yelmox")
    write(joinpath(kr, "input", "kryos_defaults.nml"), "k")
    write(joinpath(kr, "input", "shared.nml"), "kryos")
    write(joinpath(fu, "datasrc_defaults.nml"), "d")
    dest = joinpath(d, "run", "input")
    links = input_overlay(dest; kryos=kr, yelmox=yx)
    @test Set(keys(links)) == Set(["yelmo_defaults.nml", "shared.nml", "esm", "kryos_defaults.nml", "datasrc_defaults.nml"])
    @test read(joinpath(dest, "shared.nml"), String) == "kryos"          # kryos wins
    @test read(joinpath(dest, "datasrc_defaults.nml"), String) == "d"
    @test islink(joinpath(dest, "esm")) && isdir(joinpath(dest, "esm"))
    # again: links are replaced, not duplicated
    @test input_overlay(dest; kryos=kr, yelmox=yx) == links
    # a real file in the way is an error
    write(joinpath(dest, "alien.nml"), "x")
    @test isfile(joinpath(dest, "alien.nml"))
    # the kryos input of this checkout has its templates
    @test isfile(joinpath(KRYOS_ROOT, "input", "kryos_defaults.nml"))
end

@testset "Gate: instances against a par file" begin
    d = mktempdir()
    par = joinpath(d, "p.nml")
    write(par, """
    &globe
        forcing = "forcing"
        geo     = "geo_glob"
    /
    &Greenland
        ice = "ice_grl"
        geo = "geo_grl"
        climate = "climate_grl"
        surface = "surface_grl"
        shelf = "shelf_grl"
        geothermal = "ghf_grl"
        sediments = "sed_grl"
    /
    """ * join(("&$n\n    backend = \"x\"\n/\n" for n in ("forcing", "geo_glob", "ice_grl", "geo_grl", "climate_grl", "surface_grl", "shelf_grl", "ghf_grl", "sed_grl")), ""))
    g = load_gate("grl")
    @test isempty(check_instances(g, par))
    @test !isempty(check_instances(load_gate("ant"), par))
    write(par, replace(read(par, String), "ice = \"ice_grl\"" => "ice = \"ice_x\""))
    bad = check_instances(g, par)
    @test any(occursin("ice_x", b) for b in bad)
    # the gate files against the committed gate par files
    for n in ("grl", "ant")
        f = joinpath(KRYOS_ROOT, "par", "gate", "$n.nml")
        isfile(f) && @test isempty(check_instances(load_gate(n), f))
    end
end

@testset "Gate: docs/gate.md row" begin
    row = gate_row(; date="2026-10-25", gate="grl", kryos="5e46158abc", yelmox="74f3a2bac", stack="ab12cd34" => "S2",
                   flags="-O2 -fp-model precise", length=1000.0, outcome="PASS", report="/g/reports/x.md")
    cols = split(strip(row, ['|', ' ']), " | ")
    @test length(cols) == 11
    @test cols[1] == "2026-10-25" && cols[3] == "5e46158" && cols[4] == "74f3a2b" && cols[5] == "ab12cd34 (S2)"
    @test cols[7] == "1000 yr" && cols[8] == "PASS" && cols[9] == "none" && cols[10] == "none"
    @test length(split(GATE_ROW_HEADER, '\n')) == 2 && occursin("---", GATE_ROW_HEADER)
    # from reports: a failing one names the first divergence
    rep = Compare.Report("a", "b", "identity", "identity")
    rep.first_divergence = Dict{String,Any}("time" => 20.0, "static" => false,
                                            "variables" => [Dict("key" => "a.nc:x", "comp" => "", "run" => "a.nc:x", "max_abs" => 1.0)])
    push!(rep.results, Compare.VarResult(key="a.nc:x", ref_file="a.nc", ref_vars=["x"], class="E", status="fail"))
    ok = Compare.Report("a", "b", "identity", "identity")
    r2 = gate_row([ok, rep]; date="2026-10-25", gate="ant", kryos="5e46158", yelmox="74f3a2b", stack="s", flags="f", length=1000, report="r")
    @test occursin("FAIL", r2) && occursin("t = 20: a.nc:x", r2)
    @test occursin("PASS", gate_row([ok]; date="d", gate="g", kryos="5e46158", yelmox="74f3a2b", stack="s", flags="f", length=1, report="r"))
    # appended to a new file
    f = joinpath(mktempdir(), "docs", "gate.md")
    append_gate_row(f, row); append_gate_row(f, r2)
    lines = readlines(f)
    @test count(startswith("| 2026"), lines) == 2 && any(l -> l == split(GATE_ROW_HEADER, '\n')[1], lines)
end
