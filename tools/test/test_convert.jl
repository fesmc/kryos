using KryosTools.Convert
import FesmUtils.Namelist as NL

const KRYOS = normpath(joinpath(@__DIR__, "..", ".."))

mini(; kw...) = mini_yelmox(mktempdir(); kw...)

function convert_mini(root; overrides=String[], kw...)
    convert_yelmox(; yelmox=root, par="yelmox/yelmox_Test.nml", overrides=overrides, kryos=KRYOS, kw...)
end

refused(f) = try
    f(); ""
catch e
    e isa ConversionError || rethrow()
    e.msg
end

kv(groups, g) = Dict(k.name => NL.parse_value(k.raw) for k in NL.group_dict(groups)[g].keys)
namelist(conv) = NL.Group[NL.Group(og.name, [NL.Key(k.name, k.raw) for k in og.keys]) for og in conv.groups]
kvg(G, g) = Dict(k.name => NL.parse_value(k.raw) for k in G[g].keys)

@testset "Convert: miniature yelmox" begin
    root = mini()
    conv = convert_mini(root; overrides=["ctrl.time_end=1000", "tm_2D.dt=100", "ctrl.write_lossless=True"])
    out = joinpath(mktempdir(), "test.nml")
    write_par(out, conv)
    G = NL.read_namelist(out)
    names = [g.name for g in G]
    @test length(unique(names)) == length(names)

    @testset "core groups" begin
        @test kv(G, "kryos")["domains"] == "Testland"
        @test kv(G, "kryos")["restart"] == "None"
        @test kv(G, "kryos")["write_lossless"] === true
        @test kv(G, "kryos")["library"] == ""
        @test kv(G, "kryos")["phys_const_file"] == "input/phys_const.nml"
        c = kv(G, "clock")
        @test c["time_end"] === 1000.0 && c["dtt"] === 10.0 && c["method"] == "const"
        @test haskey(c, "time_const") && !haskey(c, "time_const_cal")
        @test kv(G, "out_2D")["dt"] === 100.0
        @test kv(G, "out_1D")["dt"] === 10.0               # tm_1D dt = 1 raised to dtt
        @test any(occursin("tm_1D.dt = 1.0 is below dtt", n) for n in conv.notes)
        d = kv(G, "Testland")
        @test d["ice"] == "ice_tes" && d["sediments"] == "sed_tes" && d["grid_hub"] == "TST-16KM"
        @test d["out_2D"] == "out_2D" && d["basins"] == "Zwally2012"
        @test d["experiments"] == ""                       # sim.relax = sim.opt = False: no experiment
        @test kv(G, "globe") == Dict("forcing" => "forcing", "geo" => "geo_glob")
    end

    @testset "instances" begin
        ice = kv(G, "ice_tes")
        @test ice["backend"] == "yelmo" && ice["grid"] == "TST-32KM" && ice["init_equil_time"] === 10.0
        @test ice["init_time_thrm"] === 0.0
        @test ice["out_yelmo"] == "out_2D" && ice["out_yelmo_sm"] == "out_2Dsm" && ice["out_yelmo_ts"] == "out_1D"
        geo = kv(G, "geo_tes")
        @test geo["dt"] === 20.0 && geo["isos"] == "isos" && !haskey(geo, "bsl")
        @test geo["out_2D"] == "none" && geo["out_1D"] == "none"        # write_isos = False
        gg = kv(G, "geo_glob")
        @test gg["bsl"] == "barysealevel" && !haskey(gg, "isos") && gg["clock"] == "forcing"
        cl = kv(G, "climate_tes")
        @test cl["backend"] == "snapshots" && cl["dt"] === 10.0 && cl["terms"] == "dT" && cl["reference"] == "pd"
        @test kv(G, "forcing")["coefs"] == ["alpha", "dT"]
        @test kv(G, "surface_tes")["backend"] == "chion" && kv(G, "surface_tes")["orbit"] == "orbit"
        @test kv(G, "shelf_tes")["marine_shelf"] == "marine_shelf"
    end

    @testset "data instances" begin
        g = kv(G, "ghf_tes")
        @test g["backend"] == "data" && g["dt"] === -1.0 && g["Q_geo"] == "Q_geo_tes" && g["Q_geo_min"] === 1e-4
        f = kv(G, "Q_geo_tes")
        @test f["mode"] == "static" && f["name"] == "ghf" && f["fill"] == "nearest" && f["remap"] == "auto"
        @test occursin("GHF.nc", f["filename"])
        s = kv(G, "H_sed_tes")
        @test s["mode"] == "const" && s["value"] === 0.0
    end

    @testset "library groups copied, physics groups expanded" begin
        for g in ("coef_alpha", "coef_dT", "term_dT", "snapshot_pd", "tas_pd", "chion_const", "itm_chion", "ytrc", "ytopo")
            @test g in names
        end
        @test !("tas_lgm" in names)                     # not named by any key
        @test !("relax" in names) && !("opt" in names) && !("relax_tes" in names) && !("smbpal" in names)
        @test kv(G, "snapshots")["f_p"] === 0.06 && kv(G, "snapshots")["grid_global"] == "GLOBAL-1DEG"
        @test kv(G, "ytopo") == Dict("a" => 2.0, "b" => "x")
        @test kv(G, "yelmo")["restart"] == "None" && !haskey(kv(G, "yelmo"), "phys_const")
        @test kv(G, "marine_shelf")["restart"] == "None"          # was "none"
        # the basins and thermal-forcing corrections belong to the shelf component's group
        mk = keys(kv(G, "marine_shelf"))
        @test !any(k -> k in ("basins", "corr_method", "corr_names", "corr_values", "corr_b1", "tf_correction", "tf_path", "tf_name"), mk)
        sc = kv(G, "shelf_corr_tes")
        @test sc["basins"] == "Zwally2012" && sc["corr_method"] == "tf" && sc["corr_names"] == ["b1", "b2"]
        @test sc["corr_values"] == [0.3, -1.5] && sc["corr_b2"] == "Zwally2012:2" && sc["tf_correction"] === false
        @test sc["tf_path"] == "ice_data/{domain}/{grid_name}/OCEAN.nc" && sc["tf_name"] == "dT_nl"
        @test kv(G, "shelf_tes")["corr"] == "shelf_corr_tes"
        @test findfirst(==("shelf_corr_tes"), names) > findfirst(==("marine_shelf"), names)
        @test kv(G, "isos")["restart"] == "None"
        @test kv(G, "chion")["model"] == "itm"
        @test kv(G, "itm_chion")["alb_ice"] === 0.31
    end

    @testset "accounting" begin
        s = account_summary(conv)
        @test s.total == s.mapped + s.default + s.dropped
        @test s.dropped > 0 && s.mapped > 0 && s.default > 0
        keyof(a) = "$(a.group).$(a.key)"
        d = Dict(keyof(a) => a for a in conv.acc)
        @test d["ctrl.time_end"].status == :mapped && d["ctrl.time_end"].dst == "clock.time_end"
        @test d["ctrl.run_step"].status == :dropped && !isempty(d["ctrl.run_step"].reason)
        @test d["comps.with_shelf"].status == :dropped
        @test d["sim.init_method"].status == :mapped
        @test d["ghf.ghf_const"].status == :dropped
        @test d["ghf.obs_path"].dst == "Q_geo_tes.filename"
        @test d["relax.tau1"].status == :dropped && occursin("sim.relax = False", d["relax.tau1"].reason)
        @test d["opt.opt_cf"].status == :dropped && d["opt.tf_max"].status == :dropped
        @test d["yelmo.phys_const"].status == :dropped
        @test d["marine_shelf.corr_values"].dst == "shelf_corr_tes.corr_values" && d["marine_shelf.basins"].dst == "shelf_corr_tes.basins"
        @test d["tm_1D.dt"].dst == "out_1D.dt"
        @test d["output.write_isos"].status == :default            # off: nothing to set
        # every key of every resolved group is accounted exactly once
        @test length(conv.acc) == length(unique(keyof.(conv.acc)))
        # every dropped key has a reason
        @test all(!isempty(a.reason) for a in conv.acc if a.status == :dropped)
    end

    @testset "provenance and checks" begin
        files = Dict(conv.files)
        @test haskey(files, "input/yelmo_phys_const.nml") && length(files["input/yelmo_phys_const.nml"]) == 64
        @test haskey(files, "input/alpha.dat") && haskey(files, "input/timeout.txt")
        h = Convert.par_header(conv)
        @test occursin("yelmo_phys_const.nml", h) && occursin("-p: ctrl.time_end=1000", h)
        @test any(startswith(c, "round trip:") for c in conv.checks)
        # a report is written
        rep = joinpath(dirname(out), "report.md")
        write_report(rep, conv)
        txt = read(rep, String)
        @test occursin("## Dropped", txt) && occursin("`ctrl.run_step`", txt) && occursin("## Files a run directory needs", txt)
    end
end

@testset "Convert: -p semantics" begin
    root = mini()
    conv = convert_mini(root; overrides=["ctrl.time_end=500", "forcing.coefs==alpha", "snapshots.f_p=0.07"])
    G = Dict(g.name => g for g in namelist(conv))
    @test kvg(G, "clock")["time_end"] === 500.0
    @test kvg(G, "forcing")["coefs"] == "alpha"
    @test kvg(G, "snapshots")["f_p"] === 0.07
    # init_time_thrm is the yelmo key of the ice instance
    conv4 = convert_mini(root; overrides=["sim.init_time_thrm=5.0"])
    G4 = Dict(g.name => g for g in namelist(conv4))
    @test kvg(G4, "ice_tes")["init_time_thrm"] === 5.0
    # an alias of info.json: ghfx -> ghf
    conv2 = convert_mini(root; overrides=["ghfx.use_obs=False"])
    G2 = Dict(g.name => g for g in namelist(conv2))
    @test kvg(G2, "Q_geo_tes")["mode"] == "const" && kvg(G2, "Q_geo_tes")["value"] === 0.05
    @test kvg(G2, "ghf_tes")["Q_geo_min"] === 0.0
    # a group the par file does not have is copied from the library when a -p names a key of it
    conv3 = convert_mini(root; overrides=["term_dT.fields==tas"])
    G3 = Dict(g.name => g for g in namelist(conv3))
    @test kvg(G3, "term_dT")["fields"] == "tas"
    # a key neither in the group nor its defaults is an error of the resolution
    @test occursin("not in", refused(() -> convert_mini(root; overrides=["ctrl.nokey=1"])))
end


@testset "Convert: experiments" begin
    ovr = ["sim.relax=True", "sim.opt=True"]
    conv = convert_mini(mini(); overrides=ovr)
    out = joinpath(mktempdir(), "x.nml")
    write_par(out, conv)
    G = NL.read_namelist(out)
    @test kv(G, "Testland")["experiments"] == ["relax_tes", "opt_cf_tes"]
    r = kv(G, "relax_tes")
    @test r["kind"] == "relax" && r["topo_rel"] === 4 && r["tau1"] === 100.0 && r["time1"] === 2000.0 && r["m"] === 2.0
    o = kv(G, "opt_cf_tes")
    @test o["kind"] == "opt_cb_ref" && o["method"] == "L21" && o["cf_init"] === 0.01 && o["sigma_err"] === 50e3
    @test o["fill_method"] == "cf_min" && o["cf_time_end"] === 15e3
    @test !haskey(o, "tf_min") && !haskey(o, "opt_cf")
    @test !("relax" in [g.name for g in G]) && !("opt" in [g.name for g in G])
    keyof(a) = "$(a.group).$(a.key)"
    d = Dict(keyof(a) => a for a in conv.acc)
    @test d["relax.tau2"].dst == "relax_tes.tau2" && d["opt.opt_cf"].dst == "opt_cf_tes.method"
    @test d["opt.opt_tf"].status == :dropped && occursin("no experiment", d["opt.opt_tf"].reason)   # opt_tf: WP3
    @test d["opt.tf_max"].status == :dropped
    @test d["sim.relax"].status == :dropped && d["sim.opt"].status == :dropped
    @test length(conv.acc) == length(unique(keyof.(conv.acc)))
    @test any(occursin("experiments of the domain", n) for n in conv.notes)
    # only relax
    G1 = Dict(g.name => g for g in namelist(convert_mini(mini(); overrides=["sim.relax=True"])))
    @test kvg(G1, "Testland")["experiments"] == "relax_tes" && !haskey(G1, "opt_cf_tes")
    # opt_cf = none keeps the experiment (cb_ref from outside, no optimisation)
    G2 = Dict(g.name => g for g in namelist(convert_mini(mini(); overrides=["sim.opt=True", "opt.opt_cf=none"])))
    @test kvg(G2, "opt_cf_tes")["method"] == "none" && kvg(G2, "Testland")["experiments"] == "opt_cf_tes"
    # a base that selects an experiment its par file does not have
    msg = refused(() -> convert_mini(mini(edits=["relax             = False" => "relax             = True", "&relax" => "&relax_gone"])))
    @test occursin("no [relax] group", msg)
end

@testset "Convert: refusals name the key and the phase" begin
    cases = [
        (["sim.lim_pd_ice=True"] => "sim.lim_pd_ice = true"),
        (["sim.opt=True", "opt.opt_tf=L21"] => "opt.opt_tf"),
        (["sim.init_method=recon"] => "sim.init_method"),
        (["comps.climate=esm"] => "comps.climate"),
        (["comps.with_isostasy=False"] => "comps.with_isostasy"),
        (["comps.surface_method=smbpal"] => "comps.surface_method"),
        (["ctrl.restart=restart-1.000-kyr"] => "ctrl.restart"),
        (["marine_shelf.bmb_method=pico"] => "marine_shelf.bmb_method"),
        (["marine_shelf.bmb_method=quad-slope"] => "marine_shelf.bmb_method"),
        (["output.write_cmip=True"] => "output.write_cmip"),
        (["ghf.f_stdev=1.0"] => "ghf.f_stdev"),
        (["sim.use_negis=True"] => "sim.use_negis"),
    ]
    for (ovr, key) in cases
        root = mini()
        msg = refused(() -> convert_mini(root; overrides=ovr))
        @test !isempty(msg)
        @test occursin(lowercase(key), lowercase(msg)) && occursin("phase", msg)
    end
    # a refusal lists every reason, not only the first
    msg = refused(() -> convert_mini(mini(); overrides=["sim.lim_pd_ice=True", "sim.use_negis=True"]))
    @test occursin("sim.lim_pd_ice", msg) && occursin("sim.use_negis", msg)
    # opt_tf is refused only when sim.opt selects it
    @test convert_mini(mini(); overrides=["opt.opt_tf=L21"]) isa Conversion
    # a group of a rembo par file, and a suffixed group
    msg = refused(() -> convert_mini(mini(edits=["&relax" => "&rembo\n    a = 1.0\n/\n&relax"])))
    @test occursin("&rembo", msg) && occursin("phase 2", msg)
    msg = refused(() -> convert_mini(mini(edits=["&relax" => "&ctrl_north\n    a = 1.0\n/\n&relax"])))
    @test occursin("multi-domain", msg) && occursin("phase 3", msg)
end

@testset "Convert: unaccounted keys and groups stop the conversion" begin
    msg = refused(() -> convert_mini(mini(edits=["write_lossless  = False" => "write_lossless  = False\n    newkey = 1"])))
    @test occursin("ctrl.newkey", msg) && occursin("not accounted", msg)
    # a group without a rule is not read by yelmox: a base group nothing names is dropped with a reason
    root = mini(edits=["&relax" => "&coef_unused\n    method = \"const\"\n/\n&relax"])
    conv = convert_mini(root)
    a = [x for x in conv.acc if x.group == "coef_unused"]
    @test !isempty(a) && all(x.status == :dropped && occursin("not read by yelmox", x.reason) for x in a)
end

@testset "Convert: validation against the kryos templates" begin
    conv = convert_mini(mini())
    @test validate_conversion(conv) isa Vector{String}
    # an undeclared key in a core group
    og = first(g for g in conv.groups if g.name == "ice_tes")
    push!(og.keys, Convert.OutKey("no_such_key", "1", 1, "test", "", nothing, true, nothing))
    @test occursin("no_such_key", refused(() -> validate_conversion(conv)))
    pop!(og.keys)
    # a schedule dt that is not a multiple of dtt
    sg = first(g for g in conv.groups if g.name == "out_2D")
    k = first(x for x in sg.keys if x.name == "dt")
    old = k.raw; k.raw = "15.0"; k.value = 15.0
    @test occursin("multiple of dtt", refused(() -> validate_conversion(conv)))
    k.raw = old; k.value = parse(Float64, old)
    # a named group that is not in the file
    ig = first(g for g in conv.groups if g.name == "ice_tes")
    y = first(x for x in ig.keys if x.name == "yelmo")
    y.raw = "\"nowhere\""; y.value = "nowhere"
    @test occursin("nowhere", refused(() -> validate_conversion(conv)))
end

@testset "Convert: written file reads back" begin
    conv = convert_mini(mini())
    out = joinpath(mktempdir(), "t.nml")
    write_par(out, conv)
    @test Convert.check_roundtrip(conv, out)
    # a doctored file is caught
    txt = read(out, String)
    write(out, replace(txt, "15e3" => "99.0"; count=1))
    @test occursin("round trip", refused(() -> Convert.check_roundtrip(conv, out)))
end

@testset "Convert: --keep-library keeps the references" begin
    conv = convert_mini(mini(); keep_library=true)
    names = [g.name for g in conv.groups]
    @test !("coef_dT" in names) && !("tas_pd" in names) && "coef_alpha" in names
    @test conv.groups[1].name == "kryos"
    lib = first(k for k in conv.groups[1].keys if k.name == "library")
    @test lib.value == ["input/clim_db.nml", "input/forcing_db.nml"]
end

@testset "Convert: output switches" begin
    # the geo stream keys are not declared by the kryos templates yet: left out with a warning
    conv = convert_mini(mini(edits=["write_isos    = False" => "write_isos    = True"]))
    geo = first(g for g in conv.groups if g.name == "geo_tes")
    @test first(k for k in geo.keys if k.name == "out_2D").value == "out_2D"
    declared = any(k.name == "out_isos" for k in geo.keys)
    @test declared || any(occursin("out_isos", w) for w in conv.warnings)
end

# the real yelmox bases, when a yelmox checkout is at hand (make check-tools sets YELMOXROOT)
const YELMOXROOT = get(ENV, "YELMOXROOT", expanduser("~/models/yelmox/.claude/worktrees/topo-prep"))

if isfile(joinpath(YELMOXROOT, ".runme", "info.json"))
    @testset "Convert: the yelmox bases ($YELMOXROOT)" begin
        base = joinpath(YELMOXROOT, "yelmox")
        gate_p = ["ctrl.time_end=1000", "ctrl.write_lossless=True", "tm_2D.dt=100", "tm_2Dsm.dt=10", "tm_1D.dt=10", "tm_rst.dt=500"]
        for (par, extra) in (("yelmox_Greenland.nml", String[]), ("yelmox_Antarctica.nml", ["sim.relax=False", "sim.opt=False"]),
                             ("yelmox_Antarctica.nml", ["opt.opt_tf=none"]))
            isfile(joinpath(base, par)) || continue
            conv = convert_yelmox(; yelmox=YELMOXROOT, par="yelmox/$par", overrides=vcat(gate_p, extra), kryos=KRYOS)
            s = account_summary(conv)
            @test s.total == s.mapped + s.default + s.dropped
            out = joinpath(mktempdir(), "x.nml")
            write_par(out, conv)
            G = NL.read_namelist(out)
            @test kv(G, "clock")["time_end"] === 1000.0
            @test kv(G, "kryos")["write_lossless"] === true
            for g in ("yelmo", "isos", "barysealevel", "marine_shelf", "chion")
                @test kv(G, g)["restart"] == "None"
            end
        end
        # the others are refused, naming a phase
        others = ["yelmox_GreenlandPaleo.nml" => "bias.active", "yelmox_LIS.nml" => "comps.surface_method",
                  "yelmox_North.nml" => "sim.init_method", "yelmox_Pyrenees.nml" => "comps.with_shelf",
                  "yelmox_SRG.nml" => "comps.with_isostasy", "yelmox_esm_Antarctica.nml" => "comps.climate",
                  "yelmox_esm_Greenland.nml" => "comps.climate", "yelmox_rembo_Greenland.nml" => "comps.climate",
                  "rembo_Greenland.nml" => "&rembo"]
        for (par, what) in others
            isfile(joinpath(base, par)) || continue
            msg = refused(() -> convert_yelmox(; yelmox=YELMOXROOT, par="yelmox/$par", kryos=KRYOS))
            @test occursin(what, msg) && occursin("phase", msg)
        end
        # the shipped Antarctica base needs the plugins: refused without the gate's overrides
        if isfile(joinpath(base, "yelmox_Antarctica.nml"))
            msg = refused(() -> convert_yelmox(; yelmox=YELMOXROOT, par="yelmox/yelmox_Antarctica.nml", kryos=KRYOS))
            @test occursin("opt.opt_tf", msg) && occursin("phase", msg)
            # the interim gate ant-cf: relax and the cb_ref optimisation, without opt_tf
            conv = convert_yelmox(; yelmox=YELMOXROOT, par="yelmox/yelmox_Antarctica.nml", overrides=vcat(gate_p, ["opt.opt_tf=none"]), kryos=KRYOS)
            out = joinpath(mktempdir(), "ant.nml")
            write_par(out, conv)
            G = NL.read_namelist(out)
            @test kv(G, "Antarctica")["experiments"] == ["relax_ant", "opt_cf_ant"]
            @test kv(G, "relax_ant")["time2"] === 6000.0 && kv(G, "opt_cf_ant")["cf_init"] === 0.01
        end
    end
else
    @info "no yelmox checkout at $YELMOXROOT (set YELMOXROOT): the real bases are not converted"
end
