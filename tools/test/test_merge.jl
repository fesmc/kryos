using KryosTools.Convert
import FesmUtils.Namelist as NL

const MKRYOS = normpath(joinpath(@__DIR__, "..", ".."))

# a small converted-style single-domain file
function mini_domain(dir, name, sfx; coefs="", tc="0.0", extra=String[], dom_extra="")
    path = joinpath(dir, "$name.nml")
    open(path, "w") do io
        println(io, "! header of $name\n! source: yelmox/yelmox_$name.nml\n! -p: (none)")
        println(io, """
&kryos
    domains = "$name"
    clock   = "clock"
/

&clock
    time_end   = 100.0
    time_const = $tc
/

&globe
    forcing = "forcing"
    geo     = "geo_glob"
/

&forcing
    backend = "tsgen"
    coefs   = $coefs
/

&geo_glob
    backend = "fastisostasy"
    bsl     = "barysealevel"
/

&barysealevel
    method = "fastiso"
/

&out_2D
    dt = 10.0
/

&$name
    grid_hub   = "HUB"
    topo       = "topo"
    regions    = "regions"
    ice        = "ice_$sfx"
    geo        = "geo_$sfx"
    climate    = "climate_$sfx"
$dom_extra/

&ice_$sfx
    backend = "yelmo"
    yelmo   = "yelmo"
/

&geo_$sfx
    backend = "fastisostasy"
    isos    = "isos"
/

&climate_$sfx
    backend   = "snapshots"
    snapshots = "snapshots"
    reference = "pd"
    terms     = "lgm" "ocean"
/

&topo
    path = "topo_$sfx.nc"
/

&regions
    layers = "basin"
/

&yelmo
    nml_ytopo = "ytopo"
    nz_aa     = $(sfx == "aaa" ? 5 : 10)
/

&ytopo
    solver = "$(sfx == "aaa" ? "expl" : "impl-lis")"
/

&isos
    nl = 1
/

&snapshots
    f_p = 0.05
/

&snapshot_pd
    tas = "tas_$sfx"
/

&snapshot_lgm
    tas = "tas_lgm"
/

&term_lgm
    coef  = "alpha"
    snap0 = "pd"
    snap1 = "lgm"
    pattern = "snapshots"
/

&term_ocean
    coef = "dTo"
/

&tas_$sfx
    name = "t2m_$sfx"
/
""")
        for l in extra
            println(io, l)
        end
    end
    return path
end

widegroup(n) = n in ("kryos", "clock", "forcing", "barysealevel", "globe", "geo_glob") || startswith(n, "out_") || startswith(n, "coef_")
mkey(groups, g, k) = NL.parse_value(only(x for x in NL.group_dict(groups)[g].keys if x.name == k).raw)
gnames(groups) = [g.name for g in groups]
failure(f) = try
    f(); ""
catch e
    e isa ConversionError || rethrow()
    e.msg
end

@testset "Convert: merge of two single-domain files" begin
    @testset "synthetic: renaming, repointing, run-wide groups" begin
        d = mktempdir()
        a = mini_domain(d, "Aland", "aaa"; coefs="\"alpha\" \"dTo\"", tc="-21e3")
        b = mini_domain(d, "Bland", "bbb"; coefs="\"\"", extra=["&term_other", "    coef = \"\"", "/"])
        res = Dict("clock.time_const" => "0.0")
        m = merge_par([a, b]; resolve=res)
        G = m.groups
        @test m.domains == ["Aland", "Bland"] && m.suffixes == ["aaa", "bbb"]
        @test length(unique(gnames(G))) == length(G)
        @test mkey(G, "kryos", "domains") == ["Aland", "Bland"]
        @test mkey(G, "forcing", "coefs") == ["alpha", "dTo"]                  # the union
        @test mkey(G, "clock", "time_const") === 0.0
        @test mkey(G, "clock", "time_end") === 100.0
        # shared domain groups renamed in both files, repointed
        for n in ("topo", "regions", "yelmo", "ytopo", "isos", "snapshots", "snapshot_pd", "snapshot_lgm", "term_lgm", "term_ocean")
            @test n ∉ gnames(G) && "$(n)_aaa" in gnames(G) && "$(n)_bbb" in gnames(G)
        end
        @test mkey(G, "Aland", "topo") == "topo_aaa" && mkey(G, "Bland", "topo") == "topo_bbb"
        @test mkey(G, "ice_aaa", "yelmo") == "yelmo_aaa" && mkey(G, "yelmo_bbb", "nml_ytopo") == "ytopo_bbb"
        @test mkey(G, "climate_aaa", "reference") == "pd_aaa"                  # snapshot_pd_aaa
        @test mkey(G, "climate_bbb", "terms") == ["lgm_bbb", "ocean_bbb"]
        @test mkey(G, "term_lgm_aaa", "snap0") == "pd_aaa" && mkey(G, "term_lgm_bbb", "snap1") == "lgm_bbb"
        @test mkey(G, "term_lgm_aaa", "pattern") == "snapshots"                # a keyword, not a group
        @test mkey(G, "ytopo_aaa", "solver") == "expl" && mkey(G, "ytopo_bbb", "solver") == "impl-lis"
        # a group only one file has keeps its name
        @test "term_other" in gnames(G) && "tas_aaa" in gnames(G) && "tas_bbb" in gnames(G)
        @test "tas_aaa" ∉ [r[2] for r in m.renamed]
        # the instance groups are the single-domain names; the run-wide groups are not renamed
        @test all(n in gnames(G) for n in ("ice_aaa", "geo_bbb", "climate_bbb", "barysealevel", "geo_glob", "globe", "out_2D"))
        @test length(m.renamed) == 2 * 10
        # every key of a domain that names a group reaches a group of that domain
        for i in 1:2
            S = NL.group_dict(split_merged(m, i))
            @test haskey(S, m.domains[i]) && haskey(S, "yelmo") && haskey(S, "snapshot_pd")
        end
        # the inverse gives back the original file's groups
        for (i, f) in enumerate((a, b))
            orig = NL.group_dict(NL.read_namelist(f))
            back = NL.group_dict(split_merged(m, i))
            for g in values(back)
                widegroup(g.name) && continue
                @test haskey(orig, g.name)
                haskey(orig, g.name) && @test [(k.name, NL.parse_value(k.raw)) for k in g.keys] ==
                                              [(k.name, NL.parse_value(k.raw)) for k in orig[g.name].keys]
            end
            reached = Set(keys(back))
            @test all(n in reached for n in keys(orig) if !widegroup(n) && n != "term_other")
        end
        # file order is the domain order
        mr = merge_par([b, a]; resolve=res)
        @test mr.domains == ["Bland", "Aland"] && mkey(mr.groups, "kryos", "domains") == ["Bland", "Aland"]
        @test sort(gnames(mr.groups)) == sort(gnames(G))
        @test mkey(mr.groups, "forcing", "coefs") == ["alpha", "dTo"]
        # written file reads back
        out = joinpath(d, "merged.nml")
        write_merged(out, m)
        @test gnames(NL.read_namelist(out)) == gnames(G)
        io = IOBuffer(); merge_report(io, m, out)
        @test occursin("`&topo` -> `&topo_aaa`", String(take!(io)))
    end

    @testset "stops" begin
        d = mktempdir()
        a = mini_domain(d, "Aland", "aaa"; tc="-21e3")
        b = mini_domain(d, "Bland", "bbb")
        msg = failure(() -> merge_par([a, b]))
        @test occursin("&clock time_const", msg) && occursin("-21e3", msg) && occursin("Aland.nml", msg) && occursin("Bland.nml", msg)
        @test occursin("clock.time_const", failure(() -> merge_par([a, b]; resolve=Dict("clock.time_end" => "1.0"))))
        @test occursin("agree", failure(() -> merge_par([a, b]; resolve=Dict("clock.time_const" => "0.0", "clock.time_end" => "1.0"))))
        a2 = mini_domain(d, "Cland", "ccc")
        b2 = mini_domain(d, "Dland", "ccc")
        @test occursin("same instance suffix", failure(() -> merge_par([a2, b2])))
        @test occursin("same domain", failure(() -> merge_par([a2, a2])))
        # a file that is already two-domain
        two = joinpath(d, "two.nml")
        m = merge_par([a2, mini_domain(d, "Eland", "eee")])
        write_merged(two, m)
        @test occursin("single-domain", failure(() -> merge_par([two, b])))
        # a renamed group no key names
        c = mini_domain(d, "Fland", "fff"; extra=["&orphan", "    x = 1", "/"])
        g = mini_domain(d, "Gland", "ggg"; extra=["&orphan", "    x = 2", "/"])
        @test occursin("named by no key", failure(() -> merge_par([c, g])))
        @test occursin("at least two", failure(() -> merge_par([a])))
    end

    @testset "the science case: NH + ANT" begin
        res = Dict("clock.time_const" => "0.0", "barysealevel.sl_path" => "\"none\"")
        nor = joinpath(MKRYOS, "par", "gate", "nor.nml"); ant = joinpath(MKRYOS, "par", "gate", "ant.nml")
        m = merge_par([nor, ant]; resolve=res)
        G = m.groups
        gd = NL.group_dict(G)
        @test m.domains == ["North", "Antarctica"] && m.suffixes == ["nor", "ant"]
        @test length(unique(gnames(G))) == length(G)
        @test mkey(G, "kryos", "domains") == ["North", "Antarctica"]
        @test mkey(G, "forcing", "coefs") == ["alpha", "dTo"]
        @test mkey(G, "clock", "time_const") === 0.0 && mkey(G, "barysealevel", "sl_path") == "none"
        @test mkey(G, "barysealevel", "method") == "fastiso"
        @test mkey(G, "North", "grid_hub") == "NH-32KM" && mkey(G, "Antarctica", "grid_hub") == "ANT-16KM"
        @test mkey(G, "North", "experiments") == "recon_nor" && mkey(G, "Antarctica", "experiments") == ""
        @test mkey(G, "ice_nor", "yelmo") == "yelmo_nor" && mkey(G, "ice_ant", "yelmo") == "yelmo_ant"
        @test mkey(G, "ycalv_nor", "calv_flt_method") == "vm-m16" && mkey(G, "ycalv_ant", "calv_flt_method") == "equil"
        @test mkey(G, "climate_nor", "reference") == "pd_nor" && mkey(G, "climate_ant", "reference") == "pd_ant"
        @test mkey(G, "climate_nor", "terms") == ["lgm", "ocean"]                   # their groups exist once
        @test mkey(G, "snapshot_pd_nor", "tas") == "era5_1981_2010_tas" && mkey(G, "snapshot_pd_ant", "tas") == "racmo23p2_1981_2010_tas"
        @test mkey(G, "isos_ant", "lithosphere") == "rheology_file" && mkey(G, "isos_nor", "lithosphere") == "uniform"
        @test mkey(G, "geo_nor", "isos") == "isos_nor" && mkey(G, "geo_ant", "isos") == "isos_ant"
        @test mkey(G, "shelf_nor", "corr") == "shelf_corr_nor" && haskey(gd, "shelf_corr_ant")
        @test mkey(G, "yelmo_nor", "nml_masks") == "yelmo_masks_nor"
        # nothing of one domain named in the other's keys: split back to the original files
        for (i, f) in enumerate((nor, ant))
            orig = NL.group_dict(NL.read_namelist(f))
            back = NL.group_dict(split_merged(m, i))
            for (n, g) in back
                widegroup(n) && continue
                @test haskey(orig, n)
                haskey(orig, n) && @test [(k.name, k.raw) for k in g.keys] == [(k.name, k.raw) for k in orig[n].keys]
            end
            @test all(n in keys(back) for n in keys(orig) if !widegroup(n))
        end
        # the merged file in par/ is this merge of the two base conversions
        base = merge_par([joinpath(MKRYOS, "par", "kryos_North.nml"), joinpath(MKRYOS, "par", "kryos_Antarctica.nml")];
                         resolve=Dict("clock.time_end" => "15e3", "clock.time_const" => "0.0",
                                      "barysealevel.sl_path" => "\"none\""))
        out = joinpath(mktempdir(), "kryos_NH-ANT.nml")
        write_merged(out, base)
        @test read(out, String) == read(joinpath(MKRYOS, "par", "kryos_NH-ANT.nml"), String)
    end
end
