using KryosTools.Compare
using KryosTools.Gate
using NCDatasets
import TOML

# Synthetic run folders of a two-domain gate: files named as kryos writes them (flat instance files, ts files,
# bundles with a folder per domain), enough for the checks.
function kk_write_nc(path, vars::Pair...; time=[0.0, 10.0, 20.0, 30.0])
    mkpath(dirname(path))
    NCDataset(path, "c") do ds
        defDim(ds, "time", length(time))
        v = defVar(ds, "time", Float64, ("time",)); v[:] = time
        for (n, a) in vars
            if ndims(a) == 1
                x = defVar(ds, n, Float32, ("time",)); x[:] = a
            else
                n1, n2 = size(a, 1), size(a, 2)
                haskey(ds.dim, "x") || (defDim(ds, "x", n1); defDim(ds, "y", n2))
                x = defVar(ds, n, Float32, ("x", "y", "time")); x[:, :, :] = a
            end
        end
    end
end

function kk_run(dir; doms=("Aland",), sfx=("aaa",), bump=0.0f0, budget_err=0.0f0, mice=1.0f0)
    A = Float32[3.6f14, 3.6f14, 3.6f14, 3.6f14]
    dv = Dict(s => Float32[0, 1f12 * (i + 1), 2f12 * (i + 1), 3f12 * (i + 1)] for (i, s) in enumerate(sfx))
    bsl = zeros(Float32, 4)
    for k in 2:4
        tot = sum(Float32[dv[s][k] for s in sfx])
        bsl[k] = bsl[k-1] - tot / A[k]
    end
    bsl[4] += budget_err
    kk_write_nc(joinpath(dir, "geo_glob_ts.nc"), "bsl" => bsl, "A_ocn" => A)
    kk_write_nc(joinpath(dir, "forcing_ts.nc"), "M_ice" => Float32[mice, 2mice, 3mice, 4mice])
    for (i, (d, s)) in enumerate(zip(doms, sfx))
        H = fill(Float32(i), 3, 3, 4); H[1, 1, 3:4] .+= bump
        kk_write_nc(joinpath(dir, "ice_$s.nc"), "H" => H)
        kk_write_nc(joinpath(dir, "geo_$(s)_ts.nc"), "dV_bsl" => dv[s])
        for t in (20, 30)
            b = joinpath(dir, "restart-$t")
            kk_write_nc(joinpath(b, d, "ice_$(s)_restart.nc"), "H" => H[:, :, 1:1]; time=[Float64(t)])
            kk_write_nc(joinpath(b, "forcing_imports.nc"), "M_ice" => Float32[mice]; time=[Float64(t)])
            kk_write_nc(joinpath(b, "clock_restart.nc"), "ticks" => Float32[t]; time=[Float64(t)])
        end
    end
end

function kk_gate(; checks)
    g = Dict{String,Any}("gate" => "kk-test", "kind" => "kryos", "chunks" => [0.0, 20.0, 30.0], "restart_times" => [20.0, 30.0],
                         "domain_order" => ["Aland", "Bland"], "skip" => ["forcing_ts.nc:*", "*forcing_imports.nc:*", "geo_glob_ts.nc:*"],
                         "domains" => Dict("Aland" => Dict("instances" => Dict("domain" => "Aland", "ice" => "ice_aaa", "geo" => "geo_aaa")),
                                           "Bland" => Dict("instances" => Dict("domain" => "Bland", "ice" => "ice_bbb", "geo" => "geo_bbb"))),
                         "run" => [Dict("name" => "cont", "par" => "sci", "kp" => ["barysealevel.method=\"const\""]),
                                   Dict("name" => "swap", "par" => "sci", "swap" => true)],
                         "check" => checks)
    return g
end

@testset "Gate: kryos against kryos" begin
    d = mktempdir()
    kk_run(joinpath(d, "a"); doms=("Aland",), sfx=("aaa",))
    kk_run(joinpath(d, "ab"); doms=("Aland", "Bland"), sfx=("aaa", "bbb"), mice=2.0f0)       # the forcing sees both domains
    kk_run(joinpath(d, "ab2"); doms=("Aland", "Bland"), sfx=("aaa", "bbb"), mice=2.0f0)
    kk_run(joinpath(d, "abx"); doms=("Aland", "Bland"), sfx=("aaa", "bbb"), mice=2.0f0, bump=1f-7)
    kk_run(joinpath(d, "abb"); doms=("Aland", "Bland"), sfx=("aaa", "bbb"), mice=2.0f0, budget_err=1f-6)

    @testset "stage" begin
        g = kk_gate(; checks=[])
        @test kk_stage(g, "cont") == ("sci", ["barysealevel.method=\"const\""])
        @test kk_stage(g, "swap") == ("sci", ["kryos.domains=\"Bland\" \"Aland\""])
        @test_throws ErrorException kk_stage(g, "nothing")
    end

    @testset "per-domain equality, aggregates skipped" begin
        g = kk_gate(; checks=[Dict("type" => "domain", "domain" => "Aland", "ref" => "a", "run" => "ab")])
        r = only(run_kk_checks(g, d))
        @test r.outcome == "PASS"
        rep = r.report
        @test count(x -> x.status == "skipped", rep.results) > 0                  # M_ice of the forcing, file and bundle
        @test any(f -> startswith(f, "ice_bbb"), rep.new_files)                  # the other domain's files: N
        # the other domain's folder of the bundle is not compared, this domain's is
        @test any(x -> x.key == "restart-20/Aland/ice_aaa_restart.nc:H" && x.status == "pass", rep.results)
        # a difference in the domain fails, and the report names it
        g2 = kk_gate(; checks=[Dict("type" => "domain", "domain" => "Aland", "ref" => "a", "run" => "abx")])
        r2 = only(run_kk_checks(g2, d))
        @test r2.outcome == "FAIL"
        @test any(x -> x.status == "fail" && occursin("ice_aaa", x.key), r2.report.results)
        # without the skip the forcing aggregate fails
        g3 = kk_gate(; checks=[Dict("type" => "domain", "domain" => "Aland", "ref" => "a", "run" => "ab")])
        delete!(g3, "skip")
        @test only(run_kk_checks(g3, d)).outcome == "FAIL"
    end

    @testset "identical, budget" begin
        g = kk_gate(; checks=[Dict("type" => "identical", "ref" => "ab", "run" => "ab2"),
                              Dict("type" => "budget", "run" => "ab"),
                              Dict("type" => "budget", "run" => "abb"),
                              Dict("type" => "budget", "run" => "nope"),
                              Dict("type" => "bogus")])
        rs = run_kk_checks(g, d)
        @test [r.outcome for r in rs] == ["PASS", "PASS", "FAIL", "ERROR", "ERROR"]
        @test occursin("closes" , "closes") && occursin("records", rs[2].detail)
        b = bsl_budget(joinpath(d, "ab"), ["geo_aaa", "geo_bbb"])
        @test b.ok && length(b.times) == 3 && b.max_ratio <= b.tol
        # one domain's dV left out: the budget does not close
        b1 = bsl_budget(joinpath(d, "ab"), ["geo_aaa"])
        @test !b1.ok
        # a domain with a slower cadence counts at its own times only
        bd = bsl_budget(joinpath(d, "ab"), ["geo_aaa", "geo_bbb"]; dt=Dict("geo_bbb" => 20.0))
        @test !bd.ok
        @test kk_outcome(rs) == "ERROR" && kk_outcome(rs[1:2]) == "PASS" && kk_outcome(rs[1:3]) == "FAIL"
        io = IOBuffer(); kk_markdown(io, g, rs; runs=d)
        @test occursin("| budget ab | PASS |", String(take!(io)))
    end

    @testset "chunks" begin
        # cont = a run to 30; the chunk run c00 = 0..20, c01 = 20..30 with the records of its interval
        full = joinpath(d, "cont"); kk_run(full)
        c00 = joinpath(d, "chunk", "c00"); c01 = joinpath(d, "chunk", "c01")
        kk_run(c00); kk_run(c01)
        # kk_run writes four records and bundles at 20 and 30 in both: equal to cont's, as a chunk run's would be
        g = kk_gate(; checks=[Dict("type" => "chunks", "ref" => "cont", "run" => "chunk")])
        @test only(run_kk_checks(g, d)).outcome == "PASS"
        rm(joinpath(d, "chunk", "c01"); recursive=true); kk_run(c01; bump=1f-7)
        @test only(run_kk_checks(g, d)).outcome == "FAIL"
    end

    @testset "the gate files of the repository" begin
        gd = joinpath(@__DIR__, "..", "gate")
        for f in readdir(gd)
            endswith(f, ".toml") || continue
            g = TOML.parsefile(joinpath(gd, f))
            get(g, "kind", "") == "kryos" || continue
            @test haskey(g, "domains") && haskey(g, "domain_order") && !isempty(kk_checks(g))
            runs = Set(r["name"] for r in kk_runs(g))
            for c in kk_checks(g), k in ("ref", "run")
                haskey(c, k) && @test c[k] in runs
            end
            for r in kk_runs(g)
                par, p = kk_stage(g, r["name"])
                @test isfile(joinpath(@__DIR__, "..", "..", "par", "gate", par * ".nml"))
            end
        end
    end
end
