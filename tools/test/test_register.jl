using Test

@testset "register_k3.sh" begin
    script = normpath(joinpath(@__DIR__, "..", "levante", "register_k3.sh"))
    G = mktempdir()
    sha, stack = "abc1234", "S5b"
    function mkrun(g)
        d = joinpath(G, "runs", g, "kx-$sha", stack, "cont")
        mkpath(joinpath(d, "restart-1000", "North")); mkpath(joinpath(d, "input"))
        write(joinpath(d, "ice_nor.nc"), "x"); write(joinpath(d, "restart-1000", "North", "ice_nor_restart.nc"), "y")
        write(joinpath(d, "run.log"), "log"); write(joinpath(d, "kryos.x"), "bin"); write(joinpath(d, "par.nml"), "&kryos\n/\n")
        write(joinpath(d, "input", "link.nml"), "z")
    end
    run_(args...) = (p = run(ignorestatus(pipeline(setenv(`bash $script $args`, "G" => G); stdout=devnull, stderr=devnull))); p.exitcode)
    mkrun("sci"); mkrun("sci-const")
    @test run_(sha, stack) == 1                                           # no PASS in the reports: refuses
    mkpath(joinpath(G, "reports", "KK"))
    write(joinpath(G, "reports", "KK", "summary.txt"), "sci | $sha | $stack | # sci (kryos against kryos): PASS | exit 0\n")
    @test run_(sha, stack) == 1                                           # one gate only
    open(joinpath(G, "reports", "KK", "summary.txt"), "a") do io
        println(io, "sci-const | $sha | $stack | # sci-const (kryos against kryos): FAIL | exit 1")
    end
    @test run_(sha, stack) == 1
    open(joinpath(G, "reports", "KK", "summary.txt"), "a") do io
        println(io, "sci-const | $sha | $stack | # sci-const (kryos against kryos): PASS | exit 0")
    end
    @test run_(sha, stack) == 0
    dest = joinpath(G, "ref", "K3-sci", "kx-$sha", stack)
    @test isfile(joinpath(dest, "sci", "cont", "restart-1000", "North", "ice_nor_restart.nc"))
    @test isfile(joinpath(dest, "sci-const", "cont", "par.nml")) && isfile(joinpath(dest, "registered.txt"))
    @test !ispath(joinpath(dest, "sci", "cont", "kryos.x")) && !ispath(joinpath(dest, "sci", "cont", "run.log")) &&
          !ispath(joinpath(dest, "sci", "cont", "input"))
    @test !iswritable(joinpath(dest, "sci", "cont", "ice_nor.nc"))
    @test run_(sha, stack) == 1                                           # exists: never overwritten
    chmod(dest, 0o755; recursive=true)                                    # let the temp dir be removed
end
