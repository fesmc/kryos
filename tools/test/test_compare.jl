using KryosTools.Compare
using NCDatasets
using TOML

# a pair of run directories with one file: `x` (Float32, 3x2, three records); `bump` moves one value of the last record by 1 ulp
function make_run(dir; bump=false, name="a.nc", var="x", scale=1.0, units=nothing, T=Float32)
    mkpath(dir)
    NCDataset(joinpath(dir, name), "c") do ds
        defDim(ds, "xc", 3); defDim(ds, "yc", 2); defDim(ds, "time", Inf)
        t = defVar(ds, "time", Float64, ("time",)); xc = defVar(ds, "xc", Float64, ("xc",)); yc = defVar(ds, "yc", Float64, ("yc",))
        x = defVar(ds, var, T, ("xc", "yc", "time"))
        units === nothing || (x.attrib["units"] = units)
        xc[:] = [1.0, 2.0, 3.0]; yc[:] = [1.0, 2.0]
        for k in 1:3
            t[k] = 10.0 * (k - 1)
            v = T.(reshape(1:6, 3, 2)) .* T(k)
            v = v .* T(scale)
            bump && k == 3 && (v[2, 1] = nextfloat(v[2, 1]))
            x[:, :, k] = v
        end
    end
    return dir
end

@testset "Compare: identity" begin
    d = mktempdir()
    ref = make_run(joinpath(d, "ref")); same = make_run(joinpath(d, "same")); off = make_run(joinpath(d, "off"); bump=true)
    r = compare(ref, same)
    @test exit_code(r) == 0 && r.first_divergence === nothing
    r = compare(ref, off)
    @test exit_code(r) == 1
    @test r.first_divergence["time"] == 20.0
    @test r.first_divergence["variables"][1]["key"] == "a.nc:x"
    io = IOBuffer(); report_markdown(io, r)
    @test occursin("FAIL", String(take!(io)))
    @test occursin("\"exit\":1", report_json(r))
    # a missing run directory is an error (exit 2)
    @test exit_code(compare(ref, joinpath(d, "nowhere"))) == 2
end

@testset "Compare: mapped, class C" begin
    d = mktempdir()
    ref = make_run(joinpath(d, "ref"); name="yelmo.nc", var="smb")
    good = make_run(joinpath(d, "good"); name="ice.nc", var="smb", scale=910.0, units="kg m-2 yr-1", T=Float64)
    # the run holds Float64(x) * 910, as the converter of the adapters writes it
    NCDataset(joinpath(good, "ice.nc"), "a") do ds
        x = Float64.(ds["smb"][:, :, :]) ./ 910.0
        ds["smb"][:, :, :] = Float64.(Float32.(x)) .* 910.0
    end
    map = joinpath(d, "map.toml")
    write(map, """
    version = 1
    [constants]
    rho_ice = 910.0
    [[row]]
    ref = "yelmo.nc"
    run = "ice.nc"
    class = "C"
    expr = "Float64(x) * rho_ice"
    units = "kg m-2 yr-1"
    comp = "ice"
    """)
    r = compare(ref, good; map=map)
    @test exit_code(r) == 0
    bad = make_run(joinpath(d, "bad"); name="ice.nc", var="smb", scale=910.0, units="kg m-2", T=Float64)
    NCDataset(joinpath(bad, "ice.nc"), "a") do ds
        x = Float64.(ds["smb"][:, :, :]) ./ 910.0
        ds["smb"][:, :, :] = Float64.(Float32.(x)) .* 910.0
    end
    @test exit_code(compare(ref, bad; map=map)) == 1                 # wrong units attribute
end

@testset "Compare: the gate files" begin
    for n in ("grl", "ant")
        g = load_gate(n)
        @test g["gate"] == n && !isempty(g["instances"]) && g["chunks"][end] == g["time_end"]
    end
    m = load_map(joinpath(@__DIR__, "..", "gate", "map_phase1.toml"))
    @test !isempty(m.rows) && !isempty(m.restart_rows)
    @test isempty(load_deviations(joinpath(@__DIR__, "..", "gate", "deviations.toml")))
end
