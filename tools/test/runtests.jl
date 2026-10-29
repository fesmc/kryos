using Test
using KryosTools

const KRYOS_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include("mini_yelmox.jl")

@testset "KryosTools" begin
    include("test_compare.jl")
    include("test_convert.jl")
    include("test_merge.jl")
    include("test_gate.jl")
    include("test_globe_ice.jl")
end
