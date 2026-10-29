# KryosTools.GlobeIce: the merge formulas on a hand-made map (two source cells into one target cell
# that they cover entirely) and on a partly covered target.
using KryosTools.GlobeIce

@testset "globe_ice offline merge" begin
    m = ScripMap([1, 2], [1, 1], [1.0, 1.0], 2, 1, [1.0, 1.0], [1.0, 1.0], [2.0], [2.0])
    s = OfflineSource(m, [3.0, 5.0], [10.0], [2.0, 4.0], [1.0, 1.0])
    var, F = merge_offline([s], [100.0])
    @test F ≈ [1.0]
    @test var ≈ [(3.0 * 2.0 + 5.0 * 4.0) / 10.0]                  # an amount per target area, background gone

    s2 = OfflineSource(m, [3.0, 5.0], [10.0], [2.0, 4.0], [1.0, 0.5])      # the second cell half owned
    var2, F2 = merge_offline([s2], [100.0])
    @test F2 ≈ [0.75]
    @test var2 ≈ [0.25 * 100.0 + (3.0 * 2.0 + 0.5 * 5.0 * 4.0) / 10.0]

    vi, _ = merge_offline([s], [100.0]; kind = :intensive)           # coverage x mean
    @test vi ≈ [3.0]
end
