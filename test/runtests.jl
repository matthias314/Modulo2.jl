using Modulo2, Test

@testset "broadcast" begin
    u = ZZ2Vector([1,1,0])
    v = ZZ2Vector([1,0,1])
    w = ZZ2Vector([0,1,1])
    @test v .+ w == v + w
    @test u .+ v .+ w == u + v + w
    @test u .+ (v .+ w) == u + v + w
    @test (u .+ v) .+ w == u + v + w
    @test ZZ2(1) .* v == v
    @test ZZ2(0) .* v == zero(v)
    @test ZZ2(1) .* (v .+ w) == v + w
    @test ZZ2(0) .* (v .+ w) == zero(v)
    @test ZZ2(1) .* v .+ w == v + w
    @test ZZ2(0) .* v .+ w == w
    @test v .+ ZZ2(1) .* w == v + w
    @test v .+ ZZ2(0) .* w == v
    @test ZZ2(0) .* v .+ ZZ2(1) .* w == w
    @test ZZ2(1) .* v .+ ZZ2(0) .* w == v
    @test ZZ2(1) .* (v .+ ZZ2(0) .* w) == v
end
