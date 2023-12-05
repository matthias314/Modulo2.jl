using Modulo2, Test

using Modulo2: randomarray

@testset "ZZ2" begin
    @test false == @inferred ZZ2(0)
    @test true == @inferred ZZ2(1)
    @test ZZ2(0) == @inferred zero(ZZ2)
    @test ZZ2(1) == @inferred one(ZZ2)
    @test @inferred ZZ2(0) == ZZ2(0)
    @test ZZ2(1) == ZZ2(1)
    @test ZZ2(0) != ZZ2(1)

    @test @inferred iszero(ZZ2(0))
    @test !iszero(ZZ2(1))
    @test !isone(ZZ2(0))
    @test @inferred isone(ZZ2(1))

    for x in (ZZ2(0), ZZ2(1))
        f = @inferred iseven(x)
        @test f == iszero(x)
        f = @inferred isodd(x)
        @test f == isone(x)

        y = @inferred Base.literal_pow(^, x, Val(0))
        @test isone(y)
        y = @inferred Base.literal_pow(^, x, Val(1))
        @test y == x
        y = x^4
        @test y == x
        if iszero(x)
            @test_throws Exception x^(-3)
        else
            y = x^(-3)
            @test y == x
        end

        k = 0
        y = @inferred x^k
        @test isone(y)
        k = 1
        y = @inferred x^k
        @test y == x
        k = 4
        y = x^k
        @test y == x
        k = -3
        if iszero(x)
            @test_throws Exception x^k
        else
            y = x^k
            @test y == x
        end
    end

    for x1 in (false, true), x2 in (false, true)
        y1 = ZZ2(x1)
        y2 = ZZ2(x2)
        z = @inferred y1+y2
        @test z == xor(x1, x2)
        z = @inferred y1-y2
        @test z == xor(x1, x2)
        z = @inferred y1*y2
        @test z == (x1 && x2)
        if iszero(y2)
            @test_throws Exception y1/y2
        else
            z = @inferred y1/y2
            @test z == y1
        end
    end

    @test_throws Exception inv(ZZ2(0))
    @test ZZ2(1) == inv(ZZ2(1))
end

@testset "ZZ2 convert" begin
    for x in (ZZ2(0), ZZ2(1))
        y = @inferred Bool(x)
        @test y == isone(x)
        y = @inferred convert(Bool, x)
        @test y == isone(x)
        z = @inferred ZZ2(y)
        @test z == x
        z = @inferred convert(ZZ2, y)
        @test z == x
    end

    for T in (Int, UInt, Int8, UInt16, Float32, Float64), _ in 1:5
        x = round(T, 100*rand())
        y1 = @inferred ZZ2(x)
        y2 = @inferred convert(ZZ2, x)
        @test y1 == y2 == isodd(x)
    end

    for x in (1.5, 'a')
        @test_throws Exception ZZ2(x)
        @test_throws Exception convert(ZZ2, x)
    end
end

const maxn = 20000

@testset "ZZ2Array add" begin
    for d in 1:4, _ in 1:3
        n = round(Int, maxn^(1/d))
        dims = rand(1:n, d)
        a1 = rand(ZZ2, dims...)
        a2 = rand(ZZ2, dims...)
        b1 = ZZ2Array(a1)
        b2 = ZZ2Array(a2)
        c1 = @inferred b1+b2
        c2 = @inferred b1-b2
        @test b1 == a1 && b2 == a2 && c1 == c2 == a1+b2
    end
end

@testset "ZZ2Array mul" begin
    for d in 1:4, _ in 1:3
        n = round(Int, maxn^(1/d))
        dims = rand(1:n, d)
        a = ZZ2Array(rand(ZZ2, dims...))
        c0 = @inferred ZZ2(0)*a
        c1 = @inferred ZZ2(1)*a
        @test iszero(c0)
        @test c1 == a
    end
end

@testset "det and inv" begin
    for i in 1:8
        n = round(Int, sqrt(maxn))
        dim = i <= 3 ? i : rand(1:n)
        a = randomarray(dim, dim)
        d = @inferred det(a)
        if iszero(d)
            @test_throws Exception inv(a)
        else
            b = @inferred inv(a)
            @test isone(Matrix(a)*Matrix(b)) && isone(Matrix(b)*Matrix(a))
        end
    end
end

@testset "ZZ2Array == and copyto!" begin
    for d in 1:4, i in 1:5
        if i <= 2
            dims = fill(1, d)
            dims[i == 1 ? 1 : d] = i == d == 1 ? 1 : 255
        else
            n = round(Int, maxn^(1/d))
            dims = rand(1:n, d)
        end

        b = rand(ZZ2, dims...)
        c = ZZ2Array(b)
        @test b == c

        a1 = similar(c, dims...)
        @inferred copyto!(a1, b)
        a2 = similar(a1)
        @inferred copyto!(a2, c)
        @test b == a1 == a2

        a1 = similar(c, reverse(dims)...)
        @inferred copyto!(a1, b)
        a2 = similar(a1)
        @inferred copyto!(a2, c)
        @test a1 == a2
    end
end

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
