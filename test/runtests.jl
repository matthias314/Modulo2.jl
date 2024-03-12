using Modulo2, Test

using Modulo2: randomarray

using BitIntegers

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
        x = T(rand(0:100))
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

@testset "ZZ2Array" begin
    for n in 0:3
        a = ZZ2Array(undef, ntuple(Returns(0), n)...)
        @test all(iszero, size(a))
        a = zeros(ZZ2, ntuple(i -> 5*i, n)...)
        @test iszero(a)
        a = ones(ZZ2, ntuple(i -> 5*i, n)...)
        @test all(isone, a)
    end
end

@testset "resize!" begin
    a = ones(ZZ2, 300)
    b = Vector(a)
    resize!(a, 10)
    @test a == resize!(b, 10)
    @test length(a.data) == Modulo2.M
    for i in 1:10
        a[i] = ZZ2(0)
    end
    @test iszero(a.data)
end

@testset "ZZ2Array add" begin
    b1 = fill!(ZZ2Array{0}(undef), ZZ2(0))
    b2 = fill!(ZZ2Array{0}(undef), ZZ2(1))
    c1 = @inferred b1+b2
    c2 = @inferred b1-b2
    @test c1 == c2
    @test c1[] == c2[] == b1[]+b2[]

    for d in 1:4, _ in 1:3
        n = round(Int, maxn^(1/d))
        dims = rand(1:n, d)
        a1 = rand(ZZ2, dims...)
        a2 = rand(ZZ2, dims...)
        b1 = ZZ2Array(a1)
        b2 = ZZ2Array(a2)
        c1 = @inferred b1+b2
        c2 = @inferred b1-b2
        @test b1 == a1 && b2 == a2 && c1 == c2 == a1+a2
    end
end

@testset "ZZ2Array mul" begin
    a = fill!(ZZ2Array{0}(undef), ZZ2(1))
    c0 = @inferred ZZ2(0)*a
    c1 = @inferred ZZ2(1)*a
    @test iszero(c0)
    @test c1 == a

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

@testset "ZZ2Matrix mat-vec mul" begin
    for n in (0, 1, 100, 257)
        a = randomarray(n, 2*n)
        b = randomarray(2*n)
        c = ones(ZZ2, n)

        ab = Matrix(a) * Vector(b)
        @test a*b == ab

        mul!(c, a, b, ZZ2(0), ZZ2(1))
        @test c == ones(ZZ2, n)
        mul!(c, a, b, ZZ2(0), ZZ2(0))
        @test iszero(c)
        mul!(c, a, b)
        @test c == ab
        fill!(c, ZZ2(1))
        mul!(c, a, b, 1, 1)
        @test c == ab + ones(ZZ2, n)
    end
end

@testset "ZZ2Matrix mat-mat mul" begin
    for n in (0, 1, 100, 255)
        a = randomarray(n, n+5)
        b = randomarray(n+5, n+7)
        c = ones(ZZ2, n, n+7)

        ab = Matrix(a) * Matrix(b)
        @test a*b == ab

        mul!(c, a, b, ZZ2(0), ZZ2(1))
        @test c == ones(ZZ2, n, n+7)
        mul!(c, a, b, ZZ2(0), ZZ2(0))
        @test iszero(c)
        mul!(c, a, b)
        @test c == ab
        fill!(c, ZZ2(1))
        mul!(c, a, b, 1, 1)
        @test c == ab + ones(ZZ2, n, n+7)
    end
end

@testset "det and inv" begin
    for i in 1:8
        n = round(Int, sqrt(maxn))
        dim = i <= 3 ? i-1 : rand(1:n)
        a = randomarray(dim, dim)
        d = @inferred det(a)
        dim == 0 && @test isone(d)
        if iszero(d)
            @test_throws Exception inv(a)
        else
            b = @inferred inv(a)
            @test isone(a*b)
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

        a1 = similar(b, dims...)
        @inferred copyto!(a1, c)
        @test c == a1

        a1 = similar(b, reverse(dims)...)
        @inferred copyto!(a1, c)
        a2 = similar(a1)
        copyto!(a2, b)
        @test a1 == a2
    end
end

@testset "broadcast" begin
    for d in 0:4
        dims = rand(2:4, d)
        u = randomarray(dims...)
        v = randomarray(dims...)
        w = randomarray(dims...)
        @test v .+ w == v + w
        @test u .+ v .+ w == u + v + w
        @test u .+ (v .+ w) == u + v + w
        @test (u .+ v) .+ w == u + v + w
        @test ZZ2(1) .* v == v
        @test ZZ2(0) .* v == zero(v)
        @test v .* ZZ2(1) == v
        @test v .* ZZ2(0) == zero(v)
        @test ZZ2(1) .* (v .+ w) == v + w
        @test ZZ2(0) .* (v .+ w) == zero(v)
        @test (v .+ w) .* ZZ2(1) == v + w
        @test (v .+ w) .* ZZ2(0)  == zero(v)
        @test ZZ2(1) .* v .+ w == v + w
        @test ZZ2(0) .* v .+ w == w
        @test v .+ ZZ2(1) .* w == v + w
        @test v .+ ZZ2(0) .* w == v
        @test ZZ2(0) .* v .+ ZZ2(1) .* w == w
        @test ZZ2(1) .* v .+ ZZ2(0) .* w == v
        @test ZZ2(1) .* (v .+ ZZ2(0) .* w) == v

        x = similar(u)
        x .= v .+ w
        @test x == v + w
        x .= u .+ v .+ w
        @test x == u + v + w
        x .= u .+ (v .+ w)
        @test x == u + v + w
        x .= (u .+ v) .+ w
        @test x == u + v + w
        x .= ZZ2(1) .* v
        @test x == v
        x .= ZZ2(0) .* v
        @test x == zero(v)
        x .= v .* ZZ2(1)
        @test x == v
        x .= v .* ZZ2(0)
        @test x == zero(v)
        x .= ZZ2(1) .* (v .+ w)
        @test x == v + w
        x .= ZZ2(0) .* (v .+ w)
        @test x == zero(v)
        x .= ZZ2(1) .* v .+ w
        @test x == v + w
        x .= ZZ2(0) .* v .+ w
        @test x == w
        x .= v .+ ZZ2(1) .* w
        @test x == v + w
        x .= v .+ ZZ2(0) .* w
        @test x == v
        x .= ZZ2(0) .* v .+ ZZ2(1) .* w
        @test x == w
        x .= ZZ2(1) .* v .+ ZZ2(0) .* w
        @test x == v
        x .= ZZ2(1) .* (v .+ ZZ2(0) .* w)
        @test x == v

        x = zero(u)
        x .+= u
        @test x == u
        x = zero(u)
        x .-= u
        @test x == u
        x = copy(u)
        x .*= ZZ2(0)
        @test iszero(x)
        x = copy(u)
        x .*= ZZ2(1)
        @test x == u
    end
end

@testset "ZZ2 count 1/0" begin
    for d in 0:4, _ in 1:4
        dims = rand(0:129, d)
        x = randomarray(dims...)
        m = @inferred count_ones(x)
        @test m == sum(count_ones, x; init = 0)
        @test m == sum(BitArray(x))
        @test m + count_zeros(x) == length(x)
    end
end

BitIntegers.@define_integers 80
BitIntegers.@define_integers 800

@testset "zz2vector" begin
    for T in (UInt8, Int16, UInt64, Int80, UInt256, Int800)
        n = rand(T)
        a = @inferred zz2vector(n)
        m = 8*sizeof(T)
        @test length(a) == m
        @test all(a[i+1] == isodd(n >> i) for i in 0:m-1)
        b = randomarray(111)
        c = zz2vector!(b, n)
        @test c === b == a
    end
end
