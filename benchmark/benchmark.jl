using Chairmarks

function prettytime(sample)
    io = IOBuffer()
    Chairmarks.print_time(io, sample.time)
    String(take!(io))
end

using Modulo2

using AbstractAlgebra, Nemo, Mods, LinearAlgebraX

const R2 = Mod{2}
const R3 = AbstractAlgebra.GF(2)
const R4 = Nemo.GF(2)

skipped = "*skipped*"

function bm_mul(v)
s = """
### Matrix multiplication

`Matrix{ZZ2}` refers to Julia's standard matrix multiplication with coefficients in `ZZ2`.

| size | Modulo2 | Matrix{ZZ2} | AbstractAlgebra | Nemo |
| ---: | ---: | ---: | ---: | ---: |
"""
for n in v
    a0 = rand(Bool, n, n)
    b0 = rand(Bool, n, n)

    t1 = let a1 = ZZ2Matrix(a0), b1 = ZZ2Matrix(b0)
        prettytime(@b $a1 * $b1)
    end
    if n <= 1000
        t2 = let a2 = map(ZZ2, a0), b2 = map(ZZ2, b0)
            prettytime(@b $a2 * $b2)
        end
    else
        t2 = skipped
    end
    if n <= 500
        t3 = let a3 = matrix(R3, a0), b3 = matrix(R3, b0)
            prettytime(@b $a3 * $b3)
        end
    else
        t3 = skipped
    end
    t4 = let a4 = matrix(R4, a0), b4 = matrix(R4, b0)
        prettytime(@b $a4 * $b4)
    end

    s0 = "| $n | $t1 | $t2 | $t3 | $t4 |\n"

    print(stderr, s0)
    s *= s0
end
return s
end

function bm_rank(v)
s = """
### Rank

| size | Modulo2 | LinearAlgebraX | AbstractAlgebra | Nemo |
| ---: | ---: | ---: | ---: | ---: |
"""
for n in v
    a0 = rand(Bool, n, n)

    t1 = let a1 = ZZ2Matrix(a0)
        prettytime(@b Modulo2.rank($a1))
    end
    if n <= 500
        t2 = let a2 = map(R2, a0)
            prettytime(@b rankx($a2))
        end
    else
        t2 = skipped
    end
    if n <= 500
        t3 = let a3 = matrix(R3, a0)
            prettytime(@b AbstractAlgebra.rank($a3))
        end
    else
        t3 = skipped
    end
    t4 = let a4 = matrix(R4, a0)
        prettytime(@b AbstractAlgebra.rank($a4))
    end

    s0 = "| $n | $t1 | $t2 | $t3 | $t4 |\n"

    print(stderr, s0)
    s *= s0
end
return s
end

function bm_det(v)
s = """
### Determinant

| size | Modulo2 | LinearAlgebraX | AbstractAlgebra | Nemo |
| ---: | ---: | ---: | ---: | ---: |
"""
for n in v
    a0 = rand(Bool, n, n)

    a1 = ZZ2Matrix(a0)
    a2 = map(R2, a0)
    a3 = matrix(R3, a0)
    a4 = matrix(R4, a0)

    t1 = prettytime(@b det($a1))
    if n <= 1000
        t2 = prettytime(@b detx($a2))
    else
        t2 = skipped
    end
    if n <= 500
        t3 = prettytime(@b det($a3))
    else
        t3 = skipped
    end
    t4 = prettytime(@b det($a4))

    s0 = "| $n | $t1 | $t2 | $t3 | $t4 |\n"

    print(stderr, s0)
    s *= s0
end
return s
end

function bm_inv(v)
s = """
### Matrix inverse

| size | Modulo2 | LinearAlgebraX | AbstractAlgebra | Nemo |
| ---: | ---: | ---: | ---: | ---: |
"""
for n in v
    local a0, a1

    while true
        a0 = rand(Bool, n, n)
        a1 = ZZ2Matrix(a0)
        iszero(det(a1)) || break
    end

    a2 = map(R2, a0)
    a3 = matrix(R3, a0)
    a4 = matrix(R4, a0)

    t1 = prettytime(@b inv($a1))
    if n <= 1000
        t2 = prettytime(@b detx($a2))
    else
        t2 = skipped
    end
    if n <= 500
        t3 = prettytime(@b inv($a3))
    else
        t3 = skipped
    end
    t4 = prettytime(@b inv($a4))

    s0 = "| $n | $t1 | $t2 | $t3 | $t4 |\n"

    print(stderr, s0)
    s *= s0
end
return s
end

println("""
Mods: $(pkgversion(Mods))
LinearAlgebraX: $(pkgversion(LinearAlgebraX))
AbstractAlgebra: $(pkgversion(AbstractAlgebra))
Nemo: $(pkgversion(Nemo))
""")

v = [125, 250, 500, 1000, 2000, 4000, 8000]

println(bm_mul(v))
println(bm_rank(v))
println(bm_det(v))
println(bm_inv(v))
