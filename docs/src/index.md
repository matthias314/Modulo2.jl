```@meta
DocTestSetup = quote
        using Modulo2
    # for jldoctest outside of docstrings
    end
```

# Modulo2.jl

```@docs
Modulo2
```

## Type `ZZ2`

```@docs
ZZ2
count_zeros(::ZZ2)
count_ones(::ZZ2)
```

## Type `ZZ2Array`

```@docs
ZZ2Array
resize!
zeros
ones
count_zeros(::ZZ2Array)
count_ones(::ZZ2Array)
zz2vector
zz2vector!
Modulo2.randomarray
mul!
rcef
rcef!
rref
rank
rank!
det
det!
inv
inv!
```

## Broadcasting

Broadcasting is partially implemented for the type `ZZ2Array`,
namely for addition, subtraction and scalar multiplication as well
as for assignments.
```jldoctest
julia> a, b = ZZ2Matrix([1 0; 1 1]), ZZ2Matrix([0 1; 1 0])
(ZZ2[1 0; 1 1], ZZ2[0 1; 1 0])

julia> c = ZZ2(1)
1

julia> a .+ c .* b
2×2 ZZ2Matrix:
 1  1
 0  1

julia> a .*= c
2×2 ZZ2Matrix:
 1  0
 1  1

julia> a .= c .* b
2×2 ZZ2Matrix:
 0  1
 1  0

julia> a .+= c .* b
2×2 ZZ2Matrix:
 0  0
 0  0
```

## Internal functions

```@docs
Modulo2.zeropad!
Modulo2.gauss!
```
