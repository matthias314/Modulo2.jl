# Modulo2.jl

This package provides types and functions for linear algebra modulo 2.
It defines a type `ZZ2` for integers mod 2 and `ZZ2Array` (`ZZ2Vector`, `ZZ2Matrix`)
for arrays (vectors, matrices) with elements in `ZZ2`.

## `ZZ2`

The type `ZZ2 <: Number` represents integers modulo 2. It is simply a wrapper around `Bool`
that allows to use the usual arithmetic operations `+` and `*` instead of the logical
operations `xor` and `&`.

Elements can be created from `Bool` or any other `Integer` type or via the functions `zero` and `one`.
Similarly, `Integer` types are promoted to `ZZ2`.

### Examples
```julia-repl
julia> ZZ2(1) == one(ZZ2)
true

julia> iszero(ZZ2(4))
true

julia> ZZ2(1) + 3
0

julia> typeof(ans)
ZZ2
```

## `ZZ2Array`

The type `ZZ2Array{N} <: AbstractArray{ZZ2,N}` represents `N`-dimensional arrays with elements of type `ZZ2`.
There are the abbreviations `ZZ2Vector <: AbstractVector{ZZ2}` and `ZZ2Matrix <: AbstractMatrix{ZZ2}`.

The internal representation is packed, meaning that each element only uses one bit.
However, columns are internally padded to a length that is a multiple of 256.

A `ZZ2Array` can be created from any `AbstractArray` whose elements can be converted to `ZZ2`.
One can also leave the elements undefined by using the `undef` argument.

### Examples
```julia-repl
julia> ZZ2Matrix([1 2 3; 4 5 6])
2×3 ZZ2Matrix:
 1  0  1
 0  1  0

julia> v = ZZ2Vector(undef, 2); v[1] = true; v[2] = 2.0; v
2-element ZZ2Vector:
 1
 0
```

Broadcasting is supported to some extent:
```julia-repl
julia> a = ZZ2Matrix([1 0 1; 1 1 0])
2×3 ZZ2Matrix:
 1  0  1
 1  1  0

julia> b = ones(ZZ2, 2, 3)
2×3 ZZ2Matrix:
 1  1  1
 1  1  1

julia> a .+= b
2×3 ZZ2Matrix:
 0  1  0
 0  0  1
```

## Functions related to Gaussian elimination

The package defines functions for (reduced) column echelon form (currently called `rref` instead of `rcef`),
`rank`, determinant (`det`) and inverses (`inv`). There are also the argument-modifying counterparts
`rref!`, `rank!`, `det!` and `inv!`. See the docstrings.
