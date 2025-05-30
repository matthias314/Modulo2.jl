"""
    Modulo2

This package provides types and functions for linear algebra modulo 2.
It defines a type `ZZ2` for integers mod 2 and
`ZZ2Array` (`ZZ2Vector`, `ZZ2Matrix`) for arrays (vectors, matrices)
with elements in `ZZ2`.

See also [`ZZ2`](@ref), [`ZZ2Array`](@ref).
"""
module Modulo2

import Base: show, ==, +, -, *, /, ^, inv, literal_pow,
    zero, one, iszero, isone, iseven, isodd, convert, rand, promote_rule,
    size, zeros, ones, getindex, setindex!, copy, Bool,
    count_ones, count_zeros

using Base: @propagate_inbounds, BitInteger

using Random: rand!, AbstractRNG, SamplerType

# otherwise "throw" makes code slower
@noinline throw(e) = Core.throw(e)
throw_dim(s) = throw(DimensionMismatch(s))



#
# ZZ2
#

export ZZ2

"""
    ZZ2 <: Number

A type representing integers modulo 2.

Elements can be created via the functions `zero` and `one` or
from an argument of type `Bool` or any other `Integer` type.
More generally, any argument accepted by `isodd` can be converted to `ZZ2`.
In algebraic computations with `ZZ2`, `Integer` types are promoted to `ZZ2`.
Conversely, elements of type `ZZ2` can be converted to `Bool`.

See also `Base.zero`, `Base.one`, `Base.iszero`, `Base.isone`, `Base.isodd`.

# Examples
```jldoctest
julia> ZZ2(1) == one(ZZ2)
true

julia> iszero(ZZ2(4.0))
true

julia> x = ZZ2(1) + 3
0

julia> typeof(x)
ZZ2

julia> Bool(x), convert(Bool, x)
(false, false)
```
"""
struct ZZ2 <: Number
    m::Bool
    ZZ2(m::Bool) = new(m)  # this avoids infinite recursion
end

Base.hash(a::ZZ2, h::UInt) = hash(a.m, h)

ZZ2(a::ZZ2) = a
ZZ2(x) = isinteger(x) ? ZZ2(isodd(x)) : error("cannot convert non-integer value to ZZ2")

Bool(a::ZZ2) = a.m

convert(::Type{ZZ2}, x::Number) = ZZ2(x)
convert(::Type{ZZ2}, x) = ZZ2(x)

show(io::IO, a::ZZ2) = print(io, isone(a) ? '1' : '0')

zero(::Type{ZZ2}) = ZZ2(false)
one(::Type{ZZ2}) = ZZ2(true)

iszero(a::ZZ2) = !a.m
isone(a::ZZ2) = a.m

iseven(a::ZZ2) = iszero(a)
isodd(a::ZZ2) = isone(a)

"""
    count_ones(a::ZZ2) -> Integer

Return `1` if `a` equals `ZZ2(1)` and `0` otherwise.

See also [`count_zeros(::ZZ2)`](@ref), [`count_ones(::ZZ2Array)`](@ref), `count_ones(::Integer)`.
"""
count_ones(a::ZZ2) = Int(isodd(a))

"""
    count_zeros(a::ZZ2) -> Integer

Return `1` if `a` equals `ZZ2(0)` and `0` otherwise.

See also [`count_ones(::ZZ2)`](@ref), [`count_zeros(::ZZ2Array)`](@ref), `count_zeros(::Integer)`.
"""
count_zeros(a::ZZ2) = Int(iseven(a))

+(a::ZZ2) = a
+(a::ZZ2, b::ZZ2) = ZZ2(xor(a.m, b.m))
-(a::ZZ2) = a
-(a::ZZ2, b::ZZ2) = a + b
*(a::ZZ2, b::ZZ2) = ZZ2(a.m & b.m)
/(a::ZZ2, b::ZZ2) = iszero(b) ? error("division by zero") : a
inv(a::ZZ2) = one(ZZ2)/a

function ^(a::ZZ2, n::Integer)
    if n > 0
        a
    elseif iszero(n)
        one(ZZ2)
    else
        inv(a)
    end
end

literal_pow(::typeof(^), a::ZZ2, ::Val{N}) where N = a^N

rand(rng::AbstractRNG, ::SamplerType{ZZ2}) = ZZ2(rand(rng, Bool))

promote_rule(::Type{<:Integer}, ::Type{ZZ2}) = ZZ2
promote_rule(::Type{Bool}, ::Type{ZZ2}) = ZZ2   # necessary to avoid ambiguities



#
# ZZ2Array
#

export ZZ2Array, ZZ2Vector, ZZ2Matrix,
    addcol!, swapcols!, rcef!, rcef, rref, rank, rank!,
    identity_matrix, dot, det, det!, inv!, mul!,
    zz2vector, zz2vector!

import Base: resize!, copyto!, similar, fill!, inv

using BitIntegers: UInt256, AbstractBitSigned, AbstractBitUnsigned

const GeneralBitInteger = Union{BitInteger,AbstractBitSigned,AbstractBitUnsigned}

const TA = UInt64
const TB = UInt256
const M = sizeof(TB)÷sizeof(TA)
const BA = 8*sizeof(TA)
const BB = 8*sizeof(TB)
const L = trailing_zeros(BA)
const LB = trailing_zeros(BB)

import LinearAlgebra: dot, rank, det, mul!

"""
    ZZ2Vector <: AbstractVector{ZZ2}
    ZZ2Matrix <: AbstractMatrix{ZZ2}
    ZZ2Array{N} <: AbstractArray{ZZ2,N}

An abstract vector / matrix / array type with elements of type `ZZ2`.

The internal representation is packed, meaning that each element only uses one bit.
However, columns are internally padded to a length that is a multiple of $BB.

A `ZZ2Array` can be created from any `AbstractArray` whose elements can be converted to `ZZ2`.
One can also leave the elements undefined by using the `undef` argument.

See also [`zeros`](@ref), [`ones`](@ref), [`zz2vector`](@ref), [`resize!`](@ref),
[`det`](@ref), [`inv`](@ref), [`rank`](@ref), [`rcef`](@ref), [`rref`](@ref).

# Examples
```jldoctest
julia> ZZ2Matrix([1 2 3; 4 5 6])
2×3 ZZ2Matrix:
 1  0  1
 0  1  0

julia> v = ZZ2Vector(undef, 2); v[1] = true; v[2] = 2.0; v
2-element ZZ2Vector:
 1
 0
```
"""
mutable struct ZZ2Array{N} <: AbstractArray{ZZ2,N}
    i1::Int
    const data::Array{TA,N}
    ZZ2Array{N}(i1::Integer, data::Array{TA,N}) where N = new(i1, data)
    # this avoids confusing error messages
end

function zeropad!(a::ZZ2Array{0})
    a.data[] &= TA(1)
    a
end

"""
    $(@__MODULE__).zeropad!(a::ZZ2Array{N}) where N -> a

Set the padding bits in the underlying array `a.data` to zero and return `a`.
The visible bits are not modified.

This is an internal function of the module.
"""
function zeropad!(a::ZZ2Array{N}) where N
    i1 = a.i1 & (BB-1)
    i1 == 0 && return a
    m = TB(1) << i1 - TB(1)
    data = @view reinterpret(TB, a.data)[end, ntuple(Returns(:), N-1)...]
    data .&= m
    a
end

const ZZ2Vector = ZZ2Array{1}
const ZZ2Matrix = ZZ2Array{2}

ZZ2Array{0}(::UndefInitializer, ii::Tuple{}; init = true) = ZZ2Array{0}(-1, zeros(TA))

function ZZ2Array{N}(::UndefInitializer, ii::NTuple{N,Integer}; init = true) where N
    any(<(0), ii) && error("invalid array dimensions")
    i1 = M * ((ii[1] + 1 << LB - 1) >> LB)
    data = Array{TA, N}(undef, i1, ii[2:end]...)
    a = ZZ2Array{N}(ii[1], data)
    init ? zeropad!(a) : a
end

ZZ2Array{N}(::UndefInitializer, ii::Integer...; init = true) where N = ZZ2Array{N}(undef, ii; init)

ZZ2Array(::UndefInitializer, ii::NTuple{N,Integer}; init = true) where N = ZZ2Array{N}(undef, ii; init)

ZZ2Array(::UndefInitializer, ii::Integer...; init = true) = ZZ2Array(undef, ii; init)

function ZZ2Array{N}(a::AbstractArray{T,N}) where {T,N}
    b = ZZ2Array{N}(undef, size(a))
    for i in eachindex(a, b)
        @inbounds b[i] = a[i]
    end
    b
end

ZZ2Array{N}(a::ZZ2Array{N}) where N = copy(a)

ZZ2Array(a::AbstractArray{T,N}) where {T,N} = ZZ2Array{N}(a)

similar(::Type{<:ZZ2Array}, ::Type{ZZ2}, ii::Dims) = ZZ2Array(undef, ii)
similar(::ZZ2Array, ::Type{ZZ2}, ii::Dims) = ZZ2Array(undef, ii)

"""
    resize!(a::ZZ2Vector, n::Integer) -> a

Change the length of `a` to `n` and return `a`.
In case `a` is enlarged, the new entries are undefined.
"""
function resize!(a::ZZ2Vector, n::Integer; init = true)
    n >= 0 || error("new length must be ≥ 0")
    i1 = M * ((n + 1 << LB - 1) >> LB)
    i1 == length(a.data) || resize!(a.data, i1)
    a.i1 = n
    init ? zeropad!(a) : a
end

function fill!(a::ZZ2Array, c)
    c = ZZ2(c)
    fill!(a.data, iszero(c) ? TA(0) : ~TA(0))
    isone(c) && zeropad!(a)
    a
end

# TODO: add zero_matrix ?
"""
    zeros(ZZ2, ii::NTuple{N,Integer}) where N

Return a `ZZ2Array` of size `ii` with zero entries.
"""
function zeros(::Type{ZZ2}, ii::NTuple{N,Integer}) where N
    a = ZZ2Array{N}(undef, ii; init = false)
    fill!(a, ZZ2(0))
end

zeros(::Type{ZZ2}, ::Tuple{}) = fill!(ZZ2Array{0}(undef; init = false), ZZ2(0))

"""
    ones(ZZ2, ii::NTuple{N,Integer}) where N

Return a `ZZ2Array` of size `ii` with entries `ZZ2(1)`.
"""
function ones(::Type{ZZ2}, ii::NTuple{N,Integer}) where N
    a = ZZ2Array{N}(undef, ii; init = false)
    fill!(a, ZZ2(1))
end

ones(::Type{ZZ2}, ::Tuple{}) = fill!(ZZ2Array{0}(undef; init = false), ZZ2(1))

# TODO: could probably be done more efficiently
function identity_matrix(::Type{ZZ2}, i1::Integer, i2::Integer = i1)
    a = zeros(ZZ2, i1, i2)
    for k in 1:min(i1, i2)
        a[k, k] = one(ZZ2)
    end
    a
end

zero(a::ZZ2Array) = zeros(ZZ2, size(a))

one(a::ZZ2Matrix) = identity_matrix(ZZ2, size(a)...)

size(a::ZZ2Array{0}, d) = error("dimension out of range")
size(a::ZZ2Array, d) = d == 1 ? a.i1 : size(a.data, d)

size(a::ZZ2Array{0}) = ()
size(a::ZZ2Array) = (a.i1, size(a.data)[2:end]...)

copy(a::ZZ2Array{N}) where N = ZZ2Array{N}(a.i1, copy(a.data))

function copyto!(a::ZZ2Array, b::ZZ2Array)
    if a.i1 == b.i1 || (a.i1 == BA*size(a.data, 1) && b.i1 == BA*size(b.data, 1))
        copyto!(a.data, b.data)
        a
    else
        invoke(copyto!, Tuple{AbstractArray,AbstractArray}, a, b)
    end
end

convert(::Type{ZZ2Array{N}}, a::ZZ2Array{N}) where N = a

convert(::Type{ZZ2Array{N}}, a::AbstractArray{T,N}) where {T,N} = ZZ2Array{N}(a)

convert(::Type{ZZ2Array}, a::AbstractArray{T,N}) where {T,N} = convert(ZZ2Array{N}, a)

# conversion to and from BitVector

function ZZ2Vector(v::BitVector)
    m = length(v)
    i1 = M * ((m + 1 << LB - 1) >> LB)
    w = Vector{TA}(undef, i1)
    unsafe_copyto!(w, 1, v.chunks, 1, length(v.chunks))
    zeropad!(ZZ2Vector(m, w))
end

function Base.BitVector(v::ZZ2Vector)
    m = length(v)
    w = BitVector(undef, m)
    unsafe_copyto!(w.chunks, 1, v.data, 1, Base.num_bit_chunks(m))
    w
end

"""
    zz2vector(n) -> ZZ2Vector

Return a `ZZ2Vector` containing the bit representation of `n`.
Here `n` must be a `Base.BitInteger` or a bit integer defined by the package `BitIntegers.jl`.

See also [`zz2vector!`](@ref), `Base.BitInteger`,
`BitIntegers.AbstractBitSigned`, `BitIntegers.AbstractBitUnsigned`.

# Example
```jldoctest
julia> zz2vector(Int8(35))
8-element ZZ2Vector:
 1
 1
 0
 0
 0
 1
 0
 0
```
"""
zz2vector(n::T) where T <: GeneralBitInteger = zz2vector!(ZZ2Vector(undef, 8*sizeof(T); init = false), n)

"""
    zz2vector!(a::ZZ2Vector, n) -> a

Fill `a` with the bit representation of `n` and return `a`.
Here `n` must be a `Base.BitInteger` or a bit integer defined by the package `BitIntegers.jl`.
The length of `a` is set to the bit length of `n`.

See also [`zz2vector`](@ref), `Base.BitInteger`,
`BitIntegers.AbstractBitSigned`, `BitIntegers.AbstractBitUnsigned`.
"""
function zz2vector!(a::ZZ2Vector, n::T) where T <: GeneralBitInteger
    resize!(a, 8*sizeof(T); init = false)
    for i in 1:length(a.data)
        @inbounds a.data[i] = n & ~zero(TA)
        n >>= BA
    end
    a
end

getindex(a::ZZ2Array{0}) = ZZ2(a.data[])

@inline function getindex(a::ZZ2Array{N}, ii::Vararg{Int,N}) where N
    @boundscheck checkbounds(a, ii...)
    ii1 = ii[1]-1
    i1 = (ii1 >> L) + 1
    i0 = ii1 & (1 << L - 1)
    @inbounds ZZ2(a.data[i1, ii[2:end]...] >> i0)
end

function setindex!(a::ZZ2Array{0}, x)
    a.data[] = Bool(ZZ2(x))
    a
end

@inline function setindex!(a::ZZ2Array{N}, x, ii::Vararg{Int,N}) where N
    @boundscheck checkbounds(a, ii...)
    ii1 = ii[1]-1
    i1 = (ii1 >> L) + 1
    i0 = ii1 & (1 << L - 1)
    m = TA(1) << i0
    if iszero(ZZ2(x))
        @inbounds a.data[i1, ii[2:end]...] &= ~m
    else
        @inbounds a.data[i1, ii[2:end]...] |= m
    end
    a
end

==(a::ZZ2Array, b::ZZ2Array) = a.i1 == b.i1 && a.data == b.data

"""
    count_ones(a::ZZ2Array) -> Integer

Return the number of elements of `a` equal to `ZZ2(1)`.

See also [`count_zeros(::ZZ2Array)`](@ref), [`count_ones(::ZZ2)`](@ref), `count_ones(::Integer)`.
"""
count_ones(a::ZZ2Array) = sum(count_ones, a.data; init = 0)

"""
    count_zeros(a::ZZ2Array) -> Integer

Return the number of elements of `a` equal to `ZZ2(0)`.

See also [`count_ones(::ZZ2Array)`](@ref), [`count_zeros(::ZZ2)`](@ref), `count_zeros(::Integer)`.
"""
count_zeros(a::ZZ2Array) = length(a) - count_ones(a)   # we can't sum over a.data because of the padding

function +(a::ZZ2Array{N}, b::ZZ2Array{N}) where N
    ii = size(a)
    jj = size(b)
    ii == jj || throw_dim("first array has dimensions $ii, second array has dimensions $jj")
    ZZ2Array{N}(a.i1, map(⊻, a.data, b.data))
end

# without the following methods for +, - and * one gets errors in broadcast_preserving_zero_d
+(a::ZZ2Array) = copy(a)  # see julia#58295
-(a::ZZ2Array) = copy(a)

-(a::ZZ2Array{N}, b::ZZ2Array{N}) where N = a+b

*(c::Number, a::ZZ2Array) = iszero(ZZ2(c)) ? zero(a) : copy(a)
# end of the list

*(a::ZZ2Matrix, b::AbstractVector{<:Number}) = mul!(ZZ2Vector(undef, size(a, 1); init = false), a, b)

"""
    mul!(c::ZZ2Vector, a::ZZ2Matrix, b::AbstractVector{<:Number}, α::Number = ZZ2(1), β::Number = ZZ2(0)) -> c

Store the combined matrix-vector multiply-add `α a*b + β c` in `c` and return `c`.
With the default values of `α` and `β` the product `a*b` is computed.
The elements of `b` as well as `α` and `β` must be convertible to `ZZ2`.

This function is re-exported from the module `LinearAlgebra`.
"""
function mul!(c::ZZ2Vector, a::ZZ2Matrix, b::AbstractVector{<:Number}, α::Number = ZZ2(1), β::Number = ZZ2(0))
    i1, i2 = size(a)
    j1 = size(b, 1)
    (i2 == j1 && size(c, 1) == i1) || throw_dim("matrix has dimensions ($i1, $i2), vector has length $j1")
    iszero(ZZ2(β)) && fill!(c, ZZ2(0))   # needed since we add columns to c
    if isone(ZZ2(α))
        for k in 1:i2
            @inbounds isone(ZZ2(b[k])) && addcol!(c, 1, a, k)
        end
    end
    c
end

*(a::ZZ2Matrix, b::AbstractMatrix{<:Number}) = mul!(ZZ2Array(undef, size(a, 1), size(b, 2); init = false), a, b)

"""
    mul!(c::ZZ2Matrix, a::ZZ2Matrix, b::AbstractMatrix{<:Number}, α::Number = ZZ2(1), β::Number = ZZ2(0)) -> c

Store the combined matrix-matrix multiply-add `α a*b + β c` in `c` and return `c`.
With the default values of `α` and `β` the product `a*b` is computed.
The elements of `b` as well as `α` and `β` must be convertible to `ZZ2`.

This function is re-exported from the module `LinearAlgebra`.
"""
function mul!(c::ZZ2Matrix, a::ZZ2Matrix, b::AbstractMatrix{<:Number}, α::Number = ZZ2(1), β::Number = ZZ2(0))
    i1, i2 = size(a)
    j1, j2 = size(b)
    i2 == j1 || throw_dim("first matrix has dimensions ($i1, $i2), second matrix has dimensions ($j1, $j2)")
    iszero(ZZ2(β)) && fill!(c, ZZ2(0))   # needed since we add columns to c
    if isone(ZZ2(α))
        for k2 in 1:j2, k1 in 1:j1
            @inbounds isone(ZZ2(b[k1, k2])) && addcol!(c, k2, a, k1)
        end
    end
    c
end

function dot(a::ZZ2Array, b::ZZ2Array)
    size(a) == size(b) || throw_dim("vectors/arrays must have the same size")
    s = TA(0)
    for j in eachindex(a.data)
        @inbounds s ⊻= a.data[j] & b.data[j]
    end
    ZZ2(count_ones(s))
end

const xor_ir = VERSION > v"1.12-" ? """
    %q1 = getelementptr inbounds <$M x i$BA>, ptr %0, i64 %1
    %v1 = load <$M x i$BA>, ptr %q1, align 8
    %q2 = getelementptr inbounds <$M x i$BA>, ptr %2, i64 %3
    %v2 = load <$M x i$BA>, ptr %q2, align 8
    %vr = xor <$M x i$BA> %v1, %v2
    store <$M x i$BA> %vr, ptr %q1, align 8
    ret void
""" : """
    %p1 = inttoptr i64 %0 to <$M x i$BA>*
    %q1 = getelementptr inbounds <$M x i$BA>, <$M x i$BA>* %p1, i64 %1
    %v1 = load <$M x i$BA>, <$M x i$BA>* %q1, align 8
    %p2 = inttoptr i64 %2 to <$M x i$BA>*
    %q2 = getelementptr inbounds <$M x i$BA>, <$M x i$BA>* %p2, i64 %3
    %v2 = load <$M x i$BA>, <$M x i$BA>* %q2, align 8
    %vr = xor <$M x i$BA> %v1, %v2
    store <$M x i$BA> %vr, <$M x i$BA>* %q1, align 8
    ret void
"""

function xor!(p1::Ptr{TA}, j1::Int, p2::Ptr{TA}, j2::Int)
    Base.llvmcall(xor_ir, Cvoid, Tuple{Ptr{TA},Int64,Ptr{TA},Int64}, p1, j1, p2, j2)
end

@inline function addcol!(a::ZZ2Array, k0::Integer, b::ZZ2Array, k1::Integer, range::AbstractUnitRange = axes(a.data, 1))
    i1 = size(a.data, 1)
    j1 = ((k0-1)*i1+first(range)-1) ÷ M
    j2 = ((k1-1)*i1+first(range)-1) ÷ M
    l = (length(range)+M-1) ÷ M
    for _ in 0:l-1
        xor!(pointer(a.data), j1, pointer(b.data), j2)
        j1 += 1
        j2 += 1
    end
    a
end

function swapcols!(a::ZZ2Array, k0::Integer, k1::Integer, range::AbstractUnitRange = axes(a.data, 1))
    c = a.data
    @inbounds for j in range
        c[j, k0], c[j, k1] = c[j, k1], c[j, k0]
    end
    a
end

"""
    $(@__MODULE__).gauss!(a::ZZ2Matrix, ::Val{mode}) where mode

* If `mode == :rcef`, return `(r, b)` where `r` is the rank of `a`
  and `b` its reduced column echelon form.
* If `mode == :cef`, return `(r, b)` where `r` is the rank of `a`
  and `b` a column echelon form of it.
* If `mode == :det`, return the determinant of `a`.
* If `mode == :inv`, return the inverse of `a`.

In all cases, the matrix `a` may be modified during the computation.
This is an internal function of the module.
"""
function gauss!(b::ZZ2Matrix, ::Val{mode}) where mode
# modes:
#  :rcef = reduced column echelon form (zeros left of leading ones)
#  :cef  = column echelon form (no zeros left of leading ones)
#  :det  = determinant
#  :inv  = inverse
    full = mode == :rcef || mode == :inv
    i1, i2 = size(b)
    ii1 = size(b.data, 1)
    k = 1
    if mode == :inv
        bi = identity_matrix(ZZ2, i1, i2)
    end
    for j1 in 1:i1
        flag = true
        for j2 in k:i2
            @inbounds if isone(b[j1, j2])
                if j2 != k
                    jj = (j1-1) >> L + 1
                    swapcols!(b, j2, k, jj:ii1)
                    mode == :inv && swapcols!(bi, j2, k)
                end
                for l in (full ? 1 : k+1):i2
                    if (!full || l != k) && isone(b[j1, l])
                        jj = (j1-1) >> L + 1
                        addcol!(b, l, b, k, jj:ii1)
                        mode == :inv && addcol!(bi, l, bi, k)
                    end
                end
                k += 1
                flag = false
                break
            end
        end
        if mode == :det && flag
            return ZZ2(0)
        elseif mode == :inv && flag
            error("matrix not invertible")
        end
    end
    if mode == :det
        return ZZ2(1)
    elseif mode == :inv
        return bi
    else
        return (k-1, b)
    end
end

"""
    rcef!(b::ZZ2Matrix; reduced = true) -> ZZ2Matrix

Return the tuple `(r, c)` where `r` is the rank of the matrix `b` and `c` a column echelon form of it.
If `reduced` is `true`, then the reduced column echelon form is computed.
The argument `b` may be modified during the computation, which avoids the allocation of a new matrix.

See also [`rcef`](@ref).
"""
rcef!(b::ZZ2Matrix; reduced = true) = gauss!(b, Val(reduced ? :rcef : :cef))

"""
    rcef(b::ZZ2Matrix; reduced = true) -> ZZ2Matrix

Return the tuple `(r, c)` where `r` is the rank of the matrix `b` and `c` a column echelon form of it.
If `reduced` is `true`, then the reduced column echelon form is computed.

See also [`rref`](@ref), [`rcef!`](@ref).

# Examples
```jldoctest
julia> a = ZZ2Matrix([1 0 0; 1 1 1])
2×3 ZZ2Matrix:
 1  0  0
 1  1  1

julia> rcef(a)
(2, ZZ2[1 0 0; 0 1 0])

julia> rcef(a; reduced = false)
(2, ZZ2[1 0 0; 1 1 0])
```
"""
rcef(b::ZZ2Matrix; kw...) = rcef!(copy(b); kw...)

"""
    rref(b::ZZ2Matrix; reduced = true) -> ZZ2Matrix

Return the tuple `(r, c)` where `r` is the rank of the matrix `b` and `c` a row echelon form of it.
If `reduced` is `true`, then the reduced row echelon form is computed.

Note that it is more efficient to compute a (reduced) *column* echelon form via `rcef`.

See also [`rcef`](@ref).

# Examples
```jldoctest
julia> a = ZZ2Matrix([1 1; 0 1; 0 1])
3×2 ZZ2Matrix:
 1  1
 0  1
 0  1

julia> rref(a)
(2, ZZ2[1 0; 0 1; 0 0])

julia> rref(a; reduced = false)
(2, ZZ2[1 1; 0 1; 0 0])
```
"""
function rref(b::ZZ2Matrix; kw...)
    r, c = rcef!(ZZ2Matrix(transpose(b)); kw...)
    r, transpose(c)
end

"""
    rank!(b::ZZ2Matrix) -> Int

Return the rank of the matrix `b`.
The argument `b` may be modified during the computation, which avoids the allocation of a new matrix.

See also [`rank`](@ref).
"""
rank!(b::ZZ2Matrix) = gauss!(b, Val(:cef))[1]

"""
    rank(b::ZZ2Matrix) -> Int

Return the rank of the matrix `b`.

See also [`rank!`](@ref).
"""
rank(b::ZZ2Matrix) = rank!(copy(b))

"""
    det!(b::ZZ2Matrix) -> ZZ2

Return the determinant of the matrix `b`.
The argument `b` may be modified during the computation, which avoids the allocation of a new matrix.

See also [`det`](@ref).
"""
function det!(b::ZZ2Matrix)
    ==(size(b)...) || throw_dim("matrix is not square")
    gauss!(b, Val(:det))
end

"""
    det(b::ZZ2Matrix) -> ZZ2

Return the determinant of the matrix `b`.

See also [`det!`](@ref).
"""
det(b::ZZ2Matrix) = det!(copy(b))

"""
    inv!(b::ZZ2Matrix) -> ZZ2Matrix

Return the inverse of the matrix `b`, which must be invertible.
The argument `b` may be modified during the computation, which avoids the allocation of a new matrix.

See also [`inv`](@ref).
"""
function inv!(b::ZZ2Matrix)
    ==(size(b)...) || throw_dim("matrix is not square")
    gauss!(b, Val(:inv))
end

"""
    inv(b::ZZ2Matrix) -> ZZ2Matrix

Return the inverse of the matrix `b`, which must be invertible.

See also [`inv!`](@ref).
"""
inv(b::ZZ2Matrix) = inv!(copy(b))

function randommatrix(i1, i2, k)
    a = zeros(ZZ2, i1, i2)
    for j in 1:k
        j1 = rand(1:i1)
        j2 = rand(1:i2)
        a[j1, j2] = 1
    end
    a
end

"""
    $(@__MODULE__).randomarray(ii...) -> ZZ2Array

Return a `ZZ2Array` of size `ii` with random entries.
"""
function randomarray(ii...)
    a = ZZ2Array(undef, ii; init = false)
    rand!(a.data)
    zeropad!(a)
end



#
# broadcasting
#

import Base: copy, copyto!

using Base.Broadcast: AbstractArrayStyle, DefaultArrayStyle, Broadcasted
import Base.Broadcast: BroadcastStyle

struct ZZ2ArrayStyle{N} <: AbstractArrayStyle{N} end

(::Type{<:ZZ2ArrayStyle})(::Val{N}) where N = ZZ2ArrayStyle{N}()

BroadcastStyle(::Type{ZZ2Array{N}}) where N = ZZ2ArrayStyle{N}()

similar(bc::Broadcasted{ZZ2ArrayStyle{N}}, ::Type{ZZ2}, dims) where N = similar(ZZ2Array{N}, dims)

function add!(a::ZZ2Array, b::ZZ2Array)
    if size(a) != size(b)
        throw_dim("first array has dimensions $(size(a)), second array has dimensions $(size(b))")
    end
    a.data .⊻= b.data
    a
end

function add!(a::ZZ2Array{N}, b::AbstractArray{T,N}) where {T,N}
    for ii in eachindex(a, b)
        @inbounds a[ii] += b[ii]
    end
    a
end

add!(a::ZZ2Array{N}, bc::Broadcasted{ZZ2ArrayStyle{N}, <:Any, <:Union{typeof(+), typeof(-)}}) where N =
    foldl(add!, bc.args; init = a)

function add!(a::ZZ2Array{N}, bc::Broadcasted{ZZ2ArrayStyle{N}, <:Any, typeof(*)}) where N
    a1, a2 = bc.args
    iszero(ZZ2(a1)) ? a : add!(a, a2)
end

copy_convert(a::AbstractArray{T,N}) where {T,N} = convert(ZZ2Array{N}, a)
copy_convert(a::Union{ZZ2Array,Broadcasted}) = copy(a)

function copy(bc::Broadcasted{<:ZZ2ArrayStyle, <:Any, <:Union{typeof(+), typeof(-)}})
    if bc.args isa Tuple{ZZ2Array,ZZ2Array}
        +(bc.args...)
    else
        a1, as... = bc.args
        foldl(add!, as; init = copy_convert(a1))
    end
end

function copy(bc::Broadcasted{ZZ2ArrayStyle{N}, <:Any, typeof(*)}) where N
    c, b = bc.args[1] isa Number ? bc.args : reverse(bc.args)
    iszero(ZZ2(c)) ? fill!(similar(b, ZZ2), ZZ2(0)) : copy(b)
end

function copyto!(a::ZZ2Array{N}, bc::Broadcasted{ZZ2ArrayStyle{N}, <:Any, <:Union{typeof(+), typeof(-)}}) where N
    a1, as... = bc.args
    foldl(add!, as; init = a === a1 ? a : copyto!(a, a1))
end

function copyto!(a::ZZ2Array{N}, bc::Broadcasted{ZZ2ArrayStyle{N}, <:Any, typeof(*)}) where N
    c, b = bc.args[1] isa Number ? bc.args : reverse(bc.args)
    if iszero(ZZ2(c))
        fill!(a, ZZ2(0))
    elseif a !== b
        copyto!(a, b)
    end
    a
end



#
# precompilation
#

for i in (:rcef, :cef, :det, :inv)
    precompile(gauss!, (ZZ2Matrix, Val{i}))
end

end
