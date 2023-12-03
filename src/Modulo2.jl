module Modulo2

export ZZ2, ZZ2Array, ZZ2Vector, ZZ2Matrix, addcol!, swapcols!,
    rref!, rref, rank, rank!, det, det!, inv, inv!, dot, identity_matrix

import Base: show, ==, +, -, *, /, ^, inv, literal_pow,
    zero, one, iszero, isone, iseven, isodd, convert, rand, promote_rule,
    size, zeros, ones, getindex, setindex!, copy, inv, Bool

using Base: @propagate_inbounds

using Random: AbstractRNG, SamplerType

import LinearAlgebra: dot, det

# otherwise "throw" makes code slower
@noinline throw(e) = Core.throw(e)



#
# ZZ2
#

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

+(a::ZZ2) = a
+(a::ZZ2, b::ZZ2) = ZZ2(xor(a.m, b.m))
-(a::ZZ2) = a
-(a::ZZ2, b::ZZ2) = a + b
*(a::ZZ2, b::ZZ2) = ZZ2(a.m & b.m)
/(a::ZZ2, b::ZZ2) = iszero(b) ? error("division by zero") : a
inv(a::ZZ2) = one(ZZ2)/a
^(a::ZZ2, n::Integer) = n >= 0 ? a : inv(a)

function literal_pow(::Type{^}, a::ZZ2, ::Val{N}) where N
    if N > 0
        a
    elseif N == 0
        one(ZZ2)
    else
        inv(a)
    end
end

rand(rng::AbstractRNG, ::SamplerType{ZZ2}) = ZZ2(rand(rng, Bool))

promote_rule(::Type{<:Integer}, ::Type{ZZ2}) = ZZ2
promote_rule(::Type{Bool}, ::Type{ZZ2}) = ZZ2   # necessary to avoid ambiguities



#
# ZZ2Array
#

using Base: OneTo
import Base: similar, fill!

struct ZZ2Array{N} <: AbstractArray{ZZ2,N}
    i1::Int
    data::Array{UInt,N}
    ZZ2Array{N}(i1::Int, data::Array{UInt,N}) where N = new(i1, data)
    # this avoids confusing error messages
end

const ZZ2Vector = ZZ2Array{1}
const ZZ2Matrix = ZZ2Array{2}

function ZZ2Array{N}(::UndefInitializer, ii::NTuple{N,Integer}) where N
    i1 = 4 * ((ii[1] + 1 << 8 - 1) >> 8)
    ZZ2Array{N}(ii[1], Array{UInt, N}(undef, i1, ii[2:end]...))
end

ZZ2Array{N}(::UndefInitializer, ii::Integer...) where N = ZZ2Array{N}(undef, ii)

ZZ2Array(::UndefInitializer, ii::NTuple{N,Integer}) where N = ZZ2Array{N}(undef, ii)

ZZ2Array(::UndefInitializer, ii::Integer...) = ZZ2Array(undef, ii)

function ZZ2Array{N}(a::AbstractArray{T,N}) where {T,N}
    b = ZZ2Array(undef, size(a)...)
    for i in eachindex(a)
        b[i] = a[i]
    end
    b
end

ZZ2Array(a::AbstractArray{T,N}) where {T,N} = ZZ2Array{N}(a)

function similar(::Type{ZZ2Array{N}},
        ii::Tuple{Union{Integer, OneTo}, Vararg{Union{Integer,OneTo},N}}) where N
    ZZ2Array{N}(undef, map(last, ii))
end

similar(a::A, ::Type{ZZ2}, dims::Union{Int,Tuple{Int,Vararg{Int}}} = size(a)) where A <: ZZ2Array = similar(a, dims)
similar(a::A, dim::Integer = length(a)) where A <: ZZ2Vector  = similar(A,  (dims,))
similar(a::A, dims::Tuple = size(a)) where A <: ZZ2Array  = similar(A, dims isa Integer ? (dims,) : dims)

# ZZ2Vector(a::AbstractVector) = ZZ2Array(a)
# ZZ2Matrix(a::AbstractMatrix) = ZZ2Array(a)

function fill!(a::ZZ2Array, c)
    fill!(a.data, iszero(ZZ2(c)) ? UInt(0) : ~UInt(0))
    a
end

zeros(::Type{ZZ2}, ii::Integer...) = fill!(ZZ2Array(undef, ii...), ZZ2(0))
# TODO: add zero_matrix ?

ones(::Type{ZZ2}, ii::Integer...) = fill!(ZZ2Array(undef, ii...), ZZ2(1))

# TODO: could probably be done more efficiently
function identity_matrix(::Type{ZZ2}, i1::Integer, i2::Integer = i1)
    a = zeros(ZZ2, i1, i2)
    for k in 1:min(i1, i2)
        a[k, k] = one(ZZ2)
    end
    a
end

zero(a::ZZ2Array) = zeros(ZZ2, size(a)...)

size(a::ZZ2Array, d) = d == 1 ? a.i1 : size(a.data, d)

size(a::ZZ2Array) = (a.i1, size(a.data)[2:end]...)

copy(a::ZZ2Array{N}) where N = ZZ2Array{N}(a.i1, copy(a.data))

@inline function _getindex(a::ZZ2Array, ii...)
    @boundscheck checkbounds(a, ii...)
    ii1 = ii[1]-1
    i1 = (ii1 >> 6) + 1
    i0 = ii1 & (1 << 6 - 1)
    @inbounds ZZ2(a.data[i1, ii[2:end]...] >> i0)
end

@propagate_inbounds function getindex(a::ZZ2Array{N}, ii::Vararg{Int,N}) where N
    @boundscheck checkbounds(a, ii...)
    _getindex(a, ii...)
end

@propagate_inbounds getindex(a::ZZ2Vector, i::Int) = _getindex(a, i)

function adjust_index(a::ZZ2Array, i)
    (j2, j1) = divrem(i-1, a.i1)
    j2 * (size(a.data, 1) << 6) + j1 + 1
end

@propagate_inbounds function getindex(a::ZZ2Array, i::Int)
# method for linear indexing
    @boundscheck checkbounds(a, i)
    @inbounds _getindex(a, adjust_index(a, i))
end

function _setindex!(a::ZZ2Array, x, ii...)
    ii1 = ii[1]-1
    i1 = (ii1 >> 6) + 1
    i0 = ii1 & (1 << 6 - 1)
    x = ZZ2(x)
    if iszero(x)
        a.data[i1, ii[2:end]...] &= bitrotate(~UInt(1), i0)
    else
        a.data[i1, ii[2:end]...] |= UInt(1) << i0
    end
    x
end

setindex!(a::ZZ2Array{N}, x, ii::Vararg{Int,N}) where N = _setindex!(a, x, ii...)

setindex!(a::ZZ2Vector, x, i::Int) = _setindex!(a, x, i)

function setindex!(a::ZZ2Array, x, i::Int)
# method for linear indexing
    _setindex!(a, x, adjust_index(a, i))
end

function ==(a::ZZ2Array, b::ZZ2Array)
    s = size(a)
    s == size(b) || return false
    ii1, ii... = s
    l1 = (ii1-1) >> 6 + 1
    i0 = ii1 & (1 << 6 - 1)
    m = UInt(1) << (i0+1) - UInt(1)
    ca = a.data
    cb = b.data
    i1 = size(cb, 1)
    # @show i1 ii
    for k in 0:prod(ii)-1
        for k1 in 1:l1-1
            # @show k1 k
            @inbounds ca[i1*k+k1] == cb[i1*k+k1] || return false
        end
        # @show k
        @inbounds ca[i1*k+l1] & m == cb[i1*k+l1] & m || return false
    end
    return true
end

+(a::ZZ2Array) = copy(a)

function +(a::ZZ2Array{N}, b::ZZ2Array{N}) where N
    ii = size(a)
    jj = size(b)
    ==(ii, jj) || throw(DimensionMismatch("first matrix has dimensions $ii, second matrix has dimensions $jj"))
    data = Array{UInt}(undef, size(a.data))
    data .= a.data .⊻ b.data
    ZZ2Array{N}(a.i1, data)
end

-(a::ZZ2Array) = +(a)

-(a::ZZ2Array, b::ZZ2Array) = a + b

*(c::ZZ2, a::ZZ2Array) = iszero(c) ? zero(a) : copy(a)

function *(a::ZZ2Matrix, b::ZZ2Vector)
    i1, i2 = size(a)
    j1 = size(b, 1)
    i2 == j1 || throw(DimensionMismatch("matrix has dimensions ($i1, $i2), vector has length $j1"))
    c = zeros(ZZ2, i1)
    for k in 1:i2
        if isone(b[k])
            addcol!(c, 1, a, k)
        end
    end
    c
end

function *(a::ZZ2Matrix, b::ZZ2Matrix)
    i1, i2 = size(a)
    j1, j2 = size(b)
    i2 == j1 || throw(DimensionMismatch("first matrix has dimensions ($i1, $i2), second matrix has dimensions ($j1, $j2)"))
    c = zeros(ZZ2, i1, j2)
    @inbounds for k2 in 1:j2, k1 in 1:j1
        if isone(b[k1, k2])
            addcol!(c, k2, a, k1)
        end
    end
    c
end

@inline function dot(v::ZZ2Vector, w::ZZ2Vector)
    @boundscheck size(v) != size(w) && throw(DimensionMismatch("vectors must have same length"))
    s = UInt(0)
    vc = v.data
    wc = w.data
    i1 = size(v, 1)
    n = (i1 + 1 << 6 - 1) >> 6
    # TODO: @turbo doesn't work with xor
    @inbounds for j in 1:n-1
        s ⊻= vc[j] & wc[j]
    end
    mask = ~(~UInt(0) << (i1 & (1 << 6 - 1)))
    # mask = ~UInt(0) >> ((1 << 6) - i1 & (1 << 6 - 1))
    # @show mask
    @inbounds s ⊻= vc[n] & wc[n] & mask
    ZZ2(count_ones(s))
end

@generated function xor!(p1::Ptr{UInt64}, j1::Int, p2::Ptr{UInt64}, j2::Int, ::Val{W}) where W
    ir = """
        %pp1 = inttoptr i64 %0 to <$W x i64>*
        %p1 = getelementptr inbounds <$W x i64>, <$W x i64>* %pp1, i64 %1
        %v1 = load <$W x i64>, <$W x i64>* %p1, align 8
        %pp2 = inttoptr i64 %2 to <$W x i64>*
        %p2 = getelementptr inbounds <$W x i64>, <$W x i64>* %pp2, i64 %3
        %v2 = load <$W x i64>, <$W x i64>* %p2, align 8
        %vr = xor <$W x i64> %v1, %v2
        store <$W x i64> %vr, <$W x i64>* %p1, align 8
        ret void
        """
    quote
        Base.llvmcall($ir, Cvoid, Tuple{Ptr{UInt64},Int64,Ptr{UInt64},Int64}, p1, j1, p2, j2)
    end
end

function addcol!(a::ZZ2Array, k0::Integer, b::ZZ2Array, k1::Integer, range::AbstractUnitRange = axes(a.data, 1))
    c = a.data
    i1 = size(c, 1)
    j1 = ((k0-1)*i1+first(range)-1) >> 2
    j2 = ((k1-1)*i1+first(range)-1) >> 2
    l = (last(range)-first(range)+4) >> 2   # this works because we always go to the end of the column
    # TODO: think about comment!
    for _ in 0:l-1
        xor!(pointer(c), j1, pointer(b.data), j2, Val(4))
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

function _rref!(b::ZZ2Matrix, ::Val{mode}) where mode
# mode == 0: full rref
# mode == 1: no zeros above leading ones
# mode == 2: determinant
# mode == 3: inverse
    full = mode == 0 || mode == 3
    i1, i2 = size(b)
    ii1 = size(b.data, 1)
    k = 1
    if mode == 3
        bi = identity_matrix(ZZ2, i1, i2)
    end
    for j1 in 1:i1
        flag = true
        for j2 in k:i2
            @inbounds if isone(b[j1, j2])
                if j2 != k
                    jj = (j1-1) >> 6 + 1
                    swapcols!(b, j2, k, jj:ii1)
                    mode == 3 && swapcols!(bi, j2, k)
                end
                for l in (full ? 1 : k+1):i2
                    if (!full || l != k) && isone(b[j1, l])
                        jj = (j1-1) >> 6 + 1
                        addcol!(b, l, b, k, jj:ii1)
                        mode == 3 && addcol!(bi, l, bi, k)
                    end
                end
                k += 1
                flag = false
                break
            end
        end
        if mode == 2 && flag
            return zero(ZZ2)
        elseif mode == 3 && flag
            error("matrix not invertible")
        end
    end
    if mode == 2
        return one(ZZ2)
    elseif mode == 3
        return bi
    else
        return (k-1, b)
    end
end

rref!(b::ZZ2Matrix; full = true) = _rref!(b, Val(full ? 0 : 1))

rref(b::ZZ2Matrix; kw...) = rref!(b; kw...)

rank!(b::ZZ2Matrix) = _rref!(b, Val(1))[1]

rank(b::ZZ2Matrix) = rank!(copy(b))

function det!(b::ZZ2Matrix)
    ==(size(b)...) || error("matrix not square")
    _rref!(b, Val(2))
end

det(b::ZZ2Matrix) = det!(copy(b))

function inv!(b::ZZ2Matrix)
    ==(size(b)...) || error("matrix not square")
    _rref!(b, Val(3))
end

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

function randomarray(ii...)
    N = length(ii)
    i1 = 4 * ((ii[1] + 1 << 8 - 1) >> 8)
    ZZ2Array{N}(ii[1], rand(UInt, i1, ii[2:end]...))
end



#
# precompilation
#

for i in 0:3
    precompile(_rref!, (ZZ2Matrix, Val{i}))
end

end
