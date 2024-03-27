var documenterSearchIndex = {"docs":
[{"location":"","page":"Modulo2.jl","title":"Modulo2.jl","text":"DocTestSetup = quote\n        using Modulo2\n    # for jldoctest outside of docstrings\n    end","category":"page"},{"location":"#Modulo2.jl","page":"Modulo2.jl","title":"Modulo2.jl","text":"","category":"section"},{"location":"","page":"Modulo2.jl","title":"Modulo2.jl","text":"Modulo2","category":"page"},{"location":"#Modulo2","page":"Modulo2.jl","title":"Modulo2","text":"Modulo2\n\nThis package provides types and functions for linear algebra modulo 2. It defines a type ZZ2 for integers mod 2 and ZZ2Array (ZZ2Vector, ZZ2Matrix) for arrays (vectors, matrices) with elements in ZZ2.\n\nSee also ZZ2, ZZ2Array.\n\n\n\n\n\n","category":"module"},{"location":"#Type-ZZ2","page":"Modulo2.jl","title":"Type ZZ2","text":"","category":"section"},{"location":"","page":"Modulo2.jl","title":"Modulo2.jl","text":"ZZ2\ncount_zeros(::ZZ2)\ncount_ones(::ZZ2)","category":"page"},{"location":"#Modulo2.ZZ2","page":"Modulo2.jl","title":"Modulo2.ZZ2","text":"ZZ2 <: Number\n\nA type representing integers modulo 2.\n\nElements can be created via the functions zero and one or from an argument of type Bool or any other Integer type. More generally, any argument accepted by isodd can be converted to ZZ2. In algebraic computations with ZZ2, Integer types are promoted to ZZ2. Conversely, elements of type ZZ2 can be converted to Bool.\n\nSee also Base.zero, Base.one, Base.iszero, Base.isone, Base.isodd.\n\nExamples\n\njulia> ZZ2(1) == one(ZZ2)\ntrue\n\njulia> iszero(ZZ2(4.0))\ntrue\n\njulia> x = ZZ2(1) + 3\n0\n\njulia> typeof(x)\nZZ2\n\njulia> Bool(x), convert(Bool, x)\n(false, false)\n\n\n\n\n\n","category":"type"},{"location":"#Base.count_zeros-Tuple{ZZ2}","page":"Modulo2.jl","title":"Base.count_zeros","text":"count_zeros(a::ZZ2) -> Integer\n\nReturn 1 if a equals ZZ2(0) and 0 otherwise.\n\nSee also count_ones(::ZZ2), count_zeros(::ZZ2Array), count_zeros(::Integer).\n\n\n\n\n\n","category":"method"},{"location":"#Base.count_ones-Tuple{ZZ2}","page":"Modulo2.jl","title":"Base.count_ones","text":"count_ones(a::ZZ2) -> Integer\n\nReturn 1 if a equals ZZ2(1) and 0 otherwise.\n\nSee also count_zeros(::ZZ2), count_ones(::ZZ2Array), count_ones(::Integer).\n\n\n\n\n\n","category":"method"},{"location":"#Type-ZZ2Array","page":"Modulo2.jl","title":"Type ZZ2Array","text":"","category":"section"},{"location":"","page":"Modulo2.jl","title":"Modulo2.jl","text":"ZZ2Array\nresize!\nzeros\nones\ncount_zeros(::ZZ2Array)\ncount_ones(::ZZ2Array)\nzz2vector\nzz2vector!\nModulo2.randomarray\nmul!\nrcef\nrcef!\nrref\nrank\nrank!\ndet\ndet!\ninv\ninv!","category":"page"},{"location":"#Modulo2.ZZ2Array","page":"Modulo2.jl","title":"Modulo2.ZZ2Array","text":"ZZ2Vector <: AbstractVector{ZZ2}\nZZ2Matrix <: AbstractMatrix{ZZ2}\nZZ2Array{N} <: AbstractArray{ZZ2,N}\n\nAn abstract vector / matrix / array type with elements of type ZZ2.\n\nThe internal representation is packed, meaning that each element only uses one bit. However, columns are internally padded to a length that is a multiple of 256.\n\nA ZZ2Array can be created from any AbstractArray whose elements can be converted to ZZ2. One can also leave the elements undefined by using the undef argument.\n\nSee also zeros, ones, zz2vector, resize!, det, inv, rank, rcef, rref.\n\nExamples\n\njulia> ZZ2Matrix([1 2 3; 4 5 6])\n2×3 ZZ2Matrix:\n 1  0  1\n 0  1  0\n\njulia> v = ZZ2Vector(undef, 2); v[1] = true; v[2] = 2.0; v\n2-element ZZ2Vector:\n 1\n 0\n\n\n\n\n\n","category":"type"},{"location":"#Base.resize!","page":"Modulo2.jl","title":"Base.resize!","text":"resize!(a::ZZ2Vector, n::Integer) -> a\n\nChange the length of a to n and return a. In case a is enlarged, the new entries are undefined.\n\n\n\n\n\n","category":"function"},{"location":"#Base.zeros","page":"Modulo2.jl","title":"Base.zeros","text":"zeros(ZZ2, ii::NTuple{N,Integer}) where N\n\nReturn a ZZ2Array of size ii with zero entries.\n\n\n\n\n\n","category":"function"},{"location":"#Base.ones","page":"Modulo2.jl","title":"Base.ones","text":"ones(ZZ2, ii::NTuple{N,Integer}) where N\n\nReturn a ZZ2Array of size ii with entries ZZ2(1).\n\n\n\n\n\n","category":"function"},{"location":"#Base.count_zeros-Tuple{ZZ2Array}","page":"Modulo2.jl","title":"Base.count_zeros","text":"count_zeros(a::ZZ2Array) -> Integer\n\nReturn the number of elements of a equal to ZZ2(0).\n\nSee also count_ones(::ZZ2Array), count_zeros(::ZZ2), count_zeros(::Integer).\n\n\n\n\n\n","category":"method"},{"location":"#Base.count_ones-Tuple{ZZ2Array}","page":"Modulo2.jl","title":"Base.count_ones","text":"count_ones(a::ZZ2Array) -> Integer\n\nReturn the number of elements of a equal to ZZ2(1).\n\nSee also count_zeros(::ZZ2Array), count_ones(::ZZ2), count_ones(::Integer).\n\n\n\n\n\n","category":"method"},{"location":"#Modulo2.zz2vector","page":"Modulo2.jl","title":"Modulo2.zz2vector","text":"zz2vector(n) -> ZZ2Vector\n\nReturn a ZZ2Vector containing the bit representation of n. Here n must be a Base.BitInteger or a bit integer defined by the package BitIntegers.jl.\n\nSee also zz2vector!, Base.BitInteger, BitIntegers.AbstractBitSigned, BitIntegers.AbstractBitUnsigned.\n\nExample\n\njulia> zz2vector(Int8(35))\n8-element ZZ2Vector:\n 1\n 1\n 0\n 0\n 0\n 1\n 0\n 0\n\n\n\n\n\n","category":"function"},{"location":"#Modulo2.zz2vector!","page":"Modulo2.jl","title":"Modulo2.zz2vector!","text":"zz2vector!(a::ZZ2Vector, n) -> a\n\nFill a with the bit representation of n and return a. Here n must be a Base.BitInteger or a bit integer defined by the package BitIntegers.jl. The length of a is set to the bit length of n.\n\nSee also zz2vector, Base.BitInteger, BitIntegers.AbstractBitSigned, BitIntegers.AbstractBitUnsigned.\n\n\n\n\n\n","category":"function"},{"location":"#Modulo2.randomarray","page":"Modulo2.jl","title":"Modulo2.randomarray","text":"Modulo2.randomarray(ii...) -> ZZ2Array\n\nReturn a ZZ2Array of size ii with random entries.\n\n\n\n\n\n","category":"function"},{"location":"#LinearAlgebra.mul!","page":"Modulo2.jl","title":"LinearAlgebra.mul!","text":"mul!(c::ZZ2Vector, a::ZZ2Matrix, b::AbstractVector{<:Number}, α::Number = ZZ2(1), β::Number = ZZ2(0)) -> c\n\nStore the combined matrix-vector multiply-add α a*b + β c in c and return c. With the default values of α and β the product a*b is computed. The elements of b as well as α and β must be convertible to ZZ2.\n\nThis function is re-exported from the module LinearAlgebra.\n\n\n\n\n\nmul!(c::ZZ2Matrix, a::ZZ2Matrix, b::AbstractMatrix{<:Number}, α::Number = ZZ2(1), β::Number = ZZ2(0)) -> c\n\nStore the combined matrix-matrix multiply-add α a*b + β c in c and return c. With the default values of α and β the product a*b is computed. The elements of b as well as α and β must be convertible to ZZ2.\n\nThis function is re-exported from the module LinearAlgebra.\n\n\n\n\n\n","category":"function"},{"location":"#Modulo2.rcef","page":"Modulo2.jl","title":"Modulo2.rcef","text":"rcef(b::ZZ2Matrix; reduced = true) -> ZZ2Matrix\n\nReturn the tuple (r, c) where r is the rank of the matrix b and c a column echelon form of it. If reduced is true, then the reduced column echelon form is computed.\n\nSee also rref, rcef!.\n\nExamples\n\njulia> a = ZZ2Matrix([1 0 0; 1 1 1])\n2×3 ZZ2Matrix:\n 1  0  0\n 1  1  1\n\njulia> rcef(a)\n(2, ZZ2[1 0 0; 0 1 0])\n\njulia> rcef(a; reduced = false)\n(2, ZZ2[1 0 0; 1 1 0])\n\n\n\n\n\n","category":"function"},{"location":"#Modulo2.rcef!","page":"Modulo2.jl","title":"Modulo2.rcef!","text":"rcef!(b::ZZ2Matrix; reduced = true) -> ZZ2Matrix\n\nReturn the tuple (r, c) where r is the rank of the matrix b and c a column echelon form of it. If reduced is true, then the reduced column echelon form is computed. The argument b may be modified during the computation, which avoids the allocation of a new matrix.\n\nSee also rcef.\n\n\n\n\n\n","category":"function"},{"location":"#Modulo2.rref","page":"Modulo2.jl","title":"Modulo2.rref","text":"rref(b::ZZ2Matrix; reduced = true) -> ZZ2Matrix\n\nReturn the tuple (r, c) where r is the rank of the matrix b and c a row echelon form of it. If reduced is true, then the reduced row echelon form is computed.\n\nNote that it is more efficient to compute a (reduced) column echelon form via rcef.\n\nSee also rcef.\n\nExamples\n\njulia> a = ZZ2Matrix([1 1; 0 1; 0 1])\n3×2 ZZ2Matrix:\n 1  1\n 0  1\n 0  1\n\njulia> rref(a)\n(2, ZZ2[1 0; 0 1; 0 0])\n\njulia> rref(a; reduced = false)\n(2, ZZ2[1 1; 0 1; 0 0])\n\n\n\n\n\n","category":"function"},{"location":"#LinearAlgebra.rank","page":"Modulo2.jl","title":"LinearAlgebra.rank","text":"rank(b::ZZ2Matrix) -> Int\n\nReturn the rank of the matrix b.\n\nSee also rank!.\n\n\n\n\n\n","category":"function"},{"location":"#Modulo2.rank!","page":"Modulo2.jl","title":"Modulo2.rank!","text":"rank!(b::ZZ2Matrix) -> Int\n\nReturn the rank of the matrix b. The argument b may be modified during the computation, which avoids the allocation of a new matrix.\n\nSee also rank.\n\n\n\n\n\n","category":"function"},{"location":"#LinearAlgebra.det","page":"Modulo2.jl","title":"LinearAlgebra.det","text":"det(b::ZZ2Matrix) -> ZZ2\n\nReturn the determinant of the matrix b.\n\nSee also det!.\n\n\n\n\n\n","category":"function"},{"location":"#Modulo2.det!","page":"Modulo2.jl","title":"Modulo2.det!","text":"det!(b::ZZ2Matrix) -> ZZ2\n\nReturn the determinant of the matrix b. The argument b may be modified during the computation, which avoids the allocation of a new matrix.\n\nSee also det.\n\n\n\n\n\n","category":"function"},{"location":"#Base.inv","page":"Modulo2.jl","title":"Base.inv","text":"inv(b::ZZ2Matrix) -> ZZ2Matrix\n\nReturn the inverse of the matrix b, which must be invertible.\n\nSee also inv!.\n\n\n\n\n\n","category":"function"},{"location":"#Modulo2.inv!","page":"Modulo2.jl","title":"Modulo2.inv!","text":"inv!(b::ZZ2Matrix) -> ZZ2Matrix\n\nReturn the inverse of the matrix b, which must be invertible. The argument b may be modified during the computation, which avoids the allocation of a new matrix.\n\nSee also inv.\n\n\n\n\n\n","category":"function"},{"location":"#Broadcasting","page":"Modulo2.jl","title":"Broadcasting","text":"","category":"section"},{"location":"","page":"Modulo2.jl","title":"Modulo2.jl","text":"Broadcasting is partially implemented for the type ZZ2Array, namely for addition, subtraction and scalar multiplication as well as for assignments.","category":"page"},{"location":"","page":"Modulo2.jl","title":"Modulo2.jl","text":"julia> a, b = ZZ2Matrix([1 0; 1 1]), ZZ2Matrix([0 1; 1 0])\n(ZZ2[1 0; 1 1], ZZ2[0 1; 1 0])\n\njulia> c = ZZ2(1)\n1\n\njulia> a .+ c .* b\n2×2 ZZ2Matrix:\n 1  1\n 0  1\n\njulia> a .*= c\n2×2 ZZ2Matrix:\n 1  0\n 1  1\n\njulia> a .= c .* b\n2×2 ZZ2Matrix:\n 0  1\n 1  0\n\njulia> a .+= c .* b\n2×2 ZZ2Matrix:\n 0  0\n 0  0","category":"page"},{"location":"#Internal-functions","page":"Modulo2.jl","title":"Internal functions","text":"","category":"section"},{"location":"","page":"Modulo2.jl","title":"Modulo2.jl","text":"Modulo2.zeropad!\nModulo2.gauss!","category":"page"},{"location":"#Modulo2.zeropad!","page":"Modulo2.jl","title":"Modulo2.zeropad!","text":"Modulo2.zeropad!(a::ZZ2Array{N}) where N -> a\n\nSet the padding bits in the underlying array a.data to zero and return a. The visible bits are not modified.\n\nThis is an internal function of the module.\n\n\n\n\n\n","category":"function"},{"location":"#Modulo2.gauss!","page":"Modulo2.jl","title":"Modulo2.gauss!","text":"Modulo2.gauss!(a::ZZ2Matrix, ::Val{mode}) where mode\n\nIf mode == :rcef, return (r, b) where r is the rank of a and b its reduced column echelon form.\nIf mode == :cef, return (r, b) where r is the rank of a and b a column echelon form of it.\nIf mode == :det, return the determinant of a.\nIf mode == :inv, return the inverse of a.\n\nIn all cases, the matrix a may be modified during the computation. This is an internal function of the module.\n\n\n\n\n\n","category":"function"}]
}