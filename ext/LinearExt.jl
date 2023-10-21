module LinearExt

using Modulo2: ZZ2
using LinearCombinations: LinearCombinations

Base.convert(::Type{ZZ2}, ::LinearCombinations.Sign) = ZZ2(1)

LinearCombinations.is_domain(::Type{ZZ2}) = true

LinearCombinations.has_char2(::Type{ZZ2}) = true

end
