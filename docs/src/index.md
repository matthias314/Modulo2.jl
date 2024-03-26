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

## Internal functions

```@docs
Modulo2.zeropad!
Modulo2.gauss!
```
