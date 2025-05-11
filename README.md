# Modulo2.jl

This package provides types and functions for linear algebra modulo 2.
It defines a type `ZZ2` for integers mod 2 and `ZZ2Array` (`ZZ2Vector`, `ZZ2Matrix`)
for arrays (vectors, matrices) with elements in `ZZ2`.

See the [documentation](https://matthias314.github.io/Modulo2.jl/stable/) for details.

## Benchmarks

Below we compare Modulo2.jl to the Julia packages
[LinearAlgebraX.jl](https://github.com/scheinerman/LinearAlgebraX.jl)
(with [Mods.jl](https://github.com/scheinerman/Mods.jl) for integers mod 2),
[AbstractAlgebra.jl](https://github.com/Nemocas/AbstractAlgebra.jl)
and
[Nemo.jl](https://github.com/Nemocas/Nemo.jl).
Note that unlike Modulo2.jl these packages can deal with more general coefficients than integers mod 2.
Also keep in mind that the precise results depend on the specific processor and memory.

The timings are for square matrices of the size indicated.

### Matrix multiplication

`Matrix{ZZ2}` refers to Julia's standard matrix multiplication with coefficients in `ZZ2`.

| size | Modulo2 | Matrix{ZZ2} | AbstractAlgebra | Nemo |
| ---: | ---: | ---: | ---: | ---: |
| 125 | 96.809 μs | 868.790 μs | 30.176 ms | 137.154 μs |
| 250 | 417.656 μs | 6.820 ms | 248.557 ms | 685.672 μs |
| 500 | 1.918 ms | 52.957 ms | 1.920 s | 4.383 ms |
| 1000 | 8.780 ms | 415.689 ms | *skipped* | 25.891 ms |
| 2000 | 43.421 ms | *skipped* | *skipped* | 240.183 ms |
| 4000 | 264.777 ms | *skipped* | *skipped* | 1.391 s |
| 8000 | 3.340 s | *skipped* | *skipped* | 10.282 s |

### Rank

| size | Modulo2 | LinearAlgebraX | AbstractAlgebra | Nemo |
| ---: | ---: | ---: | ---: | ---: |
| 125 | 49.495 μs | 10.139 ms | 12.187 ms | 1.162 ms |
| 250 | 227.735 μs | 125.813 ms | 86.832 ms | 4.902 ms |
| 500 | 990.822 μs | 944.049 ms | 689.586 ms | 22.287 ms |
| 1000 | 4.363 ms | *skipped* | *skipped* | 183.620 ms |
| 2000 | 20.575 ms | *skipped* | *skipped* | 821.456 ms |
| 4000 | 114.753 ms | *skipped* | *skipped* | 2.917 s |
| 8000 | 1.475 s | *skipped* | *skipped* | 14.506 s |

### Determinant

| size | Modulo2 | LinearAlgebraX | AbstractAlgebra | Nemo |
| ---: | ---: | ---: | ---: | ---: |
| 125 | 44.803 μs | 4.838 ms | 29.456 ms | 1.538 ms |
| 250 | 214.386 μs | 38.796 ms | 230.214 ms | 4.907 ms |
| 500 | 913.450 μs | 310.935 ms | 1.860 s | 28.394 ms |
| 1000 | 4.172 ms | 2.519 s | *skipped* | 250.680 ms |
| 2000 | 20.911 ms | *skipped* | *skipped* | 969.185 ms |
| 4000 | 116.159 ms | *skipped* | *skipped* | 3.172 s |
| 8000 | 1.557 s | *skipped* | *skipped* | 14.594 s |

### Matrix inverse

| size | Modulo2 | LinearAlgebraX | AbstractAlgebra | Nemo |
| ---: | ---: | ---: | ---: | ---: |
| 125 | 140.391 μs | 4.838 ms | 54.038 ms | 1.484 ms |
| 250 | 575.657 μs | 38.790 ms | 443.860 ms | 7.113 ms |
| 500 | 2.594 ms | 311.347 ms | 3.717 s | 35.020 ms |
| 1000 | 12.623 ms | 2.568 s | *skipped* | 310.629 ms |
| 2000 | 67.567 ms | *skipped* | *skipped* | 841.042 ms |
| 4000 | 509.087 ms | *skipped* | *skipped* | 5.103 s |
| 8000 | 7.770 s | *skipped* | *skipped* | 27.526 s |

Package versions:
Modulo2: 0.2.3,
Mods: 2.2.6,
LinearAlgebraX: 0.2.10,
AbstractAlgebra: 0.44.13,
Nemo: 0.49.2;
Julia 1.11.5

Computer: Intel Core i3-10110U CPU @ 2.10GHz with 8GB RAM

The benchmark code can be found in the `benchmark` directory.
