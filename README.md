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
| 125 | 931.770 μs | 12.569 ms | 291.501 ms | 1.340 ms |
| 250 | 3.933 ms | 101.226 ms | 2.287 s | 6.821 ms |
| 500 | 16.353 ms | 809.666 ms | 18.969 s | 40.489 ms |
| 1000 | 72.019 ms | 6.562 s | *skipped* | 260.374 ms |
| 2000 | 338.997 ms | *skipped* | *skipped* | 1.909 s |
| 4000 | 1.810 s | *skipped* | *skipped* | 14.918 s |
| 8000 | 25.469 s | *skipped* | *skipped* | 111.866 s |

### Rank

| size | Modulo2 | LinearAlgebraX | AbstractAlgebra | Nemo |
| ---: | ---: | ---: | ---: | ---: |
| 125 | 361.280 μs | 73.607 ms | 100.184 ms | 9.690 ms |
| 250 | 1.692 ms | 983.601 ms | 796.266 ms | 41.400 ms |
| 500 | 7.012 ms | 7.683 s | 6.530 s | 188.034 ms |
| 1000 | 32.872 ms | *skipped* | *skipped* | 1.586 s |
| 2000 | 160.160 ms | *skipped* | *skipped* | 8.824 s |
| 4000 | 905.736 ms | *skipped* | *skipped* | 46.453 s |
| 8000 | 12.606 s | *skipped* | *skipped* | 218.142 s |

### Determinant

| size | Modulo2 | LinearAlgebraX | AbstractAlgebra | Nemo |
| ---: | ---: | ---: | ---: | ---: |
| 125 | 382.600 μs | 40.103 ms | 256.972 ms | 9.643 ms |
| 250 | 1.798 ms | 328.040 ms | 2.129 s | 41.668 ms |
| 500 | 7.620 ms | 2.679 s | 17.296 s | 188.500 ms |
| 1000 | 33.899 ms | 21.416 s | *skipped* | 1.491 s |
| 2000 | 160.628 ms | *skipped* | *skipped* | 8.750 s |
| 4000 | 904.472 ms | *skipped* | *skipped* | 43.444 s |
| 8000 | 12.870 s | *skipped* | *skipped* | 199.009 s |

### Matrix inverse

| size | Modulo2 | LinearAlgebraX | AbstractAlgebra | Nemo |
| ---: | ---: | ---: | ---: | ---: |
| 125 | 1.249 ms | 40.871 ms | 523.685 ms | 11.649 ms |
| 250 | 5.062 ms | 327.477 ms | 4.521 s | 56.150 ms |
| 500 | 22.127 ms | 2.669 s | 36.821 s | 297.618 ms |
| 1000 | 98.185 ms | 21.390 s | *skipped* | 1.479 s |
| 2000 | 551.627 ms | *skipped* | *skipped* | 7.842 s |
| 4000 | 3.549 s | *skipped* | *skipped* | 43.986 s |
| 8000 | 70.131 s | *skipped* | *skipped* | 277.870 s |

Package versions:
Mods v2.2.4,
LinearAlgebraX v0.2.7,
AbstractAlgebra v0.40.1,
Nemo v0.43.1

Computer: Intel Core i3-10110U CPU @ 2.10GHz with 8GB RAM

The benchmark code can be found in the `benchmark` directory.
