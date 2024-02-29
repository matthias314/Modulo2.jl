Mods: 2.2.4
LinearAlgebraX: 0.2.7
AbstractAlgebra: 0.40.1
Nemo: 0.43.1

### Matrix Multiplication

| size | Modulo2 | Matrix{ZZ2} | AbstractAlgebra | Nemo |
| --- | --- | --- | --- | --- |
| 125 | 931.770 μs | 12.569 ms | 291.501 ms | 1.340 ms |
| 250 | 3.933 ms | 101.226 ms | 2.287 s | 6.821 ms |
| 500 | 16.353 ms | 809.666 ms | 18.969 s | 40.489 ms |
| 1000 | 72.019 ms | 6.562 s | skipped | 260.374 ms |
| 2000 | 338.997 ms | skipped | skipped | 1.909 s |
| 4000 | 1.810 s | skipped | skipped | 14.918 s |
| 8000 | 25.469 s | skipped | skipped | 111.866 s |

### Rank

| size | Modulo2 | LinearAlgebraX | AbstractAlgebra | Nemo |
| --- | --- | --- | --- | --- |
| 125 | 361.280 μs | 73.607 ms | 100.184 ms | 9.690 ms |
| 250 | 1.692 ms | 983.601 ms | 796.266 ms | 41.400 ms |
| 500 | 7.012 ms | 7.683 s | 6.530 s | 188.034 ms |
| 1000 | 32.872 ms | skipped | skipped | 1.586 s |
| 2000 | 160.160 ms | skipped | skipped | 8.824 s |
| 4000 | 905.736 ms | skipped | skipped | 46.453 s |
| 8000 | 12.606 s | skipped | skipped | 218.142 s |

### Determinant

| size | Modulo2 | LinearAlgebraX | AbstractAlgebra | Nemo |
| --- | --- | --- | --- | --- |
| 125 | 382.600 μs | 40.103 ms | 256.972 ms | 9.643 ms |
| 250 | 1.798 ms | 328.040 ms | 2.129 s | 41.668 ms |
| 500 | 7.620 ms | 2.679 s | 17.296 s | 188.500 ms |
| 1000 | 33.899 ms | 21.416 s | skipped | 1.491 s |
| 2000 | 160.628 ms | skipped | skipped | 8.750 s |
| 4000 | 904.472 ms | skipped | skipped | 43.444 s |
| 8000 | 12.870 s | skipped | skipped | 199.009 s |

### Inverse

| size | Modulo2 | LinearAlgebraX | AbstractAlgebra | Nemo |
| --- | --- | --- | --- | --- |
| 125 | 1.249 ms | 40.871 ms | 523.685 ms | 11.649 ms |
| 250 | 5.062 ms | 327.477 ms | 4.521 s | 56.150 ms |
| 500 | 22.127 ms | 2.669 s | 36.821 s | 297.618 ms |
| 1000 | 98.185 ms | 21.390 s | skipped | 1.479 s |
| 2000 | 551.627 ms | skipped | skipped | 7.842 s |
| 4000 | 3.549 s | skipped | skipped | 43.986 s |
| 8000 | 70.131 s | skipped | skipped | 277.870 s |

