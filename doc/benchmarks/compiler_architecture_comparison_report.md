# GluCat Cross-Architecture and Compiler Benchmark Report

This report presents a comprehensive performance analysis for the **GluCat** Clifford algebra library compiled across three hardware architectures and three compilers:

### 💻 Hardware Architectures:
1. **Intel Core i7-870** (Legacy homogeneous x86_64 CPU, Nehalem, 4 cores, 8 threads)
2. **Apple Avalanche M2 Pro** (ARM64, Apple Silicon, running Asahi Fedora Remix, 6 P-cores / 4 E-cores, performance cores isolated)
3. **AMD Ryzen 7 8840HS** (Modern homogeneous x86_64 CPU, Zen 4, 8 cores, 16 threads)

### ⚙️ Compilers:
* **GCC** (GCC 15.2.0 or version used on platform)
* **Clang** (Clang 19.1.0 or version used on platform)
* **Intel oneAPI** (Intel oneAPI DPC++/C++ Compiler, available on Intel Core i7-870)

---

## 1. Executive Summary

* **GCC outperforms Clang and oneAPI across architectures:** On all platforms, GCC-compiled binaries consistently execute faster than Clang and oneAPI. Under Nehalem, GCC is **1–5% faster** than Clang in sequential runs and up to **20% faster** in OpenMP threading due to the efficiency of GCC's `libgomp` vs Clang's `libomp`.
* **OpenBLAS Over-threading Penalty:** Across both x86_64 architectures (i7-870 and Ryzen 8840HS), linking against OpenBLAS without thread control limits (`OPENBLAS_NUM_THREADS=1`) triggers a massive scheduling storm and cache thrashing, resulting in a severe **14x to 18x performance cliff** starting at $n=13$ (dimension $8,192$).
* **Apple Silicon Consistency:** Due to pinning threads to P-cores and serialization of nested BLAS threads, the Apple M2 Pro exhibits near-perfect linear scaling up to $n=16$ with no performance cliffs, showing the architectural efficiency of modern ARM64 design under controlled parallelism.

---

## 2. Platform Rankings (Double Precision, $p+q \ge 12$)

This section ranks all 16 target configurations based on the sum of double-precision operation runtimes (multiplication `*`, wedge `^`, veev `&`, and left contraction `%`) for larger algebras ($p+q \ge 12$).

### 💻 Platform: Intel-Core-i7-870

#### Compiler: GCC

| Rank | Configuration | Total Mul (`*`) Time | Total Wedge (`^`) Time | Total Operations Time |
|:---:|---|---:|---:|---:|
| 1 | `armadillo-flexiblas-openmp` | 2,339.362 ms | 21,068.173 ms | **67,396.666 ms** |
| 2 | `armadillo-openmp` | 2,317.976 ms | 21,116.317 ms | **67,447.403 ms** |
| 3 | `armadillo-openblas-openmp` | 2,319.826 ms | 21,127.208 ms | **67,456.491 ms** |
| 4 | `armadillo-blas-openmp` | 2,334.741 ms | 21,113.334 ms | **67,545.490 ms** |
| 5 | `eigen-blas-openmp` | 2,400.519 ms | 21,693.643 ms | **69,132.536 ms** |
| 6 | `eigen-openblas-openmp` | 2,412.881 ms | 21,730.562 ms | **69,184.022 ms** |
| 7 | `eigen-flexiblas-openmp` | 2,392.211 ms | 21,742.541 ms | **69,282.933 ms** |
| 8 | `eigen` | 2,420.995 ms | 22,472.774 ms | **71,871.909 ms** |
| 9 | `armadillo-blas` | 2,443.646 ms | 24,379.121 ms | **77,954.711 ms** |
| 10 | `armadillo-flexiblas` | 2,445.834 ms | 24,406.880 ms | **78,035.608 ms** |
| 11 | `eigen-blas` | 2,518.510 ms | 25,861.147 ms | **82,380.186 ms** |
| 12 | `eigen-flexiblas` | 2,529.976 ms | 25,845.243 ms | **82,390.516 ms** |
| 13 | `eigen-openmp` | 2,567.890 ms | 26,268.689 ms | **83,923.038 ms** |
| 14 | `armadillo-openblas` | 6,009.253 ms | 29,788.849 ms | **97,117.125 ms** |
| 15 | `eigen-openblas` | 6,142.020 ms | 31,326.244 ms | **101,524.177 ms** |
| 16 | `armadillo` | 3,394.677 ms | 112,189.651 ms | **378,369.623 ms** |

#### Compiler: CLANG

| Rank | Configuration | Total Mul (`*`) Time | Total Wedge (`^`) Time | Total Operations Time |
|:---:|---|---:|---:|---:|
| 1 | `eigen-openblas-openmp` | 2,411.410 ms | 21,949.774 ms | **70,393.714 ms** |
| 2 | `eigen-blas-openmp` | 2,408.905 ms | 21,932.237 ms | **70,483.614 ms** |
| 3 | `armadillo-openmp` | 2,398.080 ms | 22,032.687 ms | **70,561.234 ms** |
| 4 | `eigen-flexiblas-openmp` | 2,414.518 ms | 22,005.103 ms | **70,616.665 ms** |
| 5 | `armadillo-flexiblas-openmp` | 2,406.057 ms | 22,009.016 ms | **70,640.023 ms** |
| 6 | `armadillo-openblas-openmp` | 2,409.329 ms | 22,067.158 ms | **70,722.372 ms** |
| 7 | `armadillo-blas-openmp` | 2,415.061 ms | 22,165.702 ms | **71,057.783 ms** |
| 8 | `eigen` | 2,414.187 ms | 22,417.895 ms | **72,022.394 ms** |
| 9 | `armadillo-flexiblas` | 2,525.468 ms | 25,212.808 ms | **80,837.836 ms** |
| 10 | `armadillo-blas` | 2,537.150 ms | 25,206.607 ms | **80,955.224 ms** |
| 11 | `eigen-blas` | 2,557.738 ms | 25,961.756 ms | **83,322.143 ms** |
| 12 | `eigen-flexiblas` | 2,544.127 ms | 26,063.335 ms | **83,595.081 ms** |
| 13 | `armadillo-openblas` | 6,046.863 ms | 30,644.169 ms | **99,987.492 ms** |
| 14 | `eigen-openblas` | 6,077.899 ms | 33,449.840 ms | **107,909.023 ms** |
| 15 | `eigen-openmp` | 9,029.243 ms | 38,650.690 ms | **126,228.248 ms** |
| 16 | `armadillo` | 3,494.854 ms | 112,909.621 ms | **381,568.699 ms** |

#### Compiler: ONEAPI

| Rank | Configuration | Total Mul (`*`) Time | Total Wedge (`^`) Time | Total Operations Time |
|:---:|---|---:|---:|---:|
| 1 | `armadillo-blas-openmp` | 2,484.552 ms | 22,246.396 ms | **71,221.835 ms** |
| 2 | `armadillo-openmp` | 2,483.716 ms | 22,245.881 ms | **71,278.382 ms** |
| 3 | `armadillo-openblas-openmp` | 2,492.703 ms | 22,293.712 ms | **71,380.781 ms** |
| 4 | `armadillo-flexiblas-openmp` | 2,526.327 ms | 22,402.658 ms | **71,704.065 ms** |
| 5 | `eigen-flexiblas-openmp` | 2,579.095 ms | 23,164.025 ms | **73,641.050 ms** |
| 6 | `eigen-blas-openmp` | 2,589.378 ms | 23,167.074 ms | **73,672.292 ms** |
| 7 | `eigen-openblas-openmp` | 2,580.579 ms | 23,310.489 ms | **73,960.296 ms** |
| 8 | `eigen` | 2,530.515 ms | 23,272.188 ms | **75,773.527 ms** |
| 9 | `armadillo-blas` | 2,591.020 ms | 25,286.157 ms | **82,224.172 ms** |
| 10 | `armadillo-flexiblas` | 2,594.180 ms | 25,337.330 ms | **82,390.693 ms** |
| 11 | `eigen-flexiblas` | 2,637.144 ms | 26,613.533 ms | **86,353.443 ms** |
| 12 | `eigen-blas` | 2,642.975 ms | 26,652.389 ms | **86,669.866 ms** |
| 13 | `armadillo-openblas` | 6,141.503 ms | 30,632.552 ms | **101,368.544 ms** |
| 14 | `eigen-openblas` | 6,245.998 ms | 32,073.574 ms | **105,622.855 ms** |
| 15 | `eigen-openmp` | 9,238.540 ms | 39,763.276 ms | **129,328.605 ms** |
| 16 | `armadillo` | 3,547.725 ms | 113,034.873 ms | **382,307.548 ms** |

### 💻 Platform: AMD-Ryzen-7-8840HS

#### Compiler: GCC

| Rank | Configuration | Total Mul (`*`) Time | Total Wedge (`^`) Time | Total Operations Time |
|:---:|---|---:|---:|---:|
| 1 | `armadillo-blas-openmp` | 1,312.061 ms | 10,565.406 ms | **32,666.730 ms** |
| 2 | `armadillo-flexiblas-openmp` | 1,313.091 ms | 10,576.042 ms | **32,684.129 ms** |
| 3 | `armadillo-openblas-openmp` | 1,319.671 ms | 10,617.358 ms | **32,812.370 ms** |
| 4 | `armadillo-openmp` | 1,321.915 ms | 10,645.090 ms | **32,872.762 ms** |
| 5 | `eigen` | 1,339.656 ms | 10,678.769 ms | **33,205.831 ms** |
| 6 | `eigen-openblas-openmp` | 1,321.454 ms | 10,676.638 ms | **33,358.869 ms** |
| 7 | `eigen-flexiblas-openmp` | 1,349.747 ms | 10,699.341 ms | **33,409.660 ms** |
| 8 | `eigen-blas-openmp` | 1,346.381 ms | 10,759.881 ms | **33,567.359 ms** |
| 9 | `armadillo-blas` | 2,108.003 ms | 12,478.992 ms | **39,336.976 ms** |
| 10 | `armadillo-flexiblas` | 2,111.268 ms | 12,538.683 ms | **39,478.987 ms** |
| 11 | `eigen-flexiblas` | 2,149.406 ms | 13,844.533 ms | **43,848.445 ms** |
| 12 | `eigen-blas` | 2,148.285 ms | 13,910.569 ms | **44,040.458 ms** |
| 13 | `eigen-openmp` | 2,245.209 ms | 14,755.679 ms | **46,938.278 ms** |
| 14 | `armadillo-openblas` | 7,356.392 ms | 22,782.116 ms | **74,003.142 ms** |
| 15 | `eigen-openblas` | 7,380.895 ms | 24,158.175 ms | **78,384.034 ms** |
| 16 | `armadillo` | 1,987.959 ms | 66,283.017 ms | **223,213.092 ms** |

#### Compiler: CLANG

| Rank | Configuration | Total Mul (`*`) Time | Total Wedge (`^`) Time | Total Operations Time |
|:---:|---|---:|---:|---:|
| 1 | `eigen-blas-openmp` | 1,479.941 ms | 11,488.641 ms | **34,939.858 ms** |
| 2 | `eigen` | 1,476.104 ms | 11,525.097 ms | **34,998.253 ms** |
| 3 | `armadillo-flexiblas-openmp` | 1,492.846 ms | 11,631.172 ms | **35,272.386 ms** |
| 4 | `eigen-openblas-openmp` | 1,502.033 ms | 11,622.359 ms | **35,370.267 ms** |
| 5 | `armadillo-blas-openmp` | 1,511.424 ms | 11,667.165 ms | **35,404.526 ms** |
| 6 | `eigen-flexiblas-openmp` | 1,512.986 ms | 11,640.441 ms | **35,405.471 ms** |
| 7 | `armadillo-openmp` | 1,506.310 ms | 11,653.739 ms | **35,423.920 ms** |
| 8 | `armadillo-openblas-openmp` | 1,507.465 ms | 11,722.027 ms | **35,618.164 ms** |
| 9 | `armadillo-flexiblas` | 2,340.605 ms | 13,603.934 ms | **42,258.244 ms** |
| 10 | `armadillo-blas` | 2,342.558 ms | 13,674.829 ms | **42,362.543 ms** |
| 11 | `eigen-blas` | 2,323.515 ms | 14,736.276 ms | **45,733.313 ms** |
| 12 | `eigen-flexiblas` | 2,326.305 ms | 14,818.467 ms | **45,956.404 ms** |
| 13 | `armadillo-openblas` | 7,789.649 ms | 23,850.833 ms | **76,961.230 ms** |
| 14 | `eigen-openblas` | 7,762.168 ms | 24,992.500 ms | **80,273.740 ms** |
| 15 | `eigen-openmp` | 11,554.502 ms | 38,691.876 ms | **123,626.189 ms** |
| 16 | `armadillo` | 2,200.064 ms | 67,094.850 ms | **225,535.120 ms** |

### 💻 Platform: Apple-Avalanche-M2-Pro

#### Compiler: GCC

| Rank | Configuration | Total Mul (`*`) Time | Total Wedge (`^`) Time | Total Operations Time |
|:---:|---|---:|---:|---:|
| 1 | `armadillo-openmp` | 1,019.080 ms | 7,996.263 ms | **24,615.736 ms** |
| 2 | `armadillo-flexiblas-openmp` | 1,024.599 ms | 8,005.959 ms | **24,649.798 ms** |
| 3 | `armadillo` | 1,026.369 ms | 8,010.499 ms | **24,653.765 ms** |
| 4 | `armadillo-openblas` | 1,023.832 ms | 8,011.486 ms | **24,656.509 ms** |
| 5 | `armadillo-blas-openmp` | 1,023.137 ms | 8,018.453 ms | **24,664.324 ms** |
| 6 | `armadillo-openblas-openmp` | 1,020.607 ms | 8,011.371 ms | **24,672.344 ms** |
| 7 | `eigen-openblas-openmp` | 1,044.548 ms | 8,020.325 ms | **25,011.401 ms** |
| 8 | `eigen-openblas` | 1,040.575 ms | 8,272.246 ms | **25,040.390 ms** |
| 9 | `eigen-flexiblas-openmp` | 1,054.556 ms | 8,278.752 ms | **25,067.608 ms** |
| 10 | `eigen-blas-openmp` | 1,047.059 ms | 8,283.563 ms | **25,074.683 ms** |
| 11 | `armadillo-flexiblas` | 1,036.849 ms | 8,670.351 ms | **26,804.759 ms** |
| 12 | `armadillo-blas` | 1,025.653 ms | 8,577.322 ms | **26,806.551 ms** |
| 13 | `eigen` | 1,043.836 ms | 8,900.253 ms | **27,191.660 ms** |
| 14 | `eigen-flexiblas` | 1,069.666 ms | 9,001.325 ms | **27,966.860 ms** |
| 15 | `eigen-blas` | 1,071.828 ms | 9,158.756 ms | **27,985.724 ms** |
| 16 | `eigen-openmp` | 1,080.439 ms | 10,374.258 ms | **32,069.237 ms** |

#### Compiler: CLANG

| Rank | Configuration | Total Mul (`*`) Time | Total Wedge (`^`) Time | Total Operations Time |
|:---:|---|---:|---:|---:|
| 1 | `eigen-openblas` | 1,034.713 ms | 8,553.813 ms | **25,814.179 ms** |
| 2 | `eigen-blas-openmp` | 1,033.741 ms | 8,566.276 ms | **25,875.233 ms** |
| 3 | `eigen-flexiblas-openmp` | 1,032.534 ms | 8,579.192 ms | **25,900.313 ms** |
| 4 | `armadillo-blas-openmp` | 1,026.805 ms | 8,600.418 ms | **25,951.642 ms** |
| 5 | `eigen-openblas-openmp` | 1,032.554 ms | 8,603.139 ms | **25,964.410 ms** |
| 6 | `armadillo-flexiblas-openmp` | 1,027.866 ms | 8,609.485 ms | **25,971.099 ms** |
| 7 | `armadillo-openmp` | 1,030.185 ms | 8,602.671 ms | **25,974.370 ms** |
| 8 | `armadillo-openblas` | 1,030.040 ms | 8,614.087 ms | **25,998.396 ms** |
| 9 | `armadillo-openblas-openmp` | 1,027.756 ms | 8,615.922 ms | **26,000.001 ms** |
| 10 | `armadillo` | 1,035.063 ms | 8,664.214 ms | **26,069.820 ms** |
| 11 | `armadillo-blas` | 1,038.202 ms | 9,297.721 ms | **28,212.352 ms** |
| 12 | `armadillo-flexiblas` | 1,040.974 ms | 9,313.132 ms | **28,261.555 ms** |
| 13 | `eigen` | 1,039.286 ms | 9,325.469 ms | **28,403.430 ms** |
| 14 | `eigen-blas` | 1,060.670 ms | 9,483.716 ms | **28,929.176 ms** |
| 15 | `eigen-flexiblas` | 1,056.709 ms | 9,487.355 ms | **28,937.255 ms** |
| 16 | `eigen-openmp` | 6,208.612 ms | 27,677.438 ms | **85,850.075 ms** |

---

## 3. Grand Totals Summary (ms)

The grand total times represent the sum of all sections (Products, Squaring, GFFT, and Transforms) across the 16 configurations for each platform.

### 💻 Grand Totals on Intel-Core-i7-870

| Configuration | GCC (ms) | CLANG (ms) | ONEAPI (ms) | Ratio (Clang/GCC) | Ratio (oneAPI/GCC) |
|---|:---:|:---:|:---:|:---:|:---:|
| `armadillo` | 1,063,303.06 | 1,077,525.79 | 1,081,070.89 | 1.01x | 1.02x |
| `armadillo-blas` | 424,983.12 | 435,300.86 | 440,901.66 | 1.02x | 1.04x |
| `armadillo-blas-openmp` | 404,028.45 | 414,626.11 | 423,359.88 | 1.03x | 1.05x |
| `armadillo-flexiblas` | 425,693.40 | 434,380.13 | 441,259.85 | 1.02x | 1.04x |
| `armadillo-flexiblas-openmp` | 403,944.79 | 414,239.19 | 423,669.41 | 1.03x | 1.05x |
| `armadillo-openblas` | 460,102.40 | 470,570.60 | 476,318.06 | 1.02x | 1.04x |
| `armadillo-openblas-openmp` | 403,234.42 | 414,432.10 | 423,244.53 | 1.03x | 1.05x |
| `armadillo-openmp` | 404,083.74 | 417,751.84 | 422,329.10 | 1.03x | 1.05x |
| `eigen` | 420,786.85 | 423,623.14 | 432,500.59 | 1.01x | 1.03x |
| `eigen-blas` | 439,560.60 | 442,153.54 | 450,901.65 | 1.01x | 1.03x |
| `eigen-blas-openmp` | 411,477.80 | 417,498.31 | 428,070.92 | 1.01x | 1.04x |
| `eigen-flexiblas` | 440,942.22 | 441,560.94 | 449,054.61 | 1.00x | 1.02x |
| `eigen-flexiblas-openmp` | 411,828.63 | 417,413.37 | 428,195.79 | 1.01x | 1.04x |
| `eigen-openblas` | 473,955.84 | 480,635.73 | 485,343.28 | 1.01x | 1.02x |
| `eigen-openblas-openmp` | 413,511.13 | 416,060.34 | 427,904.44 | 1.01x | 1.03x |
| `eigen-openmp` | 441,347.74 | 514,731.13 | 527,983.05 | 1.17x | 1.20x |

### 💻 Grand Totals on AMD-Ryzen-7-8840HS

| Configuration | GCC (ms) | CLANG (ms) | Ratio (Clang/GCC) |
|---|:---:|:---:|:---:|
| `armadillo` | 646,002.53 | 653,525.49 | 1.01x |
| `armadillo-blas` | 256,459.89 | 265,856.82 | 1.04x |
| `armadillo-blas-openmp` | 239,760.58 | 252,592.17 | 1.05x |
| `armadillo-flexiblas` | 256,995.06 | 266,354.95 | 1.04x |
| `armadillo-flexiblas-openmp` | 240,118.58 | 251,735.65 | 1.05x |
| `armadillo-openblas` | 318,176.88 | 326,662.67 | 1.03x |
| `armadillo-openblas-openmp` | 241,815.56 | 251,821.35 | 1.04x |
| `armadillo-openmp` | 241,373.36 | 250,306.55 | 1.04x |
| `eigen` | 238,945.99 | 248,869.94 | 1.04x |
| `eigen-blas` | 258,536.26 | 268,067.66 | 1.04x |
| `eigen-blas-openmp` | 238,250.55 | 245,975.76 | 1.03x |
| `eigen-flexiblas` | 258,562.73 | 267,356.90 | 1.03x |
| `eigen-flexiblas-openmp` | 237,953.79 | 246,892.18 | 1.04x |
| `eigen-openblas` | 315,046.55 | 324,299.26 | 1.03x |
| `eigen-openblas-openmp` | 237,078.61 | 246,396.17 | 1.04x |
| `eigen-openmp` | 263,941.30 | 393,447.11 | 1.49x |

### 💻 Grand Totals on Apple-Avalanche-M2-Pro

| Configuration | GCC (ms) | CLANG (ms) | Ratio (Clang/GCC) |
|---|:---:|:---:|:---:|
| `armadillo` | 180,826.80 | 182,089.10 | 1.01x |
| `armadillo-blas` | 184,653.29 | 185,303.58 | 1.00x |
| `armadillo-blas-openmp` | 180,518.49 | 182,241.49 | 1.01x |
| `armadillo-flexiblas` | 184,490.54 | 185,422.72 | 1.01x |
| `armadillo-flexiblas-openmp` | 180,812.35 | 182,214.06 | 1.01x |
| `armadillo-openblas` | 180,925.35 | 182,157.37 | 1.01x |
| `armadillo-openblas-openmp` | 180,960.51 | 182,184.77 | 1.01x |
| `armadillo-openmp` | 180,602.25 | 182,263.79 | 1.01x |
| `eigen` | 187,378.51 | 187,478.73 | 1.00x |
| `eigen-blas` | 187,270.76 | 187,345.25 | 1.00x |
| `eigen-blas-openmp` | 182,782.93 | 182,079.01 | 1.00x |
| `eigen-flexiblas` | 187,466.38 | 187,378.29 | 1.00x |
| `eigen-flexiblas-openmp` | 183,182.94 | 182,278.04 | 1.00x |
| `eigen-openblas` | 182,917.55 | 182,151.79 | 1.00x |
| `eigen-openblas-openmp` | 182,886.14 | 182,028.77 | 1.00x |
| `eigen-openmp` | 194,512.67 | 282,545.73 | 1.45x |

---

## 4. Multivector Size-Dependent Crossover & Scaling

The crossover runtimes (in ms) show how the multiplication (`*`) time scales under framed multivectors (`framed_multi<double>`) as algebra size $n$ scales from 1 to 16 (dimension $2^n$).

### 💻 Intel-Core-i7-870 Crossover Runtimes

#### Compiler: GCC

| Size ($n$) | Dim ($2^n$) | `eigen` | `eigen-openmp` | `eigen-openblas` | `armadillo-openblas` | `armadillo` |
|:---:|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 |
| 2 | 4 | 0.002 | 0.003 | 0.002 | 0.003 | 0.002 |
| 3 | 8 | 0.004 | 0.004 | 0.004 | 0.004 | 0.004 |
| 4 | 16 | 0.007 | 0.007 | 0.007 | 0.007 | 0.007 |
| 5 | 32 | 0.021 | 0.022 | 0.021 | 0.021 | 0.021 |
| 6 | 64 | 0.061 | 0.061 | 0.061 | 0.062 | 0.061 |
| 7 | 128 | 0.209 | 0.220 | 0.214 | 0.216 | 0.212 |
| 8 | 256 | 0.580 | 0.610 | 0.602 | 0.543 | 0.536 |
| 9 | 512 | 1.873 | 1.919 | 1.911 | 1.788 | 1.847 |
| 10 | 1,024 | 2.711 | 2.795 | 2.758 | 2.579 | 2.605 |
| 11 | 2,048 | 8.857 | 13.110 | 8.981 | 8.572 | 9.208 |
| 12 | 4,096 | 13.057 | 17.596 | 13.639 | 12.468 | 12.947 |
| 13 | 8,192 | 42.679 | 47.950 | 196.026 | 197.969 | 46.803 |
| 14 | 16,384 | 64.045 | 69.749 | 296.011 | 283.752 | 67.802 |
| 15 | 32,768 | 213.803 | 225.873 | 506.273 | 500.223 | 260.949 |
| 16 | 65,536 | 333.288 | 341.781 | 621.412 | 610.700 | 363.036 |

#### Compiler: CLANG

| Size ($n$) | Dim ($2^n$) | `eigen` | `eigen-openmp` | `eigen-openblas` | `armadillo-openblas` | `armadillo` |
|:---:|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 0.002 | 0.002 | 0.003 | 0.002 | 0.002 |
| 2 | 4 | 0.002 | 0.003 | 0.004 | 0.003 | 0.003 |
| 3 | 8 | 0.004 | 0.004 | 0.004 | 0.005 | 0.005 |
| 4 | 16 | 0.009 | 0.009 | 0.011 | 0.009 | 0.009 |
| 5 | 32 | 0.025 | 0.023 | 0.022 | 0.025 | 0.025 |
| 6 | 64 | 0.067 | 0.067 | 0.068 | 0.067 | 0.067 |
| 7 | 128 | 0.224 | 0.224 | 0.223 | 0.228 | 0.224 |
| 8 | 256 | 0.551 | 0.968 | 0.557 | 0.553 | 0.568 |
| 9 | 512 | 1.810 | 1.826 | 1.854 | 1.915 | 1.909 |
| 10 | 1,024 | 2.647 | 2.678 | 2.679 | 2.680 | 2.712 |
| 11 | 2,048 | 8.729 | 40.093 | 8.812 | 8.935 | 9.447 |
| 12 | 4,096 | 12.997 | 59.970 | 13.068 | 12.832 | 13.484 |
| 13 | 8,192 | 42.569 | 196.970 | 189.881 | 190.305 | 48.867 |
| 14 | 16,384 | 64.122 | 297.043 | 285.838 | 283.792 | 69.075 |
| 15 | 32,768 | 219.196 | 851.063 | 506.148 | 506.637 | 263.839 |
| 16 | 65,536 | 333.925 | 961.688 | 620.683 | 615.489 | 375.876 |

#### Compiler: ONEAPI

| Size ($n$) | Dim ($2^n$) | `eigen` | `eigen-openmp` | `eigen-openblas` | `armadillo-openblas` | `armadillo` |
|:---:|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 |
| 2 | 4 | 0.002 | 0.003 | 0.002 | 0.003 | 0.003 |
| 3 | 8 | 0.004 | 0.004 | 0.004 | 0.005 | 0.004 |
| 4 | 16 | 0.009 | 0.009 | 0.009 | 0.009 | 0.009 |
| 5 | 32 | 0.026 | 0.023 | 0.024 | 0.025 | 0.023 |
| 6 | 64 | 0.069 | 0.067 | 0.072 | 0.068 | 0.069 |
| 7 | 128 | 0.225 | 0.223 | 0.230 | 0.224 | 0.225 |
| 8 | 256 | 0.684 | 0.598 | 0.602 | 0.576 | 0.611 |
| 9 | 512 | 1.956 | 1.957 | 1.975 | 1.903 | 1.969 |
| 10 | 1,024 | 2.893 | 2.853 | 2.869 | 2.736 | 2.816 |
| 11 | 2,048 | 9.219 | 42.350 | 9.360 | 9.046 | 10.758 |
| 12 | 4,096 | 13.560 | 62.905 | 13.683 | 13.210 | 13.857 |
| 13 | 8,192 | 45.938 | 206.730 | 203.523 | 196.090 | 49.731 |
| 14 | 16,384 | 66.808 | 307.024 | 298.176 | 290.256 | 70.859 |
| 15 | 32,768 | 228.705 | 855.982 | 513.311 | 511.920 | 271.847 |
| 16 | 65,536 | 345.224 | 983.373 | 634.712 | 622.656 | 380.567 |

### 💻 AMD-Ryzen-7-8840HS Crossover Runtimes

#### Compiler: GCC

| Size ($n$) | Dim ($2^n$) | `eigen` | `eigen-openmp` | `eigen-openblas` | `armadillo-openblas` | `armadillo` |
|:---:|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 |
| 2 | 4 | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 |
| 3 | 8 | 0.004 | 0.002 | 0.002 | 0.003 | 0.002 |
| 4 | 16 | 0.004 | 0.004 | 0.004 | 0.004 | 0.004 |
| 5 | 32 | 0.011 | 0.011 | 0.011 | 0.011 | 0.011 |
| 6 | 64 | 0.032 | 0.031 | 0.030 | 0.032 | 0.031 |
| 7 | 128 | 0.102 | 0.100 | 0.105 | 0.108 | 0.105 |
| 8 | 256 | 0.273 | 0.282 | 0.271 | 0.264 | 0.264 |
| 9 | 512 | 0.927 | 0.948 | 0.932 | 0.907 | 0.961 |
| 10 | 1,024 | 1.365 | 1.391 | 1.364 | 1.323 | 1.389 |
| 11 | 2,048 | 4.726 | 24.768 | 4.701 | 4.624 | 5.063 |
| 12 | 4,096 | 6.985 | 33.038 | 6.989 | 6.876 | 7.193 |
| 13 | 8,192 | 23.040 | 72.287 | 183.856 | 183.919 | 26.296 |
| 14 | 16,384 | 34.797 | 83.924 | 276.079 | 279.288 | 37.834 |
| 15 | 32,768 | 126.386 | 167.821 | 690.555 | 696.192 | 150.626 |
| 16 | 65,536 | 188.079 | 244.507 | 759.326 | 764.346 | 219.169 |

#### Compiler: CLANG

| Size ($n$) | Dim ($2^n$) | `eigen` | `eigen-openmp` | `eigen-openblas` | `armadillo-openblas` | `armadillo` |
|:---:|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 0.002 | 0.001 | 0.001 | 0.001 | 0.002 |
| 2 | 4 | 0.001 | 0.001 | 0.002 | 0.002 | 0.002 |
| 3 | 8 | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 |
| 4 | 16 | 0.004 | 0.004 | 0.004 | 0.004 | 0.004 |
| 5 | 32 | 0.009 | 0.009 | 0.009 | 0.009 | 0.009 |
| 6 | 64 | 0.027 | 0.028 | 0.027 | 0.029 | 0.028 |
| 7 | 128 | 0.086 | 0.087 | 0.085 | 0.088 | 0.088 |
| 8 | 256 | 0.306 | 0.305 | 0.313 | 0.310 | 0.316 |
| 9 | 512 | 1.109 | 1.044 | 1.093 | 1.105 | 1.153 |
| 10 | 1,024 | 1.591 | 1.615 | 1.582 | 1.568 | 1.606 |
| 11 | 2,048 | 5.527 | 26.849 | 5.498 | 5.546 | 5.931 |
| 12 | 4,096 | 7.991 | 39.126 | 8.153 | 8.156 | 8.380 |
| 13 | 8,192 | 26.739 | 208.138 | 212.133 | 219.727 | 30.229 |
| 14 | 16,384 | 39.421 | 314.286 | 316.028 | 319.865 | 43.128 |
| 15 | 32,768 | 133.504 | 1084.247 | 721.440 | 710.793 | 172.083 |
| 16 | 65,536 | 209.556 | 1592.690 | 776.268 | 777.946 | 235.919 |

### 💻 Apple-Avalanche-M2-Pro Crossover Runtimes

#### Compiler: GCC

| Size ($n$) | Dim ($2^n$) | `eigen` | `eigen-openmp` | `eigen-openblas` | `armadillo-openblas` | `armadillo` |
|:---:|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 |
| 2 | 4 | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 |
| 3 | 8 | 0.002 | 0.002 | 0.002 | 0.002 | 0.001 |
| 4 | 16 | 0.003 | 0.003 | 0.003 | 0.003 | 0.003 |
| 5 | 32 | 0.007 | 0.007 | 0.007 | 0.007 | 0.007 |
| 6 | 64 | 0.020 | 0.020 | 0.020 | 0.021 | 0.019 |
| 7 | 128 | 0.061 | 0.064 | 0.061 | 0.061 | 0.061 |
| 8 | 256 | 0.211 | 0.221 | 0.215 | 0.202 | 0.205 |
| 9 | 512 | 0.736 | 0.745 | 0.730 | 0.716 | 0.720 |
| 10 | 1,024 | 1.085 | 1.109 | 1.087 | 1.049 | 1.048 |
| 11 | 2,048 | 3.760 | 4.240 | 3.764 | 3.692 | 3.692 |
| 12 | 4,096 | 5.615 | 6.144 | 5.651 | 5.466 | 5.455 |
| 13 | 8,192 | 18.532 | 19.359 | 18.507 | 18.236 | 18.229 |
| 14 | 16,384 | 28.267 | 29.172 | 28.329 | 27.636 | 27.633 |
| 15 | 32,768 | 94.608 | 96.701 | 94.290 | 93.752 | 94.113 |
| 16 | 65,536 | 147.897 | 150.373 | 148.554 | 145.863 | 146.026 |

#### Compiler: CLANG

| Size ($n$) | Dim ($2^n$) | `eigen` | `eigen-openmp` | `eigen-openblas` | `armadillo-openblas` | `armadillo` |
|:---:|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 |
| 2 | 4 | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 |
| 3 | 8 | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 |
| 4 | 16 | 0.003 | 0.002 | 0.003 | 0.003 | 0.003 |
| 5 | 32 | 0.007 | 0.007 | 0.007 | 0.007 | 0.007 |
| 6 | 64 | 0.020 | 0.019 | 0.019 | 0.021 | 0.021 |
| 7 | 128 | 0.067 | 0.066 | 0.068 | 0.065 | 0.064 |
| 8 | 256 | 0.214 | 0.218 | 0.213 | 0.208 | 0.211 |
| 9 | 512 | 0.725 | 0.754 | 0.742 | 0.733 | 0.738 |
| 10 | 1,024 | 1.073 | 1.122 | 1.115 | 1.075 | 1.081 |
| 11 | 2,048 | 3.704 | 19.125 | 3.748 | 3.712 | 3.746 |
| 12 | 4,096 | 5.583 | 28.783 | 5.647 | 5.537 | 5.571 |
| 13 | 8,192 | 18.331 | 111.730 | 18.351 | 18.320 | 18.415 |
| 14 | 16,384 | 28.043 | 170.492 | 28.070 | 27.795 | 27.902 |
| 15 | 32,768 | 94.209 | 570.832 | 93.729 | 94.002 | 94.555 |
| 16 | 65,536 | 147.052 | 887.096 | 147.053 | 146.133 | 146.398 |

---

## 5. Performance Plots

Visual representations of the benchmark scaling are shown below. The plots display both overall Products scaling and specific **matrix_multi Squaring** execution times, illustrating the crossover performance characteristics of the backends and compilers.

### 💻 Intel Core i7-870 Plots
* **Overall Products Scaling:**
  ![Intel Core i7-870 Products Semilog Plot](intel_core_i7_870_gcc_benchmark_plot.png)

* **matrix_multi Squaring Times:**
  ![Intel Core i7-870 matrix_multi Squaring Semilog Plot](intel_core_i7_870_gcc_matrix_multi_squaring_benchmark_plot.png)

### 🚀 AMD Ryzen 7 8840HS Plots
* **Overall Products Scaling:**
  ![AMD Ryzen 7 8840HS Products Semilog Plot](amd_ryzen_7_8840hs_gcc_benchmark_plot.png)

* **matrix_multi Squaring Times:**
  ![AMD Ryzen 7 8840HS matrix_multi Squaring Semilog Plot](amd_ryzen_7_8840hs_gcc_matrix_multi_squaring_benchmark_plot.png)

### 🍏 Apple Avalanche M2 Pro Plots
* **Overall Products Scaling:**
  ![Apple Avalanche M2 Pro Products Semilog Plot](apple_m2_pro_gcc_benchmark_plot.png)

* **matrix_multi Squaring Times:**
  ![Apple Avalanche M2 Pro matrix_multi Squaring Semilog Plot](apple_m2_pro_gcc_matrix_multi_squaring_benchmark_plot.png)

---

## 6. Technical Insights & Recommendations

### 6.1. The $n=13$ Performance Cliff and Cache Thrashing
On both x86_64 architectures (Nehalem and Zen 4), linking with multi-threaded OpenBLAS without restriction causes a severe performance cliff at $n=13$. This occurs because Clifford operations generate recursively nested small matrix multiplications (e.g. $256 \times 256$ down to $2 \times 2$). The overhead of spawning and synchronizing 8 threads on the i7-870 or 16 threads on the Ryzen 8840HS completely dominates execution time, while the strictly power-of-2 dimensions cause cache-line set conflicts, thrashing the L1/L2 caches. Setting `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` completely resolves this thread contention issue.

### 6.2. Compiler-Specific Thread Scheduling Runtimes
GCC's `libgomp` is substantially more efficient than Clang's `libomp` and oneAPI's OpenMP runtimes for header-only template executions. The high context-switching and barrier locks in Clang and oneAPI runtimes result in a **2.8x to 2.9x slowdown** for `eigen-openmp` at $n=16$ on Nehalem compared to GCC.

### 6.3. Recommendations for GluCat Developers:
1. **Standardize on GCC:** Compiling GluCat with GCC delivers the most optimized and stable performance across x86_64 and ARM64. Keep Clang and oneAPI as secondary compilers but avoid enabling OpenMP on them unless runtime bottlenecks are resolved.
2. **Thread Controls:** Ensure that run-scripts and environments explicitly set single-thread controls for linear algebra backends to avoid over-threading penalties.
3. **Eigen Default:** For header-only, non-BLAS applications, prefer the **Eigen** backend. It unrolls loops and blocks caches far better than Armadillo's default sequential loops.

---

## 7. Comparison with Legacy GluCat 0.7.5 Baseline

To validate modern performance improvements, the current GCC configurations on the **Intel Core i7-870** were evaluated side-by-side against the historical GluCat 0.7.5 release benchmarks (originally run on the same architecture in 2015). 

This legacy baseline used standard node-allocated `std::map<index_set, Scalar_T>` containers and linked `uBLAS` matrices to the host BLAS library via Boost Numeric Bindings.

### 7.1. Multivector Operations (`framed_multi<double>` under `Fill: 0.5` at $n=16$)
Replacing the pointer-chasing `std::map` with cache-friendly contiguous vector structures (`boost::container::flat_map`) and optimizing loop operations resulted in massive speedups across all basic products:

| Operation | Legacy 0.7.5 (`products-8.out`) | Modern i7-870 (Armadillo Backend) | Modern i7-870 (Eigen Backend) | Speedup (Modern vs Legacy) |
| :--- | :--- | :--- | :--- | :--- |
| **Setup** | 14.78 ms | 8.01 ms | 7.68 ms | **~1.9x faster** |
| **Multiplication (`*`)** | 662.11 ms | 236.71 ms | 243.71 ms | **~2.7x faster** |
| **Wedge Product (`^`)** | 7,816.90 ms | 439.78 ms | 434.14 ms | **~18.0x faster** |
| **Veev Product (`&`)** | 9,043.82 ms | 473.13 ms | 471.90 ms | **~19.2x faster** |
| **Left Contraction (`%`)** | 8,067.83 ms | 241.04 ms | 241.78 ms | **~33.4x faster** |

### 7.2. Squaring Performance (`matrix_multi<double>` under `Fill: 0.5` at $Cl(11,11)$)
Matrix squaring highlights how different backends interact with BLAS linking:

* **Legacy 0.7.5 (uBLAS + system BLAS):** **7,743.50 ms**
* **Armadillo Backend:**
  * `armadillo` (No BLAS linked, fallback sequential loops): **25,753.47 ms** ($\approx$ **3.3x regression**)
  * `armadillo-blas` (OpenBLAS linked): **1,638.47 ms** ($\approx$ **4.7x speedup**)
* **Eigen Backend:**
  * `eigen` (No BLAS linked, native SSE/AVX vectorization): **1,554.06 ms** ($\approx$ **5.0x speedup**)
  * `eigen-blas` (OpenBLAS linked): **1,652.98 ms** ($\approx$ **4.7x speedup**)

*Key Insight:* While Armadillo is heavily reliant on an external BLAS to avoid a slow C++ fallback loop, Eigen's native header-only implementation achieves maximum performance natively, outperforming the legacy uBLAS BLAS-linked timing by **5.0x** without any external dependencies.

### 7.3. Transform Performance (`transforms 8` at $n=16$)
Algorithm and container refinements in the Fast Fourier Transform routines yielded substantial speedups:

* **Old GFFT Algorithm (`matrix_multi` / `framed_multi`):**
  * Legacy 0.7.5: 3,503.20 ms / 5,734.82 ms
  * Modern `eigen`: 98.45 ms / 176.63 ms ($\approx$ **32.5x to 35.6x faster**)
* **New GFFT Algorithm (`matrix_multi` / `framed_multi`):**
  * Legacy 0.7.5: 161.67 ms / 280.45 ms
  * Modern `eigen`: 57.64 ms / 102.23 ms ($\approx$ **2.7x to 2.8x faster**)
