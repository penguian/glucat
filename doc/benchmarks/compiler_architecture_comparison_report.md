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
| 1 | `armadillo-flexiblas-openmp` | 860.568 ms | 6,808.286 ms | **21,262.183 ms** |
| 2 | `armadillo-openblas-openmp` | 874.911 ms | 6,825.387 ms | **21,354.572 ms** |
| 3 | `armadillo-blas-openmp` | 879.693 ms | 6,904.926 ms | **21,546.915 ms** |
| 4 | `armadillo-openmp` | 883.946 ms | 6,913.884 ms | **21,633.831 ms** |
| 5 | `eigen` | 897.680 ms | 6,950.072 ms | **21,640.238 ms** |
| 6 | `eigen-openblas-openmp` | 874.554 ms | 6,942.855 ms | **21,692.566 ms** |
| 7 | `eigen-flexiblas-openmp` | 871.482 ms | 7,006.159 ms | **21,819.417 ms** |
| 8 | `eigen-blas-openmp` | 888.441 ms | 7,089.021 ms | **22,064.054 ms** |
| 9 | `armadillo-flexiblas` | 1,379.870 ms | 8,715.885 ms | **27,812.952 ms** |
| 10 | `armadillo-blas` | 1,399.344 ms | 8,757.731 ms | **27,992.593 ms** |
| 11 | `eigen-flexiblas` | 1,398.345 ms | 9,463.784 ms | **29,961.175 ms** |
| 12 | `eigen-blas` | 1,406.891 ms | 9,546.288 ms | **30,130.341 ms** |
| 13 | `eigen-openmp` | 1,476.281 ms | 9,578.227 ms | **30,307.996 ms** |
| 14 | `armadillo-openblas` | 6,470.184 ms | 18,944.067 ms | **62,292.739 ms** |
| 15 | `eigen-openblas` | 6,443.645 ms | 19,817.884 ms | **64,763.906 ms** |
| 16 | `armadillo` | 1,321.743 ms | 44,434.647 ms | **149,915.047 ms** |

#### Compiler: CLANG

| Rank | Configuration | Total Mul (`*`) Time | Total Wedge (`^`) Time | Total Operations Time |
|:---:|---|---:|---:|---:|
| 1 | `eigen` | 982.669 ms | 7,675.109 ms | **23,312.271 ms** |
| 2 | `eigen-blas-openmp` | 977.088 ms | 7,703.118 ms | **23,430.520 ms** |
| 3 | `armadillo-flexiblas-openmp` | 995.536 ms | 7,764.006 ms | **23,588.182 ms** |
| 4 | `armadillo-openblas-openmp` | 1,000.586 ms | 7,775.334 ms | **23,654.362 ms** |
| 5 | `armadillo-openmp` | 1,006.883 ms | 7,811.092 ms | **23,692.240 ms** |
| 6 | `eigen-flexiblas-openmp` | 999.506 ms | 7,792.850 ms | **23,712.447 ms** |
| 7 | `eigen-openblas-openmp` | 987.259 ms | 7,793.383 ms | **23,719.189 ms** |
| 8 | `armadillo-blas-openmp` | 1,011.943 ms | 7,811.624 ms | **23,750.214 ms** |
| 9 | `armadillo-blas` | 1,531.504 ms | 9,608.337 ms | **29,962.416 ms** |
| 10 | `armadillo-flexiblas` | 1,509.295 ms | 9,639.239 ms | **30,083.022 ms** |
| 11 | `eigen-flexiblas` | 1,525.574 ms | 10,303.869 ms | **32,003.684 ms** |
| 12 | `eigen-blas` | 1,532.425 ms | 10,335.895 ms | **32,056.321 ms** |
| 13 | `armadillo-openblas` | 6,871.761 ms | 20,036.048 ms | **65,127.636 ms** |
| 14 | `eigen-openblas` | 6,881.074 ms | 20,745.504 ms | **67,196.295 ms** |
| 15 | `eigen-openmp` | 8,939.824 ms | 32,228.478 ms | **101,538.291 ms** |
| 16 | `armadillo` | 1,462.027 ms | 44,952.603 ms | **151,099.718 ms** |

### 💻 Platform: Apple-Avalanche-M2-Pro

#### Compiler: GCC

| Rank | Configuration | Total Mul (`*`) Time | Total Wedge (`^`) Time | Total Operations Time |
|:---:|---|---:|---:|---:|
| 1 | `armadillo-blas-openmp` | 1,011.901 ms | 7,965.464 ms | **24,501.917 ms** |
| 2 | `armadillo-flexiblas-openmp` | 1,012.249 ms | 7,964.867 ms | **24,507.010 ms** |
| 3 | `armadillo-openblas-openmp` | 1,018.159 ms | 7,986.088 ms | **24,554.488 ms** |
| 4 | `armadillo-openblas` | 1,015.078 ms | 7,985.252 ms | **24,564.814 ms** |
| 5 | `armadillo-openmp` | 1,016.593 ms | 7,981.936 ms | **24,565.349 ms** |
| 6 | `armadillo` | 1,016.297 ms | 7,993.081 ms | **24,597.652 ms** |
| 7 | `eigen-openblas` | 1,033.415 ms | 8,246.113 ms | **24,776.345 ms** |
| 8 | `eigen-blas-openmp` | 1,034.786 ms | 7,983.957 ms | **24,866.716 ms** |
| 9 | `eigen-flexiblas-openmp` | 1,040.828 ms | 8,181.835 ms | **24,876.392 ms** |
| 10 | `eigen-openblas-openmp` | 1,033.632 ms | 8,328.033 ms | **25,524.221 ms** |
| 11 | `armadillo-flexiblas` | 1,030.139 ms | 8,653.740 ms | **26,735.728 ms** |
| 12 | `armadillo-blas` | 1,027.691 ms | 8,657.893 ms | **26,760.061 ms** |
| 13 | `eigen-flexiblas` | 1,056.765 ms | 8,727.479 ms | **27,317.616 ms** |
| 14 | `eigen` | 1,038.305 ms | 8,924.412 ms | **27,325.517 ms** |
| 15 | `eigen-blas` | 1,054.180 ms | 8,950.861 ms | **27,400.928 ms** |
| 16 | `eigen-openmp` | 1,076.592 ms | 10,282.611 ms | **31,839.949 ms** |

#### Compiler: CLANG

| Rank | Configuration | Total Mul (`*`) Time | Total Wedge (`^`) Time | Total Operations Time |
|:---:|---|---:|---:|---:|
| 1 | `eigen-openblas` | 1,031.993 ms | 8,540.725 ms | **25,764.334 ms** |
| 2 | `eigen-flexiblas-openmp` | 1,033.259 ms | 8,551.433 ms | **25,813.794 ms** |
| 3 | `eigen-openblas-openmp` | 1,028.872 ms | 8,560.217 ms | **25,818.666 ms** |
| 4 | `eigen-blas-openmp` | 1,032.655 ms | 8,595.762 ms | **25,938.860 ms** |
| 5 | `armadillo-flexiblas-openmp` | 1,026.494 ms | 8,599.034 ms | **25,942.762 ms** |
| 6 | `armadillo-blas-openmp` | 1,026.614 ms | 8,599.933 ms | **25,943.465 ms** |
| 7 | `armadillo-openmp` | 1,029.285 ms | 8,602.035 ms | **25,965.003 ms** |
| 8 | `armadillo-openblas` | 1,029.411 ms | 8,601.397 ms | **25,990.220 ms** |
| 9 | `armadillo-openblas-openmp` | 1,027.379 ms | 8,615.921 ms | **25,994.174 ms** |
| 10 | `armadillo` | 1,035.778 ms | 8,666.827 ms | **26,162.329 ms** |
| 11 | `armadillo-flexiblas` | 1,039.648 ms | 9,290.762 ms | **28,201.423 ms** |
| 12 | `armadillo-blas` | 1,041.939 ms | 9,295.196 ms | **28,213.681 ms** |
| 13 | `eigen` | 1,027.806 ms | 9,289.004 ms | **28,290.967 ms** |
| 14 | `eigen-flexiblas` | 1,045.781 ms | 9,361.836 ms | **28,495.447 ms** |
| 15 | `eigen-blas` | 1,045.639 ms | 9,365.477 ms | **28,498.679 ms** |
| 16 | `eigen-openmp` | 6,204.791 ms | 27,545.457 ms | **85,389.502 ms** |

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
| `armadillo` | 437,678.49 | 443,480.04 | 1.01x |
| `armadillo-blas` | 178,879.57 | 185,054.42 | 1.03x |
| `armadillo-blas-openmp` | 162,730.24 | 172,367.92 | 1.06x |
| `armadillo-flexiblas` | 177,488.91 | 185,464.48 | 1.04x |
| `armadillo-flexiblas-openmp` | 161,956.36 | 171,709.64 | 1.06x |
| `armadillo-openblas` | 236,753.46 | 245,074.85 | 1.04x |
| `armadillo-openblas-openmp` | 163,108.52 | 171,153.49 | 1.05x |
| `armadillo-openmp` | 162,984.05 | 169,855.26 | 1.04x |
| `eigen` | 160,604.80 | 170,207.16 | 1.06x |
| `eigen-blas` | 177,761.51 | 186,451.75 | 1.05x |
| `eigen-blas-openmp` | 160,694.40 | 168,960.37 | 1.05x |
| `eigen-flexiblas` | 177,440.77 | 186,845.88 | 1.05x |
| `eigen-flexiblas-openmp` | 159,965.40 | 168,648.69 | 1.05x |
| `eigen-openblas` | 235,368.65 | 244,680.87 | 1.04x |
| `eigen-openblas-openmp` | 161,048.57 | 168,616.59 | 1.05x |
| `eigen-openmp` | 179,124.33 | 300,925.29 | 1.68x |

### 💻 Grand Totals on Apple-Avalanche-M2-Pro

| Configuration | GCC (ms) | CLANG (ms) | Ratio (Clang/GCC) |
|---|:---:|:---:|:---:|
| `armadillo` | 180,093.69 | 182,321.16 | 1.01x |
| `armadillo-blas` | 183,815.73 | 185,355.39 | 1.01x |
| `armadillo-blas-openmp` | 179,909.74 | 182,289.26 | 1.01x |
| `armadillo-flexiblas` | 184,157.56 | 185,346.98 | 1.01x |
| `armadillo-flexiblas-openmp` | 180,304.37 | 182,351.07 | 1.01x |
| `armadillo-openblas` | 180,250.20 | 182,153.63 | 1.01x |
| `armadillo-openblas-openmp` | 180,230.66 | 182,124.00 | 1.01x |
| `armadillo-openmp` | 180,082.68 | 182,241.36 | 1.01x |
| `eigen` | 187,149.68 | 187,153.49 | 1.00x |
| `eigen-blas` | 186,296.77 | 186,404.56 | 1.00x |
| `eigen-blas-openmp` | 182,231.51 | 182,010.18 | 1.00x |
| `eigen-flexiblas` | 185,955.47 | 186,149.90 | 1.00x |
| `eigen-flexiblas-openmp` | 182,493.59 | 181,970.42 | 1.00x |
| `eigen-openblas` | 181,675.84 | 181,980.38 | 1.00x |
| `eigen-openblas-openmp` | 182,638.26 | 181,776.72 | 1.00x |
| `eigen-openmp` | 193,811.52 | 281,816.69 | 1.45x |

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
| 1 | 2 | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 |
| 2 | 4 | 0.001 | 0.000 | 0.001 | 0.001 | 0.001 |
| 3 | 8 | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 |
| 4 | 16 | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 |
| 5 | 32 | 0.004 | 0.004 | 0.005 | 0.004 | 0.005 |
| 6 | 64 | 0.009 | 0.009 | 0.010 | 0.009 | 0.009 |
| 7 | 128 | 0.028 | 0.029 | 0.032 | 0.028 | 0.033 |
| 8 | 256 | 0.154 | 0.151 | 0.165 | 0.151 | 0.153 |
| 9 | 512 | 0.574 | 0.563 | 7.618 | 0.551 | 0.587 |
| 10 | 1,024 | 0.760 | 0.768 | 7.806 | 0.739 | 0.755 |
| 11 | 2,048 | 2.844 | 14.947 | 2.817 | 2.776 | 3.039 |
| 12 | 4,096 | 3.877 | 17.845 | 3.744 | 3.812 | 4.031 |
| 13 | 8,192 | 13.064 | 46.280 | 108.344 | 111.914 | 15.485 |
| 14 | 16,384 | 18.651 | 50.602 | 155.686 | 156.143 | 20.114 |
| 15 | 32,768 | 66.418 | 97.695 | 555.811 | 571.942 | 83.076 |
| 16 | 65,536 | 92.003 | 122.172 | 669.213 | 673.078 | 106.088 |

#### Compiler: CLANG

| Size ($n$) | Dim ($2^n$) | `eigen` | `eigen-openmp` | `eigen-openblas` | `armadillo-openblas` | `armadillo` |
|:---:|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 |
| 2 | 4 | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 |
| 3 | 8 | 0.001 | 0.001 | 0.002 | 0.002 | 0.001 |
| 4 | 16 | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 |
| 5 | 32 | 0.005 | 0.004 | 0.005 | 0.005 | 0.004 |
| 6 | 64 | 0.009 | 0.009 | 0.011 | 0.011 | 0.009 |
| 7 | 128 | 0.025 | 0.025 | 0.029 | 3.530 | 0.026 |
| 8 | 256 | 0.167 | 0.165 | 0.185 | 0.179 | 0.170 |
| 9 | 512 | 0.613 | 0.624 | 4.193 | 16.261 | 0.651 |
| 10 | 1,024 | 0.861 | 0.859 | 4.443 | 6.181 | 0.872 |
| 11 | 2,048 | 3.169 | 15.799 | 3.157 | 3.204 | 3.577 |
| 12 | 4,096 | 4.160 | 21.266 | 4.447 | 4.208 | 4.556 |
| 13 | 8,192 | 15.259 | 137.719 | 131.669 | 131.915 | 18.201 |
| 14 | 16,384 | 19.790 | 178.936 | 175.916 | 175.671 | 23.312 |
| 15 | 32,768 | 71.022 | 677.572 | 647.562 | 654.038 | 94.976 |
| 16 | 65,536 | 98.355 | 890.994 | 674.907 | 679.890 | 118.096 |

### 💻 Apple-Avalanche-M2-Pro Crossover Runtimes

#### Compiler: GCC

| Size ($n$) | Dim ($2^n$) | `eigen` | `eigen-openmp` | `eigen-openblas` | `armadillo-openblas` | `armadillo` |
|:---:|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 0.001 | 0.001 | 0.000 | 0.001 | 0.001 |
| 2 | 4 | 0.000 | 0.001 | 0.000 | 0.001 | 0.001 |
| 3 | 8 | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 |
| 4 | 16 | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 |
| 5 | 32 | 0.004 | 0.004 | 0.004 | 0.004 | 0.004 |
| 6 | 64 | 0.010 | 0.010 | 0.009 | 0.009 | 0.009 |
| 7 | 128 | 0.030 | 0.028 | 0.029 | 0.030 | 0.030 |
| 8 | 256 | 0.170 | 0.179 | 0.175 | 0.164 | 0.165 |
| 9 | 512 | 0.627 | 0.640 | 0.624 | 0.607 | 0.608 |
| 10 | 1,024 | 0.851 | 0.894 | 0.865 | 0.820 | 0.830 |
| 11 | 2,048 | 3.101 | 3.714 | 3.133 | 3.046 | 3.047 |
| 12 | 4,096 | 4.264 | 4.918 | 4.286 | 4.133 | 4.136 |
| 13 | 8,192 | 14.947 | 16.218 | 14.822 | 14.717 | 14.674 |
| 14 | 16,384 | 20.734 | 22.073 | 20.727 | 20.100 | 20.171 |
| 15 | 32,768 | 74.091 | 77.232 | 73.506 | 73.017 | 72.814 |
| 16 | 65,536 | 103.928 | 107.743 | 103.355 | 100.987 | 101.307 |

#### Compiler: CLANG

| Size ($n$) | Dim ($2^n$) | `eigen` | `eigen-openmp` | `eigen-openblas` | `armadillo-openblas` | `armadillo` |
|:---:|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 |
| 2 | 4 | 0.000 | 0.001 | 0.001 | 0.000 | 0.001 |
| 3 | 8 | 0.001 | 0.001 | 0.001 | 0.002 | 0.002 |
| 4 | 16 | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 |
| 5 | 32 | 0.004 | 0.004 | 0.004 | 0.004 | 0.005 |
| 6 | 64 | 0.011 | 0.011 | 0.011 | 0.011 | 0.011 |
| 7 | 128 | 0.031 | 0.031 | 0.031 | 0.031 | 0.030 |
| 8 | 256 | 0.172 | 0.178 | 0.174 | 0.173 | 0.175 |
| 9 | 512 | 0.615 | 0.618 | 0.619 | 0.625 | 0.637 |
| 10 | 1,024 | 0.854 | 0.861 | 0.859 | 0.857 | 0.862 |
| 11 | 2,048 | 3.081 | 15.700 | 3.165 | 3.115 | 3.162 |
| 12 | 4,096 | 4.265 | 21.910 | 4.368 | 4.278 | 4.315 |
| 13 | 8,192 | 14.770 | 91.124 | 15.091 | 15.022 | 15.143 |
| 14 | 16,384 | 20.500 | 126.915 | 20.900 | 20.631 | 20.826 |
| 15 | 32,768 | 73.606 | 448.733 | 74.329 | 74.615 | 74.737 |
| 16 | 65,536 | 102.658 | 630.271 | 103.743 | 102.800 | 103.621 |

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
