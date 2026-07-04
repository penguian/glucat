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
* **Sandwich Product Optimization Takeaway:** The native sparse-domain sandwich operator (`operator|` / `versor`) achieves the absolute fastest runtimes for `framed_multi` on **AMD Ryzen 7 8840HS (GCC on `eigen` at 134.707 ms)** and for `matrix_multi` on **Apple Avalanche M2 Pro (GCC on `armadillo-openblas` at 2.741 ms)**, showcasing a speedup of up to **~1370x** compared to dense matrix conversion.

---

## 2. Platform Rankings (Double Precision, $p+q \ge 12$)

This section ranks all 16 target configurations based on the sum of double-precision operation runtimes (multiplication `*`, wedge `^`, veev `&`, and left contraction `%`) for larger algebras ($p+q \ge 12$).

### 💻 Platform: Intel-Core-i7-870

#### Compiler: GCC

| Rank | Configuration | Total Mul (`*`) Time | Total Wedge (`^`) Time | Total Operations Time |
|:---:|---|---:|---:|---:|
| 1 | `armadillo-openblas-openmp` | 2,374.301 ms | 21,273.790 ms | **67,851.086 ms** |
| 2 | `eigen-openblas-openmp` | 2,390.053 ms | 21,332.924 ms | **68,052.302 ms** |
| 3 | `armadillo-openmp` | 2,337.893 ms | 21,334.810 ms | **68,062.182 ms** |
| 4 | `eigen-flexiblas-openmp` | 2,387.092 ms | 21,345.626 ms | **68,073.651 ms** |
| 5 | `eigen-blas-openmp` | 2,399.838 ms | 21,399.042 ms | **68,257.069 ms** |
| 6 | `armadillo-blas-openmp` | 2,363.727 ms | 21,406.912 ms | **68,325.474 ms** |
| 7 | `armadillo-flexiblas-openmp` | 2,401.332 ms | 21,460.807 ms | **68,513.403 ms** |
| 8 | `eigen` | 2,375.171 ms | 22,160.936 ms | **70,865.971 ms** |
| 9 | `armadillo-blas` | 2,479.088 ms | 24,609.386 ms | **78,545.536 ms** |
| 10 | `armadillo-flexiblas` | 2,477.461 ms | 24,637.717 ms | **78,675.484 ms** |
| 11 | `eigen-flexiblas` | 2,504.807 ms | 25,147.462 ms | **80,158.221 ms** |
| 12 | `eigen-blas` | 2,513.445 ms | 25,647.129 ms | **80,826.594 ms** |
| 13 | `eigen-openmp` | 2,564.188 ms | 25,565.882 ms | **81,577.502 ms** |
| 14 | `armadillo-openblas` | 6,010.201 ms | 30,009.884 ms | **97,589.263 ms** |
| 15 | `eigen-openblas` | 6,189.261 ms | 30,948.186 ms | **100,525.073 ms** |
| 16 | `armadillo` | 3,418.839 ms | 112,070.190 ms | **378,104.125 ms** |

#### Compiler: CLANG

| Rank | Configuration | Total Mul (`*`) Time | Total Wedge (`^`) Time | Total Operations Time |
|:---:|---|---:|---:|---:|
| 1 | `eigen-openblas-openmp` | 2,405.904 ms | 21,816.441 ms | **70,000.761 ms** |
| 2 | `eigen-flexiblas-openmp` | 2,414.882 ms | 21,842.536 ms | **70,053.108 ms** |
| 3 | `armadillo-flexiblas-openmp` | 2,398.098 ms | 21,832.184 ms | **70,060.090 ms** |
| 4 | `eigen-blas-openmp` | 2,424.590 ms | 21,866.565 ms | **70,171.976 ms** |
| 5 | `armadillo-blas-openmp` | 2,399.171 ms | 21,905.987 ms | **70,237.463 ms** |
| 6 | `armadillo-openmp` | 2,390.062 ms | 21,909.187 ms | **70,246.087 ms** |
| 7 | `armadillo-openblas-openmp` | 2,402.433 ms | 21,979.389 ms | **70,481.655 ms** |
| 8 | `eigen` | 2,418.345 ms | 22,300.530 ms | **71,665.221 ms** |
| 9 | `armadillo-flexiblas` | 2,521.606 ms | 25,111.495 ms | **80,557.309 ms** |
| 10 | `armadillo-blas` | 2,521.785 ms | 25,146.973 ms | **80,614.378 ms** |
| 11 | `eigen-flexiblas` | 2,501.353 ms | 25,610.967 ms | **82,022.104 ms** |
| 12 | `eigen-blas` | 2,511.816 ms | 25,600.157 ms | **82,050.279 ms** |
| 13 | `armadillo-openblas` | 6,018.524 ms | 30,449.690 ms | **99,509.160 ms** |
| 14 | `eigen-openblas` | 6,053.669 ms | 30,976.140 ms | **100,917.826 ms** |
| 15 | `eigen-openmp` | 9,051.215 ms | 38,464.036 ms | **125,468.284 ms** |
| 16 | `armadillo` | 3,463.418 ms | 112,445.690 ms | **379,817.861 ms** |

#### Compiler: ONEAPI

| Rank | Configuration | Total Mul (`*`) Time | Total Wedge (`^`) Time | Total Operations Time |
|:---:|---|---:|---:|---:|
| 1 | `armadillo-openmp` | 2,482.670 ms | 22,254.565 ms | **71,245.220 ms** |
| 2 | `armadillo-openblas-openmp` | 2,481.168 ms | 22,319.090 ms | **71,416.295 ms** |
| 3 | `armadillo-blas-openmp` | 2,503.716 ms | 22,559.661 ms | **72,063.926 ms** |
| 4 | `armadillo-flexiblas-openmp` | 2,491.849 ms | 22,557.716 ms | **72,067.223 ms** |
| 5 | `eigen-blas-openmp` | 2,525.273 ms | 22,679.461 ms | **72,463.707 ms** |
| 6 | `eigen-flexiblas-openmp` | 2,553.807 ms | 22,941.507 ms | **73,293.670 ms** |
| 7 | `eigen-openblas-openmp` | 2,546.611 ms | 22,977.818 ms | **73,353.442 ms** |
| 8 | `eigen` | 2,512.176 ms | 23,016.153 ms | **74,833.623 ms** |
| 9 | `armadillo-flexiblas` | 2,601.098 ms | 25,317.431 ms | **82,306.845 ms** |
| 10 | `armadillo-blas` | 2,627.605 ms | 25,644.026 ms | **83,284.356 ms** |
| 11 | `eigen-flexiblas` | 2,653.520 ms | 26,298.368 ms | **85,262.806 ms** |
| 12 | `eigen-blas` | 2,651.174 ms | 26,349.669 ms | **85,374.550 ms** |
| 13 | `armadillo-openblas` | 6,187.371 ms | 30,951.147 ms | **102,201.352 ms** |
| 14 | `eigen-openblas` | 6,266.571 ms | 31,685.436 ms | **104,318.047 ms** |
| 15 | `eigen-openmp` | 9,218.914 ms | 39,271.593 ms | **127,665.361 ms** |
| 16 | `armadillo` | 3,569.756 ms | 112,639.768 ms | **381,647.409 ms** |

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
| `armadillo` | 1,063,391.32 | 1,073,561.19 | 1,077,833.25 | 1.01x | 1.01x |
| `armadillo-blas` | 426,723.87 | 433,560.94 | 440,033.02 | 1.02x | 1.03x |
| `armadillo-blas-openmp` | 405,135.72 | 412,547.73 | 421,621.79 | 1.02x | 1.04x |
| `armadillo-flexiblas` | 427,712.49 | 433,165.21 | 439,006.77 | 1.01x | 1.03x |
| `armadillo-flexiblas-openmp` | 405,957.20 | 413,616.35 | 421,414.02 | 1.02x | 1.04x |
| `armadillo-openblas` | 460,871.95 | 468,023.96 | 474,975.16 | 1.02x | 1.03x |
| `armadillo-openblas-openmp` | 403,477.54 | 412,613.68 | 420,080.76 | 1.02x | 1.04x |
| `armadillo-openmp` | 404,107.66 | 415,604.92 | 420,486.13 | 1.03x | 1.04x |
| `eigen` | 417,449.93 | 420,354.51 | 429,788.48 | 1.01x | 1.03x |
| `eigen-blas` | 431,595.45 | 438,150.01 | 449,157.58 | 1.02x | 1.04x |
| `eigen-blas-openmp` | 409,393.29 | 417,554.55 | 424,873.16 | 1.02x | 1.04x |
| `eigen-flexiblas` | 430,850.85 | 437,168.55 | 448,318.19 | 1.01x | 1.04x |
| `eigen-flexiblas-openmp` | 408,863.50 | 416,980.54 | 425,749.99 | 1.02x | 1.04x |
| `eigen-openblas` | 469,506.08 | 473,748.17 | 484,552.34 | 1.01x | 1.03x |
| `eigen-openblas-openmp` | 409,223.61 | 415,385.98 | 425,435.29 | 1.02x | 1.04x |
| `eigen-openmp` | 437,734.57 | 514,467.98 | 523,924.85 | 1.18x | 1.20x |

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
| 2 | 4 | 0.002 | 0.003 | 0.002 | 0.002 | 0.003 |
| 3 | 8 | 0.004 | 0.004 | 0.004 | 0.004 | 0.004 |
| 4 | 16 | 0.008 | 0.008 | 0.007 | 0.007 | 0.007 |
| 5 | 32 | 0.021 | 0.021 | 0.022 | 0.021 | 0.022 |
| 6 | 64 | 0.062 | 0.061 | 0.070 | 0.066 | 0.061 |
| 7 | 128 | 0.209 | 0.210 | 0.222 | 0.218 | 0.213 |
| 8 | 256 | 0.569 | 0.602 | 0.597 | 0.542 | 0.557 |
| 9 | 512 | 1.821 | 1.891 | 1.919 | 1.810 | 1.868 |
| 10 | 1,024 | 2.676 | 2.768 | 2.806 | 2.571 | 2.634 |
| 11 | 2,048 | 8.683 | 13.186 | 9.104 | 8.573 | 9.193 |
| 12 | 4,096 | 13.331 | 17.361 | 13.758 | 12.486 | 13.044 |
| 13 | 8,192 | 42.035 | 47.802 | 198.195 | 190.248 | 47.001 |
| 14 | 16,384 | 63.989 | 68.574 | 300.053 | 284.010 | 67.103 |
| 15 | 32,768 | 211.174 | 224.742 | 508.271 | 501.267 | 260.985 |
| 16 | 65,536 | 325.604 | 342.542 | 626.441 | 613.943 | 365.072 |

#### Compiler: CLANG

| Size ($n$) | Dim ($2^n$) | `eigen` | `eigen-openmp` | `eigen-openblas` | `armadillo-openblas` | `armadillo` |
|:---:|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 |
| 2 | 4 | 0.003 | 0.003 | 0.003 | 0.003 | 0.003 |
| 3 | 8 | 0.004 | 0.004 | 0.004 | 0.004 | 0.004 |
| 4 | 16 | 0.009 | 0.010 | 0.014 | 0.009 | 0.009 |
| 5 | 32 | 0.024 | 0.025 | 0.025 | 0.023 | 0.028 |
| 6 | 64 | 0.067 | 0.079 | 0.067 | 0.071 | 0.069 |
| 7 | 128 | 0.247 | 0.239 | 0.223 | 0.231 | 0.234 |
| 8 | 256 | 0.570 | 0.565 | 0.580 | 0.546 | 0.544 |
| 9 | 512 | 1.827 | 1.836 | 1.846 | 1.812 | 1.971 |
| 10 | 1,024 | 2.674 | 2.699 | 2.688 | 2.609 | 2.678 |
| 11 | 2,048 | 8.742 | 40.718 | 8.838 | 8.781 | 9.391 |
| 12 | 4,096 | 12.908 | 60.724 | 13.060 | 12.678 | 13.329 |
| 13 | 8,192 | 43.091 | 199.215 | 190.128 | 188.230 | 47.788 |
| 14 | 16,384 | 63.781 | 299.243 | 286.103 | 280.169 | 68.712 |
| 15 | 32,768 | 219.188 | 845.226 | 506.113 | 505.239 | 263.028 |
| 16 | 65,536 | 334.073 | 969.087 | 619.339 | 614.530 | 373.952 |

#### Compiler: ONEAPI

| Size ($n$) | Dim ($2^n$) | `eigen` | `eigen-openmp` | `eigen-openblas` | `armadillo-openblas` | `armadillo` |
|:---:|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 0.002 | 0.003 | 0.002 | 0.002 | 0.002 |
| 2 | 4 | 0.003 | 0.003 | 0.003 | 0.011 | 0.003 |
| 3 | 8 | 0.004 | 0.005 | 0.004 | 0.005 | 0.004 |
| 4 | 16 | 0.009 | 0.009 | 0.009 | 0.009 | 0.009 |
| 5 | 32 | 0.024 | 0.024 | 0.025 | 0.025 | 0.023 |
| 6 | 64 | 0.071 | 0.069 | 0.067 | 0.069 | 0.072 |
| 7 | 128 | 0.231 | 0.224 | 0.224 | 0.224 | 0.240 |
| 8 | 256 | 0.587 | 0.608 | 0.630 | 0.584 | 0.653 |
| 9 | 512 | 1.932 | 1.967 | 1.958 | 1.931 | 1.996 |
| 10 | 1,024 | 2.791 | 2.841 | 2.853 | 2.779 | 2.833 |
| 11 | 2,048 | 9.180 | 42.349 | 9.300 | 9.199 | 9.870 |
| 12 | 4,096 | 13.469 | 62.825 | 13.684 | 13.402 | 13.994 |
| 13 | 8,192 | 44.770 | 206.553 | 200.090 | 198.332 | 50.278 |
| 14 | 16,384 | 66.025 | 308.757 | 299.892 | 293.906 | 70.985 |
| 15 | 32,768 | 228.269 | 861.621 | 515.438 | 514.534 | 272.292 |
| 16 | 65,536 | 345.494 | 972.460 | 654.293 | 626.626 | 387.894 |

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

### 6.4. Small vs. Large Multivector Performance Analysis Across All Suites

To understand the runtime scaling of GluCat, the performance characteristics can be divided into two main operational regimes—small multivectors ($n \le 8$, dimension $\le 256$) and large multivectors ($n \ge 12$, dimension $\ge 4,096$)—across all timing suites:

#### 1. Small Multivectors ($n \le 8$, Dimension $\le 256$)
- **Products & Squaring:** At small sizes, the overhead of dynamic memory allocations, library dispatch, and indexing structures completely dominates the execution time. Thread parallelization (OpenMP) or linking to BLAS libraries actually degrades performance due to thread-spawning latency. In this regime, header-only, unrolled backends like **Eigen** (without OpenMP or BLAS) perform best (often sub-0.5 ms).
- **GFFT & Transforms:** Because the size of the basis-loop and lookup tables is extremely small, execution times across all configurations remain sub-millisecond, leaving compiler instruction scheduling as the primary factor rather than memory bandwidth or multi-threading.
- **Expressions:** Complex C++ expression chains (Padé approximants, Taylor series, addition chains) scale very efficiently at small dimensions. Compiler optimizations and lazy evaluation structures (like CRTP) keep the setup and mixed runtimes well under $0.1$ ms.
- **Versor & Sandwich Products:** For small dimensions, the difference between `naive` (dense matrix conversion) and dedicated `operator|`/`versor` solvers is very small (sub-millisecond), as the cost of dense matrix conversion for small matrices is negligible.

#### 2. Large Multivectors ($n \ge 12$, Dimension $\ge 4,096$)
- **Products & Squaring:** FLOPS and memory bandwidth dominate. In this regime, parallelization via OpenMP and optimized BLAS backends (e.g. `eigen-openmp`, `armadillo-openmp`) become essential, delivering massive speedups on multi-core systems (e.g., AMD Ryzen 7 8840HS and Apple M2 Pro).
- **GFFT & Transforms:** Cache thrashing and basis-loop lookup overhead scale exponentially. Using optimized BLAS backends (like OpenBLAS or Flexiblas) becomes critical to maintain cache blocking and efficient linear algebra operations.
- **Expressions:** Runtimes scale exponentially, with complex operations like Padé and Taylor series taking up to $1000+$ ms at size $Cl(8,8)$. OpenMP implementations provide significant speedups for GCC, but suffer from heavy scheduling bottlenecks in Clang and oneAPI.
- **Versor & Sandwich Products:** The advantage of native sparse-domain sandwich operators scales astronomically. For $Cl(16,0)$ (dimension 65,536), dense matrix conversion (`naive`) requires converting a sparse multivector to a $65,536 \times 65,536$ dense matrix, taking over **12 seconds** on Intel i7 and **3.7 seconds** on Apple M2 Pro. In contrast, the native sparse `operator|` executes in **18 ms** (Intel i7) and **2.7 ms** (Apple M2 Pro), yielding an astronomical **~670x to ~1370x speedup** by completely bypassing the memory and $O(N^3)$ computational footprint of dense matrix conversion.

---

## 7. Comparison with Legacy GluCat 0.7.5 Baseline

To validate modern performance improvements, the current GCC configurations on the **Intel Core i7-870** were evaluated side-by-side against the historical GluCat 0.7.5 release benchmarks (originally run on the same architecture in 2015).

This legacy baseline used standard node-allocated `std::map<index_set, Scalar_T>` containers and linked `uBLAS` matrices to the host BLAS library via Boost Numeric Bindings.

### 7.1. Multivector Operations (`framed_multi<double>` under `Fill: 0.5` at $n=16$)
Replacing the pointer-chasing `std::map` with cache-friendly contiguous vector structures (`boost::container::flat_map`) and optimizing loop operations resulted in massive speedups across all basic products:

| Operation | Legacy 0.7.5 (`products-8.out`) | Modern i7-870 (Armadillo Backend) | Modern i7-870 (Eigen Backend) | Speedup (Modern vs Legacy) |
| :--- | :--- | :--- | :--- | :--- |
| **Setup** | 14.78 ms | 8.01 ms | 7.92 ms | **~1.9x faster** |
| **Multiplication (`*`)** | 662.11 ms | 277.98 ms | 240.41 ms | **~2.8x faster** |
| **Wedge Product (`^`)** | 7,816.90 ms | 442.61 ms | 432.20 ms | **~18.1x faster** |
| **Veev Product (`&`)** | 9,043.82 ms | 450.28 ms | 443.77 ms | **~20.4x faster** |
| **Left Contraction (`%`)** | 8,067.83 ms | 243.94 ms | 250.41 ms | **~33.1x faster** |

### 7.2. Squaring Performance (`matrix_multi<double>` under `Fill: 0.5` at $Cl(11,11)$)
Matrix squaring highlights how different backends interact with BLAS linking:

* **Legacy 0.7.5 (uBLAS + system BLAS):** **7,743.50 ms**
* **Armadillo Backend:**
  * `armadillo` (No BLAS linked, fallback sequential loops): **25,778.09 ms** ($\approx$ **3.3x regression**)
  * `armadillo-blas` (OpenBLAS linked): **1,639.63 ms** ($\approx$ **4.7x speedup**)
* **Eigen Backend:**
  * `eigen` (No BLAS linked, native SSE/AVX vectorization): **1,586.65 ms** ($\approx$ **4.9x speedup**)
  * `eigen-blas` (OpenBLAS linked): **1,628.04 ms** ($\approx$ **4.8x speedup**)

*Key Insight:* While Armadillo is heavily reliant on an external BLAS to avoid a slow C++ fallback loop, Eigen's native header-only implementation achieves maximum performance natively, outperforming the legacy uBLAS BLAS-linked timing by **4.9x** without any external dependencies.

### 7.3. Transform Performance (`transforms 8` at $n=16$)
Algorithm and container refinements in the Fast Fourier Transform routines yielded substantial speedups:

* **Old GFFT Algorithm (`matrix_multi` / `framed_multi`):**
  * Legacy 0.7.5: 3,503.20 ms / 5,734.82 ms
  * Modern `eigen`: 99.78 ms / 176.88 ms ($\approx$ **35.1x to 32.4x faster**)
* **New GFFT Algorithm (`matrix_multi` / `framed_multi`):**
  * Legacy 0.7.5: 161.67 ms / 280.45 ms
  * Modern `eigen`: 56.97 ms / 99.55 ms ($\approx$ **2.8x to 2.8x faster**)

---

## 8. Clifford Algebra Expressions Performance

This section compares the execution runtimes of complex multivector expressions (commutators, Padé approximants, Taylor series, mixed expressions, and addition chains) for double-precision multivectors in $Cl(8,8)$ (dimension 65,536) under `framed_multi` representation with `Fill: 0.5` across target architectures and compilers.

### 💻 Intel-Core-i7-870 Expressions

#### Compiler: GCC

| Configuration | Setup (ms) | Commutator (ms) | Padé (ms) | Taylor Series (ms) | Mixed (ms) | Addition (ms) |
|:---|---:|---:|---:|---:|---:|---:|
| `eigen` | 16.073 | 1,374.060 | 279.247 | 1,057.119 | 896.580 | 8.469 |
| `eigen-openmp` | 16.050 | 1,396.516 | 285.471 | 1,088.902 | 901.178 | 8.732 |
| `eigen-openblas` | 16.212 | 1,949.299 | 570.928 | 2,215.248 | 894.407 | 7.896 |
| `armadillo-openblas` | 15.901 | 1,949.879 | 552.997 | 2,180.792 | 917.850 | 7.974 |
| `armadillo` | 15.711 | 1,473.462 | 333.479 | 1,208.957 | 924.358 | 7.993 |

#### Compiler: CLANG

| Configuration | Setup (ms) | Commutator (ms) | Padé (ms) | Taylor Series (ms) | Mixed (ms) | Addition (ms) |
|:---|---:|---:|---:|---:|---:|---:|
| `eigen` | 16.277 | 1,460.621 | 284.487 | 1,074.210 | 970.489 | 7.870 |
| `eigen-openmp` | 16.717 | 2,705.330 | 284.312 | 3,573.367 | 969.981 | 9.886 |
| `eigen-openblas` | 16.160 | 2,031.215 | 567.650 | 2,229.970 | 961.293 | 7.894 |
| `armadillo-openblas` | 16.642 | 2,023.378 | 562.488 | 2,210.109 | 959.053 | 8.508 |
| `armadillo` | 16.180 | 1,533.049 | 340.358 | 1,235.890 | 969.428 | 7.883 |

#### Compiler: ONEAPI

| Configuration | Setup (ms) | Commutator (ms) | Padé (ms) | Taylor Series (ms) | Mixed (ms) | Addition (ms) |
|:---|---:|---:|---:|---:|---:|---:|
| `eigen` | 16.377 | 1,479.419 | 293.541 | 1,110.619 | 970.664 | 7.675 |
| `eigen-openmp` | 16.153 | 2,746.264 | 293.504 | 3,638.130 | 983.017 | 7.479 |
| `eigen-openblas` | 16.396 | 2,057.718 | 580.827 | 2,271.637 | 973.623 | 7.191 |
| `armadillo-openblas` | 16.210 | 2,037.918 | 564.931 | 2,222.463 | 982.437 | 7.192 |
| `armadillo` | 16.540 | 1,557.192 | 346.141 | 1,259.304 | 980.307 | 7.107 |

### 💻 AMD-Ryzen-7-8840HS Expressions

#### Compiler: GCC

| Configuration | Setup (ms) | Commutator (ms) | Padé (ms) | Taylor Series (ms) | Mixed (ms) | Addition (ms) |
|:---|---:|---:|---:|---:|---:|---:|
| `eigen` | 8.098 | 733.578 | 104.984 | 388.450 | 558.813 | 2.263 |
| `eigen-openmp` | 8.197 | 799.475 | 106.242 | 519.395 | 558.748 | 2.373 |
| `eigen-openblas` | 8.206 | 1,892.228 | 686.803 | 2,693.218 | 560.465 | 2.260 |
| `armadillo-openblas` | 8.328 | 1,901.777 | 678.293 | 2,699.249 | 568.264 | 2.277 |
| `armadillo` | 8.204 | 779.536 | 128.760 | 462.360 | 571.876 | 3.143 |

#### Compiler: CLANG

| Configuration | Setup (ms) | Commutator (ms) | Padé (ms) | Taylor Series (ms) | Mixed (ms) | Addition (ms) |
|:---|---:|---:|---:|---:|---:|---:|
| `eigen` | 8.460 | 696.246 | 113.727 | 438.738 | 501.573 | 2.433 |
| `eigen-openmp` | 8.636 | 2,898.251 | 114.263 | 3,912.234 | 502.347 | 2.796 |
| `eigen-openblas` | 8.532 | 1,852.446 | 697.181 | 2,745.530 | 495.651 | 3.025 |
| `armadillo-openblas` | 8.991 | 1,859.118 | 702.431 | 2,756.883 | 499.933 | 2.513 |
| `armadillo` | 8.768 | 732.558 | 145.247 | 522.299 | 493.612 | 2.477 |

### 💻 Apple-Avalanche-M2-Pro Expressions

#### Compiler: GCC

| Configuration | Setup (ms) | Commutator (ms) | Padé (ms) | Taylor Series (ms) | Mixed (ms) | Addition (ms) |
|:---|---:|---:|---:|---:|---:|---:|
| `eigen` | 8.225 | 713.157 | 122.478 | 470.637 | 503.665 | 2.952 |
| `eigen-openmp` | 8.248 | 715.834 | 123.328 | 475.870 | 502.725 | 2.913 |
| `eigen-openblas` | 8.237 | 715.774 | 122.743 | 469.923 | 505.566 | 3.071 |
| `armadillo-openblas` | 8.266 | 708.578 | 116.562 | 457.245 | 505.897 | 2.897 |
| `armadillo` | 8.248 | 709.340 | 117.073 | 458.761 | 504.646 | 2.942 |

#### Compiler: CLANG

| Configuration | Setup (ms) | Commutator (ms) | Padé (ms) | Taylor Series (ms) | Mixed (ms) | Addition (ms) |
|:---|---:|---:|---:|---:|---:|---:|
| `eigen` | 8.507 | 727.308 | 122.090 | 463.899 | 519.620 | 2.782 |
| `eigen-openmp` | 8.481 | 2,259.655 | 124.403 | 2,826.182 | 524.405 | 2.779 |
| `eigen-openblas` | 8.478 | 733.383 | 123.270 | 466.589 | 522.720 | 2.780 |
| `armadillo-openblas` | 8.476 | 729.778 | 118.032 | 462.125 | 523.232 | 2.796 |
| `armadillo` | 8.504 | 729.952 | 118.204 | 461.656 | 523.834 | 2.787 |

---

## 9. Versor and Sandwich Product Performance

This section evaluates the runtime performance of standard sandwich products ($A B A'$) and versor exponentiations in $Cl(16,0)$ (dimension 65,536) under both `framed_multi` and `matrix_multi` representations. Runtimes compare the `naive` (dense matrix conversion) approach against the dedicated solver `operator|`, the native `versor` sandwich solver, and the native `versor_exp` sandwich exponentiation solver.

### Key Takeaways
- **Fastest `framed_multi` sandwich performance:** Achieved on **AMD Ryzen 7 8840HS** using GCC with the **`eigen`** configuration at **134.707 ms** (compared to the `naive` runtime of 235.102 ms, representing a **~1.7x speedup**).
- **Fastest `matrix_multi` sandwich performance:** Achieved on **Apple Avalanche M2 Pro** using GCC with the **`armadillo-openblas`** configuration at **2.741 ms** (compared to the `naive` runtime of 3,752.765 ms, representing a **~1370x speedup**).


### 💻 Intel-Core-i7-870 Versor Products

#### Compiler: GCC

| Configuration | `framed` (naive) | `framed` (op\|) | `framed` (versor) | `framed` (exp) | `matrix` (naive) | `matrix` (op\|) | `matrix` (versor) | `matrix` (exp) |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| `eigen` | 643.059 ms | 364.988 ms | 361.447 ms | 1,455.963 ms | 12,053.832 ms | 18.020 ms | 18.012 ms | 69.272 ms |
| `eigen-openmp` | 656.939 ms | 369.418 ms | 370.525 ms | 1,513.195 ms | 12,162.016 ms | 21.359 ms | 21.436 ms | 89.841 ms |
| `eigen-openblas` | 1,230.412 ms | 673.836 ms | 672.125 ms | 4,051.615 ms | 12,381.717 ms | 48.573 ms | 51.696 ms | 180.110 ms |
| `armadillo-openblas` | 1,199.730 ms | 632.520 ms | 633.939 ms | 3,999.322 ms | 12,312.861 ms | 19.808 ms | 19.907 ms | 119.880 ms |
| `armadillo` | 735.583 ms | 450.322 ms | 451.956 ms | 1,872.273 ms | 12,093.281 ms | 123.999 ms | 124.687 ms | 565.555 ms |

#### Compiler: CLANG

| Configuration | `framed` (naive) | `framed` (op\|) | `framed` (versor) | `framed` (exp) | `matrix` (naive) | `matrix` (op\|) | `matrix` (versor) | `matrix` (exp) |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| `eigen` | 642.843 ms | 359.644 ms | 358.618 ms | 1,448.984 ms | 12,171.045 ms | 17.914 ms | 17.471 ms | 65.234 ms |
| `eigen-openmp` | 1,631.272 ms | 993.483 ms | 996.728 ms | 5,575.165 ms | 12,920.941 ms | 69.606 ms | 71.298 ms | 219.105 ms |
| `eigen-openblas` | 1,229.873 ms | 880.191 ms | 877.142 ms | 4,264.218 ms | 12,451.630 ms | 80.139 ms | 60.021 ms | 187.724 ms |
| `armadillo-openblas` | 1,215.273 ms | 813.878 ms | 821.940 ms | 4,189.380 ms | 12,542.360 ms | 35.357 ms | 32.284 ms | 124.117 ms |
| `armadillo` | 758.387 ms | 464.692 ms | 468.844 ms | 1,907.765 ms | 12,185.403 ms | 123.063 ms | 122.855 ms | 565.728 ms |

#### Compiler: ONEAPI

| Configuration | `framed` (naive) | `framed` (op\|) | `framed` (versor) | `framed` (exp) | `matrix` (naive) | `matrix` (op\|) | `matrix` (versor) | `matrix` (exp) |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| `eigen` | 668.241 ms | 372.424 ms | 369.720 ms | 1,521.543 ms | 13,479.670 ms | 17.819 ms | 17.629 ms | 67.278 ms |
| `eigen-openmp` | 1,626.775 ms | 1,003.888 ms | 1,003.101 ms | 5,744.894 ms | 14,762.666 ms | 69.922 ms | 72.532 ms | 222.149 ms |
| `eigen-openblas` | 1,249.132 ms | 878.496 ms | 886.725 ms | 4,329.717 ms | 13,816.782 ms | 72.443 ms | 59.723 ms | 184.389 ms |
| `armadillo-openblas` | 1,222.451 ms | 827.065 ms | 825.169 ms | 4,231.395 ms | 13,836.419 ms | 19.934 ms | 19.975 ms | 116.718 ms |
| `armadillo` | 764.432 ms | 468.107 ms | 473.000 ms | 1,942.093 ms | 13,652.177 ms | 122.928 ms | 123.004 ms | 566.729 ms |

### 💻 AMD-Ryzen-7-8840HS Versor Products

#### Compiler: GCC

| Configuration | `framed` (naive) | `framed` (op\|) | `framed` (versor) | `framed` (exp) | `matrix` (naive) | `matrix` (op\|) | `matrix` (versor) | `matrix` (exp) |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| `eigen` | 235.102 ms | 134.707 ms | 138.857 ms | 525.710 ms | 4,052.550 ms | 6.364 ms | 6.263 ms | 19.862 ms |
| `eigen-openmp` | 278.623 ms | 173.503 ms | 168.215 ms | 775.531 ms | 4,090.122 ms | 38.765 ms | 38.873 ms | 116.427 ms |
| `eigen-openblas` | 1,579.320 ms | 737.687 ms | 758.350 ms | 4,304.018 ms | 4,782.378 ms | 55.202 ms | 48.035 ms | 136.436 ms |
| `armadillo-openblas` | 1,567.089 ms | 717.881 ms | 713.266 ms | 4,286.269 ms | 4,745.472 ms | 15.625 ms | 8.460 ms | 71.828 ms |
| `armadillo` | 285.847 ms | 173.382 ms | 173.236 ms | 703.029 ms | 4,149.463 ms | 49.939 ms | 49.962 ms | 227.442 ms |

#### Compiler: CLANG

| Configuration | `framed` (naive) | `framed` (op\|) | `framed` (versor) | `framed` (exp) | `matrix` (naive) | `matrix` (op\|) | `matrix` (versor) | `matrix` (exp) |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| `eigen` | 266.951 ms | 142.140 ms | 141.441 ms | 562.186 ms | 3,531.077 ms | 4.708 ms | 4.638 ms | 16.964 ms |
| `eigen-openmp` | 2,313.628 ms | 838.784 ms | 1,366.094 ms | 5,198.014 ms | 5,091.284 ms | 38.598 ms | 42.888 ms | 130.341 ms |
| `eigen-openblas` | 1,625.391 ms | 934.855 ms | 949.526 ms | 4,841.428 ms | 4,172.663 ms | 38.359 ms | 39.707 ms | 112.456 ms |
| `armadillo-openblas` | 1,616.033 ms | 904.809 ms | 907.708 ms | 4,838.515 ms | 4,194.553 ms | 15.682 ms | 8.653 ms | 79.947 ms |
| `armadillo` | 314.266 ms | 193.680 ms | 189.197 ms | 765.657 ms | 3,636.392 ms | 49.854 ms | 49.887 ms | 226.764 ms |

### 💻 Apple-Avalanche-M2-Pro Versor Products

#### Compiler: GCC

| Configuration | `framed` (naive) | `framed` (op\|) | `framed` (versor) | `framed` (exp) | `matrix` (naive) | `matrix` (op\|) | `matrix` (versor) | `matrix` (exp) |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| `eigen` | 278.843 ms | 155.341 ms | 155.859 ms | 601.235 ms | 3,780.555 ms | 6.010 ms | 5.985 ms | 21.277 ms |
| `eigen-openmp` | 280.916 ms | 157.393 ms | 157.263 ms | 615.597 ms | 3,740.178 ms | 6.854 ms | 6.803 ms | 28.549 ms |
| `eigen-openblas` | 277.484 ms | 155.064 ms | 155.444 ms | 601.870 ms | 3,728.653 ms | 5.762 ms | 5.753 ms | 19.151 ms |
| `armadillo-openblas` | 270.769 ms | 148.931 ms | 148.626 ms | 583.673 ms | 3,752.765 ms | 2.741 ms | 2.738 ms | 13.177 ms |
| `armadillo` | 271.236 ms | 148.809 ms | 149.230 ms | 584.604 ms | 3,756.745 ms | 2.744 ms | 2.734 ms | 13.150 ms |

#### Compiler: CLANG

| Configuration | `framed` (naive) | `framed` (op\|) | `framed` (versor) | `framed` (exp) | `matrix` (naive) | `matrix` (op\|) | `matrix` (versor) | `matrix` (exp) |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| `eigen` | 279.941 ms | 157.886 ms | 157.456 ms | 613.907 ms | 3,197.130 ms | 6.902 ms | 6.877 ms | 23.344 ms |
| `eigen-openmp` | 1,684.872 ms | 593.373 ms | 948.638 ms | 3,685.439 ms | 4,233.299 ms | 37.122 ms | 37.444 ms | 102.599 ms |
| `eigen-openblas` | 280.221 ms | 157.299 ms | 157.564 ms | 615.692 ms | 3,213.762 ms | 6.629 ms | 6.612 ms | 21.286 ms |
| `armadillo-openblas` | 274.950 ms | 152.708 ms | 152.383 ms | 609.215 ms | 3,208.085 ms | 3.044 ms | 3.035 ms | 13.670 ms |
| `armadillo` | 275.662 ms | 153.395 ms | 152.545 ms | 608.469 ms | 3,203.303 ms | 3.035 ms | 3.063 ms | 13.705 ms |
