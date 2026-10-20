# AMD-lite vs SuiteSparse AMD

Comparison of the lightweight implementation in this directory with the
original SuiteSparse AMD package.  The original AMD source tree is not bundled
with AMD-lite; the benchmark script takes its path as an argument.

## Methodology

- Compiler: GCC 13.3.0, x86-64, `-O3 -std=c99 -DNDEBUG`.
- Original AMD was built locally from a SuiteSparse AMD source tree supplied
  to `bench/run_compare.sh`, using a minimal `SuiteSparse_config` stub
  (`SuiteSparse_malloc` / `SuiteSparse_free` mapped to `malloc` / `free`).
  This isolates the AMD code itself from the separate SuiteSparse
  configuration library.
- Both libraries were given identical CSR input.  For symmetric matrices,
  CSR and CSC describe the same pattern; `A+A'` is the same either way.
- Benchmark times are the best of 3–7 runs on a pinned CPU
  (`taskset -c 0`).  Wall times fluctuate a little with hardware and system
  load, but the ratios are stable.
- All reported `lnz` values are the libraries' own predicted nonzeros in `L`
  (`info.lnz` for AMD-lite, `Info[AMD_LNZ]` for the original).

## Code size

Library source only, excluding tests, demo, and the benchmark harness:

| | lines | shared `.so` text | shared `.so` file | static `.a` file |
|---|---:|---:|---:|---:|
| **AMD-lite** | **1,632** | **26.1 KB** | **36.7 KB** | **28.6 KB** |
| SuiteSparse AMD, full | 4,194 | 35.3 KB | 50.3 KB | 78.7 KB |
| SuiteSparse AMD, core* | 3,771 | 29.8 KB | 41.6 KB | 56.0 KB |

\* core excludes `amd_control`, `amd_info`, `amd_dump`, and `amd_version`.

- AMD-lite is about **39%** of the full original source size
  (1,632 / 4,194 lines).
- The shared library text is about **74%** of the full original
  (26.1 / 35.3 KB), or about **88%** of the core-only original.
- Both implement int32 and int64 entry points.  AMD-lite has no external
  library dependency; the original requires `SuiteSparse_config`.

The static-archive difference is larger because the original archive contains
many small objects, optional user-callable routines, and duplicated int64
wrappers.  The shared-library and text-size figures are the more meaningful
binary comparison.

## Performance

Representative results:

| matrix | n | nnz | AMD-lite | original AMD | ratio | peak workspace (lite/orig) |
|---|---:|---:|---:|---:|---:|---:|
| 2-D grid 100×100 | 10,000 | 39,600 | 2.41 ms | 2.30 ms | 1.05× | 0.6 / 0.5 MB |
| 2-D grid 500×500 | 250,000 | 998,000 | 82.61 ms | 78.73 ms | 1.05× | 14.1 / 13.2 MB |
| 3-D grid 20³ | 8,000 | 45,600 | 4.50 ms | 4.61 ms | 0.98× | 0.5 / 0.5 MB |
| 3-D grid 50³ | 125,000 | 735,000 | 113.91 ms | 112.59 ms | 1.01× | 8.1 / 7.7 MB |
| banded n=100k, b=5 | 100,000 | 999,970 | 27.75 ms | 27.18 ms | 1.02× | 8.4 / 8.0 MB |
| banded n=100k, b=50 | 100,000 | 9,997,450 | 389.52 ms | 373.88 ms | 1.04× | 49.6 / 49.2 MB |
| random n=100k, degree 4 | 100,000 | 400,000 | 281.93 ms | 278.32 ms | 1.01× | 5.6 / 5.3 MB |
| random n=100k, degree 10 | 100,000 | 1,000,000 | 1062.73 ms | 1061.75 ms | 1.00× | 8.4 / 8.0 MB |
| random n=100k, degree 20 | 100,000 | 2,000,000 | 1178.87 ms | 1172.11 ms | 1.01× | 13.0 / 12.6 MB |
| star n=100k | 100,000 | 199,998 | 3.33 ms | 3.11 ms | 1.07× | 4.7 / 4.3 MB |
| 2-D grid 300×300, unsorted rows | 90,000 | 358,800 | 17.60 ms | 17.01 ms | 1.03× | 5.1 / 6.4 MB |

Summary:

- Across all benchmark classes, AMD-lite ranges from **0.98× to 1.07×** the
  original runtime.
- The original's `lnz` and AMD-lite's `info.lnz` matched **exactly** on every
  benchmark matrix, so the ordering quality is the same in these tests.
- Peak workspace is close to the original.  For clean sorted input the new
  O(nnz) `AMD_aat`/`AMD_1`-style construction keeps memory within a few
  percent.  For unsorted input AMD-lite uses less memory because it
  normalizes once rather than keeping additional preprocessing copies.

### Full benchmark log

```
grid2d_100x100         n=   10000 nz=    39600  lite=   2.412 ms  orig=   2.297 ms  ratio= 1.05x  lnz lite=     196332 orig=     196332  mem lite=   0.6MB orig=   0.5MB  ok=1/1 st=0/0
grid2d_200x200         n=   40000 nz=   159200  lite=   9.281 ms  orig=   9.315 ms  ratio= 1.00x  lnz lite=    1041911 orig=    1041911  mem lite=   2.3MB orig=   2.1MB  ok=1/1 st=0/0
grid2d_300x300         n=   90000 nz=   358800  lite=  24.265 ms  orig=  23.937 ms  ratio= 1.01x  lnz lite=    2838059 orig=    2838059  mem lite=   5.1MB orig=   4.7MB  ok=1/1 st=0/0
grid2d_400x400         n=  160000 nz=   638400  lite=  48.849 ms  orig=  49.477 ms  ratio= 0.99x  lnz lite=    5503298 orig=    5503298  mem lite=   9.0MB orig=   8.4MB  ok=1/1 st=0/0
grid2d_500x500         n=  250000 nz=   998000  lite=  82.607 ms  orig=  78.730 ms  ratio= 1.05x  lnz lite=    8966158 orig=    8966158  mem lite=  14.1MB orig=  13.2MB  ok=1/1 st=0/0
grid3d_20^3            n=    8000 nz=    45600  lite=   4.504 ms  orig=   4.612 ms  ratio= 0.98x  lnz lite=     834282 orig=     834282  mem lite=   0.5MB orig=   0.5MB  ok=1/1 st=0/0
grid3d_30^3            n=   27000 nz=   156600  lite=  18.219 ms  orig=  18.149 ms  ratio= 1.00x  lnz lite=    5578774 orig=    5578774  mem lite=   1.7MB orig=   1.6MB  ok=1/1 st=0/0
grid3d_40^3            n=   64000 nz=   374400  lite=  51.977 ms  orig=  48.472 ms  ratio= 1.07x  lnz lite=   20550676 orig=   20550676  mem lite=   4.2MB orig=   3.9MB  ok=1/1 st=0/0
grid3d_50^3            n=  125000 nz=   735000  lite= 113.912 ms  orig= 112.588 ms  ratio= 1.01x  lnz lite=   61473756 orig=   61473756  mem lite=   8.1MB orig=   7.7MB  ok=1/1 st=0/0
banded_n100k_b5        n=  100000 nz=   999970  lite=  27.745 ms  orig=  27.182 ms  ratio= 1.02x  lnz lite=     499988 orig=     499988  mem lite=   8.4MB orig=   8.0MB  ok=1/1 st=0/0
banded_n100k_b20       n=  100000 nz=  3999580  lite= 119.449 ms  orig= 112.727 ms  ratio= 1.06x  lnz lite=    1999835 orig=    1999835  mem lite=  22.1MB orig=  21.7MB  ok=1/1 st=0/0
banded_n100k_b50       n=  100000 nz=  9997450  lite= 389.516 ms  orig= 373.884 ms  ratio= 1.04x  lnz lite=    4999025 orig=    4999025  mem lite=  49.6MB orig=  49.2MB  ok=1/1 st=0/0
random_n100k_d4        n=  100000 nz=   400000  lite= 281.930 ms  orig= 278.316 ms  ratio= 1.01x  lnz lite=  250141713 orig=  250141713  mem lite=   5.6MB orig=   5.3MB  ok=1/1 st=0/0
random_n100k_d10       n=  100000 nz=  1000000  lite=1062.725 ms  orig=1061.746 ms  ratio= 1.00x  lnz lite= 1280502208 orig= 1280502208  mem lite=   8.4MB orig=   8.0MB  ok=1/1 st=0/0
random_n100k_d20       n=  100000 nz=  2000000  lite=1178.865 ms  orig=1172.114 ms  ratio= 1.01x  lnz lite= 2333379851 orig= 2333379851  mem lite=  13.0MB orig=  12.6MB  ok=1/1 st=0/0
star_n100k             n=  100000 nz=   199998  lite=   3.334 ms  orig=   3.109 ms  ratio= 1.07x  lnz lite=      99999 orig=      99999  mem lite=   4.7MB orig=   4.3MB  ok=1/1 st=0/0
grid2d_300_unsorted    n=   90000 nz=   358800  lite=  17.597 ms  orig=  17.012 ms  ratio= 1.03x  lnz lite=    2838059 orig=    2838059  mem lite=   5.1MB orig=   6.4MB  ok=1/1 st=0/1
```

`ok` is a permutation check; `st` is the return status.  The original returns
`AMD_OK_BUT_JUMBLED` (1) for the unsorted grid, after preprocessing it.

## Reproduce

```sh
make
make test
./bench/run_compare.sh /path/to/SuiteSparse/AMD
```

The path may point at the AMD directory itself (containing `Source/` and
`Include/`) or at a SuiteSparse root (containing `AMD/Source/` and
`AMD/Include/`).  It can also be supplied through `AMD_SRC`:

```sh
AMD_SRC=/path/to/SuiteSparse/AMD ./bench/run_compare.sh
```

`bench/run_compare.sh` builds the original AMD with a minimal SuiteSparse
configuration stub, then compiles and runs `bench/compare_amd.c` against both
libraries.
