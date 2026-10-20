# AMD-lite

A lightweight, dependency-free C implementation of the **approximate minimum
degree (AMD)** fill-reducing ordering algorithm for sparse matrices.

`AMD-lite` keeps the parts that matter for ordering -- the quotient-graph AMD
kernel with supervariable detection, mass elimination, element absorption,
aggressive absorption, dense-row removal, and assembly-tree post-ordering --
and drops the SuiteSparse configuration layer, the Fortran interface, the
MATLAB interface, and the duplicated int64 source files.

Input is **CSR by default**.  Internally the pattern of `A+A'` is formed, so
the routine also works for unsymmetric matrices.  Both `int32_t` and `int64_t`
entry points are provided.

## Algorithm

The method is the one described in

> P. R. Amestoy, T. A. Davis, and I. S. Duff,
> *An approximate minimum degree ordering algorithm*,
> SIAM J. Matrix Anal. Appl., 17(4):886-905, 1996.

At each step the pivot is the supervariable with the smallest approximate
external degree.  The implementation uses a quotient graph, element
absorption, mass elimination, optional aggressive absorption, and a
depth-first post-ordering of the assembly tree.

## Build

```sh
make                # build/libamd_lite.a and build/libamd_lite.so
make test           # build and run the test suite
make demo           # order a 20-by-20 five-point grid
```

or with CMake:

```sh
cmake -S . -B build-cmake
cmake --build build-cmake
ctest --test-dir build-cmake
```

The library has no dependencies other than the C standard library and `libm`.

## Usage

```c
#include "amd_lite.h"

/* A is n-by-n in CSR form: rowptr[0..n], colidx[0..nnz-1].
 * Rows need not be sorted; duplicate entries are allowed. */
int32_t n = ...;
int32_t *rowptr = ...;
int32_t *colidx = ...;
int32_t *perm = malloc((size_t) n * sizeof(int32_t));

amd_lite_options opt;
amd_lite_defaults(&opt);       /* dense = 10.0, aggressive = 1, CSR */
/* opt.dense = -1.0;           turn off dense-row handling */
/* opt.aggressive = 0;         use AMDBAR-style absorption */

amd_lite_info info;
int status = amd_lite_order(n, rowptr, colidx, perm, &opt, &info);

if (status == AMD_LITE_OK) {
    /* perm[k] is the k-th pivot row/column. */
}
free(perm);
```

For `int64_t` indices use `amd_lite_l_order`; the API is otherwise identical.

`options` and `info` may be `NULL`.  Passing `options = NULL` selects the
defaults.  The input arrays are not modified.

### Return codes

| code | meaning |
| --- | --- |
| `AMD_LITE_OK` | success |
| `AMD_LITE_INVALID` | invalid matrix or arguments |
| `AMD_LITE_OUT_OF_MEMORY` | allocation failed or problem too large |

### Options

| field | default | meaning |
| --- | --- | --- |
| `dense` | `10.0` | a row/column is dense if its degree in `A+A'` exceeds `dense*sqrt(n)` (degrees `<= 16` are never dense); set `< 0` to disable |
| `aggressive` | `1` | nonzero enables aggressive element absorption |
| `format` | `AMD_LITE_CSR` | `AMD_LITE_CSC` may be used for column-oriented input; the ordering is invariant under transposition because only the pattern of `A+A'` matters |

### Statistics

`amd_lite_info` reports `n`, `nz`, `nzdiag`, `nz_a_plus_at`, `symmetry`,
`ndense`, `ncmpa`, `dmax`, approximate `lnz`, `ndiv`, multiply-subtract
counts for LDL' and LU, and the number of bytes allocated by the call.

## Memory

For a problem with `nzA = |pattern(A+A')|` off-diagonal entries and an index
size of `S` bytes, the ordering kernel allocates approximately

```
(1.2 * nzA + 9 * n) * S  bytes
```

The pattern-construction phase temporarily needs one normalized copy of the
input (about `nnz` indices) plus up to `2*nnz` entries for the union pattern.
For symmetric input the union is about `nnz` entries.  The symmetry statistic
is only computed when an `amd_lite_info` output is requested.

## Differences from SuiteSparse AMD

* CSR-first public API (`csr` row pointers and column indices).
* Rows may be unsorted and may contain duplicates; they are normalized
  internally.
* No `SuiteSparse_config`, Fortran, or MATLAB dependencies.
* Approximate degree and fill statistics are reported; some debug/control
  routines of the original package are omitted.
* The core ordering logic is the quotient-graph AMD algorithm from the
  original package, specialized to 32-bit and 64-bit indices.

## License

BSD 3-clause.  This is a derivative of the SuiteSparse AMD package by
Timothy A. Davis, Patrick R. Amestoy, and Iain S. Duff.  See `LICENSE`.

## Code size and performance

See [`comparison.md`](comparison.md) for the full comparison and complete
benchmark log.  Summary:

Library source, excluding tests and the demo:

| | lines | shared-library text | shared-library file |
| --- | ---: | ---: | ---: |
| AMD-lite | 1,632 | 26.1 KB | 36.7 KB |
| SuiteSparse AMD, full | 4,194 | 35.3 KB | 50.3 KB |
| SuiteSparse AMD, core* | 3,771 | 29.8 KB | 41.6 KB |

\* core = no `amd_control`, `amd_info`, `amd_dump`, or `amd_version`.
AMD-lite is about 39% of the full original source size and about 74% of its
shared-library text size.  Both were built with `gcc -O3` and include the
32-bit and 64-bit entry points.

The benchmark in `bench/compare_amd.c` uses identical CSR input for both
libraries.  `bench/run_compare.sh` builds a SuiteSparse AMD source tree
supplied by the caller (it is not bundled with AMD-lite) using a minimal
`SuiteSparse_config` stub.  Representative results on a pinned CPU, best of
3-7 runs:

| matrix | AMD-lite | original AMD | ratio | peak workspace (lite/orig) |
| --- | ---: | ---: | ---: | ---: |
| 2-D grid, 500x500 | 82.6 ms | 78.7 ms | 1.05x | 14.1 / 13.2 MB |
| 3-D grid, 50^3 | 113.9 ms | 112.6 ms | 1.01x | 8.1 / 7.7 MB |
| banded, n=100k, b=50 | 389.5 ms | 373.9 ms | 1.04x | 49.6 / 49.2 MB |
| random, n=100k, degree 10 | 1062.7 ms | 1061.7 ms | 1.00x | 8.4 / 8.0 MB |
| star, n=100k | 3.33 ms | 3.11 ms | 1.07x | 4.7 / 4.3 MB |
| 2-D grid, unsorted rows | 17.6 ms | 17.0 ms | 1.03x | 5.1 / 6.4 MB |

Across all benchmark classes the ratio stayed between 0.98x and 1.07x, and
`info.lnz` matched the original `Info[AMD_LNZ]` exactly.  The unsorted-input
case also shows the benefit of AMD-lite's lower preprocessing memory.

To reproduce, point the benchmark at a SuiteSparse AMD checkout:

```sh
./bench/run_compare.sh /path/to/SuiteSparse/AMD
# or
AMD_SRC=/path/to/SuiteSparse/AMD ./bench/run_compare.sh
```
