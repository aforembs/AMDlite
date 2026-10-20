/* amd_lite.h: lightweight approximate-minimum-degree ordering
 *
 * A compact, dependency-free implementation of the AMD ordering algorithm
 * for symmetric or unsymmetric sparse matrices supplied in CSR form.
 *
 * Algorithm: P. R. Amestoy, T. A. Davis, and I. S. Duff,
 * "An approximate minimum degree ordering algorithm", SIAM J. Matrix Anal.
 * Appl., 17(4):886-905, 1996.
 *
 * This file is a lightweight derivative of the SuiteSparse AMD package.
 * AMD, Copyright (c) 1996-2024, Timothy A. Davis, Patrick R. Amestoy, and
 * Iain S. Duff.  All Rights Reserved.  SPDX-License-Identifier: BSD-3-clause
 */

#ifndef AMD_LITE_H
#define AMD_LITE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define AMD_LITE_VERSION_MAJOR 1
#define AMD_LITE_VERSION_MINOR 0
#define AMD_LITE_VERSION_PATCH 0

/* return codes */
#define AMD_LITE_OK            0
#define AMD_LITE_INVALID      -1
#define AMD_LITE_OUT_OF_MEMORY -2

/* input formats.  The pattern of A+A' (and hence the ordering) is identical
 * for CSR and CSC input, so this flag is informational.  CSR is the default. */
#define AMD_LITE_CSR 0
#define AMD_LITE_CSC 1

/* default control values */
#define AMD_LITE_DEFAULT_DENSE      10.0
#define AMD_LITE_DEFAULT_AGGRESSIVE 1

/* ------------------------------------------------------------------------- */
/* options and statistics                                                    */
/* ------------------------------------------------------------------------- */

typedef struct
{
    /* A row/column is treated as "dense" when its degree in A+A' exceeds
     * dense * sqrt(n), but never for degrees <= 16.  Set dense < 0 to disable
     * dense-row handling.  Default: AMD_LITE_DEFAULT_DENSE (10.0). */
    double dense;

    /* Nonzero to use aggressive element absorption.  Default: 1. */
    int aggressive;

    /* AMD_LITE_CSR (0) or AMD_LITE_CSC (1).  Default: AMD_LITE_CSR. */
    int format;

} amd_lite_options;

typedef struct
{
    int      status;         /* one of the AMD_LITE_* return codes          */
    int64_t  n;              /* order of the matrix                         */
    int64_t  nz;             /* number of entries in the input              */
    int64_t  nzdiag;         /* number of diagonal entries in the input     */
    int64_t  nz_a_plus_at;   /* off-diagonal entries in pattern(A+A')       */
    int64_t  ndense;         /* number of dense rows/columns removed        */
    int64_t  ncmpa;          /* number of garbage collections               */
    int64_t  dmax;           /* max nonzeros in a column of L, incl. diag.  */
    double   symmetry;       /* matched off-diagonal fraction, in [0,1],
                                computed after duplicates are removed       */
    double   lnz;            /* approx. nonzeros in L, excluding diagonal   */
    double   ndiv;           /* divides for subsequent LDL'/LU factorization*/
    double   nmultsubs_ldl;  /* multiply-subtract pairs, LDL'               */
    double   nmultsubs_lu;   /* multiply-subtract pairs, LU                 */
    double   memory;         /* bytes allocated for the ordering workspace  */
} amd_lite_info;

/* ------------------------------------------------------------------------- */
/* API                                                                       */
/* ------------------------------------------------------------------------- */

/* Fill an options struct with the defaults. */
void amd_lite_defaults (amd_lite_options *options);

/* Compute an AMD ordering of the pattern of A+A'.
 *
 *   n:       order of the matrix (n >= 0).
 *   rowptr:  array of n+1 pointers into colidx.  rowptr[0] must be 0 and
 *            rowptr must be nondecreasing.
 *   colidx:  array of nnz column indices in [0,n).  Rows need not be sorted
 *            and duplicate entries are allowed.
 *   perm:    output permutation of size n.  Row/column i is the k-th pivot
 *            when perm[k] = i.  In MATLAB notation the reordered matrix is
 *            A(perm+1,perm+1).
 *   options: optional control parameters (NULL means defaults).
 *   info:    optional statistics (NULL is allowed).
 *
 * The input arrays are not modified.  Returns AMD_LITE_OK, AMD_LITE_INVALID,
 * or AMD_LITE_OUT_OF_MEMORY.
 */
int amd_lite_order
(
    int32_t              n,
    const int32_t       *rowptr,
    const int32_t       *colidx,
    int32_t             *perm,
    const amd_lite_options *options,
    amd_lite_info           *info
) ;

/* int64_t variant of amd_lite_order. */
int amd_lite_l_order
(
    int64_t              n,
    const int64_t       *rowptr,
    const int64_t       *colidx,
    int64_t             *perm,
    const amd_lite_options *options,
    amd_lite_info           *info
) ;

#ifdef __cplusplus
}
#endif

#endif /* AMD_LITE_H */
