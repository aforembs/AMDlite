/* test_amd_lite.c: correctness and quality tests for amd_lite. */

#include "amd_lite.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int tests_run = 0 ;
static int tests_failed = 0 ;

#define CHECK(cond) do { \
    tests_run++ ; \
    if (!(cond)) { \
        tests_failed++ ; \
        fprintf (stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond) ; \
    } \
} while (0)

/* ------------------------------------------------------------------------- */
/* exact symbolic fill for small matrices                                    */
/* ------------------------------------------------------------------------- */

typedef uint64_t u64 ;

static void exact_fill
(
    int32_t n,
    const int32_t *rp,
    const int32_t *ci,
    const int32_t *perm,
    long *fill_out,
    long *lnz_out
)
{
    int W = (n + 63) / 64 ;
    size_t aw = (size_t) n * (size_t) W ;
    u64 *adj = (u64 *) calloc (aw ? aw : 1, sizeof (u64)) ;
    u64 *rem = (u64 *) calloc ((size_t) W ? (size_t) W : 1, sizeof (u64)) ;
    u64 *nb  = (u64 *) malloc (((size_t) W ? (size_t) W : 1) * sizeof (u64)) ;
    long fill = 0, lnz = 0 ;
    int32_t i, k ;

    for (i = 0 ; i < n ; i++)
    {
        int32_t p ;
        for (p = rp [i] ; p < rp [i + 1] ; p++)
        {
            int32_t j = ci [p] ;
            if (j == i) continue ;
            adj [(size_t) i * W + j / 64] |= ((u64) 1 << (j % 64)) ;
            adj [(size_t) j * W + i / 64] |= ((u64) 1 << (i % 64)) ;
        }
        rem [i / 64] |= ((u64) 1 << (i % 64)) ;
    }

    for (k = 0 ; k < n ; k++)
    {
        int32_t p = perm [k] ;
        int deg = 0, w, a, b, c = 0 ;
        int32_t *list ;
        u64 x ;

        rem [p / 64] &= ~((u64) 1 << (p % 64)) ;
        for (w = 0 ; w < W ; w++) nb [w] = adj [(size_t) p * W + w] & rem [w] ;
        for (w = 0 ; w < W ; w++) deg += __builtin_popcountll (nb [w]) ;
        lnz += deg ;
        list = (int32_t *) malloc ((size_t) (deg ? deg : 1) * sizeof (int32_t)) ;
        for (w = 0 ; w < W ; w++)
        {
            x = nb [w] ;
            while (x)
            {
                int bit = __builtin_ctzll (x) ;
                x &= x - 1 ;
                list [c++] = (int32_t) (w * 64 + bit) ;
            }
        }
        for (a = 0 ; a < deg ; a++)
        {
            for (b = a + 1 ; b < deg ; b++)
            {
                int32_t u = list [a], v = list [b] ;
                if (!(adj [(size_t) u * W + v / 64] & ((u64) 1 << (v % 64))))
                {
                    adj [(size_t) u * W + v / 64] |= ((u64) 1 << (v % 64)) ;
                    adj [(size_t) v * W + u / 64] |= ((u64) 1 << (u % 64)) ;
                    fill++ ;
                }
            }
        }
        for (a = 0 ; a < deg ; a++)
        {
            int32_t u = list [a] ;
            adj [(size_t) u * W + p / 64] &= ~((u64) 1 << (p % 64)) ;
        }
        for (w = 0 ; w < W ; w++) adj [(size_t) p * W + w] = 0 ;
        free (list) ;
    }

    free (adj) ;
    free (rem) ;
    free (nb) ;
    *fill_out = fill ;
    *lnz_out = lnz ;
}

static int valid_perm (int32_t n, const int32_t *perm)
{
    char *seen = (char *) calloc ((size_t) (n ? n : 1), 1) ;
    int32_t i ;
    int ok = 1 ;
    for (i = 0 ; i < n ; i++)
    {
        if (perm [i] < 0 || perm [i] >= n || seen [perm [i]]) ok = 0 ;
        else seen [perm [i]] = 1 ;
    }
    free (seen) ;
    return ok ;
}

/* ------------------------------------------------------------------------- */
/* matrix generators                                                         */
/* ------------------------------------------------------------------------- */

static void make_grid (int m, int32_t **rp_out, int32_t **ci_out, int32_t *n_out)
{
    int32_t n = (int32_t) m * m ;
    int32_t *rp = (int32_t *) malloc (((size_t) n + 1) * sizeof (int32_t)) ;
    int32_t *ci = (int32_t *) malloc ((size_t) n * 4 * sizeof (int32_t)) ;
    int32_t nz = 0, i, j ;
    rp [0] = 0 ;
    for (i = 0 ; i < m ; i++)
    {
        for (j = 0 ; j < m ; j++)
        {
            int32_t u = i * m + j ;
            rp [u] = nz ;
            if (i > 0)      ci [nz++] = (i - 1) * m + j ;
            if (i < m - 1)  ci [nz++] = (i + 1) * m + j ;
            if (j > 0)      ci [nz++] = i * m + j - 1 ;
            if (j < m - 1)  ci [nz++] = i * m + j + 1 ;
        }
    }
    rp [n] = nz ;
    *rp_out = rp ;
    *ci_out = ci ;
    *n_out = n ;
}

static void make_random
(
    int32_t n,
    int deg,
    uint64_t seed,
    int32_t **rp_out,
    int32_t **ci_out
)
{
    int32_t *rp = (int32_t *) malloc (((size_t) n + 1) * sizeof (int32_t)) ;
    int32_t *ci = (int32_t *) malloc (((size_t) n * (size_t) deg * 2 + 1)
        * sizeof (int32_t)) ;
    int32_t nz = 0, i, d ;
    rp [0] = 0 ;
    for (i = 0 ; i < n ; i++)
    {
        rp [i] = nz ;
        for (d = 0 ; d < deg ; d++)
        {
            int32_t j ;
            seed = seed * 6364136223846793005ULL + 1 ;
            j = (int32_t) ((seed >> 33) % (uint64_t) n) ;
            if (j != i) ci [nz++] = j ;
        }
    }
    rp [n] = nz ;
    *rp_out = rp ;
    *ci_out = ci ;
}

static void make_star (int32_t n, int32_t **rp_out, int32_t **ci_out)
{
    int32_t *rp = (int32_t *) malloc (((size_t) n + 1) * sizeof (int32_t)) ;
    int32_t *ci = (int32_t *) malloc ((size_t) 2 * (n - 1) * sizeof (int32_t)) ;
    int32_t i, nz = 0 ;
    rp [0] = 0 ;
    for (i = 0 ; i < n ; i++)
    {
        rp [i] = nz ;
        if (i == 0)
        {
            int32_t j ;
            for (j = 1 ; j < n ; j++) ci [nz++] = j ;
        }
        else
        {
            ci [nz++] = 0 ;
        }
    }
    rp [n] = nz ;
    *rp_out = rp ;
    *ci_out = ci ;
}

/* ------------------------------------------------------------------------- */
/* tests                                                                     */
/* ------------------------------------------------------------------------- */

static void test_edge_cases (void)
{
    int32_t perm [1] = { 0 } ;
    amd_lite_info info ;

    /* n = 0: no arrays needed */
    CHECK (amd_lite_order (0, NULL, NULL, NULL, NULL, &info) == AMD_LITE_OK) ;

    /* one isolated row */
    {
        int32_t rp [2] = { 0, 0 } ;
        CHECK (amd_lite_order (1, rp, NULL, perm, NULL, &info) == AMD_LITE_OK) ;
        CHECK (perm [0] == 0) ;
    }

    /* one diagonal entry */
    {
        int32_t rp [2] = { 0, 1 } ;
        int32_t ci [1] = { 0 } ;
        CHECK (amd_lite_order (1, rp, ci, perm, NULL, &info) == AMD_LITE_OK) ;
        CHECK (perm [0] == 0) ;
        CHECK (info.nzdiag == 1) ;
        CHECK (info.nz_a_plus_at == 0) ;
    }

    /* empty matrix with n = 6 */
    {
        int32_t n = 6, i ;
        int32_t rp [7] = { 0, 0, 0, 0, 0, 0, 0 } ;
        int32_t p [6] ;
        CHECK (amd_lite_order (n, rp, NULL, p, NULL, &info) == AMD_LITE_OK) ;
        for (i = 0 ; i < n ; i++) CHECK (p [i] == i) ;
    }
}

static void test_path_and_star (void)
{
    /* path 0-1-2-3-4 */
    {
        int32_t rp [6] = { 0, 1, 3, 5, 7, 8 } ;
        int32_t ci [8] = { 1, 0, 2, 1, 3, 2, 4, 3 } ;
        int32_t p [5] ;
        amd_lite_info info ;
        long fill, lnz, nfill, nlnz ;
        int32_t nat [5] = { 0, 1, 2, 3, 4 } ;
        CHECK (amd_lite_order (5, rp, ci, p, NULL, &info) == AMD_LITE_OK) ;
        CHECK (valid_perm (5, p)) ;
        exact_fill (5, rp, ci, p, &fill, &lnz) ;
        exact_fill (5, rp, ci, nat, &nfill, &nlnz) ;
        CHECK (fill == 0) ;
        CHECK (lnz == 4) ;
        CHECK (info.lnz == (double) lnz) ;
    }

    /* star: center 0 connected to all others */
    {
        int32_t n = 50 ;
        int32_t *rp, *ci, *p, *nat ;
        amd_lite_info info ;
        long fill, lnz, nfill, nlnz ;
        int32_t i ;
        make_star (n, &rp, &ci) ;
        p = (int32_t *) malloc ((size_t) n * sizeof (int32_t)) ;
        nat = (int32_t *) malloc ((size_t) n * sizeof (int32_t)) ;
        for (i = 0 ; i < n ; i++) nat [i] = i ;
        CHECK (amd_lite_order (n, rp, ci, p, NULL, &info) == AMD_LITE_OK) ;
        CHECK (valid_perm (n, p)) ;
        CHECK (p [n - 1] == 0) ;
        exact_fill (n, rp, ci, p, &fill, &lnz) ;
        exact_fill (n, rp, ci, nat, &nfill, &nlnz) ;
        CHECK (fill == 0) ;
        CHECK (nfill > fill) ;
        free (rp) ;
        free (ci) ;
        free (p) ;
        free (nat) ;
    }
}

static void test_grid_quality (void)
{
    int m ;
    for (m = 5 ; m <= 20 ; m += 5)
    {
        int32_t *rp, *ci, *p, *nat ;
        int32_t n, i ;
        long fill, lnz, nfill, nlnz ;
        amd_lite_info info ;
        make_grid (m, &rp, &ci, &n) ;
        p = (int32_t *) malloc ((size_t) n * sizeof (int32_t)) ;
        nat = (int32_t *) malloc ((size_t) n * sizeof (int32_t)) ;
        for (i = 0 ; i < n ; i++) nat [i] = i ;
        CHECK (amd_lite_order (n, rp, ci, p, NULL, &info) == AMD_LITE_OK) ;
        CHECK (valid_perm (n, p)) ;
        exact_fill (n, rp, ci, p, &fill, &lnz) ;
        exact_fill (n, rp, ci, nat, &nfill, &nlnz) ;
        CHECK (fill < nfill) ;
        CHECK (info.lnz == (double) lnz) ;
        free (rp) ;
        free (ci) ;
        free (p) ;
        free (nat) ;
    }
}

static void test_random (void)
{
    int n, d ;
    for (n = 20 ; n <= 160 ; n += 20)
    {
        for (d = 2 ; d <= 6 ; d += 2)
        {
            int32_t *rp, *ci, *p ;
            long fill, lnz ;
            amd_lite_info info ;
            make_random ((int32_t) n, d, (uint64_t) n * 131 + (uint64_t) d,
                &rp, &ci) ;
            p = (int32_t *) malloc ((size_t) n * sizeof (int32_t)) ;
            CHECK (amd_lite_order ((int32_t) n, rp, ci, p, NULL, &info)
                == AMD_LITE_OK) ;
            CHECK (valid_perm ((int32_t) n, p)) ;
            exact_fill ((int32_t) n, rp, ci, p, &fill, &lnz) ;
            CHECK (lnz >= 0) ;
            free (rp) ;
            free (ci) ;
            free (p) ;
        }
    }
}

static void test_options (void)
{
    int32_t n = 120 ;
    int32_t *rp, *ci, *p1, *p2, *p3 ;
    amd_lite_options opt ;
    int32_t i ;

    make_random (n, 4, 987654321, &rp, &ci) ;
    p1 = (int32_t *) malloc ((size_t) n * sizeof (int32_t)) ;
    p2 = (int32_t *) malloc ((size_t) n * sizeof (int32_t)) ;
    p3 = (int32_t *) malloc ((size_t) n * sizeof (int32_t)) ;

    CHECK (amd_lite_order (n, rp, ci, p1, NULL, NULL) == AMD_LITE_OK) ;
    CHECK (valid_perm (n, p1)) ;

    amd_lite_defaults (&opt) ;
    opt.aggressive = 0 ;
    CHECK (amd_lite_order (n, rp, ci, p2, &opt, NULL) == AMD_LITE_OK) ;
    CHECK (valid_perm (n, p2)) ;

    opt.dense = -1 ;
    CHECK (amd_lite_order (n, rp, ci, p3, &opt, NULL) == AMD_LITE_OK) ;
    CHECK (valid_perm (n, p3)) ;

    /* CRC and CSC give the same graph, and therefore the same ordering */
    amd_lite_defaults (&opt) ;
    opt.format = AMD_LITE_CSC ;
    CHECK (amd_lite_order (n, rp, ci, p2, &opt, NULL) == AMD_LITE_OK) ;
    for (i = 0 ; i < n ; i++) CHECK (p1 [i] == p2 [i]) ;

    free (rp) ;
    free (ci) ;
    free (p1) ;
    free (p2) ;
    free (p3) ;
}

static void test_duplicates_and_unsorted (void)
{
    int32_t n = 8 ;
    /* unsorted, with duplicate off-diagonal entries */
    int32_t rp [9] = { 0, 3, 6, 9, 12, 15, 18, 21, 24 } ;
    int32_t ci [24] = {
        5, 1, 1,  0, 2, 0,  1, 3, 1,  2, 4, 2,
        3, 5, 3,  4, 6, 6,  5, 7, 5,  6, 0, 6
    } ;
    int32_t *p1 = (int32_t *) malloc ((size_t) n * sizeof (int32_t)) ;
    int32_t *p2 = (int32_t *) malloc ((size_t) n * sizeof (int32_t)) ;
    amd_lite_info info ;

    CHECK (amd_lite_order (n, rp, ci, p1, NULL, &info) == AMD_LITE_OK) ;
    CHECK (valid_perm (n, p1)) ;
    /* on the normalized pattern, 14 of 16 directed entries have a reverse */
    CHECK (info.symmetry > 14.0 / 16.0 - 1e-12) ;
    CHECK (info.symmetry < 14.0 / 16.0 + 1e-12) ;

    /* the same graph without duplicates, sorted */
    {
        int32_t rp2 [9] = { 0, 3, 5, 7, 9, 11, 14, 16, 18 } ;
        int32_t ci2 [18] = {
            1, 5, 7,  0, 2,  1, 3,  2, 4,
            3, 5,  0, 4, 6,  5, 7,  0, 6
        } ;
        CHECK (amd_lite_order (n, rp2, ci2, p2, NULL, NULL) == AMD_LITE_OK) ;
        for (int32_t i = 0 ; i < n ; i++) CHECK (p1 [i] == p2 [i]) ;
    }

    free (p1) ;
    free (p2) ;
}

static void test_dense (void)
{
    int32_t n = 500 ;
    int32_t *rp, *ci, *p ;
    amd_lite_info info ;
    long fill, lnz ;

    make_star (n, &rp, &ci) ;
    p = (int32_t *) malloc ((size_t) n * sizeof (int32_t)) ;
    CHECK (amd_lite_order (n, rp, ci, p, NULL, &info) == AMD_LITE_OK) ;
    CHECK (valid_perm (n, p)) ;
    CHECK (info.ndense >= 1) ;
    CHECK (p [n - 1] == 0) ;
    exact_fill (n, rp, ci, p, &fill, &lnz) ;
    CHECK (fill == 0) ;
    free (rp) ;
    free (ci) ;
    free (p) ;
}

static void test_int64 (void)
{
    int32_t n = 60 ;
    int32_t *rp32, *ci32, *p32 ;
    int64_t *rp64, *ci64, *p64 ;
    int64_t i ;
    int64_t nz ;

    make_random (n, 4, 424242, &rp32, &ci32) ;
    nz = rp32 [n] ;
    rp64 = (int64_t *) malloc (((size_t) n + 1) * sizeof (int64_t)) ;
    ci64 = (int64_t *) malloc ((size_t) nz * sizeof (int64_t)) ;
    p32 = (int32_t *) malloc ((size_t) n * sizeof (int32_t)) ;
    p64 = (int64_t *) malloc ((size_t) n * sizeof (int64_t)) ;

    for (i = 0 ; i <= n ; i++) rp64 [i] = rp32 [i] ;
    for (i = 0 ; i < nz ; i++) ci64 [i] = ci32 [i] ;

    CHECK (amd_lite_order (n, rp32, ci32, p32, NULL, NULL) == AMD_LITE_OK) ;
    CHECK (amd_lite_l_order ((int64_t) n, rp64, ci64, p64, NULL, NULL)
        == AMD_LITE_OK) ;
    for (i = 0 ; i < n ; i++) CHECK (p32 [i] == p64 [i]) ;

    free (rp32) ;
    free (ci32) ;
    free (p32) ;
    free (rp64) ;
    free (ci64) ;
    free (p64) ;
}

static void test_invalid (void)
{
    int32_t p [4] = { 0, 0, 0, 0 } ;
    amd_lite_info info ;

    /* decreasing pointer */
    {
        int32_t rp [4] = { 0, 2, 1, 3 } ;
        int32_t ci [3] = { 0, 1, 2 } ;
        CHECK (amd_lite_order (3, rp, ci, p, NULL, &info)
            == AMD_LITE_INVALID) ;
    }
    /* row index out of range */
    {
        int32_t rp [3] = { 0, 1, 2 } ;
        int32_t ci [2] = { 0, 7 } ;
        CHECK (amd_lite_order (2, rp, ci, p, NULL, &info)
            == AMD_LITE_INVALID) ;
    }
    /* negative row index */
    {
        int32_t rp [3] = { 0, 1, 2 } ;
        int32_t ci [2] = { 0, -1 } ;
        CHECK (amd_lite_order (2, rp, ci, p, NULL, &info)
            == AMD_LITE_INVALID) ;
    }
    /* first pointer not zero */
    {
        int32_t rp [2] = { 1, 1 } ;
        int32_t ci [1] = { 0 } ;
        CHECK (amd_lite_order (1, rp, ci, p, NULL, &info)
            == AMD_LITE_INVALID) ;
    }
    /* missing permutation */
    {
        int32_t rp [2] = { 0, 0 } ;
        CHECK (amd_lite_order (1, rp, NULL, NULL, NULL, &info)
            == AMD_LITE_INVALID) ;
    }
    /* bad option value */
    {
        int32_t rp [2] = { 0, 0 } ;
        amd_lite_options opt ;
        amd_lite_defaults (&opt) ;
        opt.format = 99 ;
        CHECK (amd_lite_order (1, rp, NULL, p, &opt, &info)
            == AMD_LITE_INVALID) ;
    }
}

int main (void)
{
    test_edge_cases () ;
    test_path_and_star () ;
    test_grid_quality () ;
    test_random () ;
    test_options () ;
    test_duplicates_and_unsorted () ;
    test_dense () ;
    test_int64 () ;
    test_invalid () ;

    printf ("%d checks, %d failed\n", tests_run, tests_failed) ;
    return (tests_failed == 0) ? 0 : 1 ;
}
