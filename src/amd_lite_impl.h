/* amd_lite_impl.h: type-generic implementation of the AMD ordering kernel.
 *
 * This file is included twice by amd_lite.c, once for int32_t and once for
 * int64_t indices.  The including file must define:
 *
 *   AMDL_IDX       index type (int32_t or int64_t)
 *   AMDL_UINT      corresponding unsigned type
 *   AMDL_IDX_MAX   largest value of AMDL_IDX
 *   AMDL_PREFIX    prefix for internal symbols
 *   AMDL_PUBLIC    name of the public ordering function
 *
 * AMD, Copyright (c) 1996-2024, Timothy A. Davis, Patrick R. Amestoy, and
 * Iain S. Duff.  All Rights Reserved.  SPDX-License-Identifier: BSD-3-clause
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>

#define AMDL_EMPTY  ((AMDL_IDX) -1)
#define AMDL_FLIP(i) ((AMDL_IDX) (-(i) - 2))

/* ------------------------------------------------------------------------- */
/* small helpers                                                             */
/* ------------------------------------------------------------------------- */

static void *AMDL_FN(xmalloc) (size_t nmemb, size_t size)
{
    if (nmemb == 0) nmemb = 1;
    if (size != 0 && nmemb > SIZE_MAX / size) return (NULL);
    return (malloc (nmemb * size)) ;
}

static int AMDL_FN(cmp_index) (const void *pa, const void *pb)
{
    AMDL_IDX a = *(const AMDL_IDX *) pa ;
    AMDL_IDX b = *(const AMDL_IDX *) pb ;
    return ((a > b) - (a < b)) ;
}

/* ------------------------------------------------------------------------- */
/* input validation                                                          */
/* ------------------------------------------------------------------------- */

static int AMDL_FN(validate)
(
    AMDL_IDX n,
    const AMDL_IDX *rowptr,
    const AMDL_IDX *colidx,
    int64_t *nz_out,
    int64_t *nzdiag_out,
    int *clean_out
)
{
    AMDL_IDX i ;
    int64_t nz, nzdiag = 0 ;
    int clean = 1 ;

    *nz_out = 0 ;
    *nzdiag_out = 0 ;
    *clean_out = 1 ;
    if (n < 0) return (AMD_LITE_INVALID) ;
    if (n == 0) return (AMD_LITE_OK) ;
    if (rowptr == NULL) return (AMD_LITE_INVALID) ;
    if (rowptr [0] != 0) return (AMD_LITE_INVALID) ;
    if (rowptr [n] > 0 && colidx == NULL) return (AMD_LITE_INVALID) ;

    for (i = 0 ; i < n ; i++)
    {
        AMDL_IDX p1 = rowptr [i] ;
        AMDL_IDX p2 = rowptr [i + 1] ;
        AMDL_IDX p, prev = -1 ;
        if (p1 < 0 || p2 < p1) return (AMD_LITE_INVALID) ;
        for (p = p1 ; p < p2 ; p++)
        {
            AMDL_IDX j = colidx [p] ;
            if (j < 0 || j >= n) return (AMD_LITE_INVALID) ;
            if (j <= prev) clean = 0 ;
            prev = j ;
            if (j == i) nzdiag++ ;
        }
    }
    nz = (int64_t) rowptr [n] ;
    if (nz < 0) return (AMD_LITE_INVALID) ;
    *nz_out = nz ;
    *nzdiag_out = nzdiag ;
    *clean_out = clean ;
    return (AMD_LITE_OK) ;
}

/* sort a single row; insertion sort is faster for the small rows that
 * dominate typical sparse matrices */
static void AMDL_FN(sort_row) (AMDL_IDX *a, AMDL_IDX n)
{
    AMDL_IDX i, j, v ;
    if (n < 2) return ;
    if (n <= 32)
    {
        for (i = 1 ; i < n ; i++)
        {
            v = a [i] ;
            j = i ;
            while (j > 0 && a [j - 1] > v)
            {
                a [j] = a [j - 1] ;
                j-- ;
            }
            a [j] = v ;
        }
    }
    else
    {
        qsort (a, (size_t) n, sizeof (AMDL_IDX), AMDL_FN(cmp_index)) ;
    }
}

static int AMDL_FN(build_pattern)
(
    AMDL_IDX n,
    const AMDL_IDX *rowptr,
    const AMDL_IDX *colidx,
    int64_t nz,
    int clean_input,
    AMDL_IDX **out_rp,
    AMDL_IDX **out_rj,
    int64_t *out_nzaat,
    double *out_symmetry,
    int need_symmetry
)
{
    const AMDL_IDX *Cp = NULL, *Ci = NULL ;
    AMDL_IDX *Cp_own = NULL, *Ci_own = NULL, *Clen = NULL ;
    AMDL_IDX *Tp = NULL, *Len = NULL, *rp = NULL, *rj = NULL ;
    AMDL_IDX i, j, k, p, pj ;
    uint64_t sum, nz_norm ;
    int64_t nzdiag = 0, nzboth = 0 ;
    double symmetry = 1.0 ;

    *out_rp = NULL ;
    *out_rj = NULL ;
    *out_nzaat = 0 ;
    *out_symmetry = 1.0 ;

    /* n+1 and nz must be representable as allocation sizes */
    if ((uint64_t) n > (uint64_t) (SIZE_MAX / sizeof (AMDL_IDX)) - 2 ||
        (uint64_t) nz > (uint64_t) (SIZE_MAX / sizeof (AMDL_IDX)))
    {
        return (AMD_LITE_OUT_OF_MEMORY) ;
    }

    if (nz == 0)
    {
        rp = AMDL_FN(xmalloc) ((size_t) n + 1, sizeof (AMDL_IDX)) ;
        if (rp == NULL) return (AMD_LITE_OUT_OF_MEMORY) ;
        for (i = 0 ; i <= n ; i++) rp [i] = 0 ;
        *out_rp = rp ;
        return (AMD_LITE_OK) ;
    }

    /* --------------------------------------------------------------------- */
    /* normalize the input rows: sort and remove duplicates if necessary     */
    /* --------------------------------------------------------------------- */

    if (clean_input)
    {
        Cp = rowptr ;
        Ci = colidx ;
        nz_norm = (uint64_t) nz ;
    }
    else
    {
        uint64_t dst ;

        Cp_own = AMDL_FN(xmalloc) ((size_t) n + 1, sizeof (AMDL_IDX)) ;
        Ci_own = AMDL_FN(xmalloc) ((size_t) nz, sizeof (AMDL_IDX)) ;
        Clen = AMDL_FN(xmalloc) ((size_t) n, sizeof (AMDL_IDX)) ;
        if (Cp_own == NULL || Ci_own == NULL || Clen == NULL) goto done ;

        memcpy (Cp_own, rowptr, ((size_t) n + 1) * sizeof (AMDL_IDX)) ;
        memcpy (Ci_own, colidx, (size_t) nz * sizeof (AMDL_IDX)) ;

        for (i = 0 ; i < n ; i++)
        {
            AMDL_IDX start = Cp_own [i] ;
            AMDL_IDX len = Cp_own [i + 1] - start ;
            AMDL_IDX w = 0, q ;
            int row_clean = 1 ;

            for (q = 1 ; q < len ; q++)
            {
                if (Ci_own [start + q - 1] >= Ci_own [start + q])
                {
                    row_clean = 0 ;
                    break ;
                }
            }
            if (row_clean)
            {
                Clen [i] = len ;
                continue ;
            }

            AMDL_FN(sort_row) (Ci_own + start, len) ;
            for (q = 0 ; q < len ; q++)
            {
                AMDL_IDX v = Ci_own [start + q] ;
                if (w == 0 || v != Ci_own [start + w - 1])
                {
                    Ci_own [start + w++] = v ;
                }
            }
            Clen [i] = w ;
        }

        /* compact the normalized rows into contiguous storage */
        dst = 0 ;
        for (i = 0 ; i < n ; i++)
        {
            AMDL_IDX start = Cp_own [i] ;
            AMDL_IDX len = Clen [i] ;
            if (dst != (uint64_t) start && len > 0)
            {
                memmove (Ci_own + dst, Ci_own + start,
                    (size_t) len * sizeof (AMDL_IDX)) ;
            }
            Cp_own [i] = (AMDL_IDX) dst ;
            dst += (uint64_t) len ;
            if (dst > (uint64_t) AMDL_IDX_MAX) goto done ;
        }
        Cp_own [n] = (AMDL_IDX) dst ;
        nz_norm = dst ;
        free (Clen) ;
        Clen = NULL ;
        Cp = Cp_own ;
        Ci = Ci_own ;
    }

    /* --------------------------------------------------------------------- */
    /* AMD_aat-style scan: union degrees and symmetry in O(nnz)              */
    /* --------------------------------------------------------------------- */

    Tp = AMDL_FN(xmalloc) ((size_t) n, sizeof (AMDL_IDX)) ;
    Len = AMDL_FN(xmalloc) ((size_t) n, sizeof (AMDL_IDX)) ;
    if (Tp == NULL || Len == NULL) goto done ;
    for (i = 0 ; i < n ; i++) Len [i] = 0 ;

    for (k = 0 ; k < n ; k++)
    {
        AMDL_IDX p1 = Cp [k] ;
        AMDL_IDX p2 = Cp [k + 1] ;

        for (p = p1 ; p < p2 ; )
        {
            AMDL_IDX j = Ci [p] ;
            if (j < k)
            {
                Len [j]++ ;
                Len [k]++ ;
                p++ ;
            }
            else if (j == k)
            {
                p++ ;
                nzdiag++ ;
                break ;
            }
            else
            {
                break ;
            }

            /* scan the lower part of column j up to row k */
            {
                AMDL_IDX pj2 = Cp [j + 1] ;
                for (pj = Tp [j] ; pj < pj2 ; )
                {
                    AMDL_IDX t = Ci [pj] ;
                    if (t < k)
                    {
                        Len [t]++ ;
                        Len [j]++ ;
                        pj++ ;
                    }
                    else if (t == k)
                    {
                        pj++ ;
                        nzboth++ ;
                        break ;
                    }
                    else
                    {
                        break ;
                    }
                }
                Tp [j] = pj ;
            }
        }
        Tp [k] = p ;
    }

    /* remaining lower-triangular entries */
    for (j = 0 ; j < n ; j++)
    {
        for (pj = Tp [j] ; pj < Cp [j + 1] ; pj++)
        {
            AMDL_IDX t = Ci [pj] ;
            Len [t]++ ;
            Len [j]++ ;
        }
    }

    /* --------------------------------------------------------------------- */
    /* prefix the union degrees and fill A+A' in O(nnz)                      */
    /* --------------------------------------------------------------------- */

    rp = AMDL_FN(xmalloc) ((size_t) n + 1, sizeof (AMDL_IDX)) ;
    if (rp == NULL) goto done ;

    sum = 0 ;
    for (i = 0 ; i < n ; i++)
    {
        uint64_t c = (uint64_t) Len [i] ;
        rp [i] = (AMDL_IDX) sum ;
        sum += c ;
        if (sum > (uint64_t) AMDL_IDX_MAX) goto done ;
    }
    rp [n] = (AMDL_IDX) sum ;

    if (sum > 0)
    {
        rj = AMDL_FN(xmalloc) ((size_t) sum, sizeof (AMDL_IDX)) ;
        if (rj == NULL) goto done ;
    }

    /* Len is now the per-row write cursor */
    for (i = 0 ; i < n ; i++) Len [i] = rp [i] ;

    for (k = 0 ; k < n ; k++)
    {
        AMDL_IDX p1 = Cp [k] ;
        AMDL_IDX p2 = Cp [k + 1] ;

        for (p = p1 ; p < p2 ; )
        {
            AMDL_IDX j = Ci [p] ;
            if (j < k)
            {
                rj [Len [j]++] = k ;
                rj [Len [k]++] = j ;
                p++ ;
            }
            else if (j == k)
            {
                p++ ;
                break ;
            }
            else
            {
                break ;
            }

            {
                AMDL_IDX pj2 = Cp [j + 1] ;
                for (pj = Tp [j] ; pj < pj2 ; )
                {
                    AMDL_IDX t = Ci [pj] ;
                    if (t < k)
                    {
                        rj [Len [t]++] = j ;
                        rj [Len [j]++] = t ;
                        pj++ ;
                    }
                    else if (t == k)
                    {
                        pj++ ;
                        break ;
                    }
                    else
                    {
                        break ;
                    }
                }
                Tp [j] = pj ;
            }
        }
        Tp [k] = p ;
    }

    for (j = 0 ; j < n ; j++)
    {
        for (pj = Tp [j] ; pj < Cp [j + 1] ; pj++)
        {
            AMDL_IDX t = Ci [pj] ;
            rj [Len [t]++] = j ;
            rj [Len [j]++] = t ;
        }
    }

    /* --------------------------------------------------------------------- */
    /* symmetry of the normalized pattern                                    */
    /* --------------------------------------------------------------------- */

    if (need_symmetry)
    {
        int64_t total = (int64_t) nz_norm - nzdiag ;
        if (total > 0)
        {
            symmetry = (2.0 * (double) nzboth) / (double) total ;
            if (symmetry > 1.0) symmetry = 1.0 ;
            if (symmetry < 0.0) symmetry = 0.0 ;
        }
    }

    *out_rp = rp ;
    *out_rj = rj ;
    *out_nzaat = (int64_t) sum ;
    *out_symmetry = symmetry ;
    free (Tp) ;
    free (Len) ;
    free (Cp_own) ;
    free (Ci_own) ;
    return (AMD_LITE_OK) ;

done:
    free (Tp) ;
    free (Len) ;
    free (Cp_own) ;
    free (Ci_own) ;
    free (Clen) ;
    free (rp) ;
    free (rj) ;
    return (AMD_LITE_OUT_OF_MEMORY) ;
}

/* ------------------------------------------------------------------------- */
/* flag handling for the W workspace                                         */
/* ------------------------------------------------------------------------- */

static AMDL_IDX AMDL_FN(clear_flag)
(
    AMDL_IDX wflg,
    AMDL_IDX wbig,
    AMDL_IDX *W,
    AMDL_IDX n
)
{
    AMDL_IDX x ;
    if (wflg < 2 || wflg >= wbig)
    {
        for (x = 0 ; x < n ; x++)
        {
            if (W [x] != 0) W [x] = 1 ;
        }
        wflg = 2 ;
    }
    return (wflg) ;
}

/* ------------------------------------------------------------------------- */
/* post-ordering of the assembly tree                                        */
/* ------------------------------------------------------------------------- */

static AMDL_IDX AMDL_FN(post_tree)
(
    AMDL_IDX root,
    AMDL_IDX k,
    AMDL_IDX *Child,
    const AMDL_IDX *Sibling,
    AMDL_IDX *Order,
    AMDL_IDX *Stack
)
{
    AMDL_IDX f, head, h, i ;
    head = 0 ;
    Stack [0] = root ;

    while (head >= 0)
    {
        i = Stack [head] ;
        if (Child [i] != AMDL_EMPTY)
        {
            /* push children so that the first child is popped first */
            for (f = Child [i] ; f != AMDL_EMPTY ; f = Sibling [f]) head++ ;
            h = head ;
            for (f = Child [i] ; f != AMDL_EMPTY ; f = Sibling [f])
            {
                Stack [h--] = f ;
            }
            Child [i] = AMDL_EMPTY ;
        }
        else
        {
            head-- ;
            Order [i] = k++ ;
        }
    }
    return (k) ;
}

static void AMDL_FN(postorder)
(
    AMDL_IDX nn,
    AMDL_IDX *Parent,
    AMDL_IDX *Nv,
    AMDL_IDX *Fsize,
    AMDL_IDX *Order,
    AMDL_IDX *Child,
    AMDL_IDX *Sibling,
    AMDL_IDX *Stack
)
{
    AMDL_IDX i, j, k, parent, frsize, f, fprev, maxfrsize, bigfprev,
        bigf, fnext ;

    for (j = 0 ; j < nn ; j++)
    {
        Child [j] = AMDL_EMPTY ;
        Sibling [j] = AMDL_EMPTY ;
    }

    /* build child lists; larger elements tend to end up at the end */
    for (j = nn - 1 ; j >= 0 ; j--)
    {
        if (Nv [j] > 0)
        {
            parent = Parent [j] ;
            if (parent != AMDL_EMPTY)
            {
                Sibling [j] = Child [parent] ;
                Child [parent] = j ;
            }
        }
    }

    /* place the largest child last in each child list */
    for (i = 0 ; i < nn ; i++)
    {
        if (Nv [i] > 0 && Child [i] != AMDL_EMPTY)
        {
            fprev = AMDL_EMPTY ;
            maxfrsize = AMDL_EMPTY ;
            bigfprev = AMDL_EMPTY ;
            bigf = AMDL_EMPTY ;

            for (f = Child [i] ; f != AMDL_EMPTY ; f = Sibling [f])
            {
                frsize = Fsize [f] ;
                if (frsize >= maxfrsize)
                {
                    maxfrsize = frsize ;
                    bigfprev = fprev ;
                    bigf = f ;
                }
                fprev = f ;
            }

            fnext = Sibling [bigf] ;
            if (fnext != AMDL_EMPTY)
            {
                if (bigfprev == AMDL_EMPTY)
                {
                    Child [i] = fnext ;
                }
                else
                {
                    Sibling [bigfprev] = fnext ;
                }
                Sibling [bigf] = AMDL_EMPTY ;
                Sibling [fprev] = bigf ;
            }
        }
    }

    for (i = 0 ; i < nn ; i++) Order [i] = AMDL_EMPTY ;

    k = 0 ;
    for (i = 0 ; i < nn ; i++)
    {
        if (Parent [i] == AMDL_EMPTY && Nv [i] > 0)
        {
            k = AMDL_FN(post_tree) (i, k, Child, Sibling, Order, Stack) ;
        }
    }
}

/* ------------------------------------------------------------------------- */
/* AMD ordering kernel (quotient graph, mass elimination, absorption)        */
/* ------------------------------------------------------------------------- */

static int AMDL_FN(kernel)
(
    AMDL_IDX n,
    AMDL_IDX *Pe,
    AMDL_IDX *Iw,
    AMDL_IDX *Len,
    AMDL_IDX iwlen,
    AMDL_IDX pfree,
    AMDL_IDX *Nv,
    AMDL_IDX *Next,
    AMDL_IDX *Last,
    AMDL_IDX *Head,
    AMDL_IDX *Elen,
    AMDL_IDX *Degree,
    AMDL_IDX *W,
    AMDL_IDX *P,
    double alpha,
    int aggressive,
    amd_lite_info *info
)
{
    AMDL_IDX deg, degme, e, elenme, eln, i, ilast, inext, j, jlast,
        jnext, k, knt1, knt2, knt3, lenj, ln, me, mindeg, nel, nleft,
        nvi, nvj, nvpiv, slenme, wbig, we, wflg, wnvi, ndense, ncmpa,
        dense, lemax, p, p1, p2, p3, p4, pdst, pend, pj, pme, pme1,
        pme2, pn, psrc ;
    AMDL_UINT hash ;

    double f, r, ndiv, s, nms_lu, nms_ldl, dmax, lnz, lnzme ;

    /* --------------------------------------------------------------------- */
    /* initialization                                                        */
    /* --------------------------------------------------------------------- */

    lnz = 0 ;
    ndiv = 0 ;
    nms_lu = 0 ;
    nms_ldl = 0 ;
    dmax = 1 ;
    me = AMDL_EMPTY ;
    mindeg = 0 ;
    ncmpa = 0 ;
    nel = 0 ;
    lemax = 0 ;

    if (!(alpha >= 0))
    {
        /* NaN: use the default; negative: handled below */
        if (alpha != alpha) alpha = AMD_LITE_DEFAULT_DENSE ;
    }
    if (alpha < 0)
    {
        dense = (n > 2) ? n - 2 : n ;
    }
    else
    {
        double d = alpha * sqrt ((double) n) ;
        if (d >= (double) n) dense = n ;
        else if (d <= 16.0) dense = 16 ;
        else dense = (AMDL_IDX) d ;
    }
    if (dense < 16) dense = 16 ;
    if (dense > n) dense = n ;

    for (i = 0 ; i < n ; i++)
    {
        Last [i] = AMDL_EMPTY ;
        Head [i] = AMDL_EMPTY ;
        Next [i] = AMDL_EMPTY ;
        Nv [i] = 1 ;
        W [i] = 1 ;
        Elen [i] = 0 ;
        Degree [i] = Len [i] ;
    }

    wbig = AMDL_IDX_MAX - n ;
    wflg = AMDL_FN(clear_flag) (0, wbig, W, n) ;

    /* --------------------------------------------------------------------- */
    /* initialize degree lists and eliminate empty/dense rows                */
    /* --------------------------------------------------------------------- */

    ndense = 0 ;
    for (i = 0 ; i < n ; i++)
    {
        deg = Degree [i] ;
        if (deg == 0)
        {
            Elen [i] = AMDL_FLIP (1) ;
            nel++ ;
            Pe [i] = AMDL_EMPTY ;
            W [i] = 0 ;
        }
        else if (deg > dense)
        {
            ndense++ ;
            Nv [i] = 0 ;
            Elen [i] = AMDL_EMPTY ;
            nel++ ;
            Pe [i] = AMDL_EMPTY ;
        }
        else
        {
            inext = Head [deg] ;
            if (inext != AMDL_EMPTY) Last [inext] = i ;
            Next [i] = inext ;
            Head [deg] = i ;
        }
    }

    /* --------------------------------------------------------------------- */
    /* main pivot loop                                                       */
    /* --------------------------------------------------------------------- */

    while (nel < n)
    {
        /* ----------------------------------------------------------------- */
        /* find the next supervariable of minimum degree                     */
        /* ----------------------------------------------------------------- */

        for (deg = mindeg ; deg < n ; deg++)
        {
            me = Head [deg] ;
            if (me != AMDL_EMPTY) break ;
        }
        if (deg >= n || me == AMDL_EMPTY) return (AMD_LITE_INVALID) ;
        mindeg = deg ;

        inext = Next [me] ;
        if (inext != AMDL_EMPTY) Last [inext] = AMDL_EMPTY ;
        Head [deg] = inext ;

        elenme = Elen [me] ;
        nvpiv = Nv [me] ;
        nel += nvpiv ;

        /* ----------------------------------------------------------------- */
        /* construct the new element                                         */
        /* ----------------------------------------------------------------- */

        Nv [me] = -nvpiv ;
        degme = 0 ;

        if (elenme == 0)
        {
            /* construct the new element in place */
            pme1 = Pe [me] ;
            pme2 = pme1 - 1 ;

            for (p = pme1 ; p <= pme1 + Len [me] - 1 ; p++)
            {
                i = Iw [p] ;
                nvi = Nv [i] ;
                if (nvi > 0)
                {
                    degme += nvi ;
                    Nv [i] = -nvi ;
                    Iw [++pme2] = i ;

                    /* remove i from its degree list */
                    ilast = Last [i] ;
                    inext = Next [i] ;
                    if (inext != AMDL_EMPTY) Last [inext] = ilast ;
                    if (ilast != AMDL_EMPTY)
                    {
                        Next [ilast] = inext ;
                    }
                    else
                    {
                        Head [Degree [i]] = inext ;
                    }
                }
            }
        }
        else
        {
            /* construct the new element in empty space */
            p = Pe [me] ;
            pme1 = pfree ;
            slenme = Len [me] - elenme ;

            for (knt1 = 1 ; knt1 <= elenme + 1 ; knt1++)
            {
                if (knt1 > elenme)
                {
                    e = me ;
                    pj = p ;
                    ln = slenme ;
                }
                else
                {
                    e = Iw [p++] ;
                    pj = Pe [e] ;
                    ln = Len [e] ;
                }

                for (knt2 = 1 ; knt2 <= ln ; knt2++)
                {
                    i = Iw [pj++] ;
                    nvi = Nv [i] ;

                    if (nvi > 0)
                    {
                        /* ------------------------------------------------- */
                        /* garbage collection, if needed                     */
                        /* ------------------------------------------------- */
                        if (pfree >= iwlen)
                        {
                            Pe [me] = p ;
                            Len [me] -= knt1 ;
                            if (Len [me] == 0) Pe [me] = AMDL_EMPTY ;
                            Pe [e] = pj ;
                            Len [e] = ln - knt2 ;
                            if (Len [e] == 0) Pe [e] = AMDL_EMPTY ;

                            ncmpa++ ;

                            /* store first entry of each object in Pe and */
                            /* flip it as a marker                       */
                            for (j = 0 ; j < n ; j++)
                            {
                                pn = Pe [j] ;
                                if (pn >= 0)
                                {
                                    Pe [j] = Iw [pn] ;
                                    Iw [pn] = AMDL_FLIP (j) ;
                                }
                            }

                            psrc = 0 ;
                            pdst = 0 ;
                            pend = pme1 - 1 ;

                            while (psrc <= pend)
                            {
                                j = AMDL_FLIP (Iw [psrc++]) ;
                                if (j >= 0)
                                {
                                    Iw [pdst] = Pe [j] ;
                                    Pe [j] = pdst++ ;
                                    lenj = Len [j] ;
                                    for (knt3 = 0 ; knt3 <= lenj - 2 ; knt3++)
                                    {
                                        Iw [pdst++] = Iw [psrc++] ;
                                    }
                                }
                            }

                            /* move the partially constructed element */
                            p1 = pdst ;
                            for (psrc = pme1 ; psrc <= pfree - 1 ; psrc++)
                            {
                                Iw [pdst++] = Iw [psrc] ;
                            }
                            pme1 = p1 ;
                            pfree = pdst ;
                            pj = Pe [e] ;
                            p = Pe [me] ;
                        }

                        /* ------------------------------------------------- */
                        /* add i to the new element                          */
                        /* ------------------------------------------------- */
                        degme += nvi ;
                        Nv [i] = -nvi ;
                        Iw [pfree++] = i ;

                        /* remove i from its degree list */
                        ilast = Last [i] ;
                        inext = Next [i] ;
                        if (inext != AMDL_EMPTY) Last [inext] = ilast ;
                        if (ilast != AMDL_EMPTY)
                        {
                            Next [ilast] = inext ;
                        }
                        else
                        {
                            Head [Degree [i]] = inext ;
                        }
                    }
                }

                if (e != me)
                {
                    /* absorb element e into the new element me */
                    Pe [e] = AMDL_FLIP (me) ;
                    W [e] = 0 ;
                }
            }

            pme2 = pfree - 1 ;
        }

        /* me is now an element in Iw [pme1..pme2] */
        Degree [me] = degme ;
        Pe [me] = pme1 ;
        Len [me] = pme2 - pme1 + 1 ;
        Elen [me] = AMDL_FLIP (nvpiv + degme) ;

        wflg = AMDL_FN(clear_flag) (wflg, wbig, W, n) ;

        /* ----------------------------------------------------------------- */
        /* scan 1: |Le \ Lme| for all elements in Lme                        */
        /* ----------------------------------------------------------------- */

        for (pme = pme1 ; pme <= pme2 ; pme++)
        {
            i = Iw [pme] ;
            eln = Elen [i] ;
            if (eln > 0)
            {
                nvi = -Nv [i] ;
                wnvi = wflg - nvi ;
                for (p = Pe [i] ; p <= Pe [i] + eln - 1 ; p++)
                {
                    e = Iw [p] ;
                    we = W [e] ;
                    if (we >= wflg)
                    {
                        we -= nvi ;
                    }
                    else if (we != 0)
                    {
                        we = Degree [e] + wnvi ;
                    }
                    W [e] = we ;
                }
            }
        }

        /* ----------------------------------------------------------------- */
        /* scan 2: degree update, element absorption, and hash buckets       */
        /* ----------------------------------------------------------------- */

        for (pme = pme1 ; pme <= pme2 ; pme++)
        {
            i = Iw [pme] ;
            p1 = Pe [i] ;
            p2 = p1 + Elen [i] - 1 ;
            pn = p1 ;
            hash = 0 ;
            deg = 0 ;

            /* element list of i */
            if (aggressive)
            {
                for (p = p1 ; p <= p2 ; p++)
                {
                    e = Iw [p] ;
                    we = W [e] ;
                    if (we != 0)
                    {
                        AMDL_IDX dext = we - wflg ;
                        if (dext > 0)
                        {
                            deg += dext ;
                            Iw [pn++] = e ;
                            hash += (AMDL_UINT) e ;
                        }
                        else
                        {
                            /* external degree is zero: absorb e into me */
                            Pe [e] = AMDL_FLIP (me) ;
                            W [e] = 0 ;
                        }
                    }
                }
            }
            else
            {
                for (p = p1 ; p <= p2 ; p++)
                {
                    e = Iw [p] ;
                    we = W [e] ;
                    if (we != 0)
                    {
                        AMDL_IDX dext = we - wflg ;
                        deg += dext ;
                        Iw [pn++] = e ;
                        hash += (AMDL_UINT) e ;
                    }
                }
            }

            Elen [i] = pn - p1 + 1 ;

            /* variable list of i */
            p3 = pn ;
            p4 = p1 + Len [i] ;
            for (p = p2 + 1 ; p < p4 ; p++)
            {
                j = Iw [p] ;
                nvj = Nv [j] ;
                if (nvj > 0)
                {
                    deg += nvj ;
                    Iw [pn++] = j ;
                    hash += (AMDL_UINT) j ;
                }
            }

            if (Elen [i] == 1 && p3 == pn)
            {
                /* --------------------------------------------------------- */
                /* mass elimination: i is only adjacent to me                */
                /* --------------------------------------------------------- */
                Pe [i] = AMDL_FLIP (me) ;
                nvi = -Nv [i] ;
                degme -= nvi ;
                nvpiv += nvi ;
                nel += nvi ;
                Nv [i] = 0 ;
                Elen [i] = AMDL_EMPTY ;
            }
            else
            {
                /* update the approximate degree */
                Degree [i] = (Degree [i] < deg) ? Degree [i] : deg ;

                /* move first supervariable to end of list, first element to */
                /* end of element part, and add me to the front            */
                Iw [pn] = Iw [p3] ;
                Iw [p3] = Iw [p1] ;
                Iw [p1] = me ;
                Len [i] = pn - p1 + 1 ;

                hash = hash % (AMDL_UINT) n ;
                j = Head [hash] ;
                if (j <= AMDL_EMPTY)
                {
                    Next [i] = AMDL_FLIP (j) ;
                    Head [hash] = AMDL_FLIP (i) ;
                }
                else
                {
                    Next [i] = Last [j] ;
                    Last [j] = i ;
                }
                Last [i] = (AMDL_IDX) hash ;
            }
        }

        Degree [me] = degme ;

        /* clear the flag array; adding the largest element size seen so
         * far keeps all existing W values strictly below wflg */
        if (degme > lemax) lemax = degme ;
        wflg += lemax ;
        wflg = AMDL_FN(clear_flag) (wflg, wbig, W, n) ;

        /* ----------------------------------------------------------------- */
        /* supervariable detection                                           */
        /* ----------------------------------------------------------------- */

        for (pme = pme1 ; pme <= pme2 ; pme++)
        {
            i = Iw [pme] ;
            if (Nv [i] < 0)
            {
                /* let i be the head of a hash bucket and empty the bucket */
                hash = (AMDL_UINT) Last [i] ;
                j = Head [hash] ;
                if (j == AMDL_EMPTY)
                {
                    i = AMDL_EMPTY ;
                }
                else if (j < AMDL_EMPTY)
                {
                    i = AMDL_FLIP (j) ;
                    Head [hash] = AMDL_EMPTY ;
                }
                else
                {
                    i = Last [j] ;
                    Last [j] = AMDL_EMPTY ;
                }

                while (i != AMDL_EMPTY && Next [i] != AMDL_EMPTY)
                {
                    ln = Len [i] ;
                    eln = Elen [i] ;

                    /* scatter i into W, skipping the first element (me) */
                    for (p = Pe [i] + 1 ; p <= Pe [i] + ln - 1 ; p++)
                    {
                        W [Iw [p]] = wflg ;
                    }

                    jlast = i ;
                    j = Next [i] ;
                    while (j != AMDL_EMPTY)
                    {
                        int same = (Len [j] == ln) && (Elen [j] == eln) ;
                        for (p = Pe [j] + 1 ; same && p <= Pe [j] + ln - 1 ; p++)
                        {
                            if (W [Iw [p]] != wflg) same = 0 ;
                        }
                        if (same)
                        {
                            /* j has the same pattern as i: merge j into i */
                            Pe [j] = AMDL_FLIP (i) ;
                            Nv [i] += Nv [j] ;
                            Nv [j] = 0 ;
                            Elen [j] = AMDL_EMPTY ;
                            j = Next [j] ;
                            Next [jlast] = j ;
                        }
                        else
                        {
                            jlast = j ;
                            j = Next [j] ;
                        }
                    }
                    wflg++ ;
                    i = Next [i] ;
                }
            }
        }

        /* ----------------------------------------------------------------- */
        /* restore degree lists and remove nonprincipal supervariables       */
        /* ----------------------------------------------------------------- */

        p = pme1 ;
        nleft = n - nel ;
        for (pme = pme1 ; pme <= pme2 ; pme++)
        {
            i = Iw [pme] ;
            nvi = -Nv [i] ;
            if (nvi > 0)
            {
                Nv [i] = nvi ;

                deg = Degree [i] + degme - nvi ;
                {
                    AMDL_IDX limit = nleft - nvi ;
                    if (deg > limit) deg = limit ;
                }

                inext = Head [deg] ;
                if (inext != AMDL_EMPTY) Last [inext] = i ;
                Next [i] = inext ;
                Last [i] = AMDL_EMPTY ;
                Head [deg] = i ;

                if (deg < mindeg) mindeg = deg ;
                Degree [i] = deg ;

                Iw [p++] = i ;
            }
        }

        /* ----------------------------------------------------------------- */
        /* finalize the new element                                          */
        /* ----------------------------------------------------------------- */

        Nv [me] = nvpiv ;
        Len [me] = p - pme1 ;
        if (Len [me] == 0)
        {
            Pe [me] = AMDL_EMPTY ;
            W [me] = 0 ;
        }
        if (elenme != 0)
        {
            pfree = p ;
        }

        if (info != NULL)
        {
            f = (double) nvpiv ;
            r = (double) degme + (double) ndense ;
            if (f + r > dmax) dmax = f + r ;

            lnzme = f * r + (f - 1) * f / 2 ;
            lnz += lnzme ;
            ndiv += lnzme ;

            s = f * r * r + r * (f - 1) * f + (f - 1) * f * (2 * f - 1) / 6 ;
            nms_lu += s ;
            nms_ldl += (s + lnzme) / 2 ;
        }
    }

    /* --------------------------------------------------------------------- */
    /* final statistics                                                      */
    /* --------------------------------------------------------------------- */

    if (info != NULL)
    {
        f = (double) ndense ;
        if (f > dmax) dmax = f ;

        lnzme = (f - 1) * f / 2 ;
        lnz += lnzme ;
        ndiv += lnzme ;

        s = (f - 1) * f * (2 * f - 1) / 6 ;
        nms_lu += s ;
        nms_ldl += (s + lnzme) / 2 ;

        info->ndense = (int64_t) ndense ;
        info->ncmpa = (int64_t) ncmpa ;
        info->dmax = (int64_t) dmax ;
        info->lnz = lnz ;
        info->ndiv = ndiv ;
        info->nmultsubs_ldl = nms_ldl ;
        info->nmultsubs_lu = nms_lu ;
    }

    /* --------------------------------------------------------------------- */
    /* restore the assembly tree and post-order it                           */
    /* --------------------------------------------------------------------- */

    for (i = 0 ; i < n ; i++) Pe [i] = AMDL_FLIP (Pe [i]) ;
    for (i = 0 ; i < n ; i++) Elen [i] = AMDL_FLIP (Elen [i]) ;

    /* compress paths of variables that were merged into elements */
    for (i = 0 ; i < n ; i++)
    {
        if (Nv [i] == 0)
        {
            j = Pe [i] ;
            if (j == AMDL_EMPTY) continue ;
            while (Nv [j] == 0) j = Pe [j] ;
            e = j ;
            j = i ;
            while (Nv [j] == 0)
            {
                jnext = Pe [j] ;
                Pe [j] = e ;
                j = jnext ;
            }
        }
    }

    /* postorder the assembly tree: W receives the element order */
    AMDL_FN(postorder) (n, Pe, Nv, Elen, W, Head, Next, Last) ;

    /* --------------------------------------------------------------------- */
    /* build the output permutation                                          */
    /* --------------------------------------------------------------------- */

    for (k = 0 ; k < n ; k++)
    {
        Head [k] = AMDL_EMPTY ;
        Next [k] = AMDL_EMPTY ;
    }
    for (e = 0 ; e < n ; e++)
    {
        k = W [e] ;
        if (k != AMDL_EMPTY) Head [k] = e ;
    }

    nel = 0 ;
    for (k = 0 ; k < n ; k++)
    {
        e = Head [k] ;
        if (e == AMDL_EMPTY) break ;
        Next [e] = nel ;
        nel += Nv [e] ;
    }

    /* nonprincipal variables go just before their element */
    for (i = 0 ; i < n ; i++)
    {
        if (Nv [i] == 0)
        {
            e = Pe [i] ;
            if (e != AMDL_EMPTY)
            {
                Next [i] = Next [e] ;
                Next [e]++ ;
            }
            else
            {
                Next [i] = nel++ ;
            }
        }
    }

    for (i = 0 ; i < n ; i++)
    {
        k = Next [i] ;
        if (k < 0 || k >= n) return (AMD_LITE_INVALID) ;
        P [k] = i ;
    }

    return (AMD_LITE_OK) ;
}

/* ------------------------------------------------------------------------- */
/* public entry point                                                        */
/* ------------------------------------------------------------------------- */

int AMDL_PUBLIC
(
    AMDL_IDX n,
    const AMDL_IDX *rowptr,
    const AMDL_IDX *colidx,
    AMDL_IDX *perm,
    const amd_lite_options *options,
    amd_lite_info *info
)
{
    amd_lite_options opt ;
    AMDL_IDX *Rp = NULL, *Rj = NULL, *work = NULL ;
    AMDL_IDX *Pe, *Len, *Nv, *Next, *Last, *Head, *Elen, *Degree, *W, *Iw ;
    int64_t nz = 0, nzdiag = 0, nzaat = 0 ;
    int clean_input = 0 ;
    double symmetry = 1.0 ;
    uint64_t nzaat_u, iwlen_u, work_u, bytes_u ;
    uint64_t nn ;
    int status ;

    if (info != NULL)
    {
        memset (info, 0, sizeof (amd_lite_info)) ;
        info->status = AMD_LITE_OK ;
        info->symmetry = 1.0 ;
    }

    if (options != NULL)
    {
        opt = *options ;
    }
    else
    {
        amd_lite_defaults (&opt) ;
    }
    if (opt.format != AMD_LITE_CSR && opt.format != AMD_LITE_CSC)
    {
        if (info != NULL) info->status = AMD_LITE_INVALID ;
        return (AMD_LITE_INVALID) ;
    }

    status = AMDL_FN(validate) (n, rowptr, colidx, &nz, &nzdiag,
        &clean_input) ;
    if (status != AMD_LITE_OK)
    {
        if (info != NULL) info->status = status ;
        return (status) ;
    }

    if (info != NULL)
    {
        info->n = (int64_t) n ;
        info->nz = nz ;
        info->nzdiag = nzdiag ;
    }

    if (n == 0) return (AMD_LITE_OK) ;
    if (perm == NULL)
    {
        if (info != NULL) info->status = AMD_LITE_INVALID ;
        return (AMD_LITE_INVALID) ;
    }
    if ((uint64_t) n > (uint64_t) (AMDL_IDX_MAX / 16))
    {
        if (info != NULL) info->status = AMD_LITE_OUT_OF_MEMORY ;
        return (AMD_LITE_OUT_OF_MEMORY) ;
    }

    status = AMDL_FN(build_pattern) (n, rowptr, colidx, nz, clean_input,
        &Rp, &Rj, &nzaat, &symmetry, info != NULL) ;
    if (status != AMD_LITE_OK)
    {
        if (info != NULL) info->status = status ;
        return (status) ;
    }

    nzaat_u = (uint64_t) nzaat ;
    nn = (uint64_t) n ;

    /* workspace: 9 size-n vectors plus the matrix and elbow room */
    iwlen_u = nzaat_u + nzaat_u / 5 + nn ;
    if (iwlen_u < nzaat_u + nn || iwlen_u > (uint64_t) AMDL_IDX_MAX)
    {
        free (Rp) ;
        free (Rj) ;
        if (info != NULL) info->status = AMD_LITE_OUT_OF_MEMORY ;
        return (AMD_LITE_OUT_OF_MEMORY) ;
    }
    work_u = 9 * nn + iwlen_u ;
    if (work_u < 9 * nn || work_u > SIZE_MAX / sizeof (AMDL_IDX))
    {
        free (Rp) ;
        free (Rj) ;
        if (info != NULL) info->status = AMD_LITE_OUT_OF_MEMORY ;
        return (AMD_LITE_OUT_OF_MEMORY) ;
    }

    work = AMDL_FN(xmalloc) ((size_t) work_u, sizeof (AMDL_IDX)) ;
    if (work == NULL)
    {
        free (Rp) ;
        free (Rj) ;
        if (info != NULL) info->status = AMD_LITE_OUT_OF_MEMORY ;
        return (AMD_LITE_OUT_OF_MEMORY) ;
    }

    Pe     = work ;
    Len    = Pe     + n ;
    Nv     = Len    + n ;
    Next   = Nv     + n ;
    Last   = Next   + n ;
    Head   = Last   + n ;
    Elen   = Head   + n ;
    Degree = Elen   + n ;
    W      = Degree + n ;
    Iw     = W      + n ;

    /* copy the symmetric pattern into Iw */
    {
        AMDL_IDX pfree = 0, i ;
        for (i = 0 ; i < n ; i++)
        {
            AMDL_IDX len = Rp [i + 1] - Rp [i] ;
            Pe [i] = pfree ;
            Len [i] = len ;
            if (len > 0)
            {
                memcpy (Iw + pfree, Rj + Rp [i],
                    (size_t) len * sizeof (AMDL_IDX)) ;
                pfree += len ;
            }
        }
        free (Rp) ;
        Rp = NULL ;
        free (Rj) ;
        Rj = NULL ;

        status = AMDL_FN(kernel) (n, Pe, Iw, Len, (AMDL_IDX) iwlen_u,
            pfree, Nv, Next, Last, Head, Elen, Degree, W, perm,
            opt.dense, opt.aggressive != 0, info) ;
    }

    if (info != NULL)
    {
        info->nz_a_plus_at = nzaat ;
        info->symmetry = symmetry ;
        info->status = status ;
        bytes_u = work_u * sizeof (AMDL_IDX) ;
        info->memory = (double) bytes_u ;
    }

    free (work) ;
    return (status) ;
}

#undef AMDL_EMPTY
#undef AMDL_FLIP
