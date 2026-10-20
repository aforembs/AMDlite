/* amd_lite_demo.c: order a 2-D 5-point grid and print statistics.
 *
 * usage: amd_lite_demo [grid_size]
 */

#include "amd_lite.h"

#include <stdio.h>
#include <stdlib.h>

int main (int argc, char **argv)
{
    int m = (argc > 1) ? atoi (argv [1]) : 20 ;
    int32_t n, i, j, nz = 0 ;
    int32_t *rp, *ci, *p ;
    amd_lite_info info ;
    int status ;

    if (m < 1) m = 1 ;
    n = (int32_t) m * (int32_t) m ;
    rp = (int32_t *) malloc (((size_t) n + 1) * sizeof (int32_t)) ;
    ci = (int32_t *) malloc ((size_t) n * 4 * sizeof (int32_t)) ;
    p  = (int32_t *) malloc ((size_t) n * sizeof (int32_t)) ;
    if (rp == NULL || ci == NULL || p == NULL)
    {
        fprintf (stderr, "out of memory\n") ;
        free (rp) ;
        free (ci) ;
        free (p) ;
        return (1) ;
    }

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

    status = amd_lite_order (n, rp, ci, p, NULL, &info) ;
    if (status != AMD_LITE_OK)
    {
        fprintf (stderr, "amd_lite_order failed: %d\n", status) ;
        free (rp) ;
        free (ci) ;
        free (p) ;
        return (1) ;
    }

    printf ("AMD-lite demo: %d-by-%d grid (%d rows, %d entries)\n",
        m, m, (int) n, (int) nz) ;
    printf ("status          %d\n", info.status) ;
    printf ("symmetry        %.3f\n", info.symmetry) ;
    printf ("nz(A+A')        %lld\n", (long long) info.nz_a_plus_at) ;
    printf ("dense rows      %lld\n", (long long) info.ndense) ;
    printf ("garbage coll.   %lld\n", (long long) info.ncmpa) ;
    printf ("approx nnz(L)   %.0f\n", info.lnz) ;
    printf ("divides         %.0f\n", info.ndiv) ;
    printf ("mult-sub (LDL') %.0f\n", info.nmultsubs_ldl) ;
    printf ("mult-sub (LU)   %.0f\n", info.nmultsubs_lu) ;
    printf ("max L column    %lld\n", (long long) info.dmax) ;
    printf ("memory          %.0f bytes\n", info.memory) ;

    if (n <= 100)
    {
        printf ("permutation:") ;
        for (i = 0 ; i < n ; i++) printf (" %d", (int) p [i]) ;
        printf ("\n") ;
    }

    free (rp) ;
    free (ci) ;
    free (p) ;
    return (0) ;
}
