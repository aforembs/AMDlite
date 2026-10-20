#define _POSIX_C_SOURCE 199309L
#include "amd_lite.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

/* SuiteSparse AMD, int32, CSC input; for symmetric matrices CSR == CSC */
extern int amd_order(int32_t n, const int32_t *Ap, const int32_t *Ai,
                     int32_t *P, double *Control, double *Info);

typedef struct { int32_t n; int32_t *rp, *ci; int32_t nz; } Matrix;

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static int cmp_i32(const void *a, const void *b)
{
    int32_t x = *(const int32_t *) a, y = *(const int32_t *) b;
    return (x > y) - (x < y);
}

static Matrix mat_alloc(int32_t n, int32_t cap)
{
    Matrix M;
    M.n = n;
    M.rp = (int32_t *) malloc(((size_t) n + 1) * sizeof(int32_t));
    M.ci = (int32_t *) malloc((size_t) (cap ? cap : 1) * sizeof(int32_t));
    M.nz = 0;
    M.rp[0] = 0;
    return M;
}

static void mat_free(Matrix *M)
{
    free(M->rp);
    free(M->ci);
    M->rp = M->ci = NULL;
}

static void sort_rows(Matrix *M)
{
    for (int32_t i = 0; i < M->n; i++)
    {
        int32_t p = M->rp[i], q = M->rp[i + 1], k = 0;
        if (q - p > 1) qsort(M->ci + p, (size_t) (q - p), sizeof(int32_t), cmp_i32);
        for (int32_t r = p; r < q; r++)
            if (k == 0 || M->ci[r] != M->ci[p + k - 1]) M->ci[p + k++] = M->ci[r];
        M->rp[i + 1] = p + k;
    }
    M->nz = M->rp[M->n];
}

static Matrix make_grid2d(int m)
{
    int32_t n = (int32_t) m * m;
    Matrix M = mat_alloc(n, n * 4);
    int32_t nz = 0;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            int32_t u = (int32_t) i * m + j;
            M.rp[u] = nz;
            if (i > 0)     M.ci[nz++] = (int32_t) (i - 1) * m + j;
            if (j > 0)     M.ci[nz++] = u - 1;
            if (j < m - 1) M.ci[nz++] = u + 1;
            if (i < m - 1) M.ci[nz++] = (int32_t) (i + 1) * m + j;
        }
    }
    M.rp[n] = nz;
    M.nz = nz;
    return M;
}

static Matrix make_grid3d(int d)
{
    int32_t n = (int32_t) d * d * d;
    Matrix M = mat_alloc(n, n * 6);
    int32_t nz = 0;
    for (int z = 0; z < d; z++)
    for (int y = 0; y < d; y++)
    for (int x = 0; x < d; x++)
    {
        int32_t u = (int32_t) (z * d + y) * d + x;
        int32_t nb[6], k = 0;
        M.rp[u] = nz;
        if (x > 0) nb[k++] = u - 1;
        if (x < d - 1) nb[k++] = u + 1;
        if (y > 0) nb[k++] = u - d;
        if (y < d - 1) nb[k++] = u + d;
        if (z > 0) nb[k++] = u - (int32_t) d * d;
        if (z < d - 1) nb[k++] = u + (int32_t) d * d;
        for (int a = 1; a < k; a++) { int32_t v = nb[a], b = a; while (b > 0 && nb[b-1] > v) { nb[b] = nb[b-1]; b--; } nb[b] = v; }
        for (int a = 0; a < k; a++) M.ci[nz++] = nb[a];
    }
    M.rp[n] = nz;
    M.nz = nz;
    return M;
}

static Matrix make_banded(int32_t n, int b)
{
    Matrix M = mat_alloc(n, (int32_t) n * (2 * b));
    int32_t nz = 0;
    for (int32_t i = 0; i < n; i++)
    {
        M.rp[i] = nz;
        for (int32_t j = i - b; j <= i + b; j++)
            if (j >= 0 && j < n && j != i) M.ci[nz++] = j;
    }
    M.rp[n] = nz;
    M.nz = nz;
    return M;
}

static Matrix make_star(int32_t n)
{
    Matrix M = mat_alloc(n, 2 * (n - 1));
    int32_t nz = 0;
    for (int32_t i = 0; i < n; i++)
    {
        M.rp[i] = nz;
        if (i == 0) { for (int32_t j = 1; j < n; j++) M.ci[nz++] = j; }
        else M.ci[nz++] = 0;
    }
    M.rp[n] = nz;
    M.nz = nz;
    return M;
}

static uint64_t rng_state = 88172645463325252ULL;
static uint32_t rng(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (uint32_t) (rng_state >> 11);
}

static Matrix make_random_sym(int32_t n, int avgdeg)
{
    int32_t *deg = (int32_t *) calloc((size_t) n, sizeof(int32_t));
    int64_t e = (int64_t) n * avgdeg / 2;
    int32_t *eu = (int32_t *) malloc((size_t) e * sizeof(int32_t));
    int32_t *ev = (int32_t *) malloc((size_t) e * sizeof(int32_t));
    int64_t k = 0;
    for (int64_t t = 0; t < e; t++)
    {
        int32_t u = (int32_t) (rng() % (uint32_t) n);
        int32_t v = (int32_t) (rng() % (uint32_t) n);
        if (u == v) { t--; continue; }
        eu[k] = u; ev[k] = v; k++;
        deg[u]++; deg[v]++;
    }
    int32_t cap = 0;
    for (int32_t i = 0; i < n; i++) cap += deg[i];
    Matrix M = mat_alloc(n, cap);
    M.rp[0] = 0;
    for (int32_t i = 0; i < n; i++) M.rp[i + 1] = M.rp[i] + deg[i];
    int32_t *cur = (int32_t *) malloc((size_t) n * sizeof(int32_t));
    for (int32_t i = 0; i < n; i++) cur[i] = M.rp[i];
    for (int64_t t = 0; t < k; t++)
    {
        M.ci[cur[eu[t]]++] = ev[t];
        M.ci[cur[ev[t]]++] = eu[t];
    }
    free(cur); free(deg); free(eu); free(ev);
    sort_rows(&M);
    return M;
}

static void reverse_rows(Matrix *M)
{
    for (int32_t i = 0; i < M->n; i++)
    {
        int32_t p = M->rp[i], q = M->rp[i + 1] - 1;
        while (p < q) { int32_t t = M->ci[p]; M->ci[p++] = M->ci[q]; M->ci[q--] = t; }
    }
}

static int check_perm(int32_t n, const int32_t *p)
{
    unsigned char *seen = (unsigned char *) calloc((size_t) (n ? n : 1), 1);
    int ok = 1;
    for (int32_t i = 0; i < n; i++)
    {
        if (p[i] < 0 || p[i] >= n || seen[p[i]]) { ok = 0; break; }
        seen[p[i]] = 1;
    }
    free(seen);
    return ok;
}

static void run(const char *name, Matrix *M, int reps)
{
    int32_t n = M->n;
    int32_t *p1 = (int32_t *) malloc((size_t) n * sizeof(int32_t));
    int32_t *p2 = (int32_t *) malloc((size_t) n * sizeof(int32_t));
    amd_lite_options opt;
    amd_lite_info info;
    double control[5] = { 10.0, 1.0, 0.0, 0.0, 0.0 };
    double oinfo[20];
    double best1 = 1e300, best2 = 1e300;
    int st1 = 0, st2 = 0;

    amd_lite_defaults(&opt);
    opt.dense = 10.0;
    opt.aggressive = 1;

    /* warm-up */
    st1 = amd_lite_order(n, M->rp, M->ci, p1, &opt, &info);
    st2 = amd_order(n, M->rp, M->ci, p2, control, oinfo);

    for (int r = 0; r < reps; r++)
    {
        double t;
        t = now(); amd_lite_order(n, M->rp, M->ci, p1, &opt, &info); t = now() - t;
        if (t < best1) best1 = t;
        t = now(); amd_order(n, M->rp, M->ci, p2, control, oinfo); t = now() - t;
        if (t < best2) best2 = t;
    }

    int ok1 = check_perm(n, p1);
    int ok2 = check_perm(n, p2);

    /* theoretical peak index counts for the current implementation */
    uint64_t nzaat = (uint64_t) info.nz_a_plus_at;
    int clean = 1;
    for (int32_t ii = 0; ii < n && clean; ii++)
        for (int32_t pp = M->rp[ii] + 1; pp < M->rp[ii + 1]; pp++)
            if (M->ci[pp] <= M->ci[pp - 1]) { clean = 0; break; }
    uint64_t setup_units;
    if (clean)
        setup_units = nzaat + 3ULL * (uint64_t) n + 1;
    else
    {
        uint64_t dirty = (uint64_t) M->nz + 2ULL * (uint64_t) n + 1;
        uint64_t clean_path = nzaat + 3ULL * (uint64_t) n + 1;
        setup_units = dirty > clean_path ? dirty : clean_path;
    }
    uint64_t kernel_units = 9ULL * (uint64_t) n + nzaat + nzaat / 5 + (uint64_t) n;
    uint64_t peak_units = setup_units > kernel_units ? setup_units : kernel_units;
    double mem_ours = (double) peak_units * (double) sizeof(int32_t);
    double mem_orig = oinfo[7];

    printf("%-22s n=%8d nz=%9d  lite=%8.3f ms  orig=%8.3f ms  ratio=%5.2fx  "
           "lnz lite=%11.0f orig=%11.0f  mem lite=%6.1fMB orig=%6.1fMB  ok=%d/%d st=%d/%d\n",
           name, n, M->nz, 1e3 * best1, 1e3 * best2,
           best1 / best2, info.lnz, oinfo[9],
           mem_ours / 1048576.0, mem_orig / 1048576.0,
           ok1, ok2, st1, st2);

    free(p1); free(p2);
}

int main(void)
{
    printf("AMD-lite vs SuiteSparse AMD (original built locally, -O3, same input)\n");
    printf("times are best of 3-7 runs; mem is approximate peak workspace\n\n");

    const int g2[] = { 100, 200, 300, 400, 500 };
    for (unsigned i = 0; i < sizeof(g2)/sizeof(g2[0]); i++)
    {
        Matrix M = make_grid2d(g2[i]);
        char name[64];
        snprintf(name, sizeof name, "grid2d_%dx%d", g2[i], g2[i]);
        run(name, &M, g2[i] <= 200 ? 7 : 5);
        mat_free(&M);
    }

    const int g3[] = { 20, 30, 40, 50 };
    for (unsigned i = 0; i < sizeof(g3)/sizeof(g3[0]); i++)
    {
        Matrix M = make_grid3d(g3[i]);
        char name[64];
        snprintf(name, sizeof name, "grid3d_%d^3", g3[i]);
        run(name, &M, 7);
        mat_free(&M);
    }

    const int bands[] = { 5, 20, 50 };
    for (unsigned i = 0; i < sizeof(bands)/sizeof(bands[0]); i++)
    {
        Matrix M = make_banded(100000, bands[i]);
        char name[64];
        snprintf(name, sizeof name, "banded_n100k_b%d", bands[i]);
        run(name, &M, 5);
        mat_free(&M);
    }

    const int degs[] = { 4, 10, 20 };
    for (unsigned i = 0; i < sizeof(degs)/sizeof(degs[0]); i++)
    {
        Matrix M = make_random_sym(100000, degs[i]);
        char name[64];
        snprintf(name, sizeof name, "random_n100k_d%d", degs[i]);
        run(name, &M, 5);
        mat_free(&M);
    }

    {
        Matrix M = make_star(100000);
        run("star_n100k", &M, 5);
        mat_free(&M);
    }

    {
        Matrix M = make_grid2d(300);
        reverse_rows(&M);
        run("grid2d_300_unsorted", &M, 5);
        mat_free(&M);
    }

    return 0;
}
