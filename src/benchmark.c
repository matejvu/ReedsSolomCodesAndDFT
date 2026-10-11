#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "utils.h"
#include "encoding.h"
#include "welch_berlekamp.h"
#include "fft.h"

/* Current CPU time in seconds */
static double now(void)
{
    return (double)clock() / CLOCKS_PER_SEC;
}

/* One benchmark row: N, n, alpha, number of repetitions,
   and whether to also run the slow (Task 3) decoder. */
void benchmark(uint64_t N, uint32_t n, ff_elem alpha, int repetitions, int run_slow)
{
    uint32_t p = N + 1;
    uint32_t e = (N - n) / 2;

    ff_elem *w   = malloc(n * sizeof(ff_elem));   /* original message */
    ff_elem *c   = malloc(N * sizeof(ff_elem));   /* codeword */
    ff_elem *y   = malloc(N * sizeof(ff_elem));   /* codeword with errors */
    ff_elem *w_d = malloc(n * sizeof(ff_elem));   /* decoded message */
    char *corrupted = malloc(N);                  /* is this position already corrupted? */

    /* evaluation points a_i = alpha^i, only needed by the slow decoder */
    ff_elem *a = NULL;
    if (run_slow) {
        a = malloc(N * sizeof(ff_elem));
        for (uint64_t i = 0; i < N; i++) a[i] = ff_pow(alpha, i, p);
    }

    double t_fast = 0, t_slow = 0;
    int ok_fast = 1, ok_slow = 1;

    for (int r = 0; r < repetitions; r++) {

        /* 1) random message */
        for (uint32_t i = 0; i < n; i++) {
            w[i] = rand() % p;
        }

        for (uint64_t i = 0; i < N; i++) {
            if(i<n) c[i]=w[i];
            else c[i]=0;
        }
        fft(c ,N, alpha,p);

        /* 3) exactly e errors at distinct random positions */
        for (uint64_t i = 0; i < N; i++) {
            y[i] = c[i];
            corrupted[i] = 0;
        }
        uint32_t num_errors = 0;
        while (num_errors < e) {
            uint64_t pos = rand() % N;
            if (corrupted[pos]) continue;                        /* already used -> pick again */
            corrupted[pos] = 1;
            y[pos] = ff_add(y[pos], 1 + rand() % (p - 1), p);    /* add a non-zero value */
            num_errors++;
        }

        /* 4) fast decoder: measure time and check the result */
        double t0 = now();
        int ret = wb_decode_fft(y,N,n,alpha,p,w_d);
        t_fast += now() - t0;

        if (ret != 0) ok_fast = 0;
        for (uint32_t i = 0; i < n; i++) {
            if(w_d[i]!=w[i]) ok_fast=0;
        }

        /* 5) slow decoder (only if run_slow) */
        if (run_slow) {
            double t1 = now();
            int ret_slow = wb_decode(y, a, N, n, p, w_d);
            t_slow += now() - t1;

            if (ret_slow != 0) ok_slow = 0;
            for (uint32_t i = 0; i < n; i++) {
                if (w_d[i] != w[i]) ok_slow = 0;
            }
        }
    }

    printf("N=%6llu n=%6u e=%4u | fast: %.6f s (%s)",
           (unsigned long long)N, n, e, t_fast / repetitions, ok_fast ? "OK" : "ERROR");
    if (run_slow)
        printf(" | slow: %.6f s (%s) | speedup: %.1fx",
               t_slow / repetitions, ok_slow ? "OK" : "ERROR", t_slow / t_fast);
    printf("\n");

    free(w); free(c); free(y); free(w_d); free(corrupted);
    if (a) free(a);
}

int main(void)
{
    srand(time(NULL));
    benchmark(16,    8,     3, 1000, 1);   /* both decoders, 1000 repetitions */
    benchmark(256,   64,    3, 10,   1);   /* both decoders, 10 repetitions */
    benchmark(65536, 64512, 3, 1,    0);   /* fast decoder only, once */
    return 0;
}