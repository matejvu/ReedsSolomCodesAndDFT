#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

ff_elem** alloc_matrix(uint64_t rows, uint64_t cols)
{
    ff_elem** m = malloc(rows * sizeof(ff_elem*));
    if (!m) return NULL;
    m[0] = calloc(rows * cols, sizeof(ff_elem));
    if (!m[0]) { free(m); return NULL; }
    for (uint64_t i = 1; i < rows; i++) m[i] = m[0] + i * cols;
    return m;
}

void free_matrix(ff_elem** m)
{
    if (!m) return;
    free(m[0]);
    free(m);
}

void print_matrix(ff_elem** m, uint64_t rows, uint64_t cols)
{
    for (uint64_t i = 0; i < rows; i++) {
        for (uint64_t j = 0; j < cols; j++) printf("%u ", m[i][j]);
        printf("\n");
    }
}

int gaussian_elimination(ff_elem** A, ff_elem* b, ff_elem* x, uint64_t n, uint32_t p)
{
    for (uint64_t c = 0; c < n; c++) {
        uint64_t piv = c;
        while (piv < n && A[piv][c] == 0) piv++;
        if (piv == n) return -1;
        if (piv != c) {
            for (uint64_t j = c; j < n; j++) {
                ff_elem t = A[c][j]; A[c][j] = A[piv][j]; A[piv][j] = t;
            }
            ff_elem t = b[c]; b[c] = b[piv]; b[piv] = t;
        }
        ff_elem inv = ff_inv(A[c][c], p);
        for (uint64_t j = c; j < n; j++) A[c][j] = ff_mul(A[c][j], inv, p);
        b[c] = ff_mul(b[c], inv, p);
        for (uint64_t r = 0; r < n; r++) {
            if (r == c || A[r][c] == 0) continue;
            ff_elem f = A[r][c];
            for (uint64_t j = c; j < n; j++) A[r][j] = ff_add(A[r][j], ff_neg(ff_mul(f, A[c][j], p), p), p);
            b[r] = ff_add(b[r], ff_neg(ff_mul(f, b[c], p), p), p);
        }
    }
    for (uint64_t i = 0; i < n; i++) x[i] = b[i];
    return 0;
}
