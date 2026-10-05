#include "encoding.h"

void encode_message(uint32_t p, ff_elem *w, uint32_t n, ff_elem *a_list, uint64_t N, ff_elem **c_ptr)
{
    //Forming the matrix F, based on interpolation points a.
    ff_elem** F = alloc_matrix(N, n);
    for(int i=0; i<N; i++){
        for(int j=0; j<n; j++)
            F[i][j] = ff_pow(a_list[i], j, p); //a[i]^j mod p
    }

    *c_ptr = (ff_elem*)malloc(sizeof(ff_elem) * N);
    ff_elem* c = *c_ptr;

    //Compute coding by matrix vector multiplication

    for(int i=0; i<N; i++){
        c[i] = ff_mul(c[i], 0, p);
        for(int k=0; k<n; k++){
            ff_elem term = ff_mul(F[i][k], w[k], p);
            c[i] = ff_add(c[i], term, p);
        }
    }

    free_matrix(F);
}

void decode_message(uint32_t p, ff_elem **w_ptr, uint32_t n, ff_elem *a_list, uint64_t N, ff_elem *c)
{
    //Forming the matrix F, based on first n interpolation points a
    ff_elem** F = alloc_matrix(n, n);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++)
            F[i][j] = ff_pow(a_list[i], j, p); //a[i]^j mod p
    }

    *w_ptr = (ff_elem*)malloc(sizeof(ff_elem) * n);
    ff_elem* w = *w_ptr;

    //Compute original message by gaussian elimination

    gaussian_elimination(F, c, w, n, p);

    free_matrix(F);
}
