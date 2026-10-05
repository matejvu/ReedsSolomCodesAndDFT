#include "qr_decomposition.h"

void gram_schmidt_procedure(ff_elem **A, ff_elem **Q, ff_elem **R, uint64_t m, uint64_t n)
{
    //Reduced versions
    Q = alloc_matrix(m, n);
    R = alloc_matrix(n, n);

    for(int i = 0; i<n; i++){
        ff_elem* ui = take_column(A, m, n, i);
        
    }
}

ff_elem* take_column(ff_elem** A, uint64_t n, uint64_t m, uint64_t col){
    
    ff_elem* column_vect = (ff_elem*)malloc(sizeof(ff_elem) * n);
    for(int k = 0; k<n; k)
        column_vect[k] = A[k][col];
}

ff_elem norm(ff_elem* vect, uint64_t n, uint32_t p){
    ff_elem total = 0;
    for(int i=0; i<n; i++){
        ff_elem term = vect[i];
        term = ff_mul(term, term, p);
        total = ff_add(total, term, p);
    }
    // total = ff_pow()/
}