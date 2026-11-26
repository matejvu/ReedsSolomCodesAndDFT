#include <stdlib.h>
#include "fft.h"

void fft(ff_elem *v, uint64_t N, ff_elem omega, uint32_t p){
    if (N ==1) return;

    uint64_t h = N/2;
    ff_elem *par = malloc(h * sizeof(ff_elem));
    ff_elem *nepar = malloc(h * sizeof(ff_elem));

    for(uint64_t i=0; i<h;i++){
        par[i] = v[2 *i];
        nepar[i] = v[2*i+1];
    }

    ff_elem omega2 = ff_pow(omega,2,p);
    fft(par,h,omega2,p);
    fft(nepar,h,omega2,p);

    ff_elem w = 1;                    
    for (uint64_t k = 0; k < h; k++) {
        ff_elem t = ff_mul(w,nepar[k],p);
        v[k]= ff_add(par[k],t,p);
        v[k + h] = ff_add(par[k],ff_neg(t,p),p);
        w = ff_mul(w,omega,p);
    }

    free(par);
    free(nepar);
}

void ifft(ff_elem *v, uint64_t N, ff_elem alpha, uint32_t p){
    fft(v,N,ff_inv(alpha,p),p);


    ff_elem n_inv = ff_inv(N,p);
    for (uint64_t i=0;i < N;i++){
        v[i]=ff_mul(v[i],n_inv,p);
    }
}