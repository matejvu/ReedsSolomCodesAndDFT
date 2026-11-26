


#include <stdlib.h>
#include "welch_berlekamp.h"
#include "fft.h"

void build_wb_system(ff_elem **M, ff_elem *rhs,
                     const ff_elem *y, const ff_elem *a_list,
                     uint64_t N, uint32_t n, uint32_t e, uint32_t p)
{
    for(uint64_t i=0 ; i < N ;i++){
        ff_elem pw=1;
        for (uint64_t j=0; j<n+e; j++){
            M[i][j]= pw;
            pw = ff_mul(pw, a_list[i] , p);

        }

        pw=1;
        for (uint64_t k=0; k<e; k++){
            M[i][n+e+k]= ff_neg(ff_mul(y[i],pw,p),p);
            pw = ff_mul(pw, a_list[i], p);
        }
        rhs[i]=ff_mul(y[i],pw,p);

    }

}


/* Deli G sa moničnim E (stepena e, E[e] == 1).
   G ima lenG koeficijenata, Q dobija lenG - e koeficijenata.
   Vraća 0 ako je ostatak 0, inače -1. */
int poly_div_monic(const ff_elem *G, uint32_t lenG,
                   const ff_elem *E, uint32_t e,
                   ff_elem *Q, uint32_t p)
{
    ff_elem *rem = malloc(lenG * sizeof(ff_elem));
    
    for( uint32_t i =0 ; i< lenG;i++){
        rem[i]=G[i];
    }

    uint32_t lenQ = lenG - e;

    for (int64_t k = (int64_t)lenQ - 1; k >= 0; k--) {
        /* TODO: Q[k] = ... (vodeći koeficijent ostatka) */
        Q[k] = rem[k+e];
        for (uint32_t j = 0; j <= e; j++) {
            rem[k+j] = ff_add(rem[k+j],  ff_neg(ff_mul(Q[k],E[j],p),p),p);
        }
    }
    for(uint32_t i =0 ; i< e ; i++){
        if(rem[i] != 0){
            free(rem);
            return -1;
        }   
    }
    
    free(rem);
    return 0;
}


/* Welch-Berlekamp dekoder.
   y: primljena reč (N brojeva, moguće sa greškama)
   w_out: ovde se upisuje dekodirana poruka (n brojeva)
   Vraća 0 ako je uspelo, -1 ako nije (previše grešaka). */
int wb_decode(const ff_elem *y, const ff_elem *a_list,
              uint64_t N, uint32_t n, uint32_t p, ff_elem *w_out)
{
    uint32_t e = (N - n) / 2; 
    uint32_t cols = n + 2 * e;

    ff_elem **M = alloc_matrix(N, cols);
    ff_elem *rhs = malloc(N * sizeof(ff_elem));
    ff_elem *x = malloc(cols * sizeof(ff_elem));


    build_wb_system( M, rhs, y, a_list, N,  n,  e ,  p);


    int ret = gaussian_elimination_general(M, rhs, x , N ,cols, p);

    /* 4) ako je Gaus uspeo: napravi E i podeli */
    if (ret == 0) {
        ff_elem *E = malloc((e + 1) * sizeof(ff_elem));
        for(uint32_t i =0 ; i < e; i++){
            E[i]=x[n+e+i];
        }
        E[e]=1;

        ret = poly_div_monic(x, n + e, E, e, w_out, p);

        free(E);
    }

    free_matrix(M);
    free(rhs);
    free(x);

    return ret;
}

/* Brzi Welch-Berlekamp dekoder (Task 4.3).
   Tačke su a_i = alpha^i, N stepen dvojke, p = N + 1.
   Vraća 0 ako je uspelo, -1 ako nije. */
int wb_decode_fft(const ff_elem *y, uint64_t N, uint32_t n,
                  ff_elem alpha, uint32_t p, ff_elem *w_out)
{
    uint32_t e = (N - n) / 2;

    /* KORAK 1: Y = inverzni DFT od y.
       ifft radi "u mestu", a y ne smemo da menjamo -> prvo kopija. */
    ff_elem *Y = malloc(N * sizeof(ff_elem));
    for (uint64_t i=0; i<N;i++){
        Y[i]=y[i];
    }

    ifft(Y,N,alpha,p);
    

    /* KORAK 2: napravi mali sistem B * E = desna  (e x e) */
    ff_elem **B = alloc_matrix(e, e);
    ff_elem *desna = malloc(e * sizeof(ff_elem));
    ff_elem *Ex = malloc(e * sizeof(ff_elem));      

    for (uint32_t r = 0; r < e; r++) {
        uint64_t j = n + e + r;
        for (uint32_t k = 0; k < e; k++) {
            B[r][k] = Y[j-k];
        }
        desna[r]=ff_neg(Y[j-e],p);
    }

    /* KORAK 3: reši mali sistem */
    int ret= gaussian_elimination_general(B, desna,Ex, e,e,p );
    /* KORAK 4: napravi E (dužine e+1, E[e] = 1) */
    if (ret == 0) {
        ff_elem *E = malloc((e + 1) * sizeof(ff_elem));
        for(uint32_t i =0 ; i < e; i++){
            E[i]=Ex[i];
        }
        E[e]=1;

        ff_elem *G= malloc((n+e)*sizeof(ff_elem));
        for(uint64_t j=0; j< n+e;j++){
            ff_elem s= 0;
            for(uint32_t k=0; k<=e;k++){
                uint64_t idx= (j+N-k)%N;
                s= ff_add(s,ff_mul(E[k],Y[idx],p),p);
            }
            G[j]=s;
        }

        ret = poly_div_monic(G, n + e, E, e, w_out, p);




        free(E);
        free(G);
    }

    free(Y);
  
    free_matrix(B);
    free(desna);
    free(Ex);
    return ret;
}