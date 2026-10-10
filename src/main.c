#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "encoding.h"
#include "welch_berlekamp.h"


int main(){

    ff_elem** A = alloc_matrix(3, 3);
    ff_elem* b = malloc(3 * sizeof(ff_elem));
    ff_elem* x = malloc(3 * sizeof(ff_elem));

    // Initialize A and b
    A[0][0] = 1; A[0][1] = 2; A[0][2] = 3; b[0] = 1;
    A[1][0] = 0; A[1][1] = 1; A[1][2] = 4; b[1] = 2;
    A[2][0] = 5; A[2][1] = 6; A[2][2] = 0; b[2] = 3;

    // Perform Gaussian elimination
    if (gaussian_elimination(A, b, x, 3, 7) == 0) {
        printf("Solution found:\n");
        for (uint32_t i = 0; i < 3; i++) {
            printf("x[%u] = %u\n", i, x[i]);
        }
    } else {
        printf("No solution exists.\n");
    }

    // Clean up
    free_matrix(A);
    free(b);
    free(x);

    //Test encoding

    uint32_t p = 7; // Example prime
    ff_elem* w = (ff_elem*)malloc(sizeof(ff_elem) * 3);
    ff_elem* a_list = (ff_elem*)malloc(sizeof(ff_elem) * 4);
    ff_elem* c; //= (ff_elem*)malloc(sizeof(ff_elem) * 3);

    // Initialize w and a_list
    w[0] = 1; w[1] = 5; w[2] = 4;
    a_list[0] = 1; a_list[1] = 2; a_list[2] = 3; a_list[3] = 4;

    encode_message(p, w, 3, a_list, 4, &c);

    printf("Encoded message:\n");
    for (uint32_t i = 0; i < 4; i++) {
        printf("c[%u] = %u\n", i, c[i]);
    }
    free(w);

    decode_message(p, &w, 3, a_list, 4, c);

    printf("Decoded message:\n");
    for (uint32_t i = 0; i < 3; i++) {
        printf("w[%u] = %u\n", i, w[i]);
    }

    free(w);
    free(a_list);
    free(c);

    uint32_t p2 = 17, n2 = 8, e2 = 4;
    uint64_t N2 = 16;
    ff_elem a2[16];
    ff_elem w2[8] = {1, 2, 3, 4, 5, 6, 7, 8};
    ff_elem *c2;

    // tačke a_i = 3^i mod 17
    for (int i = 0; i < 16; i++) a2[i] = ff_pow(3, i, p2);

    // kodiramo poruku w2 -> kodna reč c2 (16 brojeva)
    encode_message(p2, w2, n2, a2, N2, &c2);

    // pravimo sistem M * [g | E] = rhs
    ff_elem **M = alloc_matrix(N2, n2 + 2 * e2);
    ff_elem *rhs = malloc(N2 * sizeof(ff_elem));
    build_wb_system(M, rhs, c2, a2, N2, n2, e2, p2);

    printf("Matrica M:\n");
    print_matrix(M, N2, n2 + 2 * e2);

    free_matrix(M);
    free(rhs);
    free(c2);


    // ===== Test novog Gausa =====
    ff_elem rhs1[3] = {1, 0, 1};   // ima rešenje
    ff_elem rhs2[3] = {1, 0, 2};   // nema rešenje
    ff_elem *testovi[2] = {rhs1, rhs2};

    for (int t = 0; t < 2; t++) {
        ff_elem **G = alloc_matrix(3, 3);
        ff_elem xg[3];
        G[0][0] = 1; G[0][1] = 2; G[0][2] = 3;
        G[1][0] = 2; G[1][1] = 4; G[1][2] = 1;
        G[2][0] = 3; G[2][1] = 6; G[2][2] = 4;

        int ret = gaussian_elimination_general(G, testovi[t], xg, 3, 3, 7);
        printf("Test %d: povratna vrednost = %d\n", t + 1, ret);
        if (ret == 0) printf("x = (%u, %u, %u)\n", xg[0], xg[1], xg[2]);
        free_matrix(G);
    }

    // ===== Test deljenja polinoma =====
    ff_elem Gt[4] = {2, 0, 5, 1};   // x^3 + 5x^2 + 2
    ff_elem Et[2] = {2, 1};         // x + 2
    ff_elem Qt[3];
    int rd = poly_div_monic(Gt, 4, Et, 1, Qt, 7);
    printf("Deljenje: povratna = %d, Q = (%u, %u, %u)\n", rd, Qt[0], Qt[1], Qt[2]);


    // ===== Test Welch-Berlekamp dekodera =====
    uint32_t p3 = 17, n3 = 8;
    uint64_t N3 = 16;
    ff_elem a3[16];
    ff_elem w3[8] = {1, 2, 3, 4, 5, 6, 7, 8};
    ff_elem *c3;
    for (int i = 0; i < 16; i++) a3[i] = ff_pow(3, i, p3);
    encode_message(p3, w3, n3, a3, N3, &c3);

    int pozicije[5] = {1, 4, 7, 10, 13};   // gde kvarimo

    for (int ne = 0; ne <= 5; ne++) {
        ff_elem y[16], wd[8] = {0};
        for (int i = 0; i < 16; i++) y[i] = c3[i];                  // kopija kodne reči
        for (int k = 0; k < ne; k++)
            y[pozicije[k]] = ff_add(y[pozicije[k]], 5, p3);         // pokvari ne slova

        int r = wb_decode(y, a3, N3, n3, p3, wd);
        printf("Greske: %d -> povratna = %d, w = ", ne, r);
        for (int i = 0; i < 8; i++) printf("%u ", wd[i]);
        printf("\n");
    }
    free(c3);

    return 0;
}