#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "encoding.h"


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

    return 0;
}