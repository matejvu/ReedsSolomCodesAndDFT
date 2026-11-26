#ifndef UTILS_H
#define UTILS_H

#include "finite_field_arithmetic.h"

ff_elem** alloc_matrix(uint64_t rows, uint64_t cols);
void free_matrix(ff_elem** m);
void print_matrix(ff_elem** m, uint64_t rows, uint64_t cols);
int gaussian_elimination(ff_elem** A, ff_elem* b, ff_elem* x, uint64_t n, uint32_t p);
int gaussian_elimination_general(ff_elem** A, ff_elem* b, ff_elem* x, uint64_t rows, uint64_t cols, uint32_t p);
#endif /* UTILS_H */
