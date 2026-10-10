#ifndef WELCH_BERLEKAMP_H
#define WELCH_BERLEKAMP_H

#include <stdint.h>
#include "finite_field_arithmetic.h"
#include "utils.h"

void build_wb_system(ff_elem **M, ff_elem *rhs,
                     const ff_elem *y, const ff_elem *a_list,
                     uint64_t N, uint32_t n, uint32_t e, uint32_t p);
int poly_div_monic(const ff_elem *G, uint32_t lenG, const ff_elem *E, uint32_t e, ff_elem *Q, uint32_t p);

int wb_decode(const ff_elem *y, const ff_elem *a_list, uint64_t N, uint32_t n, uint32_t p, ff_elem *w_out);
#endif