#ifndef QR_DECOMPOSITION_H
#define QR_DECOMPOSITION_H

#include "finite_field_arithmetic.h"
#include "utils.h"

void gram_schmidt_procedure(ff_elem** A, ff_elem** Q, ff_elem** R, uint64_t m, uint64_t n);

#endif /* QR_DECOMPOSITION_H */
