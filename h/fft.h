#ifndef FFT_H
#define FFT_H

#include <stdint.h>
#include "finite_field_arithmetic.h"

/* v <- DFT od v. N mora biti stepen dvojke. */
void fft(ff_elem *v, uint64_t N, ff_elem omega, uint32_t p);

/* v <- inverzni DFT od v  */
void ifft(ff_elem *v, uint64_t N, ff_elem alpha, uint32_t p);

#endif