/*
 * finite_field_arithmetic.h
 *
 * Basic arithmetic in the prime field K := Z/pZ, for primes p < 2^32.
 *
 * - The prime p is passed as the last argument of every function.
 *   It must be a prime with 2 <= p < 2^32.
 * - Inputs can be any integers, including negative ones or values larger
 *   than p: they are reduced modulo p automatically.
 * - Results are always returned in the range [0, p-1].
 *
 * Example:
 *   #include "finite_field_arithmetic.h"
 *
 *   uint32_t p = 101;
 *   ff_elem c = ff_mul(45, 77, p);   // 45 * 77 mod 101
 *   ff_elem d = ff_div(c, 77, p);    // back to 45
 *   ff_elem e = ff_neg(1, p);        // 100, i.e. -1 mod 101
 *
 * Compile together with your own code, e.g.:
 *   gcc -Wall -O2 main.c finite_field_arithmetic.c -o main
 */

#ifndef FINITE_FIELD_ARITHMETIC_H
#define FINITE_FIELD_ARITHMETIC_H

#include <stdint.h>

/* A field element, in the range [0, p-1]. */
typedef uint32_t ff_elem;

/* Returns a + b in K. */
ff_elem ff_add(int64_t a, int64_t b, uint32_t p);

/* Returns -a in K. */
ff_elem ff_neg(int64_t a, uint32_t p);

/* Returns a * b in K. */
ff_elem ff_mul(int64_t a, int64_t b, uint32_t p);

/* Returns a^e in K, for a non-negative exponent e.
 * By convention, a^0 = 1 (including 0^0). */
ff_elem ff_pow(int64_t a, uint64_t e, uint32_t p);

/* Returns the inverse of a in K, i.e. the x with a * x = 1 in K.
 * a must be non-zero modulo p (checked with assert). */
ff_elem ff_inv(int64_t a, uint32_t p);

/* Returns (a / b) in K, i.e. a * b^(-1) in K.
 * b must be non-zero modulo p (checked with assert). */
ff_elem ff_div(int64_t a, int64_t b, uint32_t p);

#endif /* FINITE_FIELD_ARITHMETIC_H */
