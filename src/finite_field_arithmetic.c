/*
 * finite_field_arithmetic.c
 *
 * Implementation of basic arithmetic in F_p for primes p < 2^32.
 * See finite_field_arithmetic.h for documentation.
 */

#include <assert.h>
#include "finite_field_arithmetic.h"

/* Reduces any signed 64-bit integer to the range [0, p-1]. */
static ff_elem reduce(int64_t x, uint32_t p)
{
    int64_t r = x % (int64_t)p;   /* in (-p, p), same sign as x */
    if (r < 0)
        r += p;
    return (ff_elem)r;
}

ff_elem ff_add(int64_t a, int64_t b, uint32_t p)
{
    /* Sum of two reduced values is < 2^33: no overflow on 64 bits. */
    uint64_t s = (uint64_t)reduce(a, p) + reduce(b, p);
    if (s >= p)
        s -= p;
    return (ff_elem)s;
}

ff_elem ff_neg(int64_t a, uint32_t p)
{
    ff_elem r = reduce(a, p);
    return (r == 0) ? 0 : p - r;
}

ff_elem ff_mul(int64_t a, int64_t b, uint32_t p)
{
    /* Product of two reduced values is < 2^64: no overflow on 64 bits. */
    return (ff_elem)(((uint64_t)reduce(a, p) * reduce(b, p)) % p);
}

ff_elem ff_pow(int64_t a, uint64_t e, uint32_t p)
{
    ff_elem base = reduce(a, p);
    ff_elem result = 1;

    /* Square-and-multiply. */
    while (e > 0) {
        if (e & 1)
            result = ff_mul(result, base, p);
        base = ff_mul(base, base, p);
        e >>= 1;
    }
    return result;
}

ff_elem ff_inv(int64_t a, uint32_t p)
{
    ff_elem x = reduce(a, p);
    assert(x != 0 && "ff_inv: 0 has no inverse");

    /* Extended Euclidean algorithm on (p, x), tracking only the
     * coefficient of x. Invariant: t * x = r (mod p). */
    int64_t t = 0, newt = 1;
    int64_t r = p, newr = x;

    while (newr != 0) {
        int64_t q = r / newr;
        int64_t tmp;

        tmp = t - q * newt;  t = newt;  newt = tmp;
        tmp = r - q * newr;  r = newr;  newr = tmp;
    }
    /* Now r = gcd(x, p) = 1, so t * x = 1 (mod p). */
    return reduce(t, p);
}

ff_elem ff_div(int64_t a, int64_t b, uint32_t p)
{
    return ff_mul(a, ff_inv(b, p), p);
}
