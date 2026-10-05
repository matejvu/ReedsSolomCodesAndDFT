#ifndef ENCODING_H
#define ENCODING_H

#include <stdint.h>
#include <stdlib.h>
#include "finite_field_arithmetic.h"
#include "utils.h"

// Encodes the message w given in module p with interpolation points from a_list length N.
// Returns c, the encoded message.
void encode_message(uint32_t p, ff_elem* w, uint32_t n, ff_elem* a_list, uint64_t N, ff_elem** c);
//Decodes the message c(w) given in module p with interpolation points from a_list length N.
//Return w, the original message.
void decode_message(uint32_t p, ff_elem** w, uint32_t n, ff_elem* a_list, uint64_t N, ff_elem* c);
 

#endif /* ENCODING_H */