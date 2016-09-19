#ifndef __SORT128_H_
#define __SORT128_H_

#include <stdint.h>

#include "uint128.h"

void msb_radix_quicksort128( uint128 *data, uint64_t length, unsigned low_bit, unsigned hi_bit );

#endif
