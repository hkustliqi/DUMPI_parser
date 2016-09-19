#ifndef __UINT128_H_
#define __UINT128_H_

#include <stdint.h>

#define SWAP128(X,Y) {uint128 t; t = (X); (X) = (Y); (Y) = t;}

typedef union {
  struct {
    uint64_t word[2]; 
  } as64;
  // Architecture dependent vector types here.
} uint128;

#endif
