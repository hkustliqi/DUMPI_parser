#include <stdio.h>
#include <assert.h>
#include <stdint.h>

#include "uint128.h"

int mask_bit_present( uint128 x, uint128 m ) {
  uint128 guard;

  guard.as64.word[0] = x.as64.word[0] & m.as64.word[0];
  guard.as64.word[1] = x.as64.word[1] & m.as64.word[1];
  
  return guard.as64.word[0] || guard.as64.word[1];
}

int gt128( uint128 x, uint128 y ) {
  if (x.as64.word[0] != y.as64.word[0])
    return x.as64.word[0] > y.as64.word[0];
  else
    return x.as64.word[1] > y.as64.word[1];
}
  
void _insertion_sort128( uint128 *a, uint64_t n )
{
  uint64_t i, j;
  uint128 look;

  if( n == 1 )
    return;

  for( i = 1; i < n; i++ ) {
    look = a[i];
    
    j = i - 1;
    while( j >= 0 && gt128( a[j], look ) ) {
      a[j+1] = a[j];
      j--;
    }
    a[j+1] = look;
  }
  return;
}

void _msb_radix_qsort128( uint128 *data, int64_t left, int64_t right, unsigned low_bit, unsigned high_bit ) {
  uint64_t b, i, j;
  uint128 mask;

  mask.as64.word[0] = 0x8000000000000000L;
  mask.as64.word[1] = 0x8000000000000000L;

  if( right <= left )
    return;

  // Switch to insertion sort for small arrays.
  /* if( right - left <= 16 ) { */
  /*   _insertion_sort128( &data[left], right + 1 - left ); */
  /* } */

  if( high_bit >= 64 ) {
    mask.as64.word[0] >>= (128 - high_bit - 1); 
    mask.as64.word[1] = 0;
  }
  else {
    mask.as64.word[0] = 0;
    mask.as64.word[1] >>= (64 - high_bit - 1);
  }

  i = left;
  j = right;

  while ( i != j ) {
    while ( !mask_bit_present( data[i], mask ) && i < j )
      i++; 

    while ( mask_bit_present( data[j], mask ) && j > i ) 
      j--;
   
    SWAP128( data[j], data[i] );
  }
  
  if ( !mask_bit_present( data[right], mask ) )
    j++;

  _msb_radix_qsort128( data, left, j-1, low_bit, high_bit-1);
  _msb_radix_qsort128( data, j,  right, low_bit, high_bit-1);
}

void msb_radix_quicksort128( uint128 *data, uint64_t size, unsigned low_bit, unsigned high_bit ) {
  _msb_radix_qsort128( data, 0, size - 1, low_bit, high_bit );
}

