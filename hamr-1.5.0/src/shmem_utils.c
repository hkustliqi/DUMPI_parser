#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
//#include <mpp/shmem.h>
#include <mpi.h>

#include "shmem_utils.h"

#define _SHMEM_REDUCE_SYNC_SIZE 32
#define _SHMEM_BCAST_SYNC_SIZE 32
#define _SHMEM_REDUCE_MIN_WRKDATA_SIZE 32
#define _SHMEM_SYNC_VALUE 1


//#define  RANK       shmem_my_pe()
//#define  NPES       shmem_n_pes()
//#define  BARRIER()  shmem_barrier_all()

 

#define  MIN(A,B)   (((A)<(B))?(A):(B))
#define  MAX(A,B)   (((A)>(B))?(A):(B))


int shmem_check_shptrs( void *s, char *str, int RANK, int NPES) {
  int                  i;
  static uint64_t shptrs;
  long            reduced;
  
  // Place pointer in variable depending on id
  shptrs = (long)s;

  // Global OR reduction to all, then verify local pointer matches
  // the reduced value.
  reduced = mpp_accum_or_long( shptrs, RANK, NPES );
  
  if( reduced != shptrs ) 
    return 0;
 
  return 1;
}

void * shmalloc_safe( size_t size, char *file, int line, int RANK, int NPES ) {

#define _printf  if (RANK == 0) printf 
#define _fprintf if (RANK == 0) fprintf


  void *buffer;
  //buffer = shmalloc(size);
  buffer = malloc(size);

  if( !buffer ) {
    fprintf(stderr, "ERROR: shmalloc failed from allocation at %s:%d.\n", 
	    file, line);
    exit(1);
  }
  if( !shmem_check_shptrs( buffer, file, RANK, NPES ) ) {
    //fprintf(stderr, "ERROR: Inconsistent shmem pointers "
    //    "from allocation at %s:%d\n", file, line);
    //exit(1);
  }
  return buffer;
}

/* void * shmalloc_aligned_safe( void *aligned, size_t size, int alignment, char *file, int line ) { */
/*   uintptr_t mask = ~(uintptr_t)(alignment - 1); */
/*   void * mem = shmalloc_safe ( size + alignment, file, line ); */
/*   aligned = (void *)(((uintptr_t)mem + alignment - 1) & mask); */
/*   printf("aligned Returning %lx\n", (long) mem); */
/*   return mem; */
/* } */

long  mpp_accum_long( long val, int RANK, int NPES  ) {
  int i;
  static long dst, src;
  static long Sync[_SHMEM_REDUCE_SYNC_SIZE];
  static long Work[MAX(4, _SHMEM_REDUCE_MIN_WRKDATA_SIZE)];
  static int init=0;

  if (! init) {
    for (i = 0; i < _SHMEM_REDUCE_SYNC_SIZE; i++)
      Sync[i] = _SHMEM_SYNC_VALUE;
    init = 1;
  }

  src = val;
  //BARRIER();
  MPI_Barrier(MPI_COMM_WORLD);
  //shmem_long_sum_to_all (&dst, &src, 1, 0, 0, NPES, Work, Sync);
  MPI_Allreduce (&src, &dst, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

  return dst;
}

long  mpp_accum_or_long( long val, int RANK, int NPES ) {
  int i;
  static long dst, src;
  static long Sync[_SHMEM_REDUCE_SYNC_SIZE];
  static long Work[MAX(4, _SHMEM_REDUCE_MIN_WRKDATA_SIZE)];
  static int init=0;

  if (! init) {
    for (i = 0; i < _SHMEM_REDUCE_SYNC_SIZE; i++)
      Sync[i] = _SHMEM_SYNC_VALUE;
    init = 1;
  }
  
  src = val;
  //BARRIER();
  MPI_Barrier(MPI_COMM_WORLD);
  //shmem_long_or_to_all (&dst, &src, 1, 0, 0, NPES, Work, Sync);
  MPI_Allreduce (&src, &dst, 1, MPI_LONG, MPI_LOR, MPI_COMM_WORLD);

  return dst;
}

void mpp_broadcast_long( long *val, int root, int RANK, int NPES  ) {
  int i;
  static long dst, src;
  static long Sync[_SHMEM_BCAST_SYNC_SIZE];
  static int init=0;
  
  if (! init) {
    for (i = 0; i < _SHMEM_BCAST_SYNC_SIZE; i++)
      Sync[i] = _SHMEM_SYNC_VALUE;
    init = 1;
  }

  src = *val;
  //BARRIER();
  MPI_Barrier(MPI_COMM_WORLD);
  //shmem_broadcast64 (&dst, &src, 1, root, 0, 0, NPES, Sync);
  MPI_Bcast (&src, 1, MPI_UINT64_T, root, MPI_COMM_WORLD); 

  if( RANK != root )
    *val = dst;
}

double  mpp_accum_double( double val, int RANK, int NPES ) {
  int i;
  static double dst, src;
  static long Sync[_SHMEM_REDUCE_SYNC_SIZE];
  static double Work[MAX(4, _SHMEM_REDUCE_MIN_WRKDATA_SIZE)];
  static int init=0;

  if (! init) {
    for (i = 0; i < _SHMEM_REDUCE_SYNC_SIZE; i++)
      Sync[i] = _SHMEM_SYNC_VALUE;
    init = 1;
  }

  src = val;
  //BARRIER();
  MPI_Barrier(MPI_COMM_WORLD);
  //shmem_double_sum_to_all (&dst, &src, 1, 0, 0, NPES, Work, Sync);
  MPI_Allreduce (&src, &dst, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return dst;
}

double  mpp_max_double( double val, int RANK, int NPES  ) {
  int i;
  static double dst, src;
  static long Sync[_SHMEM_REDUCE_SYNC_SIZE];
  static double Work[MAX(4, _SHMEM_REDUCE_MIN_WRKDATA_SIZE)];
  static int init=0;

  if (! init) {
    for (i = 0; i < _SHMEM_REDUCE_SYNC_SIZE; i++)
      Sync[i] = _SHMEM_SYNC_VALUE;
    init = 1;
  }

  src = val;
  //BARRIER();
  MPI_Barrier(MPI_COMM_WORLD);
  //shmem_double_max_to_all (&dst, &src, 1, 0, 0, NPES, Work, Sync);
  MPI_Allreduce (&src, &dst, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  return dst;
}

double  mpp_min_double( double val, int RANK, int NPES ) {
  int i;
  static double dst, src;
  static long Sync[_SHMEM_REDUCE_SYNC_SIZE];
  static double Work[MAX(4, _SHMEM_REDUCE_MIN_WRKDATA_SIZE)];
  static int init=0;

  if (! init) {
    for (i = 0; i < _SHMEM_REDUCE_SYNC_SIZE; i++)
      Sync[i] = _SHMEM_SYNC_VALUE;
    init = 1;
  }

  src = val;
  //BARRIER();
  MPI_Barrier(MPI_COMM_WORLD);
  //shmem_double_min_to_all (&dst, &src, 1, 0, 0, NPES, Work, Sync);
  MPI_Allreduce (&src, &dst, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  return dst;
}
