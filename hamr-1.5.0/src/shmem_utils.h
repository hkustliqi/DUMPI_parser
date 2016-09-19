#ifndef __SHMEM_UTILS_H
#define __SHMEM_UTILS_H

//#include <mpp/shmem.h>
#include <mpi.h>


//#define  RANK       shmem_my_pe()
//#define  NPES       shmem_n_pes()
//#define  BARRIER()  shmem_barrier_all()



//#define _printf  if (RANK == 0) printf 
//#define _fprintf if (RANK == 0) fprintf 

/**
 * Compares SHMEM symmetric memory pointers to determine consistency.  If not 
 * consistent then exits here. 
 *
 * @param rank Rank of the calling process or thread.
 * @param npes Total number of processes.
 * @param s    Address to be checked.
 * @param str  String holding the name of the shmalloc'd array to be checked.
 */
void check_shptrs( void *s, char *str, int RANK, int NPES );

/**
 * Wrapper for shmalloc that also checks to determine if the allocation
 * is consistent across PEs.  Exits immediately if not and displays the
 * file and line number of the failed allocation.
 * 
 * @param size Number of bytes to allocate.
 * @param file Normally the __FILE__ constant should be passed in.
 * @param line Normally the __LINE__ constant should be passed in.
 */
void * shmalloc_safe( size_t size, char *file, int line, int RANK, int NPES );

/**
 * Broadcasts val from root to all PEs.  The broadcasted value is returned.
 */
void  mpp_broadcast_long( long *val, int root, int RANK, int NPES );
long  mpp_accum_long( long val, int RANK, int NPES );
long  mpp_accum_or_long( long val, int RANK, int NPES );

double  mpp_accum_double( double val, int RANK, int NPES );
double  mpp_min_double( double val, int RANK, int NPES );
double  mpp_max_double( double val, int RANK, int NPES );
#endif 
