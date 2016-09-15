#ifndef __NETSTRESS_UTILS_H
#define __NETSTRESS_UTILS_H

/**
 * @file netstress_utils.h
 *
 * @section DESCRIPTION
 * These are utility functions that are used within netstress.
 *
 */

#include <getopt.h>
#include <float.h>

#include "netstress.h"

typedef struct bench_t{
  long long min_size; 
  long long max_size; 
  int iter; 
  long seed; 
  int ngrps; 
  int quiet;
  int very_quiet;
  int dry_run;
  int alltoall;
  double timeout;
}bench_t;

/**
 * Handle the user options that are used by the benchmark.
 * 
 * @param rank  [in]  Rank of calling PE.
 * @param bench [out] Benchmark options struct.
 */
void handle_options(int argc, char *argv[], int rank, bench_t *bench);

/**
 * Determine if PEs are collocated on a single node.
 *
 * @param me  [in]  Block to be inverted.
 * @param ppg [in]  The block one row under x.
 * @return
 *       0 if all PEs in group are collocated on single node.
 *       1 if there is a PE in a group that is not collocated on the node.
 */
int check_collocation(int me, int ppg);

/**
 * Determine units of a memory value.
 *
 * @param bytes [in]  Size in bytes.
 * @param out   [out] Units.
 * @return 
 *       Converted value of 'bytes' in 'out' units.
 */
double units_of(uint64_t bytes, char *out);

/**
 * Set seed as the same on all PEs.
 *
 * @param me           [in]  Rank of calling PE.
 * @param npes         [in]  Total number of calling PEs.
 * @param initial_seed [out] Default seed (epoch time).
 * @return 
 *       Seed, which is now the same on all PEs.
 */
long set_seed(int me, int npes, unsigned long initial_seed);

#endif

