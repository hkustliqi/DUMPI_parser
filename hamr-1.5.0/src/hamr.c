#include <errno.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <assert.h>

//#include <mpp/shmem.h>
#include <mpi.h>

#include "config.h"
#include "uint128.h"
#include "parse_opts.h"
#include "timer.h"
#include "rand64.h"
#include "shmem_utils.h"
#include "sort128.h"

#define  GIBI  1073741824L
#define  MEBI  1048576



static int g_npes;              ///< Number of PEs
static int g_pe_bits;           ///< ceil(log_2( npes ))


#define  PE_OF(X)         ((X) >> (64 - g_pe_bits))
#define  IS_POW_OF_TWO(X) (!((X) & ((X) - 1)))

void
randomize(uint128 * A, uint64_t len, int high_only)
{
  uint64_t i;
  uint64_t nibble;
  uint64_t r;

  /*
   * Generating a random 64 bit integer for both the high and low words
   * starts makes the benchmark take too long for large problem sizes, 
   * so we only randomize the high word after each iteration.
   */
  for (i = 0; i < len; i++) {
    A[i].as64.word[0] = rand64();
    if (!high_only)
      A[i].as64.word[1] = rand64();
  }

  /*
   * If we don't have a power of two number of PEs then take the low 
   * g_pe_bits of the high word of each key, manipulate it to be in the 
   * range 0 to g_npes, then OR it into the high bits of the high word.
   */
  if (!IS_POW_OF_TWO(g_npes)) {
    for (i = 0; i < len; i++) {
      nibble = A[i].as64.word[0] % g_npes;
      A[i].as64.word[0] = A[i].as64.word[0] >> g_pe_bits;
      A[i].as64.word[0] |= nibble << (64 - g_pe_bits);
    }
  }
}

void
verify_exchange(const uint128 * table, uint64_t fill_count,
                uint64_t recv_count, int quiet, int RANK, int NPES)
{
#define _printf  if (RANK == 0) printf 
#define _fprintf if (RANK == 0) fprintf


  uint64_t incorrect_count = 0;
  uint64_t total_keys = 0;
  uint64_t ii;

  for (ii = 0; ii < recv_count; ii++) {
    if (RANK != PE_OF(table[ii].as64.word[0])) {
      incorrect_count++;
    }
  }
  if (incorrect_count) {
    fprintf(stderr, "(%2d) incorrect_count = %ld/%ld (%.2lf%%)\n",
            RANK, incorrect_count, recv_count,
            100 * (double) incorrect_count / recv_count);
  }

  //BARRIER();
  MPI_Barrier(MPI_COMM_WORLD);

  incorrect_count = mpp_accum_long(incorrect_count, RANK, NPES);
  total_keys = mpp_accum_long(recv_count, RANK, NPES);

  /*
   * Each PE starts with fill_count keys so we check to make sure that the 
   * sum of all recv_counts is fill_count * g_npes.
   */
  if (RANK == 0 && !quiet && (total_keys != fill_count * g_npes)) {
    fprintf(stderr, "WARNING: final table size is missing/has extra elements! "
            "(%ld at end vs %ld to start). "
            "Try reducing the block size (see -h).\n",
            total_keys, fill_count * g_npes);
  }

  if (!quiet) {
    _printf("Error rate = %ld/%ld (%lf%%)\n", incorrect_count, total_keys,
            100 * (double) incorrect_count / total_keys);
  }
}

void
put_buckets(int dest, const uint128 * restrict send_buckets,
            uint64_t *restrict send_counts, uint64_t bucket_len,
            uint128 * restrict recv_buff, uint64_t *recv_count)
{
  long dest_index;

  if (send_counts[dest] > 0) {
    //dest_index = shmem_long_fadd((long *) recv_count, send_counts[dest], dest);
 /*   MPI_Win win;
    MPI_Win_create((long *) recv_count, send_counts[dest], sizeof(MPI_LONG), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);
    MPI_Fetch_and_op(&send_counts[dest], &dest_index, MPI_LONG, dest, 0, MPI_SUM, win);
    MPI_Win_fence(0, win);
    MPI_Win_free(&win);
*/
//Qi: the following MPI_GET was causing segamentation fault. Disabling it for now.
//    long * recvData;

/*    MPI_Win_create(&recvData, 1, sizeof(MPI_LONG), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);
    MPI_Get((long *) recv_count, 1, MPI_LONG, dest, 0, 1, MPI_LONG, win);
    MPI_Win_fence(0, win);
    MPI_Win_free(&win);    
    dest_index = send_counts[dest] + *recvData;
*/
    //dest_index = send_counts[dest];

    //shmem_put128(&recv_buff[dest_index], &send_buckets[dest * bucket_len],
    //             send_counts[dest], dest);
    MPI_Win win;
    MPI_Win_create(&recv_buff[dest_index], send_counts[dest], sizeof(MPI_LONG_DOUBLE), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);
    MPI_Put(&send_buckets[dest * bucket_len], send_counts[dest], MPI_LONG_DOUBLE, dest, 0, send_counts[dest], MPI_LONG_DOUBLE, win);
    MPI_Win_fence(0, win);
    MPI_Win_free(&win);


    /*MPI_Request req;
    MPI_Irecv(&recv_buff[0], send_counts[dest], MPI_BYTE,
		dest, 0, MPI_COMM_WORLD, &req);
    MPI_Send(&send_buckets[dest * bucket_len], send_counts[dest], MPI_BYTE,
	       dest, 0, MPI_COMM_WORLD);
    */
    send_counts[dest] = 0;
  }
}

void
header(float GiB_mem, uint64_t pe_table_len, uint64_t fill_count,
       uint64_t scratch_len, uint64_t seed, float throttle,
       uint64_t max_message_size, int RANK, int NPES )
{

#define _printf  if (RANK == 0) printf 
#define _fprintf if (RANK == 0) fprintf


  _printf("*********************************************************\n");
  //_printf("* %s\n", PACKAGE_STRING);
  _printf("* GiB of mem                       = %.2f\n", GiB_mem);
  _printf("* Number of PEs                    = %d\n", g_npes);
  _printf("* Table size per PE (elements)     = %ld\n", pe_table_len);
  _printf("* Table size per PE (GiB)          = %.3f\n",
          (float) pe_table_len * sizeof(uint128) / GIBI);
  _printf("* Max bucket overhead per PE (GiB) = %.3f\n",
          (float) g_npes * max_message_size / (GIBI));
  _printf("* Initial fill count per PE        = %ld\n", fill_count);
  _printf("* Total size (elements)            = %ld\n", fill_count * g_npes);
  _printf("* Scratch size (elements)          = %ld\n", scratch_len);
  _printf("* Seed value                       = %ld\n", seed);
  _printf("* Sync throttle                    = %.2f%%\n", throttle * 100.0);
  _printf("*********************************************************\n");
}

uint64_t
key_exchange(uint128 * restrict table, uint64_t sort_count,
             uint128 * restrict scratch_table, uint64_t scratch_len,
             uint64_t bucket_len, uint64_t block_size, int RANK)
{
  static uint64_t recv_count;

  uint64_t ii;

  uint64_t *send_counts;

  uint128 *send_buckets;

  // calloc
  send_counts = calloc(g_npes, sizeof(uint64_t));
  send_buckets = calloc(g_npes, bucket_len * sizeof(uint128));

  //memset( send_counts,  0, g_npes * sizeof( uint64_t ) );
  //memset( send_buckets, 0, g_npes * bucket_len * sizeof(uint128) );

  recv_count = 0;

  //BARRIER();                    // Correctness of recv_count
  MPI_Barrier(MPI_COMM_WORLD);

  // Start draining to buckets at the end of scratch space.  Other PEs 
  // will begin to fill at the start of table concurrently.
  for (ii = scratch_len; ii < sort_count; ii++) {
    uint64_t dest = PE_OF(table[ii].as64.word[0]);
    uint64_t count = send_counts[dest];
    send_buckets[dest * bucket_len + count] = table[ii];
    send_counts[dest]++;
    if (send_counts[dest] >= bucket_len) {
      put_buckets(dest, send_buckets, send_counts, bucket_len, table,
                  (uint64_t *) &recv_count);
    }

    if (ii % block_size == 0)
      //BARRIER();
      MPI_Barrier(MPI_COMM_WORLD);
  }
printf("exchange 1 \n");

  // Take care of the keys saved in the scratch space
  for (ii = 0; ii < scratch_len; ii++) {
    uint64_t dest = PE_OF(scratch_table[ii].as64.word[0]);
    uint64_t count = send_counts[dest];
    send_buckets[dest * bucket_len + count] = scratch_table[ii];
    send_counts[dest]++;
    if (send_counts[dest] >= bucket_len) {
      put_buckets(dest, send_buckets, send_counts, bucket_len, table,
                  (uint64_t *) &recv_count);
    }

    if (ii % block_size == 0)
      //BARRIER();
      MPI_Barrier(MPI_COMM_WORLD);
  }
printf("exchange 2 \n");
  // Flush the remaining partially filled buckets.
  for (ii = 0; ii < g_npes; ii++) {
    int dest = (ii + RANK) % g_npes;
    put_buckets(dest, send_buckets, send_counts, bucket_len, table,
                (uint64_t *) &recv_count);
  }

  free(send_counts);
  free(send_buckets);

  //BARRIER();                    // Correctness of recv_count 
  MPI_Barrier(MPI_COMM_WORLD);

  return recv_count;
}

int
main(int argc, char *argv[])
{
  int ii;
  int retval = 0;               ///< Value to return to the shell.
  int req_arg_index = 1;        ///< Index in argv of the required argument
  int min_msg_size;
  int max_msg_size;
  
  float GiB_mem;
  float fill_fraction = 0.9;

  uint64_t seed = 0;            ///< Seed for random number generator.
  uint64_t pe_table_len = 0;    ///< Number of 128 bit records in A.
  uint64_t fill_count;
  uint64_t block_size;
  uint64_t bucket_len;
  uint64_t scratch_len;

  timer t_exch;
  timer t_sort;
  timer t_init;
  
  uint64_t recv_count;

  uint128 *scratch_table;
  uint128 *table;

                  /* d  q  t  h  v */
  opts_t options = { 0, 0, 0, 0, 0,
                     0,         // time init
                     0,         // do final sort 
                     0.9,       // Throttle
                     1,         // c
                     4096,      // m
                     1048576,   // M
                     2.0,       // step
                     1L << 24,  // w
                     0,         // s
                     1 };       // i

  timer_clear(&t_init);
  timer_start(&t_init);
  //start_pes(0);
  MPI_Init(&argc, &argv);

  int id_int;
  int W_int;
  MPI_Comm_rank(MPI_COMM_WORLD, &id_int);
  MPI_Comm_size(MPI_COMM_WORLD, &W_int);
  #define RANK id_int
  #define NPES W_int

#define _printf  if (RANK == 0) printf 
#define _fprintf if (RANK == 0) fprintf

  //BARRIER();
  MPI_Barrier(MPI_COMM_WORLD);
  timer_stop(&t_init);
  
  timer_clear(&t_exch);
  timer_clear(&t_sort);

  if (!NPES) {
    fprintf(stderr, "Error determining the number of PEs. Compile with "
            "either -DCRAY or -DSGI\n");
    exit(1);
  }

  g_npes = NPES;
  g_pe_bits = (int) ceil((log(g_npes) / log(2)));

  if (argc < 2) {
    usage(RANK, NPES);
    exit(1);
  }

  req_arg_index = parse_opts(argc, argv, &options, RANK, NPES);
  
 
  if (g_npes == 1) {
    fprintf(stderr, "ERROR: Run with more than a single PE.\n");
    exit(1);
  }

  if (req_arg_index == argc) {
    _fprintf(stderr, "ERROR: A memory size is required.\n");
    usage(RANK, NPES);
    exit(1);
  }

  if (options.opt_sync_throttle <= 0.0000) {
    _fprintf(stderr, "ERROR: Throttle setting must be greater than 0.\n");
    exit(1);
  }

  if (options.opt_min_message_size < 16) {
    _fprintf(stderr,
             "WARNING: Smallest message size is 16 bytes (setting to 16).\n");
    options.opt_min_message_size = 16;
  }

  if (options.opt_min_message_size > options.opt_max_message_size) {
    _fprintf(stderr,
             "ERROR: Max message size must be greater that the minimum.\n");
    exit(1);
  }

  if (options.opt_step <= 1.0) {
    _fprintf(stderr, "ERROR: Step size must be > 1.0\n");
    exit(1);
  }
  
  if (options.opt_seed) {
    seed = options.opt_seed;
  }
  else {
    if (RANK == 0) {
      seed = time(0);
    }
    mpp_broadcast_long((long *) &seed, 0, RANK, NPES);
  }

  //BARRIER();
  MPI_Barrier(MPI_COMM_WORLD);

  srand64((uint64_t) seed + RANK);

  GiB_mem = atof(argv[req_arg_index]);
  pe_table_len = GiB_mem * GIBI / (g_npes * sizeof(uint128));
  fill_count = pe_table_len * fill_fraction;

  scratch_len = options.opt_scratch_size / sizeof(uint128);


  if (!options.flag_quiet) {
    header(GiB_mem, pe_table_len, fill_count, scratch_len, seed,
           options.opt_sync_throttle, options.opt_max_message_size, RANK, NPES);
  }
  
  if (options.flag_time_init) {
    _printf("Init time: %.4f seconds\n", t_init.accum_wall);
  }
 
  if (scratch_len > fill_count) {
    _fprintf(stderr, "ERROR: Fill count must be greater than scratch size.\n");
    exit(1);
  }

  if (options.flag_dry_run) {
    goto cleanup;
  }

  if (!options.flag_quiet) {
    _printf("%10s\t%5s\t%23s\t%9s\n", "Message", "Wall",
            "B/W per PE (MiB/s)", "(GiB/s)");
    _printf("%10s\t%5s\t%7s\t%7s\t%7s\t%9s\n", "size (B)",
            "time", "MIN", "AVG", "MAX", "TOTAL");

  }

  scratch_table = malloc(scratch_len * sizeof(uint128));
  table = shmalloc_safe(pe_table_len * sizeof(uint128), __FILE__, __LINE__, RANK, NPES);

  /*
   * block_size is the "throttle".  The key exchange is syncrhonized every
   * block, so as block_size --> 1 more barrier calls are inserted.
   */
  block_size = scratch_len * options.opt_sync_throttle;

  uint64_t m;
  min_msg_size = options.opt_min_message_size;
  max_msg_size = options.opt_max_message_size;
printf("Rank = %u, min= %u, max = %u, step = %f \n", RANK, min_msg_size, max_msg_size, options.opt_step);

  for (m = min_msg_size; m <= max_msg_size; m *= options.opt_step) {
    for (ii = 0; ii < options.opt_iterations; ii++) {
      timer_clear(&t_exch);
      randomize(table, fill_count, ii);

      /*
       * Copy scratch_len keys to the scratch space so that we can start 
       * filling in exchanged values at the beginning of table.  The more
       * asynchronous we are the larger scratch space we'll need to ensure that
       * no PE drains from table at a rate lower than it is being filled by the
       * other PEs. 
       */
      timer_start(&t_exch);
      memcpy(scratch_table, table, scratch_len * sizeof(uint128));

      //BARRIER();
      MPI_Barrier(MPI_COMM_WORLD);

      bucket_len = m / sizeof(uint128);
printf("start of key exchange \n");
      recv_count = key_exchange(table, fill_count, scratch_table, scratch_len,
                                bucket_len, block_size, RANK);
      timer_stop(&t_exch);
printf("end of key exchange \n");
      if (options.opt_check_level >= 2) {
        verify_exchange(table, fill_count, recv_count, options.flag_quiet, RANK, NPES);
      }

      {
        double my_bandwidth = fill_count * sizeof(uint128) / t_exch.accum_wall;
        double avg_bandwidth = mpp_accum_double(my_bandwidth, RANK, NPES) / g_npes;
        double min_bandwidth = mpp_min_double(my_bandwidth, RANK, NPES);
        double max_bandwidth = mpp_max_double(my_bandwidth, RANK, NPES);
        double time_average = mpp_accum_double(t_exch.accum_wall, RANK, NPES) / g_npes;

        _printf("%10ld\t%5.1f\t%7.1f\t%7.1f\t%7.1f\t%9.3f\n",
                bucket_len * sizeof(uint128), time_average,
                min_bandwidth / MEBI, avg_bandwidth / MEBI,
                max_bandwidth / MEBI, avg_bandwidth * g_npes / GIBI);
      }
    }
    //BARRIER();  // All PEs start the next iteration together.
    MPI_Barrier(MPI_COMM_WORLD);
  }

  if (options.opt_check_level >= 1) {
    if (!options.flag_quiet) {
      _printf("Verification...\n");
    }
    verify_exchange(table, fill_count, recv_count, options.flag_quiet, RANK, NPES);
  }

  if (options.flag_final_sort) {
    double sort_time_average;
    if (!options.flag_quiet) {
      _printf("Sorting...\n");
    }

    timer_start(&t_sort);
    msb_radix_quicksort128(table, recv_count, 0, 127 - g_pe_bits);
    timer_stop(&t_sort);
    sort_time_average = mpp_accum_double(t_sort.accum_wall, RANK, NPES) / g_npes;
    //BARRIER();
    MPI_Barrier(MPI_COMM_WORLD);
    _printf("Local sort (average): %5.2fs (wall), " 
        "%6.2f MiQuadwords/s per PE\n",
         sort_time_average, fill_count / (sort_time_average * MEBI));
  }

#ifdef DEBUG
  int fd;
  char fname[20];
  sprintf(fname, "%d.out", RANK);
  fd = open(fname, O_WRONLY | O_CREAT | O_TRUNC, 00600);
  write(fd, table, recv_count * sizeof(uint128));
  close(fd);
#endif

  //shfree(table);
  free(table);

  free(scratch_table);

cleanup:
#ifndef _SGI_SOURCE
  //shmem_finalize();
  MPI_Finalize();
#endif

  return (retval);
}
