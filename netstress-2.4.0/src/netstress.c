/* netstress.c
 *
 * Compilation:
 *      netstress is designed to run in the UPC, SHMEM, and MPI environments.  
 * 
 */


#include "netstress.h"

#define INCR      2
#define MIN_SIZE  (1L<<3)  // 8 bytes
#define MAX_SIZE  (1L<<23) // 8 mebibytes

#ifdef _MPI2
#error Netstress MPI2 version is not yet implemented.
#endif

/* ========================================================================== */

#ifdef  _SHMEM
static long pSync[_SHMEM_BARRIER_SYNC_SIZE];
static long pSyncBcast[_SHMEM_BCAST_SYNC_SIZE];
#define barrier_group() shmem_barrier(firstpe, 0, pepg, pSync)
unsigned long * src;
unsigned long * dst;
#endif

#if defined(_MPI) || defined(_MPI2)
#define barrier_group() MPI_Barrier(grp_comm)
unsigned long * src;
unsigned long * dst;
#endif

#ifdef _UPC
#define barrier_group() upc_barrier
shared unsigned long * src;
shared unsigned long * dst;
#endif

/* ========================================================================== */
// Use reduction to compute min, max, and avg over PE groups

#define  _MIN(A,B)   (((A)<(B))?(A):(B))
#define  _MAX(A,B)   (((A)>(B))?(A):(B))

#ifdef  _SHMEM
static long rSync[_SHMEM_REDUCE_SYNC_SIZE];
static double Work[_MAX(4,_SHMEM_REDUCE_MIN_WRKDATA_SIZE)];
#endif

#ifdef _UPC
// UPC reduction over a group of PEs
static double all_mask_reduce(upc_op_t op, double value, double time, int me, int start, int total)
{
  int i;
  shared double *R, *G;
  double r;

  G = upc_all_alloc(THREADS, THREADS*sizeof(double)); // group
  R = upc_all_alloc(THREADS, THREADS*sizeof(double)); // reduction

  upc_barrier;
  for (i = 0; i < THREADS; i++) {
    G[i+me*THREADS] = (i>=start && i<(start+total)) ? time : value;
  }
  upc_all_prefix_reduceD(R, G, op, THREADS*THREADS, 1, NULL, UPC_IN_ALLSYNC | UPC_OUT_ALLSYNC);
  r = R[(start+total)*THREADS-1];
  
  if (op == UPC_ADD) {
    r = (r - R[start*THREADS-1])/total;
  }
  
  upc_barrier;
  if (me == 0) {
    upc_free(G);
    upc_free(R);
  }
  
  return r;
}
#endif

#ifdef _MPI
static void min_max_avg_mpi(double time, int start, int total, double stats[3], MPI_Comm comm)
{
  double t_min, t_max, t_sum;
  double time_s;

  time_s = time;
  MPI_Allreduce(&time_s, &t_min, 1, MPI_DOUBLE, MPI_MIN, comm);
  MPI_Allreduce(&time_s, &t_max, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(&time_s, &t_sum, 1, MPI_DOUBLE, MPI_SUM, comm);

  stats[0] = t_min;       // minimum time
  stats[1] = t_max;       // maximum time
  stats[2] = t_sum/total; // average time
}
#endif

#ifdef _SHMEM
static void min_max_avg_shmem(double time, int start, int total, double stats[3])
{
  static double t_min, t_max, t_sum;

  int i; 
  static double time_s;
  
  for (i = 0; i < _SHMEM_REDUCE_SYNC_SIZE; i++){
    rSync[i] = _SHMEM_SYNC_VALUE;
  }
  barrier_all();
  
  time_s = time;
  shmem_double_min_to_all(&t_min, &time_s, 1, start, 0, total, Work, rSync);
  shmem_double_max_to_all(&t_max, &time_s, 1, start, 0, total, Work, rSync);
  shmem_double_sum_to_all(&t_sum, &time_s, 1, start, 0, total, Work, rSync);
   
  stats[0] = t_min;       // minimum time
  stats[1] = t_max;       // maximum time
  stats[2] = t_sum/total; // average time
}
#endif

#ifdef _UPC
static void min_max_avg_upc(double time, int start, int total, double stats[3], int me)
{
  static double t_min, t_max, t_sum;

  int i; 
  static double time_s;
 
  t_min = all_mask_reduce(UPC_MIN, DBL_MAX, time, me, start, total);
  t_max = all_mask_reduce(UPC_MAX, 0.0,     time, me, start, total);
  t_sum = all_mask_reduce(UPC_ADD, 0.0,     time, me, start, total);
  
  stats[0] = t_min;       // minimum time
  stats[1] = t_max;       // maximum time
  stats[2] = t_sum/total; // average time
}
#endif


/* ========================================================================== */

// Bandwidth of a single all-to-all
static double bandwidth(int iter, int pepg, int numints, double time)
{
  return (iter*pepg*sizeof(unsigned long)*numints)/(time*MEGABYTE);
}

int main(int argc, char **argv) {
  int i, cnt;
  int next, numints;

  int me;               // rank of PE
  int npes;             // total number of PEs
  int pepg;             // number of PEs per group
  int firstpe;          // first PE in group

#ifdef DEBUG
  unsigned long my_checksum;
  unsigned long recvd_checksum;
  int j;
#endif
                /* m,        M,        i,  s, g, q, Q, d, a  t*/
  bench_t bench = {MIN_SIZE, MAX_SIZE, 30, 0, 1, 0, 0, 0, 2, 0.0}; 

  double size;
  int fact;
  char order[STRSIZE], hostname[STRSIZE], langname[STRSIZE];

  static double all_stats[3];
  static double grp_stats[3];

  // Declarations
#ifdef _SHMEM
  unsigned long *lsrc, *ldst;
  timer init;
  strcpy(langname, "SHMEM");
#endif
#ifdef _UPC  
  static shared timer init;
  unsigned long *lsrc, *ldst;
  strcpy(langname, "UPC");
#endif
#if defined(_MPI) || defined(_MPI2)
  timer init;
  unsigned long *lsrc, *ldst;
  MPI_Request req; MPI_Status sts;
  strcpy(langname, "MPI");
#endif

  // Initialization
  timer_clear((timer*)&init);
  timer_start((timer*)&init);
#ifdef _SHMEM
  start_pes(0);
#if defined(CRAY) || defined (_CRAYC)
  npes = shmem_n_pes();
  me   = shmem_my_pe();
#else
  me = _my_pe();
  npes = _num_pes();
#endif
  for (i = 0; i < _SHMEM_BARRIER_SYNC_SIZE; i++) {
    pSync[i] = _SHMEM_SYNC_VALUE;
  }
  for (i = 0; i < _SHMEM_BCAST_SYNC_SIZE; i++) {
    pSyncBcast[i] = _SHMEM_SYNC_VALUE;
  }
  barrier_all();
  timer t[npes];
  double start_time[npes];
#endif
  
#ifdef _UPC
  me   = MYTHREAD;
  npes = THREADS;
  static shared timer t[THREADS];
  static shared double start_time[THREADS];
#endif
 
#if defined(_MPI) || defined(_MPI2)
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  timer t[npes];
  double start_time[npes];
#endif
  timer_stop((timer*)&init);  
  _printf("# Total initialization time: %.2fs\n\n", init.accum_wall);
  start_time[0] = 0.0;

  bench.seed = set_seed(me, npes, (unsigned long)time(NULL));

  handle_options(argc, argv, me, &bench);

  barrier_all();

  // Everyone seeds their generators with the same value, which will
  // allow each PE to calculate a local checksum to compare directly
  // with the checksum of the incoming data.
  srand64(bench.seed);
  
  // Error checks
  if (npes % bench.ngrps) {
    _fprintf(stderr, "ERROR: number of groups must divide number of PEs.\n");
    exit(1);
  } else if (bench.min_size < MIN_SIZE) {
    _fprintf(stderr, "ERROR: min_size must be at least %ld bytes.\n", MIN_SIZE);
    exit(1);
  } else if (bench.max_size < bench.min_size) {
    _fprintf(stderr, "ERROR: min_size must not exceed max_size.\n");
    exit(1);
  } else if (bench.alltoall > 2) {
    _fprintf(stderr, "ERROR: all-to-all must be either \'get\' or \'put\'.\n");
    exit(1);
  }

  pepg = npes/bench.ngrps;
  firstpe = (me/pepg)*pepg;

  if (bench.very_quiet || bench.ngrps == 1) {
    bench.quiet = 1;
  }

  // Print header
  if (!bench.very_quiet) {
    _printf("===================================\n");
    _printf(" Netstress-%s\n"          , langname);
    _printf(" All-to-all type  = ");
    switch (bench.alltoall) {
    case 0:
      _printf("%10s\n", "get"); break;
    case 1:
      _printf("%10s\n", "put"); break;
    default:
      _printf("%10s\n", "put & get"); break;
    }
    _printf(" Total PEs        = %10d\n", npes);
    _printf(" Number of groups = %10d\n", bench.ngrps);
    _printf(" PEs per group    = %10d\n", npes/bench.ngrps);
    size = units_of(bench.min_size, order);
    _printf(" Min message size = %10.1lf %sB\n", size, order);
    size = units_of(bench.max_size, order);
    _printf(" Max message size = %10.1lf %sB\n", size, order);
    _printf(" Loop iterations  = %10d\n", bench.iter);
    //    _printf(" Seed value       = %10ld\n", bench.seed);
#ifdef DEBUG
      _printf(" Data verification enabled\n");
#endif
    _printf("===================================\n");
  }
  fflush(stdout);

  if (bench.dry_run) {
    exit(1);
  }

  if (!bench.quiet && !bench.very_quiet) {
    _printf("%6s %12s %12s\n", "Group", "PEs", "Host");
    fflush(stdout);

    if (gethostname(hostname, sizeof(hostname)) == 1) {
      strcpy(hostname, "host_error");
    }
    
    if (check_collocation(me, pepg)) { // Not currently implemented in MPI
      strcpy(hostname, "multiple\0");
    }

    barrier_all();
    _group_printf("%5d: %4d -> %4d %12s\n",
		  me/pepg, firstpe, firstpe+pepg-1, hostname);
    fflush(stdout);
    
/* #ifdef DEBUG */
/*       barrier_all(); // ensure header has finished printing */
/*       printf("Rank %4d running on host: %s\n", me, hostname); */
/* #endif */
  }

  if (bench.very_quiet) {
    _printf("# Number of PEs: %d --- PEs per Group: %d\n", npes, pepg);
  }

  // Allocate message buffer
#ifdef _SHMEM
  src = shmalloc(bench.max_size);
  dst = shmalloc(bench.max_size);
  lsrc = &src[0];
  ldst = &dst[0];
#endif

#ifdef _UPC
  src = upc_all_alloc(THREADS, bench.max_size);
  dst = upc_all_alloc(THREADS, bench.max_size);
  ldst = (unsigned long *) &(dst[MYTHREAD]);
  lsrc = (unsigned long *) &(src[MYTHREAD]);
#endif

#if defined(_MPI) || defined(_MPI2)
  src  = malloc(bench.max_size);
  dst  = malloc(bench.max_size);
  lsrc = &src[0];
  ldst = &dst[0];
#endif

  if (!src || !dst) {
    _fprintf(stderr, "Allocation failure.\n");
    exit(1);
  }

  // Fill message buffer with random 64-bit integers
  for (i = 0; i < bench.max_size/sizeof(unsigned long); i++) {
    ldst[i] = rand64();
    lsrc[i] = rand64();
  }

  if (!bench.very_quiet) {
    _printf("\n%6s %13s %10s %10s %10s %10s \n", "Group", "Size",
	   "Max Time", "Min B/W", "Avg B/W", "Max B/W");
    _printf("%6s %13s %10s %10s %10s %10s\n", "", "       ",
	   "(sec)", "(MB/s)", "(MB/s)", "(MB/s)");
  }
  barrier_all();
  fflush(stdout);

#ifdef _MPI
  MPI_Comm grp_comm;
  int new_me = me%pepg;
  int my_grp = me/pepg;

  firstpe = 0; // for MPI each group has PEs 0 thru pepg-1

  MPI_Comm_split(MPI_COMM_WORLD, my_grp, me, &grp_comm);
#endif
  
  // Warm up the interconnect by performing an all-to-all
  // 1000 times with a single integer.
  for (cnt = 1; cnt <= 1000; cnt++) {
    for (i = firstpe; i < firstpe+pepg; i++) {
      next = ((me+i) % pepg) + firstpe;
#ifdef _SHMEM
      shmem_put64(dst, src, MIN_SIZE, next);
#endif
      
#ifdef _UPC
      upc_memput(&(dst[next]), lsrc, MIN_SIZE*sizeof(unsigned long));
#endif

#ifdef _MPI
      next %= pepg;
      MPI_Irecv(dst, MIN_SIZE * sizeof(unsigned long), MPI_BYTE,
		(pepg+new_me-i)%pepg+firstpe, i, grp_comm, &req);
      MPI_Send(src, MIN_SIZE * sizeof(unsigned long), MPI_BYTE,
	       next, i, grp_comm);
      MPI_Wait(&req, &sts);
#endif

      //barrier_group();
    }
  }
  barrier_all();
  
  /***********************/
  /*** BEGIN stressing ***/
  /***********************/
  start_time[me] = wall();
  for (numints = bench.min_size/sizeof(unsigned long); numints <= (bench.max_size/sizeof(unsigned long)); numints *= INCR) {

    if (me%pepg == 0) {
      if (bench.very_quiet) {
	size = sizeof(unsigned long)*numints;
	strcpy(order, "");
      } else {
	size = units_of(sizeof(unsigned long)*numints, order);
      }
    }
    _printf("%5s: %10ld%2sB ", "all", (unsigned long)size, order);
    fflush(stdout);
    barrier_all();

    timer_clear((timer*)&t[me]);
    timer_start((timer*)&t[me]);
    if (bench.alltoall == 0 || bench.alltoall == 2) { // begin get
      for (cnt = 1; cnt <= bench.iter; cnt++) {
	for (i = firstpe; i < firstpe+pepg; i++) {
	  next = ((me+i) % pepg) + firstpe;
#ifdef _MPI
	  next %= pepg;
#endif
	  
#ifdef DEBUG
	  if (i == firstpe) {
	    my_checksum    = 0;
	    recvd_checksum = 0;
	  }
	  for (j = numints/INCR; j < numints; j++) {
	    my_checksum ^= lsrc[j];
	  }
	  barrier_group();
#endif

	  
#ifdef _SHMEM
	  shmem_get64(dst, src, numints, next);
#endif
	  
#ifdef _UPC
	  upc_memget(ldst, &(src[next]), numints*sizeof(unsigned long));
#endif

#ifdef _MPI
	MPI_Isend(src, numints * sizeof(unsigned long), MPI_BYTE,
		  ((pepg+new_me-i)%pepg)+firstpe, i, grp_comm, &req);
	MPI_Recv(dst, numints * sizeof(unsigned long), MPI_BYTE,
		 next, i, grp_comm, &sts);
	MPI_Wait(&req, &sts);
       // to generate traffic matrix
       // _printf("%6s %6s %10s \n", "Src", "Dst", "Size");
       // _group_printf("%5d %10ld %10.4lf \n",
	//      	    me/pepg, ((pepg+new_me-i)%pepg)+firstpe, numints * sizeof(unsigned long));
       // fflush(stdout);
#endif

#ifdef DEBUG
	  barrier_group();
	  
	  for (j = numints/INCR; j < numints; j++) {
	    recvd_checksum ^= ldst[j];
	  }
	    
	  if (my_checksum != recvd_checksum ) {
	    fprintf(stderr,
		    "\nERROR: rank %d (%s) get from rank %d checksum failure, "
		    "size: %d.\n",  me, hostname, next,
		    numints*sizeof(unsigned long));
	  }
#endif
	  
	  //  barrier_group();
	} // end single get
	// Check for stonewall
	if ((int)bench.timeout > 0) {
	  barrier_all(); // need to barrier so everyone has the same time
	  if (wall() - start_time[me] >= bench.timeout) {
	    _printf("\n# Time limit of %.2fs reached.\n", bench.timeout);
	    goto cleanup;
	  }
	}
      } // end iter for get
    } // end all gets
    
    if (bench.alltoall == 2) {
      barrier_group();
    }
    
    if (bench.alltoall == 1 || bench.alltoall == 2) { // begin put   
      for (cnt = 1; cnt <= bench.iter; cnt++) {
	for (i = firstpe; i < firstpe+pepg; i++) {
	  next = ((me+i) % pepg) + firstpe;               
#ifdef _MPI
	  next %= pepg;
#endif

#ifdef DEBUG	  
	  if (i == firstpe) {
	    my_checksum    = 0;
	    recvd_checksum = 0;
	  }
	  for (j = numints/INCR; j < numints; j++ ) {
	    my_checksum ^= lsrc[j];
	  }
	  barrier_group();
#endif
	  
	  
#ifdef _SHMEM
	  shmem_put64(dst, src, numints, next);
#endif
	  
#ifdef _UPC
	  //#pragma _CRI pgas_defer_sync;
	  upc_memput(&(dst[next]), lsrc, numints*sizeof(unsigned long));
	  //	  upc_fence;
#endif	 

#ifdef _MPI
	MPI_Irecv(dst, numints * sizeof(unsigned long), MPI_BYTE,
		  (pepg+new_me-i)%pepg+firstpe, i, grp_comm, &req);
	MPI_Send(src, numints * sizeof(unsigned long), MPI_BYTE,
		 next, i, grp_comm);
	MPI_Wait(&req, &sts);
#endif

#ifdef DEBUG 
	  barrier_group();
	  
	  for (j = numints/INCR; j < numints; j++) {
	    recvd_checksum ^= ldst[j];
	  }
	  
	  if (my_checksum != recvd_checksum ) {
	    fprintf(stderr,
		    "\nERROR: rank %d (%s) get from rank %d checksum failure, "
		    "size: %d.\n",  me, hostname, next, 
		    numints*sizeof(unsigned long));
	  }
#endif
	  
	  // barrier_group();
	} // end single put
	// Check for stonewall
	if ((int)bench.timeout > 0) {
	  barrier_all(); // need to barrier so everyone has the same time
	  if (wall() - start_time[me] >= bench.timeout) {
	    _printf("\n# Time limit of %.2fs reached.\n", bench.timeout);
	    goto cleanup;
	  }
	}
      } // end iter for put      
    } // end all puts
    
    timer_stop((timer*)&t[me]);
    barrier_all();
    
    /***********************/
    /***  END stressing  ***/
    /***********************/

    fact = (bench.alltoall == 2) ? 2 : 1; // number of all-to-alls
    
    // Report overall stats
#ifdef _MPI
    min_max_avg_mpi(t[me].accum_wall, 0, npes, all_stats, MPI_COMM_WORLD);
#endif
#ifdef _SHMEM  
    min_max_avg_shmem(t[me].accum_wall, 0, npes, all_stats);
#endif
#ifdef _UPC
    min_max_avg_upc(t[me].accum_wall, 0, npes, all_stats, me);
#endif
    
    _printf("%10.4lf %10.2lf %10.2lf %10.2lf \n",
	    all_stats[1],
	    fact*bandwidth(bench.iter, pepg, numints, all_stats[1]),
	    fact*bandwidth(bench.iter, pepg, numints, all_stats[2]),
	    fact*bandwidth(bench.iter, pepg, numints, all_stats[0]) );
    fflush(stdout);
    barrier_all();
    
    // Report group stats
    if (!bench.quiet) {
      barrier_group();
#ifdef _MPI
      min_max_avg_mpi(t[me].accum_wall, firstpe, pepg, grp_stats, grp_comm);
#endif
#ifdef _SHMEM  
      min_max_avg_shmem(t[me].accum_wall, firstpe, pepg, grp_stats);
#endif
#ifdef _UPC
      min_max_avg_upc(t[me].accum_wall, firstpe, pepg, grp_stats, me);
#endif
      _group_printf("%5d: %10ld%2sB %10.4lf %10.2lf %10.2lf %10.2lf \n",
		    me/pepg, (int)size, order, grp_stats[1],
		    fact*bandwidth(bench.iter, pepg, numints, grp_stats[1]),
		    fact*bandwidth(bench.iter, pepg, numints, grp_stats[2]),
		    fact*bandwidth(bench.iter, pepg, numints, grp_stats[0]) );
      fflush(stdout);
      barrier_all();
      _printf("\n");
      fflush(stdout);
    }
    
  }
  
 cleanup:
  if ((int)start_time[0] > 0)
    _printf("\n# Total run time: %.2fs\n", wall() - start_time[me]);
  barrier_all();

#ifdef _SHMEM
  shfree(dst);
  shfree(src);
#if defined(CRAY) || defined(_CRAYC)
  shmem_finalize();
#endif
#endif

#ifdef _UPC
  if (me == 0) {
    upc_free(src);
    upc_free(dst);
  }
#endif

#if defined ( _MPI ) || defined ( _MPI2 )
  MPI_Finalize();
  free(src);
  free(dst);
#endif

  return(0);
}

