/* netstress_utils.c */

#include "netstress_utils.h"

/* ========================================================================== */

static void usage(){
  printf("\n"
	 "NAME\n"
	 "\n"
	 "\tnetstress - network stressing benchmark\n"
	 "\n"
	 "SYNOPSIS\n"
	 "\n"
	 "\tnetstress [OPTIONS]\n"
	 "\n"
	 "DESCRIPTION\n"
	 "\n"
	 "\tNetstress performs asynchronous all-to-alls using either puts,\n"
	 "\tgets, or both, with iteratively doubling message sizes.\n"
	 "\n"
	 "OPTIONS\n\n"
	 "\t-h, --help\n"
	 "\t\tDisplay this menu.\n\n"
	 "\t-a <get|put>, --all-to-all=<get|put>\n"
	 "\t\tAll-to-all type. get or put. If the option is not used, the\n"
	 "\t\tdefault behavior is to both get and put.\n\n"
	 "\t-d, --dry-run\n"
	 "\t\tPrint the header and then exit.\n\n"
	 "\t-g INT, --groups=INT\n"
	 "\t\tPartition PEs into GRPS groups. Default is 1.\n\n"
	 "\t-i INT, --iterations=INT\n"
	 "\t\tRun ITER loop iterations for each message size.\n"
	 "\t\tDefault is 30.\n\n"
	 "\t-m SIZE, --min-size=SIZE\n"
	 "\t\tBegin looping with buffer size SIZE. Default is 8 bytes.\n\n"
	 "\t-M SIZE, --max-size=SIZE\n"
	 "\t\tLoop until reaching a buffer of size SIZE. Default is 8 MB.\n\n"
	 "\t-q, --quiet\n"
	 "\t\tQuiet output. Do not print group stats.\n\n"
	 "\t-Q, --very-quiet\n"
	 "\t\tJust print summary stats and a header suitable for GNUplot.\n\n"
	 //	 "\t-s INT, --seed=INT\n"
	 //	 "\t\tUse SEED as the random number seed.\n\n"
	 "\t-t DOUBLE --timeout=DOUBLE\n"
	 "\t\tRun until this number of seconds is reached, then abort.\n"
	 "\t\tThe abort will wait until current iteration completes. Note\n"
	 "\t\tthat there is a barrier between each iteration when running\n"
	 "\t\tin this mode.\n\n"
	 "\t-v, --version\n"
	 "\t\tPrint version number.\n"
	 "\n");
}

static long long parse_size(char *optarg)
{
  int j;
  char units[8] = "KMGTPEZY";
  double size = atof(optarg);
  
  for (j = 0; j < strlen(units); j++) {
    size *= 1024.0;
    if (strchr(optarg, units[j]) || strchr(optarg, tolower(units[j])))
      return (long long)size;
  }
  if (j == strlen(units))
    return (long long) atof(optarg);

  return -1;
}

// Command line options handling
void handle_options(int argc, char *argv[], int rank, bench_t *bench){

  while (1){
    static struct option long_options[] =
      {
	{"all-to-all", required_argument, 0, 'a'},
	{"dry-run",    no_argument,       0, 'd'},
	{"groups",     required_argument, 0, 'g'},
	{"help",       no_argument,       0, 'h'},
	{"iterations", required_argument, 0, 'i'},
	{"min-size",   required_argument, 0, 'm'},
	{"max-size",   required_argument, 0, 'M'},
	{"quiet",      no_argument,       0, 'q'},
	{"very-quiet", no_argument,       0, 'Q'},
	{"seed",       required_argument, 0, 's'},
	{"timeout",    required_argument, 0, 't'},
	{"version",    required_argument, 0, 'v'},
	{0, 0, 0, 0}
      };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    
    int c = getopt_long (argc, argv, "hi:m:M:g:s:qQda:vt:",
			 long_options, &option_index);
    
    /* Detect the end of the options. */
    if (c == -1)
      break;
    
    switch (c)
      {
      case 'a':
	if (strcmp(optarg, "get") == 0)
	  bench->alltoall = 0;
	else if (strcmp(optarg, "put") == 0)
	  bench->alltoall = 1;
	else
	  bench->alltoall = 3; // error
	break;

      case 'd':
	bench->dry_run = 1;
	break;

      case 'g':
	bench->ngrps = atoi(optarg);
	break;

      case 'h':
	if (rank == 0) 
	  usage();
	barrier_all();
	exit(0);

      case 'i':
	bench->iter = atoi(optarg);
	if(bench->iter < 1)
	  bench->iter = INT_MAX;
	break;

      case 'm': 
	bench->min_size = parse_size(optarg);
       	break; 
	
      case 'M': 
	bench->max_size = parse_size(optarg);
       	break;        

      case 'q':
	bench->quiet = 1;
	break;
	
      case 'Q':
	bench->very_quiet = 1;
	break;

      case 's':
	bench->seed = atoi(optarg);
	break;    

      case 't':
	bench->timeout = atof(optarg);
	break;  

      case 'v':
	if (rank == 0) 
	  printf("\n  Netstress v%s\n\n", NETSTRESS_VERSION);
	barrier_all();
	exit(0);

      case '?':
	/* getopt_long already printed an error message. */
	if (rank == 0) 
	  usage();
	barrier_all();
	exit(1);
	break;
	
      default:
	abort ();
      }
  }
}

/* ========================================================================== */

// Check if groups of PEs are located on a single node
int check_collocation(int me, int ppg)
{
  int flag = 0, i;
  char *currhost;  

#ifdef _MPI
  char *srchost;
  MPI_Status sts;
  MPI_Request req;  
  MPI_Comm grp_comm;
  int new_me = me%ppg;
  int my_grp = me/ppg;
  char prevhost[STRSIZE]; 

  MPI_Comm_split(MPI_COMM_WORLD, my_grp, me, &grp_comm);

  currhost = malloc(STRSIZE*sizeof(char));
  srchost  = malloc(STRSIZE*sizeof(char));
  memset(srchost, 0, STRSIZE);

  if (gethostname(srchost, sizeof(srchost)) == 1) {
    strcpy(srchost, "host_error");
  }
  
  MPI_Barrier(grp_comm);
  strcpy(prevhost, srchost);
  
  for (i = 0; i < ppg; i++) { 
    MPI_Barrier(grp_comm);
    MPI_Isend(srchost, STRSIZE, MPI_CHAR,(ppg+new_me-i)%ppg, i, grp_comm, &req);
    MPI_Recv(currhost, STRSIZE, MPI_CHAR,(new_me+i) % ppg, i, grp_comm, &sts);
    MPI_Wait(&req, &sts);
   
    if (strncmp(currhost, prevhost, strlen(currhost)) != 0) { flag = 1; }

    strncpy(prevhost, currhost, STRSIZE);
  }
  MPI_Barrier(grp_comm);
  
  free(currhost);
  free(srchost);
#endif

#ifdef _SHMEM
  char *srchost;
  int size;
  currhost = shmalloc(STRSIZE*sizeof(char));
  srchost  = shmalloc(STRSIZE*sizeof(char));
  
  if (gethostname(srchost, sizeof(srchost)) == 1) {
    strcpy(srchost, "host_error");
  }

  size = strlen(srchost);

  if (me%ppg == 0) {
    char prevhost[size];
    strcpy(prevhost, srchost);
    for (i = (me/ppg)*ppg+1; i < (me/ppg)*ppg+ppg; i++) {
      shmem_getmem(currhost, srchost, size, i);
      if (strcmp(currhost, prevhost) != 0) { flag = 1; break; }
      memset(prevhost, 0, size);
      strcpy(prevhost, currhost);
    }
  }
  barrier_all();
  shfree(currhost);
  shfree(srchost);
#endif

#ifdef _UPC
  shared char * srchost;
  char *lcurrhost, *lsrchost;
  int size;
  currhost = (char*) malloc(STRSIZE*sizeof(char));
  srchost  = (shared char*) upc_all_alloc(THREADS*STRSIZE, sizeof(char));
  lsrchost  = (char *) &(srchost[MYTHREAD]);

  if (gethostname(&lsrchost[0], sizeof(&lsrchost[0])) == 1) {
    strcpy(&lsrchost[0], "host_error");
  }

  size = strlen(&lsrchost[0]);

  if (me%ppg == 0) {
    char prevhost[size];
    strcpy(prevhost, &lsrchost[0]);
    for (i = (me/ppg)*ppg+1; i < (me/ppg)*ppg+ppg; i++) {
      upc_memget(currhost, &(srchost[i]), size*sizeof(char));
      if (strcmp(currhost, prevhost) != 0) { flag = 1; break; }
      memset(prevhost, 0, size);
      strcpy(prevhost, currhost);
    }
  }
  barrier_all();
  
  free(currhost);
  if (me == 0) {
    upc_free(srchost);
  }
#endif
  
  return flag;
}

/* ========================================================================== */

// Determine order of bytes
double units_of(uint64_t bytes, char *out) {
  double partial = bytes;
  const char *units[9] = {"", "K", "M", "G", "T", "P", "E", "Z", "Y" };
  int i = 0;
  
  while ( i < 9 ) {
    if ( partial >= 1024.0 ) {
      partial /= 1024.0;
      i++;
    }
    else
      break;
  }
  
  strcpy(out, units[i]);
  return partial;
}

/* ========================================================================== */

// Make sure everyone has the same time as their initial seed
long set_seed(int me, int npes, unsigned long initial_seed)
{
#ifdef _SHMEM
  long my_seed[npes];
  long *get_seed;
  get_seed = shmalloc(npes*sizeof(long));
  get_seed[me] = initial_seed;
  barrier_all();
  shmem_get64(&get_seed[me], &get_seed[0], 1, 0);
  my_seed[me] = get_seed[me];
  barrier_all();
  shfree(get_seed);
#endif

#ifdef _UPC
  static shared long my_seed[THREADS];
  static shared long get_seed[THREADS];
  get_seed[me] = initial_seed;
  barrier_all();
  upc_all_broadcast(get_seed, &get_seed[0], sizeof(int), 
                    UPC_IN_NOSYNC | UPC_OUT_NOSYNC);
  my_seed[me] = get_seed[me];
  barrier_all();
#endif

#ifdef _MPI
  long my_seed[npes];
  long get_seed = initial_seed;
  barrier_all();
  MPI_Bcast(&get_seed, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  my_seed[me] = get_seed;
  barrier_all();
#endif

  return my_seed[me];
}

/* ========================================================================== */

