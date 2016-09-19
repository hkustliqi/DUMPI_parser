/** 
 * parse_opts.c 
 *
 * Processes user-provided command line options. 
 * The -h option can be used to view all options.
 */

#include <unistd.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "parse_opts.h"
#include "shmem_utils.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef PACKAGE_VERSION
#error "Definition for PACKAGE_VERSION missing."
#endif

static char *optmap = "0:dqhvtnc:m:M:s:i:";

static struct option long_options[] =
{
  {"dry-run",         no_argument,       0, 'd'},
  {"quiet",           no_argument,       0, 'q'},
  {"time-detail",     no_argument,       0, 't'}, 
  {"help",            no_argument,       0, 'h'},
  {"version",         no_argument,       0, 'v'},
  {"time-init",       no_argument,       0,  0 },
  {"do-final-sort",   no_argument,       0,  0 },
        
  {"check-level",     required_argument, 0, 'c'},
  {"throttle",        required_argument, 0,  0 },
  {"max-message-size",required_argument, 0, 'M'},
  {"min-message-size",required_argument, 0, 'm'},
  {"step"            ,required_argument, 0,  0 },
     
  {"log-scratch-size",required_argument, 0,  0 },
  {"seed",            required_argument, 0, 's'},
  {"iterations",      required_argument, 0, 'i'},

  {0, 0, 0, 0}
};

void usage(int RANK, int NPES) 
{

#define _printf  if (RANK == 0) printf 
#define _fprintf if (RANK == 0) fprintf

  _fprintf(stdout,
"\n"
"NAME\n"
"\n"
"\tHAMR - Huge Asynchronous 128-bit MSB Radix Sort\n"
"\n"
"SYNOPSIS\n"
"\n"
"\thamr [OPTIONS] GiB_OF_MEMORY\n"
"\n"
"DESCRIPTION\n"
"\n"
"\tHAMR fills GiB_OF_MEMORY worth of random 128-bit unsigned integers\n"
"\tand sorts them in parallel. The sort is in-place (with an arbitrarily\n"
"\tsmall buffer used for scratch space) and will use a variety of message\n"
"\tsizes for the key exchange.  Note: it is very easy to oversubscribe\n"
"\tmemory with this benchmark (see notes for guidance).  This benchmark\n"
"\tshould perform well on systems with excellent all-to-all and atomic\n"
"\tfetch and increment performance.\n"
"\n"
"\tIf running with a number of PEs that is not a power of two, then HAMR\n"
"\twill generate the random data in such a way as to ensure that the high\n"
"\tceil(log_2(N)) bits will be in the range 0 to N-1 for N PEs.\n"
"\n"
"OPTIONS\n\n"
"\t-h --help\n"
"\t\tDisplay this help menu and then quit.\n"
"\n"
"\t-v --version\n"
"\t\tDisplay version information and then quit.\n"
"\n"
"\t-d --dry-run\n" 
"\t\tPrint out problem size and memory information, then exit.\n"
"\n"
"\t-s SEED --seed=SEED\n"
"\t\tUse SEED to seed the random number generator.\n"
"\t\tThe default seed is the Unix epoch time.\n"
"\n"
"\t-q --quiet\n"
"\t\tPrint timing information only.  Useful for redirection\n"
"\t\tto GNUplot for graphing results.\n"
"\n"
"\t--do-final-sort\n"
"\t\tPerform the final serial sort on each PE after the key\n"
"\t\texchange.  Default is no.\n"
"\n"
"\t-c LEVEL --check=LEVEL\n"
"\t\tVerification of the key exchange phase.\n"
"\t\tLEVEL=0 : No verifcation is performed.\n"
"\t\tLEVEL=1 : Verification is performed only after last\n"
"\t\t          message size. (Default)\n"
"\t\tLEVEL=2 : Verification after EVERY message size.\n"
"\n"
"\t-m SIZE --min-message-size=SIZE\n"
"\t\tStart by using messages of SIZE bytes.  SIZE may be in\n"
"\t\tthe form x^y.  Default is 1024 (1 KiB).\n"
"\n"
"\t-M SIZE --max-message-size=SIZE\n"
"\t\tUse a maximum message size of SIZE bytes.\n"
"\t\tSIZE may be in the form x^y.  Default is 2^20 (1 MiB).\n"
"\n"
"\t--step=FACTOR\n"
"\t\tMultiply by FACTOR when iterating over message sizes.\n"
"\t\tMust be > 1.0.  Default is 2.0.\n"
"\n"
"\t-i N --iterations=N\n"
"\t\tFor each message size, perform N iterations of key\n"
"\t\texchange.  New random data is generated between\n"
"\t\teach iteration.  Default is N=1.\n"
"\n"
"\t--time-init\n"
"\t\tPrint out the amount of time spent in start_pes.\n"
"\n"
"ADVANCED OPTIONS\n\n"
"\t--throttle=THROTTLE\n"
"\t\tSets the synchronization throttle position.  Synchronization\n"
"\t\thappens when each PE has transmitted a block of\n"
"\t\tWORKSPACE * THROTTLE bytes.  A value of 0 inserts the maximum\n"
"\t\tnumber of barriers.  Default is 0.9.  Values greater than 1.0\n"
"\t\tare possible, but not recommended because it can introduce\n"
"\t\terrors if there is a large difference in the running times\n"
"\t\tof the PEs. High throttle values would make it possible, for\n"
"\t\tinstance, for a particularly fast PE to send all of its data\n"
"\t\tto a slow remote PE before the remote PE has had a chance\n"
"\t\tto create adequate space in its data array. In general, large\n"
"\t\tthrottle values may require larger scratch sizes in order to\n"
"\t\tensure correctness.  The particular balance will be system\n"
"\t\tdependent.\n"
"\n"
"\t--log-scratch-size=LOG_SIZE\n"
"\t\tEach PE will allocate WORKSPACE=2^LOG_SIZE bytes of scratch\n"
"\t\tspace to use for the sort.  Default is 24.  This value controls\n"
"\t\tthe \"in-placeness\" of the sort.  As the scratch space increases\n"
"\t\tthe amount of data to sort decreases (for a fixed memory size).\n"
"\t\tLarger scratch spaces require less synchronization.\n"  
"\n"
"USAGE AND OUTPUT NOTES\n"
"\n"
"\tThis benchmark consists of two phases: an asynchronous parallel key\n"
"\texchange and a local MSB radix sort.\n\n"
"\tIn phase 1, each PE allocates GiB_OF_MEMORY on the symmetric heap, \n"
"\tthen fills 90%% of that memory with random data.  It also allocates\n"
"\tmemory for its local scratch space on the private heap, then performs\n"
"\tthe following loop:\n\n"
"\tfor each m from min-message-size to max-message-size\n\tdo\n"
"\t\t1. Allocate NPES * m bytes on the private heap to serve as\n"
"\t\t   sending buckets.\n"
"\t\t2. Randomize the data to be sorted.\n"
"\t\t3. Copy the first WORKSPACE bytes from the unsorted data to the\n"
"\t\t   scratch space.\n"
"\t\t4. Iterate over the unsorted data moving elements into the\n"
"\t\t   destination buckets.\n"
"\t\t5. When any bucket is full send it.\n"
"\t\t6. If the amount sent is a multiple of WORKSPACE * THROTTLE\n"
"\t\t   then synchronize.\n"
"\tdone\n\n"
"\tIt is easy to oversubscribe memory for larger message sizes if it is\n"
"\talso desired to fill memory with data to be sorted.  Use the -m, -M\n"
"\toptions to ensure that\n"
"\tGiB_OF_MEMORY*2^30 + NPES * NPES * message_size + WORKSPACE is less\n"
"\tthan the sum of memory on all nodes in use for all message sizes.\n"
"\tUse the --log-scratch-size option to control the size of WORKSPACE.\n"
"\n"
"\tThe Wall time column in the output consists only of the communication\n"
"\ttime and does not include any time generating the random numbers at\n"
"\teach iteration.  For large problem sizes the random data generation\n"
"\ttime will not be insignificant.\n"
 );
  //BARRIER();
  MPI_Barrier(MPI_COMM_WORLD);
} /* end: void usage() */

/* ======================================================================== */
long parse_arg(char *arg)
{
  long val = (strrchr(arg,'^')) ? 
    (long)pow(atof(arg),atof(&arg[strcspn(arg,"^")+1])) : (long)atof(arg);

  return val;
}

int parse_opts(int argc, char *argv[], opts_t *options, int RANK, int NPES)
{

#define _printf  if (RANK == 0) printf 
#define _fprintf if (RANK == 0) fprintf


  int c;
  int option_index = 0;

  while ((c=getopt_long(argc,argv,optmap,long_options,&option_index)) != -1) {

    switch (c) {
    case 0:
      if (!strcmp(long_options[option_index].name, "log-scratch-size")) {
        options->opt_scratch_size     = pow(2,atoi(optarg));
        break;
      }
      else if (!strcmp(long_options[option_index].name, "do-final-sort")) {
        options->flag_final_sort = 1;
        break;
      }
      else if (!strcmp(long_options[option_index].name, "throttle")) {
        options->opt_sync_throttle    = atof(optarg);
        break;
      }
      else if (!strcmp(long_options[option_index].name, "step")) {
        options->opt_step = atof(optarg);
        break;
      }
      else if (!strcmp(long_options[option_index].name, "time-init")) {
        options->flag_time_init = 1;
        break;
      }
      break;
      
    case 'n':  // Recognize -n for backward compatability.
      options->flag_final_sort = 0;  
      break;
      
    case 'd': 
      options->flag_dry_run    = 1;
      break;
      
    case 'q':
      options->flag_quiet      = 1;
      break;
      
    case 'c':
      options->opt_check_level = atoi(optarg);
      break;

    case 's': 
      options->opt_seed        = atoi(optarg);
      break;

    case 'i':
      options->opt_iterations  = atoi(optarg);
      break;

    case 'm':
      options->opt_min_message_size = parse_arg(optarg);
      break;

    case 'M':
      options->opt_max_message_size =  parse_arg(optarg);
      break;

    case 'v': /* version information */
      _printf("hamr version %s\n", PACKAGE_VERSION);
      exit(0);
      break;
      
    default:
      usage(RANK, NPES);
      exit(0);
      break;
    } 
  } 
  
  return optind;
}

