#ifndef __PARSE_OPTS_H_
#define __PARSE_OPTS_H_

#include <getopt.h>

typedef struct {
  // Flags
  int flag_dry_run;             // -d 
  int flag_quiet;               // -q 
  int flag_time_detail;         // -t 
  int flag_help;                // -h
  int flag_version;             // -v 
  int flag_time_init;           // --time-init
  int flag_final_sort;          // --do-final-sort

  // Req'd args
  double opt_sync_throttle;     // --throttle

  int opt_check_level;          // -c 
  int opt_min_message_size;     // -m
  int opt_max_message_size;     // -M
  double opt_step;              // --step
  int opt_scratch_size;         // -w 
  int opt_seed;                 // -s 
  int opt_iterations;           // -i
} opts_t;

int  parse_opts(int argc, char *argv[], opts_t *options, int RANK, int NPES);
void usage(int RANK, int NPES);
#endif
