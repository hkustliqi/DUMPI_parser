/***      UNCLASSIFIED//FOR OFFICIAL USE ONLY      ***/

#ifndef _simple_timer_h_
#define _simple_timer_h_

#include <inttypes.h>
#include <sys/time.h>

typedef  int64_t simple_timer_t;

static int64_t simple_timer( simple_timer_t *T)
{
  struct timeval t;

  gettimeofday(&t,NULL);

  int64_t rtn = (t.tv_sec*1000000 + t.tv_usec) - *T;
  *T          = t.tv_sec*1000000 + t.tv_usec;

  return rtn;
}

#endif

/***      UNCLASSIFIED//FOR OFFICIAL USE ONLY      ***/
