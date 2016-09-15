#ifndef __NETSTRESS_H
#define __NETSTRESS_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <limits.h>
#include <sys/time.h>
#include <string.h>
#include <unistd.h>
//#include <math.h>

#include "rand64.h"
#include "netstress_utils.h"
#include "timer.h"
//#include "shmem_utils.h"

#define _XOPEN_SOURCE 500
#define STRSIZE        32

#define KILOBYTE (1L << 10)
#define MEGABYTE (1L << 20)
#define GIGABYTE (1L << 30)

#ifdef _SHMEM
#include <mpp/shmem.h>
#define barrier_all()   shmem_barrier_all()
#endif

#if defined(_MPI) || defined(_MPI2)
#include <mpi.h>
#define barrier_all() MPI_Barrier(MPI_COMM_WORLD)
#endif

#ifdef _UPC
#include <upc.h>
#include <upc_collective.h>
#define barrier_all() upc_barrier
#endif     

/* printing macros */
#ifndef _printf
#define _printf  if (me == 0) printf 
#endif

#ifndef _fprintf
#define _fprintf if (me == 0) fprintf 
#endif

#ifndef _group_printf
#define _group_printf if (me%pepg == 0) printf
#endif

#ifndef _group_fprintf
#define _group_fprintf if (me%pepg == 0) fprintf
#endif

#endif
