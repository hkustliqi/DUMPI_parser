/*****************************************************************
 *
 * CONFIDENCE
 *
 * Variability Metrics for System Interconnects
 *
 * Version: 2010.06.18
 *
 * Copyright (c) 2008, 2009, 2010 UT-Battelle, LLC.
 *
 * All Rights Reserved. See the file 'LICENSING.pdf' for the
 * License Agreement.
 *
 * Author: Jeffery A. Kuehn
 *         Extreme Scale Systems Center
 *         Oak Ridge National Laboratory 
 *         PO Box 2008 MS6164
 *         Oak Ridge TN 37831-6164
 *         kuehn@ornl.gov
 *
 * The author also wishes to acknowledge and thank the following
 * collaborators for our many conversations and their critiques,
 * questions, and suggestions:
 *
 *         Stephen Poole
 *         Stephen Hodson
 *         Blair Sullivan
 *         Christopher Groer
 *         Richard Barrett
 *         [PARTNER_STAFF_1]
 *         [PARTNER_STAFF_2]
 *
 * Oak Ridge National Laboratory is operated by UT-Battelle, LLC
 * for the U.S. Department of Energy.
 *
 * This work was supported by the United States Department of
 * Defense and used resources of the Extreme Scale Systems Center
 * at Oak Ridge National Laboratory.
 *
 ****************************************************************/
#include <math.h>
#ifndef HAVE_CONFIDENCE_H
#define HAVE_CONFIDENCE_H

/**************************************************************
 * TYPES
 **************************************************************/
typedef unsigned long long ul64;

struct HISTOGRAM {
	double min0;		/* minimum (mins=1) */
	double mod0, mods;	/* mode: 0,scaled */
	double med0, meds;	/* median: 0,scaled */
	double max0, maxs;	/* maximum: 0,scaled */
	double m10, m1m, m1s;	/* 1st moment: 0,min,scaled */
	double m20, m2m, m2s;	/* 2nd moment: 0,min,scaled */
	double m30, m3m, m3s;	/* 3rd moment: 0,min,scaled */
	double m40, m4m, m4s;	/* 4th moment: 0,min,scaled */
	ul64 nsamples;		/* number of samples in the distribution */
	ul64 *dist;		/* pointer to histogram array: dist[nbins] */
	char *label;		/* label */
};

struct MEASUREMENT {		/* group communications and timer measurements */
	struct HISTOGRAM tim;	/* timer */
	struct HISTOGRAM onNodeOnesided;	/* local communication */
	struct HISTOGRAM onNodePairwise;	/* local communication */
	struct HISTOGRAM onNodeOnesidedMinimum;	/* local communication */
	struct HISTOGRAM onNodePairwiseMinimum;	/* local communication */
	struct HISTOGRAM offNodeOnesided;	/* remote communication */
	struct HISTOGRAM offNodePairwise;	/* remote communication */
	struct HISTOGRAM offNodeOnesidedMinimum;	/* remote communication */
	struct HISTOGRAM offNodePairwiseMinimum;	/* remote communication */
	int nbins;		/* number of bins in the histograms */
	double binwidth;	/* bin interval : [n*binwidth,(n+1)*binwidth( */
	int buflen;		/* message size in bytes */
	char *label;
};

struct COMMBUFF {
	void *b;		/* buffer data */
	int l;			/* buffer length in bytes */
};

typedef struct HISTOGRAM histogram_t;
typedef struct MEASUREMENT measurement_t;
typedef struct COMMBUFF commbuff_t;

typedef struct HISTOGRAM *histogram_p;
typedef struct MEASUREMENT *measurement_p;
typedef struct COMMBUFF *commbuff_p;

/**************************************************************
 * FUNCTION PROTOTYPES
 **************************************************************/
void *measurement_create(char *label, int buflen, int numbins, double binsize);
void *measurement_destroy(void *m);
void measurement_collect(void *m);
void measurement_bin(measurement_t * p, double *t, double *cos, double *cpw, int LOCAL);
void measurement_analyze(void *m, double scale);
void measurement_cryweird(void *normal, void *me);
void measurement_aggregate(void *g, void *l);
void measurement_serialize(void *m, int WRITINGPE);
void measurement_fmthist(histogram_p h, char *label);

double measurement_moment(int mom, ul64 * dist, ul64 center, ul64 nsamples, ul64 nbins, double binwidth);
ul64 measurement_samplecount(ul64 * dist, int nbins);

void comm_initialize(int *argc, char **argv[]);	/* generic interface */
void *comm_newbuffer(int nbytes);
void comm_test1(measurement_p p);
void comm_test2(measurement_p p);
void comm_test3(measurement_p p);
void comm_aggregate(measurement_p g, measurement_p l);
void comm_showmapping();
void comm_finalize();

void comm_MPI_initialize(int *argc, char **argv[]);	/* MPI internal */
void *comm_MPI_newbuffer(int nbytes);
void comm_MPI_test1(measurement_p p);
void comm_MPI_test2(measurement_p p);
void comm_MPI_test3(measurement_p p);
void comm_MPI_finalize();

void comm_SHMEM_initialize(int *argc, char **argv[]);	/* SHMEM internal */
void *comm_SHMEM_newbuffer(int nbytes);
void comm_SHMEM_test1(measurement_p p);
void comm_SHMEM_test2(measurement_p p);
void comm_SHMEM_test3(measurement_p p);
void comm_SHMEM_finalize();

ul64 comm_getnodeid();
int comm_ceil2(int n);
void getoptions(int argc, char *argv[], char **envp);
extern inline int time2bin(double t);
extern inline double bin2time(int b);
extern inline double bin2midtime(int b);


/**************************************************************
 * COMM PACKAGE MACROS
 * NAMEBUFFSIZE -- POSIX and OpenMPI both require at least 256
 *                 while UNIX requires 8-14
 * ROOTONLY     -- readabililty macro for that serial stuff
 **************************************************************/
#define NAMEBUFFSIZE 256
#define ROOTONLY   if(ThisRankID==RootRankID)

/**************************************************************
 * GLOBAL VARIABLES
 **************************************************************/
#ifndef MAIN
#define EXTERN(x) extern x
#define EXTERNINIT(x,y) extern x
#else
#define EXTERN(x) x
#define EXTERNINIT(x,y) x = y
#endif

EXTERN(int TestType);
EXTERN(int ThisRankID);
EXTERN(int RootRankID);
EXTERN(int NumRanks);
EXTERN(int NumStages);
EXTERN(int NumWarmUp);
EXTERN(int NumMessages);
EXTERN(int NumCycles);
EXTERN(ul64 NumMessagesGlobal);
EXTERN(long *NodeID);
EXTERN(int BufLen);		/* bytes */
EXTERN(int NumBins);
EXTERN(double BinSize);		/* seconds */
EXTERN(double MaxHistTime);	/* seconds */
EXTERN(double HistScale);	/* seconds */
EXTERN(int LogarithmicBinning);
EXTERN(int RankMapping);
EXTERN(char CaseName[NAMEBUFFSIZE]);
EXTERNINIT(char *COPYRIGHT,
"########################################################################\n" \
"#                                                                      #\n" \
"# CONFIDENCE --  VARIABILITY METRICS FOR SYSTEM INTERCONNECTS          #\n" \
"#                                                                      #\n" \
"# Copyright 2008, 2009, 2010 UT-Battelle, LLC.                         #\n" \
"#         This program and the results it generates represent          #\n" \
"#         unpublished proprietary research owned by UT-Battelle, LLC   #\n" \
"#         Explicit permission is required for all uses.                #\n" \
"#                                                                      #\n" \
"# Author:                                                              #\n" \
"#         Jeffery A. Kuehn                                             #\n" \
"#         Extreme Scale Systems Center                                 #\n" \
"#         Oak Ridge National Laboratory                                #\n" \
"#         kuehn@ornl.gov                                               #\n" \
"#                                                                      #\n" \
"# Acknowledgements:                                                    #\n" \
"#         This work was supported by the United States Department of   #\n" \
"#         Defense & used resources of the Extreme Scale Systems Center #\n" \
"#         at Oak Ridge National Laboratory.                            #\n" \
"#                                                                      #\n" \
"########################################################################\n" );

#ifndef MAIN
#undef EXTERN
#endif

#endif				/* HAVE_CONFIDENCE_H */
