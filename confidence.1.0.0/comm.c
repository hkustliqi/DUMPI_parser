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
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <errno.h>
#include <mpi.h>
#include "config.h"
#include "confidence.h"
#include "orbtimer.h"

/*
 * routines in this file abstract the communication layer.
 *
 * so far, only MPI is implemented
 *
 */

void comm_aggregate(measurement_p g, measurement_p l) {
	/* collects local measurments into a global measurement */
	int ierr;
	assert(l->nbins == g->nbins);
	ierr = 0;
#ifdef SHMEM
	/* unimplemented */
#else				/* MPI case */
	ierr += MPI_Allreduce(l->onNodeOnesided.dist,
			      g->onNodeOnesided.dist,
			      l->nbins,
			      MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD);
	ierr += MPI_Allreduce(l->onNodePairwise.dist,
			      g->onNodePairwise.dist,
			      l->nbins,
			      MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD);
	ierr += MPI_Allreduce(l->onNodeOnesidedMinimum.dist,
			      g->onNodeOnesidedMinimum.dist,
			      l->nbins,
			      MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD);
	ierr += MPI_Allreduce(l->onNodePairwiseMinimum.dist,
			      g->onNodePairwiseMinimum.dist,
			      l->nbins,
			      MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD);
	ierr += MPI_Allreduce(l->offNodeOnesided.dist,
			      g->offNodeOnesided.dist,
			      l->nbins,
			      MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD);
	ierr += MPI_Allreduce(l->offNodePairwise.dist,
			      g->offNodePairwise.dist,
			      l->nbins,
			      MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD);
	ierr += MPI_Allreduce(l->offNodeOnesidedMinimum.dist,
			      g->offNodeOnesidedMinimum.dist,
			      l->nbins,
			      MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD);
	ierr += MPI_Allreduce(l->offNodePairwiseMinimum.dist,
			      g->offNodePairwiseMinimum.dist,
			      l->nbins,
			      MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD);
	ierr += MPI_Allreduce(l->tim.dist,
			      g->tim.dist,
			      l->nbins,
			      MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD);
	assert(ierr == 0);
#endif
	return;
}

void comm_SHMEM_initialize(int *argc, char **argv[]) {
}

void comm_SHMEM_finalize() {
}

void *comm_newbuffer_SHMEM(int nbytes) {
	commbuff_p p;
	return (void *)p;
}

unsigned long long comm_getnodeid() {
	/* 
	 * Assuming the system's nodes are "numbered" after some fashion, this
	 * function returns a unique integer for each hardware 'node'.  The
	 * nodeid's of two MPI ranks are compared and the communication is
	 * assumed to be local (ie. at a latency potentially less than that of
	 * the network) if the nodeid's are the same.  If something other than
	 * digits are used to differentiate the nodes this function should be
	 * replaced with a function which generates a useful nodeid for the
	 * target machine. If none of the NODEID_<TYPE> macros are set, we
	 * punt, and return nodeid's suggesting every rank on a separate node.
	 */
	unsigned long long id;
	int tmp, ierr;
#if defined(NODEID_GETHOSTNAME) || defined(NODEID_MPI) || defined(NODEID_SLURM)
	int len;
	char *nodename;
	char namebuff[NAMEBUFFSIZE];
	char nid[NAMEBUFFSIZE];
	char *pn, *pp, *limit;
#    if defined(NODEID_MPI)
	ierr = MPI_Get_processor_name(namebuff, &len);
	nodename = namebuff;
#    elif defined(NODEID_GETHOSTNAME)
	gethostname(namebuff, NAMEBUFFSIZE);
	nodename = namebuff;
#    elif defined(NODEID_SLURM)
	nodename = getenv("SLURM_NODEID");
#    endif
	len = strlen(nodename);	/* MPI should have been okay */
	limit = nodename + strlen(nodename);
	pp = nodename - 1;
	pn = nid;
	/* copy all (and only) the digits to a new buffer */
	while (++pp < nodename + strlen(nodename))
		if (isdigit(*pp))
			*pn++ = *pp;
	*pn = (char)0;
	errno = 0;
	id = strtoull(nid, NULL, 0);
	ierr = errno;
	if (ierr == ERANGE || ierr == EINVAL) {
		perror("strtol");
		exit(ierr);
	}
#else				/* NODEID_<TYPE> not defined. Punt. One MPI rank per node. */
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &tmp);
	assert(ierr == 0);
	id = tmp;
#endif				/* NODEID_<TYPE> */
	return id;
}

void comm_initialize_MPI(int *argc, char **argv[]) {
	int i;
	unsigned long long mynodeid;
	int ierr = 0;
	ierr += MPI_Init(argc, argv);
	ierr += MPI_Comm_rank(MPI_COMM_WORLD, &ThisRankID);
	ierr += MPI_Comm_size(MPI_COMM_WORLD, &NumRanks);
	assert(ierr == 0);
	RootRankID = 0;
	NumStages = comm_ceil2(NumRanks);
	NodeID = (long *)malloc(NumRanks * sizeof(long));
	assert(NodeID != 0);
	NodeID[ThisRankID] = mynodeid = comm_getnodeid();
	ierr += MPI_Allgather(&mynodeid, sizeof(unsigned long long), MPI_BYTE,
			      NodeID, sizeof(unsigned long long), MPI_BYTE, MPI_COMM_WORLD);
	assert(ierr == 0);
	return;
}

void comm_finalize_MPI() {
	int ierr = 0;
	ierr = MPI_Finalize();
	return;
}

void *comm_newbuffer_MPI(int nbytes) {
	commbuff_p p;
	p = (commbuff_p) malloc(sizeof(commbuff_t));
	assert(p != NULL);
	p->l = nbytes;
	p->b = (void *)malloc((size_t) nbytes);
	assert(p->b != NULL);
	return (void *)p;
}

void comm_initialize(int *argc, char **argv[]) {
#ifdef SHMEM
	comm_initialize_SHMEM(argc, argv);
#else
	comm_initialize_MPI(argc, argv);
#endif
	return;
}

void comm_finalize() {
#ifdef SHMEM
	comm_finalize_SHMEM();
#else
	comm_finalize_MPI();
#endif
	return;
}

void *comm_newbuffer(int nbytes) {
#ifdef SHMEM
	return comm_newbuffer_SHMEM(nbytes);
#else
	return comm_newbuffer_MPI(nbytes);
#endif
}

int comm_ceil2(int n) {
	int c = 1;
	while (c < n)
		c *= 2;
	return c;
}

void comm_showmapping() {
	int i;
	FILE *Fmapping;
	char fname[512];
	ROOTONLY {
		if (RankMapping ==1) {
			snprintf(fname, 512, "%s/%s.RankToNodeMapping.%dof%d", CaseName, "global", ThisRankID, NumRanks);
			Fmapping = fopen(fname, "w");
			fprintf(Fmapping, "#Task Rank, NodeIDs, (from RootTask)\n");
			for (i = 0; i < NumRanks; i++)
				fprintf(Fmapping,"%d, %d\n", i, NodeID[i]);
			fclose(Fmapping);
		}
	}
	return;
}
