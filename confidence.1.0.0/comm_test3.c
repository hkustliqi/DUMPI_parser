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
/*
 * Tests a single (simultaneous, bidirectional) exchange, but extracts
 * both one-sided and pairwise variability.
 *
 * Pros: 
 * - Assesses both one-sided and pairwise variability with minimal averaging
 * - Provides a Least Upper Bound of the network's minimum latency
 * - Quantifies network topology effects 
 * - Provides a baseline minimum for comparison
 *
 * Cons:
 * - Requires additional storage
 */

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

void comm_test3(measurement_p p) {
#ifdef SHMEM
	comm_SHMEM_test3(p);
#else
	comm_MPI_test3(p);
#endif
	return;
}

void comm_SHMEM_test3(measurement_p p) {
	return;
}

void comm_MPI_test3(measurement_p p) {
	commbuff_t *sbuf, *rbuf;
	ul64 *ptimd, *posd, *ppwd, *posmd, *ppwmd;
	double *cos, *cpw, *t;
	double cosmin, cpwmin, binwidth, maxhisttime;
	int i, icycle, istage, ibucket, ierr, thatRankID, nbins;
	ORB_t t1, t2, t3;
	MPI_Status mpistatus;
	sbuf = comm_newbuffer(p->buflen);	/* exchange buffers */
	rbuf = comm_newbuffer(p->buflen);
	cos = (double *)malloc(NumMessages * sizeof(double));	/* array for onesided kernel timings */
	assert(cos != NULL);
	cpw = (double *)malloc(NumMessages * sizeof(double));	/* array for pairwise kernel timings */
	assert(cpw != NULL);
	t = (double *)malloc(NumMessages * sizeof(double));	/* array for timer overehad timings */
	assert(t != NULL);
	ierr = MPI_Barrier(MPI_COMM_WORLD);	/* pre-synchronize all tasks */
	/*****************************************************************************
	 * A full set of samples for this task consists of message exchanges with each
	 * possible partner. The innermost loop below exchanges some number of messages
	 * between a particular pairing of partners. The middle loop steps through the
	 * possible partners. While the outmost allows us to aggregate multiple sets
	 * of samples to increase the total number of samples.
	 *****************************************************************************/
	for (icycle = 0; icycle < NumCycles; icycle++) {
		for (istage = 0; istage < NumStages; istage++) {	/* step through the stage schedule */
			thatRankID = ThisRankID ^ istage;	/* who's my buddy for this stage? */
			if ((thatRankID < NumRanks) && (thatRankID != ThisRankID)) {	/* valid pairing */
				/* valid pair, proceed with test */
				ierr = 0;
				/***************************************/
				/* warm-up / pre-synchronize this pair */
				/***************************************/
				for (i = 0; i < NumWarmUp; i++) {
					ORB_read(t1);
					ORB_read(t2);
					ierr += MPI_Sendrecv(sbuf->b, p->buflen, MPI_BYTE, thatRankID, 0,
							     rbuf->b, p->buflen, MPI_BYTE, thatRankID, 0,
							     MPI_COMM_WORLD, &mpistatus);
					ORB_read(t3);
				}
				assert(ierr == 0);
				/************************************************************/
				/* BEGIN PERFORMANCE KERNEL -- gather samples for this pair */
				/************************************************************/
				for (i = 0; i < NumMessages; i++) {
					ORB_read(t1);	/* for timer overhead estimate */
					ORB_read(t2);
					/***************************************/
					/* begin timed communication primitive */
					/***************************************/
					ierr += MPI_Sendrecv(sbuf->b, p->buflen, MPI_BYTE, thatRankID, 0,
							     rbuf->b, p->buflen, MPI_BYTE, thatRankID, 0,
							     MPI_COMM_WORLD, &mpistatus);
					/*************************************/
					/* end timed communication primitive */
					/*************************************/
					ORB_read(t3);
					/* save the timings */
					t[i] = ORB_seconds(t2, t1);
					cos[i] = ORB_seconds(t3, t2);
				}
				/************************************************************/
				/* END PERFORMANCE KERNEL -- samples gathered for this pair */
				/************************************************************/
				assert(ierr == 0);
				/* exchange array of local timings with partner */
				ierr += MPI_Sendrecv(cos, NumMessages, MPI_DOUBLE, thatRankID, 0,
						     cpw, NumMessages, MPI_DOUBLE, thatRankID, 0,
						     MPI_COMM_WORLD, &mpistatus);
				assert(ierr == 0);
				for (i = 0; i < NumMessages; i++) {	/* pairwise as average, comparable to one-sided 
									 */
					cpw[i] = (cpw[i] + cos[i]) / 2.0;
				}
				/* bin the t, cos, and cpw results for this cycle of this pair in p */
				measurement_bin(p, t, cos, cpw, (NodeID[ThisRankID] == NodeID[thatRankID]));
			}	/* if valid pairing */
		}		/* for istage */
	}			/* for icycle */
	free(t);
	free(cpw);
	free(cos);
	free(rbuf);
	free(sbuf);
	return;
}
