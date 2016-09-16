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
 * Tests a single (simultaneous, bidirectional) exchange.
 *
 * Pros: 
 * Emphasizes the network variability with minimal impact from averaging.
 *
 * Cons:
 * Can significantly underestimate the network's actual minimum latency.
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

void comm_test1(measurement_p p) {
#ifdef SHMEM
	comm_SHMEM_test1(p);
#else
	comm_MPI_test1(p);
#endif
	return;
}

void comm_MPI_test1(measurement_p p) {
	commbuff_t *sbuf, *rbuf;
	double delta, binwidth, maxhisttime;
	int i, istage, ibucket, ierr, thatRankID, nbins;
	ORB_t t1, t2, t3;
	MPI_Status mpistatus;
	nbins = p->nbins;
	binwidth = p->binwidth;
	sbuf = comm_newbuffer(p->buflen);	/* allocate the exchange buffers */
	rbuf = comm_newbuffer(p->buflen);
	ierr = MPI_Barrier(MPI_COMM_WORLD);	/* pre-synchronize all tasks */
	for (istage = 0; istage < NumStages; istage++) {	/* step through the stage schedule */
		thatRankID = ThisRankID ^ istage;
		if ((thatRankID < NumRanks) && (thatRankID != ThisRankID)) {
			ierr = 0;
			for (i = 0; i < NumWarmUp; i++) {	/* warm-up / pre-synchronize this pair */
				ORB_read(t1);
				ORB_read(t2);
				ierr += MPI_Sendrecv(sbuf->b, p->buflen, MPI_BYTE, thatRankID, 0,
						     rbuf->b, p->buflen, MPI_BYTE, thatRankID, 0,
						     MPI_COMM_WORLD, &mpistatus);
				ORB_read(t3);
				assert(ierr == 0);
			}
			for (i = 0; i < NumMessages; i++) {	/* gather histogram data for this pair */
				ORB_read(t1);	/* for timer overhead estimate */
				ORB_read(t2);	/***************** begin timed communication primitive */
				ierr += MPI_Sendrecv(sbuf->b, p->buflen, MPI_BYTE, thatRankID, 0,
						     rbuf->b, p->buflen, MPI_BYTE, thatRankID, 0,
						     MPI_COMM_WORLD, &mpistatus);
				ORB_read(t3);	/******************* end timed communication primitive */
				assert(ierr == 0);
				/* bin the communication cost */
				delta = ORB_seconds(t3, t2);
				/* Note for future: */
				/* could exchange this delta with thatRankID to estimate L1+L2 */
				/* for a least upper bound approximation to L0 vs Lmin==Lsr */
				/* then bin the LUB estimates for L0 for each pair to obtain a */
				/* stepped topology representation */
				ibucket = time2bin(delta);
				if (ibucket >= 0) {
					if (ibucket > (nbins - 1))
						ibucket = nbins - 1;
					if (NodeID[ThisRankID] == NodeID[thatRankID]) {
						(p->onNodeOnesided.dist)[ibucket]++;	/* bin as local */
					} else {
						(p->offNodeOnesided.dist)[ibucket]++;	/* bin as remote */
					}
				}
				/* bin the timer cost */
				delta = ORB_seconds(t2, t1);
				ibucket = time2bin(delta);
				if (ibucket >= 0) {
					if (ibucket > (nbins - 1))
						ibucket = nbins - 1;
					(p->tim.dist)[ibucket]++;
				}
			}
		}
	}
	free(rbuf);
	free(sbuf);
	return;
}
