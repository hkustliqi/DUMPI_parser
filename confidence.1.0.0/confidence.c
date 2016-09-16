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
#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "confidence.h"

int main(int argc, char *argv[], char **envp) {
	void *l, *g;

	comm_initialize(&argc, &argv);
	ROOTONLY printf("%s\n", COPYRIGHT);
	getoptions(argc, argv, envp);
	ROOTONLY printf("Confidence: calibrating timer...\n");
	ORB_calibrate();	/* assumes consistent clock behavior across cluster */
	l = measurement_create("local", BufLen, NumBins, BinSize);
	g = measurement_create("global", BufLen, NumBins, BinSize);

	ROOTONLY printf("Confidence: testing...\n");
	measurement_collect(l);
	ROOTONLY printf("Confidence: local analysis...\n");
	measurement_analyze(l, -1.0);
	ROOTONLY printf("Confidence: remote analysis\n");
	measurement_aggregate(g, l);
	measurement_analyze(g, -1.0);
	measurement_cryweird(g, l);
	ROOTONLY printf("Confidence: saving results\n");
	ROOTONLY measurement_serialize(g, RootRankID);
	/* measurement_serialize(l, ThisRankID); */

	ROOTONLY {

		/* for the folks who want a simple score rather than the full data in the histogram and summary files */
#define NODIVIDEBYZERO(_N_) ( (_N_ == 0) ? (1) : (_N_))
		printf("Experimental Figure-of-Merit for off-node communication:  %6.3f%%\n",
		       (((measurement_p) g)->offNodePairwiseMinimum.min0 -
			((measurement_p) g)->offNodeOnesidedMinimum.min0)
		       / NODIVIDEBYZERO(sqrt(sqrt(((measurement_p) g)->offNodePairwise.m40))) * 100.0);

		/* 
		   printf("OneSided Score (Topology Neutral): %6.3f%%\n",
		   ((measurement_p)g)->offNodeOnesidedMinimum.min0 /
		   sqrt(sqrt(((measurement_p)g)->offNodeOnesided.m40))*100.0); printf("Pairwise Score (Topology
		   Neutral): %6.3f%%\n", ((measurement_p)g)->offNodePairwiseMinimum.min0 /
		   sqrt(sqrt(((measurement_p)g)->offNodePairwise.m40))*100.0); printf("OneSided Score (Topology
		   Adjusted): %6.3f%%\n", sqrt(sqrt(((measurement_p)g)->offNodeOnesidedMinimum.m40 /
		   ((measurement_p)g)->offNodeOnesided.m40))*100.0); printf("Pairwise Score (Topology Adjusted):
		   %6.3f%%\n", sqrt(sqrt(((measurement_p)g)->offNodePairwiseMinimum.m40 /
		   ((measurement_p)g)->offNodePairwise.m40))*100.0); */
	}

	g = measurement_destroy(g);
	l = measurement_destroy(l);
	comm_finalize();
	return 0;
}
