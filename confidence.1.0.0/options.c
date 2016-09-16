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
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "config.h"
#include "confidence.h"

/*
 * defaults for command line options
 */
void setdefaults() {
	NumMessages = 100000;	/* messages per cycle */
	NumCycles = 10;		/* cycles -- future: time limit the data collection */
	NumWarmUp = 100;	/* keep this < 1% of NumMessages */
	NumMessagesGlobal = (ul64) NumCycles *(ul64) NumMessages *(ul64) (NumRanks-1);
	BufLen = 1;		/* small message */
	NumBins = 1000;		/* with log binning, don't need much more */
	BinSize = 50.0e-9;	/* 50ns works well with x86_64 assm timers */
	LogarithmicBinning = 0;
	RankMapping = 0;
	TestType = 3;		/* new style */
	strcpy(CaseName, "OUTPUT_DIRECTORY");	/* user should replace */
	return;
}

/*
 * parse command line options without seatbelts/sanity/consistency checks.
 */
void getoptions(int argc, char *argv[], char **envp) {
	int ierr, opt;
	extern char *optarg;
	extern int optind, opterr, optopt;
	setdefaults();
	ierr = 0;
	while ((opt = getopt(argc, argv, "B:C:G:M:W:m:n:w:T:N:lr")) != -1) {
		switch (opt) {
		case 'B':
			BufLen = strtol(optarg, NULL, 0);
			if (NumCycles == 0)
				ierr++;
			break;
		case 'C':
			NumCycles = strtol(optarg, NULL, 0);
			if (NumCycles == 0)
				ierr++;
			break;
		case 'G':
			NumMessagesGlobal = strtoull(optarg, NULL, 0);
			if (NumMessagesGlobal == 0)
				ierr++;
			break;
		case 'M':
			NumMessages = strtol(optarg, NULL, 0);
			if (NumMessages == 0)
				ierr++;
			break;
		case 'W':
			NumWarmUp = strtol(optarg, NULL, 0);
			if (NumWarmUp == 0)
				ierr++;
			break;
		case 'l':
			LogarithmicBinning = 1;
			break;
		case 'r':
			RankMapping = 1;
			break;
		case 'n':
			NumBins = strtol(optarg, NULL, 0);
			if (NumBins == 0)
				ierr++;
			break;
		case 'w':
			BinSize = strtod(optarg, NULL);
			if (BinSize <= 0.0)
				ierr++;
			break;
		case 'm':
			MaxHistTime = strtod(optarg, NULL);
			if (MaxHistTime <= 0.0)
				ierr++;
			break;
		case 'T':
			TestType = strtol(optarg, NULL, 0);
			if (TestType < 1 || TestType > 3)
				ierr++;
			break;
		case 'N':
			strncpy(CaseName, optarg, NAMEBUFFSIZE);
			if (strlen(CaseName) == 0)
				ierr++;
			break;
		default:	/* ? */
			ierr++;
		}
	}
	if (ierr != 0) {
		ROOTONLY {
			setdefaults();
			fprintf(stderr,
				"Usage: %s [-N s] [-T 1|2|3] [-B n] [-C n] [[-G n]|[-M n]] [-W n] [-l [-m n]] [-n n] [-r]\n",
				argv[0]);
			fprintf(stderr, "\t -N <casename> \t name directory for output (default: %s)\n", CaseName);
			fprintf(stderr, "\t -r            \t save the rank-to-node mapping\n");
			/* fprintf(stderr, "\t -T [1|2|3]    \t which test to run (default: %d)\n", TestType); */ /* tests 1 & 2 are old and undocumented */
			fprintf(stderr, "\t -B <buflen>   \t buffer length for message tests in bytes (default: %d)\n", BufLen);
			fprintf(stderr, "\t -C <cycles>   \t number of cycles of all-pairs collections (default: %d)\n", NumCycles);
			fprintf(stderr, "\t -G <messages> \t total number of global messages to be exchanged\n");
			fprintf(stderr, "\t -M <messages> \t number of messages to exchange per pair (default: %d)\n", NumMessages);
			fprintf(stderr, "\t -W <warmup>   \t number of warm-up messages before timing (default: %d)\n", NumWarmUp);
			fprintf(stderr, "\t -l            \t switch from (default) linear binning to logarithmic binning\n");
			fprintf(stderr, "\t -w <binwidth> \t width of FIRST histogram bin in seconds (default: %g)\n", BinSize);
			fprintf(stderr, "\t -m <time>     \t reset maximum message time to bin (log binning only)\n");
			fprintf(stderr, "\t -n <bins>     \t number of bins in histograms (default: %d)\n", NumBins);
		}
		exit(1);
	}
	if (LogarithmicBinning == 1) {
		MaxHistTime = 1.0;
		HistScale = ((double)NumBins) / log(MaxHistTime / BinSize);
	} else {		/* LINEAR binning */
		HistScale = NumBins * BinSize;
		MaxHistTime = NumBins * BinSize;
	}
}
