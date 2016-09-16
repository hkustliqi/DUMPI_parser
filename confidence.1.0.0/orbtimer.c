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
/**************************************************************************
 * ORBtimer -- The Oak Ridge Benchmarking Timer Library
 *
 * Author: Jeffery A. Kuehn
 * Copyright (C) 2009 UT-Battelle, LLC.
 * All Rights Reserved.
 *
 **************************************************************************/

#include <unistd.h>
#include <sys/time.h>
#include "config.h"

#define  ORBTIMER_LIBRARY
#include "orbtimer.h"

static ORB_tick_t Csum = 0;
static ORB_tick_t Gsum = 0;
static ORB_tick_t Dsum = 0;
static ORB_tick_t Fsum = 0;
static ORB_tick_t nsamples = 0;
static ORB_tick_t ndummy = 0;

void ORB_calibrate() {
	int i, j;
	double seconds;
	struct timeval tv1, tv2;
	ORB_t t0, t1, t2, t3, t4, cycles;
	ORB_tick_t nsam;
	ORB_tick_t cmin, gmin, csum, gsum, c21, c32;

#if defined(ORB_IS_FIXEDFREQUENCY)
	ORB_ref_freq = ORB_IS_FIXEDFREQUENCY;
#else
#if defined(ORB_IS_FLOATINGPOINT)
#   error .........................................................
#   error .... Calibration algorithm does not support floating ....
#   error .... point timers without a known fixed frequency    ....
#   error .........................................................
#endif /* ORB_IS_FLOATINGPOINT  */
#endif /* ORB_IS_FIXEDFREQUENCY */
	cmin = ORB_min_lat_cyc;
	gmin = GTD_min_lat_cyc + ORB_avg_lat_cyc;
	for (j = 0; j < 4; j++) {
		nsam = csum = gsum = 0;	/* keep only last */
		for (i = 0; i < 1000000; i++) {	/* sample */
			ORB_read(t1);
			ORB_read(t2);
			gettimeofday(&tv1, 0);
			ORB_read(t3);
			c21 = ORB_cycles_u(t2, t1);
			c32 = ORB_cycles_u(t3, t2);
			if ((c21 >= 0) && (c32 >= 0)) {	/* GTD timers aren't monotonic */
				if (c21 < cmin)
					cmin = c21;
				if (c32 < gmin)
					gmin = c32;
				csum += c21;
				gsum += c32;
				nsam++;
			}
		}
		ndummy += nsam + csum + gsum;
	}
	Csum += csum;
	Gsum += gsum;
	nsamples += nsam;
	ORB_avg_lat_cyc = (Csum + (nsamples >> 1)) / nsamples;
	ORB_min_lat_cyc = cmin;
	GTD_avg_lat_cyc = (Gsum - Csum + (nsamples >> 1)) / nsamples;
	GTD_min_lat_cyc = gmin - ORB_avg_lat_cyc;
#if !defined(ORB_IS_FIXEDFREQUENCY)	/* discover frequency */
	gettimeofday(&tv1, 0);
	ORB_read(t1);
	sleep(5);		/* region = sleep()+2*ORB()+(0.5+0.5)*gtd() */
	ORB_read(t2);
	gettimeofday(&tv2, 0);
	seconds = ((double)(tv2.tv_sec - tv1.tv_sec)) + ((double)(tv2.tv_usec - tv1.tv_usec)) / 1.0e+6;
	cycles = ORB_cycles_u(t2, t1) + ORB_avg_lat_cyc + GTD_avg_lat_cyc;
	ORB_ref_freq = ((double)(cycles)) / seconds;
#endif				/* ORB_IS_FIXEDFREQUENCY */
	ORB_avg_lat_sec = ORB_avg_lat_cyc / ORB_ref_freq;
	ORB_min_lat_sec = ORB_min_lat_cyc / ORB_ref_freq;
	GTD_avg_lat_sec = GTD_avg_lat_cyc / ORB_ref_freq;
	GTD_min_lat_sec = GTD_min_lat_cyc / ORB_ref_freq;
}
