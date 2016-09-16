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
#include <stdio.h>
#include <assert.h>
#include "config.h"
#include "confidence.h"

/**********************************************
 * Convert a time (in seconds) to a bin number
 **********************************************/
inline int time2bin(double t) {
	int b;
	if (t < MaxHistTime) {
		if (LogarithmicBinning == 0) {	/* LINEAR binning */
			b = (int)(t / BinSize);
		} else {	/* LOGARITHMIC BINNING */
			b = (int)(HistScale * log(t / BinSize));
			b = ((b > 0) ? b : 0);	/* drop the small ones in the zeroeth bin */
		}
	} else {		/* drop the big ones in the last bin */
		b = NumBins - 1;
	}
	return b;
}

/**************************************************
 * Convert a bin to it's corresponding bottom time
 **************************************************/
inline double bin2time(int b) {
	double t;
	if (b == 0)
		return 0.0;	/* first bin is special... goes to zero */
	if (LogarithmicBinning == 0) {	/* LINEAR binning */
		t = BinSize * ((double)b);
	} else {		/* LOGARITHMIC binning */
		t = BinSize * exp((((double)b) / HistScale));
	}
	return t;
}

/***************************************************
 * Convert a bin to it's corresponding midpoint time
 ***************************************************/
inline double bin2midtime(int b) {
	double t;
	if (LogarithmicBinning == 0) {	/* LINEAR binning */
		t = BinSize * ((double)b + 0.5);
	} else {		/* LOGARITHMIC binning */
		t = BinSize * exp((((double)b + 0.5) / HistScale));
	}
	return t;
}

/**********************************************
 * Compute the moments for a histogram
 **********************************************/
void measurement_moments(histogram_p h, double center, ul64 nbins, double binwidth, double *m1, double *m2, double *m3,
			 double *m4) {
	/* calculate four moments of a distribution about a particular bin, given nbins, and binwidth */
	int i;
	double x;
	*m1 = 0.0;
	*m2 = 0.0;
	*m3 = 0.0;
	*m4 = 0.0;
	if (h->nsamples != 0) {
		for (i = 0; i < nbins; i++) {
			x = bin2midtime(i) - center;
			*m1 += (h->dist)[i] * x;
			*m2 += (h->dist)[i] * x * x;
			*m3 += (h->dist)[i] * x * x * x;
			*m4 += (h->dist)[i] * x * x * x * x;
		}
		*m1 /= h->nsamples;
		*m2 /= h->nsamples;
		*m3 /= h->nsamples;
		*m4 /= h->nsamples;
	}
	return;
}

/**********************************************
 * Count the samples in a histogram
 **********************************************/
ul64 measurement_samplecount(ul64 * dist, int nbins) {
	int i;
	ul64 nsamples;
	nsamples = 0;
	for (i = 0; i < nbins; i++) {
		nsamples += dist[i];
	}
	return nsamples;
}

/**********************************************
 * Compute the statistics on a histogram
 **********************************************/
void measurement_histogram(histogram_p h, int nbins, double binwidth, double scale) {
	/* compute stats for an individual histogram */
	double s;
	int i, j;
	ul64 tmp;
	tmp = 0;
	h->nsamples = measurement_samplecount(h->dist, nbins);	/* samples */
	i = -1;
	while ((h->dist)[++i] == 0) ;	/* minimum */
	h->min0 = bin2midtime(i);
	j = 0;
	for (i = 0; i < nbins; i++)
		if ((h->dist)[i] > (h->dist)[j])
			j = i;		/* mode */
	h->mod0 = bin2midtime(j);
	i = -1;
	while ((tmp += (h->dist)[++i]) < (h->nsamples) / 2) ;	/* median */
	h->med0 = bin2midtime(i);
	i = nbins;
	while ((h->dist)[--i] == 0) ;	/* maximum */
	h->max0 = bin2midtime(i);
	/* compute moments */
	measurement_moments(h, 0.0, nbins, binwidth, &(h->m10), &(h->m20), &(h->m30), &(h->m40));
	measurement_moments(h, (h->min0), nbins, binwidth, &(h->m1m), &(h->m2m), &(h->m3m), &(h->m4m));
	if (scale <= 0.0) {
		s = (h->min0);
	} else {
		s = scale;
	}
	h->mods = (h->mod0) / s;
	h->meds = (h->med0) / s;
	h->maxs = (h->max0) / s;
	h->m1s = (h->m10) / s;	/* 1st moment: scaled */
	h->m2s = (h->m20) / s / s;	/* 2nd moment: scaled */
	h->m3s = (h->m30) / s / s / s;	/* 3rd moment: scaled */
	h->m4s = (h->m40) / s / s / s / s;	/* 4th moment: scaled */

	/* old code... this form of the scaled metrics worked well, but confused novices */
	/* h->m1s = (h->m1m) / s; /* 1st moment: scaled */
	/* h->m2s = (h->m2m) / s / s; /* 2nd moment: scaled */
	/* h->m3s = (h->m3m) / s / s / s; /* 3rd moment: scaled */
	/* h->m4s = (h->m4m) / s / s / s / s; /* 4th moment: scaled */
}

/**********************************************
 * Call the analysis routines on each histogram
 **********************************************/
void measurement_analyze(void *m, double scale) {
	/* analyze raw measurement data with summary statistics */
	measurement_p p;
	p = (measurement_p) m;
	measurement_histogram(&(p->onNodeOnesided), p->nbins, p->binwidth, scale);
	measurement_histogram(&(p->onNodePairwise), p->nbins, p->binwidth, scale);
	measurement_histogram(&(p->onNodeOnesidedMinimum), p->nbins, p->binwidth, scale);
	measurement_histogram(&(p->onNodePairwiseMinimum), p->nbins, p->binwidth, scale);
	measurement_histogram(&(p->offNodeOnesided), p->nbins, p->binwidth, scale);
	measurement_histogram(&(p->offNodePairwise), p->nbins, p->binwidth, scale);
	measurement_histogram(&(p->offNodeOnesidedMinimum), p->nbins, p->binwidth, scale);
	measurement_histogram(&(p->offNodePairwiseMinimum), p->nbins, p->binwidth, scale);
	measurement_histogram(&(p->tim), p->nbins, p->binwidth, scale);
}

/**********************************************
 * constructor routine for MEASUREMENT
 **********************************************/
void *measurement_create(char *label, int buflen, int numbins, double binsize) {
	/* create a container to hold a collected dataset */
	measurement_p p;
	p = (measurement_p) malloc(sizeof(measurement_t));
	assert(p != NULL);

	p->tim.dist = (ul64 *) calloc((size_t) numbins, (size_t) sizeof(ul64));
	assert(p->tim.dist != NULL);
	p->tim.label = "timer";

	p->onNodeOnesided.dist = (ul64 *) calloc((size_t) numbins, (size_t) sizeof(ul64));
	assert(p->onNodeOnesided.dist != NULL);
	p->onNodeOnesided.label = "onNodeOnesided";

	p->onNodeOnesidedMinimum.dist = (ul64 *) calloc((size_t) numbins, (size_t) sizeof(ul64));
	assert(p->onNodeOnesidedMinimum.dist != NULL);
	p->onNodeOnesidedMinimum.label = "onNodeOnesidedMinimum";

	p->onNodePairwise.dist = (ul64 *) calloc((size_t) numbins, (size_t) sizeof(ul64));
	assert(p->onNodePairwise.dist != NULL);
	p->onNodePairwise.label = "onNodePairwise";

	p->onNodePairwiseMinimum.dist = (ul64 *) calloc((size_t) numbins, (size_t) sizeof(ul64));
	assert(p->onNodePairwiseMinimum.dist != NULL);
	p->onNodePairwiseMinimum.label = "onNodePairwiseMinimum";

	p->offNodeOnesided.dist = (ul64 *) calloc((size_t) numbins, (size_t) sizeof(ul64));
	assert(p->offNodeOnesided.dist != NULL);
	p->offNodeOnesided.label = "offNodeOnesided";

	p->offNodeOnesidedMinimum.dist = (ul64 *) calloc((size_t) numbins, (size_t) sizeof(ul64));
	assert(p->offNodeOnesidedMinimum.dist != NULL);
	p->offNodeOnesidedMinimum.label = "offNodeOnesidedMinimum";

	p->offNodePairwise.dist = (ul64 *) calloc((size_t) numbins, (size_t) sizeof(ul64));
	assert(p->offNodePairwise.dist != NULL);
	p->offNodePairwise.label = "offNodePairwise";

	p->offNodePairwiseMinimum.dist = (ul64 *) calloc((size_t) numbins, (size_t) sizeof(ul64));
	assert(p->offNodePairwiseMinimum.dist != NULL);
	p->offNodePairwiseMinimum.label = "offNodePairwiseMinimum";

	p->nbins = numbins;
	p->binwidth = binsize;
	p->buflen = buflen;
	p->label = (char *)malloc(strlen(label) + 2);
	assert(p->label != NULL);
	strcpy(p->label, label);
	return (void *)p;
}

/**********************************************
 * Destructor routine for MEASUREMENT
 **********************************************/
void *measurement_destroy(void *m) {
	/* destroy a container holding a collected dataset */
	measurement_p p;
	p = (measurement_p) m;
	free(p->label);
	free(p->offNodePairwiseMinimum.dist);
	free(p->offNodePairwise.dist);
	free(p->offNodeOnesidedMinimum.dist);
	free(p->offNodeOnesided.dist);
	free(p->onNodePairwiseMinimum.dist);
	free(p->onNodePairwise.dist);
	free(p->onNodeOnesidedMinimum.dist);
	free(p->onNodeOnesided.dist);
	free(p->tim.dist);
	free(p);
	return (void *)NULL;
}

/**********************************************
 * Select among comm_test?() based on cmdline
 **********************************************/
void measurement_collect(void *m) {
	/* capture timings to compose raw data into a measurement */
	/* eventually select among several tests */
	switch (TestType) {
	case 1:
		comm_test1((measurement_p) m);
		break;
	case 2:
		comm_test2((measurement_p) m);
		break;
	case 3:
	default:
		comm_test3((measurement_p) m);
		break;
	}
}

/**********************************************
 * Bin a block of measurements between a pair
 **********************************************/
void measurement_bin(measurement_t * p, double *t, double *cos, double *cpw, int LOCAL) {
	double cosmin, cpwmin, binwidth, maxhisttime;
	ul64 *ptimd, *posd, *ppwd, *posmd, *ppwmd;
	int i, nbins;
	if (LOCAL) {		/* bin these values as local communication */
		posd = p->onNodeOnesided.dist;
		ppwd = p->onNodePairwise.dist;
		posmd = p->onNodeOnesidedMinimum.dist;
		ppwmd = p->onNodePairwiseMinimum.dist;
	} else {		/* bin these values as remote communication */
		posd = p->offNodeOnesided.dist;
		ppwd = p->offNodePairwise.dist;
		posmd = p->offNodeOnesidedMinimum.dist;
		ppwmd = p->offNodePairwiseMinimum.dist;
	}
	ptimd = p->tim.dist;
	cosmin = cpwmin = 1.0e+16;
	for (i = 0; i < NumMessages; i++) {
		/* bin the individual results */
		if (t[i] >= 0.0)
			ptimd[time2bin(t[i])]++;
		if (cos[i] >= 0.0)
			posd[time2bin(cos[i])]++;
		if (cpw[i] >= 0.0)
			ppwd[time2bin(cpw[i])]++;
		/* save the minimums for now */
		if ((cos[i] > 0.0) && (cos[i] < cosmin))
			cosmin = cos[i];
		if ((cpw[i] > 0.0) && (cpw[i] < cpwmin))
			cpwmin = cpw[i];
	}
	/* now bin the minimums for this communications pair */
	if (cosmin > 0.0)
		posmd[time2bin(cosmin)]++;
	if (cpwmin > 0.0)
		ppwmd[time2bin(cpwmin)]++;
}

/****************************************************************************
 * FUTURE HOOK: compare local measurements to global "norms" to find outliers
 ****************************************************************************/
void measurement_cryweird(void *normal, void *me) {
}

/**********************************************
 * call through to comm_ layer aggregation
 **********************************************/
void measurement_aggregate(void *g, void *l) {
	/* aggregate raw data from measurements across the system into a single raw measurement */
	comm_aggregate((measurement_p) g, (measurement_p) l);
}

#define FNAMESIZE 512
/* FIXME */
/**********************************************
 * save the data to disk
 **********************************************/
void measurement_serialize(void *m, int writingRankID) {
	/* if I am the target, then I will output the raw data and statistics for visualization */
	int i;
	FILE *Fhist, *Fpdf, *Fcdf;
	char fname[FNAMESIZE];
	char dname[FNAMESIZE];
	double r, s, t, u, v, w, x, y, z;
	double binwidth, binbot, binmid, bintop, maxhisttime;
	measurement_p p;
	p = (measurement_p) m;
	/* have to switch from binwidth based descriptions here... */
	if (writingRankID == ThisRankID) {
		mkdir(CaseName, 0755);
		/* raw histogram data */
		snprintf(fname, FNAMESIZE, "%s/%s.HIST.%d", CaseName, p->label, ThisRankID);
		Fhist = fopen(fname, "w");
		fprintf(Fhist, "%s",COPYRIGHT);
		fprintf(Fhist, "# Casename:          %s\n", CaseName);
		fprintf(Fhist, "# TestType:          %d\n", TestType);
		fprintf(Fhist, "# NumRanks:          %d\n", NumRanks);
		fprintf(Fhist, "# MLabel:            %s\n", p->label);
		fprintf(Fhist, "# Message Size:      %d\n", BufLen);
		fprintf(Fhist, "# Message Pattern:   %d cycle(s) through an all-pairs schedule\n", NumCycles);
		fprintf(Fhist, "#                    of %d warmups and %d messages per pair\n", NumWarmUp, NumMessages);
		if (LogarithmicBinning == 1) {
			fprintf(Fhist, "# Binning:           Logarithmic, ending at %g seconds\n", MaxHistTime);
		} else {
			fprintf(Fhist, "# Binning:           Linear, ending at %g seconds\n", MaxHistTime);
		}
		fprintf(Fhist, "#%6s %18s %17s %15s %15s %15s %15s %15s %15s %15s %15s\n",
			"bin", " (us) to  (us)", "timer",
			"onNd-OS", "onNd-PW", "onNd-OS-min", "onNd-PW-min",
			"offNd-OS", "offNd-PW", "offNd-OS-min", "offNd-PW-min");
		for (i = 0; i < p->nbins; i++) {
			binbot = bin2time(i);
			bintop = bin2time((i + 1));
			binwidth = bintop - binbot;
			fprintf(Fhist,
				"%6d %11.4g %11.4g %15llu %15llu %15llu %15llu %15llu %15llu %15llu %15llu %15llu\n", i,
				binbot * 1.0e+6, bintop * 1.0e+6,
				(p->tim.dist)[i],
				(p->onNodeOnesided.dist)[i],
				(p->onNodePairwise.dist)[i],
				(p->onNodeOnesidedMinimum.dist)[i],
				(p->onNodePairwiseMinimum.dist)[i],
				(p->offNodeOnesided.dist)[i],
				(p->offNodePairwise.dist)[i],
				(p->offNodeOnesidedMinimum.dist)[i],
				(p->offNodePairwiseMinimum.dist)[i]);
		}
		fclose(Fhist);
		/* PDF data */
		snprintf(fname, FNAMESIZE, "%s/%s.PDF.%d", CaseName, p->label, ThisRankID);
		Fpdf = fopen(fname, "w");
		fprintf(Fpdf, "%s",COPYRIGHT);
		fprintf(Fpdf, "# Casename:          %s\n", CaseName);
		fprintf(Fpdf, "# TestType:          %d\n", TestType);
		fprintf(Fpdf, "# NumRanks:          %d\n", NumRanks);
		fprintf(Fpdf, "# MLabel:            %s\n", p->label);
		fprintf(Fpdf, "# Message Size:      %d\n", BufLen);
		fprintf(Fpdf, "# Message Pattern:   %d cycle(s) through an all-pairs schedule\n", NumCycles);
		fprintf(Fpdf, "#                    of %d warmups and %d messages per pair\n", NumWarmUp, NumMessages);
		if (LogarithmicBinning == 1) {
			fprintf(Fpdf, "# Binning:           Logarithmic, ending at %g seconds\n", MaxHistTime);
		} else {
			fprintf(Fpdf, "# Binning:           Linear, ending at %g seconds\n", MaxHistTime);
		}
		fprintf(Fpdf, "#%6s %18s %17s %15s %15s %15s %15s %15s %15s %15s %15s\n",
			"bin", " (us) to  (us)", "timer",
			"onNd-OS", "onNd-PW", "onNd-OS-min", "onNd-PW-min",
			"offNd-OS", "offNd-PW", "offNd-OS-min", "offNd-PW-min");
		for (i = 0; i < p->nbins; i++) {
			binbot = bin2time(i);
			bintop = bin2time((i + 1));
			binwidth = bintop - binbot;
			fprintf(Fpdf,
				"%6d %11.4g %11.4g %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n", i,
				binbot * 1.0e+6, bintop * 1.0e+6,
#define NODIVIDEBYZERO(_N_) ( (_N_ == 0) ? (1) : (_N_))
				(double)(p->tim.dist)[i]                    / binwidth / (double)NODIVIDEBYZERO(p->tim.nsamples),
				(double)(p->onNodeOnesided.dist)[i]         / binwidth / (double)NODIVIDEBYZERO(p->onNodeOnesided.nsamples),
				(double)(p->onNodePairwise.dist)[i]         / binwidth / (double)NODIVIDEBYZERO(p->onNodePairwise.nsamples),
				(double)(p->onNodeOnesidedMinimum.dist)[i]  / binwidth / (double)NODIVIDEBYZERO(p->onNodeOnesidedMinimum.nsamples),
				(double)(p->onNodePairwiseMinimum.dist)[i]  / binwidth / (double)NODIVIDEBYZERO(p->onNodePairwiseMinimum.nsamples),
				(double)(p->offNodeOnesided.dist)[i]        / binwidth / (double)NODIVIDEBYZERO(p->offNodeOnesided.nsamples),
				(double)(p->offNodePairwise.dist)[i]        / binwidth / (double)NODIVIDEBYZERO(p->offNodePairwise.nsamples),
				(double)(p->offNodeOnesidedMinimum.dist)[i] / binwidth / (double)NODIVIDEBYZERO(p->offNodeOnesidedMinimum.nsamples),
				(double)(p->offNodePairwiseMinimum.dist)[i] / binwidth / (double)NODIVIDEBYZERO(p->offNodePairwiseMinimum.nsamples));
		}
		fclose(Fpdf);
		/* CDF data */
		snprintf(fname, FNAMESIZE, "%s/%s.CDF.%d", CaseName, p->label, ThisRankID);
		Fcdf = fopen(fname, "w");
		fprintf(Fcdf, "%s",COPYRIGHT);
		fprintf(Fcdf, "# Casename:          %s\n", CaseName);
		fprintf(Fcdf, "# TestType:          %d\n", TestType);
		fprintf(Fcdf, "# NumRanks:          %d\n", NumRanks);
		fprintf(Fcdf, "# MLabel:            %s\n", p->label);
		fprintf(Fcdf, "# Message Size:      %d\n", BufLen);
		fprintf(Fcdf, "# Message Pattern:   %d cycle(s) through an all-pairs schedule\n", NumCycles);
		fprintf(Fcdf, "#                    of %d warmups and %d messages per pair\n", NumWarmUp, NumMessages);
		if (LogarithmicBinning == 1) {
			fprintf(Fcdf, "# Binning:           Logarithmic, ending at %g seconds\n", MaxHistTime);
		} else {
			fprintf(Fcdf, "# Binning:           Linear, ending at %g seconds\n", MaxHistTime);
		}
		fprintf(Fcdf, "%#6s %18s %17s %15s %15s %15s %15s %15s %15s %15s %15s\n",
			"bin", " (us) to  (us)", "timer",
			"onNd-OS", "onNd-PW", "onNd-OS-min", "onNd-PW-min",
			"offNd-OS", "offNd-PW", "offNd-OS-min", "offNd-PW-min");
		r = s = t = u = v = w = x = y = z = 0.0;
		for (i = 0; i < p->nbins; i++) {
			binbot = bin2time(i);
			bintop = bin2time((i + 1));
			binwidth = bintop - binbot;
			r += (double)(p->tim.dist)[i]                    / (double)NODIVIDEBYZERO(p->tim.nsamples);
			s += (double)(p->onNodeOnesided.dist)[i]         / (double)NODIVIDEBYZERO(p->onNodeOnesided.nsamples);
			t += (double)(p->onNodePairwise.dist)[i]         / (double)NODIVIDEBYZERO(p->onNodePairwise.nsamples);
			u += (double)(p->onNodeOnesidedMinimum.dist)[i]  / (double)NODIVIDEBYZERO(p->onNodeOnesidedMinimum.nsamples);
			v += (double)(p->onNodePairwiseMinimum.dist)[i]  / (double)NODIVIDEBYZERO(p->onNodePairwiseMinimum.nsamples);
			w += (double)(p->offNodeOnesided.dist)[i]        / (double)NODIVIDEBYZERO(p->offNodeOnesided.nsamples);
			x += (double)(p->offNodePairwise.dist)[i]        / (double)NODIVIDEBYZERO(p->offNodePairwise.nsamples);
			y += (double)(p->offNodeOnesidedMinimum.dist)[i] / (double)NODIVIDEBYZERO(p->offNodeOnesidedMinimum.nsamples);
			z += (double)(p->offNodePairwiseMinimum.dist)[i] / (double)NODIVIDEBYZERO(p->offNodePairwiseMinimum.nsamples);
			fprintf(Fcdf,
				"%6d %11.4g %11.4g %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n", i,
				binbot * 1.0e+6, bintop * 1.0e+6, r, s, t, u, v, w, x, y, z);
		}
		fclose(Fcdf);
		/* Statistical Summaries */
		measurement_fmthist(&(p->tim), p->label);
		measurement_fmthist(&(p->onNodeOnesided), p->label);
		measurement_fmthist(&(p->onNodeOnesidedMinimum), p->label);
		measurement_fmthist(&(p->onNodePairwise), p->label);
		measurement_fmthist(&(p->onNodePairwiseMinimum), p->label);
		measurement_fmthist(&(p->offNodeOnesided), p->label);
		measurement_fmthist(&(p->offNodeOnesidedMinimum), p->label);
		measurement_fmthist(&(p->offNodePairwise), p->label);
		measurement_fmthist(&(p->offNodePairwiseMinimum), p->label);
		comm_showmapping();
	}
}

/**********************************************
 * save the histogram summary stats to disk
 **********************************************/
void measurement_fmthist(histogram_p h, char *label) {
	char fname[FNAMESIZE];
	FILE *Fstat;
	snprintf(fname, FNAMESIZE, "%s/%s.STAT.%s.%d", CaseName, label, h->label, ThisRankID);
	Fstat = fopen(fname, "w");
	fprintf(Fstat, "%s",COPYRIGHT);
	fprintf(Fstat, "# Casename:          %s\n", CaseName);
	fprintf(Fstat, "# TestType:          %d\n", TestType);
	fprintf(Fstat, "# NumRanks:          %d\n", NumRanks);
	fprintf(Fstat, "# MLabel:            %s\n", label);
	fprintf(Fstat, "# HLabel:            %s\n", h->label);
	fprintf(Fstat, "# Message Size:      %d\n", BufLen);
	fprintf(Fstat, "# Message Pattern:   %d cycle(s) through an all-pairs schedule\n", NumCycles);
	fprintf(Fstat, "#                    of %d warmups and %d messages per pair\n", NumWarmUp, NumMessages);
	if (LogarithmicBinning == 1) {
		fprintf(Fstat, "# Binning:           Logarithmic, ending at %g seconds\n", MaxHistTime);
	} else {
		fprintf(Fstat, "# Binning:           Linear, ending at %g seconds\n", MaxHistTime);
	}
	fprintf(Fstat, "# Number of Samples: %15llu in %d bins\n", h->nsamples, NumBins);	/* number of samples in the distribution */
	fprintf(Fstat, "\n");
	fprintf(Fstat, "Minimum:      %15.2g usec     %15.2g * minLatency\n", h->min0 * 1.0e+6, 1.0);			/* minimum: 0, 0-scaled */
	fprintf(Fstat, "Mode:         %15.2g usec     %15.2g * minLatency\n", h->mod0 * 1.0e+6, h->mods);		/* mode: 0, 0-scaled */
	fprintf(Fstat, "Median:       %15.2g usec     %15.2g * minLatency\n", h->med0 * 1.0e+6, h->meds);		/* median: 0, 0-scaled */
	fprintf(Fstat, "Mean:         %15.2g usec     %15.2g * minLatency\n", h->m10 * 1.0e+6, h->m10 / h->min0);	/* 1st moment: 0, 0-scaled */
	fprintf(Fstat, "Maximum:      %15.2g usec     %15.2g * minLatency\n", h->max0 * 1.0e+6, h->maxs);		/* maximum: 0, 0-scaled */
	fprintf(Fstat, "\n");
	fprintf(Fstat, "R1(Mean):     %15.2g usec     %15.2g * minLatency\n", h->m10 * 1.0e+6, h->m1s);			/* 1st moment: 0,min-scaled */
	fprintf(Fstat, "R2(Variance): %15.2g usec     %15.2g * minLatency\n", sqrt(h->m20) * 1.0e+6, sqrt(h->m2s));	/* 2nd moment: 0,min-scaled */
	fprintf(Fstat, "R3(Skewness): %15.2g usec     %15.2g * minLatency\n", cbrt(h->m30) * 1.0e+6, cbrt(h->m3s));	/* 3rd moment: 0,min-scaled */
	fprintf(Fstat, "R4(Kurtosis): %15.2g usec     %15.2g * minLatency\n", sqrt(sqrt(h->m40)) * 1.0e+6, sqrt(sqrt(h->m4s)));	/* 4th moment: 0,min-scaled */
	fprintf(Fstat, "\n");
	fprintf(Fstat, "Mean:         %15.2g seconds     %15.2g * minLatency\n", h->m10, h->m1s);			/* 1st moment: 0,min-scaled */
	fprintf(Fstat, "Variance:     %15.2g seconds**2  %15.2g * minLatency**2\n", h->m20, h->m2s);			/* 2nd moment: 0,min-scaled */
	fprintf(Fstat, "Skewness:     %15.2g seconds**3  %15.2g * minLatency**3\n", h->m30, h->m3s);			/* 3rd moment: 0,min-scaled */
	fprintf(Fstat, "Kurtosis:     %15.2g seconds**4  %15.2g * minLatency**4\n", h->m40, h->m4s);			/* 4th moment: 0,min-scaled */
	fclose(Fstat);
}
