#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include <fftw3.h>
#include "fitsio.h"
#include "readPfits.h"
//#include "ptime.h"
#include <gsl/gsl_multimin.h>

#define NP 2048
#define K 4149.37759 // 1.0/2.41
//#define K 4148.808

typedef struct subintegration{
	char fname[128];     // name of data file

	int indexSub;        // subintegration index
	int indexChn;        // channel index
	double *p_multi;

	double *rms;    // rms for each profile
} subintegration;

int read_prof (subintegration *sub, pheader *header);

int cal_rms (double phase, double b, double *rms, params param);
int find_peak (int n0, int n, double *s, int *position);
double find_peak_value (int n, double *s);
int def_off_pulse (int nphase, double *in, double frac_off);
int off_pulse (int nphase, int index, double *in, double *out, double frac_off);
int remove_baseline (double *in, int index, double frac_off, int n, double *out);

void initialiseSub(subintegration *sub, pheader *header);
void demallocSub(subintegration *sub, pheader *phead);

