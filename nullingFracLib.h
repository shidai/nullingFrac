#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
//#include <fftw3.h>
#include "fitsio.h"
#include "readPfits.h"
#include "cpgplot.h"
//#include "ptime.h"
//#include <gsl/gsl_multimin.h>

#define NP 2048
#define K 4149.37759 // 1.0/2.41
//#define K 4148.808

typedef struct subintegration{
	char fname[128];     // name of data file

	int indexSub;        // subintegration index
	int indexChn;        // channel index
	double *p_multi;
	double *dat_scl;
	double *dat_offs;
} subintegration;

int read_prof (subintegration *sub, pheader *header);

int def_off_pulse (int nphase, double *in, double frac_off);
int def_on_pulse (int nphase, double *in, double frac_on);
int get_region (int nphase, int index, double *in, double *out, double frac_off);
int cal_pulseEnergy (double *in, double frac_on, int n, double out[2]);

void initialiseSub (subintegration *sub, pheader *header);
void demallocSub (subintegration *sub, pheader *phead);

int histogram (double *data, int n, float *x, float *val, float low, float up, int step);
int makePlot (char *fname, char *dname, int number);
float find_max_value (int n, float *s);
