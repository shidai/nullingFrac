// ptimeT library
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <fftw3.h>
//#include <gsl/gsl_multimin.h>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>

#include "nullingFracLib.h"
#include "ptimeLib.h"
//#include "T2toolkit.h"
//#include "tempo2pred.h"
#define ITMAX 100000  // Maximum allowed number of iterations.
#define EPS 1.0e-16 // Machine double floating-point precision.
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

int read_prof (subintegration *sub, pheader *header)
{  
  fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
  int status;
  int colnum;

  status = 0;

	// open psrfits
  if ( fits_open_file(&fptr, sub->fname, READONLY, &status) )          // open the file
  {
      printf( "error while openning file\n" );
  }

	// move to subint
	fits_movnam_hdu(fptr, BINARY_TBL, (char *)"SUBINT",0,&status);

  int frow;
  int felem;
  int nelem;
  int null;
  int anynull;

	/////////////////////////////////////////////////////////////////////////
	// read profile
  if ( fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status) )           // get the row number
  {
      printf( "error while getting the colnum number\n" );
	}

	frow = sub->indexSub;
	felem = 1;
	nelem = header->nbin*header->nchan*header->npol;
	null = 0;
	anynull = 0;

	fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, sub->p_multi, &anynull, &status);           // read the column

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// read weights
 
	if ( fits_get_colnum(fptr, CASEINSEN, "DAT_WTS", &colnum, &status) )           // get the colnum number
	{
		printf( "error while getting the colnum number\n" );
	}

	frow = sub->indexSub;
  felem = 1;
  nelem = header->nchan;
  null = 0;
  anynull = 0;

	fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, sub->wts, &anynull, &status);           // read the column

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// read channel frequency
 
	if ( fits_get_colnum(fptr, CASEINSEN, "DAT_FREQ", &colnum, &status) )           // get the colnum number
	{
		printf( "error while getting the colnum number\n" );
	}

	frow = sub->indexSub;
  felem = 1;
  nelem = header->nchan;
  null = 0;
  anynull = 0;

	fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, sub->freq, &anynull, &status);           // read the column

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// read channel frequency SSB
 
	double temp[nelem];
	if ( fits_get_colnum(fptr, CASEINSEN, "FREQ_SSB", &colnum, &status) )           // get the colnum number
	{
		printf( "error while getting the colnum number\n" );
	}

	frow = sub->indexSub;
  felem = 1;
  nelem = header->nchan;
  null = 0;
  anynull = 0;

	fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, temp, &anynull, &status);           // read the column

	int i;
	for (i = 0; i < nelem; i++)
	{
		sub->freqSSB[i] = temp[i]/1000000.0;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// read psr frequency at SSB
 
	if ( fits_get_colnum(fptr, CASEINSEN, "BATFREQ", &colnum, &status) )           // read pulse freq at SSB
	{
		printf( "error while getting the colnum number\n" );
	}

	frow = sub->indexSub;
  felem = 1;
  nelem = 1;
  null = 0;
  anynull = 0;

	fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, &sub->batFreq, &anynull, &status);           // read the column

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// read offs
  
	if ( fits_get_colnum(fptr, CASEINSEN, "OFFS_SUB", &colnum, &status) )           // get the colnum number
	{
		printf( "error while getting the colnum number\n" );
	}

	frow = sub->indexSub;
  felem = 1;
  nelem = 1;
  null = 0;
  anynull = 0;

	fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, &sub->offs, &anynull, &status);           // read the column

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// close psrfits file
	if ( fits_close_file(fptr, &status) )
	{
		printf( " error while closing the file " );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////

	double weight, frequency;
	frequency = 0.0;
	weight = 0.0;
	
	int z;
	for (z = 0; z < header->nchan; z++)
	{
		frequency += sub->freq[z]*sub->wts[z];
		weight += sub->wts[z];
	}
	frequency = frequency/weight;
	sub->frequency = frequency;
	
	// get the pulse period of this subintegration
	long double mjd0;  // the mjd of each subint
	T2Predictor pred;

	// get the period
	//print_t2pred(sub->fname);   // output t2pred.dat
	T2Predictor_Init(&pred);  // prepare the predictor
	
	if (T2Predictor_ReadFits(&pred,sub->fname))
	{
		printf("Error: unable to read predictor\n");
		exit(1);
	}

	// get the period at mjd0
	mjd0 = (long double)(header->imjd) + ((long double)(header->smjd) + (long double)(header->stt_offs) + (long double)(sub->offs))/86400.0L;
	//printf ("imjd is %ld \n", imjd);
	//printf ("mjd0 is %.15Lf \n", mjd0);
	sub->Cperiod = 1.0/T2Predictor_GetFrequency(&pred,mjd0,sub->frequency);

	for (z = 0; z < header->nchan; z++)
	{
		sub->period[z] = 1.0/T2Predictor_GetFrequency(&pred,mjd0,sub->freq[z]);
	}

	T2Predictor_Destroy(&pred);

	return 0;
}

int cal_rms (double phase, double b, double *rms, params param)
// calculate the rms of each subchannel  
{
	double gk;
	int i,j,n;

	gk=0.0;
	n=0;

	for (i = 0; i < param.nchn; i++)
	{
	    for (j = 0; j < param.num; j++)
	    {
		    //gk+=a_p[i]*a_p[i]+b*b*a_s[i]*a_s[i]-2.0*b*a_s[i]*a_p[i]*cos(p_p[i]-p_s[i]-(i+1)*phase)+a*a*1024.0*1024.0-2.0*a*1024.0*a_p[i]*cos(p_p[i])+2.0*a*b*1024.0*a_s[i]*cos(p_s[i]+(i+1)*phase);
		    gk+=param.a_p[i][j]*param.a_p[i][j]+b*b*param.a_s[i][j]*param.a_s[i][j]-2.0*b*param.a_s[i][j]*param.a_p[i][j]*cos(param.p_p[i][j]-param.p_s[i][j]-(j+1)*phase);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		n++;
		}
	}
	
	(*rms)=sqrt(gk/n);

	return 0;
}

int find_peak (int n0, int n, double *s, int *position)
{
	int i;
	double temp[n];
	double peak;

	for (i = n*n0; i < n*n0+n; i++)
	{
		temp[i-n*n0] = s[i];
	}

	double a, b, c;
	for (i = 0; i < n-1; i++)
	{
		//a = fabs(temp[i]);
		//b = fabs(temp[i+1]);
		a = temp[i];
		b = temp[i+1];
		c = (a >= b ? a : b);

		temp[i+1] = c;
	}
	peak = temp[n-1];

	for (i = n*n0; i < n*n0+n; i++)
	{
		if (fabs(peak-s[i]) < 1.0e-3)
		{
			(*position) = i;
		}
	}

	return 0;
}

double find_peak_value (int n, double *s)
{
	int i;
	double temp[n];

	for (i = 0; i < n; i++)
	{
		temp[i] = s[i];
	}

	double a, b, c;
	for (i = 0; i < n-1; i++)
	{
		a = temp[i];
		b = temp[i+1];
		//a = fabs(temp[i]);
		//b = fabs(temp[i+1]);
		//c = (fabs(a) >= fabs(b) ? a : b);
		c = (a >= b ? a : b);

		temp[i+1] = c;
	}

	return temp[n-1];
}

int def_off_pulse (int nphase, double *in, double frac_off)
// define the off pulse region based on I, return the starting index of off pulse region
// using frac_off to calculate the off pulse region
{
	int n = nphase;
	int num_off = (int)(n*frac_off);
	int i,j;
	double small;
	double temp;
	int index = 0;

	for (i = 0; i < n; i++)
	{
		if (i == 0)
		{
			small = 0.0;
			for(j = 0; j < num_off; j++)
			{
				small += (in[j]+30000.0)*(in[j]+30000.0);  // make all numbers positive
			}
			small = sqrt(small/num_off);
		}
			
		temp = 0.0;
		for(j = 0; j < num_off; j++)
		{
			if ((i+j) > n-1)
			{
				temp += (in[(i+j)-(n-1)]+30000.0)*(in[(i+j)-(n-1)]+30000.0);
			}
			else 
			{
				temp += (in[i+j]+30000.0)*(in[i+j]+30000.0);
			}
		}
		temp = sqrt(temp/num_off);

		small = (temp <= small ? temp : small);
		index = (temp <= small ? i : index);
		//printf ("%d %lf %lf\n", index, small, ave);
	}

	return index;
}

int off_pulse (int nphase, int index, double *in, double *out, double frac_off)
// get the off_pulse region
{
	int n = nphase;
	int num_off = (int)(n*frac_off);
	int i;

	for (i = 0; i < num_off; i++)
	{
		if ((index+i) > n-1)
		{
			out[i] = in[(index+i)-(n-1)];
		}
		else 
		{
			out[i] = in[index+i];
		}
	}

	return 0;
}

int remove_baseline (double *in, int index, double frac_off, int n, double *out)
{
	// define the off_pulse range, frac_off is the fraction of the phase
	// index is the starting point of the off_pulse range
	int num_off = (int)(n*frac_off);
	double off_0[num_off];

	off_pulse (n, index, in, off_0, frac_off);

	int i;
	double baseline = 0.0;
    for (i = 0; i < num_off; i++)
    {
        baseline += off_0[i];
        //average_s += s_off[i];
    }
	baseline = baseline/num_off;

    //printf ("the baseline of std is: %lf \n", baseline);
    //printf ("average is: %lf %lf\n", average, average_s);

	for (i = 0; i < n; i++)
	{
		out[i] = (in[i]-baseline);
		//s_norm[i] = (s[i]-baseline)/(s_peak-baseline);
	}
	
	return 0;
}

void initialiseSub(subintegration *sub, pheader *header)
{
	int nphase;
	int nchn;
	int npol;
	
	nchn = header->nchan; 
	npol = header->npol; 
	nphase = header->nbin; 	

	sub->rms = (double *)malloc(sizeof(double)*nchn);
	sub->p_multi = (double *)malloc(sizeof(double)*nchn*npol*nphase);
}

void demallocSub(subintegration *sub, pheader *phead)
{
	free(sub->rms);
	free(sub->p_multi);
	free(sub);
}
