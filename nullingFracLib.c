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

	// close psrfits file
	if ( fits_close_file(fptr, &status) )
	{
		printf( " error while closing the file " );
	}

	return 0;
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

int def_on_pulse (int nphase, double *in, double frac_on)
// define the on pulse region based on I, return the starting index of on pulse region
{
	int n = nphase;
	int num_off = (int)(n*frac_on);
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

		small = (temp >= small ? temp : small);
		index = (temp >= small ? i : index);
		//printf ("%d %lf %lf\n", index, small, ave);
	}

	return index;
}

int get_region (int nphase, int index, double *in, double *out, double frac_off)
// get the off_pulse or on_pulse region
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

int cal_pulseEnergy (double *in, double frac_on, int n, double out[2])
{
	// define the on_pulse range, frac_on is the fraction of the phase
	// define the off_pulse range, frac_off is the fraction of the phase
	double on[n];
	double off[n];
	int num_on = (int)(n*frac_on);
	int num_off = (int)(n*frac_on);
	double on_0[num_on];
	double off_0[num_off];

	int index_on;
	int index_off;
	
	int i;
	double baseline = 0.0;
	double on_energy = 0.0;
	double off_energy = 0.0;

	index_on = def_on_pulse (n, in, frac_on);
	index_off = def_off_pulse (n, in, frac_on);

	get_region (n, index_on, in, on_0, frac_on);
	get_region (n, index_off, in, off_0, frac_on);

	// remove baseline
	for (i = 0; i < num_off; i++)
	{
		baseline += off_0[i];
	}
	baseline = baseline/num_off;

	for (i = 0; i < num_on; i++)
	{
		on[i] = (on_0[i]-baseline);
		off[i] = (off_0[i]-baseline);
		//s_norm[i] = (s[i]-baseline)/(s_peak-baseline);
	}

	// calculate on pulse and off pulse energy
	for (i = 0; i < num_on; i++)
	{
		on_energy += on[i];
		off_energy += off[i];
	}

	out[0] = on_energy;
	out[1] = off_energy;
	
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

	sub->p_multi = (double *)malloc(sizeof(double)*nchn*npol*nphase);
}

void demallocSub(subintegration *sub, pheader *phead)
{
	free(sub->p_multi);
	free(sub);
}

int histogram (double *data, int n, float *x, float *val, float low, float up, int step)
{
	int i,j,count;
	float width;
	float *temp;

	temp = (float*)malloc(sizeof(float)*(step+1));

	width = (up-low)/step;
	for (i=0; i<step; i++)
	{
		x[i] = low + i*width + width/2.0;
	}

	for (i=0; i<=step; i++)
	{
		temp[i] = low + i*width;
	}

	for (i=0; i<step; i++)
	{
		count = 0;
		for (j=0; j<n; j++)
		{
			if (data[j]>=temp[i] && data[j]<temp[i+1])
			{
				count += 1;
			}
		}
		val [i] = count;
	}

	free(temp);
	return 0;
}

int makePlot (char *fname, char *dname, int number)
{
	FILE *fpt;
	double *on_energy;
	double *off_energy;
	double on, off;
	int i, j;

	on_energy = (double *)malloc(sizeof(double)*number);
	off_energy = (double *)malloc(sizeof(double)*number);

	if ((fpt = fopen(fname, "r")) == NULL)
	{
		fprintf (stdout, "Can't open file\n");
		exit(1);
	}

	i = 0;
	while (fscanf(fpt, "%lf %lf", &on, &off) == 2)
	{
		on_energy[i] = on;
		off_energy[i] = off;
		i++;
	}

	if (fclose (fpt) != 0)
		fprintf (stderr, "Error closing\n");

	/////////////////////////////////////////////////

	float *xHis_on; // x axis of the histogram
	float *val_on;  // data value of the histogram
	float *xHis_off; // x axis of the histogram
	float *val_off;  // data value of the histogram
	int step = 500; // steps in the histogram

	//char caption[1024];
	//char text[1024];

	float max;

	// make histogram
	xHis_on = (float*)malloc(sizeof(float)*step);
	val_on = (float*)malloc(sizeof(float)*step);
	xHis_off = (float*)malloc(sizeof(float)*step);
	val_off = (float*)malloc(sizeof(float)*step);

	histogram (on_energy, number, xHis_on, val_on, -5.0, 5.0, step);
	histogram (off_energy, number, xHis_off, val_off, -5.0, 5.0, step);

	// plot 
	//cpgbeg(0,"/xs",1,1);
	cpgbeg(0,dname,1,1);

      	cpgsch(2); // set character height
      	cpgscf(2); // set character font

	// find the max
	max = find_max_value(step,val_on);
      	//cpgenv(-5,5,0,4500,0,1); // set window and viewport and draw labeled frame
      	cpgenv(-5,5,0,max+0.1*max,0,1); // set window and viewport and draw labeled frame

	//sprintf(caption, "%s", "Flux density histogram");
      	cpglab("Flux (mJy)","Number","");
	cpgbin(step,xHis_on,val_on,0);
	cpgsci(2);
	cpgbin(step,xHis_off,val_off,0);
	///////////////////////////////////////////////////////
	cpgend();
	////////////////////

	free(on_energy);
	free(off_energy);
	free(xHis_on);
	free(val_on);
	free(xHis_off);
	free(val_off);

	return 0;
}

float find_max_value (int n, float *s)
{
	int i;
	float *temp;

	temp = (float *)malloc(sizeof(float)*n);

	for (i = 0; i < n; i++)
	{
		temp[i] = s[i];
	}

	float a, b, c;
	for (i = 0; i < n-1; i++)
	{
		a = temp[i];
		b = temp[i+1];
															
		c = (a >= b ? a : b);
		
		temp[i+1] = c;
	
	}
	
	c = temp[n-1];
	free(temp);

	return c;
}

float find_min_value (int n, float *s)
{
	int i;
	float *temp;

	temp = (float *)malloc(sizeof(float)*n);

	for (i = 0; i < n; i++)
	{
		temp[i] = s[i];
	}

	float a, b, c;
	for (i = 0; i < n-1; i++)
	{
		a = temp[i];
		b = temp[i+1];
															
		c = (a <= b ? a : b);
		
		temp[i+1] = c;
	
	}
	
	c = temp[n-1];
	free(temp);

	return c;
}

