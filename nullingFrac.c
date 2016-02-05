// calculate nulling fraction 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "nullingFracLib.h"

int main (int argc, char *argv[])
{
	int h,i,j,k;

	pheader *header;
	header = (pheader *)malloc(sizeof(pheader));

	fitsfile *fp;

	subintegration *sub;
	sub = (subintegration *)malloc(sizeof(subintegration));

	//////////////////////////////////////////////////////
	char oname[128];   // name of output .tim
	char dname[128];   // name of output .tim
	double on_frac; // on pulse fraction

	int index, n;
	for (i=0;i<argc;i++)
	{
		if (strcmp(argv[i],"-f") == 0)
		{
			index = i + 1;
			n = 0;
			while ( (index + n) < argc && strcmp(argv[index+n],"-o") != 0 && strcmp(argv[index+n],"-dev") != 0 && strcmp(argv[index+n],"-frac") != 0 )
			{
				n++;
			}
		}
		else if (strcmp(argv[i],"-o")==0)
		{
			strcpy(oname,argv[++i]);
		}
		else if (strcmp(argv[i],"-dev")==0)
		{
			strcpy(dname,argv[++i]);
		}
		else if (strcmp(argv[i],"-frac")==0)
		{
			on_frac = atof(argv[++i]);
			printf ("On pulse fraction is %.1f: \n", on_frac);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////
	// open file to write toa 
	FILE *fpt;
	if ((fpt = fopen(oname, "w+")) == NULL)
	{
		fprintf (stdout, "Can't open file\n");
		exit(1);
	}
	
	/////////////////////////////////////////////////////////////////////////////////
	// start to deal with different data file
	int number = 0;
	for (k = index; k < index + n; k++)
	{
		// get the data file name
		strcpy(sub->fname,argv[k]);
		printf ("%s\n", sub->fname);

		// open psrfits file
		fp = openFitsFile(sub->fname);

		// read header info
		loadPrimaryHeader(fp,header);
		closeFitsFile(fp);
		////////////////////////////////////////////////////

		int nphase;
		int nchn;
		int nsub;
		int npol;
	
		nchn = header->nchan; 
		npol = header->npol;  // npol should be one
		nsub = header->nsub; 
		nphase = header->nbin; 	

		//printf ("%d\n", nchn);
		////////////////////////////////////////////////

		double p_temp[nphase];
		double energy[2];

		initialiseSub(sub, header);

		// start to derive toa from different subint
		for (h = 1; h <= nsub; h++)
		{
			sub->indexSub = h;
			// read profiles from data file
			read_prof(sub, header);

			// get pulse period of this subintegration

			// start to derive toas for different channels
			for (i = 0; i < nchn; i++)
			{
				sub->indexChn = i;
				for (j = 0; j < nphase; j++)
				{
					p_temp[j] = sub->p_multi[i*nphase + j];
					//printf ("%lf\n", p_temp[j]);
				}
				cal_pulseEnergy (p_temp, on_frac, nphase, energy);
				fprintf (fpt, "%lf %lf\n", energy[0], energy[1]);
				number++;
			}

		}
	}

	if (fclose (fpt) != 0)
		fprintf (stderr, "Error closing\n");

	free(header);
	demallocSub(sub, header);
	
	makePlot (oname, dname, number);

	return 0;
}

