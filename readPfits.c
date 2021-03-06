#include <stdio.h>
#include <math.h>
#include <string.h>
#include "fitsio.h"
#include "readPfits.h"

// Routines for reading PSRFITS files

void closeFitsFile(fitsfile *fp)
{
  int status=0;
  fits_close_file(fp,&status);
  fits_report_error(stderr,status);
  if (status)
    {
      printf("Error closing file\n");
      exit(1);
    }

}

fitsfile * openFitsFile(char *fname)
{
  fitsfile *fp;
  int status=0;
  fits_open_file(&fp,fname,READONLY,&status);
  fits_report_error(stderr,status);
  if (status)
    {
      printf("Error opening file >%s<\n",fname);
      exit(1);
    }
  return fp;
}

void loadPrimaryHeader(fitsfile *fp,pheader *phead)
{
  int status=0;
  int nkey=-1;
  int morekeys=-1;
  int i;
  char keyname[128],val[128],comment[128];

  fits_get_hdrspace(fp,&nkey,&morekeys,&status);

  phead->nhead = nkey;
  // Allocate memory
  phead->keyname = (char **)malloc(sizeof(char *)*nkey);
  phead->val = (char **)malloc(sizeof(char *)*nkey);
  phead->comment = (char **)malloc(sizeof(char *)*nkey);

  for (i=0;i<nkey;i++)
	{
		phead->keyname[i] = (char *)malloc(sizeof(char)*128);
		phead->val[i] = (char *)malloc(sizeof(char)*128);
		phead->comment[i] = (char *)malloc(sizeof(char)*128);
	}

  // Complete allocating memory

  for (i=1;i<=nkey;i++)
	{
		fits_read_keyn(fp,i+1,phead->keyname[i-1],phead->val[i-1],phead->comment[i-1],&status);
		if (strcmp(phead->keyname[i-1],"OBSFREQ")==0)
			sscanf(phead->val[i-1],"%f",&(phead->freq));
		else if (strcmp(phead->keyname[i-1],"OBSFREQ_SSB")==0)
			sscanf(phead->val[i-1],"%f",&(phead->freqSSB));
		else if (strcmp(phead->keyname[i-1],"STT_IMJD")==0)
			sscanf(phead->val[i-1],"%d",&(phead->imjd));
		else if (strcmp(phead->keyname[i-1],"STT_SMJD")==0)
			sscanf(phead->val[i-1],"%f",&(phead->smjd));
		else if (strcmp(phead->keyname[i-1],"STT_OFFS")==0)
			sscanf(phead->val[i-1],"%f",&(phead->stt_offs));
		else if (strcmp(phead->keyname[i-1],"OBSBW")==0)
			sscanf(phead->val[i-1],"%f",&(phead->bw));
	}
  // Read specific parameters
  fits_read_key(fp,TSTRING,(char *)"OBS_MODE",phead->obsMode,NULL,&status);
  fits_read_key(fp,TSTRING,(char *)"SRC_NAME",phead->source,NULL,&status);
  if (status)
	{
		fits_report_error(stderr,status);
		exit(1);
	}

  // Now load information from the subintegration table
  fits_movnam_hdu(fp,BINARY_TBL,(char *)"SUBINT",1,&status);
  if (status)
	{
		printf("No subintegration table\n");
		status=0;
	}
  else
	{
		fits_read_key(fp,TINT,(char *)"NAXIS2",&(phead->nsub),NULL,&status);
		if (status)
		{
			printf("Reading naxis2\n");
			fits_report_error(stderr,status);
			exit(1);
		}     
		fits_read_key(fp,TINT,(char *)"NCHAN",&(phead->nchan),NULL,&status);
    
		if (status)
		{
			printf("Reading nchan\n");
			fits_report_error(stderr,status);
			exit(1);
		}
      
		/*      fits_read_key(fp,TFLOAT,(char *)"ZERO_OFF",&(phead->zeroOff),NULL,&status);
    if (status)
		{
		  printf("Reading zero_off\n");
		  fits_report_error(stderr,status);
		  phead->zeroOff = 0;
		  status=0;
		  printf("Complete reading zero_off\n");
		}  */
				/*
  	    printf("COMMENTED OUT READING ZERO_OFF AND NBITS IN PFITS.C -- PUT BACK -- PROBLEM WITH SOME FILES\n"); 
  	    fits_read_key(fp,TINT,(char *)"NBITS",&(phead->nbits),NULL,&status);
  	    if (status)
		{
		  printf("Reading nbits\n");
		  fits_report_error(stderr,status);
		  exit(1);
		}
		*/

		fits_read_key(fp,TINT,(char *)"NPOL",&(phead->npol),NULL,&status);
    
		if (status)
		{
			printf("Reading npol\n");
			fits_report_error(stderr,status);
			exit(1);
		}
      
		fits_read_key(fp,TINT,(char *)"NSBLK",&(phead->nsblk),NULL,&status);
		if (status)
		{
			printf("Reading nsblk\n");
			fits_report_error(stderr,status);
			exit(1);
		}

		fits_read_key(fp,TINT,(char *)"NBIN",&(phead->nbin),NULL,&status);
		if (status)
		{
			printf("Reading nbin\n");
			fits_report_error(stderr,status);
			exit(1);
		}

		//      printf("nbin = %d (%d)\n",phead->nbin,status);
		fits_read_key(fp,TFLOAT,(char *)"CHAN_BW",&(phead->chanbw),NULL,&status);
    
		if (phead->chanbw < 0 && phead->bw > 0)
			phead->bw*=-1;
      
		fits_read_key(fp,TFLOAT,(char *)"TBIN",&(phead->tsamp),NULL,&status);
	}
  
	fits_movnam_hdu(fp,BINARY_TBL,(char *)"PSRPARAM",1,&status);
  if (status)
	{
		printf("No PSRPARM table\n");
		status=0;
	}
  else
	{
		int len,i,colnum;
		char **line,str1[1024],str2[1024];
		char nval[128]="UNKNOWN";
		int anynul=0;
		float tt;
   
		fits_read_key(fp,TINT,(char *)"NAXIS2",&len,NULL,&status);
   
		fits_get_colnum(fp,CASEINSEN,(char *)"PARAM",&colnum,&status);
		if (status) 
		{
			printf("Unable to find data in the psrparam table in FITS file\n");
			exit(1);
		}

		line = (char **)malloc(sizeof(char *));
		line[0] = (char *)malloc(sizeof(char)*1024); 

		for (i=0;i<len;i++)
		{
			fits_read_col_str(fp,colnum,i+1,1,1,nval,line,&anynul,&status);
	  
			if (sscanf(line[0],"%s %s",str1,str2)==2)
	    {
				if (strcasecmp(str1,"DM")==0)
					sscanf(str2,"%f",&(phead->dm));
	      if (strcasecmp(str1,"F0")==0)
				{
					sscanf(str2,"%f",&tt);
					phead->period = 1.0/tt;
				}
			}
	  //	  printf("Read: %s\n",line[0]);
		}
      //      printf("Lenght = %d\n",len);
  
		free(line[0]);
		free(line);
	}
}
