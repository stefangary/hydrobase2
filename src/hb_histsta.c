/*  hb_histsta.c
***********************************************************************
*  Usage:  hb_histsta infile_list  -O<outfile>
*                                 [-S]
*                                 [-T<tmin>/<tmax>/<tinc>]
*                                 [-X<xmin>/<xmax>/<xinc> &
*                                  -Y<ymin>/<ymax>/<yinc>]
*                                 [-D<directory>] [-Eextent] [-h]
*
*  Reads HydroBase data files in the infile_list
*  and bins the stations and scans according to
*  time (-T), depth (default), and horizontal
*  position (need both -X and -Y flags).  The
*  bin flags are used to set the bin ranges
*  (_min,_max) and increments (_inc).
*
*  Histogram binning can count the number of
*  stations in each bin (default) or the number
*  of data points (scans) in each bin (if -S
*  is present).
*
*  The other optional flags specify the root
*  directory (-D) and the extension (-E) of the
*  files.
*
*  Output goes to the output file (-O), all with
*  columnar format and no headers.  For the
*  depth and time histograms, format is two cols:
*      
*      (t|z)bin_center   number_of_stations
*
*  and for the position binning, output is in
*  three columns:
*
*      lon_bin_center lat_bin_center number_of_stations
*
**************************************************************************
*/ 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hydrobase.h"

#define  PRINT_MSG   1 
#define  MAXRANGES  200 
#define    DIR     ""
#define    EXTENT   ""

/* Prototypes for locally defined functions */

/* Print usage if user error. */
void print_usage(char *);

/* Get axis information from command line. */
void get_axis(char*, double *, double *, double *);

/* Get the index of the value in the depth array. */
int find_index(double, double *, int);

main(int argc, char **argv) {

  /* Define array pointers to store min/max depth,
   * temperature, and salinity at each depth bin. */

  int    nzbin;
  double *zbinedge, *zbinc; /* Bin edges and centers. */
  int    *zhist; /* Depth histogram (on bin centers). */
  int    *zobsf; /* Flag profile for determining whether there has
		  * been an observation at a given level or not. */

  int    ntbin;
  double *tbinedge, *tbinc;
  double date;
  int    *thist;

  int    nxbin, nybin;
  double  *xbinedge, *xbinc;
  double  *ybinedge, *ybinc;
  int    **xyhist;

  /* Local counter variables:
   * nstat = number of stations read (over all files)
   * nscan = number of scans read (over all files) */
  int    nstat, nscan;
  int    status, indx, indy;
  int    n, i, j;


  /* Flags to keep track of whether the user
   * has specified -X, -Y, -Z, or -T flags */
  int    xflag;
  int    yflag;
  int    tflag;
  int    zflag;
  int    sflag;

  /* Numbers to store the inputs of 
   * bin ranges and increments. */
  double tmin, tmax;
  double xmin, xmax;
  double ymin, ymax;
  double tinc;
  double xinc, yinc;

  /* nfiles = number of files on command line
   * curfile = current file number */
  int    nfiles, curfile;
  
  /* File ID of hydrobase format input file */
  int    infile;
  
  /* File ID of output histogram files (txt files) */
  FILE   *outfile;
  
  /* Pointers to strings for directory and
   * extension roots as well as a generic
   * string for reading command line. */
  char   *dir, *extent, *st;
  
  /* Input header and data storage */
  struct HYDRO_HDR hdr;   
  struct HYDRO_DATA data;

  /* Error flag for reading command line */
  int error;

  /* Are there command line arguments? */
  if (argc < 2) {
    print_usage(argv[0]);
    exit(1);
  }

  /* Set default values... */
  dir = DIR;
  extent = EXTENT;
  outfile = (FILE *)NULL;
  infile = 0;
  nstat = nscan = nfiles = 0; /* Initialize station, scan, and file counts */
  tflag = xflag = yflag = sflag = 0; /* Initialize command line arg flags to false */
  zflag = 1;
  curfile = 1;
  nfiles = 0;
   
  /* Initialize input data storage */
  for (i = 0; i < MAXPROP; ++i) {
    data.observ[i] = (double *) NULL;
  }
  hdr.prop_id = (int *)NULL;

  /* Parse command line arguments... */
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
               case 'D':
		 /* Set the input file directory. */
		 dir = &argv[i][2];
		 break;

               case 'E':
		 /* Set the input file extension. */
		 extent = &argv[i][2];
		 break;

               case 'X':
		 xflag = 1;
		 zflag = 0;
		 get_axis(&argv[i][2],&xmin,&xmax,&xinc);
		 /*fprintf(stderr,"X: [%lf : %lf : %lf]\n",xmin,xinc,xmax);*/
		 break;

               case 'Y':
		 yflag = 1;
		 zflag = 0;
		 get_axis(&argv[i][2],&ymin,&ymax,&yinc);
		 /*fprintf(stderr,"Y: [%lf : %lf : %lf]\n",ymin,yinc,ymax);*/
		 break;

               case 'O':
		 outfile = fopen(&argv[i][2],"w");
		 if (outfile == NULL) {
		   fprintf(stderr,"\nError opening %s for output.\n",&argv[i][2]);
		   exit(1);
		 }
		 break;

               case 'S':
		 /* We output #scans rather than #stations. */
		 sflag = 1;
		 break;

               case 'T':
		 tflag = 1;
		 zflag = 0;
		 get_axis(&argv[i][2],&tmin,&tmax,&tinc);
		 /*fprintf(stderr,"T: [%lf : %lf : %lf]\n",tmin,tinc,tmax);*/
		 break;

	       case 'h':
	          print_usage(argv[0]);
		  exit(0);

               default :
                   fprintf(stderr,"\nError parsing command line");
                   fprintf(stderr,"\n in particular: %s\n", argv[i]);
                   fprintf(stderr,"Use -h for help\n");
                   exit(1);
      }  /* end switch */
    }  /* end if */
    else {
      /* If none of the switches apply, then the
       * command line argument must be one of the
       * list of infiles.  Augment the input file
       * counter. */
      ++nfiles;
    }
  }  /* end for */
  
  /* Check that the user has specifed an outfile */
  if ( outfile == NULL ) {
    print_usage(argv[0]);
    fprintf(stderr,"\n\nYou must specify the outfile with -O.\n");
    exit(1);
  }

  /* Check that the user has specified both -X and -Y
   * if position binning is desired. */
  if ( xflag != yflag) {
    print_usage(argv[0]);
    fprintf(stderr,"\n\nYou must specify both -X and -Y flags for position binning.\n");
    exit(1);
  }

  /* Check that the user has not specified more than
   * one type of binning. */
  if ( (xflag == tflag) && (xflag == 1) ) {
    print_usage(argv[0]);
    fprintf(stderr,"\n\nYou may only do one type of histogram with each invocation.\n");
    exit(1);
  }

  /* Check that the user has specified at least one infile. */
  if ( nfiles < 1 ) {
    print_usage(argv[0]);
    fprintf(stderr,"\n\nYou must specify Hydrobase input files.\n");
    exit(1);
  }

  /* Initialize the bins depending on the type of
   * binning requested by the user. */
  /*-----------------------------------------------------*/

  /* Check if we are to do time binning */
  if ( tflag == 1 ) {
    /* Compute the number of bins - truncate round down. */
    ntbin = (int)((tmax - tmin)/tinc);

    /* Recompute the increment if the user gave an uneven inc. */
    tinc = (tmax - tmin)/((double)ntbin);

    fprintf(stderr,"T: [%lf : %lf : %lf]\n",tmin,tinc,tmax);

    tbinedge = (double *) calloc((ntbin+1), sizeof(double));
    tbinc    = (double *) calloc((ntbin), sizeof(double));
    thist    = (int *) calloc((ntbin), sizeof(int));

    /* Compute the bin edges. */
    tbinedge[0] = tmin;
    for ( i = 1; i <= ntbin; i++) {
      tbinedge[i] = tbinedge[i-1] + tinc;
      /*fprintf(stderr,"\nT bin edges are: %f",tbinedge[i]);*/
    }

    /* Compute the bin centers. */
    for ( i = 0; i < ntbin; i++) {
      tbinc[i] = (tbinedge[i] + tbinedge[i+1])/2.0;
      /*fprintf(stderr,"\nT bin centers are: %f",tbinc[i]);*/
      thist[i] = 0;
    }
  } /* End of time binning check. */

  /*-----------------------------------------------------*/

  /* Check if we are to do position binning */
  if ( xflag == 1 && yflag == 1 ) {
    /* Compute the number of bins - truncate round down. */
    nxbin = (int)((xmax - xmin)/xinc);
    nybin = (int)((ymax - ymin)/yinc);

    /* Recompute the increment if the user gave an uneven inc. */
    xinc = (xmax - xmin)/((double)nxbin);
    yinc = (ymax - ymin)/((double)nybin);

    fprintf(stderr,"X: [%lf : %lf : %lf]\n",xmin,xinc,xmax);
    fprintf(stderr,"Y: [%lf : %lf : %lf]\n",ymin,yinc,ymax);

    xbinedge = (double *) calloc((nxbin+1), sizeof(double));
    ybinedge = (double *) calloc((nybin+1), sizeof(double));
    xbinc    = (double *) calloc((nxbin), sizeof(double));
    ybinc    = (double *) calloc((nybin), sizeof(double));
    xyhist   = (int **) calloc((nxbin), sizeof(int *));
    for ( i = 0; i < nxbin; i++) {
      xyhist[i] = (int *) calloc((nybin), sizeof(int));
    }

    /* Compute the bin edges. */
    xbinedge[0] = xmin;
    for ( i = 1; i <= nxbin; i++) {
      xbinedge[i] = xbinedge[i-1] + xinc;
      /*fprintf(stderr,"\nX bin edges are: %f",xbinedge[i]);*/
    }

    ybinedge[0] = ymin;
    for ( i = 1; i <= nybin; i++) {
      ybinedge[i] = ybinedge[i-1] + yinc;
      /*fprintf(stderr,"\nX bin edges are: %f",ybinedge[i]);*/
    }

    /* Compute the bin centers. */
    for ( i = 0; i < nxbin; i++) {
      xbinc[i] = (xbinedge[i] + xbinedge[i+1])/2.0;
      /*fprintf(stderr,"\nX bin centers are: %f",xbinc[i]);*/
    }

    for ( i = 0; i < nybin; i++) {
      ybinc[i] = (ybinedge[i] + ybinedge[i+1])/2.0;
      /*fprintf(stderr,"\nY bin centers are: %f",ybinc[i]);*/
    }

    /* Initialize the position histogram with zero counts. */
    for ( i = 0; i < nxbin; i++) {
      for (j = 0; j < nybin; j++) {
	xyhist[i][j] = 0;
      }
    }
  }

  /*-----------------------------------------------------*/

  /* Otherwise, the default is to do depth binning */
  if ( zflag == 1) {

    /* Manually define the depth binning. */
    nzbin = 38;
    zbinedge = (double *) calloc((nzbin+1), sizeof(double));

    zbinedge[0] = 0.0;
    zbinedge[1] = 10.0;
    zbinedge[2] = 20.0;
    zbinedge[3] = 30.0;
    zbinedge[4] = 50.0;
    zbinedge[5] = 75.0;
    zbinedge[6] = 100.0;
    zbinedge[7] = 125.0;
    zbinedge[8] = 150.0;
    zbinedge[9] = 200.0;
    zbinedge[10] = 250.0;
    zbinedge[11] = 300.0;
    zbinedge[12] = 400.0;
    zbinedge[13] = 500.0;
    zbinedge[14] = 600.0;
    zbinedge[15] = 700.0;
    zbinedge[16] = 800.0;
    zbinedge[17] = 900.0;
    zbinedge[18] = 1000.0;
    zbinedge[19] = 1100.0;
    zbinedge[20] = 1200.0;
    zbinedge[21] = 1300.0;
    zbinedge[22] = 1400.0;
    zbinedge[23] = 1500.0;
    zbinedge[24] = 1750.0;
    zbinedge[25] = 2000.0;
    zbinedge[26] = 2500.0;
    zbinedge[27] = 3000.0;
    zbinedge[28] = 3500.0;
    zbinedge[29] = 4000.0;
    zbinedge[30] = 4500.0;
    zbinedge[31] = 5000.0;
    zbinedge[32] = 5500.0;
    zbinedge[33] = 6000.0;
    zbinedge[34] = 6500.0;
    zbinedge[35] = 7000.0;
    zbinedge[36] = 7500.0;
    zbinedge[37] = 8000.0;
    zbinedge[38] = 8500.0;

    /*
    for ( i = 0; i <= nzbin; i++ ) {
      fprintf(stderr,"\nDepth bin edges are: %f",zbinedge[i]);
    }
    */

    /* Initialize/compute depth bin centers and histogram. */
    zbinc = (double *) calloc(nzbin, sizeof(double));
    zhist = (int *) calloc(nzbin, sizeof(int));
    zobsf = (int *) calloc(nzbin, sizeof(int));

    for ( i = 0; i < nzbin; i++) {
      zbinc[i] = (zbinedge[i] + zbinedge[i+1])/2.0;
      zhist[i] = 0;
      zobsf[i] = 0;
    }
    
    /*
    for ( i = 0; i < nzbin; i++ ) {
      fprintf(stderr,"\nDepth bin centers are: %f",zbinc[i]);
    }
    */
  } /* End of depth binning check. */
  
  /*-----------------------------------------------------*/

  /* For each file in list_of_infiles, */
  do {
    
    /* Open the current infile. */
    infile = open_hydro_file(dir, argv[curfile], extent, PRINT_MSG);
    if (infile  < 0)
      /* We do not have a valid infile.x */
      goto NEXTFILE;
    
    /* Loop for each station */
    while ((status = get_station(infile, &hdr, &data)) == 0) {
      /* get_station returns 0 for sucessful read. */
      /* Input data are in structures hdr and data. */

      /* Augment the number of stations opened. */
      ++nstat;

      /* Check binning type. */
      /*-----------------------------------------------------*/

      if ( tflag == 1 ) {
	/* We are to do time binning. */

	/* Calculate the fractional date (in years) of this profile. */
	date = hdr.year + (hdr.month-1)/12 + (hdr.day-1)/365;

	/* Find the time bin that bounds the date of this profile. */
	indx = find_index(date, tbinedge, ntbin);

	if ( sflag == 0 ) {
	  /* We compute the number of stations. */
	  thist[indx]++;
	}

	if ( sflag == 1 ) {
	  /* We compute the number of measurements. */
	  thist[indx] += hdr.nobs;
	}
      } /* End of time binning check. */

      /*-----------------------------------------------------*/

      if ( xflag == 1 && yflag == 1 ) {
	/* We are to do position binning. */
	
	/*fprintf(stderr,"\nStation position: %f E, %f N",hdr.lon,hdr.lat);*/

	/* Find the x and y bins that hold this profile. */
	indx = find_index(hdr.lon, xbinedge, nxbin);
	indy = find_index(hdr.lat, ybinedge, nybin);

	/*fprintf(stderr,"\nIndex in grid: %d E, %d N",indx,indy);*/

	if ( sflag == 0 ) {
	  /* We compute the number of stations. */
	  xyhist[indx][indy]++;
	}

	if ( sflag == 1 ) {
	  /* We compute the number of measurements. */
	  xyhist[indx][indy] += hdr.nobs;
	}
      }
      
      /*-----------------------------------------------------*/

      if ( zflag == 1) {
	/* We are binning stations by depth. */

	/* For each station, reinitialize the observation flag. */
	for (i = 0; i < nzbin; i++) {
	  zobsf[i] = 0;
	}

	/* Loop for each data scan/depth in profile, */
	for (i = 0; i < hdr.nobs; ++i) {

	  /* Find the depth bin that bounds this depth */
	  indx = find_index(data.observ[(int)DE][i], zbinedge, nzbin);

	  /* Count the observations at this level. */
	  zobsf[indx]++;

	  /* If counting number of observations: */
	  if( sflag == 1 ) {
	  
	    /* Augment the counter in the histogram for that depth bin. */
	    zhist[indx] += zobsf[indx];
	  }

	} /* End of loop over all depths. */

	if( sflag == 0 ) {
	  /* Loop over all levels and test whether or
	   * not there were observation at a given level. */
	  for ( i = 0; i < nzbin; i++) { 
	    if( zobsf[i] > 0 ) {
	      zhist[i]++;
	    }
	  }
	}
      } /* End of depth binning check. */

      /*-----------------------------------------------------*/

    } /* End of loop over each station. */
    
    report_status(status, stderr);
    if ( nfiles ) close(infile);

  NEXTFILE:
    ;
  } while (curfile++ < nfiles); /* End of loop over all input files */
  
  /* Print information about total number of stations opened. */
  fprintf(stderr,"Total number of stations: %d",nstat);

  /* Write to the output file depending on binning type. */
  /*-----------------------------------------------------*/

  if ( tflag == 1 ) {
    /* We are to do time binning. */
    for ( i = 0; i < ntbin; i++) {
      fprintf(outfile,"%f %d\n",tbinc[i],thist[i]);
    }
  } /* End of time binning check*/

  /*-----------------------------------------------------*/

  if ( xflag == 1 && yflag == 1 ) {
    /* We are to do position binning. */
    for ( i = 0; i < nxbin; i++) {
      for ( j = 0; j < nybin; j++) {
	fprintf(outfile,"%f %f %d\n",xbinc[i],ybinc[j],xyhist[i][j]);
      }
    }
  } /* End of position binning check. */
      
  /*-----------------------------------------------------*/

  if ( zflag == 1) {
    /* We are binning stations by depth. */
    for ( i = 0; i < nzbin; i++) {
      fprintf(outfile,"%f %d\n",zbinc[i],zhist[i]);
    }
  } /* End of depth binning check. */

  /*-----------------------------------------------------*/

  fprintf(stderr,"\nEnd of %s\n", argv[0]);
  exit(0);
}  /* End main */

/****************************************************************************/

void print_usage(char *program)
{
  fprintf(stderr,"\n --------------------------------------------------------");
  fprintf(stderr,"\n Usage: %s list_of_infile(s) -O<outfile>",program);
  fprintf(stderr,"\n                             [-T<tmin>/<tmax>/<tinc>]");
  fprintf(stderr,"\n                             [-X<xmin>/<xmax>/<xinc> &");
  fprintf(stderr,"\n                              -Y<ymin>/<ymax>/<yinc>]");
  fprintf(stderr,"\n                             [-D<directory>] [-Eextent] [-S] [-h]");
  fprintf(stderr,"\n");
  fprintf(stderr,"\n List of infiles should be first argument(s) ");
  fprintf(stderr,"\n   -O  : specifies output file for histograms.  ");
  fprintf(stderr,"\n        ex: -Odepth_histogram.txt ");
  fprintf(stderr,"\n   -T  : specifies binning, in years, for time histogram.  ");
  fprintf(stderr,"\n        ex: -T1980/2009/5 ");
  fprintf(stderr,"\n   -X  : specifies longitude binning, in degrees E, for position histogram.  ");
  fprintf(stderr,"\n        ex: -X-85.3/0.2/1 ");
  fprintf(stderr,"\n   -Y  : specifies latitude binning, in degrees N, for position histogram.  ");
  fprintf(stderr,"\n        ex: -Y0.0/65.2/1 ");
  fprintf(stderr,"\n   -D  : specifies dir for input files (default is ./) ");
  fprintf(stderr,"\n        ex: -D../data/ ");
  fprintf(stderr,"\n   -E  : specifies file extent for input files (default is no extent)");  
  fprintf(stderr,"\n        ex: -E.dat ");
  fprintf(stderr,"\n   -S  : if listed, prints #scans, if missing, #stations.");
  fprintf(stderr,"\n   -h  : help -- prints this message. ");
  fprintf(stderr,"\n");
  fprintf(stderr,"\n This program will compute depth (default), time (-T)");
  fprintf(stderr,"\n and position (both -X and -Y required) histograms of");
  fprintf(stderr,"\n the all the hydrobase stations in the list_of_infile(s).");
  fprintf(stderr,"\n");
  fprintf(stderr,"\n --------------------------------------------------------");
  fprintf(stderr,"\n");
  return;
} /* End print_usage() */


/****************************************************************************/
/* This function uses the infromation in the
 * text string pointed to by argin to determine
 * the values of vmin, vmax, and vinc.  The
 * expected input format is:
 *
 * </>min/max/inc
 *
 * where the first slash is optional. */

void get_axis(char *argin, double *vmin, double *vmax, double *vinc) {

  int error;

  if (*argin == '/')
    ++argin;

  error = (sscanf(argin,"%lf", vmin) != 1);
  while(*(argin++) != '/')
    ;

  error = (sscanf(argin,"%lf", vmax) != 1);
  while(*(argin++) != '/')
    ;

  error = (sscanf(argin,"%lf", vinc) != 1);

  /*fprintf(stderr,"Axis: [%lf : %lf : %lf]\n",*vmin,*vinc,*vmax);*/

  if (*vmin >= *vmax) {
    fprintf(stderr,"Lower bound must be less than upper bound.\n");
    exit(1);
  }

  if ( (*vmax - *vmin) <= *vinc ) {
    fprintf(stderr,"Axis increment must be less than bound separation.\n");
    exit(1);
  }
}
/****************************************************************************/

/****************************************************************************/

int find_index(double sortme, double *binedges, int nbin) {

/* Determines the index (returned value) corresponding
 * the bin (defined by binedges) into which the value
 * of sortme falls.  Also, give the number of bins
 * (nbins) to prevent overrun.
 *
 * Returns the index number, or -1 if the value is not in
 * the range of the bins. */

   int i;

   i = -1;
   while ( ++i < nbin ) {
     /* The maximum value of i is nedges - 2 since
      * nbin = nedges - 1.  This is equal to the last
      * index of the bincenters since all arrays start
      * counting by 0. For example:
      *
      * nedges = 39
      * nbins  = 38
      *
      * Then imax = 37, which works because the bin
      * centers array is valid from [0,37], not [1,38]!
      * Simiarly, the maximum value accessed in binedges
      * is [0,38] because of the i+1 term. */
     if (sortme >= binedges[i] && sortme < binedges[i+1] )
        return (i);
   }

   /* If the bin was not found, return error flag value. */
   return (-1);
   
}  /* End find_index() */

/****************************************************************************/
