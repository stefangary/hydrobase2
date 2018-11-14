/*  hb_interpolate.c

................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth  Curry
                             Woods Hole Oceanographic Inst
                             1993
                             updated to ANSI Mar 2000
................................................................................
____________________________________________________________________________
  USAGE:  

hb_interpolate filename(_roots) -P<prs_int> [-D<dirname>] [-E<file_extent>]

  list of filenames MUST be first argument or input expected from stdin

 -P : pressure interval for output series.

 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum

____________________________________________________________________________
                                                   ____________________________________________________________________________
*/


#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"


/* input file pathnames */

#define    EXTENT   ""
#define    DIR      ""

   /* global variables to store station ... */

struct HYDRO_DATA data, newdata;
struct HYDRO_HDR hdr;
int prsint;

  /* prototypes for locally defined functions */
  
void print_usage(char *);
void get_hydro_data(int );


main (argc, argv)
int   argc;
char *argv[];
{
   int     index, nprops = 0;
   int      i, count;
   int     curfile = 1, nfiles = 0; 
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   char    *extent, *dir, *st;
   int     infile;


/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    prsint = 0;
    error = 0;
 
/* initialize these ... */

   for (i = 0; i < MAXPROP; ++i) {
      data.observ[i] = (double *) NULL;
      newdata.observ[i] = (double *) NULL;
   }
   hdr.prop_id = (int *) NULL;



/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'D':                   /* get input dir */
                        dir = &argv[i][2];
                        break;
               case 'E':                    /* get file extent */
                        extent = &argv[i][2];
                        break;

               case 'P':
                        
                       error = (sscanf(&argv[i][2],"%d", &prsint)==1) ? 0:1;
                       break;

	       case 'h':  
	               print_usage(argv[0]);
	               exit(0);

               default:
                        error = 1;

          }    /* end switch */

          if (error ) {
             fprintf(stderr,"\nError parsing command line args.\n");
             fprintf(stderr,"     in particular: '%s'\n", argv[i]);
             exit(1);
          }

       }  /* end if */

       else  {
           ++nfiles;
       }
   }  /* end for */

   if ( !prsint) {
       fprintf(stderr,"\nYou must specify input file(s) and a prs_int.\n");
       exit(1);
   }
   if (!nfiles) {
       fprintf(stderr,"\nhb_interpolate: Expecting input from stdin ... ");
       infile = STDIN;
   }

/* loop for each input file */

   do {

     if (nfiles > 0) {
       infile = open_hydro_file(dir, argv[curfile], extent, print_msg);
       if (infile < 0)
          goto NEXTFILE;
     }

     
            /* read each file completely */

 
     get_hydro_data(infile);

NEXTFILE:

     if (nfiles > 0) 
         close(infile);

   } while (curfile++ < nfiles );

   fprintf(stderr,"End of hb_interpolate.\n");
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
 fprintf(stderr,"\n%s reads each station in a HydroBase file and outputs a new file with observations at the specified pressure interval\n", program);

   fprintf(stderr,"\nUsage:  %s filename_root(s) -P<prsint> [-D<dirname>] [-E<file_extent>] [-h]", program);

   fprintf(stderr,"\n\n  List of filenames MUST be first argument");
   fprintf(stderr,"\n  or station input is expected from stdin.");
   fprintf(stderr,"\n    -P  : pressure interval for output array");
   fprintf(stderr,"\n   [-D] : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n   [-h] : help...... prints this message.");
   fprintf(stderr,"\n\n NOTE: hb_interpolate will not extrapolate beyond");
   fprintf(stderr,"\n       upper and lower bounds of input station.\n\n");
   return;
}

/*****************************************************************************/
void get_hydro_data(int file)

   /*  Reads each station in a HydroBase file and outputs observations at
      the specified pressure interval.    
   */
{
   int error, i, ii, j, npts, n;
   double *tmpp, *tmpx;
/* read each station in file ... */

/* Initialize temporary holders of data.  calloc will
 * zero out the variables as well as allocate memory.
 * Also, original HB2 had 6000 here, but for really
 * deep casts 6000m+ with high resolution pressure,
 * this limited hb_interpolate and caused segfaults.
 * Upped this number to 11000, sfg.*/
    tmpp = (double *) calloc(11000,sizeof(double));
    tmpx = (double *) calloc(11000,sizeof(double));

    while ((error = get_station(file, &hdr, &data)) == 0) {
    
      /* Loop over all possible properties in station...*/
      for (ii = 0; ii < MAXPROP; ++ii) {
	if  (data.observ[ii] != NULL ) {
	  /* If this property is present, allocate space for
	   * interpolated values. Also note that this value
	   * needed to be increased from 6000 to a larger
	   * number to avoid segfaults on deep casts at
	   * high pressure resolution.*/
	  newdata.observ[ii] = (double *) calloc(11000, sizeof (double));

	  /* For the case of pressure, set the shallowest
	   * (least) pressure to the nearest, deepest good pressure
	   * interval (sometimes -9 for pressure can throw things off)
	   * and add the interval to each successive level.*/
	  if (ii == (int)PR) {

	    /* Copy only the good pressure values to a temporary array*/
	    n = 0;
	    for (j = 0; j < hdr.nobs; ++j) {
	      if (data.observ[(int)PR][j] > -8.9) {
		tmpp[n] = data.observ[(int)PR][j];
		++n;
	      }
	    }

	    /* Older lines of code here do not use the tmpp variable
	     * (the good pressure values) which means that if there
	     * are spurious -9 or HB_MISSING in the pressure series
	     * as there are sometime from the netcdf files, then
	     * the whole pressure interpolation process can be thrown
	     * off and profiles will be cut short.
	    newdata.observ[ii][0] = NINT(data.observ[ii][0] / prsint) * prsint;
	    npts = 1;
                  
	    while ((newdata.observ[ii][npts] = newdata.observ[ii][npts-1] + prsint ) <= data.observ[ii][hdr.nobs-1])
	      ++npts;
	    */

	    /* Essentially the same operation above, but using tmpp,
	     * the HB_MISSING filtered pressure, instead.*/
	    newdata.observ[ii][0] = NINT(tmpp[0] / prsint) * prsint;
	    npts = 1;
                  
	    while ((newdata.observ[ii][npts] = newdata.observ[ii][npts-1] + prsint ) <= tmpp[n-1])
	      ++npts;
	  }

	  /* For all other properties, they are interpolated.*/
	  else {
	    n = 0;
	    for (j = 0; j < hdr.nobs; ++j) {
	      /* Note which data points are flagged as missing.  This
	       * loop originally had:
	      if (data.observ[ii][j] > -8.9)
		tmpp[n] = data.observ[(int)PR][j];
	      tmpx[n] = data.observ[ii][j];
	      ++n;
	      */
	      /* But I think it's bad.  Should be all bracketed:*/
	      if (data.observ[ii][j] > -8.9) {
		tmpp[n] = data.observ[(int)PR][j];
		tmpx[n] = data.observ[ii][j];
		++n;
	      }
	      /* Tests with the two versions verify that this was
	       * indeed a bug because otherwise the bad values are
	       * included as part of the interpolation.  We want to
	       * skip over them.*/
	    }
	    for (j = 0; j < npts; ++j) {
	      newdata.observ[ii][j] = hb_linterp(newdata.observ[(int)PR][j], tmpp, tmpx,n);
	      if (newdata.observ[ii][j] < -8.9)
		newdata.observ[ii][j] = HB_MISSING;
	    }
	    
	  }
	}
      }
	   
      /* output the station */
      hdr.nobs = npts;
      newdata.nobs = npts;
      newdata.nprops = hdr.nprops;	
      
      if (hdr.nobs > 0 )
	write_hydro_station(STDOUT, &hdr, &newdata);
	    
 
      for (i = 0; i < MAXPROP; ++i) {
	if (data.observ[i] != NULL) {
          free(data.observ[i]);
          data.observ[i] = NULL;
	}
	
	if (newdata.observ[i] != NULL) {
          free(newdata.observ[i]);
          newdata.observ[i] = NULL;
	}
       
      } 
      free(hdr.prop_id); 
      hdr.prop_id = NULL;
       
    }  /* end while */

    free(tmpp);
    free(tmpx);   
    
    if (error > 0)
      report_status(error, stderr);
    
    return;
    
} /* end get_hydro_data() */

