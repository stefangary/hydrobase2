/*   hb_statchk_ts.c
.
. USAGE:  hb_statchk_ts -Iinfile -Ooutfile -Sstatfile [-Ddeep_statfile] 
.                       -N#stderr_T-S [-M#stderr_surf/rhomax]
.                      [-Bbadfile] [-Llogfile]  [-Pplotfile_rootname]
.
.   Examines each data point in a file and eliminates points which are 
.   statistically distant from the t-s relation defined for the bin.  
.
.   Linear theta-s curves and a standard error of the estimate
.   are computed  for density (sigma) bins by the program
.   hb_statfit_ts.  These coefficients determine an acceptable
.   range of values in each bin  +- n standard errors from the
.   line.  If a t-s point falls outside this range, the entire
.   scan is eliminated (because a scan must have both to place
.   it in a density bin). A point falling in no density bin or
.   in a density bin containing less than 3 pts is eliminated.
.   A station bearing >10 (or >50%) bad scans is completely
.   eliminated.   The number of eliminated points is tallied
.   as a percentage of the total number of scans and written
.   to the logfile (or stdout).
.
. Authors
. Ruth Curry - original code
. Stefan Gary - added comments
*/

#include <stdio.h>
#include <stdlib.h>   /*sfg added to deal with exit and malloc warnings*/
#include <math.h>
#include <string.h>
#include "hydrobase.h"

#define   PRINT_MSG    1


#define   RED         2       /* pen colors */
#define   BLACK       1
#define   SQU         0       /* display syml set */
#define   TRI         2
#define   X           4

#define   FLAG   -9.9
#define   PERCENT_LIM .65     /* fraction of station allowed bad before
                                 declaring entire station bad */

/* ***********  global variables  ************ */

/* This structure will store all the
 * pertinant statistical parameters
 * about a given density bin.  This
 * information is loaded from the
 * statfile or the deepfile. */
struct coefficients {
  char   xvar;         /* identifies independent variable: T or S */
  float  m1, b1, e1;   /* slope, intercept, std error for t-s */
  int    n1;           /* number of points in density bin */
} ; 

struct coefficients *coeff;
float *sigmin, *sigmax;
float percent_lim;

/* Number of density bins from statfile. */
int  nbins;

float  rhomax, nstderr_surf, nstderr;
FILE *logfile, *tsfile;
int  plotflag;
float *plot_s, *plot_t, *plot_sg, *plot_p;
int  *tspen, *sym;
char rflag;
/************************************/

/* Prototypes for internally defined functions */

/* Print help to the user. */
void    print_usage(char *);

/* Load information from the statfile. */
int load_coeffs(FILE *);

/* Any high density (deep waters) bins
 * have their coeffs overwritten by
 * corresponding values in a separate
 * deepfile. */
int load_deep_coeffs(FILE *);

/* Check that each scan within a station
 * fits within the statistical parameters
 * loaded from the statfile and/or
 * deepfile.  Returns the good scans and
 * updates a pointer to the bad scans. */
struct HYDRO_DATA *check_scans(int, struct HYDRO_HDR *, int *,  struct HYDRO_DATA **);

/************************************/
main (int argc, char **argv)
{
   int     i;
   int     badsta, badscan, badts;
   int     totscan, totsta;
   int     status; 
   int     infile, outfile, badfile;
   char    iflag, oflag, sflag, nflag, dflag, bflag;
   char   *str, *statfilename, *deepfilename; 
   char   fname[80];
   struct HYDRO_HDR hdr;
   struct HYDRO_DATA  *newdata, *baddata;
   FILE   *statfile, *deepfile;


   /* Check that there are command line arguments. */
   if (argc < 4) {
      print_usage(argv[0]);
      exit(1);
   }

   /* Set default values and flags. */
   rflag = iflag = oflag = sflag = dflag = bflag = 0;
   plotflag = 0;
   logfile = stdout;
   percent_lim = PERCENT_LIM;

   /* Initialize list of prop IDs. */
   hdr.prop_id = (int *) NULL;
   
   /* Parse command line arguments... */
   for (i = 1; i < argc; i++) {
     if (argv[i][0] == '-') {
       switch (argv[i][1]) {
               case 'I':
		 infile = open_hydro_file("",&argv[i][2],"",PRINT_MSG);
		 if (infile < 0) {
		   exit(1);
		 }
		 iflag = 1;
		 break;

               case 'O':
		 outfile = create_hydro_file(&argv[i][2], OVERWRITE);
		 if (outfile < 0) {
		   fprintf(stderr,"\nUnable to open %s for output\n", &argv[i][2]);
		   exit(1);
		 }
		 fprintf(stderr,"\nOpened %s for output in overwrite mode", &argv[i][2]);
		 oflag = 1;
		 break;
		 
               case 'B':
		 badfile = create_hydro_file(&argv[i][2], OVERWRITE);
		 if (badfile < 0) {
		   fprintf(stderr,"\nUnable to open %s for output\n", &argv[i][2]);
		   exit(1);
		 }
		 fprintf(stderr,"\nOpened %s for output in overwrite mode", &argv[i][2]);
		 bflag = 1;
		 break;
                   
               case 'S':
		 statfile = fopen(&argv[i][2], "r");
		 if (statfile == NULL) {
		   fprintf(stderr,"\nunable to open statfile: %s\n", &argv[i][2]);
		   exit(1);
		 }
		 fprintf(stderr,"\nOpened %s \n", &argv[i][2]);
		 statfilename = &argv[i][2];
		 sflag = 1;
		 break;

               case 'D':
		 deepfile = fopen(&argv[i][2], "r");
		 if (deepfile == NULL) {
		   fprintf(stderr,"\nunable to open deep_statfile: %s\n", &argv[i][2]);
		   exit(1);
		 }
		 deepfilename = &argv[i][2];
		 fprintf(stderr,"\nOpened %s ", &argv[i][2]);
		 dflag = 1;
		 break;
		 
               case 'N':
		 if (sscanf(&argv[i][2],"%f", &nstderr) != 1) {
		   fprintf(stderr,"\nError parsing %s", argv[i]);
		   exit(1);
		 }
		 nflag = 1;
		 break;
                  
               case 'M':
		 if (sscanf(&argv[i][2],"%f/%f", &nstderr_surf, &rhomax) != 2) {
		   fprintf(stderr,"\nError parsing %s", argv[i]);
		   fprintf(stderr,"\nCorrect usage:  -M<nstderr_surf>/<rhomax>\n");
		   exit(1);
		 }
		 rflag = 1;
		 break;
		 
               case 'P':
		 plotflag = 1;
		 strcpy(fname, &argv[i][2]);
		 strncat(fname, "_ts.dat",8);
		 tsfile = fopen(fname,"w");
		 break;

               case 'Q':
		 if (sscanf(&argv[i][2],"%f", &percent_lim) != 1) {
		   fprintf(stderr,"\nError parsing %s", argv[i]);
		   exit(1);
		 }
		 break;
                  
               case 'L':
		 logfile = fopen(&argv[i][2],"w");
		 if (logfile == NULL) {
		   fprintf(stderr,"\nUnable to open logfile for writing: %s\n", &argv[i][2]);
		   exit(1);
		 }
		 fprintf(logfile,"%s  ", &argv[i][2]);
		 break;

               default :
		 print_usage(argv[0]);
		 fprintf(stderr,"\nError parsing command line");
		 fprintf(stderr,"\n in particular: %s\n", argv[i]);
		 exit(1);
       }  /* End switch */
     }  /* End check for command line flag. */
   }  /* End for loop over each command line argument. */

   /* Check for required flags. */
   if (!iflag || !oflag || !sflag || !nflag) {
       print_usage(argv[0]);
       fprintf(stderr,"\n\nYou must specify I, O, S and N options!!!\n");
       exit(1);
   }

   /* Load information from the statfile
    * and store in the coeffs structure
    * (which is also allocated in this
    * routine). */
   nbins = load_coeffs(statfile);

   /* If we have deep coefficients, load them
    * instead.  Will overwrite whatever coeffs
    * and bin information were initialized in
    * load_coeffs because there is only one
    * coeffs structure.. */
   if (dflag) {
       load_deep_coeffs(deepfile);
   }

   /* Initialize counts of stations and scans. */
   totsta = badscan = totscan = 0;
   badsta = 0;

   /* Change percent to fraction, if necessary. */
   if (percent_lim > 1.0)    
     percent_lim /= 100.0;
   
   /* Loop for each station ... */
   while ((status = read_hydro_hdr(infile, &hdr)) == 0) {

     /* Increment the number of stations and scans read. */
     ++totsta;
     totscan += hdr.nobs;

     /* Check the scans of the current station */
     newdata = check_scans(infile, &hdr, &badts, &baddata);

     /* Increment the number of bad scans based
      * on the results of check_scans. */
     badscan += badts;
       
     hdr.qual[3] = '1';

     /* Check that the station is acceptable
      * based on the test that there are non-
      * zero nobs.  If the station is rejected
      * by check_scans, nobs = 0. */
     if ((hdr.nobs = newdata->nobs) > 0)
       write_hydro_station(outfile, &hdr, newdata);
     else
       /* For a bad station, increment the
	* count of rejected stations. */
       ++badsta;

     /* If the user requests saving rejected stations,
      * write the station, if rejected, to outfile. */
     if (bflag) {
       if ((hdr.nobs = baddata->nobs) > 0)
	 write_hydro_station(badfile, &hdr, baddata);
     }
      
     /* Move to the next station. */
     get_separator(infile);

     /* Recycle memory space by clearing all
      * data for each variable.  Also, free
      * all memory associated with the data
      * structures. */
     for (i = 0; i < hdr.nprops; ++i)  {
       free(newdata->observ[hdr.prop_id[i]]);
       free(baddata->observ[hdr.prop_id[i]]);
     }
     free(newdata);
     free(baddata);
     
   }  /* End while loop over each station */ 
 
   report_status(status, stderr);
   close(infile);
   
   /* Print out tallies of eliminated data ... */
   fprintf(logfile, "%6d bad scans of %6d : (%5.1f%%)\n", badscan, totscan, (float)badscan*100./(float)totscan);
   fprintf(logfile, "%6d bad stations of %6d : (%5.1f%%)\n", badsta, totsta, (float)badsta*100./(float)totsta);
   /*   fprintf(logfile, "    : %s : %s\n", statfilename, deepfilename);*/

   fflush(logfile);

   if (plotflag) {
     fclose(tsfile);
   }

   fprintf(stderr,"\n\nEnd of %s. \n", argv[0]);
   exit(0);
}   /* End main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nUsage:  %s -Iinfile -Ooutfile -Sstatfile [-Bbad_file] [-Ddeep_statfile]  [-N#stderr_ts] [-M#stderr_surf/rho_max] [-Qpercent] [-Llogfile] [-Pfile_rootname]", program);

   fprintf(stderr," ");
   fprintf(stderr,"\n   -I  : specifies input HydroBase file");  
   fprintf(stderr,"\n        ex: -I7102_0.raw ");
   fprintf(stderr,"\n   -O  : specifies output file  ");
   fprintf(stderr,"\n        ex: -O7102_0.st1 ");
   fprintf(stderr,"\n   -S  : specifies file with statistics for each sigma bin.");
   fprintf(stderr,"\n        ex: -S7102_0.tsfit ");
   fprintf(stderr,"\n  [-B] : specifies file to which bad data will be written  ");
   fprintf(stderr,"\n        ex: -O7102_0.bad ");
   fprintf(stderr,"\n  [-D] : specifies file in which statistics for sigma bins deeper than ");
   fprintf(stderr,"\n         1000m are defined.");
   fprintf(stderr,"\n        ex: -D7102.statfit ");
   fprintf(stderr,"\n  [-N] : # of standard errors to define good data.");
   fprintf(stderr,"\n        ex: -N2.3 ");
   fprintf(stderr,"\n  [-M] : # of standard errors for surface waters defining");
   fprintf(stderr,"\n         good data and the density boundary (sigma0) to");
   fprintf(stderr,"\n         which this # of std errors is applied.");
   fprintf(stderr,"\n        ex: -N5/26.6 ");
   fprintf(stderr,"\n  [-P] : enables plotfile option and specifies root name");
   fprintf(stderr,"\n         of files to which _ts.dat will be appended.");
   fprintf(stderr,"\n        ex: -P7102_0");
   fprintf(stderr,"\n  [-Q] : quota (in percent) of obs that can be bad");
   fprintf(stderr,"\n         before declaring entire station bad. ");
   fprintf(stderr,"\n         default is %f ", PERCENT_LIM);
   fprintf(stderr,"\n\n");  
   return;
} /* end print_usage() */


/****************************************************************************/
/* Given the statfile, load information about the
 * density bins and statistics necessary for the
 * rest of the program.  The statfile is of the
 * format:
 *
 * number_of_bins
 * bin# sigmin sigmax <S|T|X> slope intercept error #pts
 * bin# sigmin sigmax <S|T|X> slope intercept error #pts
 * ...
 * */

int load_coeffs(FILE *file)
{
  /* n     = number of bins
   * i     = bin counter in loop
   * index = temp storage of bin number from file */
  int n, i, index;

  /* Read the statfile header; a single value
   * specifying the number of density bins. */
  fscanf(file, "%d", &n);

  /* Allocate space for the profile of
   * min/max bin densities + coeffs */
  sigmin = (float *) malloc(n * sizeof(float));
  sigmax = (float *) malloc(n * sizeof(float));
  coeff = (struct coefficients *) malloc(n * sizeof(struct coefficients));

  /* Read information from each bin. */
  for (i = 0; i < n; ++i) {

    fscanf(file,"%d%f%f %c", &index, &sigmin[i], &sigmax[i], &coeff[i].xvar);
    fscanf(file,"%f%f%f%d", &coeff[i].m1, &coeff[i].b1, &coeff[i].e1, &coeff[i].n1);

    /* Check that the bin number in the file
     * matches the current count of bins. */
    if (index != i) {
      fprintf(stderr,"\nMismatch of index in load_coeffs().\n");
      exit(1);
    }
  } /* End looping over all bins in statfile */

  close(file);
  return(n);

} /* End load_coeffs() */

/****************************************************************************/
int load_deep_coeffs(FILE *file)
/* Reads from a statfile and stores the values corresponding to
 * sigma values > 30 (referenced to 2000 db or 4000 db) in the array,
 * "coeff".  A little bit of checking is done to determine if the
 * sigma bins match up with the ones established by load_coeffs().
 * Any sigma values > 30 in coeff are overwritten with the values
 * provided here.
 */
{
   int n, i, index, junk;
   float s1, s2;

   /* Read the header line - a single value
    * representing the number of density bins. */
   fscanf(file, "%d", &n);

   /* Test that the number of bins in the deepfile
    * matches the number of bins in the statfile. */
   if (n != nbins) {
      fprintf(stderr,"Deep statfile has different number of sigma bins.\n");
   }

   /* Loop over each density bin, */
   for (i = 0; i < n; ++i) {

     /* Read the bin number and min/max sigma bounds. */
     fscanf(file,"%d%f%f", &index, &s1, &s2);

     /* Check that this sigma value is a deep value. */
     if (s1 > 30.) {
       /* Check that sigmin/sigmax are the same as
	* in the statfile. */
       if ((s1 < (sigmin[index]-.0001)) || (s1 > (sigmin[index]+.0001)) ) {
	 fprintf(stderr,"\nFatal Error: ");
	 fprintf(stderr," Deep sigma bins do not match regular bins.\n");
	 exit(1);
       }
       
       /* If the bins match up, then load new values
	* for these deep coefficients. */
       fscanf(file," %c%f%f%f%d", &coeff[index].xvar, &coeff[index].m1, 
	      &coeff[index].b1, &coeff[index].e1, &coeff[index].n1);
     }
     else {
       /* Skip over this line because this bin
	* is not dense enough. */
       fscanf(file," %*c%*f%*f%*f%*d", &junk);
     } /* End of high density check. */
   } /* End for loop over all bins */
   
   return (n);
   
} /* End load_deep_coeffs() */

/****************************************************************************/

struct HYDRO_DATA *check_scans(int fd, struct HYDRO_HDR *hptr, int *badts_ptr, struct HYDRO_DATA **badaddr)
/*     
         fd   file descriptor for already open HydroBase file 
       hptr   ptr to HydroBase2 header info
  badts_ptr   returns # of eliminated scans here 
    badaddr   address to return buffer of bad scans
*/
   /* Returns a pointer to a struct HYDRO_DATA containing
    * all the values which pass the statistical criteria.
    * If the entire station is eliminated, the field .nobs
    * will have the value 0.
    *
    * Also returns the bad scans in a struct of format
    * HYDRO_DATA called badaddr.  Memory is allocated here
    * for both the bad and good scans and must be freed in
    * the main program. */
{
  /* Number of standard deviations defining
   * acceptable envelope around T-S curve. */
  float nstde;

  /* Reference pressure/depth, temp
   * wrt pref, surface pressure. */
  double  pref, tref, p0 = 0.0;

  /* Temp storage of single (scan values)
   * of pressure, depth, temp, potential
   * temp, salt, potential density, x,
   * and y. */
  double  p, d, t, th, s, sig, x, y;

  /* Array (one entry for each variable
   * in the station) stores one line from
   * each scan. */
  double *scan;

  int bad;
  int n, nb, i, j, k, bin, found;

  /* Data structures, all identical in size,
   * initialized here to hold the good, the
   * bad, and the buffered data. */
  struct HYDRO_DATA *dataptr, *badptr, *stabuffer;
  
  /* Initialize counters... */
  bad = 0;     
  n = nb = 0;
  
  /* Allocate memory for scan info - array length of nprops. */
  scan = (double *) calloc((size_t)hptr->nprops, sizeof(double));
  
  /* Allocate memory for good, bad, and buffer data storage. */
  dataptr = (struct HYDRO_DATA *) calloc((size_t)1, sizeof(struct HYDRO_DATA));
  badptr = (struct HYDRO_DATA *) calloc((size_t)1, sizeof(struct HYDRO_DATA));
  stabuffer = (struct HYDRO_DATA *) calloc((size_t)1, sizeof(struct HYDRO_DATA));
  
  /* For each variable in the station, */
  for (i = 0; i < hptr->nprops; ++i) {
    /* Allocate memory for the actual storage of the data. */
    dataptr->observ[hptr->prop_id[i]] = (double *) calloc((size_t)hptr->nobs, sizeof(double));
    badptr->observ[hptr->prop_id[i]] = (double *) calloc((size_t)hptr->nobs, sizeof(double));
    stabuffer->observ[hptr->prop_id[i]] = (double *) calloc((size_t)hptr->nobs, sizeof(double));
  } /* End of loop over each station variable */
  
  /* Initialize the numbers of levels/obs and
   * variables/props in each of the possible
   * data storage structures.  It is important
   * that the nobs start explicitly at zero
   * because these values will be used as checks
   * for the existance of good/bad obs back
   * in the main program. */
  dataptr->nobs = 0;
  dataptr->nprops = hptr->nprops;
  badptr->nobs = 0;
  badptr->nprops = hptr->nprops;
  stabuffer->nprops = hptr->nprops;
  stabuffer->nobs = hptr->nobs;

  /* If the user requests a plot file,
   * allocate arrays for storing plot info. */
  if (plotflag) {
    plot_s = (float *) calloc((size_t)hptr->nobs, sizeof(float));
    plot_t = (float *) calloc((size_t)hptr->nobs, sizeof(float));
    plot_sg = (float *) calloc((size_t)hptr->nobs, sizeof(float));
    plot_p = (float *) calloc((size_t)hptr->nobs, sizeof(float));
    sym = (int *) calloc((size_t)hptr->nobs, sizeof(int));
    tspen = (int *) calloc((size_t)hptr->nobs, sizeof(int));
  }

  /* Loop for each observation depth, */
  for (j = 0; j < hptr->nobs; ++j) {
    
    /* Read each scan into the variable scan. */
    if (get_data_scan(fd, scan, hptr->nprops, hptr->prop_id) > 0) {
      fprintf(stderr, "ERROR reading data in get_data_scan()!");
      return(dataptr);
    }
       
    /* Buffer each scan - copy the read scan
     * into the buffer data storage array. A
     * complete copy of the station grows
     * as each scan is read in. */
    for (k = 0; k < hptr->nprops; ++k) {
      stabuffer->observ[hptr->prop_id[k]][j] = scan[k];
    }

    /* Find the prop_id (location in scan)
     * of p, d, t, and s.  Start by lifting
     * flags at each variable. */
    p = d = t = s = FLAG;

    /* Next, cycle through all variables in
     * the station/scan and assign values
     * from the scan to pressure, depth,
     * temp, and salinity. */
    for (i = 0; i < hptr->nprops; ++i) {
      switch ((enum property) hptr->prop_id[i]) {
      case PR:
	p = scan[i];
	break;

      case DE:
	d = scan[i];
	break;

      case TE:
	t = scan[i];
	break;

      case SA:
	s = scan[i];
	break;

      default:
	;
      }  /* End switch */
    } /* End of finding values in scan */
    
    /* Compute potential temperature (ref to surface). */
    th = hb_theta(s, t, p, p0);

    /* Determine ref pressure by assigning any
     * values above 1000m to zero ref. p., all
     * values below 3000m to 4000m reference
     * level and all values in between to a
     * 2000m reference level.  Same method
     * for determining pref (and sigma, later)
     * as in hb_statfit_ts.c */
    pref = 2000.0;
    if (p <= 1000.0)
      pref = 0.0;
    if (p > 3000.0)
      pref = 4000.0;

    /* Compute the potential temperature
     * referenced to the reference depth. */
    tref = hb_theta(s, t, p, pref);
    
    /* Compute the potential density with
     * respect to the reference depth. */
    hb_svan(s, tref, pref, &sig);
    
    /* If the user requests plot information,
     * write this to the plot storage vars. */
    if (plotflag) {
      plot_sg[j] = sig;
      plot_t[j] = th;
      plot_s[j] = s;
      plot_p[j] = p;
      sym[j] = (int) pref/1000;
      tspen[j] = BLACK;
    }
          
    /* Store the number of acceptable
     * standard dev. from the mean T-S
     * curve. */
    nstde = nstderr;
 
    /* Check if this scan is in
     * the "surface" waters and
     * use a special surface value
     * for the standard deviation. */
    if (rflag) {
      if (sig < rhomax)
	nstde = nstderr_surf;
    }

    /* Search for density bin for this scan. */
    found = 0;  /* Found flag */
    bin = -1;   /* Initial value for bin index */

    /* Increment the bin index; first loop
     * starts a 0.  Keep checking successive bins
     * until we find the correct bin.  However,
     * in order to end the loop, we need to go
     * through the ++bin operation one more time,
     * hence the --bin operation at the very end
     * in order to store the correct index of the
     * bin where the found flag was set to 1. */
    while ((++bin < nbins) && !found) {
      if ((sig >= sigmin[bin]) && (sig < sigmax[bin]))
	/* Scan is within a density bin */
	found = 1;
    }
    --bin;
    
    /* Check whether the scan is out
     * of statistical range. */
    if (!found) {
      /* For the case when the scan is not in
       * sigma bins, augment number of bad scans. */
      ++bad;
      
      /* Copy over all variables in the scan to
       * the bad data. */
      for (i = 0; i < hptr->nprops; ++i) {
	badptr->observ[hptr->prop_id[i]][nb] = scan[i];
      }
      ++nb;
      if (plotflag) {
	tspen[j] = RED;
      }
      
    }
    else if (coeff[bin].n1 < 3) {
      /* For the case when there are not enough
       * pts in the data in each sigma bin (as
       * determined by hb_statfit_ts.c), augment
       * the number of bad scans. */
      ++bad;

      /* Copy over all variables in the scan to
       * the bad data. */
      for (i = 0; i < hptr->nprops; ++i) {
	badptr->observ[hptr->prop_id[i]][nb] = scan[i];
      }
      ++nb;
      if (plotflag) {
	tspen[j] = RED;
      }
    }
    else  {
      /* This could be a good scans, check T-S */
      x = th;
      y = s;

      /* Independent variable check, as determined
       * in hb_statfit_ts.c */
      if (coeff[bin].xvar == 'S') {
	x = s;
	y = th;
      }
      
      /* Check if the T-S point lies within the
       * envelope of the mean T-S curve by
       * computing the difference between the
       * regression's expected value and the actual
       * value to the number of allowable standard
       * deviations from the mean curve. */
      if (ABS(y - (x * coeff[bin].m1 + coeff[bin].b1)) > (nstde *coeff[bin].e1)) {
	/* Outside the envelope => bad scan.
	 * Augment the bad scan counter. */
	++bad;

	/* Copy over all variables in the scan to
	 * the bad data. */
	for (i = 0; i < hptr->nprops; ++i) {
	  badptr->observ[hptr->prop_id[i]][nb] = scan[i];
	}
	++nb;
	if (plotflag) {
	  tspen[j] = RED;
	}
      }
      else {
	/* Inside the envelope => good scan.
	 * Copy good scan to good data structure
	 * and augment the number of good points. */
	for (i = 0; i < hptr->nprops; ++i) {
	  dataptr->observ[hptr->prop_id[i]][n] = scan[i];
	}
	++n;
      } /* End of good scan check */
    } /* End of bin suitability check (#pts, within a density bin) */
  } /* End for loop over each scan */
  
  
  /* Check the tallies.  If # of bad
   * pts > 65%, eliminate the entire station. */
  if (bad >  hptr->nobs * percent_lim) {
    /* Update bad station nobs. */
    *badts_ptr = hptr->nobs;
    
    /* Set output data nobs to zero so
     * main program knows that the station
     * is rejected. */
    dataptr->nobs = 0;
    
    free(badptr);
    
    /* Replace all the bad scans with
     * the buffered station */
    badptr = stabuffer;
  }
  else {
    /* We keep the station but update the
     * number of good (bad) scans to the
     * good (bad) data structures. */
    dataptr->nobs = n;
    badptr->nobs = nb;

    /* Free the buffer if station has good scans */
    free(stabuffer);       
    *badts_ptr = bad;
  } /* End of reject station check */ 
  
    /* If the user requests plotting information, */
  if (plotflag) {
    
    if (bad > hptr->nobs * percent_lim) {
      for (i = 0; i < hptr->nobs; ++i) {
	tspen[i] = RED;
      }
    }
    for (i = 0; i < hptr->nobs; ++i) {

      /* Write to GMT xy scatter plot file for T-S
       * diagram with colors based on depth and
       * a plus symbol if rejected or dot if accepted. */

      if(tspen[i] == RED) {
	/* Outcast scan */
	fprintf(tsfile, " %6.3f %6.3f %6.3f +\n", plot_s[i], plot_t[i], 
	      plot_p[i]);
      }
      else{
	/* Accepted scan */
	fprintf(tsfile, " %6.3f %6.3f %6.3f c\n", plot_s[i], plot_t[i], 
	      plot_p[i]);
      }
    }
    free(plot_s);
    free(plot_t);
    free(plot_sg);
    free(plot_p);
    free(tspen);
    free(sym);
  }  /* end if plotflag */
  
  free(scan);
  *badaddr = badptr;
  return(dataptr);
  
} /* end check_scans() */

/****************************************************************************/
