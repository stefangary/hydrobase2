/*  hb_statfit_ts.c
................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth G Curry
                             Woods Hole Oceanographic Instition
                             updated 2000 to ANSI standards

			     Commenting by Stefan Gary
................................................................................
-------------------------------------------------------------------------------
.  USAGE:  hb_statfit  list_of_files -Ooutfile -Ssigbin_file [-Ddir] [-Eextent] [-Mmaxslope] [-Xminpts] [-Pplotfiles_root_name]
_______________________________________________________________________________
.   Computes linear fits to the t,s data for each sigma bin specified  
.   in the sigbin_file.  The theta-salt curve is approximated by connecting
.   the mean t,s point in each density bin.  The line tangent to this curve
.   at the mean t,s point then approximates the t-s relation in each bin. 
.   We determine the slope, y-intercept, and standard error for each bin.
.   Initially, we assume theta as the independent variable. If the slope of
.   the linear fit exceeds a critical value (maxslope), we designate salt
.   as the independent variable and compute the statistics for theta = F(salt)
.   and oxygen = F(salt).  This helps to identify situations where the 
.   F(theta) fit for a particular density bin does not do a good job of
.   specifying a particular density bin, i.e. F(theta) does not approximate the 
.   overall t-s relationship, F(salt) is better.
.
.     
.   The output file contains 1 header line consisting of the number of bins
.   followed by a line for each density bin (X => empty bin, T|S signify the
.   independent variable).  Values are whitespace delimited.
.
.   number_of_bins
.   bin# sigmin sigmax <T|S|X> slope intercept error n_data_pts_in_bin
.   bin# sigmin sigmax <T|S|X> slope intercept error n_data_pts_in_bin
_______________________________________________________________________________
*/
#include <stdio.h>
#include <stdlib.h>   /*sfg added to deal with exit and malloc warnings*/
#include <string.h>   /*sfg added to deal with strncat warnings*/
#include <math.h>
#include "hb_sigbins.h"
#include "hydrobase.h"

#define    PRINT_MSG    1
#define    MAXSLOPE     1.0  /* (delta-Salt/delta-Theta) value to determine
                                which should be the x-variable in the fit */
#define    MINPTS     3     /* min # of points desired in a density bin */
                                 
#define    EXTENT    ""
#define    DIR       ""

/*******  Global variables *******/

/* Each binrec struct contains (defined
 * in hb_sigbins.h):
 * nxy
 * xsum
 * ysum
 * xysum
 * xxsum
 * yysum
 *
 * Pointers to binrec allow us to access
 * a profile of binrec structures. */
struct binrec  *FofT_bins, *FofS_bins;

/* Min/max pot. dens. of each bin - vertical profiles */
float  *sigmin, *sigmax;

/* Reference level for each bin - vertical profile */
int  *reflev;

/* Number of bins in bin profiles */
int  nbins;

double  maxslope;

FILE *slopefile, *stddevfile, *meanfile;

int  plotflag;  

/*******************************/
/* Prototypes for locally defined functions */

/* Gives user info in case of mistake. */
void   print_usage(char *);

/* Reads density bins from file
 * and allocate memory for bin
 * profiles/structures. */
int    load_bins(char *);

/* For the given binrec profile
 * with n points, set all values
 * of each structure in the
 * profile to zero. */
void   initialize(struct binrec *, int);

/* Sort each scan of a profile into its
 * density bin and compute running intermediate
 * values for the linear regression. */
void   insertdata(struct HYDRO_HDR *,struct HYDRO_DATA *);

/* Compute the linear regression at each level (bin). */
double xyslope(struct binrec *, int, int, int *);

/**/
void   plot_slope(int, double, double, double, double);

main (int argc, char **argv)
{
   int     n, status, t_is_x, check_x; 
   int     i, nfiles = 0, curfile = 1; 
   int     infile, iflag, oflag, sflag; 
   int     minpts;
   FILE   *outfile; 
   char   *binfilename, fname[80];
   char   *dir, *extent;
   double  m, b, E, x, ybar, xbar;
   struct HYDRO_HDR hdr;
   struct HYDRO_DATA data;
   struct binrec *bin;

   /* Are there command line arguments? */
   if (argc < 4) {
      print_usage(argv[0]);
      exit(1);
   }

   /* Set default values... */
   oflag = sflag = 0;
   plotflag = 0;
   minpts = MINPTS;
   maxslope = MAXSLOPE;
   dir = DIR;
   extent = EXTENT;

   /* Initialize some variables */
   for (i = 0; i < MAXPROP; ++i) {
     data.observ[i] = (double *) NULL;
   }
   hdr.prop_id = (int *)NULL;

   /* Parse command line arguments... */
   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'D':
                        dir = &argv[i][2];
                        break;
               case 'E':
                        extent = &argv[i][2];
                        break;
               case 'O':
                   outfile = fopen(&argv[i][2],"w");
                   if (outfile == NULL) {
                      fprintf(stderr,"\nUnable to open %s for output. \n", &argv[i][2]);
                      exit(1);
                   }
                   oflag = 1;
                   break;

               case 'S':
                  binfilename = &argv[i][2];
                  sflag = 1;
                  break;

               case 'M':
                  if (sscanf(&argv[i][2],"%lf", &maxslope) != 1) {
                    fprintf(stderr,"\nError parsing %s", argv[i]);
                    exit(1);
                  }
                  break;

               case 'X':
                  if (sscanf(&argv[i][2],"%d", &minpts) != 1) {
                    fprintf(stderr,"\nError parsing %s", argv[i]);
                    exit(1);
                  }
                  break;

               case 'P':
                  plotflag = 1;
                  strcpy(fname,&argv[i][2]);
                  strncat(fname,".mean.ts",9);
                  meanfile = fopen(fname,"w");
                  if (meanfile == NULL) {
                     fprintf(stderr,"\nUnable to open %s\n", fname);
                     exit(1);
                  }

                  strcpy(fname,&argv[i][2]);
                  strncat(fname,".slope.ts",10);
                  slopefile = fopen(fname,"w");
                  if (slopefile == NULL) {
                     fprintf(stderr,"\nUnable to open %s\n", fname);
                     exit(1);
                  }
                   strcpy(fname,&argv[i][2]);
                  strncat(fname,".stddev.ts",10);
                  stddevfile = fopen(fname,"w");
                  if (stddevfile == NULL) {
                     fprintf(stderr,"\nUnable to open %s\n", fname);
                     exit(1);
                  }
                 
                  break;

               case 'h':
                 print_usage(argv[0]);
                 exit(0);
                 
               default:
                 fprintf(stderr,"\nError parsing command line");
                 fprintf(stderr,"\n in particular: %s\n", argv[i]);
                 exit(1);
            }  /* end switch */
       }  /* end if */
       else  {
           ++nfiles;
       }

   }  /* End looping over command line arguments */

   /* Check that command line contains sufficient
    * information to proceed. */
   if (!nfiles || !oflag || !sflag) {
       print_usage(argv[0]);
       fprintf(stderr,"\n\nYou must specify input, output, and sigma bin files!!!\n");
       exit(1);
   }

   /* Define sigma bins.  Reference level
    * (pressure) is assigned based on the
    * values of the bin edges. */
   if ((nbins = load_bins(binfilename)) < 0) {
       exit(1);
   }

   /* Allocate space for the structures that
    * store statistical properties of each bin. */
   initialize(FofT_bins, nbins);
   initialize(FofS_bins, nbins);

   /* Log the number of bins in header of output file. */
   fprintf(outfile,"%d", nbins);

   /* Loop for each file */
   do {
      infile = open_hydro_file(dir, argv[curfile], extent, PRINT_MSG);
      if (infile  < 0) 
	goto NEXTFILE;
       
      /* Loop for each station */
      while ((status = get_station(infile, &hdr, &data)) == 0) {  

	/* Calculate (running values) of
	 * intermediate values for the linear
	 * regression, at all depths. */
	insertdata(&hdr, &data);

      }  /* End while loop over each station */ 
 
      report_status(status, stderr);
      close(infile);
   NEXTFILE:
      ;
   } while (curfile++ < nfiles); /* End loop over each file */

   /*  Compute slope, intercept and standard deviation
    *  and write to outfile.  Check each bin for
    *  existing pts; there must be at least 3 pts to
    *  compute the standard error ... */
   fprintf(stderr,"\n computing slope and variance ...\n");

   /* Loop over each bin */
   for (i = 0; i < nbins; ++i) {

     /* Write bin number and density bin
      * bounds to output file. */
     fprintf(outfile,"\n%4d %6.3f %6.3f", i, sigmin[i], sigmax[i]);

     /* Check that we have at least three
      * points in the current bin. */
     if ((n = FofT_bins[i].nxy) >= 3) {

       /* We require that xyslope check
	* whether or not to use t or s as
	* the independent variable if
	* check_x = 1.  For check_x =0, we
	* require t always to be the
	* independent var, regardless of
	* the slope. */
       check_x = 1;

       /* Compute the slope:
	* FofT_bins = bins to compute slope on
	* i         = level in bins to compute slope
	* check_x   = check whether t is independent variable
	* t_is_x    = returned exit flag:
	*                1 => t is indep. var., m = dt/ds
	*                0 => s is indep. var., m = ds/dt
	*               -1 => not enough data,  m = 0 */
       m = xyslope(FofT_bins, i, check_x, &t_is_x);

       if (t_is_x >= 0) {
	 /* We have a sucessful calculation of the
	  * slope (regardless t or s as indep. var.) */
	 
	 /* For t_is_x = 1, use FofT_bins
	  * For t_is_x = 0, use FofS_bins */
	 bin = t_is_x ? FofT_bins : FofS_bins;

	 /* Compute the mean value of x
	  * and y in the bin */
	 ybar = bin[i].ysum / n;
	 xbar = bin[i].xsum / n;

	 /* Compute the y-intercept.  Note that the
	  * actual line lies only between the centroids
	  * of the two adjacent bins and does not
	  * necessarily go though the centroid of the
	  * current bin.  However, since we're computing
	  * and estimate of the slope of the line at
	  * this bin, we are justified in determining
	  * the line with the slope m that must go
	  * through mean value (and hence pinning down
	  * intercept b) but this is NOT the same
	  * as the line computed in xyslope. */
	 b = ybar - m * xbar;

	 /* Compute an estimate for the variance
	  * by computing the estimated standard error
	  * of the regession model.  According to
	  * McClave and Dietrich (intro stats book,
	  * p. 682),
	  *
	  * s^2 = SSE/(n-2)
	  *
	  * where SSE  = SSyy - mSSxy
	  *
	  *       SSxx = xxsum - n*xbar*xbar
	  *
	  *       SSyy = yysum - n*ybar*ybar
	  *
	  *       SSxy = xysum - n*xbar*ybar
	  *
	  * Also, they say that for a typical linear
	  * regression,
	  *
	  *          m = SSxy/SSxx
	  *
	  * which in our case means that SSxy/SSxx is
	  * an estimate for m because we really have
	  * an m computed from different density bins
	  * and not directly from the data in the
	  * central bin.
	  *
	  * The upshot of all this is that:
	  *
	  * -m*SSxy = -m*SSxy + (SSxy)^2/SSxx - (SSxy)^2/SSxx
	  *
	  *         = -m*SSxy + SSxx*(SSxy/SSxx)^2 -m*SSxy
	  * 
	  *         = -2*m*SSxy + (m^2)*SSxx
	  *
	  * which means that:
	  *
	  * s^2 = (SSyy - 2*m*SSxy + (m^2)*SSxx)/(n-2),
	  *
	  * which is the same as the expression for E,
	  * below, assuming, of course, that the slope
	  * computed in xyslope is a good estimate for
	  * the actual regression slope.
	  *
	  * In hb_statchk_ts.c, we use this error, E, to
	  * build a confidence interval between a value
	  * on the regression line and an actual data
	  * point. E by itself is NOT strictly the correct
	  * estimate for the standard deviation of the
	  * prediction error.  Rather, the correct value
	  * is (M&D, p. 707):
	  *
	  * s_(y-yhat) = s*sqrt(1 + 1/n + (x_p - xbar)^2/SSxx)
	  *
	  *
	  * The following argument tries to place an upper
	  * bound on the sqrt() term.
	  *
	  * (1) For the case when (x_p - xbar) ~ xbar
	  *     (worse case scenario than when x_p ~ xbar,
	  *     valid because density bins are pretty well
	  *     coalesced), the last term is:
	  *
	  *     xbar^2/SSxx = xbar^2/(xxsum - n*xbar^2)
	  *
	  *                 = 1/(xxsum/xbar^2 - n)
	  *
	  *                 = 1/(xxsum*n^2/xsum^2 - n)
	  *
	  *                 = 1/(n[xxsum*n/xsum^2 - 1])
	  *
	  *     Since xxsum > xsum^2 for any data set, we
	  *     have, roughly.
	  *
	  *     xbar^2/SSxx < 1/(n^2)
	  *
	  * (2) The square root terms are subject to the
	  *     upper limit:
	  *
	  *     sqrt(1 + 1/n + 1/n^2)
	  *
	  *     For the case with the minimum number of
	  *     observations (thus making the sqrt term
	  *     the largest it can be), we have n = 3.
	  *
	  *     sqrt(1 + 1/3 + 1/9) = 1.201850
	  *
	  *     This means that if we set 2 standard
	  *     deviations, we should really expect
	  *     an envelope as wide as ~2.4*stddev.  In
	  *     Lozier et al., 1995, they choose a 2.3
	  *     std. deviations about the mean T-S curve
	  *     which would be the upper bound of a
	  *     2 standard deviation interval - 98%
	  *     confidence versus 95% confidence.
	  */
	 E = ABS((bin[i].yysum - n*ybar*ybar - 2*m*(bin[i].xysum  
						    - n*ybar*xbar) + m*m*(bin[i].xxsum - n*xbar*xbar))/(n-2));

	 /* Compute the standard deviation from the variance. */
	 E = sqrt(E);
	 
	 /* Print the independent variable to output. */
	 if (t_is_x)
	   fprintf(outfile," T");
	 else 
	   fprintf(outfile," S");
	 
	 /* Print the slope, intercept, std., and
          * number of points in bin to output file. */
	 fprintf(outfile," %10.5lf %7.3lf %8.4lf %5d",  m, b, E, n);
	 
	 /* Print plotting information to plot file. */
	 if (plotflag) {
	   plot_slope( t_is_x, bin[i].xsum/bin[i].nxy, m, b, E);
	 } /* End of plotting check */
       }
       else {
	 /* Unsucessful exit from xyslope, line added by sfg. */
	 fprintf(outfile," X  0 0 0 %5d",n);
       }/* End xyslope calculation check.*/
     }
     else {
       /* We do not have sufficient numbers of points
	* in this density bin. */
       fprintf(outfile," X  0 0 0 0");
     } /* End bin number of points check */
   }  /* End for loop over each bin */

   exit(0);

}  /* End of main() */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nUsage:  %s list_of_filenames -Ooutfile -Ssigma_binfile [-Ddirname] [-Eextent] [-Mmaxslope] [-Xminpts] [-Pplotfile]", program);

   fprintf(stderr," ");
   fprintf(stderr,"\n   -O  : specifies output file  ");
   fprintf(stderr,"\n          ex: -O7102.ts_fit ");
   fprintf(stderr,"\n   -S  : specifies file in which sigma bins are defined.");
   fprintf(stderr,"\n          ex: -Ssigbins.natl ");
   fprintf(stderr,"\n  [-D] : specifies directory of input files (default is ./) ");
   fprintf(stderr,"\n        ex: -D../data/ ");
   fprintf(stderr,"\n  [-E]  : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n  [-M] : maximum slope (delta-s/delta-t) for which theta");
   fprintf(stderr,"\n         will be used as x variable in linear fit");
   fprintf(stderr,"\n  [-X] : min # of points in density bin to keep data");
   fprintf(stderr,"\n         (default = 3)");
   fprintf(stderr,"\n          ex: -X3 ");
   fprintf(stderr,"\n  [-P] : turns on plotfile option and specifies fileroot ");
   fprintf(stderr,"\n          to which will be appended mean.ts, stddev.ts and slope.ts");
   fprintf(stderr,"\n          ex: -P7207");
   fprintf(stderr,"\n\n");  
   return;
} /* end print_usage() */

/****************************************************************************/
/* Reads the binfile and loads the bin definitions
 * into the global variables.  Returns the number
 * of bins or -1 if an error occurs.  In case of
 * an error, an appropriate message is written to
 * the stderr device.  The binfile is of the
 * format (one line only):
 * sigmin_1 sigmax_1 sigmin_2 sigmax_2 ... sigmin_n sigmax_n
 */
int load_bins(char *name) {

   /* n = store number of bins
    * i = bin loop counter */
   int n, i;

   float r;

   FILE *binfile;

   /* Open the binfile for reading. */
   binfile = fopen(name,"r");

   /* Check for binfile existance. */
   if (binfile == NULL) {
       fprintf(stderr,"\nUnable to open Sigma Bin File: %s \n", name);
       return(-1);
   }

   /* Determine number of bins defined in binfile. */
   n = 0;
   while (fscanf(binfile,"%*f %f", &r) == 1)
     /* Note that the first real read is skipped
      * due to the * and only the second real is
      * assigned to &r.*/
        ++n;
   close(binfile);

   /* Allocate profile memory. */
   sigmin = (float *) malloc(n * sizeof(float));
   sigmax = (float *) malloc(n * sizeof(float));
   reflev = (int *) malloc(n * sizeof(int));
   FofT_bins = (struct binrec *) malloc(n * sizeof(struct binrec));
   FofS_bins = (struct binrec *) malloc(n * sizeof(struct binrec));

   /* Reopen binfile and read in density bins.
    * NOTE: The following code assigns ref.
    * pressures according to the value of the
    * given sigma according to:
    * sigma < 30      => ref. p. = 0
    * sigma = [30,40] => ref. p. = 2000
    * sigma > 40      => ref. p. = 4000
    *
    * This is a valid assumption since in the
    * T-S diagrams over a reasonable range of
    * values T = [0,30], S = [34,38], we see
    * almost all the values at
    * P = 0    => sigma = 20s
    * P = 2000 => sigma = 30s
    * P = 4000 => sigma = 40s
    *
    * and deviations from this trend occur in
    * regions of the T-S diagram that don't
    * apply to the real ocean at that depth.
    * For example, at reference pressure of
    * 4000db, we can't have 20oC water!
    */
   binfile = fopen(name,"r");
   i = 0;
   while (fscanf(binfile,"%f %f", &sigmin[i], &sigmax[i]) == 2) {

     /* Default reference level is sigma_2. */
     reflev[i] = 2000;

     /* Modify reference level depending on MIN density. */
     if (sigmin[i] < 30)
       reflev[i] = 0;
     if (sigmin[i] > 40)
       reflev[i] = 4000;

     /* Augment counter to go to next bin. */
     ++i;
   } /* End of loop through all bins. */

   close(binfile);

   /* Check that the number of bins read in is the same
    * as the number of bins counted. */
   if (n != i) {
      fprintf(stderr, "\nNumber of bins defined: %d differs from what was allocated for: %d", i, n);
      return(-1);
   }

   /* If all goes well, return the number of bins counted,
    * which is the same as the number of bins read in. */
   return (n);

}  /* End load_bins() */

/****************************************************************************/

void initialize(struct binrec *b, int n)

/* For the given binrec profile (b is a pointer
 * start of array of binrecs) with n levels
 * (n = number of binrecs), set all values
 * of each structure in the profile to zero. */
{
  int i; /* Level counter in binrec profile */
  
  for ( i = 0; i < n; ++i) {
    b[i].nxy = 0;
    b[i].xsum = 0.0;
    b[i].ysum = 0.0;
    b[i].xysum = 0.0;
    b[i].xxsum = 0.0;
    b[i].yysum = 0.0;
  }
  return;
} /* end initialize() */

/****************************************************************************/
/* For a data in a given profile (in
 * hydrobase struct format), determine
 * which density bin the data reside
 * and then augment running calculations
 * of intermediate calculations of
 * regressions accordingly. */
void insertdata(struct HYDRO_HDR *hptr, struct HYDRO_DATA *dptr) {

  /* Temporary single values of pot. temp (th),
   * referenced pot. temp (tref), and
   * pot. density (sigma).  p0 is always
   * the zero surface reference pressure. */
  double    th, pref, tref, sigma, p0 = 0.0;

  /* Temporary single values of temp (t),
   * salt (s), pressure (p), and depth(d)
   * at a single scan. */
  double    t, s, p, d;

  register  i, j;
  char      found;
  
  /* Check that we have pressure, temp,
   * and salinity in this station.  If
   * not, skip the station. */
  if (!available((int)PR, hptr) || 
      !available((int)TE, hptr) || 
      !available((int)SA, hptr) )
    return;
  
  /* Loop for each observation level. */
  for (j = 0; j < hptr->nobs; ++j) {
    
    t = dptr->observ[(int)TE][j];
    s = dptr->observ[(int)SA][j];
    p = dptr->observ[(int)PR][j];
    d = dptr->observ[(int)DE][j];
    
    /* Test if sigma can be computed.
     * For weird values of t and s,
     * skip to the next level with
     * the continue command. Originally,
     * the lower bound on s was set to
     * -8, but sfg changed this to 0.0.
     * This will also also statfit to
     * ignore any HB_MISSING values
     * (set to -9). */
    if ((t < -8) || (s < 0.0)) {
      continue;
    }
     
    /* Determine ref pressure:
     * Depths less than 1km = 0 dbar
     * Depths 1km to 3km = 2000 dbar
     * Depths greater 3km = 4000 dbar */
    pref = 2000.;
    if (d <= 1000.)
      pref = 0.;
    if (d > 3000.)
      pref = 4000.;
    
    /* Compute potential temperature (wrt surface) */
    th = hb_theta(s, t, p, p0);
    
    /* Compute potential temp wrt reference depth */
    tref = hb_theta(s, t, p, pref);
    
    /* Compute sigma, referenced to pref in both t and p */
    hb_svan(s, tref, pref, &sigma);   
     
    /* Search for bin that holds sigma */
    i =  0;  
    found = 0;
    do {
      if ((sigma >= sigmin[i]) && (sigma < sigmax[i]))
	found = 1;
    } while ( (++i < nbins) &&  !found );
    --i; 
    
    if (!found) {
      /* The current density value does not
       * fit into one of the defined bins -
       * register this and skip the value. */
      fprintf(stderr,
	      "\nSigma out of range: %9.3lf :  d=%5.0lf t=%6.2lf s=%6.2lf", 
	      sigma, d, t, s);
    }
    else {
      /* The current density value fits in
       * one of the defined bins, add its
       * values to the running accumulation
       * of online statistics. */

      /* T is independent variable for FofT
       * so x = th and y = s.  T is the
       * dependent variable for FofS so
       * x = s and y = th. */

      FofT_bins[i].xsum += th;
      FofS_bins[i].ysum += th;

      FofT_bins[i].xxsum += th * th;
      FofS_bins[i].yysum += th * th;

      ++FofT_bins[i].nxy;
      ++FofS_bins[i].nxy;

      FofT_bins[i].ysum += s;
      FofS_bins[i].xsum += s;

      FofT_bins[i].yysum += s * s;
      FofS_bins[i].xxsum += s * s;

      FofT_bins[i].xysum += th * s;
      FofS_bins[i].xysum += th * s;

    } /* End of in density bin check */ 
    
  }  /* End for loop over all levels in profile x*/
  
  return;
  
}  /* End insertdata() */

/****************************************************************************/

double xyslope(struct binrec *bin, int j, int determine_indep_var, int *flag_ptr)
/*  
 *                 bin   pointer to array of binrecs
 *                   j   index to bin array being considered  
 * determine_indep_var   0 or 1 
 *            flag_ptr   returned value to indicate success 
 *
 *      Returns the slope of the line tangent to the xy-curve at the 
 *      mean (x,y) point in the specified density bin. If this bin is an
 *      endpoint of the xy-curve, then the slope is simply the slope of the
 *      line connecting to its nearest neighbor.
 *
 *      If determine_indep_var is zero, the x-variable in the binrec
 *      is fixed as the independent variable.  For a non-zero value
 *      of determine_indep_var, the magnitude of the slope and the value
 *      of the global variable, maxslope, determine whether x or y is
 *      used as the independent variable.
 *
 *      The value at flag_ptr is set according to the following:
 *
 *      flag_ptr :      1 = use x as independent variable;
 *                      0 = use y  "                 "
 *                     -1 = not enough points in bin or in adjacent bins
 *                          to determine the slope.
 *
 *      If *flag_ptr = -1, xyslope returns the value 0.0.
 */
{
   double x2, y2, x1, y1;
   double slope;
   struct binrec *bin1, *bin2, *binj;
   int b1, b2;
   int i;

   /* Focus on this particular bin. */
   binj = &bin[j];
       
   /* Check that there are enough points in this bin. */
   if (binj->nxy < MINPTS) {
     /* Too few points, exit with error flag, slope = 0. */
      *flag_ptr = -1;
      return (0.0);
   }

   /* Find nearest bins with enough points.
    * Look no further than 3 bins away and
    * don't mix reference levels */

   /* Initialize bin indices */
   b1 = b2 = 0;

   /* Search for the closest backward
    * (upper layer) bin by decrementing
    * the bin index. */
   i = j;
   /* Check that:
    *  (1) i is not too low
    *  (2) we have not yet found the upper bin
    *  (3) we are not too far from the center bin
    *  (4) we are still in the same reference level */
   while ((--i >= 0) && !b1 && ((j-i) <= 3) && (reflev[i] == reflev[j])) {
     /* Check the number of points */
     if (bin[i].nxy >= MINPTS) {
       /* We have enough points, set the upper bin
        * structure and index. */
       b1 = 1;
       bin1 = &bin[i];
     }
   }

   /* Search for the closest forward
    * (lower) bin by incrementing the
    * bin index. */
   i = j;
   /* Check that:
    *  (1) i is not too high
    *  (2) we have not yet found the upper bin
    *  (3) we are not too far from the center bin
    *  (4) we are still in the same reference level */
   while ((++i < nbins) && !b2 && ((i-j) <= 3) && (reflev[i] == reflev[j])) {
     /* Check the number of points */
      if (bin[i].nxy >= MINPTS) {
	/* We have enough points, set the lower bin
	 * structure and index. */
	b2 = 1;
	bin2 = &bin[i];
      }
   }

   /* If no neighbors (close enough, in same ref.
    * level, enough points) were found, try relaxing
    * the restrictions about ref level and bin
    * separation.  */
   if ( !b1 && !b2 ) {
     
     i = j;
     while ((--i >= 0) && !b1) {
       if (bin[i].nxy >= MINPTS) {
	 b1 = 1;
	 bin1 = &bin[i];
       }
     }   

     i = j;
     while ((++i < nbins) && !b2) {
       if (bin[i].nxy >= MINPTS) {
	 b2 = 1;
	 bin2 = &bin[i];
       }
     }
   }

   /* If no neighbors have enough points, then this
    * is the only bin and we cannot compute a slope. */  
   if ( !b1 && !b2) {
       *flag_ptr = -1;
       return (0.0);
   }

   
   /* If unable to find an appropriate bin on one
    * side, use this bin as an endpoint in
    * computing the slope. */   
   if (!b1) {
        bin1 = binj;
   }

   if (!b2) {
        bin2 = binj;
   }

   /* Compute mean x-y points for the 2 selected
    * bins and the slope of the line connecting
    * them.  This is like a centered difference
    * scheme. */
   x1 = bin1->xsum / bin1->nxy;
   y1 = bin1->ysum / bin1->nxy;
   x2 = bin2->xsum / bin2->nxy;
   y2 = bin2->ysum / bin2->nxy;

   /* Infinite slope check */
   if (x2 != x1) 
     slope = (y2 - y1) / (x2 - x1);
   else
     slope = 999999.0;

   if (plotflag ) {
     /* Write mean values to meanfile for plotting */
     fprintf(meanfile,"%.3lf %.3lf \n", (binj->ysum/binj->nxy), (binj->xsum/binj->nxy));
   }
   
   /* Check that the calculated slope is not larger
    * than the maxslope - if it is, make y the indep. var. */
   if ( determine_indep_var && (ABS(slope) > maxslope))  {
     *flag_ptr = 0;
     return (1/slope);
   }
   
   /* Return values with x as the indep. var. */
   *flag_ptr = 1;
   return (slope);
   
}  /* End xyslope() */


/****************************************************************************/

void plot_slope(int xfirst, double x, double m, double b, double e)
/*     
    xfirst     1 if x-variable is theta; 0 if salt 
         x     mean x for bin 
         m     slope 
         b     y-intercept 
         e     y-standard error 
*/
{
  if (xfirst) {
    fprintf(slopefile,"%.3lf %.3lf \n",  m*(x-1)+b, x-1);
    fprintf(slopefile,"%.3lf %.3lf  \n", (m*x+b), x);
    fprintf(slopefile,"%.3lf %.3lf \n", (m*(x+1)+b), x+1);
    fprintf(slopefile,"> \n");
    fprintf(stddevfile,"%.3lf %.3lf \n",  m*(x-1)+b+2*e, x-1);
    fprintf(stddevfile,"%.3lf %.3lf \n", (m*x+b+2*e), x);
    fprintf(stddevfile,"%.3lf %.3lf \n", (m*(x+1)+b+2*e), x+1);
    fprintf(stddevfile,"> \n");
    fprintf(stddevfile,"%.3lf %.3lf\n", m*(x-1)+b-2*e, x-1);
    fprintf(stddevfile,"%.3lf %.3lf \n", (m*x+b-2*e), x);
    fprintf(stddevfile,"%.3lf %.3lf \n", (m*(x+1)+b-2*e),x+1);
    fprintf(stddevfile,"> \n");
  }
  else {
    fprintf(slopefile,"%.3lf %.3lf \n",  x-.1, m*(x-.1)+b);
    fprintf(slopefile,"%.3lf %.3lf \n", x, m*x+b );
    fprintf(slopefile,"%.3lf %.3lf \n", x+.1, m*(x+.1)+b );
    fprintf(slopefile,"> \n");
    fprintf(stddevfile,"%.3lf %.3lf\n", x-.1, m*(x-.1)+b+2*e);
    fprintf(stddevfile,"%.3lf %.3lf\n", x, m*x+b+2*e);
    fprintf(stddevfile,"%.3lf %.3lf\n",x+.1, m*(x+.1)+b+2*e );
    fprintf(stddevfile,"> \n");
    fprintf(stddevfile,"%.3lf %.3lf \n", x-.1, m*(x-.1)+b-2*e);
    fprintf(stddevfile,"%.3lf %.3lf \n", x, m*x+b-2*e );
    fprintf(stddevfile,"%.3lf %.3lf \n", x+.1, m*(x+.1)+b-2*e );
    fprintf(stddevfile,"> \n");
  }
 
  return;
}  /* End plot_slope() */

/****************************************************************************/
