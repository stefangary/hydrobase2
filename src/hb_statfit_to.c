/*  hb_statfit_to.c
................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth G Curry
                             Woods Hole Oceanographic Instition
                             updated 2000 to ANSI standards
................................................................................
-------------------------------------------------------------------------------
.  USAGE:  hb_statfit  list_of_files -Ooutfile -Ssigbin_file [-Ddir] [-Eextent] [-Mmaxslope] [-Xminpts] [-Pplotfiles_root_name]
_______________________________________________________________________________
.   Computes linear fits to the t,s data for each sigma bin specified  
.   in the sigbin file.  The theta-salt curve is approximated by connecting
.   the mean t,s point in each density bin.  The line tangent to this curve
.   at the mean t,s point then approximates the t-s relation in each bin. 
.   We determine the slope, y-intercept, and standard error for each bin.
.   Initially, we assume theta as the independent variable. If the slope of
.   the linear fit exceeds a critical value (maxslope), we designate salt
.   as the independent variable and compute the statistics for theta = F(salt)
.   and oxygen = F(salt).  This helps to identify situations where the 
.   F(theta) fit for a particular density bin does not approximate the 
.   overall t-s relationship.
.
.     
.   The output file contains 1 header line followed by a line for each density
.   bin:
.
.   number of bins
.   bin #, T or S,  slope, intercept, error, n  
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

/* ******  global variables ****** */

struct binrec  *FofT_bins, *FofS_bins;
        float  *sigmin, *sigmax;
	  int  *reflev;
          int  nbins;
       double  maxslope;
         FILE *slopefile, *meanfile;
          int  plotflag;  
/*******************************/

/* prototypes for locally defined functions */

void print_usage(char *);
int load_bins(char *);
void initialize(struct binrec *, int);
void insertdata(struct HYDRO_HDR *,struct HYDRO_DATA *);
double xyslope(struct binrec *, int, int, int *);
void plot_slope(FILE *, int, double, double, double, double);

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

/* are there command line arguments? */

   if (argc < 4) {
      print_usage(argv[0]);
      exit(1);
   }

/*  set these default values... */

   oflag = sflag = 0;
   plotflag = 0;
   minpts = MINPTS;
   maxslope = MAXSLOPE;
   dir = DIR;
   extent = EXTENT;

/* initialize some variables */
   for (i = 0; i < MAXPROP; ++i) {
     data.observ[i] = (double *) NULL;
   }
   hdr.prop_id = (int *)NULL;

/* parse command line arguments... */

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
                  strncat(fname,".mean.to",9);
                  meanfile = fopen(fname,"w");
                  if (meanfile == NULL) {
                     fprintf(stderr,"\nUnable to open %s\n", fname);
                     exit(1);
                  }

                  strcpy(fname,&argv[i][2]);
                  strncat(fname,".slope.to",10);
                  slopefile = fopen(fname,"w");
                  if (slopefile == NULL) {
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

   }  /* end for */


   if (!nfiles || !oflag || !sflag) {
       print_usage(argv[0]);
       fprintf(stderr,"\n\nYou must specify input, output, and sigma bin files!!!\n");
       exit(1);
   }


/* define sigma bins ... */

   if ((nbins = load_bins(binfilename)) < 0) {
       exit(1);
   }

   initialize(FofT_bins, nbins);
   initialize(FofS_bins, nbins);

   fprintf(outfile,"%d", nbins);


/* loop for each file */

   do {
      infile = open_hydro_file(dir, argv[curfile], extent, PRINT_MSG);
      if (infile  < 0) 
       goto NEXTFILE;
       
   /* loop for each station */

      while ((status = get_station(infile, &hdr, &data)) == 0) {  

          insertdata(&hdr, &data);
       
      }  /*end while */ 
 
      report_status(status, stderr);
      close(infile);
NEXTFILE:
      ;
   } while (curfile++ < nfiles);


/*  compute slope, intercept and standard deviation and write to outfile ... 
    check each bin for existing pts; there must be at least 3 pts to compute
    the standard error ... */

     fprintf(stderr,"\n computing slope and variance ...\n");

     for (i = 0; i < nbins; ++i) {
       fprintf(outfile,"\n%4d %6.3f %6.3f", i, sigmin[i], sigmax[i]);

       if ((n = FofT_bins[i].nxy) >= 3) {

          check_x = 1;
          m = xyslope(FofT_bins, i, check_x, &t_is_x);

          if (t_is_x >= 0) {

             bin = t_is_x ? FofT_bins : FofS_bins;
             ybar = bin[i].ysum / n;
             xbar = bin[i].xsum / n;
             b = ybar - m * xbar;
             E = ABS((bin[i].yysum - n*ybar*ybar - 2*m*(bin[i].xysum  
                - n*ybar*xbar) + m*m*(bin[i].xxsum - n*xbar*xbar)) / (n-2));
             E = sqrt(E);
          
             if (t_is_x)
                 fprintf(outfile," T");
             else 
                 fprintf(outfile," O");

             fprintf(outfile," %10.5lf %7.3lf %8.4lf %5d",  m, b, E, n);

             if (plotflag) {
               plot_slope(slopefile, t_is_x, bin[i].xsum/bin[i].nxy, m, b, E);
             }


          } /* end if */
       }   /* end if */
       else {
         fprintf(outfile," X  0 0 0 0");
       }
     }  /* end for */


   exit(0);

}  /* end of main() */

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
   fprintf(stderr,"\n          to which will be appended mean.ts and slope.ts");
   fprintf(stderr,"\n          ex: -P7207");
   fprintf(stderr,"\n\n");  
   return;
} /* end print_usage() */


/****************************************************************************/
int load_bins(char *name)

   /* reads the binfile and loads the bin definitions into the global 
      variables.  Returns the number of bins or -1 if an
      error occurs.  In case of an error, an appropriate message is 
      written to the stderr device. */
{
   int n, i;
   float r;
   FILE *binfile;


   binfile = fopen(name,"r");
   if (binfile == NULL) {
       fprintf(stderr,"\nUnable to open Sigma Bin File: %s \n", name);
       return(-1);
   }

/* determine number of bins defined in file and set up structures ... */

   n = 0;
   while (fscanf(binfile,"%*f %f", &r) == 1) 
        ++n;

   close(binfile);

   sigmin = (float *) malloc(n * sizeof(float));
   sigmax = (float *) malloc(n * sizeof(float));
   reflev = (int *) malloc(n * sizeof(int));
   FofT_bins = (struct binrec *) malloc(n * sizeof(struct binrec));
   FofS_bins = (struct binrec *) malloc(n * sizeof(struct binrec));

   binfile = fopen(name,"r");
   i = 0;
   while (fscanf(binfile,"%f %f", &sigmin[i], &sigmax[i]) == 2) {
       reflev[i] = 2000;
       if (sigmin[i] < 30)
	        reflev[i] = 0;
       if (sigmin[i] > 40)
	        reflev[i] = 4000;
       ++i;
   }

   close(binfile);
   if (n != i) {
      fprintf(stderr, "\nNumber of bins defined: %d differs from what was allocated for: %d", i, n);
      return(-1);
   }
   return (n);

}  /* end load_bins() */
/****************************************************************************/
void initialize(struct binrec *b, int n)
/*       b    pointer to start of array of binrecs 
         n    number of binrecs */
{
   int i;

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

void insertdata(struct HYDRO_HDR *hptr, struct HYDRO_DATA *dptr)
{
   double    th, pref, tref, sigma, p0 = 0.0;
   double    t, s, p, d, o;
   register  i, j;
   char      found;

   if (!available((int)PR, hptr) || 
       !available((int)TE, hptr) || 
       !available((int)OX, hptr) )
       return;


/* loop for each observation level... */

     for (j = 0; j < hptr->nobs; ++j) {

        t = dptr->observ[(int)TE][j];
        s = dptr->observ[(int)SA][j];
        p = dptr->observ[(int)PR][j];
        d = dptr->observ[(int)DE][j];
        o = dptr->observ[(int)OX][j];
	
	if (o < 0)
	   continue;

        if ((t < -8) || (s < -8)) {          /* can sigma be computed? */
           continue;
        }

        pref = 2000.;
        if (d <= 1000.)  /* determine ref pressure */
              pref = 0.;
        if (d > 3000.)
              pref = 4000.;

        th = hb_theta(s, t, p, p0);
        tref = hb_theta(s, t, p, pref);   /* compute temp, referenced */
        hb_svan(s, tref, pref, &sigma);   /* compute sigma, referenced */

        i =  0;  
        found = 0;
        do {                           /* search for correct bin */
            if ((sigma >= sigmin[i]) && (sigma < sigmax[i]))
               found = 1;
        } while ( (++i < nbins) &&  !found );
        --i; 

        if (!found) {                  /* write message to stderr device */
            fprintf(stderr,
           "\nSigma out of range: %9.3lf :  d=%5.0lf t=%6.2lf s=%6.2lf", 
            sigma, d, t, o);
        }
        else {                         /* increment counters */
            FofT_bins[i].xsum += th;
            FofS_bins[i].ysum += th;
            FofT_bins[i].xxsum += th * th;
            FofS_bins[i].yysum += th * th;
            ++FofT_bins[i].nxy;
            ++FofS_bins[i].nxy;
            FofT_bins[i].ysum += o;
            FofS_bins[i].xsum += o;
            FofT_bins[i].yysum += o * o;
            FofS_bins[i].xxsum += o * o;
            FofT_bins[i].xysum += th * o;
            FofS_bins[i].xysum += th * o;
         } 

     }  /* end for */

   return;

}  /* end insertdata() */

/****************************************************************************/
double xyslope(struct binrec *bin, int j, int determine_indep_var, int *flag_ptr)
/*  
                  bin   pointer to array of binrecs
                    j   index to bin array being considered  
  determine_indep_var   0 or 1 
             flag_ptr   returned value to indicate success 

 Returns the slope of the line tangent to the xy-curve at the 
       mean (x,y) point in the specified density bin. If this bin is an
       endpoint of the xy-curve, then the slope is simply the slope of the
       line connecting to its nearest neighbor. If determine_indep_var 
       is zero, the x-variable in the binrec is fixed as the independent
       variable.  For a non-zero value of determine_indep_var, the magnitude 
       of the slope and the value of the global variable, maxslope, 
       determine whether  x or y is used as the independent variable.
       The value at flag_ptr is set according to the following:

         flag_ptr :      1 = use x as independent variable;
                         0 = use y  "                 "
                        -1 = not enough points in bin or in adjacent bins
                             to determine the slope.

       If *flag_ptr = -1, xyslope returns the value 0.0.
*/
{
   double x2, y2, x1, y1;
   double slope;
   struct binrec *bin1, *bin2, *binj;
   int b1, b2;
   int i;



   binj = &bin[j];
       
 /*  enough points in this bin? */

   if (binj->nxy < MINPTS) {
      *flag_ptr = -1;
      return (0.0);
   }

 /*  find nearest bins with enough points...
       [look no further than 3 sig bins away and
         don't mix reference levels] */
   
   b1 = b2 = 0;

   i = j;
   while ((--i >= 0) && !b1 && ((j-i) <= 3) && (reflev[i] == reflev[j])) {
      if (bin[i].nxy >= MINPTS) {
           b1 = 1;
           bin1 = &bin[i];
      }
   }

   i = j;
   while ((++i < nbins) && !b2 && ((i-j) <= 3) && (reflev[i] == reflev[j])) {
      if (bin[i].nxy >= MINPTS) {
           b2 = 1;
           bin2 = &bin[i];
      }
   }

/*  if no neighbors, try relaxing the restrictions about ref level  */

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

/* if no neighbors have enough points, then this is the only bin and we
   cannot compute a slope. */  

   if ( !b1 && !b2) {
       *flag_ptr = -1;
       return (0.0);
   }

   
/* if unable to find an appropriate bin on one side, use this bin as an
   endpoint in computing the slope  ... */   
   if (!b1) {
        bin1 = binj;
   }

   if (!b2) {
        bin2 = binj;
   }

/* compute mean x-y points for the 2 selected bins and the slope of the line
   connecting them... */

   x1 = bin1->xsum / bin1->nxy;
   y1 = bin1->ysum / bin1->nxy;
   x2 = bin2->xsum / bin2->nxy;
   y2 = bin2->ysum / bin2->nxy;
   if (x2 != x1) 
      slope = (y2 - y1) / (x2 - x1);
   else
      slope = 999999.0;

   if (plotflag ) {
        fprintf(meanfile,"%.3lf %.3lf \n", (binj->ysum/binj->nxy), (binj->xsum/binj->nxy));
   }
   
   if ( determine_indep_var && (ABS(slope) > maxslope))  {
      *flag_ptr = 0;
      return (1/slope);
   }

   *flag_ptr = 1;
   return (slope);

}  /* end xyslope() */


/****************************************************************************/

void plot_slope(FILE *plotfile, int xfirst, double x, double m, double b, double e)
/*     
    xfirst     1 if x-variable is theta; 0 if salt 
         x     mean x for bin 
         m     slope 
         b     y-intercept 
         e     y-standard error 
*/
{
  if (xfirst) {
    fprintf(plotfile,"%.3lf %.3lf \n",  m*(x-1)+b, x-1);
    fprintf(plotfile,"%.3lf %.3lf  \n", (m*x+b), x);
    fprintf(plotfile,"%.3lf %.3lf \n", (m*(x+1)+b), x+1);
    fprintf(plotfile,"> \n");
    fprintf(plotfile,"%.3lf %.3lf \n",  m*(x-1)+b+2*e, x-1);
    fprintf(plotfile,"%.3lf %.3lf \n", (m*x+b+2*e), x);
    fprintf(plotfile,"%.3lf %.3lf \n", (m*(x+1)+b+2*e), x+1);
    fprintf(plotfile,"> \n");
    fprintf(plotfile,"%.3lf %.3lf\n", m*(x-1)+b-2*e, x-1);
    fprintf(plotfile,"%.3lf %.3lf \n", (m*x+b-2*e), x);
    fprintf(plotfile,"%.3lf %.3lf \n", (m*(x+1)+b-2*e),x+1);
    fprintf(plotfile,"> \n");
  }
  else {
    fprintf(plotfile,"%.3lf %.3lf \n",  x-.1, m*(x-.1)+b);
    fprintf(plotfile,"%.3lf %.3lf \n", x, m*x+b );
    fprintf(plotfile,"%.3lf %.3lf \n", x+.1, m*(x+.1)+b );
    fprintf(plotfile,"> \n");
    fprintf(plotfile,"%.3lf %.3lf\n", x-.1, m*(x-.1)+b+2*e);
    fprintf(plotfile,"%.3lf %.3lf\n", x, m*x+b+2*e);
    fprintf(plotfile,"%.3lf %.3lf\n",x+.1, m*(x+.1)+b+2*e );
    fprintf(plotfile,"> \n");
    fprintf(plotfile,"%.3lf %.3lf \n", x-.1, m*(x-.1)+b-2*e);
    fprintf(plotfile,"%.3lf %.3lf \n", x, m*x+b-2*e );
    fprintf(plotfile,"%.3lf %.3lf \n", x+.1, m*(x+.1)+b-2*e );
    fprintf(plotfile,"> \n");
  }
 
  return;
}  /* end plot_slope() */
