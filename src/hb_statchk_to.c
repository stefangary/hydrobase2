/*   hb_statchk_to.c
.
. USAGE:  hb_statchk_ts -Iinfile -Ooutfile -Sstatfile [-Ddeep_statfile] 
          -N#stderr_T-S [-M#stderr_surf/rhomax]
           [-Bbadfile] [-Llogfile]  [-Pplotfile_rootname]
.   Examines each data point in a file and eliminates points which are 
.   statistically distant from the t-ox relation defined for the bin.  
.   Linear theta-ox plus a standard error of estimate are computed 
.   for density (sigma) bins by the program hb_statfit_to.  These coefficients  
.   determine an acceptable range of values in each bin  -- n standard 
.   errors from the line.  If a t-ox point falls outside this range, it is 
.   assigned HB_MISSING. A point falling in no density bin or in a density bin containing less 
.   than 3 pts is given HB_MISSING.  A station bearing >10 (or >50%) bad ox observations is 
.   completely set to HB_MISSING.   The 
.   number of eliminated points is tallied as a percentage of the total number 
.   of scans and written to the logfile (or stdout).

*/

#include <stdio.h>
#include <stdlib.h>   /*sfg added to deal with exit and malloc warnings*/
#include <math.h>
#include <string.h>
#include "hydrobase.h"

#define   PRINT_MSG    1


#define   RED         2       /* pen colors */
#define   BLACK       1
#define   SQU         0       /* disspla symbol set */
#define   TRI         2
#define   X           4

#define   FLAG   -9.9
#define   PERCENT_LIM .5     /* percent of station allowed bad before
                                 declaring entire station bad */

/* ***********  global variables  ************ */

struct coefficients {
   char   xvar;         /* identifies independent variable: T or S */
   float  m1, b1, e1;   /* slope, intercept, std error for t-s */
   int    n1;
} ; 

struct coefficients *coeff;
float *sigmin, *sigmax;
float percent_lim;
int  nbins;
float  rhomax, nstderr_surf, nstderr;
FILE *logfile, *tsfile;
int  plotflag;
float *plot_s, *plot_t, *plot_sg;
int  *tspen, *sym;
char rflag;
/************************************/

/* prototypes for internally defined functions */

void    print_usage(char *);
int load_coeffs(FILE *);
int load_deep_coeffs(FILE *);
struct HYDRO_DATA *check_scans(int, struct HYDRO_HDR *, int *, int *,  struct HYDRO_DATA **);


main (int argc, char **argv)
{
   int     i;
   int     badsta, badscan, badox, goodox;
   int     totscan, totsta;
   int     status; 
   int     infile, outfile, badfile;
   char    iflag, oflag, sflag, nflag, dflag, bflag;
   char   *str, *statfilename, *deepfilename; 
   char   fname[80];
   struct HYDRO_HDR hdr;
   struct HYDRO_DATA  *newdata, *baddata;
   FILE   *statfile, *deepfile;


/* are there command line arguments? */

   if (argc < 4) {
      print_usage(argv[0]);
      exit(1);
   }

/*  set these default values... */

   rflag = iflag = oflag = sflag = dflag = bflag = 0;
   plotflag = 0;
   logfile = stdout;
   percent_lim = PERCENT_LIM;


/* initialize these... */

   hdr.prop_id = (int *) NULL;
   

/* parse command line arguments... */

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
                   fprintf(stderr,"\nOpened %s ", &argv[i][2]);
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
                  strncat(fname, "_to.dat",8);
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
            }  /* end switch */
       }  /* end if */
   }  /* end for */


   if (!iflag || !oflag || !sflag || !nflag) {
       print_usage(argv[0]);
       fprintf(stderr,"\n\nYou must specify I, O, S and N options!!!\n");
       exit(1);
   }

   nbins = load_coeffs(statfile);

   if (dflag) {
       load_deep_coeffs(deepfile);
   }
   if (percent_lim > 1.0)    /* change percent to fraction, if necessary */
      percent_lim /= 100.0;

   totsta = badscan = totscan = 0;
   badsta = 0;
   
/* loop for each station ... */

   while ((status = read_hydro_hdr(infile, &hdr)) == 0) {  
   
       newdata = check_scans(infile, &hdr, &goodox, &badox, &baddata);
       badscan += badox;
       totscan = goodox + badox;
       ++totsta;
       
       hdr.qual[2] = '1';

       write_hydro_station(outfile, &hdr, newdata);
       
       if (!goodox  && badscan )
          ++badsta;
	  
       if (bflag) {
         if ((hdr.nobs = baddata->nobs) > 0)
           write_hydro_station(badfile, &hdr, baddata);
       }
       
       get_separator(infile);

       /* recycle memory space ... */

       for (i = 0; i < hdr.nprops; ++i)  {
          free(newdata->observ[hdr.prop_id[i]]);
          free(baddata->observ[hdr.prop_id[i]]);
       }
       free(newdata);
       free(baddata);
       
   }  /*end while */ 
 
   report_status(status, stderr);
   close(infile);
   
/* print out tallies of eliminated data ... */

   fprintf(logfile, "%6d %6d (%5.1f%%)", badscan, totscan, (float)badscan*100./(float)totscan);
   fprintf(logfile, "   %6d %6d (%5.1f%%)", badsta, totsta, (float)badsta*100./(float)totsta);
   fprintf(logfile, "    : %s : %s\n", statfilename, deepfilename);

   fflush(logfile);

   if (plotflag) {
     fclose(tsfile);
   }

   fprintf(stderr,"\n\nEnd of %s. \n", argv[0]);
   exit(0);
}   /* end main */

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
int load_coeffs(FILE *file)
{
   int n, i, index;

   fscanf(file, "%d", &n);

   sigmin = (float *) malloc(n * sizeof(float));
   sigmax = (float *) malloc(n * sizeof(float));
   coeff = (struct coefficients *) malloc(n * sizeof(struct coefficients));

   for (i = 0; i < n; ++i) {
      fscanf(file,"%d%f%f %c", &index, &sigmin[i], &sigmax[i], &coeff[i].xvar);
      fscanf(file,"%f%f%f%d", &coeff[i].m1, &coeff[i].b1, &coeff[i].e1, &coeff[i].n1);
      if (index != i) {
         fprintf(stderr,"\nMismatch of index in load_coeffs().\n");
         exit(1);
      }
   }

   close(file);
   return(n);

} /* end load_coeffs() */

/****************************************************************************/
int load_deep_coeffs(FILE *file)
    /* Reads from a statfile and stores the values corresponding to
       sigma values > 30 (referenced to 2000 db or 4000 db) in the array,
       "coeff".  A little bit of checking is done to determine if the
       sigma bins match up with the ones established by load_coeffs(). */
{
   int n, i, index, junk;
   float s1, s2;

   fscanf(file, "%d", &n);
   if (n != nbins) {
      fprintf(stderr,"Deep statfile has different number of sigma bins.\n");
   }

   for (i = 0; i < n; ++i) {
      fscanf(file,"%d%f%f", &index, &s1, &s2);

      if (s1 > 30.) {   
         if ((s1 < (sigmin[index]-.0001)) || (s1 > (sigmin[index]+.0001)) ) {
            fprintf(stderr,"\nFatal Error: ");
            fprintf(stderr," Deep sigma bins do not match regular bins.\n");
            exit(1);
         }

         fscanf(file," %c%f%f%f%d", &coeff[index].xvar, &coeff[index].m1, 
                   &coeff[index].b1, &coeff[index].e1, &coeff[index].n1);
      }
      else {          /* skip over this line */
         fscanf(file," %*c%*f%*f%*f%*d", &junk);
      }
   } /* end for */

   return (n);


} /* end load_deep_coeffs() */
/****************************************************************************/

struct HYDRO_DATA *check_scans(int fd, struct HYDRO_HDR *hptr, int *ngoodox_ptr, int *nbadox_ptr, struct HYDRO_DATA **badaddr)
/*     
         fd   file descriptor for already open HydroBase file 
       hptr   ptr to HydroBase2 header info
  ngoodox_ptr   returns # of good scans here 
  nbadox_ptr   returns # of eliminated scans here 
    badaddr   address to return buffer of bad scans 
*/

   /* Returns a pointer to a struct HYDRO_DATA containing all the values
      which pass the statistical criteria.  If the entire station is
      eliminated, the field .nobs will have the value 0.     
      Also returns the bad scans in a struct HYDRO_DATA.  Memory is allocated
      here for both the bad and good scans and must be freed in the main program.
      
      */
{
   float nstde;
   double  pref, tref, p0 = 0.0;
   double  p, d, t, th, s, o, sig, x, y;
   double *scan;
   int bad;
   int n, nb, i, j, k, bin, found;
   struct HYDRO_DATA *dataptr, *badptr, *stabuffer;

/* initialize counters... */

    bad = 0;     
    n = nb = 0;

/* allocate memory for these structures ... */

    scan = (double *) calloc((size_t)hptr->nprops, sizeof(double));
    dataptr = (struct HYDRO_DATA *) calloc((size_t)1, sizeof(struct HYDRO_DATA));
    badptr = (struct HYDRO_DATA *) calloc((size_t)1, sizeof(struct HYDRO_DATA));
    stabuffer = (struct HYDRO_DATA *) calloc((size_t)1, sizeof(struct HYDRO_DATA));
    for (i = 0; i < hptr->nprops; ++i) {
       dataptr->observ[hptr->prop_id[i]] = (double *) calloc((size_t)hptr->nobs, sizeof(double));
       badptr->observ[hptr->prop_id[i]] = (double *) calloc((size_t)hptr->nobs, sizeof(double));
       stabuffer->observ[hptr->prop_id[i]] = (double *) calloc((size_t)hptr->nobs, sizeof(double));
    }
    dataptr->nobs = hptr->nobs;         /* for output */
    dataptr->nprops = hptr->nprops;
    badptr->nobs = 0;
    badptr->nprops = hptr->nprops;
    stabuffer->nprops = hptr->nprops;
    stabuffer->nobs = hptr->nobs;
    *ngoodox_ptr = 0;
    *nbadox_ptr = 0;

    if (plotflag) {
          plot_s = (float *) calloc((size_t)hptr->nobs, sizeof(float));
          plot_t = (float *) calloc((size_t)hptr->nobs, sizeof(float));
          plot_sg = (float *) calloc((size_t)hptr->nobs, sizeof(float));
          sym = (int *) calloc((size_t)hptr->nobs, sizeof(int));
          tspen = (int *) calloc((size_t)hptr->nobs, sizeof(int));
    }

/*  loop for each observation depth... */
       
    for (j = 0; j < hptr->nobs; ++j) {

       if (get_data_scan(fd, scan, hptr->nprops, hptr->prop_id) > 0) {
               fprintf(stderr, "ERROR reading data in get_data_scan()!");
               return(dataptr);
       }
       
       /* buffer each scan */
       
       for (k = 0; k < hptr->nprops; ++k) {
          stabuffer->observ[hptr->prop_id[k]][j] = scan[k];
	  dataptr->observ[hptr->prop_id[k]][j] = scan[k];
       }
       
       p = d = t = s  = o = FLAG;

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
                 case OX:
                      o = scan[i];
                      break;
                 default:
                      ;
          }  /* end switch */
       }
       
       
       if (o >= 0 ) {
          th = hb_theta(s, t, p, p0);

         /* determine ref pressure and compute density... */

          pref = 2000.0;
          if (p <= 1000.0)
                 pref = 0.0;
          if (p > 3000.0)
                 pref = 4000.0;
          
          tref = hb_theta(s, t, p, pref);
          hb_svan(s, tref, pref, &sig);

          if (plotflag) {
             plot_sg[j] = sig;
             plot_t[j] = th;
             plot_s[j] = o;
             sym[j] = (int) pref/1000;
             tspen[j] = BLACK;
          }
          
      /* search for density bin ... */

          nstde = nstderr;
          if (rflag) {   /* check if it is in the "surface" waters */
             if (sig < rhomax)
                nstde = nstderr_surf;
          }
          found = 0;
          bin = -1;
          while ((++bin < nbins) && !found) {
                if ((sig >= sigmin[bin]) && (sig < sigmax[bin]))
                    found = 1;
          }
          --bin;

    /* is it out of statistical range? */

          if (!found) {                  /* not in sigma bins */
                 ++bad;
		 dataptr->observ[(int)OX][j] = HB_MISSING;
                 for (i = 0; i < hptr->nprops; ++i) {
                     badptr->observ[hptr->prop_id[i]][nb] = scan[i];
		    
                 }
                 ++nb;
                 if (plotflag) {
                    tspen[j] = RED;
                 }
              
          }
          else if (coeff[bin].n1 < 3) {  /* not enough pts */
                 ++bad;
		 dataptr->observ[(int)OX][j] = HB_MISSING;
                 for (i = 0; i < hptr->nprops; ++i) {
                     badptr->observ[hptr->prop_id[i]][nb] = scan[i];
                 }
                 ++nb;
                 if (plotflag) {
                    tspen[j] = RED;
                 }
          }
          else  {                         /* check t-ox */
              x = th;
              y = o;
              if (coeff[bin].xvar == 'O') {
                   x = o;
                   y = th;
              }

              if (ABS(y - (x * coeff[bin].m1 + coeff[bin].b1)) > (nstde *coeff[bin].e1)) {
                 ++bad;                /* bad value */
		 dataptr->observ[(int)OX][j] = HB_MISSING;
                 for (i = 0; i < hptr->nprops; ++i) {
                  badptr->observ[hptr->prop_id[i]][nb] = scan[i];
                 }
                 ++nb;
                 if (plotflag) {
                    tspen[j] = RED;
                 }
              }
              else {                   /* good value */
               ++(*ngoodox_ptr);;
             }
          }
       } /* end if o > 0 */	  
    } /* end for */


   /* check the tallies... if # of bad pts > 50%, eliminate the
      entire station. */
    
    if (bad >  hptr->nobs * percent_lim) {
       *nbadox_ptr = hptr->nobs;
       free(badptr);
       badptr = stabuffer;    /* replace the bad scans with the buffered station */
    }
    else {
       badptr->nobs = nb;
       free(stabuffer);       /* free the buffer if station has good scans */
       *nbadox_ptr = bad;
    }

    if (plotflag) {

       if (bad > hptr->nobs * percent_lim) {
          for (i = 0; i < hptr->nobs; ++i) {
             tspen[i] = RED;
          }
       }

       for (i = 0; i < hptr->nobs; ++i) {
          fprintf(tsfile, " %6.3f %6.3f %6.3f %3d %3d \n", plot_s[i], plot_t[i], 
                     plot_sg[i], tspen[i], sym[i]);
       }

       free(plot_s);
       free(plot_t);
       free(plot_sg);
       free(tspen);
       free(sym);
    }  /* end if plotflag */
   

    free(scan);
    *badaddr = badptr;
    return(dataptr);

} /* end check_scans() */
