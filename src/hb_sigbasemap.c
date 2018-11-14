/*  hb_sigbasemap.c
***********************************************************************
*  Usage:  hb_sigbasemap -O<outfile>
*                        -S/<smin>/<smax>/<sinc>
*                        -T/<tmin>/<tmax>/<tinc>
*                        -P<refp>
*                       [-h]
*
*  This program computes an xyz file which is the
*  potential density at the given reference
*  pressure (-P) over the given salinity (-S) and
*  temperature (-T) intervals.
*
*  Output goes to the output file (-O), all with
*  columnar format and no headers.  The format is:
*
*  salinty temperature sigma
*
**************************************************************************
*/ 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hydrobase.h"

/* Prototypes for locally defined functions */

/* Print usage if user error. */
void print_usage(char *);

/* Get axis information from command line. */
void get_axis(char*, double *, double *, double *);

main(int argc, char **argv) {

  /* Define array pointers to store min/max
   * temperature and salinity as well as arrays. */

  int    ntval, nsval; /* Number of nodes */
  double *tval, *sval; /* Arrays of values */
  double tmin, smin;   /* Minima */
  double tmax, smax;   /* Maxima */
  double tinc, sinc;   /* Increments */

  /* Local counter variables */
  int    i, j;

  /* Temporary holders of single values,
   * t = temp, s = salt, tref = potential
   * temp at reference pressure, sigma =
   * specific volume anamoly (pot dens). */
  double t, s, tref, sigma, refp;

  /* Flags to keep track of whether the user
   * has specified -T, -S, and -P flags. */
  int    tflag;
  int    sflag;
  int    pflag;

  /* File ID of output file (txt xyz files) */
  FILE   *outfile;
  
  /* Pointer to a generic
   * string for reading command line. */
  char  *st;
  
  /* Error flag for reading command line */
  int error;

  /*-----------------------------------------------------*/
  /* Look at the command line. */
  /*-----------------------------------------------------*/

  /* Are there command line arguments? */
  if (argc < 2) {
    print_usage(argv[0]);
    exit(1);
  }

  /* Set default values... */
  outfile = (FILE *)NULL;
  tflag = sflag = pflag = 0; /* Initialize command line arg flags to false */
   
  /* Parse command line arguments... */
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
               case 'O':
		 outfile = fopen(&argv[i][2],"w");
		 if (outfile == NULL) {
		   fprintf(stderr,"\nError opening %s for output.\n",&argv[i][2]);
		   exit(1);
		 }
		 break;

               case 'S':
		 sflag = 1;
		 get_axis(&argv[i][2],&smin,&smax,&sinc);
		 /*fprintf(stderr,"S: [%lf : %lf : %lf]\n",smin,sinc,smax);*/
		 break;

               case 'T':
		 tflag = 1;
		 get_axis(&argv[i][2],&tmin,&tmax,&tinc);
		 /*fprintf(stderr,"T: [%lf : %lf : %lf]\n",tmin,tinc,tmax);*/
		 break;

               case 'P':
		 pflag = 1;
		 error = (sscanf(&argv[i][2], "%lf", &refp) != 1);
		 /*fprintf(stderr,"refp is: %lf\n",refp);*/
		 if (error != 0) {
		   print_usage(argv[0]);
		   fprintf(stderr,"Error reading reference pressure, refp: %lf",refp);
		   exit(0);
		 }
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
      fprintf(stderr,"\nError parsing command line");
      fprintf(stderr,"\n in particular: %s\n", argv[i]);
      fprintf(stderr,"Use -h for help\n");
      exit(1);
    }
  }  /* end for */
  
  /* Check that the user has specifed an outfile */
  if ( outfile == NULL ) {
    print_usage(argv[0]);
    fprintf(stderr,"\n\nYou must specify the outfile with -O.\n");
    exit(1);
  }

  /* Check that the user has specified -T, -S, and -P flags. */
  if ( tflag == 0 || sflag == 0 || pflag == 0 ) {
    print_usage(argv[0]);
    fprintf(stderr,"\n\nYou must specify all -T, -S, and -P flags.\n");
    exit(1);
  }

  /*-----------------------------------------------------*/
  /* Initialize the value arrays. */
  /*-----------------------------------------------------*/

  /* Compute the number of values - truncate round down.
   * Add one extra value for the endpoint, i.e.:
   * 1 2 3 . 4 . 5 . 6 . 7 8 9 10, choose min = 3, max = 7, inc = 0.5
   *     ^               ^
   * (7 - 3)/0.5 = 8
   *
   * But there are 9 numbers between (and including) 3 and 7
   * with a spacing of 0.5.*/
  ntval = (int)((tmax - tmin)/tinc) + 1;
  nsval = (int)((smax - smin)/sinc) + 1;

  /* Recompute the increment if the user gave an uneven inc. */
  tinc = (tmax - tmin)/((double)ntval);
  sinc = (smax - smin)/((double)nsval);

  /*fprintf(stderr,"T: [%lf : %lf : %lf]\n",tmin,tinc,tmax);
    fprintf(stderr,"S: [%lf : %lf : %lf]\n",smin,sinc,smax);*/

  tval = (double *) calloc(ntval, sizeof(double));
  sval = (double *) calloc(nsval, sizeof(double));

  /* Compute the values. */
  tval[0] = tmin;
  for ( i = 1; i < ntval; i++) {
    tval[i] = tval[i-1] + tinc;
    /*fprintf(stderr,"\nT values are: %lf",tval[i]);*/
  }

  sval[0] = smin;
  for ( i = 1; i < nsval; i++) {
    sval[i] = sval[i-1] + sinc;
    /*fprintf(stderr,"\nS values are: %lf",sval[i]);*/
  }

  /*-----------------------------------------------------*/
  /* Compute sigma and write to file. */
  /*-----------------------------------------------------*/

  for ( i = 0; i < ntval; i++ ) {
    
    for ( j = 0; j < nsval; j++ ) {

      /* Assign single values. */
      t = tval[i];
      s = sval[j];

      /* Test if sigma can be computed.
       * For weird values of t and s,
       * skip to the next level with
       * the continue command. */
      if ((t < -8) || (s < -8)) {
	continue;
      }
     
      /* Compute potential temp wrt reference depth */
      tref = hb_theta(s, t, refp, refp);
    
      /* Compute sigma, referenced to refp in both t and p */
      hb_svan(s, tref, refp, &sigma); 

      fprintf(outfile,"%lf %lf %lf\n",s,t,sigma);
    }
  }
  
  /*-----------------------------------------------------*/

  fprintf(stderr,"\nEnd of %s\n", argv[0]);
  exit(0);
}  /* End main */

/****************************************************************************/

void print_usage(char *program)
{
  fprintf(stderr,"\n --------------------------------------------------------");
  fprintf(stderr,"\n Usage: %s -O<outfile>",program);
  fprintf(stderr,"\n                    -T/<tmin>/<tmax>/<tinc>");
  fprintf(stderr,"\n                    -S/<smin>/<smax>/<sinc>");
  fprintf(stderr,"\n                    -P<refp>");
  fprintf(stderr,"\n                   [-h]");
  fprintf(stderr,"\n");
  fprintf(stderr,"\n   -O  : specifies output file for density map.  ");
  fprintf(stderr,"\n        ex: -Obasemap.xyz ");
  fprintf(stderr,"\n   -T  : specifies temperature input values [oC].");
  fprintf(stderr,"\n        ex: -T/0/30/0.1");
  fprintf(stderr,"\n   -S  : specifies salinity input values [PSU].");
  fprintf(stderr,"\n        ex: -S/34.4/37.6/0.01");
  fprintf(stderr,"\n   -P  : specifies reference pressure [dbar ~ m].");
  fprintf(stderr,"\n        ex: -P2000");
  fprintf(stderr,"\n  [-h] : help -- prints this message. ");
  fprintf(stderr,"\n");
  fprintf(stderr,"\n This program computes an xyz file which is the");
  fprintf(stderr,"\n potential density at the given reference");
  fprintf(stderr,"\n pressure (-P) over the given salinity (-S) and");
  fprintf(stderr,"\n temperature (-T) intervals.\n\n");
  fprintf(stderr,"\n Output goes to the output file (-O), all with");
  fprintf(stderr,"\n columnar format and no headers.  The format is:\n");
  fprintf(stderr,"\n salinty temperature sigma\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"\n --------------------------------------------------------\n");
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
