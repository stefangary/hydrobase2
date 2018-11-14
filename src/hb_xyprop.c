/* hb_xyprop.c
................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             1993
                             updated to ANSI standard May2000

			     Modified Stefan Gary April 2010
			     to add -P and -T flags.
................................................................................
................................................................................
.   Creates a list of  property pairs for each observation
.   level in a file of hydro stations.
.   -M Separates stations with a '>' to permit points to be
.      connected in GMT.
.   -D delimits a depth range for properties.
.   -F option specifies a minimum value of X,Y as a limit
.      (default is -9 to avoid outputting missing values).
.   -P option adds position info at each line (lon,lat,dep)
.   -T option adds time info at each line (year,month,day)
.  If no input files are specified, reads from STDIN; writes to STDOUT
................................................................................
................................................................................
*/

#include <stdio.h>
#include <stdlib.h>   /*sfg added to deal with exit and malloc warnings*/
#include <string.h>
#include "hydrobase.h"


/* input file pathnames */

#define    EXTENT   ""
#define    DIR      ""


double flag;

/*  prototypes for locally defined functions */

void print_usage(char *);

int main (int argc, char **argv)

{
   double  d;
   struct HYDRO_HDR hdr;
   struct HYDRO_DATA data;
   int    print_msg = 1;       /* print message in file open routines*/
   int  i, n, error;
   int nfiles, curfile;
   short dflag, xflag, yflag, mflag;
   short pflag, tflag;
   int mindepth, maxdepth;
   int scan_ok, any_scans, status;
   int infile;
   int xindex, yindex;
   enum property X, Y;
   char *s, str_format[30];
   char *dir, *extent;


/* check for command line arguments */

   if (argc < 2) {
     print_usage(argv[0]);
     exit(1);
   }

/* initialize these */

   xflag = yflag = dflag = mflag = 0;
   pflag = tflag = 0;
   flag = HB_MISSING + 0.1;
   dir = DIR;
   extent = EXTENT;
   infile = STDIN;
   for (i = 0; i < MAXPROP; ++i)
      data.observ[i] = (double *) NULL;
   hdr.prop_id = (int *) NULL;
   nfiles = 0;
   curfile = 1;
   error = 0;

/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
     if (argv[i][0] == '-') {
       switch (argv[i][1]) {
       case 'D':                   /* Get input dir */
	 dir = &argv[i][2];
	 break;
       case 'E':                   /* Get file extent */
	 extent = &argv[i][2];
	 break;
       case 'F':                   /* Set minimum acceptable value*/
	 s = &argv[i][2];  
	 error += (sscanf(s,"%lf", &flag) != 1);
	 break;
       case 'M':                   /* Specify multiple segment option */
	 mflag = 1;
	 break;
       case 'P':                   /* User requests position information */
	 pflag = 1;
	 break;
       case 'T':                   /* User requests time information */
	 tflag = 1;
	 break;
       case 'X':                   /* Get x-property */
	 xflag = 1;
	 xindex = get_prop_indx(&argv[i][2]);
	 if (xindex < 0) {
	   fprintf(stderr,"\nUnknown X property requested: %s\n", &argv[i][2]);
	   exit(1); 
	 }
	 X = (enum property) xindex;
	 break;
       case 'Y':                   /* Get y-property */
	 yflag = 1;
	 yindex = get_prop_indx(&argv[i][2]);
	 if (yindex < 0) {
	   fprintf(stderr,"\nUnknown Y property requested: %s\n", &argv[i][2]);
	   exit(1); 
	 }
	 Y = (enum property) yindex;
	 break;
       case 'Z':                   /* Get depth range */
	 dflag = 1;
	 s = &argv[i][2];  
	 error += (sscanf(s,"%d", &mindepth) != 1);
	 s = strchr(s,'/');  /* check for another delimiter*/
	 error = (s == NULL);
	 if (s != NULL) {
	   ++s;  /* move past delimiter */
	   error += (sscanf(s,"%d", &maxdepth) != 1);
	 }
	 if (error) break;
	 if (maxdepth < mindepth) {
	   fprintf(stderr,"\nmindepth exceeds maxdepth\n");
	   exit(1);
	 }
	 break;
       case 'h':                    /* Print help */
	 print_usage(argv[0]);
	 exit(0);
       default:
	 error = 1; 
       }    /* End flag parsing switch */
       
       if (error) {
	 fprintf(stderr,"\nError parsing command line args.\n");
	 fprintf(stderr,"     in particular: '%s'\n", argv[i]);
	 fprintf(stderr,"Use -h for a complete USAGE message.\n");
	 exit(1);
       }
       
     }  /* End if we encounter a '-' (a flag). */
       
     else {
          ++nfiles;
     }
   }  /* End for loop over command line tokens. */

   if (! (xflag && yflag) ) {
       fprintf(stderr,"\nYou must specify -X and -Y properties ...\n");
       exit(1);
   }
   
   /* Original line includes newline character in str_format.
   sprintf(str_format," %c%d.%dlf  %c%d.%dlf \n", '%',get_field_width(xindex),
       get_field_precis(xindex), '%', get_field_width(yindex),
       get_field_precis(yindex));
   This was moved to below so that additional information may be
   added to each line.
   */

   sprintf(str_format," %c%d.%dlf  %c%d.%dlf ", '%',get_field_width(xindex),
       get_field_precis(xindex), '%', get_field_width(yindex),
       get_field_precis(yindex));
   
   /* Loop to read each file */
   do {
   
     if (nfiles) {
       infile = open_hydro_file(dir, argv[curfile], extent, print_msg);
       if (infile < 0)
	 goto NEXTFILE;
     }
     else {
       fprintf(stderr,"Expecting station input from stdin...\n" ); 
     }
     
     /* Loop to read each station */
     while ((status = get_station(infile, &hdr, &data)) == 0) {  
       any_scans = 0;
       if (available(X, &hdr) && available(Y, &hdr)) {
	 for (n = 0; n < hdr.nobs; ++n) {
          
	   /* Check minimum value of data in scan */
	   scan_ok = (data.observ[xindex][n] > flag) && ( data.observ[yindex][n] > flag);

	   /* Check depth range of scan */
	   if (scan_ok && dflag) {
	     d = data.observ[(int)DE][n];
	     scan_ok = (d >= (double)mindepth) && (d <= (double)maxdepth);
	   }
           
	   /* If the scan is good, write to output file.*/
	   if (scan_ok) {
	     ++any_scans;
	     fprintf(stdout,str_format, data.observ[xindex][n],
		     data.observ[yindex][n]);

	     /* Check for position information request */
	     if (pflag) {
	       fprintf(stdout," %10.4lf %10.4lf %10.3lf ",hdr.lon,hdr.lat,data.observ[(int)DE][n]);
	     } /* End of position info request check.

	     /* Check for time information request */
	     if (tflag) {
	       fprintf(stdout," %4d %4d %4d ",hdr.year,hdr.month,hdr.day);
	     } /* End of time info request check.

	     /* Add a new line after all possible
	      * output for this scan is written. */
	     fprintf(stdout,"\n");

	   }   /* End check to write output */

	 }   /* End for loop over all scans */
          
	 if (any_scans && mflag) {           /* Station separator */
	   fprintf(stdout,"> \n" ); 
	 } /* End if check on whether to output a station delimitor */
	 
       }   /* End if check on availability of X and Y in this station */ 
     }   /* End while loop over each station in the current file */
       
     report_status(status, stderr);
 
   NEXTFILE:
     if (nfiles) {
	close(infile);
     }
   
   } while (curfile++ < nfiles );  /* End of loop over all files */

}  /* End main program */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nCreates a listing of x,y property pairs.  Use hb_propcalc");
   fprintf(stderr,"\nto compute properties and pipe the output to");
   fprintf(stderr,"\nthis utility.  Use the -M option to separate profiles with a >.");
   fprintf(stderr,"\nUse the -Z option to delimit a depth range.\n");
   fprintf(stderr,"\nUsage:  %s -X<x_property> -Y<y_property> [-D<dirname>] [-E<file_extent>] [-Zmindepth/maxdepth] [-M] [-P] [-T]", program);
   fprintf(stderr,"\n\n    -X  : 2-char property mnemonic to specify X");
   fprintf(stderr,"\n    -Y  : 2-char property mnemonic to specify Y");
   fprintf(stderr,"\n   OPTIONS:");
   fprintf(stderr,"\n   [-D] : specifies directory for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n   [-F] : specify flag (minimum value) for output.  default [%.2lf]", flag);
   fprintf(stderr,"\n   [-M] : separate stations with > symbol");
   fprintf(stderr,"\n   [-P] : appends lon, lat, depth to each line.");
   fprintf(stderr,"\n   [-T] : appends year, month, day to each line.");
   fprintf(stderr,"\n   [-Z] : optional depth delimiter for output");
   fprintf(stderr,"\n   [-h] : help -- prints this message");
   fprintf(stderr,"\n\n");  
   return;
}

