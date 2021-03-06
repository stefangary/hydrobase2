/*  hb_msnum.c
------------------------------------------------------------------------
* Reads lon lat from command line and returns the WMO square number.
*   
*   USAGE: hb_msnum <lon> <lat>
* -----------------------------------------------------------------------
* Originally based on wod01_convert.c, part of the original
* Hydrobase2 distribution.
* Stefan Gary, December 2020
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hydrobase.h"

#define   BLANK     0x20   /* ascii code for blank */
#define   MISSING   -9.9   /* missing value flag */
#define   DIR    ""
#define   EXTENT ""

/************************************************************

 HydroBase globally defined variables

 The structure HYDRO_HDR contains:
       char country[3], ship[3]
       char origin, instrument
       int  cruise, station, year, month, day
       float lat, lon
       int  pdr, nobs, nprops, ms10, ms1
       int *prop_id
       char qual[nqual]

The structure HYDRO_DATA contains:
       int nobs, nprops
       double *observ[MAXPROP]

*************************************************************/
  
struct HYDRO_HDR  hdr, bad_hdr;
struct HYDRO_DATA data, bad_data;

void print_usage(char *);

/***********************************************************/

/* Getting command line values at execution:
 * argc = number of parsed entries on command line
 *        with argv[0] = command_name.
 *
 * argv = pointer to pointer to char => 2D array of strings
 *        where argv[i][j] = jth char of the ith string
 */
main (int  argc, char **argv)
{
  
   /* Command line check values:
    * bopt = 0 if -B is not present, 1 if -B is present
    * topt = 0 if -T is not present, 1 if -T is present
    */
   int error, bopt, topt, vopt;

   /* Count stations read, output (out), and rejected (bad) */
   int staread, staout, stabad;

   /* Local counters */
   int  i, j, curfile = 1, nfiles = 0;

   char *dir, *extent, *st;
   FILE *shipfile;
   
/* Check for presence of command line arguments */
   if (argc < 2) {
     /* Only the command_name is at command line
      * so print usage and quit. */
     print_usage(argv[0]);
     exit(1);
   }
   
/* Set default values defined above */
   dir = DIR;  
   extent = EXTENT;

/* Nullify all data in the data storage structures */
   for (i = 0; i < MAXPROP; ++i) {
     data.observ[i] = NULL;
     bad_data.observ[i] = NULL;
   }
   hdr.prop_id = NULL;
   bad_hdr.prop_id = NULL;

/* Parse command line arguments */

/* Comment out this more general block and just read lon lat as inputs.
 * Removed the comment end lines so whole block could be commented out;
 * Will need to replace them if you want more functionality.
   for (i = 1; i < argc; i++) {
     For each command line argument, skipping argv[0] = command_name
 
     if (argv[i][0] == '-') {
       Then we have a switch, determine which type...
       switch (argv[i][1]) {

       case 'D':   /* Get input dir
	 dir = &argv[i][2];
	 break;

       case 'E':   /* Get file extent
	 extent = &argv[i][2];
	 break;

       case 'h':  
	 print_usage(argv[0]);
	 exit(0);

       case 'V':
	 vopt = 1;
	 break;

       default:
	 error = 1; 
       }    /* end switch

       if (error ) {
             fprintf(stderr,"\nError parsing command line args.\n");
             fprintf(stderr,"     in particular: '%s'\n", argv[i]);
             exit(1);
          }

     }  /* end if argument is a switch

     else  {
       /* Command line argument is not a switch, rather
	* it must be lon and lat Augment
	* number of input files.
       hdr.lon = atof(argv[i]);
       hdr.lat = atof(argv[i+1]);
       break;
     }
   }  /* end for each command line argument */

   hdr.lon = atof(argv[1]);
   hdr.lat = atof(argv[2]);

   fprintf(stdout,"%i",ms10(hdr.lat,hdr.lon,&hdr.ms1));
   
   exit(0);
  
} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program) {

  fprintf(stderr,"\nConverts a longitude latitude pair to the WMO");
  fprintf(stderr,"\nMarsden Square number.  Usage:");
  fprintf(stderr,"\nhb_msnum <lon> <lat>\n\n");
   return;
}
   
