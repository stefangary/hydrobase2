/*  hb_coriolis.c

................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             1993
................................................................................
____________________________________________________________________________
  USAGE:  

coriolis latitude

computes coriolis parameter for a given latitude.
*/

#include <stdio.h>
#include <stdlib.h>   /*sfg added to deal with exit warnings*/
#include <math.h>
#include "hydrobase.h"

main (int argc, char **argv)
{
int error;
double lat;

/* are there command line arguments? */

   if (argc < 2 || argv[1][1] == 'h') {
      fprintf(stderr,"\n %s Returns Coriolis parameter for a given latitude", argv[0]);
      fprintf(stderr,"\n\nUsage: %s <latitude>\n", argv[0]);
      exit(1);
   }

   error = (sscanf(argv[1],"%lf", &lat) == 1) ? 0 : 1;          
   if (error ) {
       fprintf(stderr,"\nError parsing latitude.\n");
       exit(1);
   }
   
   
   fprintf(stdout, "  %.8lf  ",  hb_coriol(lat));           
   exit(0);

}   /* end main */

