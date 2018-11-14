/*  hb_monochk.c

................................................................................
                          *******  HydroBase2 *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             1993
                             updated to ANSI C 1999
................................................................................
____________________________________________________________________________
  USAGE:  

 hb_saltshift filename(_roots) -S<offset> [-D<dirname>] [-E<file_extent>] 

  list of filenames MUST be first argument!

 -S : amount to shift salinity;

 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum

____________________________________________________________________________
Shifts the salinity by the amount specified with -S.                                                    ____________________________________________________________________________
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

    /* flags to indicate which properties are to be output and which are
       needed for computation */
      

   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;

int nadjust;

main (int argc, char **argv)
{
   int     index, nprops = 0;
   int     i;
   int     curfile = 1, nfiles = 0; 
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   char    *extent, *dir, *st;
   int     infile;
   void    print_usage(char *);
   void    get_hydro_data(int);


/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    error = 0;
    nadjust = 0;

/* initialize these ... */

   for (i = 0; i < MAXPROP; ++i) {
      station.observ[i] = (double *) NULL;
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

   
   infile = STDIN;

/* loop for each input file */

   do {

     if (nfiles) {
       infile = open_hydro_file(dir, argv[curfile], extent, print_msg);
       if (infile < 0)
          goto NEXTFILE;

     }
     else {
       fprintf(stderr, "\n Expecting input from stdin....");
     }
            /* read each file completely */

 
         get_hydro_data(infile);

NEXTFILE:

         close(infile);

   } while (curfile++ < nfiles );
   
   
   
   fprintf(stderr,"Number of observation depths adjusted:  %d\n", nadjust);

   fprintf(stderr,"\n\nEnd of hb_monochk\n");
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\n\n***************************************************");
   fprintf(stderr,"\n hb_monochk checks for monotonically increasing depths and adjusts by shifting decimal place.  If the resulting depths are still not monotonic, the (adjusted) station is written to stderr ");
   fprintf(stderr,"\n***************************************************");

                                                       
   fprintf(stderr,"\nUsage:  %s filename_root(s)    [-D<dirname>] [-E<file_extent>]", program);

   fprintf(stderr,"\n\n  List of filenames MUST be first argument or station file is");
   fprintf(stderr,"\n  is expected to be piped from stdin.");
  fprintf(stderr,"\n    [-D] : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");

   fprintf(stderr,"\n\n");  
   return;
}

/*****************************************************************************/
void get_hydro_data(int file)
   /*  Reads each station in a HydroBase file and computes property values
       at each standard level.  This module requires that the HydroBase
       file contains a minimum of pr, de, te, sa observations.  */
{
   int error, i;
   double dlast;

/* read each station in file ... */

    while ((error = get_station(file, &hdr, &station)) == 0) {
         dlast = 0.0;      
          for (i = 0; i < hdr.nobs; ++i) {
             if (station.observ[(int)DE][i] < (dlast - 2)  && station.observ[(int)DE][i] < 500. ) {
                station.observ[(int)DE][i]   *= 10.0;
                station.observ[(int)PR][i]   *= 10.0;
                ++nadjust;
             }
             dlast = station.observ[(int)DE][i];   
          }
        
        /* now check */
        
           dlast = 0.0;
           error = 0;
           for (i = 0; i < hdr.nobs; ++i) {
             if (station.observ[(int)DE][i] < (dlast-2)) 
                  ++error;
           }
           
       if (error)
          write_hydro_station(STDERR, &hdr, &station);
       else
          write_hydro_station(STDOUT, &hdr, &station);

       for (i = 0; i < MAXPROP; ++i) {
          if (station.observ[i] != NULL) {
             free((void *)station.observ[i]);
             station.observ[i] = NULL;
          }
       }    
   }  /* end while */

   if (error > 0)
      report_status(error, stderr);

   return;

} /* end get_hydro_data() */

