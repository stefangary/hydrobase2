/*coastal_sta.c
................................................................................
                              *  HydroBase 2 *
................................................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             1993
			     Updated to ANSI-C Feb 2000
................................................................................
................................................................................
.  separates HydroBase2 station files into coastal stations (<200 m) and
.  blue water stations (>= 200 m) 
.
................................................................................
................................................................................

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hydrobase.h"

#define    PRINT_MSG 1

void    print_usage(char *);

main (int argc, char **argv)
{
   int     i, status, error;
   int     infile, cfile, ofile;
   int     c_flag, o_flag, target;
   int    nread, ncoast, nother;
   struct HYDRO_HDR h;
   struct HYDRO_DATA data;

/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
   

/*  set these default values */

   c_flag = o_flag = 0;
   ncoast = nother = nread = 0;
   target = 200;

/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      switch (argv[i][0]) {
         case '-':
             switch (argv[i][1]) {
                case 'C':
                      c_flag = 1;
                      cfile = create_hydro_file(&argv[i][2], NOCLOBBER);
                      if (cfile < 0) {
                         fprintf(stderr,"\nUnable to open %s for output.", &argv[i][2]);
                         fprintf(stderr,"\nDoes it exist already?\n");
                         exit(1);
                      }
                      fprintf(stderr,"\nOpened %s for output", &argv[i][2]);
                      break;
                case 'D':
                      error = sscanf(&argv[i][2],"%d", &target) != 1;
                      if (error)
                         fprintf(stderr,"\nError parsing argument: %s\n", argv[i]);
                      break;
                case 'O':
                      o_flag = 1;
                      ofile = create_hydro_file(&argv[i][2], NOCLOBBER);
                      if (ofile < 0) {
                         fprintf(stderr,"\nUnable to open %s for output.", &argv[i][2]);
                         fprintf(stderr,"\nDoes it exist already?\n");
                         exit(1);
                      }
                      fprintf(stderr,"\nOpened %s for output", &argv[i][2]);
                      break;
                case 'h':
                      print_usage(argv[0]);
                      exit(0);
                default:
                      fprintf(stderr,"\nError parsing command line");
                      fprintf(stderr,"\n in particular: %s\n", argv[i]);
                      exit(1);
             }  /* end switch */
             break;
         default:
              infile = open_hydro_file("", argv[i], "", PRINT_MSG);
              if (infile  < 0) {
                 fprintf(stderr,"\nError opening %s for input", argv[i]);
                 exit(1);
              }
      }  /* end switch */
   }  /* end for */


/* initialize some variables */

   for (i = 0; i < MAXPROP; ++i) {
     data.observ[i] = (double *) NULL;
   }
   h.prop_id = (int *)NULL;


     /* loop for each station */

   while ((status = get_station(infile, &h, &data)) == 0) { 
       ++nread;
   
       if (h.pdr <= 0)
          fprintf(stderr,"\nWARNING! pdr field is %d ",h.pdr);

       if (h.pdr < target) {
         error = write_hydro_station(cfile, &h, &data);
         ++ncoast;
       } 
       else {
         error = write_hydro_station(ofile, &h, &data);
         ++nother;
       } 
       
       if (error) {
             fprintf(stderr,"\nError in write_hydro_station(): error code = %d.\n", error);
             exit(1);
       }

   }  /* end while */

   report_status(status, stderr);
   close(infile);
   
   fprintf(stderr,"\n Total stations read in: %d ", nread);
   fprintf(stderr,"\n Number of coastal stations: %d", ncoast);
   fprintf(stderr,"\n Number of blue water stations: %d", nother);


  fprintf(stderr,"\nEnd of %s.\n", argv[0]);
  exit(0);

} /* end main() */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\n%s sorts stations by pdr depth.\n", program);

   fprintf(stderr,"\nUsage: %s infile", program);
   fprintf(stderr,"  -Ccoastal_file -Oother_file [-Dpdr_depth]\n");
   fprintf(stderr,"\n    -C  : specifies file to which stations shallower than -D<pdr_depth> are written");  
   fprintf(stderr,"\n        ex: -O7301.btl.shelf.raw");
   fprintf(stderr,"\n    -O  : specifies file to which deep stations are written");  
   fprintf(stderr,"\n        ex: -O7301.btl.deep.raw");
   fprintf(stderr,"\n   [-D] : specifies depth criterion ");
   fprintf(stderr,"\n          default is 200 meters ");
   fprintf(stderr,"\n\n");  
   return;
}

