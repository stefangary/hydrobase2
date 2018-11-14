/*  hb_getpos.c
................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             1993
			     Updated Feb 2000 to ANSI standards
................................................................................
................................................................................
.  Extracts the lat, lon from HydroBase data files and writes it out to 
.  the stdout device.
.
................................................................................
*/

#include <stdio.h>
#include <stdlib.h>   /*sfg added to deal with exit and malloc warnings*/
#include <fcntl.h>
#include <string.h>
#include "hydrobase.h"


#define    EXTENT    ""
#define    DIR       ""
#define    PRINT_MSG  1               /* 0 or 1 */

main (int argc, char **argv)
{
   int     i, nfiles, curfile; 
   int     error, status, infile;
   short   bopt, latfirst, output_stano, output_time;
   float  xmin, xmax, ymin, ymax;
   int     xdateline;
   char   *dir, *extent, *st;
   struct HYDRO_HDR hdr;
   struct HYDRO_DATA data;
   void print_usage(char *);



/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    bopt = 0;
    output_stano = 0;
    output_time = 0;
    latfirst = 0;
    error = 0;
    xdateline = 0;
    nfiles = 0;
    curfile = 1;


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

               case 'B':                    /* get grid bounds */
                        bopt = 1;
                        st = &argv[i][2];
                           if (*st == '/')
                               ++st;
                        error = (sscanf(st,"%f", &xmin) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &xmax) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &ymin) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &ymax) != 1);
                        if (xmin > 0 && xmax < 0)
                           xmax += 360;
                        if (xmax > 180)
                           xdateline = 1;
                        break;
		case 'S':
		        output_stano = 1;
			break;

	        case 'T':
		        output_time = 1;
			break;
               case ':':
                        latfirst = 1;
                        break;
               case 'h':
                        print_usage(argv[0]);
			exit(0);
                        break;
               default :
                        print_usage(argv[0]);
                        fprintf(stderr,"\nError parsing command line");
                        fprintf(stderr,"\n in particular: %s\n", argv[i]);
                        exit(1);
            }  /* end switch */

            if (error ) {
                print_usage(argv[0]);
                fprintf(stderr,"\nError parsing command line args.\n");
                fprintf(stderr,"     in particular: '%s'\n", argv[i]);
                exit(1);
            }
            
       }  /* end if */
       else  {
           ++nfiles;
       }
   }  /* end for */


/* initialize these... */

   for (i = 0; i < MAXPROP; ++i)
      data.observ[i] = (double *) NULL;
   hdr.prop_id = (int *) NULL;


/* loop for each file */

   do {
   if ( !nfiles) {
      infile = STDIN;
      fprintf(stderr,"\n Expecting data from stdin....  ");
   }
   else {
      infile = open_hydro_file(dir, argv[curfile], extent, PRINT_MSG);
      if (infile  < 0) 
       goto NEXTFILE;
   }
   
     /* loop for each station */

     while ((status = get_station(infile, &hdr, &data)) == 0) { 
 
       if (xdateline && hdr.lon < 0)
          hdr.lon += 360;
          
       if (bopt) {
          if ((hdr.lon <= xmax) && (hdr.lon >= xmin) && 
           (hdr.lat <= ymax) && (hdr.lat >= ymin) ) {
	     if (output_stano)
	        fprintf(stdout, "%4d ", hdr.station);
	     if (output_time)
	        fprintf(stdout, "%4d %4d %4d ",hdr.year,hdr.month,hdr.day);
	     if (latfirst)  
               fprintf(stdout, "%8.3f %8.3f\n", hdr.lat, hdr.lon);
	     else
               fprintf(stdout, "%8.3f %8.3f\n", hdr.lon, hdr.lat);
	  }     
       }
       else {
	     if (output_stano)
	        fprintf(stdout, "%4d ", hdr.station);
	     if (output_time)
	        fprintf(stdout, "%4d %4d %4d ",hdr.year,hdr.month,hdr.day);
	     if (latfirst)  
               fprintf(stdout, "%8.3f %8.3f\n", hdr.lat, hdr.lon);
	     else
               fprintf(stdout, "%8.3f %8.3f\n", hdr.lon, hdr.lat);
       }

     }  /* end while */ 

     report_status(status, stderr);
     close(infile);

NEXTFILE:
      ;
   } while (curfile++ < nfiles);

   fprintf(stderr,"\n\n End of hb_getpos. \n");
   exit(0);

}  /* end of main() */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nExtracts latitude/longitude from HydroBase data files");
   fprintf(stderr,"\n    and writes it to the stdout device.\n");
   fprintf(stderr,"\nUsage:  %s list_of_filenames", program);

   fprintf(stderr," [-Ddirname] [-Eextent] [-Bminlon/maxlon/minlat/maxlat] [-:]");
   fprintf(stderr,"\n If no infiles are specified, input is expected from stdin");
   fprintf(stderr,"\n    -D  : specifies directory of input files (default is ./) ");
   fprintf(stderr,"\n        ex: -D../data/ ");
   fprintf(stderr,"\n    -E  : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n        ex: -E.dat ");
   fprintf(stderr,"\n    -B  : specifies optional boundaries.");
   fprintf(stderr,"\n        ex: -B-90/0/0/65 ");
   fprintf(stderr,"\n    -S  : output station # in first column.");
   fprintf(stderr,"\n    -:  : outputs latitude/longitude. ");
   fprintf(stderr,"\n        default is lon/lat ");
   fprintf(stderr,"\n    -T  : appends time information (year, month, day).");
   fprintf(stderr,"\n\n");  
   return;
}

