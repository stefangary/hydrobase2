/*  hb_decimate.c

................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth  Curry
                             Woods Hole Oceanographic Inst
                             1993
                             updated to ANSI Mar 2000
................................................................................
____________________________________________________________________________
  USAGE:  

hb_decimate filename(_roots) -P<prs_int> [-D<dirname>] [-E<file_extent>]

  list of filenames MUST be first argument or input expected from stdin

 -P : pressure interval for output series.

 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum

____________________________________________________________________________
                                                   ____________________________________________________________________________
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

   /* global variables to store station ... */

struct HYDRO_DATA data, newdata;
struct HYDRO_HDR hdr;
int prsint;
int      iflag;

  /* prototypes for locally defined functions */
  
void print_usage(char *);
void get_hydro_data(int );


main (argc, argv)
int   argc;
char *argv[];
{
   int     index, nprops = 0;
   int      i, count;
   int     curfile = 1, nfiles = 0; 
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   char    *extent, *dir, *st;
   int     infile;


/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    prsint = 0;
    error = 0;
    iflag = 0;
 
/* initialize these ... */

   for (i = 0; i < MAXPROP; ++i) {
      data.observ[i] = (double *) NULL;
      newdata.observ[i] = (double *) NULL;
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

               case 'I':
                        
                       iflag = 1;
                       break;
               case 'P':
                        
                       error = (sscanf(&argv[i][2],"%d", &prsint)==1) ? 0:1;
                       break;

	       case 'h':  
	               print_usage(argv[0]);
	               exit(0);

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

   if ( !prsint) {
       fprintf(stderr,"\nYou must specify input file(s) and a prs_int.\n");
       exit(1);
   }
   if (!nfiles) {
       fprintf(stderr,"\nExpecting input from stdin ... ");
       infile = STDIN;
   }

/* loop for each input file */

   do {

     if (nfiles > 0) {
       infile = open_hydro_file(dir, argv[curfile], extent, print_msg);
       if (infile < 0)
          goto NEXTFILE;
     }

     
            /* read each file completely */

 
     get_hydro_data(infile);

NEXTFILE:

     if (nfiles > 0) 
         close(infile);

   } while (curfile++ < nfiles );

   fprintf(stderr,"End of hb_decimate.\n");
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
 fprintf(stderr,"\n%s reads each station in a HydroBase file and outputs a new file with observations at the specified pressure interval\n", program);

   fprintf(stderr,"\nUsage:  %s filename_root(s)  -P<prsint>    [-D<dirname>] [-E<file_extent>] [-h]", program);

   fprintf(stderr,"\n\n  List of filenames MUST be first argument");
   fprintf(stderr,"\n  or station input is expected from stdin.");
   fprintf(stderr,"\n    -P  : pressure interval for output array");
   fprintf(stderr,"\n   [-I] : interpolate between levels spaced farther than pressure interval  ");
   fprintf(stderr,"\n   [-D] : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"-h help...... prints this message. \n");

   fprintf(stderr,"\n\n");  
   return;
}

/*****************************************************************************/
void get_hydro_data(int file)

   /*  Reads each station in a HydroBase file and outputs observations at
      the specified pressure interval.    
   */
{
   int error, i, ii, j, npts, deltap;
/* read each station in file ... */

    while ((error = get_station(file, &hdr, &data)) == 0) {

     if (!iflag) {
     /* decimate the arrays according to specified prsint */
     
        if (hdr.nobs > 10) {
           deltap = (int) (data.observ[(int)PR][5] - data.observ[(int)PR][4]);
   	   if (deltap < 1)
	      deltap = 1;
           npts = prsint / deltap;  /*  get # of pts per interval*/
    
           if (prsint % deltap)
              fprintf(stderr,"WARNING: prsint (%d) requested is not an even multiple of the pressure sorted ctd file %d\n", prsint, deltap);
     
           if (npts == 0) {
              npts = 1;
           }
     
           j = 0;
           for (i = 0; i < hdr.nobs; i+=npts) {
              for (ii = 0; ii < MAXPROP; ++ii) {
               if (data.observ[ii] != NULL) 
                  data.observ[ii][j] = data.observ[ii][i];
             }
             ++j;
           }
     
           /* add bottom observation  ... */
     
           if ((i-npts) != (--hdr.nobs)) {
             for (ii = 0; ii < MAXPROP; ++ii) {
               if (data.observ[ii] != NULL) 
                  data.observ[ii][j] = data.observ[ii][hdr.nobs];
             }
             ++j;
           }
     
           hdr.nobs = data.nobs = j;
	
        } /* end if hdr.nobs > 10 */
	
          if (hdr.nobs > 0 )
            write_hydro_station(STDOUT, &hdr, &data);
	
    } /* end if !iflag */
    
    if (iflag) {
    
            for (ii = 0; ii < MAXPROP; ++ii) {
		    if  (data.observ[ii] != NULL )  
	                     newdata.observ[ii] = (double *) calloc(6000, sizeof (double));
            }
    
          /* load first observation level */
            j = 0;
	    npts = 0;
            for (ii = 0; ii < MAXPROP; ++ii) {
               if (data.observ[ii] != NULL) 
	           newdata.observ[ii][npts] = data.observ[ii][j];
	    }
	   ++npts;
	   
           for (j = 1; j < hdr.nobs; ++j) {
	   
              deltap = (int) (data.observ[(int)DE][j] - data.observ[(int)DE][j-1]);
	      if (deltap > prsint) {
	      
	         newdata.observ[(int)DE][npts] = data.observ[(int)DE][j-1] + prsint;
		 
	          while (newdata.observ[(int)DE][npts] < data.observ[(int)DE][j] ) {
                     for (ii = 0; ii < MAXPROP; ++ii) {
                        if ((ii != (int)DE) &&  (data.observ[ii] != NULL) )
	                   newdata.observ[ii][npts] = hb_linterp(newdata.observ[(int)DE][npts], &data.observ[(int)DE][j-1], &data.observ[ii][j-1], 2) ;
		      }  /* end for */
		     ++npts;
	             newdata.observ[(int)DE][npts] = newdata.observ[(int)DE][npts-1] + prsint;
		  } /* end while */
		} /* end if deltap */
		
		/*load jth level */
                for (ii = 0; ii < MAXPROP; ++ii) {
		       if  (data.observ[ii] != NULL )  {
	                     newdata.observ[ii][npts] = data.observ[ii][j];
		       }
		}
		 ++npts;

	   } /* end for j */
	   
      /* output the station */
	   hdr.nobs = npts;
	   newdata.nobs = npts;
	   newdata.nprops = hdr.nprops;	
    
          if (hdr.nobs > 0 )
            write_hydro_station(STDOUT, &hdr, &newdata);
	    
    } /* end if iflag */

     for (i = 0; i < MAXPROP; ++i) {
       if (data.observ[i] != NULL) {
          free(data.observ[i]);
          data.observ[i] = NULL;
       }
       
       if (newdata.observ[i] != NULL) {
          free(newdata.observ[i]);
          newdata.observ[i] = NULL;
       }
       
     }    
     free(hdr.prop_id); 
     hdr.prop_id = NULL;
       
   }  /* end while */

   if (error > 0)
      report_status(error, stderr);

   return;

} /* end get_hydro_data() */

