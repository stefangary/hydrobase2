/*  hb_countsta.c
***********************************************************************
*  Usage:  hb_countsta infile_list [-D<directory>] [-Eextent] [-h]
*
*  Reads HydroBase data files in the infile_list
*  and simply counts the number of stations and
*  writes the integer result to stdout.
*
*  The optional flags specify the root directory
*  (-D) and the extension (-E) of the files.
*
*  The presence of the optional -S flag will output
*  the number of scans (data points) while the
*  absence of this flag will output just the number
*  of stations.
*
*  NOT YET IMPLEMENTED!!!!!!
*  The optional -M flag will count the number of
*  stations (or scans if -S is present) that have
*  at least one missing value in them.
*
**************************************************************************
*/ 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hydrobase.h"

#define  PRINT_MSG   0 /* Set to zero to suppress printing. */ 
#define  MAXRANGES  200 
#define    DIR     ""
#define    EXTENT   ""

/* Prototypes for locally defined functions */

/* Print usage if user error. */
void print_usage(char *);

main(int argc, char **argv) {

  /* Local counter variables:
   * nstat = number of stations read (over all files)
   * nscan = number of scans read (over all files) */
  int    nstat, nscan;
  int    status, indx, indy;
  int    n, i, j, s_flag;

  /* nfiles = number of files on command line
   * curfile = current file number */
  int    nfiles, curfile;
  
  /* File ID of hydrobase format input file */
  int    infile;
  
  /* Pointers to strings for directory and
   * extension roots as well as a generic
   * string for reading command line. */
  char   *dir, *extent, *st;
  
  /* Input header and data storage */
  struct HYDRO_HDR hdr;   
  struct HYDRO_DATA data;

  /* Error flag for reading command line */
  int error;

  /* Are there command line arguments? */
  if (argc < 2) {
    print_usage(argv[0]);
    exit(1);
  }

  /* Set default values... */
  dir = DIR;
  extent = EXTENT;
  infile = 0;
  nstat = nscan = nfiles = 0; /* Initialize station, scan, and file counts */
  curfile = 1;
  nfiles = 0;
  s_flag = 0;
 
  /* Initialize input data storage */
  for (i = 0; i < MAXPROP; ++i) {
    data.observ[i] = (double *) NULL;
  }
  hdr.prop_id = (int *)NULL;

  /* Parse command line arguments... */
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
               case 'D':
		 /* Set the input file directory. */
		 dir = &argv[i][2];
		 break;

               case 'E':
		 /* Set the input file extension. */
		 extent = &argv[i][2];
		 break;

               case 'S':
		 /* We seek to count scans and not stations. */
		 s_flag = 1;
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
      /* If none of the switches apply, then the
       * command line argument must be one of the
       * list of infiles.  Augment the input file
       * counter. */
      ++nfiles;
    }
  }  /* end for */
  
  /* Check that the user has specified at least one infile. */
  if ( nfiles < 1 ) {
    print_usage(argv[0]);
    fprintf(stderr,"\n\nYou must specify Hydrobase input files.\n");
    exit(1);
  }

  /*-----------------------------------------------------*/

  /* For each file in list_of_infiles, */
  do {
    
    /* Open the current infile. */
    infile = open_hydro_file(dir, argv[curfile], extent, PRINT_MSG);
    if (infile  < 0)
      /* We do not have a valid infile.x */
      goto NEXTFILE;
    
    /* Loop for each station */
    while ((status = get_station(infile, &hdr, &data)) == 0) {
      /* get_station returns 0 for sucessful read. */
      /* Input data are in structures hdr and data. */

      /* Augment the number of stations opened. */
      ++nstat;

      /* Add number of observations in station (scans)
       * to the total number of scans over all stations. */
      nscan = nscan + hdr.nobs;

    } /* End of loop over each station. */
    
    report_status(status, stderr);
    if ( nfiles ) close(infile);

  NEXTFILE:
    ;
  } while (curfile++ < nfiles); /* End of loop over all input files */
  
  /* Print information about total number of
   * stations opened or scans read. */
  if ( s_flag == 0 )
    fprintf(stdout,"%d",nstat);

  if ( s_flag == 1 )
    fprintf(stdout,"%d",nscan);

}  /* End main */

/****************************************************************************/

void print_usage(char *program)
{
  fprintf(stderr,"\n --------------------------------------------------------");
  fprintf(stderr,"\n Usage: %s list_of_infile(s) [-D<directory>] [-Eextent] [-h]",program);
  fprintf(stderr,"\n");
  fprintf(stderr,"\n List of infiles should be first argument(s) ");
  fprintf(stderr,"\n   -D  : specifies dir for input files (default is ./) ");
  fprintf(stderr,"\n        ex: -D../data/ ");
  fprintf(stderr,"\n   -E  : specifies file extent for input files (default is no extent)");  
  fprintf(stderr,"\n        ex: -E.dat ");
  fprintf(stderr,"\n   -S  : output in units of scans (data points) while");
  fprintf(stderr,"\n         ommission of flag outputs number of stations.");
  fprintf(stderr,"\n  [-h] : help -- prints this message. ");
  fprintf(stderr,"\n");
  fprintf(stderr,"\n This program will compute the number of stations in");
  fprintf(stderr,"\n the given list of infiles (all in hydrobase format).");
  fprintf(stderr,"\n");
  fprintf(stderr,"\n --------------------------------------------------------");
  fprintf(stderr,"\n");
  return;
} /* End print_usage() */

/***************************************************************************/
