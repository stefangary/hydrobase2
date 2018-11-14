/* hb_dupcheck.c
 *********************************************************
 * This program looks at the data in each file and
 * removes any scans for which any single value is
 * missing.
 *
 * Stefan Gary, Feb. 2015.
 **********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "hydrobase.h"

/* Default values for options. */
#define    DIR     ""
#define    EXTENT   ""
#define    PRINT_MSG 1

/* Prototypes for locally defined functions */

/* Print information if user error at command line. */
   void    print_usage(char *);

main (int argc, char **argv)
{

  /* File counting */
  int curfile = 1, nfiles = 0; 
  int i, j, infile;

  /* Keep track of which flags are
   * entered at command line. */
  int h_flag=0;
  
  /* Keeping track of errors during
   * command line or hydrobase reads. */
  int status;
  int error = 0;

  /* Keep track of number of stations
   * in, out, and duplicated. */
  int num_sta_read = 0;
  int num_scan_read = 0;
  int num_scan_w_missing = 0;
  int num_sta_w_missing = 0;

  /* Declare hydrobase header/data storage
   * structures for input and output of data. */
  struct HYDRO_HDR hdr_in;
  struct HYDRO_DATA data_in;
  struct HYDRO_HDR hdr_out;
  struct HYDRO_DATA data_out;

  /* Strings for directory, extension,
   * and command line input. */
  char *dir, *extent, *s;

  /* Are there command line arguments? */
  if (argc < 2) {
    print_usage(argv[0]);
    exit(1);
  }
  
  /* Set default values */
  dir = DIR;
  extent = EXTENT;
  
  /* Parse the command line arguments */
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      case 'D':
	d_flag = 1;
	dir = &argv[i][2];
	break;
	
      case 'E':
	e_flag = 1;
	extent = &argv[i][2];
	break;
	
      case 'h':
	h_flag = 1;
	print_usage(argv[0]);
	exit(0);
	break;
	
      default  :
	fprintf(stderr,"\nError parsing command line");
	fprintf(stderr,"\n in particular: %s\n", argv[i]);
	exit(1);
      }  /* end switch */
    }  /* end if */
    else  {
      ++nfiles;
    }
  }  /* end for */
  
  /* Check that the user has specified at
   * least one input file. */
  if ( !nfiles ) {
    print_usage(argv[0]);
    fprintf(stderr,"\n\nYou must specify at least one input file!\n");
    exit(1);
  }

  /* Space to store input data for a profile is
   * explicitly non-existant. */
  for (i = 0; i < MAXPROP; ++i) {
    data.observ[i] = (double *) NULL;
  }
  hdr.prop_id = (int *)NULL;
  
  /* Loop over each file listed on command line. */
  fprintf(stderr,"\nLooping over each station in each file...");
  curfile = 1;
  i = 0;
  do {
    infile = open_hydro_file(dir, argv[curfile], extent, PRINT_MSG);
    if (infile  < 0) 
       goto NEXTFILE;
    
    /* Loop for each station */
    while ((status = get_station(infile, &hdr_in, &data_in)) == 0) { 
      /* ...and load the station to memory during 
       * the process of the loop check.  Ideally, it
       * would be nice to be able to work scan-by-scan,
       * but since hdr.nobs needs to reflect the number
       * of scans BEFORE it is written, one must read
       * the whole input station first.*/

      /* Allocate space for a scan_missing_flag for
       * each scan so copying is easy, later. */

      /* Loop over each scan in the station and
       * count how many have no missing values.*/

      /* Allocate space for the output station.*/

      /* Copy scans with no missing values to the
       * output station.*/

      /* Write out a whole station. */
      error = write_hydro_station(STDOUT, &hdr_out, &data_out);
	
	/* If there was an error in the writing process, */
	if (error) {
	  fprintf(stderr,"\nError code write_hydro_station() = %d.\n", error);
	  exit(1);
	}

    } /* End of loop over each station */

    report_status(status, stderr);
    close(infile);
    
  NEXTFILE:
    ;
  } while (curfile++ < nfiles);

  /* Write stations to output, sorting
   * duplicates away from a clean single copy. */
  fprintf(stderr,"\nWriting to output files...");

  /* Initialize the file count. */
  curfile = 1;
  i = 0;

  /* Loop for each input file in the
   * same order as when they were first read. */
  do {
    infile = open_hydro_file(dir, argv[curfile], extent, PRINT_MSG);
    if (infile  < 0) 
       goto NEXTFILE3;
    
    /* Loop for each station */
    while ((status = get_station(infile, &hdr, &data)) == 0) { 
      
      /* If the station is part of the singles file, */
      if ( dupflag[i] == 0 ) {
	

	
	/* Augment the counter of output (kept) files. */
	++nout;
      }
      else {
	/* The station is a duplicate. */
	
	/* Augment the duplicate station counter. */
	++ndup;
	
	if (o_flag) {
	  /* The user wants to save duplicate data. */
	  error = write_hydro_station(dupfile, &hdr, &data);
	}
	
      } /* End of station OK check */
      
      /* Augment total station counter. */
      i++;

    }  /* End while loop over each station */
    
    report_status(status, stderr);
    close(infile);
    
  NEXTFILE3:
    ;
  } while (curfile++ < nfiles);
  
  fprintf(stderr,"\n%d stations read", num_sta_read);
  fprintf(stderr,"\n%d stations with missing values",num_sta_w_missing);
  fprintf(stderr,"\n%d total scans read",num_scan_read);
  fprintf(stderr,"\n%d total scans with missing values",num_scan_w_missing);
  fprintf(stderr,"\n\nEnd of hb_rmmissing.\n\n");
  exit(0);
  
} /* End main() */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\n%s removes scans with missing values from HB files.", program);
   fprintf(stderr,"\n");
   fprintf(stderr,"\nUsage: %s filename_root(s)", program);
   fprintf(stderr,"      [-D<dirname>] [-E<extent>] [-h]");
   fprintf(stderr,"\n");
   fprintf(stderr,"\n    -D  : specifies dirname (default is ./) ");
   fprintf(stderr,"\n        ex: -D../data/ ");
   fprintf(stderr,"\n    -E  : specifies file extent (default is no extent)");  
   fprintf(stderr,"\n        ex: -E.dat ");
   fprintf(stderr,"\n    -h  help ... prints this message.");
   fprintf(stderr,"\n");
   fprintf(stderr,"\n Output HB format files are sent to stdout.");
   fprintf(stderr,"\n Any scan for which at least one variable is");
   fprintf(stderr,"\n equal to HB_MISSING = %f is removed",HB_MISSING);
   fprintf(stderr,"\n One could extend this to search for values");
   fprintf(stderr,"\n in particular ranges and exclude them, but");
   fprintf(stderr,"\n that is not implemented here.");
   fprintf(stderr,"\n\n");
   return;
}
/**************************************************************************/
