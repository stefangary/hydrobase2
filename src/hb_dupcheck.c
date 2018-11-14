/* hb_dupcheck.c
 *********************************************************
 * This program looks at the header information in each
 * station in a hydrobase file and checks whether or not
 * it is duplicated by any other stations in that same
 * file.
 *
 * This is for use with the Hydrobase2 routines and is
 * based on hb_extract.c, by Ruth Curry.
 * If you have a large dataset, try to subdivide your
 * data set into discrete chunks in space or time with
 * hb_extract before using this routine.
 *
 * Stefan Gary, October, 2009.
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
#define    TINC     1
#define    LINC     0.01
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
  int h_flag=0, d_flag=0, e_flag=0;
  int l_flag=0, t_flag=0, o_flag=0;
  
  /* Minimum time spacing */
  int tinc;
  
  /* Minimum distance spacing */
  float linc;

  /* Keeping track of errors during
   * command line or hydrobase reads. */
  int status;
  int error = 0;

  /* Keep track of number of stations
   * in, out, and duplicated. */
  int nin = 0, nout = 0, ndup = 0;
  
  /* File ID of the duplicate files. */
  int dupfile;

  /* Declare hydrobase header/data storage
   * structures for input and output of data. */
  struct HYDRO_HDR hdr;
  struct HYDRO_DATA data;

  /* Declare an array of headers (an array
   * of structures) for storage of all header
   * information in the file. */
  struct HYDRO_HDR *hdrs;

  /* Declare an array of integers to store
   * whether a file is a duplicate or not. */
  int *dupflag;

  /* Strings for directory, extension,
   * and command line input. */
  char *dir, *extent, *s;

  /* Intermediate values */
  float datei, datej;
  float distij;

  /* Are there command line arguments? */
  if (argc < 2) {
    print_usage(argv[0]);
    exit(1);
  }
  
  /* Set default values */
  dir = DIR;
  extent = EXTENT;
  tinc = TINC;
  linc = LINC;
  
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
	
      case 'O':
	o_flag = 1;
	dupfile = create_hydro_file(&argv[i][2], NOCLOBBER);
	if (dupfile < 0) {
	  fprintf(stderr,"\nUnable to open %s for writing.", &argv[i][2]);
	  fprintf(stderr,"\nIt may exist already?\n");
	  exit(1);
	}
	break;
	
      case 'T':
	t_flag = 1;
	s = &argv[i][2];
	error = sscanf(s,"%d", &tinc) != 1;
	break;
	
      case 'L':
	l_flag = 1;
	s = &argv[i][2];
	error = sscanf(s,"%f", &linc) != 1;
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

  fprintf(stderr,"\nWe have a time increment of %d days.",tinc);
  fprintf(stderr,"\nWe have a space increment of %f deg.",linc);

  /* Space to store input data for a profile is
   * explicitly non-existant. */
  for (i = 0; i < MAXPROP; ++i) {
    data.observ[i] = (double *) NULL;
  }
  hdr.prop_id = (int *)NULL;
  
  fprintf(stderr,"\nCounting the number of stations...");
  /* Initialize the count of number of stations
   * and the file count. */
  nin = 0;
  curfile = 1;
  
  /* Loop for each input file */
  do {
    infile = open_hydro_file(dir, argv[curfile], extent, PRINT_MSG);
    if (infile  < 0) 
      goto NEXTFILE;
    
    /* Loop for each station */
     while ((status = get_station(infile, &hdr, &data)) == 0) { 
       
       /* Augment the number of stations */
       nin++;
       
     } /* End of loop over all stations. */
     
     report_status(status, stderr);
     close(infile);

  NEXTFILE:
     ;
  } while (curfile++ < nfiles);

  fprintf(stderr,"\nThere are %d input stations.",nin);

  /* Allocate memory for the header storage. */
  hdrs = (struct HYDRO_HDR *)calloc(nin, sizeof(struct HYDRO_HDR));

  /* Allocate memory for the duplicate flags
   * and initialize the flags to 0
   * (all are assumed not duplicates). */
  dupflag = (int *)calloc(nin, sizeof(int));
  for ( i = 0; i < nin; i++ ) {
    dupflag[i] = 0;
  }
  
  /* Copy header information to storage structure. */
  fprintf(stderr,"\nCopying header information to memory...");
  curfile = 1;
  i = 0;
  do {
    infile = open_hydro_file(dir, argv[curfile], extent, PRINT_MSG);
    if (infile  < 0) 
       goto NEXTFILE2;
    
    /* Loop for each station */
    while ((status = get_station(infile, &hdr, &data)) == 0) { 

      hdrs[i] = hdr;

      i++;

    } /* End of loop over each station */

    report_status(status, stderr);
    close(infile);
    
  NEXTFILE2:
    ;
  } while (curfile++ < nfiles);

  /* Check for duplicates in the header information. */
  fprintf(stderr,"\nChecking for duplicates...");

  /* For each station in the list, */
  for ( i = 0; i < nin; i++ ) {

    fprintf(stderr,"\rChecking station # %d of %d",(i+1),nin);

    /* Check that it is not already a duplicate. */
    if( dupflag[i] == 0 ) {
      /* This station is not already flagged.
       * Do nothing for flagged stations. */

      /* Cycle through all the other (non-dup)
       * stations to see if there are any
       * duplicates of this one. */
      for ( j = (i+1); j < nin; j++ ) {

	/* Check that this station is not
	 * a duplicate. */
	if( dupflag[j] == 0 ) {

	  /* Check duplication in time by
	   * computing the sumdate of each
	   * profile and the distance, in
	   * degrees, between the profiles. */
	  if ( (hdrs[i].year == hdrs[j].year) && (hdrs[i].month == hdrs[j].month) ) {

	    distij = sqrt( pow((hdrs[i].lon - hdrs[j].lon),2) + pow((hdrs[i].lat - hdrs[j].lat),2) );

	    if ( (abs(hdrs[i].day-hdrs[j].day) <= tinc) && (distij <= linc) ) {
	      /* We have isolated duplicate stations.
	       * Now, keep the station that has the most
	       * observations.*/
	      if ( hdrs[i].nobs >= hdrs[j].nobs ) {
		dupflag[j] = 1;
		/*fprintf(stderr,"\n=== Duplicate! %f < %f ===",distij,linc);*/
	      }
	      else {
		dupflag[i] = 1;
		/*fprintf(stderr,"\n=== Duplicate! %f < %f ===",distij,linc);*/
	      }
	    } /* End of duplicate assignment. */
	  } /* End of year and month check. */
	} /* End of inner dup check. */ 
      } /* End of cycling over other stations. */
    } /* End of initial dup check. */
  } /* End of cycling over all stations. */

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
	
	/* Write out the whole station. */
	error = write_hydro_station(STDOUT, &hdr, &data);
	
	/* If there was an error in the writing process, */
	if (error) {
	  fprintf(stderr,"\nError code write_hydro_station() = %d.\n", error);
	  exit(1);
	}
	
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
  
  fprintf(stderr,"\n%d stations read", nin);
  fprintf(stderr,"\n%d stations unique",nout);
  fprintf(stderr,"\n%d stations duplicates",ndup);
  fprintf(stderr,"\n\nEnd of hb_dupcheck.\n\n");
  exit(0);
  
} /* End main() */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\n%s checks for duplicate stations in HydroBase files.", program);
   fprintf(stderr,"\n");
   fprintf(stderr,"\nUsage: %s filename_root(s)", program);
   fprintf(stderr,"      [-D<dirname>] [-E<extent>] [-O<outcast_file>] ");
   fprintf(stderr,"      [-T<tinc>]    [-L<length>] [-h]");
   fprintf(stderr,"\n");
   fprintf(stderr,"\n    -D  : specifies dirname (default is ./) ");
   fprintf(stderr,"\n        ex: -D../data/ ");
   fprintf(stderr,"\n    -E  : specifies file extent (default is no extent)");  
   fprintf(stderr,"\n        ex: -E.dat ");
   fprintf(stderr,"\n    -O  : specifies file to which duplicate stations are written");  
   fprintf(stderr,"\n        ex: -OVE1982.sta");
   fprintf(stderr,"\n    -T  : specifies the minimum time separation between");
   fprintf(stderr,"\n          duplicates in days. [1 day]");
   fprintf(stderr,"\n        ex: -T1");
   fprintf(stderr,"\n    -L  : specifies the minimum distance separation");
   fprintf(stderr,"\n          duplicates in degrees.  [0.01 degrees]");
   fprintf(stderr,"\n    -h  help ... prints this message.");
   fprintf(stderr,"\n");
   fprintf(stderr,"\n Note that duplicates are defined by any two stations");
   fprintf(stderr,"\n that are within both the time and distance minimum");
   fprintf(stderr,"\n separations that are specified by -T and -L.");
   fprintf(stderr,"\n\n");
   return;
}
/**************************************************************************/
