/*  hb_countsta.c
***********************************************************************
*  Usage:  hb_thicksta infile_list [-D<directory>] [-Eextent] [-h]
*                                  [-T<tmin>] [-F<tmin>] [-O<outcast>]
*  Reads HydroBase data files in the infile_list,
*  determines the height of each station (by subtracting
*  the minimum observed value from the maximum observed
*  value) and writes the x,y,thickness triplet to
*  standard output.
*
*  The optional flags specify the root directory
*  (-D) and the extension (-E) of the files.
*
*  The presence of the optional -T<tmin> flag will output
*  only the thickness of the thickest continuous
*  interval.  The absence of this flag will return
*  the whole thickness of all scans.
*
*  The presence of the optional -F<tmin> flag will filter
*  out all stations that have a data gap of at least
*  the size tmin.  The -O<outcast_file_name> option is only
*  valid with the -F flag and allows for all removed files
*  to be written to a separate output file.
*
*  The presence of the optional -A<tmin> flag will compute
*  the thickness at each station ignoring any breaks in
*  the data greater than tmin.  This is similar to the
*  default mode of operation of the program (which simply
*  subtracts the shallowest depth from the deepest depth)
*  except than any blanks are removed from the total.
*
*  NOTE: The output with -F is a hydrobase file while the
*  output with the non-F operations is a simply .xyz file.
*  Both outputs go to stdout.
*
*  Stefan Gary, April 2013, based on hb_extract by Ruth Curry.
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
  int    n, i, j, t_flag, o_flag, f_flag, a_flag;

  /* nfiles = number of files on command line
   * curfile = current file number 
   * outcastfile = id of outcast file */
  int    nfiles, curfile, outcastfile;
  
  /* File ID of hydrobase format input file */
  int    infile;
  
  /* Pointers to strings for directory and
   * extension roots as well as a generic
   * string for reading command line. */
  char   *dir, *extent, *sss;
  
  /* Temporary values for thickness accumulators */
  double thickest, thick, netthick;
  float tmin;

  /* Input header and data storage */
  struct HYDRO_HDR hdr;   
  struct HYDRO_DATA data;

  /* Error flag for reading command line */
  int error, staOK;

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
  t_flag = 0;
  f_flag = 0;
  a_flag = 0;
  tmin = 50;

  /* Initialize input data storage */
  for (i = 0; i < MAXPROP; ++i) {
    data.observ[i] = (double *) NULL;
  }
  hdr.prop_id = (int *)NULL;

  /* Parse command line arguments... */
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
               case 'A':
		 /* We seek to filter out any blanks from
		  * the total thickness calculation. */
		 if ( t_flag == 1 ) {
		   fprintf(stderr,"\nCannot specify both -A and -T.");
		   exit(1);
		 }
		 if ( f_flag == 1 ) {
		   fprintf(stderr,"\nCannot specify both -A and -F.");
		   exit(1);
		 }
		 a_flag = 1;
		 error = (sscanf(&argv[i][2],"%f", &tmin) != 1);
		 break;

               case 'D':
		 /* Set the input file directory. */
		 dir = &argv[i][2];
		 break;

               case 'E':
		 /* Set the input file extension. */
		 extent = &argv[i][2];
		 break;

               case 'F':
		 /* We seek to filter out all stations that have a data gap. */
		 if ( t_flag == 1 ) {
		   fprintf(stderr,"\nCannot specify both -F and -T.");
		   exit(1);
		 }
		 if ( a_flag == 1 ) {
		   fprintf(stderr,"\nCannot specify both -F and -A.");
		   exit(1);
		 }
		 f_flag = 1;
		 error = (sscanf(&argv[i][2],"%f", &tmin) != 1);
		 break;

               case 'O':
		 o_flag = 1;
		 outcastfile = create_hydro_file(&argv[i][2], NOCLOBBER);
		 if (outcastfile < 0) {
		   fprintf(stderr,"\nUnable to open %s for writing.", &argv[i][2]);
		   fprintf(stderr,"\nIt may exist already?\n");
		   exit(1);
		 }
		 break;

               case 'T':
		 /* We seek to calculate continuous thickness and not total thickness. */
		 if ( f_flag == 1 ) {
		   fprintf(stderr,"\nCannot specify both -T and -F.");
		   exit(1);
		 }
		 if ( a_flag == 1 ) {
		   fprintf(stderr,"\nCannot specify both -T and -A.");
		   exit(1);
		 }
		 t_flag = 1;
		 error = (sscanf(&argv[i][2],"%f", &tmin) != 1);
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

  /* Check that if an outcast file is requested
   * the the user is also running in -F mode. */
  if ( o_flag == 1 && f_flag != 1 ) {
    print_usage(argv[0]);
    fprintf(stderr,"\n\nYou cannot specify -O without -F.\n");
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

      /* Check that there are at least 2 
       * observations in this station. 
       * Otherwise, thickness is meaningless.*/
      if ( hdr.nobs > 1 ) {

	if ( t_flag == 0 && f_flag == 0 && a_flag == 0 ) {
	  /* Compute the total thickness of the
	   * station by subtracting the depth
	   * of the upper measurement from the
	   * depth of the bottom measurement.
	   * No need to loop over each scan.
	   * Do not print out a station if only
	   * one scan is present. */
	  fprintf(stdout,"%8.3f %8.3f %9.3f\n",hdr.lon,hdr.lat,data.observ[(int)DE][(hdr.nobs-1)] - data.observ[(int)DE][0]);
	} else {
	  /* We care about the thickness of the
	   * thickest continuous swath or we need
	   * to at least filter out data blanks
	   * when computing thickness. */

	  /* Initialize the thickest segment */
	  thickest = 0.0;  /* Thickest swath encounters so far */
	  thick = 0.0;     /* Thickness of current swath */
	  netthick = 0.0;  /* Thickness of column ignoring blanks */
	  staOK = 1;

	  /* Loop over each scan, starting 
	   * with the second scan (because we
	   * want to check if the second scan
	   * is continuous after the first 
	   * scan */
	  for (j = 1; j < hdr.nobs; ++j) {
	    /* Check if the sepration between
	     * this point and the next point is
	     * greater than 20 m. */
	    if ( (data.observ[(int)DE][j] - data.observ[(int)DE][(j-1)]) > tmin ) {
	      /* We have reached a discontinuity */
	      staOK = 0;

	      /* Test if the accumulated thickness 
	       * exceeds the thickest detected
	       * layer so far.*/
	      if ( thick > thickest ) {
		thickest = thick;
	      }
		
	      /* Reset the current thickness accumulator*/
	      thick = 0.0;
	    }
	    else {
	      /* This is a continuous swath of data 
	       * so augment the thickness accumulators. */
	      thick = thick + data.observ[(int)DE][j]-data.observ[(int)DE][(j-1)];
	      netthick = netthick + data.observ[(int)DE][j]-data.observ[(int)DE][(j-1)];
	    }
	  } /* End of loop over each scan */

	  /* One last check to see if the last
	   * continuous swath was the deepest one.*/
	  if ( thick > thickest ) {
	    thickest = thick;
	  }

	  /* Write to output */
	  if ( t_flag == 1 ) {
	    /* In thickest mode, we only output an .xyz file */
	    fprintf(stdout,"%8.3f %8.3f %9.3f\n",hdr.lon,hdr.lat,thickest);
	  }
	  if ( a_flag == 1 ) {
	    /* In net thickest mode, we only output an .xyz file */
	    fprintf(stdout,"%8.3f %8.3f %9.3f\n",hdr.lon,hdr.lat,netthick);
	  }
	  if ( f_flag == 1 ) {
	    /* In filter mode, we output a hydrobase file 
	     * and, if requested, the outcast hydrobase file. */
	    if ( staOK == 1 ) {
	      error = write_hydro_station(STDOUT, &hdr, &data);
	      /* If there was an error in the writing process, */
	      if (error) {
		fprintf(stderr,"\nError in write_hydro_station(): error code = %d.\n", error);
		exit(1);
	      }
	    } else {
	      if ( o_flag == 1 ) {
		/* The user wants the whole outcast station. */
		error = write_hydro_station(outcastfile, &hdr, &data);
		/* If there was an error in the writing process, */
		if (error) {
		  fprintf(stderr,"\nError in write_hydro_station(): error code = %d.\n", error);
		  exit(1);
		}
	      }
	    } /* End of outcastfile check */
	  } /* End of data format output check */
	} /* End of -T or -F check */
      } 
      if ( hdr.nobs <= 1 && o_flag == 1) {
	/* Less than 2 observations, so only write this
	 * station to outcast file, if requested. */
	error = write_hydro_station(outcastfile, &hdr, &data);
	/* If there was an error in the writing process, */
	if (error) {
	  fprintf(stderr,"\nError in write_hydro_station(): error code = %d.\n", error);
	  exit(1);
	}
      } /* End of check for 2 obs */
    } /* End of loop over each station */
    
    report_status(status, stderr);
    if ( nfiles ) close(infile);

  NEXTFILE:
    ;
  } while (curfile++ < nfiles); /* End of loop over all input files */
  
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
  fprintf(stderr,"\n   -T  : output only the thickest continuous unit rather");
  fprintf(stderr,"\n         than the thickness encompassing all scans.  Any");
  fprintf(stderr,"\n         break in the data greater than 20m is considered");
  fprintf(stderr,"\n         a discontinuous interval.");
  fprintf(stderr,"\n  [-h] : help -- prints this message. ");
  fprintf(stderr,"\n");
  fprintf(stderr,"\n This program will compute the total height of the");
  fprintf(stderr,"\n scans in each station in the given list of infiles");
  fprintf(stderr,"\n (all in hydrobase format).  This is useful in tracking");
  fprintf(stderr,"\n the position of watermasses in conjunction with the");
  fprintf(stderr,"\n tools hb_siftsta and hb_siftlevs.");
  fprintf(stderr,"\n");
  fprintf(stderr,"\n --------------------------------------------------------");
  fprintf(stderr,"\n");
  return;
} /* End print_usage() */

/***************************************************************************/
