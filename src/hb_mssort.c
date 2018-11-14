/* hb_mssort.c
................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             1993
                             Updated March 2000
................................................................................
............................................................................

Sorts file(s) of HydroBase stations by geographic squares. The size of square
is specified by the -S parameter.  Up to MAXF files can be simultaneously
opened. Any stations that do not fit into any of the open files are
stored in a file called msextra.dat -- which can be sorted subsequently.
A summary of numbers of stations in each file is written to stdout.
____________________________________________________________________________
*/
#include <stdio.h>
#include <stdlib.h>   /*sfg added to deal with exit and malloc warnings*/
#include <string.h>
#include <errno.h>
#include "hydrobase.h"

#define   MAXF     1000     /* max # output files. Orig = 50, commets say: "cannot exceed 60", but /usr/include/linux.h and /usr/include/fs.h allow for at least 1024. */
#define   PMODE    0666   /* read & write permission for output files */
#define   DIR     ""                 
#define   EXTENT   ""
#define   PRINT_MSG  1    /* 0 or 1 */
#define   NFOUT5    4    /* # of 5 degree squares in ten degrees */
#define   NFOUT1    100    /* # of 1 degree squares in ten degrees */
#define   NFOUT2    16    /* # of 2.5 degree squares in ten degrees */


int msq1to5[100] =  { 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,  /* converts ms_1 to ms_5 */
                   0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
                   0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
                   0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
                   0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
                   2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
                   2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
                   2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
                   2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
                   2, 2, 2, 2, 2, 3, 3, 3, 3, 3  };


char *fext5[NFOUT5] = { "_0", "_1", "_2", "_3"};

char *fext1[NFOUT1] = { "_00", "_01", "_02", "_03","_04",
                      "_05", "_06", "_07", "_08","_09",
                      "_10", "_11", "_12", "_13","_14",
                      "_15", "_16", "_17", "_18","_19",
                      "_20", "_21", "_22", "_23","_24",
                      "_25", "_26", "_27", "_28","_29",
                      "_30", "_31", "_32", "_33","_34",
                      "_35", "_36", "_37", "_38","_39",
                      "_40", "_41", "_42", "_43","_44",
                      "_45", "_46", "_47", "_48","_49",
                      "_50", "_51", "_52", "_53","_54",
                      "_55", "_56", "_57", "_58","_59",
                      "_60", "_61", "_62", "_63","_64",
                      "_65", "_66", "_67", "_68","_69",
                      "_70", "_71", "_72", "_73","_74",
                      "_75", "_76", "_77", "_78","_79",
                      "_80", "_81", "_82", "_83","_84",
                      "_85", "_86", "_87", "_88","_89",
                      "_90", "_91", "_92", "_93","_94",
                      "_95", "_96", "_97", "_98","_99"};

char *fext2[NFOUT2] = { "_0h", "_1h", "_2h", "_3h",
                      "_4h", "_5h", "_6h", "_7h",
                      "_8h", "_9h", "_Ah", "_Bh",
                      "_Ch", "_Dh", "_Eh", "_Fh"};

/* Store the MS# corresponding to each open file
 * as well as the scratch file. */
int msq1[MAXF + 1];
int msq10[MAXF + 1];

int count[MAXF + 1];
int outfile[MAXF + 1];          /* File descriptors returned by open()*/
char *msname[MAXF + 1];
int nextfile;                   /* Counter of number of currently open files. */

struct HYDRO_HDR hdr;
struct HYDRO_DATA data;

  /* prototypes for locally defined functions */

void print_usage(char *);  
int openfile(char *, char *, char *, int );
int sort10(char *, char *, int);
int sort5( char *, char *, int);
int getmsq2_5(float, float);
int sort2( char *, char *, int);
int sort1( char *, char *, int);


main (int argc, char **argv) {

  int j, error;
  int  infile, mode;
  int  nfiles, i;
  int  sortdim;
  int  curfile = 1,  status;
  char *fname, name[80];
  char *outpath, *extin, *extout, *dir;


  /* Check for command line arguments. */
  if (argc < 2) {
    print_usage(argv[0]);
    exit(1);
  }

  /* Set default values. */
  dir = DIR;
  outpath = DIR;
  extin = EXTENT;
  extout = EXTENT;
  mode = NOCLOBBER;
  nfiles = 0;
  sortdim = 0;
  nextfile = 0;
  error = 0;
  
  /* Parse command line arguments. */
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
        case 'D':
	  dir = &argv[i][2];
	  break;

        case 'E':
	  extin = &argv[i][2];
	  break;

        case 'O':
	  outpath = &argv[i][2];
	  break;

        case 'N':
	  extout = &argv[i][2];
	  break;

        case 'A':
	  mode = APPEND;
	  break;
          
        case 'T':
	  mode = OVERWRITE;
	  break;

        case 'S':
	  if (argv[i][2] == '\0') {
	    fprintf(stderr,"\nOptions for size of geographic square:");
	    fprintf(stderr,"\n 1   <1 deg sq>");
	    fprintf(stderr,"\n 2   <2.5 deg sq>");
	    fprintf(stderr,"\n 5   <5 deg sq>");
	    fprintf(stderr,"\n10   <10 deg sq>\n");
	    exit(0);
	  }
	  error = (sscanf(&argv[i][2],"%d", &sortdim) != 1);
	  break;

        case 'h': 
	  print_usage(argv[0]); 
	  exit(0);
   
        default :
	  error = 1;
                        
      }  /* End command line parsing switch. */
            
      if (error ) {
	fprintf(stderr,"\nError parsing command line args.\n");
	fprintf(stderr,"     in particular: '%s'\n", argv[i]);
	fprintf(stderr,"Use -h for a complete usage statement.\n");
	exit(1);
      }
    }  /* End check of command line error. */
    else  {
      /* This command line argument must be
       * an input file. */
      ++nfiles;
    }
  }  /* End for loop over all command line arguments. */

  if (sortdim == 0) {
    fprintf(stderr,"\nYou must specify the sorting parameter with -S\n");
    exit(1);
  }
   
      
  /* Initialize square counts. */
  for (j=0; j <= MAXF; ++j){
    count[j] = 0;
    msq1[j] = -1;
    msq10[j] = -1;
    outfile[j] = -1; 
  }

  /* Explicitly nullify all data/prop IDs. */
  for (i = 0; i < MAXPROP; ++i)
    data.observ[i] = (double *) NULL;
  hdr.prop_id = (int *) NULL;

  /* Open the scratch file.  First copy
   * the outpath to the scratch file name. */
  strcpy(name, outpath);

  /* Check for easily forgotten slash in path
   * and then add the default name of the
   * scratch file. */
  j = strlen(name);
  if (j && (name[j-1] != '/')) 
    strncat(name,"/",2);
  strncat(name, "msextra.dat", 12);
  
  /* Finally, create the scratch file (last
   * file in the outfile list. */
  if ((outfile[MAXF] = create_hydro_file(name,mode)) < 0) {
    switch (mode) {
      case APPEND:
	fprintf(stderr,"\nUnable to open %s in append mode.", name);
	exit(1);

      case OVERWRITE:
	fprintf(stderr,"\nUnable to open %s in overwrite mode.", name);
	exit(1);

      default:
	fprintf(stderr,"\nUnable to open %s in noclobber node.", name);
	fprintf(stderr,"  It may already exist...\n");
	exit(1);
    } /* End switch */
  }
      
  /* Loop for each input file. */
  do {
    if (nfiles == 0) {
      infile = STDIN;
      fprintf(stderr,"\nExpecting station input from STDIN ....\n");
    }
    else {
      infile = open_hydro_file(dir, argv[curfile], extin, PRINT_MSG);
      if (infile < 0) 
	goto NEXTFILE;
    }

    /* Loop to read each station.  First
     * read the station data into structures. */
    while ((status = get_station(infile, &hdr, &data)) == 0) {

      /* Sort according to requested square size.
       * Each of the sort routines returns the
       * index number of the file into which the
       * current station fits. */
      switch (sortdim) {
        case 1:
	  j = sort1(outpath, extout, mode);
	  break;

        case 2:
           j = sort2(outpath, extout, mode);
           break;

        case 5:
           j = sort5(outpath, extout, mode);
           break;

        case 10:
           j = sort10(outpath, extout, mode);
           break;

        default:
           fprintf(stderr,"\nInvalid size of sorting square. Options: 1,2,5,10\n");
           exit(1);      
     } /* End switch test for sorting square size. */
     

      /* Write station to appropriate outfile */
      write_hydro_station(outfile[j], &hdr, &data);

      /* Augment the number of stations in this file. */
      ++count[j];

    }  /* End while */
    report_status(status, stderr);

    close(infile);
  NEXTFILE:
    ;
  }  while (curfile++ < nfiles);


  /* Write summary of distribution to stdout and close output files. */
  fprintf(stdout, "\n\n Accounting of stations written to each file: \n");
  
  for (j=0; j < nextfile; ++j) {
    fprintf(stdout, "   %s%s     %d \n", msname[j], extout, count[j]);
    close(outfile[j]);
  }

  fprintf(stdout,"  msextra.dat    %d \n", count[MAXF]);
  close(outfile[MAXF]);                            /* close scratch file */
  
  exit(0);
} /* End of main */

/****************************************************************************/
void print_usage(char *program) {
  fprintf(stderr,"\n%s sorts file(s) of HydroBase stations by geographic squares.", program);
  fprintf(stderr,"\nThe size of square is specified by the -S parameter.");
  fprintf(stderr,"\nUp to MAXF files can be simultaneously opened.");
  fprintf(stderr,"\nAny stations that do not fit into an open file");
  fprintf(stderr,"\nare stored in a file called msextra.dat -- which can be");
  fprintf(stderr,"\nsorted subsequently. An accounting of numbers of stations ");
  fprintf(stderr,"\nin each file is written to stdout.\n");
  fprintf(stderr,"\nUsage:  %s list_of_filename_root(s)", program);

  fprintf(stderr," -S<size> [-Ddir] [-Eextent][-Ooutpath] [-Nnew_extent] [-A] [-T]");
  fprintf(stderr," > logfile\n");
  fprintf(stderr,"\n    -S  : size of geographic squares for sorting ");
  fprintf(stderr,"\n        options: 1 = 1 deg square ");
  fprintf(stderr,"\n                 2 = 2.5 deg  ");
  fprintf(stderr,"\n                 5 = 5 deg  ");
  fprintf(stderr,"\n                10 = 10 deg      ex: -S10 ");
  fprintf(stderr,"\n   [-D] : specifies directory of input files (default is ./) ");
  fprintf(stderr,"\n        ex: -D../data/ ");
  fprintf(stderr,"\n   [-E] : specifies input file extent (default is no extent)");  
  fprintf(stderr,"\n        ex: -E.dat ");
  fprintf(stderr,"\n   [-O] : specifies directory of output files (default is ./) ");
  fprintf(stderr,"\n        ex: -O../natl/ ");
  fprintf(stderr,"\n   [-N] : specifies output file extent (default is no extent)");  
  fprintf(stderr,"\n        ex: -N.dat ");
  fprintf(stderr,"\n   [-A] : append to existing files (default is do not alter an existing file.)");
  fprintf(stderr,"\n   [-T] : truncate existing files (default is do not alter an existing file.)");
  fprintf(stderr,"\n   [-h] : help... prints this message.)");
  fprintf(stderr,"\n\n");  
  return;
} /* End print_usage() */

/****************************************************************************/
int sort10(char * outpath, char *extout, int mode) {

  int j; /* Internal counter */
  
  /* Check if this ms file is already open
   * by cycling through all the file counters
   * and exiting the loop if we exceed the
   * number of currently open files or if
   * we reach an ms file that this profile
   * falls into. */
  j = 0;
  while ((j < nextfile) && (msq10[j] != hdr.ms10))
    ++j;

  /* If we exceeded the number of currently open
   * files while checking above, we open another
   * output file (as long as we are within the
   * total open file limits. */
  if ((j == nextfile) && (j < MAXF)) {

    /* Allocate space for the MS# */
    msname[j] = (char *) calloc(5, sizeof(char));

    /* Write the 4-digit MS# to that variable. */
    sprintf(msname[j],"%4d", hdr.ms10);

    /* Open that file. */
    outfile[j] = openfile( outpath, msname[j], extout, mode);

    /* Update that array that stores MS# for each open file. */
    msq10[j] = hdr.ms10;

    /* Augment counter of currently open files. */
    ++nextfile;
  } 
      
  return j;

} /* End sort10() */

/****************************************************************************/
int sort5(char * outpath, char *extout, int mode) {

  int j; /* Internal counter */

  /* Check if this ms file is already open
   * and at the same time, find out if this
   * station fits into an existing open file
   * which cooresponds to a sorting box. */
  j = 0;   
  while ((j < nextfile) && !((msq1[j] == msq1to5[hdr.ms1])&&(msq10[j] == hdr.ms10))) 
    ++j; 

  /* Create a new output file if we have
   * encountered a station that falls
   * outside of the currently opened files
   * that coorespond to sorting boxes. See
   * comments in sort10 for details. */
  if ((j == nextfile) && (j < MAXF)) {
    msq1[j] = msq1to5[hdr.ms1];
    msq10[j] = hdr.ms10;
    msname[j] = (char *) calloc(7, sizeof(char));
    sprintf(msname[j],"%4d%s", hdr.ms10, fext5[msq1[j]]);
    outfile[j] = openfile(outpath, msname[j], extout, mode);
    ++nextfile;
  } 
  return j;
} /* End sort5() */

/****************************************************************************/
int getmsq2_5(float lat, float lon) {
  /* Returns the 2.5 degree Marsden
   * Square designation given coordinates. */

   int k, kk;

   k =  (int)(ABS(lat + .00001) * 10) % 100;
   kk = (int)((ABS(lon + .00001) ) * 10) % 100;
   return ((k / 25) * 4 + (kk / 25));
}  /* end msq2_5() */

/****************************************************************************/
int sort2(char * outpath, char *extout, int mode) {

  int j;

  /* Check if this ms file is already open. */
  j = 0;   
  while ((j < nextfile) && !((msq1[j] == getmsq2_5(hdr.lat, hdr.lon))&&(msq10[j] == hdr.ms10))) 
    ++j; 

  /* If we didn't sort the current station
   * into an already open file, open a new
   * output file. */
  if ((j == nextfile) && (j < MAXF)) {
    msq1[j] = getmsq2_5(hdr.lat, hdr.lon);
    msq10[j] = hdr.ms10;
    msname[j] = (char *) calloc(8, sizeof(char));
    sprintf(msname[j],"%4d%s", hdr.ms10, fext2[msq1[j]]);
    outfile[j] = openfile(outpath, msname[j], extout, mode);
    ++nextfile;
  } 
  return j;
} /* End sort2() */

/****************************************************************************/
int sort1(char * outpath, char *extout, int mode) {

  int j;

  /* Check if this ms file is already open. */
  j = 0;   
  while ((j < nextfile) && !((msq1[j] == hdr.ms1)&&(msq10[j] == hdr.ms10))) 
    ++j; 

  /* open  an output file */
  if ((j == nextfile) && (j < MAXF)) {
    msq1[j] = hdr.ms1;
    msq10[j] = hdr.ms10;
    msname[j] = (char *) calloc(8, sizeof(char));
    sprintf(msname[j],"%4d%s", hdr.ms10, fext1[hdr.ms1]);
    outfile[j] = openfile(outpath, msname[j], extout, mode);
    ++nextfile;
  } 
  return j;
} /* End sort1() */

/****************************************************************************/
int openfile(char *outpath, char *n, char *extent, int mode) {

  /* Opens an existing file or creates a new file for output.
   * The files are named: <outpath><n>.<extent>.
   *
   * function arguments:
   *       n:
   * outpath:      null terminated strings 
   *  extent:
   *    mode:     APPEND, OVERWRITE or NOCLOBBER for output file 
   */

  char fname[200];
  int  i;    

  /* Move beyond any leading dot in file extent. */
  if (*extent == '.')
    ++extent;

  /* Check pathname for easily forgotten slash */
  i = strlen(outpath);
  if (i  && (outpath[i-1] != '/')) 
    sprintf(fname, "%s/%s.%s", outpath, n, extent);  
  else
    sprintf(fname, "%s%s.%s", outpath, n, extent);  
      

  if ((i = create_hydro_file(fname,mode)) < 0)  {
    switch (mode) {
      case APPEND:
	fprintf(stderr,"\nUnable to open %s in append mode.", fname);
	exit(1);

      case OVERWRITE:
	fprintf(stderr,"\nUnable to open %s.", fname);
	exit(1);

      default:
	fprintf(stderr,"\nUnable to open %s.", fname);
	fprintf(stderr,"  It may already exist...\n");
	exit(1);
    }
  }

  switch (mode) {
    case APPEND:
      fprintf(stderr,"\nOpening or appending %s.", fname);
      return(i);

    case OVERWRITE:
      fprintf(stderr,"\nOpening or overwriting %s.", fname);
      return(i);

    default:
      fprintf(stderr,"\nOpening %s.", fname);
      return(i);
  }

}  /* End of openfile */

/****************************************************************************/
