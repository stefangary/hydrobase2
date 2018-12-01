/*  wod13_convert.c
----------------------------------------------------------------------------------
*   Reads World Ocean Database 2013/18 observed level hydrographic data files
*   extracts :    header info
*                 p,d,t,s,ox,n2,n3,si,p4,f1,f2,f3,he,tu
*   
*   USAGE: wod13_convert filelist -S<shipfile> -T<instr_type> [-B<badfile>] [-D<dirname>] [-E<extent>] [-h] [-O<outfile>]
----------------------------------------------------------------------------------
* Originally based on wod01_convert.c, part of the original
* Hydrobase2 distribution.  Additional comments and changes
* added by Stefan Gary based on the Aug. 06, 2015
* version of wodC.c from WOD website.
*
* Also, looked at the wod_convert.c distributed in HB3 and
* included the ITS-68/90 temperature scale detection here
* as well.  The wod_convert.c version in HB3 is not updated
* to WOD13 because still uses the old 40 std depth levels
* while WOD13 uses 137 std depth levels.  For compatibility
* whith HB2 (which does not use ITS-90) and with HB3 (which
* can use either one), ALL OUTPUT WILL BE IN ITS-68.
*
* Tested with a subset of WOD18 native data and appears to
* work just fine.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hydrobase.h"

#define   BLANK     0x20   /* ascii code for blank */
#define   MISSING   -9.9   /* missing value flag */
#define   DIR    ""
#define   EXTENT ""

/*******************************************************************

 Definitions and global vars from NODC's wodC.c program

********************************************************************/

/*******************************************************************

 maxtax:  Max number of taxa delimitors in a taxa set
 maxsec:  Max number of secondary headers
 maxbio:  Max nymber of biological headers
 maxparm: Max number of measured and/or calculated variables
 maxpsec: Max number of variable specific secondary headers

********************************************************************/

#define maxtax 30
#define maxsec 100
#define maxbio 50
#define maxparm 100
#define maxpsec 25 * maxparm

/**********************************************************

 *FP - WOD DATA FILE TO BE OPENED FOR READING

***********************************************************/

FILE *fp;

/**********************************************************

 Header information for each cast

 cc        - Country code of country responsible for cast
 icruise   - NODC/WOD assigned cruise number (unique w/ cc)
 ostation  - WOD unique cast number
 year      - Year of cast
 month     - Month of cast
 day       - Day of cast
 hour      - Hour of cast
 latitude  - Position of cast
 longitude - Position of cast
 levels    - Number of depth levels measured at cast
 isoor     - Set to zero for observed level data, one for std. levels
 nparm     - number of parameters (variables) in cast 
 ip2       - Variable codes for variables measured at cast
 iperror   - Error flags for each variable measured at cast

 htotfig   - Total number of figures for 1: hour, 2: latitude, 3: longitude
 hsigfig   - Number of significant figures for 1: hour, 2: latitude, 3: longitude
 hrightfig - Number of figures to the right of the decimal for
             1: hour, 2: latitude, 3: longitude

**************************************************************/

char cc[2];
int icruise=0, ostation=0, year=0, month=0, day=0;
int hour, longitude, latitude;
int levels, isoor, nparm, ip2[maxparm], iperror[maxparm];
int htotfig[3], hsigfig[3], hrightfig[3];

/*************************************************************

 Character and primary investigator data

 origcfig - Number of characters in originators cruise code
            (ZERO IF NOT GIVEN)
 origc    - Originators cruise code
 origsfig - Number of characters in orginators station code
            (ZERO IF NOT GIVEN)
 origs    - Originators station code

 npi      - Number of primary investigators
 ipip     - Variable for PI
 ipi      - PI

*****************************************************************/

int origcfig, origsfig;
char origc[30], origs[30];
int ipip[maxparm], ipi[maxparm], npi;

/*************************************************************

 Secondary header information

 nsec      - Number of secondary headers
 stotfig   - Total figures for each secondary header
 ssigfig   - Number of significant figures for each secondary header
 srightfig - Number of figures right of decimal for each sec. header
 seccode   - Secondary header code
 secval    - Secondary header value

**************************************************************/

int nsec;
int stotfig[maxsec], ssigfig[maxsec], srightfig[maxsec];
int seccode[maxsec], secval[maxsec];

/*************************************************************

 Variable specific secondary header information

 npsec      - Number of var. specific sec. headers
 pstotfig   - Total figures for each var. specific sec. header
 pssigfig   - Number of sig. figs. for each var. specific sec. header
 psrightfig - Num figs. right of decimal foreach var. specific sec. header
 psecparm   - Value of var. specific sec. header
 pseccode   - Variable specific secondary header code
 psecval    - Variable specific secondary header value

**************************************************************/

int npsec;
int pstotfig[maxpsec],pssigfig[maxpsec],psrightfig[maxpsec];
int psecparm[maxpsec],pseccode[maxpsec],psecval[maxpsec];

/*************************************************************

 Biological header information

 nbio      - Number of biological headers
 btotfig   - Total figures for each biological header
 bsigfig   - Number of significant figures for each biological header
 brightfig - Number of figures right of decimal for each bio. header
 biocode   - Biological header code
 bioval    - Biological header value

**************************************************************/

int nbio;
int btotfig[maxbio], bsigfig[maxbio], brightfig[maxbio];
int biocode[maxbio], bioval[maxbio];

/*************************************************************

 Taxa set information

 ntsets     - Number of taxa sets at cast
 ntloc      - Number of delimitors for each taxa set
 ntcode     - Taxa code for each taxa variable
 ntval      - Value for each taxa code
 nterr      - Error flag for each taxa value
 ntoerr     - Originators flag for each taxa value
 nttotfig   - Total figures for each taxa value
 ntsigfig   - Number of significant fig. for each taxa value
 ntrightfig - Number of fig. right of the decimal for each taxa value

***************************************************************/

int ntsets;
int *ntloc,*ntcode,*ntval,*nterr,*ntoerr,*nttotfig,*ntsigfig,*ntrightfig;

/*************************************************************

 Depth information

 depth     - Depth at each depth level
 zerr      - Depth error flag
 zoerr     - Depth originators flag
 ztotfig   - Total figures for each depth value
 zsigfig   - Number of sig. figs. for each depth
 zrightfig - Number of figs right of the decimal for each depth

***************************************************************/

int *depth,*zerr,*zoerr,*ztotfig,*zsigfig,*zrightfig;

/*************************************************************

 Measured variable information

 dval      - Data value at each depth level
 derr      - Error flag for each variable at each level
 doerr     - Originators flag for each variable at each level
 dtotfig   - Total figures for each data value
 dsigfig   - Number of significant figures for each data value
 drightfig - Number of figs right of the decimal for each data value

***************************************************************/

int *dataval,*derr,*doerr,*dtotfig,*dsigfig,*drightfig;

/*************************************************************

 Size delimitors

 isize     - Present array size needed to fit all measured variables
 zsize     - Present array size needed to fit all depths
 ntsetmax  - Max number of taxa sets
 isizemax  - Max array size yet encountered for measured variables
 zsizemax  - Max array size yet encountered for depths

*************************************************************/

int isize, zsize;
int ntsetsmax=0, isizemax=0, zsizemax=0;
int iVERSflag;   /* WOD13 addition */

/*************************************************************

 Internal arrays

 tenp   - Powers of ten (double array in wodC.c, 9 entries)
 sdepth - Standard depth levels - only 40 levels, so this
          variable applies only to WOD09 and earlier.

*************************************************************/

/* WOD13 added 1x10^8 */
float tenp[] = { 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000,
		 100000000 };

/* WOD13 added periods at the end of each number.*/
int sdepth[] = { 0., 10., 20., 30., 50., 75., 100., 125., 150., 200., 250.,
                  300., 400., 500., 600., 700., 800., 900., 1000., 1100.,
                  1200., 1300., 1400., 1500., 1750., 2000., 2500., 3000.,
                  3500., 4000., 4500., 5000., 5500., 6000., 6500., 7000.,
                  7500., 8000., 8500., 9000. };

/*******************************************************
 * Up to here, matches all the variable declarations
 * in the WOD13 version of wodC.c.
 *******************************************************/

/*************************************************************

 Variables for reading ship codes from external files

 maxscode  - Max number of ship codes (lines) read from file
 nodc_ship - Ship name from the shipcode file.

*************************************************************/

int maxscode;
char **nodc_ship;   /* shipcode table indexed by OCL shipcode */

/************************************************************

 HydroBase globally defined variables

 The structure HYDRO_HDR contains:
       char country[3], ship[3]
       char origin, instrument
       int  cruise, station, year, month, day
       float lat, lon
       int  pdr, nobs, nprops, ms10, ms1
       int *prop_id
       char qual[nqual]

The structure HYDRO_DATA contains:
       int nobs, nprops
       double *observ[MAXPROP]

*************************************************************/
  
struct HYDRO_HDR  hdr, bad_hdr;
struct HYDRO_DATA data, bad_data;

/************************************************************

 Prototypes for locally defined functions

*************************************************************/

void print_usage(char *);
FILE *openfile(char *, char *, char *);
int check_sta();
int check_profile_flag(int);
int check_obs_flag(int);
int oclread();
void spacer(int);
int extractc(int, int *, char * );
int extracti(int, int *, int *, int *,  int *, int);
int nocrfgetc();

/***********************************************************/

/* Getting command line values at execution:
 * argc = number of parsed entries on command line
 *        with argv[0] = command_name.
 *
 * argv = pointer to pointer to char => 2D array of strings
 *        where argv[i][j] = jth char of the ith string
 */
main (int  argc, char **argv)
{
  
   /* Command line check values:
    * bopt = 0 if -B is not present, 1 if -B is present
    * topt = 0 if -T is not present, 1 if -T is present
    */
   int error, bopt, topt, vopt;

   /* File ID of output file and rejected values file */
   int outfile, badfile;

   /* Count stations read, output (out), and rejected (bad) */
   int staread, staout, stabad;

   /* Local counters */
   int  i, j, curfile = 1, nfiles = 0;

   /* Buffer to read in shipcode data */ 
   char buf[200];

   /* Buffer to store shipindex */
   char shipindex[5];

   /* Buffer to store shipcode */
   char shipcode[5];

   /* Buffer to store shipname */
   char shipname[190];

   char *dir, *extent, *st;
   FILE *shipfile;
   
/* Check for presence of command line arguments */
   if (argc < 2) {
     /* Only the command_name is at command line
      * so print usage and quit. */
     print_usage(argv[0]);
     exit(1);
   }
   
/* Set default values defined above */
   dir = DIR;  
   extent = EXTENT;

/* Set more default values */ 
   error = 0;
   staout = staread = stabad = 0; /* initialize station counts */
   outfile = STDOUT;

/* Set check values */
   topt = bopt = vopt = 0;
   shipfile = NULL;
   maxscode = 12000;

/* Nullify all data in the data storage structures */
   for (i = 0; i < MAXPROP; ++i) {
     data.observ[i] = NULL;
     bad_data.observ[i] = NULL;
   }
   hdr.prop_id = NULL;
   bad_hdr.prop_id = NULL;

/* Parse command line arguments */

   for (i = 1; i < argc; i++) {
     /* For each command line argument, skipping argv[0] = command_name */
     if (argv[i][0] == '-') {
       /* Then we have a switch, determine which type... */
       switch (argv[i][1]) {

       case 'B':   /* Open file for rejected data */

	 badfile = create_hydro_file(&argv[i][2], NOCLOBBER);
	 if (badfile < 1) {
	   fprintf(stderr,"\nError opening output file: %s\n", &argv[i][2]);
	   fprintf(stderr,"Does file already exist?\n"); 
	   exit(1);
	 }

	 fprintf(stderr,"\nBad data will be written to: %s", &argv[i][2]);
	 bopt = 1;
	 break;

       case 'D':   /* Get input dir */
	 dir = &argv[i][2];
	 break;

       case 'E':   /* Get file extent */
	 extent = &argv[i][2];
	 break;

       case 'O':   /* Open file for converted data  */
	 outfile = create_hydro_file(&argv[i][2], OVERWRITE);
	 if (outfile < 1) {
	   fprintf(stderr,"\nError opening output file: %s\n", &argv[i][2]);
	   exit(1);
	 }
	 fprintf(stderr,"\nOpening or overwriting: %s", &argv[i][2]);
	 break;
         
       case 'S':   /* Open shipcode file  */
	 shipfile = fopen(&argv[i][2], "r");
	 if (shipfile == NULL) {
	   fprintf(stderr,"\nError opening ship file: %s\n", &argv[i][2]);
	   exit(1);
	 }
	 break;

       case 'T':   /* Get datatype  */
	 switch (argv[i][2]) {
	 case 'c':  /* fall through */
	 case 'b':
	 case 'f':
	 case 's':
	 case 'u':
	 case 'm':
	   hdr.instrument = argv[i][2];
	   bad_hdr.instrument = argv[i][2];
	   break;
	 default:
	   fprintf(stderr,"\nError in datatype specification: %s\n", argv[i]);
	   fprintf(stderr,"Legitimate choices: \n");
	   fprintf(stderr,"   [c]td  \n");
	   fprintf(stderr,"   [b]ottle \n");
	   fprintf(stderr,"   [f]loat  \n");
	   fprintf(stderr,"   [s]easoar  \n");
	   fprintf(stderr,"   [m]oored profiler  \n");
	   fprintf(stderr,"   [u]nknown  \n\n");
	   exit(1);
	 }
	 topt = 1;
	 break;

       case 'h':  
	 print_usage(argv[0]);
	 exit(0);

       case 'V':
	 vopt = 1;
	 break;

       default:
	 error = 1; 
       }    /* end switch */

       if (error ) {
             fprintf(stderr,"\nError parsing command line args.\n");
             fprintf(stderr,"     in particular: '%s'\n", argv[i]);
             exit(1);
          }

       }  /* end if argument is a switch */

       else  {
	 /* Command line argument is not a switch, rather
	  * it must be an additional file name.  Augment
	  * number of input files. */
           ++nfiles;
       }
   }  /* end for each command line argument */
   
   if (! topt) {
       fprintf(stderr,"\nYou must specify a data type: c,b,f,s,m,or u\n");
       exit(1);
   }

   if (! bopt) {
       fprintf(stderr,"\nBad stations will not be saved in a separate file\n");
   }
   
   if (! nfiles) {
       fprintf(stderr,"\nExpecting data to come from STDIN...");
       fp = stdin;
   }
 
 
   /* Read in shipcodes file, downloaded from internet
    * with commas all converted to spaces. */  
    if (shipfile == NULL) {
       fprintf(stderr,"\nMust specify path to OCL/NODC shipcode file. \n");
       exit(1);
    }

    /* Initialize shipcode table */
    nodc_ship = (char **) calloc (maxscode, sizeof(char *));
    if (nodc_ship == NULL) {
          fprintf(stderr,"Unable to allocate memory for shipcode table\n");
          exit(1);
       }

    for (i = 0; i < maxscode; ++i) {
      nodc_ship[i] = (char *) calloc(4, sizeof(char));
      if (nodc_ship[i] == NULL) {
	fprintf(stderr,"Unable to allocate memory for shipcode table\n");
	exit(1);
       }
    }

    fprintf(stderr,"\nReading shipcode file.... ");
     /* Write each line to buf */
    while ((error = fscanf(shipfile,"%[^\n]",buf)) == 1) {
      getc(shipfile);  /* move past end of line marker */
   
      if ( vopt == 1 ) {
	fprintf(stderr,"\nCurrent line is: %s",&buf);
      }

      /* Second col is index to shipcode table*/
      error = sscanf(buf,"%[^','],%[^','],%s",shipcode,shipindex);
	
      if ( vopt == 1 ) {
	fprintf(stderr,"\n Index number is: %s",&shipindex);
	fprintf(stderr,"\n Shipcode is: %s",&shipcode);
      }

      /* Convert shipindex to integer */
      i = atoi(shipindex);

      if ( vopt == 1 ) {
	fprintf(stderr,"\n Converted index number is: %d",i);
      }

      if (i >= 0  && i < maxscode) {
	/* Test that we are within acceptable range of
	 * shipcodes and then read the shipcode from the
	 * buffer (2nd column). The first two char are
	 * the country code and the last two are the
	 * platform code. */
	st = &shipcode[0];
	for (j = 0; j < 4; ++j) {
	  nodc_ship[i][j] = *st;
	  ++st;
	}
	if ( vopt == 1 ) {
	  fprintf(stderr,"\n Stored shipcode is: %s \n",nodc_ship[i]);
	}

      } /* End of test for index range */
    }   /* End of looping over lines in file */
    
/********************************************************

 Initialize WOD dynamic arrays

*********************************************************/

  spacer(1);

/********************************************************

 Load profiles from the input file(s)

********************************************************/

  /* Enter standard level depths in case we
   * have standard level data.  This is only
   * needed for WOD09 and earlier.  WOD13
   * files read standard level depths from
   * the data by default.  There is a potential
   * bug here in case a single WOD file has
   * mixed standard and obs depth levels (or
   * the multiple input files are mixed).*/
  for ( j = 0; j < 40; j++ ) 
    *(depth+j)= *(sdepth+j);
  
  /* Loop for each file.  Variable curfile is init
   * to 1 in declarations. */
  do {
    if (nfiles) {
      /*fprintf(stderr,"\n Opening file: %s",argv[curfile]);*/
      fp = openfile(dir, argv[curfile], extent);
      if (fp == NULL) {
	fprintf(stderr,"\n Skipping %s because cannot open file.",argv[curfile]);
	goto NEXTFILE;
      } /* End of null file check */
    }
    
    /* Loop for each station.  The call to oclread()
     * will read in a profile from the input file
     * with the profile data + metadata going into
     * the WOD global variables declared above.*/
    while ( oclread() >= 0 && !feof(fp) ) {

      /* Augment the read station counter */
      ++staread;
	
      if (check_sta() > 0) {
	/* The station is good, so write entire station to output */
	write_hydro_station(outfile, &hdr, &data);

	/* Augment output station counter */
	++staout;
      }
	   
      if ( bad_hdr.nobs > 0 ) {
	if ( bopt == 1 ) {
	/* If we have -B<badfile> enabled, bopt = 1
	 * and if the profile has more than zero
	 * observations, write this profile to
	 * rejected station file. */
	write_hydro_station(badfile, &bad_hdr, &bad_data);
	}

	/* Augment the rejected station counter */
	++stabad;
      }

      /* Clean up -> clear all pointers in both good
       * and bad header/data information. */

      if (hdr.prop_id != NULL) {
	free((void *)hdr.prop_id);
	hdr.prop_id = NULL;
      } 
                 
      if (bad_hdr.prop_id != NULL) {
	free((void *)bad_hdr.prop_id);
	bad_hdr.prop_id = NULL;
      } 
      
      for (i=0; i< MAXPROP; ++i) {
	if (data.observ[i] != NULL) {
	  free((void *) data.observ[i]);
	  data.observ[i] = NULL;
	}
	if (bad_data.observ[i] != NULL) {
	  free((void *) bad_data.observ[i]);
	  bad_data.observ[i] = NULL;
	}
      }  /* End of clean up looping over all properties */

    }  /* End of while looping over all stations*/
    
  NEXTFILE:
    if (nfiles) {
      fclose(fp);
    }
  } while (curfile++ < nfiles);  /* End of while looping over files */
  
  fprintf(stderr,"\nEnd of conversion.");
  fprintf(stderr,"\n  %d stations read in", staread);
  fprintf(stderr,"\n  %d stations accepted", staout);
  fprintf(stderr,"\n  %d stations rejected ", staread-staout);
  fprintf(stderr,"\n  %d stations contained bad scans \n\n", stabad - (staread-staout));
  exit(0);
  
} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program) {

   fprintf(stderr,"\nConverts .OSD, .CTD, and .PFL files (observed levels)from WOD05 to HydroBase\n");
   fprintf(stderr,"\nUsage: %s filelist -S<shipfile> -T<instr_type> [-B<badfile>] [-D<dirname>] [-E<extent>] [-h] [-O<outfile>]", program);
   fprintf(stderr,"\n\n  List of filenames must be first argument!");
   fprintf(stderr,"\n  Otherwise, input is expected from stdin.");
   fprintf(stderr,"\n-S : path/filename for OCL/NODC shipcode table.");
   fprintf(stderr,"\n-T : specifies data_type code.  Choose one:");
   fprintf(stderr,"\n   [c]td");
   fprintf(stderr,"\n   [b]ottle");
   fprintf(stderr,"\n   [f]loat");
   fprintf(stderr,"\n   [s]easoar");
   fprintf(stderr,"\n   [m]oored profiler");
   fprintf(stderr,"\n   [u]nknown  \n");
   fprintf(stderr,"\n\n   OPTIONS: ");
   fprintf(stderr,"\n-B : specifies output file for rejected stations (noclobber)");
   fprintf(stderr,"\n          defaults to stdout");
   fprintf(stderr,"\n-D : dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n-E : input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n-O : output file (default overwrites an existing file)");
   fprintf(stderr,"\n-h : help...... prints this message. ");
   fprintf(stderr,"\n-V : Verbose, mostly for debugging, to stderr.");
   fprintf(stderr,"\n\n");
   fprintf(stderr,"\nBe careful that your file names are not too long.");
   fprintf(stderr,"\nLong file names can cause stack smashing errors");
   fprintf(stderr,"\nsince there are no explicit checks for string length.");
   return;
}
   
/****************************************************************************
 * This openfile opens the file whose name is given by the inputs:
 * dir/root.extent.  A pointer to the file is returned. */

FILE *openfile(char *dir, char *root, char *extent) {

  char st[80];  /* Buffer to become file name */
  int i;
  FILE *infile;
   
  /* Copy directory info to file name */
  strcpy(st, dir);

  /* Check if dir ends in '/', if not, add '/'. */
  if ((i=strlen(dir)) != 0) {
    if (dir[i-1] != '/')
      strncat(st,"/",1);
  }

  /* Add the file root name and extentsion */
  strcat(st, root);
  strcat(st, extent);

  /* Open the file for reading only */
  infile = fopen(st,"r");

  /* Test whether the open was successful */
  if (infile == NULL)
    fprintf(stderr,"\n Unable to open %s \n", st);
  else
    fprintf(stderr,"\n   opened %s ... ", st);
  
  /* Pass the pointer to the file back to calling program */
  return(infile);
   
}  /* end openfile() */

/****************************************************************************
 * The function check_sta returns the number of observations
 * in the station (levels) if the station has acceptable
 * data, 0 if not.
 * 
 * Basically, this routine copies all the data that was
 * loaded by oclread() into the global WOD-formatted arrays
 * into the global hydrobase structures.  A direct copy is
 * put into the bad structure and an edited (only good profiles
 * and scans retained) is copied to the good structure.
 *
 * No explicit inputs are needed in the function definition
 * because global data are read from the WOD and HB variables.
 ****************************************************************************/

int check_sta() {

  /* Arrays for keeping track of the presence and
   * status/flag of each variable in the profile */
  int *prop_avail, *prop_OK;

  /* Storage of index of property within arrays */
  int *dindex;

  int **data_flagged;
  int offset, i, ii, j, nobs, delta;
  int found, index;
  float prsint;
  int l_t68;
  double t90_tmp;
  double t90_to_t68;
  
  /* Set conversion factor from T90 to T68.
   * T90 values * t90_to_t68 => T68 values
   * All output will be in T68. */
  t90_to_t68 = 1.00024;

   /* Allocate some space */
   prop_avail = (int *) calloc(MAXPROP, sizeof(int));
   prop_OK = (int *) calloc(MAXPROP, sizeof(int));
   dindex = (int *) calloc(MAXPROP, sizeof(int));
   data_flagged = (int **) calloc(MAXPROP, sizeof(int *));
   
   /* Store metadata in HydroBase station header.
    * Note that we store duplicate information in
    * both the good header and the bad header. */      
   hdr.lat = (float) latitude / tenp[*(hrightfig+1)];
   bad_hdr.lat = hdr.lat;

   hdr.lon = longitude/ tenp[ *(hrightfig+2)];
   bad_hdr.lon = hdr.lon;

   hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
   bad_hdr.ms10 = hdr.ms10;

   bad_hdr.ms1 = hdr.ms1;

   hdr.origin = '1';
   bad_hdr.origin = '1';

   /* HB3 version has slightly different
    * country assignment.  The country
    * assignment is done again based on
    * the ship codes table, so this value
    * only persists if the ship table
    * cannot be looked up properly.  Use
    * this as a double check for ship
    * table parsing. */
   strncpy(hdr.country, cc, 2);
   hdr.country[2] = '\0';
   /* Older HB2 version:
   strncpy(hdr.country, cc, 3); 
   */
   strncpy(bad_hdr.country, cc, 3);

   hdr.year = year;
   bad_hdr.year = year;

   hdr.month = month;
   bad_hdr.month = month;

   hdr.day = day;
   bad_hdr.day = day;
   
   /* Compute default values for cruise and station # */
   hdr.cruise =  icruise - ((icruise/100000) * 100000);
   hdr.station = ostation - ((ostation/10000) * 10000);
   
   /* Use originator's cruise and station #s if available.
    * Start by setting end of string delimiters. */
   *(origc+origcfig) = '\0';
   *(origs+origsfig) = '\0';

   /* If there are figures for the cruise #, set cruise #.*/
   if (origcfig > 0) {      
      if ((i = atoi(origc)) > 0)
         hdr.cruise = i - ((i/100000) * 100000);
   }

   /* If there are figures for the station #, set station #.*/
   if (origsfig > 0) {
      if ((i = atoi(origs)) > 0)
         hdr.station = i - ((i/10000) * 10000);
   }
   bad_hdr.cruise = hdr.cruise;
   bad_hdr.station = hdr.station;
   
   /*Search for a platform code in secondary header arrays  */
   found = 0;
   i = 0;
   while (!found && (i < nsec)) {
     /* secondary header code for OCL Platform
      * so if seccode[i] == 3, then found is assigned
      * a value of true, otherwise, stays at 0. */
      found = seccode[i] == 3;  

      /* Move onto the next secondary header. */
      ++i;
   }
   
   /* Default ship code */
   strncpy(hdr.ship, "XX", 3);
   if (found) {
     /* Move back one to the secondary header that
      * resulted in the correct match, above. */
      --i;
      index = secval[i] / tenp[srightfig[i]];
      if (index >= maxscode) {
         fprintf(stderr,"\nError parsing index to shipcodes.asc table: %d\n",index);
      }
      else {
	/* We have a valid ship code */
	if (index > 0 && (nodc_ship[index][2] != '\0')) {
	  hdr.ship[0] = nodc_ship[index][2]; 
	  hdr.ship[1] = nodc_ship[index][3]; 
	  hdr.ship[2] = '\0';
	}
      }
   }
   
   /* Copy ship code to both good and bad struct */
   strncpy(bad_hdr.ship, hdr.ship, 3);
   
   /*Search for a PDR depth in secondary header arrays  */
   found = 0;
   i = 0;
   while (!found && (i < nsec)) {
      found = seccode[i] == 10;    /* secondary header code for seafloor depth */
      ++i;
   }

   /* Default depth */
   hdr.pdr = bad_hdr.pdr = 0;   
   if (found) {
      --i;   
      hdr.pdr = secval[i] / tenp[srightfig[i]]; 
      bad_hdr.pdr = hdr.pdr; 
   }
      
   /* Determine temperature scale - new for WOD09, WOD13
    * There are the scale codes as listed by:
    * http://data.nodc.noaa.gov/woa/WOD/CODES/TXT/v_3_scale.txt
    * 102,Temperature: T68 (IPTS-68) -> TE in HB2
    * 103,Temperature: ITS-90        -> T90 in HB3
    * HB2 does not have ITS90 capability, so instead
    * of writing to TE or T90, we must convert
    * any explicitly noted T90 to TE and always output
    * temperature data in TE -> ITS-68.*/
   found = 0;
   ii = 0;
   l_t68 = 1;          /* Flag at 1 for TE, 0 for T90 */
   /* Loop over all variable specific secondary headers */
   while (!found  && (ii < npsec)) {
     /* psecparm == 1 is a flag saying this is a temperature header
      * pseccode is the code value, if == 3, then it's a scale b/c table 3
      * psecval is the value corresponding to this code, 102 or 103 */
     if ( (psecparm[ii] == 1) && (pseccode[ii] == 3)) {
       found = 1;
       /* OCL code for T90 == 103, for T68 == 102 */
       if (psecval[ii] == 103) {
	 /* Change from default T68 to T90 */
	 l_t68 = 0;
       }else{
	 if ( psecval[ii] == 102 ) {
	   /* All is OK - we know that this is explicitly T68 */
	 }else{
	   fprintf(stderr,"\n Unknown temperature scale! Assuming T68.");
	 }
       }
     }
     ++ii;
   }

   /* Check to see what properties are available.
    * Variable codes are determined from WOD tutorial,
    * Table 3, page 17. */

   /* Also, perform a check on the profile flags.
    * check_profile_flag will return 1 for any
    * standard WOD05 flag - so it looks like all
    * profiles will be included.  The real data
    * checking is happening at each level with
    * the data flags. */

   /* TE, SA, OX, ... etc. are declared as an
    * enum property in hydorbase.h */

   for (i = 0; i < nparm; ++i) {
     switch (ip2[i]) {
     case 1:
       /* WOD temp code = 1 */
       prop_avail[(int)TE] = 1;	         
       prop_OK[(int)TE] = check_profile_flag(iperror[i]);
       dindex[(int)TE] = i;
       break;
     case 2:
       /* WOD salt code = 2 */
       prop_avail[(int)SA] =1;
       prop_OK[(int)SA] = check_profile_flag(iperror[i]);
       dindex[(int)SA] = i;
       break;
     case 3:
       /* WOD oxygen code = 3 */
       prop_avail[(int)OX] = 1;
       prop_OK[(int)OX]  = check_profile_flag(iperror[i]);
       dindex[(int)OX] = i;
       break;
     case 4:
       /* WOD phosphate code = 4 */
       prop_avail[(int)P4] = 1;
       prop_OK[(int)P4] = check_profile_flag(iperror[i]);
       dindex[(int)P4] = i;
       break;
     case 6:
       /* WOD silicate code = 6 */
       prop_avail[(int)SI] = 1;
       prop_OK[(int)SI] = check_profile_flag(iperror[i]);
       dindex[(int)SI] = i;
       break;
     case 7:
       /* WOD code seven is undefined in WOD05! */
       fprintf(stderr,"Warning: WOD05 does not define N2 = 7!");
       prop_avail[(int)N2] = 1;
       prop_OK[(int)N2] = check_profile_flag(iperror[i]);
       dindex[(int)N2] = i;
       break;
     case 8:
       /* WOD nitrate + nitrite code = 8 */
       prop_avail[(int)N3] = 1;
       prop_OK[(int)N3] = check_profile_flag(iperror[i]);
       dindex[(int)N3] = i;
       break;
     case 25:
       /* WOD pressure code = 25 */
       prop_avail[(int)PR] = 1;
       prop_OK[(int)PR] = check_profile_flag(iperror[i]);
       dindex[(int)PR] = i;
       break;
     case 33:
       /* WOD tritium code = 33 */
       prop_avail[(int)TU] = 1;
       prop_OK[(int)TU] = check_profile_flag(iperror[i]);
       dindex[(int)TU] = i;
       break;
     case 34:
       /* WOD helium code = 34 */
       prop_avail[(int)HE] = 1;
       prop_OK[(int)HE] = check_profile_flag(iperror[i]);
       dindex[(int)HE] = i;
       break;
     case 40:
       /* WOD CFC11 code = 40 */
       prop_avail[(int)F1] = 1;
       prop_OK[(int)F1] = check_profile_flag(iperror[i]);
       dindex[(int)F1] = i;
       break;
     case 41:
       /* WOD CFC12 code = 41 */
       prop_avail[(int)F2] = 1;
       prop_OK[(int)F2] = check_profile_flag(iperror[i]);
       dindex[(int)F2] = i;
       break;
     case 42:
       /* WOD CFC113 code = 42 */
       prop_avail[(int)F3] = 1;
       prop_OK[(int)F3] = check_profile_flag(iperror[i]);
       dindex[(int)F3] = i;
       break;
     default:
       ;	   
     } /* End switch */
   } /* End for looping over all params to check availability/flags */
   
   /* Load property data into HydroBase station structures */
   hdr.nprops = 0;	
   bad_hdr.nprops = 0;	
   
   for (i = 0; i < MAXPROP; ++i) {
     if (prop_avail[i] ) {

       /* Allocate space */
       data.observ[i] = (double *) malloc(levels * sizeof(double));
       bad_data.observ[i] = (double *) malloc(levels * sizeof(double));
       data_flagged[i] = (int *) calloc(levels, sizeof(int));

       /* If we have a temperature variable and the temp
	* scale was detected as T90, then we need to
	* convert to T68 as we copy from infile to HB2. */
       if ( i == (int)TE && l_t68 == 0 ) {
	 /* Temperature variable and need to convert.  Same code
	  * as below but with temperature scale conversion hard
	  * coded in and fewer comments. */
	 for (j=0; j < levels; ++j) {
	   offset = dindex[i]*levels;
	   t90_tmp = (double) (*(dataval+offset+j) /tenp[*(drightfig+offset+j)]);
	   data.observ[i][j] = t90_tmp*t90_to_t68;
	   data_flagged[i][j] = check_obs_flag(*(derr+offset+j)) ? 0 : 1;
	   if (t90_tmp < -9.) 
	     data_flagged[i][j] = 1;
	 }
       }else{
	 /* No need to convert, could be temp or other variable. */

	 /* Copy data from dataval (read from infile) to hydrobase data */
	 for (j=0; j < levels; ++j) {

	   /* Offset for counting between columns in
	    * an array - pointer arithmetic */
	   offset = dindex[i]*levels;

	   /* Convert to real values using figs wrt decimal.*/	 
	   data.observ[i][j] = (double) (*(dataval+offset+j) /tenp[*(drightfig+offset+j)]);

	   /* Check flags on data - if check_obs_flag returns 1 = true,
	    * then the conditional operator will assign a 0 to data_flagged.
	    * Otherwise, we get a 1 to data_flagged.  Note that this
	    * operation flips the T/F result of check_obs_flag. */
	   data_flagged[i][j] = check_obs_flag(*(derr+offset+j)) ? 0 : 1;

	   /* Check missing value? */
	   if (data.observ[i][j] < -9.) 
	     data_flagged[i][j] = 1;

	 } /* End of loop over all levels */
       }  /* End of check for temperature conversion */

       /* Augment number of profiles to bad struct b/c
        * bad struct incrementally stores all data,
        * regardless of good or bad. */
       ++bad_hdr.nprops;

       if ( prop_OK[i])
	 /* Add number of valid profiles to good struct */
	 ++hdr.nprops;
     }
   } /* End of loop over all properties (variables) to load data */
	   
   /* Load depths */
   prop_avail[(int)DE] = 1;  /* By definition */
   prop_OK[(int)DE] = 1;

   /* Augment the number of values in good/bad struct b/c adding depth prop. */
   ++hdr.nprops;
   ++bad_hdr.nprops;

   /* Allocate space for depth */
   data.observ[(int)DE] = (double *)malloc (levels * sizeof(double));
   bad_data.observ[(int)DE] = (double *)malloc (levels * sizeof(double));
   data_flagged[(int)DE] = (int *) calloc(levels, sizeof(int));

   for (j=0; j < levels; ++j) {

     /* Copy depth information */
     data.observ[(int)DE][j] = (double) (depth[j] / tenp[zrightfig[j]]);

     /* Copy flag for variable.  Note T/F flip from
      * conditional operator. */
     data_flagged[(int)DE][j] = check_obs_flag(zerr[j]) ? 0 : 1;

     if (data.observ[(int)DE][j] < -9.) 
       data_flagged[(int)DE][j] = 1;
   } /* End of loop over all levels */
   
   /* If pressure wasn't a stored parameter, compute it. */
   if (!prop_avail[(int)PR]) {

     /* Allocate space for pressure */
     data.observ[(int)PR] = (double *)malloc (levels * sizeof(double));
     bad_data.observ[(int)PR] = (double *)malloc (levels * sizeof(double));
     data_flagged[(int)PR] = (int *) calloc(levels, sizeof(int));

     /* Compute the pressure at each level */
     for (j=0; j<levels; ++j) {
       data.observ[(int)PR][j] = hb_p80(data.observ[(int)DE][j], (double)hdr.lat);
       data_flagged[(int)PR][j] = data_flagged[(int)DE][j];
     } /* End of loop over all levels */

     /* Include pressure as a property in both bad
      * and good structures */
     prop_avail[(int)PR] = 1;
     prop_OK[(int)PR] = 1;  
     ++hdr.nprops;
     ++bad_hdr.nprops;
   } /* End of pressure variable presence check */

   /* Check for requisite parameters, must
    * have both T and S.  If either is missing,
    * we return 0 and go back to main program. */
   if ( ! (prop_OK[(int)TE] && prop_OK[(int)SA])) {
   
     /* Write station in bad file if missing T or S */
     for (i=0; i< MAXPROP; ++i) {
       for (j=0; j < levels; ++j) {
	 if (prop_avail[i])
	   bad_data.observ[i][j] = data.observ[i][j];
       }
     } /* End of copying station data to bad struct */
     
     /* Create a prop id for each prop */
     bad_hdr.prop_id = (int *) malloc(bad_hdr.nprops * sizeof(int));

     /* Assign an ID to each prop id */
     j = 0;
     for (i=0; i< MAXPROP; ++i) {
       /* Assign an id iff the prop is available */
       if (prop_avail[i]) {
	 bad_hdr.prop_id[j++] = i;
	 free((void *)data_flagged[i]);
       } /* End of prop available check */
     } /* End of loop over all properties to assign ID */

     /* Assign nobs and nprops */
     bad_hdr.nobs = bad_data.nobs = levels; 
     bad_data.nprops = bad_hdr.nprops;

     /* Clear memory */
     free((void *)data_flagged);
     free((void *)prop_avail);
     free((void *)prop_OK);
     free((void *)dindex);

     /* Quit this function and return to main b/c we reject this profile */
     return(0);

   } /* End of copying whole station to bad struct */
   
   
   /* Tally number of good/bad levels. A bad level = p,t,or s flagged 
      in which case the level is moved to bad structure.  */
   hdr.nobs = 0;  
   bad_hdr.nobs = 0;

   for (j=0; j<levels; ++j) {
     
     /* All accepted points have check_obs_flag = 1, but data_flagged = 0 */
     if (data_flagged[(int)DE][j] || data_flagged[(int)PR][j] || data_flagged[(int)TE][j] || data_flagged[(int)SA][j] || data.observ[(int)PR][j] < 0 ) {
       /* These are rejected data! 
        * Copy all data to the bad struct */
       for (i=0; i< MAXPROP; ++i) {
	 if (prop_avail[i])
	   bad_data.observ[i][bad_hdr.nobs] = data.observ[i][j];
       } /* End of loop over each property (variable) to copy data */

       /* Increase the counter for the bad struct */
       ++bad_hdr.nobs;
     } /* End of bad data check */
     else  {
       /* Accepted data! */
       for (i=0; i< MAXPROP; ++i) {
	 if (prop_OK[i]) {
	   /* Copy the property/variable data to good struct */
	   data.observ[i][hdr.nobs] = data.observ[i][j];

	   /* Any individual flagged values are set to missing */
	   if (data_flagged[i][j])
	     data.observ[i][hdr.nobs] = MISSING;
	 } /* End of individual data flag check */
       } /* End of good property copy */

       /* Augment the number of good observations */
       ++hdr.nobs;
     } /* End of accepted data check */
   }  /* End loop over all levels */ 

   /* Update nobs in both data and header struct */
   data.nobs = hdr.nobs ;  
   bad_data.nobs = bad_hdr.nobs; 
   
   /* Creat space for property IDs - previous prop_id
    * allocation is performed only if we are missing
    * T or S.  In that case, this section of code is
    * not reached. */
   hdr.prop_id = (int *) malloc(hdr.nprops * sizeof(int));
   bad_hdr.prop_id = (int *) malloc(bad_hdr.nprops * sizeof(int));

   /* Initialize pointer arithmetic counters */
   j = 0;
   i = 0;

   /* Assign a prop ids :
    * (a) to each OK property in the good struct
    * (b) to each available property in the bad struct */
   for (index = 0; index < MAXPROP; ++index) {

     if (prop_OK[index]) {
       /* We have a acceptable property/profile */
       hdr.prop_id[j] = index;

       /* Check that we do not exceed the number of properties.
	* This is a self consistency double-check. */
       if (++j > hdr.nprops) {
           fprintf(stderr,"\nError counting hdr.nprops in check_sta()!\n");
           exit(1);
       } /* End of # props check */
     } /* End of prop OK check */

     if (prop_avail[index]) {
       bad_hdr.prop_id[i] = index;
       free((void *)data_flagged[index]);
       if (++i > bad_hdr.nprops) {
	 fprintf(stderr,"\nError counting bad_hdr.nprops in check_sta()!\n");
           exit(1);
       } /* End of # props check */
     } /* End of prop available check */
   }
   
   /* Check interval of pressure series.
    * Decimate to 10 db if delta-p is smaller */
   if (hdr.nobs > 30) {
     /* There are enough observations to warrant decimation check */
     
     /* Determine average pressure interval from the 15th to the 30th points */
     prsint = 0; /* Initialize the running sum */
     for (i = 15; i < 30; ++i) {
       prsint += (float)(data.observ[(int)PR][i+1] - data.observ[(int)PR][i]);
     } /* End of computing the running sum */

     /* Compute the average value */
     prsint = prsint / 15;
     
     /* Check the average pressure interval */
     if (prsint < 10) {

       /* Assign a new pressure interval in index units:
	* Delta can be integers [1,10] which correspond to
        * the span [prsint close to ten, prsint much smaller than ten]. */
       delta = (int)(10 / prsint);

       /* Cycle through each of the properties/profiles in the station */
       for (i = 0; i < hdr.nprops; ++i) {
	 nobs = 0;   /* Initialize a counter */

	 /* Copy the data in strides of delta */
	 for (j = 0; j < hdr.nobs; j += delta) {
	   data.observ[hdr.prop_id[i]][nobs] = data.observ[hdr.prop_id[i]][j];
	   ++nobs;
	 } /* End of stiding loop over all observational levels */

	 /* Check whether or not the last stride went past
	  * the level closest to the bottom.  If it did, copy
	  * the bottom value into the data. */
	 if ((j - delta) < (hdr.nobs-1)) {
	   data.observ[hdr.prop_id[i]][nobs] = data.observ[hdr.prop_id[i]][hdr.nobs-1];
	    ++nobs;
	 }
       } /* End of looping over each observed level */

       /* If we have decimated the data, update the number of
        * observational levels. Note that the skipped/strided
        * over data still resides in the data array, but the only
        * in-order data is in the first section of data.observ
        * and we update nobs in the structure to basically wipe
        * knowledge of the residual data. */
       hdr.nobs = data.nobs = nobs;
     } /* End of pressure interval check */
   } /* End of decimation check

   /* Deallocate/free up memory in preparation for the next station */
   free((void *)prop_avail);
   free((void *)prop_OK);
   free((void *)dindex);
   free((void *)data_flagged);
   
   /* Return the number of good observations */
   return (hdr.nobs) ;
   
}  /* End check_sta() */

/****************************************************************************
 * The function check_profile_flag will return 1 for a single digit
 * integer input [0,9].  For any other integer input, zero is returned.
 *
 * Flag = 0 = accepted
 * Flag = 1 = failed annual standard deviation check (Ruth allows)
 * Flag = 2 = two or more density inversions (Ruth allows)
 * Flag = 3 = flagged cruise
 * Flag = 4 = failed seasonal standard deviation check
 * Flag = 5 = failed monthly standard deviation check
 * Flag = 6 = failed annual and seasonal standard deviation check
 * Flag = 7 = failed annual and monthly standard deviation check
 * Flag = 8 = failed seasonal and monthly standard deviation check
 * Flag = 9 = failed annual, seasonal and monthly standard deviation check
 *
 * Since the only WOD flags are [0,9] and 0 is the accepted cast flag,
 * if you add flags to the case list, then those values will be ignored
 * flags.
 ****************************************************************************/
int check_profile_flag(int flag) {

  switch (flag) {
   
  case 0:   /* fall through */
  case 1:   /* Ruth allows */
  case 2:   /* Ruth allows */
    /*
  case 3:   Flagged cruise!
  case 4:
  case 5:
  case 6:
  case 7:
  case 8:
  case 9:
    */
    return 1;
  default:
    return 0;
  } 
  
} /* end check_profile_flag() */

/****************************************************************************
 * The function check_obs_flag will return 1 for integer inputs listed
 * in the case n: lines and returns zero for any other integer input.
 *
 * Flag = 0 = accepted value
 * Flag = 1 = failed broad range check
 * Flag = 2 = failed inversion check
 * Flag = 3 = failed gradient check
 * Flag = 4 = level bullseye flag
 * Flag = 5 = Flags 2 + 3
 * Flag = 6 = Flags 1 + 2
 * Flag = 7 = Flags 1 + 3
 * Flag = 8 = Flag 1 + zero anamoly check
 * Flag = 9 = Flags 1 + 2 + 3
 *
 * Ultimately, the flags listed in this routine are "flags to ignore."
 * The range check and bullseye check are the flags we're watching out for.
 *
 * sfg added flag 1 b/c many PFL stations rejected otherwise.  This
 * does not seem to affect OSD and CTD station rejection.
 *
 * sfg added flag 5 b/c flag 5 is the logical combination of flag 2 + flag 3
 ***************************************************************************/
int check_obs_flag(int flag) {

   switch (flag) {
   
   case 0:   /* fall through */
   case 1:   /* sfg added */
   case 2:   /* Ruth allows */
   case 3:   /* Ruth allows */
     /*
   case 4: Bullseye!
     */
   case 5:   /* sfg added b/c 2+3, allowed above */
   case 6:   /* sfg added b/c 1+2 */
   case 7:   /* sfg added b/c 1+3 */
     /*
   case 8:
     */
   case 9:   /* sfg added b/c 1+2+3 */
     return 1;
   default:
     return 0;
   } 
   
} /* end check_obs_flag() */

/****************************************************************************
 * The function oclread will read one profile at a time from the currently
 * opened WOD05-formatted file.  All data and meta-data are stored in the
 * WOD global variables declared at the top.
 *
 * If an EOF is encountered, this function returns -1.
 *
 * This code, along with the helper functions nocrfgetc, spacer,
 * extracti, extractc, and oclread were directly imported from
 * the wesite distributing WOD05:
 * http://www.nodc.noaa.gov/OC5/SELECT/dbsearch/sysinfo.html
 *
 * SFG compared the WOD09 version of this with the updated WOD13
 * version and made changes to reflect WOD09.  My guess is that
 * these changes are still backwards compatible with WOD09 data.
 ***************************************************************************/

int oclread() {

  int i,j;

  char wodform;
  int totfig=0, sigfig=0, rightfig=0;
  int nbytet,ntypec,nbytec,nbytes,nbyteb;
  int ninfc,ntoff, doff;
  int missing=-9999;
  int npinfs=0,npinfe=0,npinf;
  int iend=0;

  /* Read in WOD format code
   * 'C' is for WOD13
   * 'B' is for WOD05,09
   * 'A' is for WOD01  */
  totfig= 1;
  if ( (iend = extractc(0,&totfig,&wodform)) == -1 ) {
    fprintf(stderr,"oclread: Cannot extract wodform!");
    return iend;
  }
  /* Added this line for WOD13, not present in WOD09 */
  if ( wodform == 'C' ) iVERSflag=2;

  /* Read in number of bytes in ASCII station, internal integer */
  if ( ( iend = extracti(0,&totfig,&sigfig,&rightfig,&nbytet,missing)) == -1 ) {
    fprintf(stderr,"oclread: Cannot extract nbytet!");
    return iend;
  }
  /* Read in OCL station number */
  if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&ostation,missing)) == -1 ) {
    fprintf(stderr,"oclread: Cannot extract ostation!");
    return iend;
  }
  /* Read in NODC country code */
  totfig= 2;
  if ( (iend = extractc(0,&totfig,cc)) == -1 ) {
    fprintf(stderr,"oclread: Cannot extract cc!");
    return iend;
  }
  /* Read in cruise number */
  if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&icruise,missing)) == -1) {
    fprintf(stderr,"oclread: Cannot extract icruise!");
    return iend;
  }
  /* Read in year, month, day, time */
  totfig=4;
  if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&year,missing)) == -1) {
    fprintf(stderr,"oclread: Cannot extract year!");
    return iend;
  }
  totfig=2;
  if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&month,missing)) == -1) {
    fprintf(stderr,"oclread: Cannot extract month!");
    return iend;
  }
  totfig=2;
  if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&day,missing)) == -1) {
    fprintf(stderr,"oclread: Cannot extract day!");
    return iend;
  }
  if ( (iend = extracti(1,htotfig,hsigfig,hrightfig,&hour, 9999)) == -1) {
    fprintf(stderr,"oclread: Cannot extract hour!");
    return iend;
  }
  /* Read in latitude and longitude */
  if ( (iend = extracti(1,(htotfig+1),(hsigfig+1),(hrightfig+1),&latitude, -9999)) == -1) {
    fprintf(stderr,"oclread: Cannot extract latitude!");
    return iend;
  }
  if ( (iend = extracti(1,(htotfig+2),(hsigfig+2),(hrightfig+2),&longitude, -99999)) == -1) {
    fprintf(stderr,"oclread: Cannot extract longitude!");
    return iend;
  }
  /* Read in number of levels */
  if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&levels,missing)) == -1) {
    fprintf(stderr,"oclread: Cannot extract levels!");
    return iend;
  }
  /* Read in observed (0) or standard (1) levels */
  totfig=1;
  if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&isoor,missing)) == -1) {
    fprintf(stderr,"oclread: Cannot extract isoor!");
    return iend;
  }
  /* Read in number of parameters at this station (not incl. bio) */
  totfig=2;
  if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&nparm,missing)) == -1) {
    fprintf(stderr,"oclread: Cannot extract nparm!");
    return iend;
  }
  /* Read in each parameter code and profile error code */
  for ( i =0; i< nparm; i++ ) {

    if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(ip2+i),
			  missing)) == -1) return iend;
    totfig=1;
    if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(iperror+i),
			  missing)) == -1) return iend;

    /* Read number of parameter specific second header variables */
    if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&npinf,missing))
	 == -1 ) return iend;
    npinfe += npinf;

    /* Read in each parameter specific second header variable */
    for ( j = npinfs; j < npinfe; j++ ) {

      *(psecparm+j) = *(ip2+i);
      if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(pseccode+j),
			    missing)) == -1) return iend;
      if ( (iend = extracti(1,(pstotfig+j),(pssigfig+j),(psrightfig+j),
			    (psecval+j), missing)) == -1) return iend;
      
    } /* End of loop over each param. spec. sec. header var. */

    npinfs += npinf;

  } /* End of loop over each param. */

  npsec = npinfe;

  /* Read in number of bytes in character and PI fields */
  if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&nbytec,
			missing)) == -1) return iend;

  /* Read in num of info types (max 3: cruise code, station code, PI) */
  origcfig= 0;
  origsfig= 0;
  npi= 0;  
  if ( nbytec > 0 ) {
    totfig=1;
    if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&ninfc,missing))
	 == -1) return iend;
    
    /* Read in type of info: 1 = cruise code, 2 = station code, 3 = PI */
    for ( i= 0; i < ninfc; i++) {

      totfig=1;
      if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&ntypec,missing))
	   == -1) return iend;

      /* Read in originators cruise code */
      if ( ntypec == 1) {

	totfig=2;
	if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&origcfig,
			      missing)) == -1) return iend;
	if ( (iend = extractc(0,&origcfig,origc)) == -1 ) return iend;
	
      } /* End of cruise code check */

      /* Read in originators station code */
      else if ( ntypec == 2 ) {

	totfig=2;
	if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&origsfig,
			      missing)) == -1) return iend;
	if ( (iend = extractc(0,&origsfig,origs)) == -1 ) return iend;
	
      } /* End of station code check */

      /* Read in PI information */
      else if ( ntypec == 3 ) {
	
	/* Read in number of PI at this station */
	totfig=2;
	if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&npi,missing))
	     == -1) return iend;

	/* Read in each param. code and PI code for that param. */
	for ( j =0; j< npi; j++ ) {
	  
	  if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(ipip+j),
				missing)) == -1) return iend;
	  if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(ipi+j),
				missing)) == -1) return iend;
	} /* End of loop over each PI code */

      } /* End of PI info check */

    } /* End of info type check (cruise, station, PI) */

  } /* End of amount of info check */

  /* Read in number of bytes in sec. header fields */
  if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&nbytes,
			missing)) == -1) return iend;
  nsec = 0;
  if ( nbytes > 0 ) {

    /* Read in number of sec. header params. */
    if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&nsec,
			  missing)) == -1) return iend;

    /* Read in each sec. header param. code and value */
    for ( i =0; i< nsec; i++ ) {
      
      if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(seccode+i),
			    missing)) == -1) return iend;
      if ( (iend = extracti(1,(stotfig+i),(ssigfig+i),(srightfig+i),
			    (secval+i), missing)) == -1) return iend;
      
    } /* End of loop over each sec. header param */
 
  } /* End of sec. header fields check */

  /* Read in number of bytes in bio header and taxa set fields */
  if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&nbyteb,
			missing)) == -1) return iend;
  nbio= 0;
  ntsets= 0;
  if ( nbyteb > 0 ) {
    
    /* Read in number of bio header param. */
    if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&nbio,
			  missing)) == -1) return iend;

    /* Read in each bio header param. code and value */
    for ( i =0; i< nbio; i++ ) {
      
      if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(biocode+i),
			    missing)) == -1) return iend;
      if ( (iend = extracti(1,(btotfig+i),(bsigfig+i),(brightfig+i),
			    (bioval+i), missing)) == -1) return iend;
      
    } /* End of loop over bio header params. */
  
    /* Read in number of taxa sets */
    if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&ntsets,
			  missing)) == -1) return iend;

    /* If # taxa sets > prev. max, reallocate space */
    if ( ntsets > ntsetsmax ) {

      ntsetsmax= ntsets;
      spacer(2);

    } /* End of size of taxa set check */

    if ( ntsets > 0 ) {

      /* Read in # entries for this taxa set, calc. read data offset */
      for ( j = 0; j < ntsets; j++) {

	if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(ntloc+j),
			      missing)) == -1) return iend;
	
	ntoff= maxtax * j;

	/* Read each taxa set param. code, value, error flag, and orig. flag */
	for ( i =0; i< *(ntloc+j); i++ ) {

	  if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(ntcode+ntoff+i),
				missing)) == -1) return iend;
	  if ( (iend = extracti(1,(nttotfig+ntoff+i),(ntsigfig+ntoff+i),
				(ntrightfig+ntoff+i), (ntval+ntoff+i), missing))
	       == -1) return iend;
	  totfig=1;
	  if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(nterr+ntoff+i),
				missing)) == -1) return iend;
	  totfig=1;
	  if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(ntoerr+ntoff+i),
				missing)) == -1) return iend;
	  
	} /* End of looping over each taxa set param code, etc. */

      } /* End of looping over each entry of taxa set. */
  
    } /* End of number of taxa set check */

  } /* End of bio header check */

  /* If # depths is greater than prev. max, reallocate space */
  zsize= levels;
  /* WOD09 commented out */
  /*if ( isoor == 0 && zsize > zsizemax ) { */
  /* WOD13 new */
  if ( zsize > zsizemax ) {
    zsizemax=zsize;
    spacer(3);
  } /* End of depths size check */

  /* If # of depth * # param. > prev. max, reallocate space */
  isize= nparm * levels;
  if ( isize > isizemax ) {
    isizemax=isize;
    spacer(4); 
  } /* End of data storage size check */

  /* Read in each depth value, error flag, and orig. flag */
  for ( j = 0; j < levels; j++ ) {

    /* WOD13 adds || wodform == 'C' */
    if ( isoor == 0 || wodform == 'C' ) {

      if ( (iend = extracti(1,(ztotfig+j),(zsigfig+j),
			    (zrightfig+j), (depth+j), missing))
	   == -1) return iend;
      totfig=1;
      if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(zerr+j),
			    missing)) == -1) return iend;
      totfig=1;
      if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(zoerr+j),
			    missing)) == -1) return iend;

    } /* End of obs or std depth level check */

    /* Read in each data value */
    for ( i =0; i< nparm; i++ ) {

      doff= i * levels;

      if ( (iend = extracti(1,(dtotfig+doff+j),(dsigfig+doff+j),
			    (drightfig+doff+j), (dataval+doff+j), missing))
	   == -1) return iend;

      if ( *(dtotfig+doff+j) > 0 ) {
	
	totfig=1;
	if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(derr+doff+j),
			      missing)) == -1) return iend;
	totfig=1;
	if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(doerr+doff+j),
			      missing)) == -1) return iend;

      } /*  */
      
      else {
	
	*(derr+doff+j)=0;
	*(doerr+doff+j)=0;
	*(drightfig+doff+j)=2;
	
      } 
      
      
    } /* End of loop over all data params (variables) */
    
  } /* End of loop over all levels */

  /* Read to end of station */
  while ( ( i = fgetc(fp)) != '\n' && !feof(fp) );
  return iend;
  
} /* End of oclread */

/************************************************************
 * The function spacer sets up original spacing for all
 * dynamic arrays.  The integer input, intime, has four
 * possible settings, each determining run mode of spacing:
 *
 * 1 => Initialize all dynamic arrays
 * 2 => Redimension taxa arrays
 * 3 => Redimension depth
 * 4 => Redimension measured parameter arrays
 *
 * NOTE: Looks like the only difference WOD09 to WOD13 is
 * the number of standard depth levels has shifted from
 * 40 to 137.  Otherwise, same.  I added the void specifier
 * in the function header.
 *************************************************************/

void spacer(int intime) {

  if ( intime == 1 ) {

/***************************************************************

 ALLOCATE SPACE FOR TAXA SET DATA

****************************************************************/

    ntsets=1;
    ntsetsmax=1;
    if ( (ntcode =calloc( ntsets * maxtax, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
    if ( (ntloc =calloc( ntsets * maxtax, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
    if ( (nttotfig =calloc( ntsets * maxtax, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
    if ( (ntsigfig =calloc( ntsets * maxtax, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
    if ( (ntrightfig =calloc( ntsets * maxtax, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
    if ( (ntval =calloc( ntsets * maxtax, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
    if ( (nterr =calloc( ntsets * maxtax, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
    if ( (ntoerr =calloc( ntsets * maxtax, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);

/***************************************************************

 ALLOCATE SPACE FOR DEPTH AND MEASURED PARAMETERS

****************************************************************/

    /* WOD09 had 40 for all sizes in this block.
     * WOD13 has 137 standard depth levels. */
    zsize= 137;
    zsizemax= 137;
    if ( (ztotfig =calloc(zsize, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
    if ( (zsigfig =calloc(zsize, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
    if ( (zrightfig =calloc(zsize, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
    if ( (depth =calloc( zsize, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
    if ( (zerr =calloc( zsize, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
    if ( (zoerr =calloc( zsize, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
 
    /* WOD09 had 40 for all sizes in this block.
     * WOD13 has 137 standard depth levels. */
    isize=137;
    isizemax= 137;
    if ( (dtotfig =calloc(isize, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA\n");
    if ( (dsigfig =calloc(isize, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA\n");
    if ( (drightfig =calloc(isize, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA\n");
    if ( (dataval =calloc( isize, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA\n");
    if ( (derr =calloc( isize, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA\n");
    if ( (doerr =calloc( isize, sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA\n"); 

  } /* End of initialize dynamic arrays check */

/***************************************************************

 REALLOCATE SPACE FOR TAXA DATA

****************************************************************/

  else if ( intime == 2 ) {

    ntsetsmax = ntsets;
    if ( (ntcode =realloc( ntcode, ntsets * maxtax * sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
    if ( (ntloc =realloc( ntloc, ntsets * maxtax * sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
    if ( (nttotfig =realloc( nttotfig, ntsets * maxtax * sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
    if ( (ntsigfig =realloc( ntsigfig, ntsets * maxtax * sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
    if ( (ntrightfig =realloc( ntrightfig, ntsets * maxtax * sizeof(int)))
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
    if ( (ntval =realloc( ntval, ntsets * maxtax * sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
    if ( (nterr =realloc( nterr, ntsets * maxtax * sizeof(int)) )
         == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
    if ( (ntoerr =realloc( ntoerr, ntsets * maxtax * sizeof(int)) )
         == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
    
  } /* End of taxa data reallocation */

/***********************************************************

 REALLOCATE SPACE FOR DEPTH

************************************************************/

  else if ( intime == 3 ) {

    if ( (ztotfig =realloc( ztotfig, zsize * sizeof(int)) )
	 == NULL )
      printf( " #1 NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
    if ( (zsigfig =realloc( zsigfig, zsize * sizeof(int)) )
	 == NULL )
      printf( " #2 NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
    if ( (zrightfig =realloc( zrightfig, zsize * sizeof(int)) )
	 == NULL )
      printf( " #3 NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
    if ( (depth =realloc( depth, zsize * sizeof(int)) )
	 == NULL )
      printf( " #4 NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
    if ( (zerr =realloc( zerr, zsize * sizeof(int)) )
       == NULL )
      printf( " #5 NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
    if ( (zoerr =realloc( zerr, zsize * sizeof(int)) )
	 == NULL )
      printf( " #6 NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
    
  } /* End of depth data reallocation */

/************************************************************

 REALLOCATE SPACE FOR MEASURED PARAMETERS

*************************************************************/

  else if ( intime == 4 ) {

    if ( (dtotfig =realloc( dtotfig, isize * sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA");
    if ( (dsigfig =realloc( dsigfig, isize * sizeof(int)) )
       == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA");
    if ( (drightfig =realloc( drightfig, isize * sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA");
    if ( (dataval =realloc( dataval, isize * sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA");
    if ( (derr =realloc( derr, isize * sizeof(int)) )
       == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA");
    if ( (doerr =realloc( derr, isize * sizeof(int)) )
	 == NULL )
      printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA");
    
  } /* End of measured parameters reallocation */
  
} /* End of spacer function */

/**********************************************************
 * The function extracti extracts an integer value from
 * a file.  There are three types of extractions:
 *
 * (0) Extracting the number of bytes in the integer value
 *     and then using that information to get the value
 *     itself.  (Internal integer value)
 *
 * (1) Extract first the number of significant figures,
 *     number of total figures, and number of figures to
 *     the right of the decimal before extracting the
 *     actual value.  If the number of sig figs is a
 *     negative sign, the data value is a missing value.
 *     (Output data value)
 *
 * (2) Extracting a value w/ input number of figures (totfig)
 *     (Output data value from given # of figures)
 *
 * The inputs are:
 *
 * int type      - Defines the type of extraction <0|1|2>
 *
 * int *totfig   - Number of figs. in ASCII representation
 *                 of value being extracted.
 *
 * int *sigfig   - Number of sig fig in value being extracted
 *
 * int *rightfig - Number of places to the right of the
 *                 decimal in value being extracted.
 *
 * int *value    - Value being extracted, as an integer
 *
 * int missing   - Missing value marker
 *
 * If the function encounters EOF, -1 is returned.
 * Otherwise, zero is returned for sucessful extration.
 *
 * NOTE: IDENTICAL WOD09 and WOD13 except for int
 * declaration in function header.
 **********************************************************/

int extracti(int type,int *totfig,int *sigfig,
             int *rightfig,int *value,int missing) {

  int sign,j,i;

  /* Skip if end of line character */
  if ( type != 2 ) i=nocrfgetc();

  /* If this is end of file, set EOF and notify main */
  if ( feof(fp) ) return -1;

  else {
    /* If this is a missing value (I = '-'),
     * set value accordingly. */
    if ( i == '-' ) {
      *value = missing;
      *totfig = 0;
      *sigfig = 0;
      *rightfig = 0;
      return 0;
    } /* End missing value check */
    
    else {

      /* Read in and/or set number of sig figs,
       * total number of figs, # figs right of
       * decimal if type == 1.  If type 0 (internal
       * integer) just set totfig. */
      if ( type == 1 ) {

	*sigfig = i - '0';

	i = nocrfgetc();
	*totfig= i - '0';
	
	i = nocrfgetc();
	*rightfig= i - '0';

      } /* End of type 1 check */

      else if ( type == 0 ) *totfig= i - '0'; 

      /* Read in value, including sign */
      
      *value= 0;  /* value is init to zero */
      for ( j = 1; j <= *totfig; j++ ) {

	i = nocrfgetc(); /* Read in each additional figure */
	
	/* Multiply the previously aquired figures
         * by ten and add the newest figure, i.e.
         * 123 = 12*10 + 3 */
	if ( j > 1 ) *value= 10 * *value + ( i - '0' );
	else {
	  /* Check for the sign of the value */
	  sign = (i == '-') ? -1 : 1;
	  if ( sign == 1 && i != ' ') *value = (i - '0');
	} /* End of sign check */

      } /* End of reading in value, figure by figure */

      /* Apply the sign to the output */
      *value *= sign;

    } /* End of missing/valid value check */
    
  } /* End of EOF check */

  return 0;
  
} /* End of function extracti */

/*************************************************************
 * The function extractc extracts character data from ocl
 * files.  This function runs in two modes determined by
 * the type input variable.  The inputs are:
 *
 * int type    - type = 0 if totfig is given, 1 if totfig tbd.
 *
 * int *totfig - number of bytes in character data
 *
 * char *cdata - array of character data
 *
 * When an EOF is encountered, -1 is returned.
 * Zero is returned otherwise after a sucessful read.
 *
 * NOTE: NO DIFFERENCE BETWEEN WOD09 and WOD13
 * VERSIONS in wodC.c except for the int declaration
 * in the function header, below.
 *************************************************************/

int extractc(int type, int *totfig, char *cdata) {

 int i,j;

 /* Check whether or not we are at EOF
  * This is done here to check that the
  * value for the number of figures that
  * is read in is a valid number and not
  * and EOF marker.*/
 if ( type !=0 ) i = nocrfgetc();
 
 if ( feof(fp) ) return -1;

 else {
   /* Compute the number of char */
   if ( type !=0 ) *totfig = i - '0';

   /* Read in each character */
   for ( j = 0; j < *totfig; j++ ) {

     i = nocrfgetc();
     *(cdata+j) = i;
     
   } /* End of reading each character */
   
 } /* End of EOF check */
 
 /* If EOF was encountered during the reading */
 if ( feof(fp) ) return -1;

 else return 0;

} /* End of function extractc */

/***********************************************
 * The function nocrfgetc (no carriage return
 * fgetc) reads ocl format skipping end of line
 * characters in a way which will work in PC
 * and UNIX environment.  The PC end of line
 * character (^M) is followed by the \n character
 * on some UNIX platforms.
 *
 * Uses: fp - global variable file id
 * 
 * Returns input character or -1 for EOF
 *
 * NOTE: No difference WOD09 and WOD13 except
 * for the int declaration in front of the
 * function name.  Retained here.
 ***********************************************/

int nocrfgetc() {

 int i;

 while ( !feof(fp) && isprint ( (i=fgetc(fp)) ) == 0 );

 if ( feof(fp) ) i = -1;

 return i;

} /* End function nocrfgetc */
