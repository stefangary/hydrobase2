/*  sgnc_convert.c
----------------------------------------------------------------------------------
*   Reads SeaGlider netCDF files created by the basestation post-processing
*   software and converts the files into HB2 format.
*
*   extracts :    header info
*                 p,d,t,s,ox (no fluorometer data!)
*   
*   USAGE: sgnc_convert filelist [-B<badfile>] [-D<dirname>]
*                                [-E<extent>] [-h] [-O<outfile>]
----------------------------------------------------------------------------------
* Based on the shell of wod09_convert.c because that
* conversion tool is the one SFG is most familiar with.
* Stefan Gary, October, 2014.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hydrobase.h"
#include "netcdf.h"

#define   BLANK     0x20   /* ascii code for blank */
#define   MISSING   -9.9   /* missing value flag */
#define   DIR    ""
#define   EXTENT ""

/*******************************************************************

 Definitions and global vars.  These variables are global because
 in C anything definied outside of a code block is global.

********************************************************************/

/*
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
int levels, isoor, nparm;
int htotfig[3], hsigfig[3], hrightfig[3];

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

 Internal arrays

 tenp   - Powers of ten (double array in wodC.c, 9 entries)
 sdepth - Standard depth levels

*************************************************************/

float tenp[] = { 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000 };
int sdepth[] = { 0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250,
                  300, 400, 500, 600, 700, 800, 900, 1000, 1100,
                  1200, 1300, 1400, 1500, 1750, 2000, 2500, 3000,
                  3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000,
                  7500, 8000, 8500, 9000 };

/*************************************************************
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
int read_and_check_sta(int, int);
int check_profile_flag(int);
int check_obs_flag(int);
void wipe_good_and_bad_sta();
int cdf_get_global_att_int(int, char *);
/*char *cdf_get_global_att_str(int, char *); for date */
double cdf_get_var_double(int, char *, int);

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
  int error, bopt, verbose;

   /* File ID of output file and rejected values file */
   int outfile, badfile;

   /* Count stations read, output (out), and rejected (bad)
    * for all stations, for dives, and for climbs separately. */
   int staread, staout, stabad;
   int divread, divout, divbad;
   int cliread, cliout, clibad;

   /* Integer locations within netcdf file to read dive
    * and climb separate from each other. */
   int dive_beg, dive_end;
   int climb_beg, climb_end;
   int dive_levels, climb_levels;

   /* Local counters */
   int  i, j, ii, curfile = 1, nfiles = 0;

   char *dir, *extent;
   char st[80]; /* Buffer to become file name */

   int ncid, nc_status;
   int glider;
   int dive_number;
   int mission;
   
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
   divout = divread = divbad = 0;
   cliout = cliread = clibad = 0;
   outfile = STDOUT;

/* Set check values */
   bopt = 0;
   verbose = 0;

/* Nullify all data in the data storage structures */
   for (i = 0; i < MAXPROP; ++i) {
     data.observ[i] = NULL;
     bad_data.observ[i] = NULL;
   }
   hdr.prop_id = NULL;
   bad_hdr.prop_id = NULL;

/**********************************************
 * Default parameters for all SeaGlider files.*
 **********************************************/

   /* WORKING HERE; move to autodetect within check_and_read_sta()
    * Need to autodetect because could change ncfile to dive file
    * if mixing multiple glider missions together. */
/* Currently, only T and S available so
 * nprops is always 2.  In the future,
 * may wish to make a command line
 * optional flag to include O2 and fluor.*/
   hdr.nprops = 2;
   bad_hdr.nprops = 2;

   data.nprops = 2;
   bad_data.nprops = 2;

/**********************************************/
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
         
       case 'h':  
	 print_usage(argv[0]);
	 exit(0);

       case 'v':
	 verbose = 1;
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
   
   if (! bopt) {
       fprintf(stderr,"\nBad stations will not be saved in a separate file\n");
   }
   
   if (! nfiles) {
       fprintf(stderr,"\n ERROR: No file names specified as input!");
       fprintf(stderr,"\n NetCDF cannot read from stdin, so user");
       fprintf(stderr,"\n must specify file names explicitly!");
       exit(1);
   }
 
  /* Loop for each file.  Variable curfile is init
   * to 1 in declarations. It is assumed that each
   * file has only one seaglider dive in it. 
   */
  do {
    if (nfiles) {

      /*========================================================*/
      /* Open the netcdf file and check validity.*/
      ncid = cdf_open(dir,argv[curfile],extent,1);

      /* Check for unopenable files and print error.  Will
       * continue to operate on remaining valid files.*/
      if ( ncid == -1 ) {
	if ( verbose ) {
	  fprintf(stderr,"\n Skipping %s because cannot open file.",argv[curfile]);
	}
	goto NEXTFILE;
      } /* End of null file check */

      /*========================================================*/
      /* Load basic header information that doesn't change for
       * dive or climb. */

      /* The SeaGlider S/N and mission number are
       * stored in the 5 digit cruise number as
       * IIIMM (first three digits are SG ID and
       * last two digits are mission number).
       * These numbers are determined from the
       * netCDF file automatically, below.
       * Since there is the potential to clash
       * with existing ship names and country
       * codes, just keep to XX values.*/
      strncpy(hdr.ship, "XX\0", 3);
      strncpy(bad_hdr.ship, "XX\0", 3);

      hdr.origin = '1';
      bad_hdr.origin = '1';

      strncpy(hdr.country, "XX\0", 3); 
      strncpy(bad_hdr.country, "XX\0", 3);

      if ( verbose ) {
	fprintf(stderr,"\n Country set to %s",hdr.country);
	fprintf(stderr,"\n Ship set to %s",hdr.ship);
	fprintf(stderr,"\n Origin set to %c",hdr.origin);
      }

      /* Get the Seaglider ID.*/
      glider = cdf_get_global_att_int(ncid,"glider");
      if ( verbose ) {fprintf(stderr,"\n Read glider ID %d",glider);}

      /* Get the  mission number.*/
      mission = cdf_get_global_att_int(ncid,"mission");
      if ( verbose ) {fprintf(stderr,"\n Read mission %d",mission);}

      /* Combine glider ID and mission number into the cruise number,
       * IIIMM (first three digits glider ID, last two are mission number.*/
      hdr.cruise = glider*100 + mission;
      bad_hdr.cruise = glider*100 + mission;
      if (verbose) {fprintf(stderr,"\n Cruise number set to %d",hdr.cruise);}

      /* Get the dive number.*/
      dive_number = cdf_get_global_att_int(ncid,"dive_number");
      if (verbose) {fprintf(stderr,"\n Read dive_number %d",dive_number);}
      hdr.station = dive_number;
      bad_hdr.station = dive_number;
      if (verbose) {fprintf(stderr,"\n Station number set to %d",hdr.station);}

      /*========================================================*/
      /* Primary goal here is to determine reading
       * interval for the climb. No specific values
       * state the end of the dive because the dive
       * is ended only when pumping begins.*/

      /* Get the first value in the time vector - start of the dive.
       * Units of seconds since 00:00Z GMT 1970.*/

      /* Get the start_of_climb_time and add it to the start of the
       * of the dive time to get the time in seconds since 00:00Z GMT
       * 1970 of when the climb starts.*/

      /* Load the whole ctd_time vector, in units of seconds since
       * 00:00Z GMT 1970.*/

      /* Scan through the ctd_time vector to find the index for which
       * the climb started.*/

      /*========================================================
       * Set basic parameters for this dive (instrument code,
       * lon, lat, date -> these last three based on GPS2, while
       * the lon, lat, and date for climb is GPSE (below).
       =========================================================*/
      if(verbose){fprintf(stderr,"\n======Now focus on dive======");}
      hdr.instrument = 'g';
      bad_hdr.instrument = 'g';
      if(verbose){fprintf(stderr,"\n Instru. code to %c",hdr.instrument);}

      /* Get GPS2, the GPS value obtained immediately before the dive.*/
      /* log_gps_lon
         log_gps_lat
	 are preparsed double arrays, that have already converted the
	 degrees, minutes, decimal minutes in degrees.decimal degrees.

	 GPS1 is [1] (don't need this one)
	 GPS2 is [2] (pos for dive)
	 GPSE is [3] (pos for climb)

	 note these are variables, not attributes.

	 one can use the GPS1_qc, GPS2_qc, and GPSE_qc variables to
	 precheck for valid locations.

	 CTD_qc is an overall flag for quality of the profile.

	 As for the date (m, d, y and not seconds since 1970), the
	 quick way is to look at the time_coverage_min and max values
	 in the global attributes.

      /* Now that we know when the climb begins, we can call a
       * single function that will get information from the dive
       * and then get information from the climb based on the
       * given index points.
       * about the position by finding the average position during
       * the dive/climb.*/
   /*hdr.lat = (float) latitude / tenp[*(hrightfig+1)];
   bad_hdr.lat = hdr.lat;

   hdr.lon = longitude/ tenp[ *(hrightfig+2)];
   bad_hdr.lon = hdr.lon;

   hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
   bad_hdr.ms10 = hdr.ms10;

   bad_hdr.ms1 = hdr.ms1;*/

      /*
      Notes about useful variables:

start_of_climb_time - time in seconds since the 

      */
      dive_beg = 0;
      dive_end = 0;
      dive_levels = 0;

      /*========================================================*/
      /* Read the dive into hydrobase database variables. */

      /* HYDRO_HDR contains:
       char instrument
       int  year, month, day
       float lat, lon
       int  pdr, nobs, nprops, ms10, ms1
       int *prop_id
       char qual[nqual]

The structure HYDRO_DATA contains:
       int nobs, nprops
       double *observ[MAXPROP]
      */
      /*
      prop_id
      qual
      */

      /* Update instrument code for the dive. */


      /* Read the dive */
      dive_levels = read_and_check_sta(dive_beg,dive_end);

      /* Write dive to output. */
      ++divread;
      ++staread;

      if (dive_levels > 0) {
	/* The station is good, so write entire station to output */

	/* Calculate the avg. lon, lat, timestamp for dive. */

	write_hydro_station(outfile, &hdr, &data);

	/* Augment output station counter */
	++staout;
	++divout;
      }

      if (bopt && (bad_hdr.nobs > 0)) {
	/* If we have -B<badfile> enabled, bopt = 1
	 * and if the profile has more than zero
	 * observations, write this profile to
	 * rejected station file. */
	write_hydro_station(badfile, &bad_hdr, &bad_data);
      
	/* Augment the rejected station counter */
	++stabad;
	++divbad;
      }

      /* Clear memory in prep for climb. */
      wipe_good_and_bad_sta();

      /*========================================================*/
      /* Read the climb into hydrobase database variables. */
      /* WORKING HERE - MIRROR FINAL FORMAT FOR DIVE */
      cliread++;
      staread++;

      climb_beg = 0;
      climb_end = 0;
      climb_levels = 0;

      /* Calculate the avg. lon, lat, timestamp for climb. */
      /* Update the header for climb. */
      if(verbose){fprintf(stderr,"\n======Now focus on climb======");}
      hdr.instrument = 'h';
      bad_hdr.instrument = 'h';
      if(verbose){fprintf(stderr,"\n Instru. code to %c",hdr.instrument);}

      /* Write climb to output. */

      /* Clear memory in prep for next file (if avail). */
      wipe_good_and_bad_sta();
      /*========================================================*/

    }
    
  NEXTFILE:
    if (nfiles) {
      /* Close netcdf file. */
      cdf_close(ncid);
    }
  } while (curfile++ < nfiles);  /* End of while looping over files */
  
  fprintf(stderr,"\nEnd of conversion.");
  fprintf(stderr,"\n  %d stations read in", staread);
  fprintf(stderr,"\n  %d dives read in", divread);
  fprintf(stderr,"\n  %d climbs read in", cliread);

  fprintf(stderr,"\n  %d stations accepted", staout);
  fprintf(stderr,"\n  %d dives accepted", divout);
  fprintf(stderr,"\n  %d climbs accepted", cliout);

  fprintf(stderr,"\n  %d stations rejected ", staread-staout);
  fprintf(stderr,"\n  %d dives rejected ", divread-divout);
  fprintf(stderr,"\n  %d climbs rejected ", cliread-cliout);

  fprintf(stderr,"\n  %d stations contained bad scans \n\n", stabad - (staread-staout));
  fprintf(stderr,"\n  %d dives contained bad scans \n\n", divbad - (divread-divout));
  fprintf(stderr,"\n  %d climbs contained bad scans \n\n", clibad - (cliread-cliout));
  exit(0);
  
} /* end main() */ 
  
/*****************************************************************************/
void wipe_good_and_bad_sta() {
  /* Clean up -> clear all pointers in both good
   * and bad header/data information. The header
   * and data structures are both global vars. */

  int i;

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
}

/*****************************************************************************/
void print_usage(char *program) {

   fprintf(stderr,"\nConverts SeaGlider netCDF files from basestation 2.08 to HydroBase\n");
   fprintf(stderr,"\nUsage: %s filelist [-B<badfile>] [-D<dirname>]",program);
   fprintf(stderr,"\n            [-E<extent>] [-h] [-O<outfile>] [-v]");
   fprintf(stderr,"\n");
   fprintf(stderr,"\n  List of filenames must be first argument!");
   fprintf(stderr,"\n  Input cannot be expected from stdin because netcdf");
   fprintf(stderr,"\n  must open each netcdf file. Routine loads T, S.");
   fprintf(stderr,"\n  Oxygen and Fluorescence has yet to be implemented.");
   fprintf(stderr,"\n");
   fprintf(stderr,"\n  This routine will read the dive files, e.g. p*.nc,");
   fprintf(stderr,"\n  and DOES NOT read the mission_profile or");
   fprintf(stderr,"\n  mission_timeseries files.  This is because the");
   fprintf(stderr,"\n  dive files have QC flags and the ultimate goal is");
   fprintf(stderr,"\n  to port the data, in as raw a format as possible,");
   fprintf(stderr,"\n  into the other HB tools for glider-ship inter-QC.");
   fprintf(stderr,"\n\n   OPTIONS: ");
   fprintf(stderr,"\n-B : Specify output file for rejected stations (noclobber)");
   fprintf(stderr,"\n          Defaults to stdout.");
   fprintf(stderr,"\n-D : dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n-E : input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.nc ");
   fprintf(stderr,"\n-O : Specify output file.  Will overwrite an");
   fprintf(stderr,"\n            existing file with same name.");
   fprintf(stderr,"\n            Not specifying this option sends");
   fprintf(stderr,"\n            results to stdout.");
   fprintf(stderr,"\n-h : help...... prints this message. ");
   fprintf(stderr,"\n-v : verbose output gives lots of extra info,");
   fprintf(stderr,"\n     mostly for development.");
   fprintf(stderr,"\n\n   ASSUMPTIONS: ");  
   fprintf(stderr,"\n  1) There is only one dive-climb pair for each dive file.");
   fprintf(stderr,"\n     The dive and the climb are separated into two");
   fprintf(stderr,"\n     HB-format stations within the same output file.");
   fprintf(stderr,"\n  2) Longitude and Latitude from GPS only (not the flight");
   fprintf(stderr,"\n     model) will be used to determine the dive and climb");
   fprintf(stderr,"\n     positions simply by assuming that the midpoint of the");
   fprintf(stderr,"\n     dive is at the mean of GPS2 and GPSE, etc.");
   fprintf(stderr,"\n  3) Same approach (in 2) applies to the date and time.");
   fprintf(stderr,"\n  4) Instrument code is g for dive and h for climb.");
   fprintf(stderr,"\n  5) XX and XX will be sent to ship and country codes.");
   fprintf(stderr,"\n  6) ...but one can still track down the glider and mission");
   fprintf(stderr,"\n     because the cruise number will be IIIMM, where III");
   fprintf(stderr,"\n     is the 3 digit cruise number and MM is the 2 digit");
   fprintf(stderr,"\n     (and zero-padded) mission number.  The cruise num");
   fprintf(stderr,"\n     is not zero-padded.");
   return;
}
   
/****************************************************************************
 * The function cdf_get_global_att_int will read the value of the global
 * attribute for the given netcdf file ID and attribute name.
 * Stefan Gary, Feb. 2015
 ****************************************************************************/
int cdf_get_global_att_int(int cdf_file_id, char *attribute_name) {
  /* Assume CDF file already opened,
   * need to find length of attribute,
   * then allocate space,
   * and then finally read the value.
   * Returns a -1 if there was an error along the way.*/

  int nc_status;
  int attribute_length;
  int *attribute_value;

  nc_status = nc_inq_attlen(cdf_file_id,NC_GLOBAL,attribute_name,&attribute_length);
  if ( nc_status != NC_NOERR ) {
    fprintf(stderr,"\n WARNING: Cannot read attribute length!");
    return(-1);
  }
  attribute_value = (int *) malloc(attribute_length * sizeof(int));
  nc_status = nc_get_att_int(cdf_file_id,NC_GLOBAL,attribute_name,attribute_value);
  if ( nc_status != NC_NOERR ) {
    fprintf(stderr,"\n WARNING: Cannot read attribute value!");
    return(-1);
  }
  return(*attribute_value);
} /* end of cdf_get_global_att_int */

/****************************************************************************
 * The function cdf_get_var_double will read the value of a single
 * double value in the given netcdf file (already open) and variable
 * name and the index value within the variable.
 * Stefan Gary, Feb. 2015
 ****************************************************************************/
double cdf_get_var_double(int cdf_file_id, char *var_name, int index) {
  /* Assume CDF file already opened,
   * need to find length of attribute,
   * then allocate space,
   * and then finally read the value.
   * Returns a -1 if there was an error along the way.
   * WORKING HERE

   */

  int nc_status;
  int attribute_length;
  int *attribute_value;

  nc_status = nc_inq_attlen(cdf_file_id,NC_GLOBAL,var_name,&attribute_length);
  if ( nc_status != NC_NOERR ) {
    fprintf(stderr,"\n WARNING: Cannot read attribute length!");
    return(-1);
  }
  attribute_value = (int *) malloc(attribute_length * sizeof(int));
  nc_status = nc_get_att_int(cdf_file_id,NC_GLOBAL,var_name,attribute_value);
  if ( nc_status != NC_NOERR ) {
    fprintf(stderr,"\n WARNING: Cannot read attribute value!");
    return(-1);
  }
  return(*attribute_value);
} /* end of cdf_get_var_double */
/****************************************************************************
 * The function read_and_check_sta returns the number of observations
 * in the station (levels) if the station has acceptable data, 0 if not.
 * 
 * Basically, this routine copies all the data that was
 * loaded by oclread() into the global WOD-formatted arrays
 * into the global hydrobase structures.  A direct copy is
 * put into the bad structure and an edited (only good profiles
 * and scans retained) is copied to the good structure.
 ****************************************************************************/

int read_and_check_sta(int start_read, int end_read) {

  /* Arrays for keeping track of the presence and
   * status/flag of each variable in the profile */
  int *prop_avail, *prop_OK;

  /* Storage of index of property within arrays */
  int *dindex;

  int **data_flagged;
  int offset, i, j, nobs, delta;
  int found, index;
  float prsint;
   
   /* Allocate some space */
   prop_avail = (int *) calloc(MAXPROP, sizeof(int));
   prop_OK = (int *) calloc(MAXPROP, sizeof(int));
   dindex = (int *) calloc(MAXPROP, sizeof(int));
   data_flagged = (int **) calloc(MAXPROP, sizeof(int *));
   
   /* Store metadata in HydroBase station header.
    * Note that we store duplicate information in
    * both the good header and the bad header. */      



   /*
   hdr.year = year;
   bad_hdr.year = year;

   hdr.month = month;
   bad_hdr.month = month;

   hdr.day = day;
   bad_hdr.day = day;
   */
   /* Compute default values for cruise and station # */
   /*hdr.cruise =  icruise - ((icruise/100000) * 100000);
   hdr.station = ostation - ((ostation/10000) * 10000);
   */
   /* Use originator's cruise and station #s if available.
    * Start by setting end of string delimiters. */
   /**(origc+origcfig) = '\0';
   *(origs+origsfig) = '\0';
   */
   /* If there are figures for the cruise #, set cruise #.*/
   /*if (origcfig > 0) {      
      if ((i = atoi(origc)) > 0)
         hdr.cruise = i - ((i/100000) * 100000);
   }
   */
   /* If there are figures for the station #, set station #.*/
   /*if (origsfig > 0) {
      if ((i = atoi(origs)) > 0)
         hdr.station = i - ((i/10000) * 10000);
   }
   bad_hdr.cruise = hdr.cruise;
   bad_hdr.station = hdr.station;
   */
	/* We have a valid ship code
	if (index > 0 && (nodc_ship[index][2] != '\0')) {
	  hdr.ship[0] = nodc_ship[index][2]; 
	  hdr.ship[1] = nodc_ship[index][3]; 
	  hdr.ship[2] = '\0';
	}
	*/

   /* Copy ship code to both good and bad struct */
   /*strncpy(bad_hdr.ship, hdr.ship, 3);*/
   
   /*Search for a PDR depth in secondary header arrays  */
   found = 0;
   i = 0;

   /*   while (!found && (i < nsec)) {*/
   /*      found = seccode[i] == 10;  */  /* secondary header code for seafloor depth */
   /*   ++i;
	}*/
   /* Default depth */
   /*   hdr.pdr = bad_hdr.pdr = 0;   
   if (found) {
      --i;   
      hdr.pdr = secval[i] / tenp[srightfig[i]]; 
      bad_hdr.pdr = hdr.pdr; 
   }
   */ 

   /* Also, perform a check on the profile flags.
    * check_profile_flag will return 1 for any
    * standard WOD05 flag - so it looks like all
    * profiles will be included.  The real data
    * checking is happening at each level with
    * the data flags. */

   /* TE, SA, OX, ... etc. are declared as an
    * enum property in hydorbase.h */
   /* Ideally auto detect what data the glider has,
    * but for now, just do it manually */
   /* The netcdf file does not provide QC flags for pressure and
    * depth.  
   prop_avail[(int)TE] = 1;	         
   prop_OK[(int)TE] = 1;
   dindex[(int)TE] = i;

   prop_avail[(int)SA] =1;
   /*prop_OK[(int)SA] = check_profile_flag(netcdf_file.CTD_qc);*/
   dindex[(int)SA] = i;
   
   prop_avail[(int)PR] = 1;
   prop_OK[(int)PR] = 1;
   dindex[(int)PR] = i;
      
   prop_avail[(int)DE] = 1;
   prop_OK[(int)DE] = 1;
   dindex[(int)DE] = i;

   /* Load property data into HydroBase station structures */
   hdr.nprops = 0;	
   bad_hdr.nprops = 0;	
   
   for (i = 0; i < MAXPROP; ++i) {
     if (prop_avail[i] ) {

       /* Allocate space */
       data.observ[i] = (double *) malloc(levels * sizeof(double));
       bad_data.observ[i] = (double *) malloc(levels * sizeof(double));
       data_flagged[i] = (int *) calloc(levels, sizeof(int));

       /* Copy data from dataval (read from infile) to hydrobase data */
       for (j=0; j < levels; ++j) {

	 /* Offset for counting between columns in
          * an array - pointer arithmetic */
	 offset = dindex[i]*levels;

	 /* Convert to real values using figs wrt decimal */
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

       /* Augment number of profiles to bad struct b/c
        * bad struct incrementally stores all data,
        * regardless of good or bad. */
       ++bad_hdr.nprops;

       if ( prop_OK[i])
	 /* Add number of valid profiles to good struct */
	 ++hdr.nprops;
     }
   } /* End of loop over all properties (variables) to load data */
	   
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
   
   /* Create space for property IDs - previous prop_id
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
   
   /* Deallocate/free up memory in preparation for the next station */
   free((void *)prop_avail);
   free((void *)prop_OK);
   free((void *)dindex);
   free((void *)data_flagged);
   
   /* Return the number of good observations */
   return (hdr.nobs) ;
   
}  /* End read_and_check_sta() */

/****************************************************************************
 * This function reflects the QC flags defined by the basestation 2.08
 * software, which is also consistent with the Argo project QC flags.
 * There is a potential to use the profile flag separately from the
 * individual data flags for a single pt., hence the two different functions
 * check_profile_flag and check_obs_flag.
 *
 * Flag = 0 = No QC performed
 * Flag = 1 = Value is OK
 * Flag = 2 = Value is probably good
 * Flag = 3 = Value is probably bad, potentially correctable
 * Flag = 4 = Untrustworthy and uncorrectable
 * Flag = 5 = Explicit manual change
 * Flag = 6 = Explicitly not sampled (versus flag 9, below)
 * Flag = 7 = Not used.
 * Flag = 8 = Interpolated value
 * Flag = 9 = Value missing, instrument timed out
 ****************************************************************************/
int check_profile_flag(int flag) {

  switch (flag) {
   
  case 0:   /* fall through */
  case 1:   
  case 2:
  case 3:
  case 5:
  case 8:
    return 1;
  default:
    return 0;
  } 
  
} /* end check_profile_flag() */

/****************************************************************************
 * The function check_obs_flag will return 1 for integer inputs listed
 * in the case n: lines and returns zero for any other integer input.
 *
 * This function reflects the QC flags defined by the basestation 2.08
 * software, which is also consistent with the Argo project QC flags.
 * Flag = 0 = No QC performed
 * Flag = 1 = Value is OK
 * Flag = 2 = Value is probably good
 * Flag = 3 = Value is probably bad, potentially correctable
 * Flag = 4 = Untrustworthy and uncorrectable
 * Flag = 5 = Explicit manual change
 * Flag = 6 = Explicitly not sampled (versus flag 9, below)
 * Flag = 7 = Not used.
 * Flag = 8 = Interpolated value
 * Flag = 9 = Value missing, instrument timed out
 *
 * Ultimately, the flags listed in this routine are "flags to ignore."
 * The range check and bullseye check are the flags we're watching out for.
 *
 ***************************************************************************/
int check_obs_flag(int flag) {

   switch (flag) {
   
   case 0:   /* fall through */
   case 1:
   case 2:
   case 3:
   case 5:
   case 8:
     return 1;
   default:
     return 0;
   } 
   
} /* end check_obs_flag() */
