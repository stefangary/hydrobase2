/*  ctd01nc_convert.c
................................................................................
.   Reads CTD NetCDF profile format files
.   extracts :    header info
.                 p,t,s
.   
.   USAGE: ctd01nc_convert infile_list  [-B<badfile>] [-D<dir>] [-E<extent>] 
                        [-O<outfile>]  [-V]
...............................................................................
*/

#define     LENNAME     1000
#define     M_LEVELS   6000
#define     MPROPS  5  /* pr de te sa ox */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "netcdf.h"
#include "hydrobase.h"


  /*  HydroBase globally defined variables */
  
struct HYDRO_HDR  hdr;
struct HYDRO_DATA data;

FILE   *listfile;
char   *buffer;
char   theFile[LENNAME];
int    fdcdf;
int    lverbose;
int    luse_list;
int    npropsout;
int    bopt, topt;
int    outfile, badfile, tonly_file;
int    i, j,  staread, staout, stabad, sta_tonly;

  /*  prototypes for locally defined functions */
void   print_usage (char *);
FILE   *openfile (char *, char *, char *);
int    open_ncfile (char *, char *, char *);
int    ctd01nc_read (void);
void   check_sta (int, int);
int    get_param_ad (char *, int, int, double *);

main (int  argc, char **argv)
{

   int    error, nfiles, curfile, status;
   char   *dir, *extent;
   char   *st;
   char   listfilename[15];
   FILE   *fp;

/* check for command line arguments */

   if (argc < 2) {
     print_usage(argv[0]);
     exit(1);
   }

/*  set these default values */

    bopt = topt = 0;
    buffer = NULL;
    curfile = 1;
    dir = "";
    error = 0;
    extent = "";
    strcpy (listfilename, "CTD_FILE_LIST\0");
    luse_list = 0;
    lverbose = 0;
    nfiles = 0;
    outfile = STDOUT;
    staout = staread = stabad = sta_tonly = 0;

     /* Set the maximum property value, as define for the system 
        in MAXPROP in the include file, and the maximum number for this
        application.
      */
    npropsout = MPROPS;
    for (i = 0; i < MAXPROP; ++i) {
       data.observ[i] = NULL;
       }
    hdr.prop_id = NULL;
    strncpy (hdr.country, "XX", 2);
    strncpy (hdr.ship, "XX", 2);
    hdr.cruise = -1;
 

/*  parse the command line arguments */
   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
        switch (argv[i][1]) {

          case 'B':                    /* open  file for output */
            badfile = create_hydro_file(&argv[i][2], OVERWRITE);
            if (badfile < 1) {
              fprintf (stderr,"\nError opening output file: %s\n", 
                       argv[i][2]);
              exit(1);
              }
            fprintf(stderr,"\nBad data will be written to: %s", &argv[i][2]);
            bopt = 1;
            break;

          case 'C':                    /* get cruise number */
            st = &argv[i][2];
            if (*st == '/' )  {
              ++st;
              }
            error = (sscanf(st,"%d", &hdr.cruise) != 1);
            break;

          case 'D':                   /* get input dir */
            dir = &argv[i][2];
            break;

          case 'E':                    /* get file extent */
            extent = &argv[i][2];
            break;

          case 'N':                    /* get country code */
            st = &argv[i][2];
            if (*st == '/' )  {
              ++st;
              }
            error = (sscanf(st,"%2s", hdr.country) != 1);
            break;

          case 'O':                    /* get output file  */
            outfile = create_hydro_file (&argv[i][2], OVERWRITE);
            if (outfile < 1) {
              fprintf (stderr,"\nError opening output file: %s\n", &argv[i][2]);
              exit(1);
              }
            fprintf(stderr,"\nOpening or overwriting: %s", &argv[i][2]);
            break;

          case 'S':                    /* get ship code */
            st = &argv[i][2];
            if (*st == '/' )  {
              ++st;
              }
            error = (sscanf(st,"%2s", hdr.ship) != 1);
            break;

          case 'T':                    /* open  file for temp_only output */
            tonly_file = create_hydro_file(&argv[i][2], OVERWRITE);
            if (badfile < 1) {
              fprintf (stderr,"\nError opening output file: %s\n", 
                       argv[i][2]);
              exit(1);
              }
            fprintf(stderr,"\nTemp-only data will be written to: %s", &argv[i][2]);
            topt = 1;
            break;

          case 'V':                    /* set verbose/debugging.  */
            lverbose = 1;
            if (strlen (argv[i]) > 2)
              sscanf (&argv[i][2], "%ld", &lverbose);
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

        }  /* end if option flag. */

       else  {
        ++nfiles;
       }
     }  /* end for each argument. */

   
   if (! bopt) {
       fprintf(stderr,"\nBad stations will not be saved in a separate file\n");
   }
   
   if (! topt) {
       fprintf(stderr,"\nTemp-only profiles will not be saved in a separate file\n");
   }
   if (! nfiles) {
       fprintf(stderr,"\nExpecting data to come from STDIN...");
       fp = stdin;
   }

   if (strncmp (hdr.country, "XX", 2) == 0)  {
     fprintf (stderr, "\nCountry code not specified in -N option.\n");
     fprintf (stderr, "  Default XX being used.\n");
     }
   if (strncmp (hdr.ship, "XX", 2) == 0)  {
     fprintf (stderr, "\nShip code not specified in -S option.\n");
     fprintf (stderr, "  Default XX being used.\n");
     }
   if (hdr.cruise == -1)  {
     fprintf (stderr, "\nCruise number must be specified in -C option.\n");
     exit (1);
     }
   


     /* A little play time. Allow the filenames to come from
        a list, for now defined to be the file identifier 
        CTD_FILE_LIST 
        and containing one complete file spec per record.
      */
if (nfiles == 1  &&
    strncmp (&argv[1][0], "CTD_FILE_LIST\0", 14) == 0)  {
  listfile = fopen (argv[1], "r");
  if (listfile == NULL)  {
    fprintf (stderr, "\nError opening %s\n", 
             listfilename);
    exit (-1);
    }
   else  {
    fprintf (stderr, "\n   using files from: %s\n",
             listfilename);
    luse_list = 1;
    }
  }
    
 /*  loop for each input file */
  do {
  
    if (nfiles) {

          /* Check for CTD_FILE_LIST. */
      if (luse_list == 1)  {
        fgets (theFile, LENNAME, listfile); 
               /* Remove the trailing \n. */
        theFile[strlen(theFile)-1] = '\0';
        if (feof (listfile) )  {
          fclose (listfile);
          nfiles = nfiles - 1;
          break;
          }
               /* Bump the nfiles so it keeps processing the do .. while. */
        nfiles ++;
        }
          /* Or just get the current argument. */
       else  {
        strcpy (theFile, argv[curfile]);
        }
      error = open_ncfile (dir, theFile, extent);
      if (lverbose > 0)  fprintf (stderr, "\n  OPENED FILE status %ld\n", error);
      if (error != NC_NOERR)  {
        goto NEXTFILE;
        }
       }

     /* loop for each station, i.e. profile. */
     
    do  {
      status = ctd01nc_read ();
        
    if (hdr.prop_id != NULL) {
      free(hdr.prop_id);
      hdr.prop_id = NULL;
      }
	          
    } while (status != EOF );
       
NEXTFILE:
  if (nfiles)  {
    error = nc_close (fdcdf);
    if (error != NC_NOERR)  {
      fprintf (stderr, "  Error closing NetCDF file %s\n",
               error, theFile);
      exit (-1);
      }
    }
  } while   (curfile++ < nfiles);       

fprintf(stderr,"\nEnd of conversion.");
fprintf(stderr,"\n  %d files processed", nfiles);
fprintf(stderr,"\n  %d stations read in", staread);
fprintf(stderr,"\n  %d stations accepted", staout);
fprintf(stderr,"\n  %d stations rejected ", staread-staout);
fprintf(stderr,"\n  %d stations contained bad scans \n\n", stabad);
exit(0);

} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"%s converts QC'd observed profile data from a CTD \n", program);
   fprintf(stderr,"cruise into HydroBase2 station format. \n");
   fprintf(stderr,"While an attempt was made at generality the specific \n");
   fprintf(stderr,"application was for two OVIDE cruises of the THALASSA. \n");
   fprintf(stderr,"\nUSAGE:  %s infile_list -Ccruise [-O<outfile>]");
   fprintf(stderr,"\n [-B<bad_file>]  [-D<dir>] [-E<extent>] ");
   fprintf(stderr,"\n [-N<country_code>]  [-S<ship_code>] ");
   fprintf(stderr,"\n infile_list: List of filenames must be first arguments.");
   fprintf(stderr,"\n   Alternatively there may be one file whose identifier");
   fprintf(stderr,"\n   begins with the characters CTD_FILE_LIST. This file");
   fprintf(stderr,"\n   must contain the file list, one per record. They will");
   fprintf(stderr,"\n   interact with the -D and -E options in the expected");
   fprintf(stderr,"\n   manner.");
   fprintf(stderr,"\n OPTIONS:     ");
   fprintf(stderr,"\n    [-B]:  output file for bad observations ");
   fprintf(stderr,"\n    [-C]:  specifies cruise number ");
   fprintf(stderr,"\n    [-D]:  directory for input files ");
   fprintf(stderr,"\n    [-E]:  extent for input files ");
   fprintf(stderr,"\n    [-N]:  2-char country code ");
   fprintf(stderr,"\n    [-O]:  specify output file for HydroBase station files ");
   fprintf(stderr,"\n    [-S]:  2-char ship code ");
   fprintf(stderr,"\n               --or stations are output to STDOUT");
   fprintf(stderr,"\n    [-V]:  verbose flag [1 2]  default for -V: 1");
   fprintf(stderr," \n");
   return;
}
   
/*****************************************************************************/
/*****************************************************************************/
FILE *openfile(char *dir, char *root, char *extent)
{
   char st[80];
   int i;
   FILE *infile;
   
   strcpy(st, dir);
   if ((i=strlen(dir)) != 0) {
      if (dir[i-1] != '/')
         strncat(st,"/",1);
   }
   strcat(st, root);
   strcat(st, extent);
   infile = fopen(st,"r");
   if (infile == NULL)
         fprintf(stderr,"\n Unable to open %s \n", st);
   else
         fprintf(stderr,"   opened %s ... \n", st);
   
   return(infile);
   
}  /* end openfile() */
   
/*****************************************************************************/
int open_ncfile (char *dir, char *root, char *extent)
{
   char st[80];
   int i, iop;
   
   strcpy(st, dir);
   if ((i=strlen(dir)) != 0) {
      if (dir[i-1] != '/')
         strncat(st,"/",1);
   }
   strcat(st, root);
   strcat(st, extent);
   iop = nc_open (st, NC_NOWRITE, &fdcdf);
   if (iop != NC_NOERR)  {
     fprintf(stderr,"\n Unable to open %s \n", st);
     }
    else
     fprintf(stderr,"   opened %s ... \n", st);
   
   return (iop);
   
}  /* end open_ncfile() */
/*****************************************************************************/
int ctd01nc_read ()
   /*  Upon entry, buffer is NULL if first attempt to read this file.
      Reads one entire profile file. 
       with all its properties and flags into global Hydrobase structures 
      Returns status = 1 for successful read, 
                       EOF for end of file 
	An error reading or parsing causes an error message and exit from program.
 */
{
   char   theVar[30], ctem[M_LEVELS], ctem_fill[M_LEVELS];
   int    status, ctd01nc_id, ctd01nc_id2;
   int    i, ii, count, badcount, propindex;
   int    ib[3], ic[3], iop;
   int    item, item_fill;
   int    idate[3], itime[3];
   int    dimid, varid;
   int    j;
   int    k_prof;
   int    k_pr_loc, k_te_loc, k_sa_loc, k_ox_loc;
   int    last_valid;
   int    n_prof, n_levels, n_param;
   int    npts;
   int    qca[10];
   double dtem, dtem_fill, dtema[M_LEVELS];
   double juldate, juldate_ref;
     /* Nominal pressure interval and desired output interval. */
   double DELTAP, prsint;

   
/* If buffer is empty, this is first line of file.
   But now we are expecting a NetCDF file and by the CTD rules
   there should only be one profile per file. But we don't trust them
   so much. I'm usually more trusting about checking for error codes
   but here we will belabor the checking.
 */   

     /* Get the number of profiles in this file group. */
   if (buffer == NULL) {
     strcpy (theVar, "N_PROF\0");
     status = nc_inq_dimid (fdcdf, theVar, &dimid);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting ctd01nc dimension id for %s\n.", theVar);
       exit (1);
       } 
     status = nc_inq_dimlen (fdcdf, dimid, &n_prof);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting ctd01nc dimension value for %s\n.", theVar);
       exit (1);
       } 
     if (n_prof < 1)  {
       fprintf(stderr,"Error checking ctd01nc dimension value for %s\n.", theVar);
       exit (1);
       } 
     if (lverbose > 0)  fprintf (stderr, "  %s  =  %ld\n",  theVar, n_prof);


     /* Get the number of pressure levels.
        This is the maximum number of the set and used for the array 
        size so each individual station will have to be tested.
      */
     strcpy (theVar, "N_LEVELS\0");
     status = nc_inq_dimid (fdcdf, theVar, &dimid);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting ctd01nc dimension id for %s\n.", theVar);
       exit (1);
       } 
     status = nc_inq_dimlen (fdcdf, dimid, &n_levels);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting ctd01nc dimension value for %s\n.", theVar);
       exit (1);
       } 
     if (n_levels < 1)  {
       fprintf(stderr,"Error checking ctd01nc dimension value for %s\n.", theVar);
       exit (1);
       } 
     if (lverbose > 0)  fprintf (stderr, "  %s  =  %ld\n",  theVar, n_levels);
     if (n_levels >= M_LEVELS)  { 
       fprintf (stderr, "  *** AMAZING. WHO WOULD HAVE THOUGHT THERE WOULD\n");
       fprintf (stderr, "      BE MORE THAN %ld LEVELS IN A PROFILE.\n", M_LEVELS);
       exit (-1);
       }


     /* Get the number of parameters in this group. 
        This is again a maximum dimension number.
      */
     strcpy (theVar, "N_PARAM\0");
     status = nc_inq_dimid (fdcdf, theVar, &dimid);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting ctd01nc dimension id for %s\n.", theVar);
       exit (1);
       } 
     status = nc_inq_dimlen (fdcdf, dimid, &n_param);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting ctd01nc dimension value for %s\n.", theVar);
       exit (1);
       } 
     if (n_param <= 1)  {
       fprintf(stderr,"Error checking ctd01nc dimension value for %s\n.", theVar);
       exit (1);
       } 
     if (lverbose > 0)  fprintf (stderr, "  %s  =  %ld\n",  theVar, n_param);


     /* We will need the reference date for later.
        As originally explained the files are grouped by day but
        they could be ordered any way in general so don't rely on the
        file name for the date.
        For the CTD format being used this is a global attribute.
      */
     strcpy (theVar, "Reference_date_time\0");
     status = nc_get_att_text (fdcdf, NC_GLOBAL, theVar, ctem);
     ctem[14] = '\0';
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting ctd01nc variable value for %s.\n", theVar);
       exit (1);
       } 

     sscanf (ctem, "%4ld%2ld%2ld%2ld%2ld%2ld",
             &idate[0], &idate[1], &idate[2],
             &itime[0], &itime[1], &itime[2]);
     gdate2dday (idate, itime, &juldate_ref);
     if (lverbose > 0)  
       fprintf (stderr, "  %s  = %lf    %s  \n",  theVar, juldate_ref, ctem);
    
     }   /* End if file initialization. */



  /* Now we just plow through to get the necessary variables for the
     HYDROBASE header. We could be clever and make general purpose retrievals
     but the CTD format is pretty well established and there are not 
     an excessive number of parameters so just repeat the basic block.
   */

   for (k_prof=0; k_prof<n_prof; k_prof++)  {
     if (lverbose > 0)  fprintf (stderr, "\n\n  Start of Station.\n");
     ++staread;
     count = 0;  
     badcount = 0;
     qca[0] = qca[1] = qca[2] = 0;
     /* Get the locations. */
   ib[0] = k_prof; ib[1] = 0; ib[2] = 0;
   ic[0] = 1;      ic[1] = 0; ic[2] = 0;
   strcpy (theVar, "LATITUDE_BEGIN\0");
   status = nc_inq_varid (fdcdf, theVar, &varid);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting ctd01nc variable id for %s.\n", theVar);
     exit (1);
     } 
   status = nc_get_att_double (fdcdf, varid, "_FillValue", &dtem_fill);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting ctd01nc fill attribute for %s.\n", theVar);
     exit (1);
     } 
   status = nc_get_var1_double (fdcdf, varid, ib, &dtem);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting ctd01nc variable value for %s.\n", theVar);
     exit (1);
     } 
   if (dtem == dtem_fill)  {
     fprintf(stderr,"Invalid/missing variable value for %s.\n", theVar);
     exit (1);
     } 
   hdr.lat = (float) dtem;
   if (lverbose > 0)
     fprintf (stderr, "  %s: %lf   %lf\n", theVar, dtem, dtem_fill);

   strcpy (theVar, "LONGITUDE_BEGIN\0");
   status = nc_inq_varid (fdcdf, theVar, &varid);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting ctd01nc variable id for %s.\n", theVar);
     exit (1);
     } 
   status = nc_get_att_double (fdcdf, varid, "_FillValue", &dtem_fill);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting ctd01nc fill attribute for %s.\n", theVar);
     exit (1);
     } 
   status = nc_get_var1_double (fdcdf, varid, ib, &dtem);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting ctd01nc variable value for %s.\n", theVar);
     exit (1);
     } 
   if (dtem == dtem_fill)  {
     fprintf(stderr,"Invalid/missing variable value for %s.\n", theVar);
     exit (1);
     } 
   hdr.lon = (float) dtem;
   if (lverbose > 0)  
     fprintf (stderr, "  %s: %lf   %lf\n", theVar, dtem, dtem_fill);

   strcpy (theVar, "BOTTOM_DEPTH\0");
   status = nc_inq_varid (fdcdf, theVar, &varid);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting ctd01nc variable id for %s.\n", theVar);
     exit (1);
     } 
   status = nc_get_att_double (fdcdf, varid, "_FillValue", &dtem_fill);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting ctd01nc fill attribute for %s.\n", theVar);
     exit (1);
     } 
   status = nc_get_var1_double (fdcdf, varid, ib, &dtem);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting ctd01nc variable value for %s.\n", theVar);
     exit (1);
     } 
   if (dtem == dtem_fill)  {
     fprintf(stderr,"Invalid/missing variable value for %s.\n", theVar);
     exit (1);
     } 
   hdr.pdr = (int) dtem;
   if (lverbose > 0)
     fprintf (stderr, "  %s: %lf   %lf\n", theVar, dtem, dtem_fill);


     /* Play the CTD date game. */
   strcpy (theVar, "JULD_BEGIN\0");
   status = nc_inq_varid (fdcdf, theVar, &varid);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting ctd01nc variable id for %s.\n", theVar);
     exit (1);
     } 
   status = nc_get_att_double (fdcdf, varid, "_FillValue", &dtem_fill);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting ctd01nc fill attribute for %s.\n", theVar);
     exit (1);
     } 
   status = nc_get_var1_double (fdcdf, varid, ib, &dtem);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting ctd01nc variable value for %s.\n", theVar);
     exit (1);
     } 
   if (dtem == dtem_fill || dtem == 0.0)  {
     fprintf(stderr,"Invalid/missing variable value for %s.\n", theVar);
       fprintf(stderr,"[%s]\n", dtem);
       /* try to use date in last header, or make a false date */
       if (hdr.year < 1990 || hdr.year > 2020) {  
           hdr.year = 9999;
	   hdr.month = 99;
	   hdr.day = 99;
	 }
    } 
   else {
         item = dtem + juldate_ref;
         jdate (item, &hdr.day, &hdr.month, &hdr.year);
   }
   if (lverbose > 0)
     fprintf (stderr, "  %s: %ld   %ld  %ld\n", theVar, hdr.year, hdr.month, hdr.day);


     /* Map CRUISE_NAME. */
   strcpy (theVar, "CRUISE_NAME\0");
   ib[0] = k_prof; ib[1] = 0; ib[2] = 0;
   ic[0] = 1;      ic[1] = 16; ic[2] = 0;
   status = nc_inq_varid (fdcdf, theVar, &varid);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting ctd01nc variable id for %s.\n", theVar);
     exit (1);
     } 
   status = nc_get_att_text (fdcdf, varid, "_FillValue", ctem_fill);
   ctem_fill[0] = '\0';
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting ctd01nc fill attribute for %s.\n", theVar);
     exit (1);
     } 
   status = nc_get_vara_text (fdcdf, varid, ib, ic, ctem);
   ctem[16] = '\0';
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting ctd01nc variable value for %s.\n", theVar);
     exit (1);
     } 
   if (strcmp (ctem, ctem_fill) == 0)  {
     fprintf(stderr,"Invalid/missing variable value for %s.\n", theVar);
     exit (1);
     } 


     /* Map STATION_NUMBER. */
   strcpy (theVar, "STATION_NUMBER\0");
   ib[0] = k_prof; ib[1] = 0; ib[2] = 0;
   ic[0] = 1;      ic[1] = 0; ic[2] = 0;
   status = nc_inq_varid (fdcdf, theVar, &varid);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting ctd01nc variable id for %s.\n", theVar);
     exit (1);
     } 
   status = nc_get_att_int (fdcdf, varid, "_FillValue", &item_fill);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting ctd01nc fill attribute for %s.\n", theVar);
     exit (1);
     } 
   status = nc_get_var1_int (fdcdf, varid, ib, &item);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting ctd01nc variable value for %s.\n", theVar);
     exit (1);
     } 
   if (item == item_fill)  {
     fprintf(stderr,"Invalid/missing variable value for %s.\n", theVar);
     exit (1);
     } 
   hdr.station = item;
   if (hdr.station > 9999) {
     hdr.station = hdr.station % 1000;
     }
   if (lverbose > 0) 
     fprintf (stderr, "  %s: %ld   %ld\n", theVar, item, item_fill);


     /* A couple of hardwired fields. */
   status = sscanf (ctem, "%ld", &hdr.cruise);
   hdr.origin = '9';         /* Currently 9 = individual PI.  */
   hdr.instrument = 'c';     /* Currently ctd. */

   if (lverbose > 0)
     fprintf (stderr, "  SHIP/CRUISE: %s   %ld\n", hdr.ship, hdr.cruise);

     
     /* Start getting the profile data. */
     
        /* alloc space in global structures */
    
    hdr.nprops = MPROPS;
    hdr.prop_id = (int *) calloc(hdr.nprops, sizeof(int));
    hdr.prop_id[0] = (int)PR; 
    hdr.prop_id[1] = (int)DE;
    hdr.prop_id[2] = (int)TE;
    hdr.prop_id[3] = (int)SA;
    hdr.prop_id[4] = (int)OX;
     
    for (i = 0; i < hdr.nprops; ++i) {
      free_and_alloc(&data.observ[hdr.prop_id[i]], 6000);
      }

   for (i=0; i<n_levels; i++)  {
     data.observ[(int)PR][i] = HB_MISSING;
     data.observ[(int)DE][i] = HB_MISSING;
     data.observ[(int)TE][i] = HB_MISSING;
     data.observ[(int)SA][i] = HB_MISSING;
     data.observ[(int)OX][i] = HB_MISSING;
     }

   /* Get the pressures. 
      We could get the basic info at the beginning but then we would have
      to store all the variables uniquely for later use. This may
      take a bit longer but keeps it general.
    */
   strcpy (theVar, "PRES\0");
   status = get_param_ad (theVar, k_prof, n_levels, data.observ[(int)PR]);
   if (status != NC_NOERR)  exit (1);
   k_pr_loc = 1;
        /* Find the last valid pressure as the array was dimensioned to
           handle the largest in the group. Technically if it were
           just one profile we could skip this but I don't trust the
           recording mechanism. 
           We may also want to put the CTD quality checks in here as well.
           */
   last_valid = -1;
   for (i=0; i<n_levels; i++)  {
     if (data.observ[(int)PR][i] != HB_MISSING)  {
       last_valid = i;
       }
     }  /* End for loop validity check.  */
   count = last_valid + 1;


   strcpy (theVar, "TEMP\0");
   status = get_param_ad (theVar, k_prof, count, data.observ[(int)TE]);
   if (status != NC_NOERR)  exit (1);
   k_te_loc = 2;


   strcpy (theVar, "PSAL\0");
   status = get_param_ad (theVar, k_prof, count, data.observ[(int)SA]);
   if (status != NC_NOERR)  exit (1);
   k_sa_loc = 3;


   strcpy (theVar, "OXYL\0");
   status = get_param_ad (theVar, k_prof, count, data.observ[(int)OX]);
   if (status != NC_NOERR)  exit (1);
   k_ox_loc = 4;


     /* At this point we are going to decimate the data.
        Following the code given in primer_ctdconvert.
      */

   DELTAP = 1.0;
   prsint = 10.0;

   npts = NINT(prsint / DELTAP);  /*  get # of pts per interval*/
   if (npts == 0) {
     fprintf(stderr,"Bad number of points per interval: %d \n", npts);
     exit(1);
     }

   j = 0;
   for (i = 0; i < count; i += npts) {
     data.observ[(int)PR][j] = data.observ[(int)PR][i];
     data.observ[(int)DE][j] = hb_depth (data.observ[(int)PR][j], hdr.lat);
     data.observ[(int)TE][j] = data.observ[(int)TE][i];
     data.observ[(int)SA][j] = data.observ[(int)SA][i];
     data.observ[(int)OX][j] = data.observ[(int)OX][i];
     ++j;
     }

     /* add bottom observation  ... */

   if ((i - npts) != (--count)) {
     data.observ[(int)PR][j] = data.observ[(int)PR][count];
     data.observ[(int)DE][j] = hb_depth (data.observ[(int)PR][count], hdr.lat);
     data.observ[(int)TE][j] = data.observ[(int)TE][count];
     data.observ[(int)SA][j] = data.observ[(int)SA][count];
     data.observ[(int)OX][j] = data.observ[(int)OX][count];
     ++j;
     }
   count = j;


     /* Set the number of observations/levels.
        This is done in CHECK_STA.
        We also want to adjust the number of parameters. By executive decision
        we could have just PR and TE with no SA. Also by executive decision
        we could not realistically expect PR and SA with no TE so if TE
        is absent the whole profile is a wash. We have moved the actual
        number of properties out to a global variable and it is reset here 
        for the call to CHECK_STA then after the call to 
        WRITE_HYDRO_STATION it will be reset to the default. This should 
        not be necessary as the freeing and allocation of space is
        currently being done just after the file initialization sequence
        in APNC_READ but outside the loop for multiple profiles. So it
        is correctly dimensioned but we'll try to keep it safe.
      */
   if (k_sa_loc < 0)  --npropsout ;
 
   check_sta (count, badcount);

     /* All data and headers are now stored.
        Is it any good. Just to make things a bit challenging we check at
        least PRES and TEMP to see if they are okay. If they are bad then
        the profile is bad. But if they are okay and the PSAL is present 
        but bad then just omit that.
      */
	    
   if (hdr.nobs > 0)  {
          /* Station is bad no matter how we cut it. */
     if ( (qca[0]+qca[1]) > 0)  {
       if (bopt)
              write_hydro_station(badfile, &hdr, &data);
       ++stabad;
       }
          /* Station  only PRES,DEPTH */
      else if (npropsout < 3)  {
       if (bopt)
            write_hydro_station(badfile, &hdr, &data);
       ++stabad;
       }
          /* Station  has good PRES,DEPTH, TEMP. */
      else if (npropsout == 3)  {
         if (topt)
            write_hydro_station(tonly_file, &hdr, &data);
         ++sta_tonly;
       }
          /* Station good as far as it goes.  Check for salinity. */
      else if (npropsout == MPROPS )  {
       if (qca[2] == 0) {        
         write_hydro_station(outfile, &hdr, &data);
         ++staout;
         }
          /* A final adjustment to the number of nprops. */
        else  {
	   if (bopt)
                 write_hydro_station(badfile, &hdr, &data);
           ++stabad;
	   if (topt) {
              --hdr.nprops;
	      --data.nprops;
              write_hydro_station(tonly_file, &hdr, &data);
              ++sta_tonly;
	    }
         }  /* End of stations with at least good PRES, TEMP. */
       }  /* End if station has saliinity. */
     }  /* End if station has some observations. */
	   


   }  /* End of loop for multiple profiles. */


   status = EOF;

   return(EOF);
}  /* end ctd01nc_read() */
/*****************************************************************************/
/* It wasn't worth it for three but now that we are getting the adjusted
   values as well it is marginal for the double array.
   Get a NetCDF double array for a parameter. 
 */
int  get_param_ad (char *theVar, int k_prof, int count, double *theArray )

{
int     ib[3], ic[3];
int     status;
int     varid;
double  dtem_fill;

if (lverbose > 0)  
  fprintf (stderr, "  LOAD %s: %ld  %ld\n", theVar, k_prof, count);

status = nc_inq_varid (fdcdf, theVar, &varid);
if (status != NC_NOERR) {
  fprintf(stderr,"Error getting ctd01nc variable id for %s.\n", theVar);
  return (status);
  } 
status = nc_get_att_double (fdcdf, varid, "_FillValue", &dtem_fill);
if (status != NC_NOERR) {
  fprintf(stderr,"Error getting ctd01nc fill attribute for %s.\n", theVar);
  return (status);
  } 

ib[0] = k_prof;   ib[1] = 0;        ib[2] = 0;
ic[0] = 1;        ic[1] = count;    ic[2] = 0;
status = nc_get_vara_double (fdcdf, varid, ib, ic, theArray);
if (status != NC_NOERR) {
  fprintf(stderr,"Error getting ctd01nc variable value for %s.\n", theVar);
  fprintf (stderr, "  ERROR STATUS: %ld\n", status);
  return (status);
  } 
for (i=0; i<count; i++)  {
  if (theArray[i] == dtem_fill)  {
    theArray[i] = HB_MISSING;
    }
   else  {
    if (lverbose > 1)
      fprintf (stderr, "    %s: %ld  %.3lf\n", theVar, i, theArray[i]);
    }
  }  /* End for loop validity check.  */

return (NC_NOERR);
}  /* End of get_param_ad. */

/*****************************************************************************/
void check_sta (int count, int badcount)
   /* Gets station ready for output.  
      This is a little different that the original check_sta in that
      the original just looked for really wild points. Then it
      appears it wrote the good levels with appropriate header to
      the hydro file and bad levels with appropriate header to the
      bad file. This is
      using the quality code and as such all data is in the normal
      hdr and data structures. It's just a case of where we write it out.
    */
{
   int    i;
   
   hdr.nobs = count;
   
   if (count) {    /* pr, te, sa -- so far */
       
       hdr.nprops = npropsout;
       data.nprops = npropsout;
       data.nobs = hdr.nobs;
 
       hdr.qual[0] = '0';
       hdr.qual[1] = '0';
       hdr.qual[2] = '1';
       hdr.qual[3] = '1';
	  
       hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
       
       /* compute depth from pressure explicitly */
       for (i = 0; i < hdr.nobs; ++i)  {
          data.observ[(int)DE][i]  = HB_MISSING;
	  if (data.observ[(int)PR][i] >= 0)
             data.observ[(int)DE][i] = hb_depth(data.observ[(int)PR][i], (double)hdr.lat);
       }
            
   }  

   return;

} /* end check_sta */


