/*  hb_propcalc.c

................................................................................
                          *******  HydroBase *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             1993
                             updated to ANSI C 1999
................................................................................
____________________________________________________________________________
  USAGE:  

 propcalc filename(_roots) -P<property_list> [-D<dirname>] [-E<file_extent>] [-W<window/w_incr>]

  list of filenames MUST be first argument!

 -P : list of properties to evaluate; ex: -Ppr/te/th/sa/ox
      -P (by itself) will print a list of the available properties;

 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum
 -W : specifies window width/increment for gradient properties (default is -W100/10)

____________________________________________________________________________
propcalc computes properties at each standard level in a station and 
outputs a (hydrobase format) ascii station file containing all the properties
specified with the -P option.                                                    ____________________________________________________________________________
*/


#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"
#include "hb_gamma.h"
#include "hb_paths.h"


/* input file pathnames */

#define    EXTENT   ""
#define    DIR      ""

    /* flags to indicate which properties are to be output and which are
       needed for computation */
      
short prop_req[MAXPROP];         /* set of props requested for output */
short prop_needed[MAXPROP];      /* set of props requested || defining a
                                    surface || needed for computation */

   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;
double s_pref;           /* ref lev for computing a non-standard sigma level */
double pe_pref;           /* ref lev for computing potential energy */
double ht_pref;           /* ref lev for computing dynamic height */
double ih_pref;           /* ref lev for computing integrated dynamic height */
int window, w_incr;     /* used to specify pr window for gradient properties*/

struct GAMMA_NC ginfo;   /* used for neutral density */

main (int argc, char **argv)
{
   short   popt;
   int     index, nprops = 0;
   int     i;
   int     curfile = 1, nfiles = 0; 
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   char    *extent, *dir, *st;
   int     infile;
   void    print_usage(char *);
   void    get_hydro_data(int);
   int     parse_prop_list(char *);


/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    infile = STDIN;
    popt = 0;
    error = 0;
    s_pref = -1;
    window = 100;
    w_incr = 10;

/* initialize these ... */

   for (i = 0; i < MAXPROP; ++i) {
       prop_req[i] = 0;
       prop_needed[i] = 0;
      station.observ[i] = (double *) NULL;
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

               case 'P':
                        popt = 1;
                        if ( argv[i][2] == '\0') {
                               print_prop_menu();
                               exit(0);
                        }
                        nprops = parse_prop_list(&argv[i][2]);
                        if (nprops <= 0)
                             error = 1;
                        break;

               case 'W':
                        error = (sscanf(&argv[i][2],"%d", &window) == 1) ? 0 : 1;
                        st = &argv[i][2];
                        while (*(st++) != '\0') {
                            if (*st == '/') {
                              ++st;
                              error = (sscanf(st,"%d", &w_incr) == 1) ? 0 : 1;
                              break;
                            }
                        }
                        
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

   if (  !popt) {
       fprintf(stderr,"\nYou must specify input file(s) and properties.\n");
       exit(1);
   }

/* loop for each input file */

   do {
   
     if (nfiles) {
       infile = open_hydro_file(dir, argv[curfile], extent, print_msg);
       if (infile < 0)
          goto NEXTFILE;
     }
     else 
          fprintf(stderr,"\nhb_propcalc: Station input expected from stdin...\n");
     
            /* read each file completely */

 
         get_hydro_data(infile);

NEXTFILE:
        if (nfiles)
         close(infile);

   } while (curfile++ < nfiles );

   fprintf(stderr,"\nEnd of %s.\n", argv[0]);
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\n\n***************************************************");
   fprintf(stderr,"\n hb_propcalc computes properties at each pressure level");
   fprintf(stderr,"\n in a station and outputs a HydroBase format"); 
   fprintf(stderr,"\n station file containing all the properties specified");
   fprintf(stderr,"\n with the -P option.");
   fprintf(stderr,"\n***************************************************");

                                                       
   fprintf(stderr,"\nUsage:  %s filename_root(s)  -P<list_of_properties>   [-W<window>[/<w_incr>]] [-D<dirname>] [-E<file_extent>]", program);

   fprintf(stderr,"\n\n  List of filename (roots) must be first arguments");
   fprintf(stderr,"\n-P : properties to list out;");
   fprintf(stderr,"\n          ex:  -Ppr/th/sa/ox/ht");
   fprintf(stderr,"\n               -P (by itself) produces a list of available properties");
   fprintf(stderr,"\n\n   OPTIONS: ");
   fprintf(stderr,"\n-W : Specifies pressure window (db) for computing ");
   fprintf(stderr,"\n          gradient properties (bvfreq, pv...) This ");
   fprintf(stderr,"\n          constitutes the range over which observations");
   fprintf(stderr,"\n          are incorporated into the gradient computation.");
   fprintf(stderr,"\n          The window is centered around each pressure ");
   fprintf(stderr,"\n          level being evaluated.");
   fprintf(stderr,"\n          w_incr specifies how finely to subdivide the");
   fprintf(stderr,"\n          window into increments(db) to create");
   fprintf(stderr,"\n          an evenly spaced pressure series over the window.");
   fprintf(stderr,"\n          defaults: -W100/10 ");
   fprintf(stderr,"\n-D : specifies directory input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n-E : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n-h help...... prints this message.");

   fprintf(stderr,"\n\n");  
   return;
}

/*****************************************************************************/
int parse_prop_list(char *st)

   /*  Parses a list of property mnemonics and sets global flags accordingly.
       Returns the number of properties.  An error will cause an exit.
        
       char *st:    the list with -P stripped off  */
{
   int index, nprops;
   char prop[4];
   double ref_val;

   nprops = 0;

   do {
      if (*st == '/')
         ++st;
      sscanf(st,"%2s", prop);
      index = get_prop_indx(prop);
      if (index < 0) {
         fprintf(stderr,"\n Unknown property requested: %s\n", prop);
         exit(1);
      }
      prop_req[index] = 1;
      prop_needed[index] = 1;
      ++nprops;
      ++st;
      ++st;


      /* !**!  Special cases for properties */

      if (((enum property)index == S_ ) || ((enum property) index == HT) || ((enum property) index == IH) || ((enum property) index == PE)) {
         if (sscanf(st, "%lf", &ref_val) != 1) {
             fprintf(stderr, "\n Specify a ref pressure for  %.2s", prop);
             fprintf(stderr,"\n   ex: -P%.2s1500/th/sa\n", prop);
             exit(1);
         }
         while (!(*st=='/' || *st=='\0' || *st==' '))
            ++st;

        switch ((enum property) index) {
           case S_:
              if (s_pref >= 0) {   /* check for previous request? */
                  fprintf(stderr,"Sorry. You have requested the property s_ with multiple reference levels.\n");
                  fprintf(stderr,"You can only use one of those prefs");
                  fprintf(stderr,"  You will probably have to do multiple runs for each different pref you want to associate with s_\n\n");
                  exit(1);
              }
              s_pref = ref_val;
              break;
           case PE:
              pe_pref = ref_val;
              break;
           case HT:
              ht_pref = ref_val;
              break;
	   case IH:
	      ih_pref = ref_val;
	      break;
           default:
              ;
        }  /* end switch */
      }

   /* end of special cases */

  } while (*st == '/');
  
 /* !**!  Special cases for properties */
  
  if (prop_req[(int)GE] && !prop_req[(int)GN]) {
     prop_req[(int)GN] = 1;
     prop_needed[(int)GN] = 1;
     ++nprops;
  }
  
  if (prop_req[(int)GE] || prop_req[(int)GN]) 
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);
     
 /* end of special cases */

  return (nprops);
}  /* end parse_prop_list() */
 

/*****************************************************************************/
void get_hydro_data(int file) {
/* Reads each station in a HydroBase file
 * and computes property values at each
 * standard level.  This module requires
 * that the HydroBase file contains a
 * minimum of pr, de, te, sa observations.  */

  int error, i, j, nreq, pts_avail, ratio_done;
  int *tmp;
  double dlat, dlon;

  /* Read each station in file ... */
  while ((error = get_station(file, &hdr, &station)) == 0) {

    /* Check if basic properties de, pr, te, and sa are available. */
    if ( ! available(PR, &hdr) && available(DE, &hdr)) {

      /* For the case when pressure is not present but
       * we do have depth, compute the pressure. */

      /* Allocate space for the pressure profile. */
      free_and_alloc(&station.observ[(int)PR], hdr.nobs);

      /* Compute the pressure. */
      for ( j=0; j < hdr.nobs; ++j) {
	station.observ[(int)PR][j] = hb_p80(station.observ[(int)DE][j], (double)hdr.lat);
      }

      /* Reset the property identifiers for this station
       * to include pressure as one of the props. */
      tmp = hdr.prop_id;
      hdr.prop_id = (int *)calloc(++hdr.nprops, sizeof(int));

      for (i=0; i < station.nprops; ++i)
	hdr.prop_id[i] = tmp[i];
      
      hdr.prop_id[i] = (int) PR;
      free(tmp);   
      ++station.nprops;

    } /* End of basic props check. */

    /* Do we have pressure, temp, and salt? */
    pts_avail = (available(PR, &hdr) && available(TE, &hdr) 
                 && available(SA, &hdr));

    ratio_done = 0;

    /* Compute appropriate properties at each level in station. */
    
    /* !**! Special cases for individual properties... */
    for (i = 0; i < MAXPROP; ++i) {
      if (prop_needed[i] && pts_avail && !available((enum property)i, &hdr)) {
	switch ((enum property) i) {
	     case DE:
               free_and_alloc(&station.observ[(int)DE], hdr.nobs);
               for ( j=0; j < hdr.nobs; ++j) {
	          station.observ[(int)DE][j] = hb_depth(station.observ[(int)PR][j], (double)hdr.lat);
	        }
	        break;

             case OX:
               if (available(O2, &hdr)) {
                    free_and_alloc(&station.observ[i], hdr.nobs);
                    for (j=0; j < hdr.nobs; ++j) {
                       station.observ[i][j] = ox_kg2l(station.observ[(int)O2][j], station.observ[(int)PR][j], station.observ[(int)TE][j],station.observ[(int)SA][j]);
                    }
               }
               break;

             case O2:
               if (available(OX, &hdr)) {
                    free_and_alloc(&station.observ[i], hdr.nobs);
                    for (j=0; j < hdr.nobs; ++j) {
                       station.observ[i][j] = ox_l2kg(station.observ[(int)OX][j], station.observ[(int)PR][j], station.observ[(int)TE][j],station.observ[(int)SA][j]);
                    }
               }   
               break;

             case TH:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_theta(hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA] );
               break;

             case TP:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_dtdp(station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], hdr.nobs, window, w_incr);
               break;

             case TZ:
               free_and_alloc(&station.observ[i], hdr.nobs);
	       dlat = (double) hdr.lat;
               compute_dtdz(station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA],hdr.nobs,window,w_incr,dlat);
               break;

             case HC:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_heat_capacity(hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA] );
               break;

             case S0:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(0., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               break;

             case S1:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(1000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               break;

             case S2:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(2000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               break;

             case S3:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(3000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               break;

             case S4:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(4000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               break;

             case S_:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(s_pref, hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               break;

	     case SD:
               free_and_alloc(&station.observ[i], hdr.nobs);
	       compute_dref_sigma(station.observ[(int)DE], hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               break;

             case HT:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_height(hdr.nobs, station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], ht_pref, station.observ[i]);
               break;

	     case IH:
	       free_and_alloc(&station.observ[i], hdr.nobs);
	       compute_htdz_over_f(hdr.nobs, station.observ[(int)DE], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], ih_pref, hdr.lat, station.observ[i]);
	       break;

             case PE:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_energy(hdr.nobs, station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], pe_pref, station.observ[i]);
               break;

             case SV:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sp_vol( hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               break;

             case VA:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_svan( hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
               break;

             case DR:
             case AL:
             case BE:
	        if (! ratio_done ) {
                    free_and_alloc(&station.observ[(int)DR], hdr.nobs);
                    free_and_alloc(&station.observ[(int)AL], hdr.nobs);
                    free_and_alloc(&station.observ[(int)BE], hdr.nobs);
                    compute_ratio( hdr.nobs, station.observ[(int)DR], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], station.observ[(int)AL], station.observ[(int)BE]);
		     ratio_done = 1;
		 }
               break;

             case VS:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sound_vel( station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], hdr.nobs);
               break;

             case PV: 
               free_and_alloc(&station.observ[i], hdr.nobs);

	       /* Get station latitude */
               dlat = (double) hdr.lat;

	       /* Compute buoyancy frequency = sqrt( -g/rho * drho/dz ) */
               buoy_freq(station.observ[i], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], hdr.nobs, window, w_incr);

	       /* Square the buoyancy frequency for input
		* into po_vort. */
               for (j = 0; j < hdr.nobs; ++j) {
                 if (station.observ[i][j] > -99.) 
                   station.observ[i][j] *= station.observ[i][j];
               }
	       /* Compute potential vorticity.  Note that the
		* 1e-12 (ms)^-1 units conversion is done inside
		* po_vort(). */
               po_vort(station.observ[i], station.observ[i], hdr.nobs, dlat);

               break;

	     case RR: 
               free_and_alloc(&station.observ[i], hdr.nobs);

	       /* Get station latitude */
               dlat = (double) hdr.lat;

	       /* Compute Rossby Radius */
	       compute_approx_rossby_radius(station.observ[i], hdr.nobs, hdr.pdr, station.observ[(int)DE], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], dlat, window, w_incr);
	       
               break;
               
             case BF:
               free_and_alloc(&station.observ[i], hdr.nobs);
               buoy_freq(station.observ[i], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], hdr.nobs, window, w_incr);
               for (j = 0; j < hdr.nobs; ++j) {
		 if (station.observ[i][j] > -9998.)
		   /* Convert the units to 1e-5 rad/s
		    * with this extra multiplication. */
		   station.observ[i][j] *= 1.0e5;
               }
               break;

             case GN:
	       if (!prop_req[(int)GE]) {
                  free_and_alloc(&station.observ[i], hdr.nobs);
                  dlat = (double) hdr.lat;
                  dlon = (double) hdr.lon;
		  compute_gamma_n(&ginfo, hdr.nobs, station.observ[i], station.observ[(int)PR],station.observ[(int)TE],station.observ[(int)SA], dlon, dlat);
               }
	       break;
	       
             case GE:
	       if (pts_avail) {
                  free_and_alloc(&station.observ[(int)GE], hdr.nobs);
                  free_and_alloc(&station.observ[(int)GN], hdr.nobs);
                  dlat = (double) hdr.lat;
                  dlon = (double) hdr.lon;
		  compute_gamma_nerr(&ginfo, hdr.nobs, station.observ[(int)GN],   
		     station.observ[(int)GE], station.observ[(int)PR],
		     station.observ[(int)TE], station.observ[(int)SA], dlon,
		     dlat);
               }
	       break;

             default:
               break;
	} /* End property calculation switch */
      } /* End property needed check */
     } /* End for loop over all properties */
     
    /* Determine availability of requested properties... */
    nreq = 0;
    for (i = 0; i < MAXPROP; ++i) {
      if (prop_req[i] && station.observ[i] != NULL )
	++nreq;
    }
    
    /* Update the station header with requested
     * properties and then output the station. */
    hdr.nprops = station.nprops = nreq;
    free((void *)hdr.prop_id);
    hdr.prop_id = (int *) calloc((size_t)nreq, sizeof(int));
    j = 0;
    for (i = 0; i < MAXPROP; ++i) {
      if (! prop_req[i]) {
	free(station.observ[i]);
	station.observ[i] = NULL;
      }
      else if (station.observ[i] != NULL)  {
	hdr.prop_id[j++] = i;
      }
    }
    write_hydro_station(STDOUT, &hdr, &station);
    
    /* Clear memory of station data. */
    for (i = 0; i < MAXPROP; ++i) {
      if (station.observ[i] != NULL) {
	free((void *)station.observ[i]);
	station.observ[i] = NULL;
      }
    } /* End of loop over all properties to clear memory. */  
  }  /* End while loop over all stations. */
  
  if (error > 0)
    report_status(error, stderr);
  
  return;

} /* end get_hydro_data() */

/***********************************************************/
