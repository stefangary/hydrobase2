/*  hb_geostroph_section.c

................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth  Curry
                             Woods Hole Oceanographic Inst
                             2001
................................................................................
____________________________________________________________________________
  USAGE:  

 hb_geostroph_section filename(_roots) -B<bin_file> -R<ref_lev_code> -V<output_velocity_file> -T<output_transport_file> -Z<std_lev_file> [-U<ref_lev_velocity>] [-L] [-N] [-O<option>] [-X<la|lo|di>] [-D<dirname>] [-E<file_extent>]

  input HydroBase2 station format file with ordered consecutive profiles.

 -R : Reference level option:  2-char property or bo (for bottom);
 -Z : file specifying pressure series for computing dynamic heights;

[-B]: file that specifies bins for computing transports (prop/min/max)
      If -B is not specified, but -T is, the bins will be the pressure series
      specified with -Z.
[-D]: specifies directory for input files (default is current dir)  
      ex: -D/d1/hbase/data
[-E]: specifies input_file extent (default is no extent) ex: -E.sum
[-L]: add list of bins to output transport file. 
[-N]: output distance in nm (km is default)
[-T]: name of output transport file.
[-U]: ref level velocity (default is 0 cm/sec). 
[-V]: name of output velocity file.
[-X]: output either lat | lon | distance. default is all 3
____________________________________________________________________________
hb_geostroph_section computes velocity and/or transports for pairs of 
consecutive stations.   
                                                    
*/


#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"
#include "hb_gamma.h"
#include "hb_paths.h"

#define  NMperKM  .539593
#define   RADperDEG 0.017453292             /* pi/180 */
#define  EarthRadius  3437.746873       /* in nm */


/* input file pathnames */

#define    EXTENT   ""
#define    DIR      ""


double *distance;
double *stdplevs;
int    min, kilo, nplevs;
float prev_lat, prev_lon;
double cumdist;
double FLAG;

struct GAMMA_NC ginfo;   /* used for neutral density */

   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;
double s_pref;           /* ref lev for computing a non-standard sigma level */

   /* prototypes for locally defined functions  ... */
   
void    print_usage(char *);
void    get_hydro_data(int);
int     parse_r_opt(char *);
int     parse_x_opt(char *);
int     get_std_levs(char *);
double  get_x_prop(int);
void get_prop(int, int);

main (int argc, char **argv)
{
   int     index, nsta;
   int     xopt, ropt, vopt, topt, zopt;
   int     i, max, sta;
   int     curfile = 1;
   int     nfiles = 0;
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   char    *extent, *dir;
   int     infile;
   float   depth;
   FILE    *velfile, *transfile, *zfile;

/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    FLAG = HB_MISSING - 0.1;
    zopt = xopt = ropt = topt = vopt = 0;
    kilo = 1;
    error = 0;
    s_pref = -1;
    infile = STDIN;
    transfile = velfile = NULL;
    cumdist = 0;

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
               case 'N':
                        kilo = 0;
                        break;
               case 'R':            /* reference level */
                        ropt = 1;
			error = parse_r_opt(&argv[i][2]);
                        break;
			
               case 'T':            /* transport file */
			
                        topt = 1;
			transfile = fopen(&argv[i][2], "w");
                        if (transfile == NULL) {
                           fprintf(stderr,"\nUnable to open %s",  &argv[i][2]);
                           exit(1);
                        }
                        fprintf(stderr,"\nOpened %s",  &argv[i][2]);
                        break;
                case 'V':            /* velocity file */
			
                        vopt = 1;
			velfile = fopen(&argv[i][2], "w");
                        if (velfile == NULL) {
                           fprintf(stderr,"\nUnable to open %s",  &argv[i][2]);
                           exit(1);
                        }
                        fprintf(stderr,"\nOpened %s",  &argv[i][2]);
                        break;

               case 'X':
                        xopt = parse_x_opt(&argv[i][2]);
			if (xopt < 0) 
			      error = 1;
			break;
               case 'Z':
	                zopt = 1;
	                nplevs = get_std_levels(&argv[i][2]);
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

   if (! (ropt && zopt)  ) {
       fprintf(stderr,"\nYou must specify -R and -Z arguments.");
       fprintf(stderr,"\nFor further help, type: %s -h   \n", argv[0]);
       exit(1);
   }

   if (ref_lev_index == (int)GN) 
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);
     
   if (nfiles > 0) {
          infile = open_hydro_file(dir, argv[curfile], extent, print_msg);
          if (infile < 0)
             exit(1);
   }
 
   error = get_station(infile, &hdr, &station);
   
   if (error) {
      report_status(error);
      exit(1);
   }
   get_dyn_ht_profile(&hdr, &station, dh1, z1, &nz1);
   
   while ((error = get_station(infile, &hdr, &station)) == 0) {
      get_dyn_ht_profile(&hdr, &station, dh2, z2, &nz2);
      dist = compute_distance(lat1, lon1, lat2, lon2);
      get_shear_profile(lat1, lon1, lat2, lon2, dh1, dh2, z1, z2, nz1, nz2, vout, zout, &nzout, &dcl);

      if (vopt) {
         /* output the avlat, avlon, cumdist, zout, vout */
      } 
      
      if (topt) {
      
      } 
      
          
      z1 = z2;
      dh1 = dh2;
      lat1 = lat2;
      lon1 = lon2;
   }

   if (error) {
      report_status(error);
      exit(1);
   }
   
   fprintf(stderr,"\n\nEnd of %s\n", argv[0]);
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nComputes geostrophic velocity and/or transports");
   fprintf(stderr,"\nalong a hydrographic section.");

    fprintf(stderr,"\n\nUSAGE:  %s filename(_roots) -R<ref_lev_type[/<value>]> -Z<std_lev_file> [-B<binfile>] [-D<infile_dir>] [-E<infile_extent>] [-L] [-N] [-T<transport_outfile>] [-U<ref_lev_velocity>] [-V<velocity_outfile>] [-X<la|lo|di>] ", program);
   fprintf(stderr,"\n\n  List of filenames must be first argument or input is expected to come from stdin");
   fprintf(stderr,"\n-R : reference level:  [bo]ttom OR 2 char mnemonic / value");  
   fprintf(stderr,"\n     ex: -Rbo        use bottom as reference level. ");
   fprintf(stderr,"\n     ex: -Rpr/1000   pressure = 1000 db for ref lev. ");
   fprintf(stderr,"\n     ex: -Rth/1.9    theta = 1.9 degC for ref lev. ");
   fprintf(stderr,"\n-Z : file of pressure series for computing dyn height profiles.");
   fprintf(stderr,"\n\n   OPTIONS: ");
   fprintf(stderr,"\n[-D] : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n[-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n[-L] : list bins in the transport output.");
   fprintf(stderr,"\n[-N] : distances are in nautical miles (default is km).");
   fprintf(stderr,"\n[-T] : name of output file for transports."); 
   fprintf(stderr,"\n[-U] : velocity at reference level OR name of file containing ref lev velocity for each pair of stations.");  
    
   fprintf(stderr,"\n[-V] : name of output file for velocities.");  
   fprintf(stderr,"\n[-X] : output x-property: [la]t,[lo]n, or [di]stance;");
   fprintf(stderr,"\n          ex:  -Xlo");
   fprintf(stderr,"\n       Default is to output all three properties.");
   
   
   fprintf(stderr,"\n-h : help...... prints this message. \n");

   fprintf(stderr,"\n\n");  
   return;
}

/*****************************************************************************/
int parse_r_opt(char *st)
    /*  Parses property mnemonics and sets global flags accordingly.
        Returns 0 if successful, or -1 for an error. */
{
   char prop[4];

/* empty argument, print out menu */
   
   if ( *st == '\0') {
      fprintf(stderr,"\nReference Level options:\n");
      fprintf(stderr,"\nbo : bottom");
      print_prop_menu();
      exit(0);
   }
			
   

   if (*st == '/')
      ++st;
      
      
/* bottom option...  */

   if (*st == 'b' || *st == 'B') {
      if (*(st+1) == 'o' || *(st+1) == 'O') {
         ref_lev_index = -1;
	 return;
      }
   }
   
/* property option...  */
      
   sscanf(st,"%2s", prop);
   ref_lev_index = get_prop_indx(prop);
   if (ref_lev_index < 0) {
         fprintf(stderr,"\n Unknown property requested: %s\n", prop);
         return(-1);
   }
   
   ++st;
   ++st;
   if (*st == '/')
      ++st;
      
      
  /* !**! Special cases for individual properties ... */

   /* properties requiring reference levels */
   
   if (((enum property)ref_lev_index == S_ ) )) {
         
	 if (sscanf(st, "%lf", &s_pref) != 1) {
             fprintf(stderr, "\n Specify a ref pressure for  %.2s", prop);
             fprintf(stderr,"\n   ex: -R%.2s1500/34.65\n", prop);
             exit(1);
         }
         while (!(*st=='/' || *st=='\0' || *st==' '))
            ++st;
   }

   if (*st == '/')
      ++st;
   if (sscanf(st,"%lf", &ref_lev_value) != 1)
      return(-1);
      
      
   return(0);
}  /* end parse_r_opt() */

/*****************************************************************************/ 
int parse_x_opt(char *st)

{
int xopt;

   if ( *st == '\0') {
      fprintf(stderr,"\nX-Axis options:\n");
      fprintf(stderr,"   la : latitude");
      fprintf(stderr,"   lo : longitude");
      fprintf(stderr,"   di : distance\n");
      exit(0);
   }
			
   switch (*st) {
       case 'l':
            switch (*(++st)) {
               case 'a':
                  xopt = 1;
                  break;
               case 'o':
                  xopt = 2;
                  break;
               default:
                  xopt = -1;
            } /* end switch */
            break;
       case 'd':
            switch (*(++st)) {
               case 'i':
                   xopt = 3;
                   break;
               default:
                   xopt = -1;
            } /* end switch */
            break;
	    
       default:
            xopt = -1;
                           
   } /* end switch */
   
   return(xopt);
}  /* end parse_x_opt() */
/*****************************************************************************/
int get_std_levels(char *st)
   /* Opens the named file, reads a list of pressure values into the global
      array stdplevs, and closes the file. Returns the number of values + 1.
   */ 
{
int error, n;
FILE *fptr;


   if ( *st == '\0') {
      fprintf(stderr,"\nNo file name for pressure levels was specified.\n");
      exit(1);
   }

   fptr = fopen(st,"r");
   if (fptr == NULL) {
      fprintf(stderr,"\nUnable to open %s for reading.\n", st);
      exit(1);
   }			
   fprintf(stderr,"\nOpened %s", st);
   
   n = 0;
   while (fscanf(fptr,"%lf", stdplevs) != EOF)
      ++n;

   if (n == 0) {
      fprintf(stderr,"\nUnable to read any pressure values\n");
      exit(1);
   }
         
   rewind(fptr);
   ++n;
   
   stdplevs = (double *) calloc(n, sizeof(double));

   n = 0;   
   while (fscanf(fptr,"%lf", &stdplevs[n]) != EOF)
      ++n;
      
   fclose(fptr);   
   fprintf(stderr,"\n   ... Read %d pressure levels", n);
   return(++n);
   
}  /* end get_std_levels() */
/*****************************************************************************/
void get_prop(int i, int main_props_avail)

      /* computes, if necessary and able, a property at each level in station ... 
      
                   i:  index corresponding to property 
    main_props_avail:  > 0 (true) if pr, te, sa are available 
    
    */
{
   int j;
   double dlat, dlon;

/* !**! Special cases for individual properties... */

       if ( !available((enum property)i, &hdr) && main_props_avail) {
          switch ((enum property) i) {
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
             case HT:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_height(hdr.nobs, station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], ht_pref, station.observ[i]);
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

             case PV: 
               free_and_alloc(&station.observ[i], hdr.nobs);
               dlat = (double) hdr.lat;
               buoy_freq(station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], hdr.nobs, window, w_incr);
               for (j = 0; j < hdr.nobs; ++j) {
                 if (station.observ[i][j] > -99.)  
                    station.observ[i][j] *= station.observ[i][j];
               }
               po_vort(station.observ[i],station.observ[i], hdr.nobs, dlat);

               break;

             case RR: 
               free_and_alloc(&station.observ[i], hdr.nobs);
               dlat = (double) hdr.lat;
	       compute_approx_rossby_radius(station.observ[i], hdr.nobs, hdr.pdr, station.observ[(int)DE], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], dlat, window, w_incr);

               break;

	       
             case BF:
               free_and_alloc(&station.observ[i], hdr.nobs);
               buoy_freq(station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], hdr.nobs, window, w_incr);
               for (j = 0; j < hdr.nobs; ++j) {
                    if (station.observ[i][j] > -99.)   /* convert the units */
                        station.observ[i][j] *= 1.0e5;
               }
               break;

	       
             case GN:
               free_and_alloc(&station.observ[i], hdr.nobs);
               dlat = (double) hdr.lat;
               dlon = (double) hdr.lon;
	       compute_gamma_n(&ginfo, hdr.nobs, station.observ[i], station.observ[(int)PR],station.observ[(int)TE],station.observ[(int)SA], dlon, dlat);
	       break;
	       
             case GE:
               free_and_alloc(&station.observ[(int)GE], hdr.nobs);
               free_and_alloc(&station.observ[(int)GN], hdr.nobs);
               dlat = (double) hdr.lat;
               dlon = (double) hdr.lon;
	       compute_gamma_nerr(&ginfo, hdr.nobs, station.observ[(int)GN],   
		     station.observ[(int)GE], station.observ[(int)PR],
		     station.observ[(int)TE], station.observ[(int)SA], dlon,
		     dlat);
	       break;

             default:
               break;
          } /* end switch */
       } /* end if */
     
  return;
} /* end get_prop() */
/*****************************************************************************/
double get_x_prop(int option)

   /* option:  1 = lat, 2 = lon, 3 = distance */
{
static int n = 0;
double dist, dx, dy;

   switch (option) {
      case 1:
         return ((double)hdr.lat);
      case 2:
         return ((double)hdr.lon);
      case 3:
         if (sopt) {
           return(distance[hdr.station-min]);
         }
         
         if (n > 0) {
           dy = (double) ABS(hdr.lat - prev_lat);
           dx = cos((double)(hdr.lat + prev_lat) *.5 * RADperDEG) * ABS(hdr.lon - prev_lon);
           dist = RADperDEG * EarthRadius * sqrt(dx * dx + dy * dy);  /* in nautical miles */
           if (kilo)
             dist /= NMperKM; 
           cumdist += dist;
         }
         prev_lat = hdr.lat;
         prev_lon = hdr.lon;
         ++n;
         return (cumdist);
      default:
        return (-1);
   } /* end switch */
}  /* end get_x_prop() */
/*****************************************************************************/
void get_hydro_data(int file)

   /*  Reads each station in a HydroBase file and computes property values
       at each standard level.  This module requires that the HydroBase
       file contains a minimum of pr, de, te, sa observations.  */
{
   int error, i, j, mainprops;
   double x;


/* read each station in file ... */

    while ((error = get_station(file, &hdr, &station)) == 0) {

  /* check that pr, de, te, and sa are available ... */
       mainprops = 1;
       if (!(available(PR, &hdr) && available(TE, &hdr) && available(SA, &hdr))) {
         mainprops = 0;
       }
       
      get_prop(yopt, mainprops);
      get_prop(zopt, mainprops);
      x = get_x_prop(xopt);

    
 /* output the station */

    if ((station.observ[zopt] != NULL) &&  (station.observ[yopt] != NULL) ) {
       for (i = 0; i < hdr.nobs; ++i) {
          if ((station.observ[zopt][i] > FLAG) && (station.observ[yopt][i] > FLAG) )
             fprintf(stdout,"%10.3lf %10.4lf %10.4lf\n", x, station.observ[yopt][i], station.observ[zopt][i]);
       }
    }
    for (i = 0; i < MAXPROP; ++i) {
       if (station.observ[i] != NULL) {
          free(station.observ[i]);
          station.observ[i] = NULL;
       }
    }    
   }  /* end while */

   if (error > 0)
      report_status(error, stderr);

   return;

} /* end get_hydro_data() */

