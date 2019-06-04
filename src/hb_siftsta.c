/*  hb_siftsta.c

................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             1993
................................................................................
  Evaluates a specified property and prints out any station which
  contains values falling between the specified minima and maxima.  
  Designed to search for errant data.
____________________________________________________________________________
  USAGE:  

 hb_siftsta filename(_roots) -P<property/min/max> [-O<outcast file>]  [-Zzmin/max] [-Ttmin/tmax] [-D<dirname>] [-E<file_extent>]

  list of filenames MUST be first argument!

 -P : property to evaluate with min/max values to search for; 
        ex: -Ps2/40/100
      -P (by itself) will print a list of the available properties;
      
 -O : specifies file for data that do NOT meet criteria 
 -S : specifies min/max sigma-0 range
 -T : specifies theta range
 -Z : specifies min/max depth range
 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum

____________________________________________________________________________
____________________________________________________________________________
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
int window, w_incr;     /* used to specify pr window for gradient properties*/
double xmin, xmax;
double tmin, tmax;    /* theta limits over which to search */
double zmin, zmax;    /* depth limits over which to search */
double sigmin, sigmax;    /* sigma limits over which to search */
int propindex;
int count;
int o_flag, t_flag, s_flag;
int outcastfile;

struct GAMMA_NC ginfo;   /* used for neutral density */

  /* prototypes for locally defined functions */
  
void print_usage(char *);
int parse_prop_list(char *);
void get_hydro_data(int);

main (int argc, char **argv)
{
   short   popt;
   int     i;
   int     curfile = 1, nfiles = 0; 
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   char    *extent, *dir, *st;
   int     infile;


/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    popt = 0;
    t_flag = s_flag = o_flag = 0;
    error = 0;
    s_pref = -1;
    window = 100;
    w_incr = 10;
    count = 0;
    zmin = 0;
    zmax = 10000;
    tmin = -3;
    tmax = 40;
    sigmin = 0;
    sigmax = 100;

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

               case 'O':
                        o_flag = 1;
                        outcastfile = create_hydro_file(&argv[i][2], NOCLOBBER);
                        if (outcastfile < 0) {
                          fprintf(stderr, "\nUnable to open %s for writing.  It may already exist?\n", &argv[i][2]);
                          exit(1);
                        }
                        break;
               case 'P':
                        popt = 1;
                        if ( argv[i][2] == '\0') {
                               print_prop_menu();
                               exit(0);
                        }
                        error = parse_prop_list(&argv[i][2]);
                        break;

                case 'S':
                        error = (sscanf(&argv[i][2],"%lf/%lf", &sigmin, &sigmax) == 2) ? 0 : 1;
                        prop_needed[(int)S0] = 1;
                        s_flag = 1;
                        break;
                      
                case 'T':
                        error = (sscanf(&argv[i][2],"%lf/%lf", &tmin, &tmax) == 2) ? 0 : 1;
                        prop_needed[(int)TH] = 1;
                        t_flag = 1;
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

                case 'Z':
                        error = (sscanf(&argv[i][2],"%lf/%lf", &zmin, &zmax) == 2) ? 0 : 1;
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
       fprintf(stderr,"\nYou must specify a property with -P.\n");
       exit(1);
   }

  if (prop_req[(int)GE] || prop_req[(int)GN]) 
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);

   prop_needed[(int)TH] = 1;
   prop_needed[(int)S0] = 1;
        
/* loop for each input file */

   do {

     infile = STDIN;
     if (nfiles > 0) {
       infile = open_hydro_file(dir, argv[curfile], extent, print_msg);
       if (infile < 0)
          goto NEXTFILE;
     }

     
            /* read each file completely */

 
         get_hydro_data(infile);
         fprintf(stderr,"\nNumber of stations that fit criteria: %d\n", count);

NEXTFILE:

     if (nfiles > 0) 
         close(infile);

   } while (curfile++ < nfiles );

   fprintf(stderr,"\n\nEnd of hb_siftsta.\n");
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\n%s determines whether a station has a value of", program);
   fprintf(stderr,"\na property which falls between specified minima/maxima.");
   fprintf(stderr,"\nIf so, it prints the station to stdout.\n");
   fprintf(stderr,"\nUsage:  %s filename_root(s)  -P<property/min/max>  [-D<dirname>] [-E<file_extent>] [-O<outcast_file>]] [-W<window>[/<w_incr>]][-Z<zmin>/<zmax>] [-S<sigmin>/<sigmax>] [-T<tmin>/<tmax>] [-h]", program);

   fprintf(stderr,"\n\n  List of filenames MUST be first argument!");
   fprintf(stderr,"\n-P : property mnemonic/min/max;");
   fprintf(stderr,"\n          ex:  -Psa/0/30");
   fprintf(stderr,"\n               -P (by itself) produces a list of available    properties");
   fprintf(stderr,"\n\n   OPTIONS: ");
   fprintf(stderr,"\n-D : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n-E : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n-O  : filename for stations which do NOT meet criteria;");
   fprintf(stderr,"\n-S : Optional sigma-0 limits min/max to search;");
   fprintf(stderr,"\n          ex:  -S27/27.8");
   fprintf(stderr,"\n-T : Optional theta limits min/max to search;");
   fprintf(stderr,"\n          ex:  -T3/4");
   fprintf(stderr,"\n-Z : Optional pressure limits min/max to search;");
   fprintf(stderr,"\n          ex:  -Z0/1000");
   fprintf(stderr,"\n-W : Specifies pressure window (db) for computing ");
   fprintf(stderr,"\n          gradient properties (bvfreq, pv...) This ");
   fprintf(stderr,"\n          constitutes the range over which observations");
   fprintf(stderr,"\n          are incorporated into the gradient computation.");
   fprintf(stderr,"\n          The window is centered around each pressure ");
   fprintf(stderr,"\n          level being evaluated.");
   fprintf(stderr,"\n          w_incr specifies how finely to subdivide the");
   fprintf(stderr,"\n          window into increments(db) to create");
   fprintf(stderr,"\n          an evenly spaced pressure series over the window.");
   fprintf(stderr,"\n          defaults: -W100/10");
   fprintf(stderr,"\n-h :  help...... prints this message. \n");

   fprintf(stderr,"\n\n");  
   return;
}

/*****************************************************************************/
int parse_prop_list(char *st)
    /*  Parses the  property mnemonic, min, max values and sets global 
        variables accordingly.
        Returns 0.  An error will cause an exit. */
{
   int  error;
   char prop[4];
   double ref_val;
  
      if (*st == '/')
         ++st;
      sscanf(st,"%2s", prop);
      propindex = get_prop_indx(prop);
      if (propindex < 0) {
         fprintf(stderr,"\n Unknown property requested: %s\n", prop);
         exit(1);
      }
      prop_req[propindex] = 1;
      prop_needed[propindex] = 1;
      ++st;
      ++st;


      /* !**!  Special cases for properties */

      if (((enum property)propindex == S_ ) || ((enum property) propindex == PE) || ((enum property) propindex == HT)) {
         if (sscanf(st, "%lf", &ref_val) != 1) {
             fprintf(stderr, "\n Specify a ref pressure for  %.2s", prop);
             fprintf(stderr,"\n   ex: -P%.2s1500/34.35/34.5\n", prop);
             exit(1);
         }
         while (!(*st=='/' || *st=='\0' || *st==' '))
            ++st;

         switch ((enum property) propindex) {
           case S_:
              s_pref = ref_val;
              break;
           case PE:
              pe_pref = ref_val;
              break;
           case HT:
              ht_pref = ref_val;
              break;
           default:
              ;
        }  /* end switch */
      }

   /* end of special cases */

   if (*st == '/')
     ++st;
   error = (sscanf(st,"%lf/%lf", &xmin, &xmax) == 2) ? 0 : 1;
   if (error) {
      fprintf(stderr,"\nUnable to parse min/max values in -P option: %s\n", st);
      exit(1);
   }

  return (0);
}  /* end parse_prop_list() */
 

/*****************************************************************************/
void get_hydro_data(int file)
   /*  Reads each station in a HydroBase file and computes property values
       at each standard level.  This module requires that the HydroBase
       file contains a minimum of pr, de, te, sa observations.  */
{
   int error, i, j, nreq, found, pts_avail;
   int input_s0, input_th;
   int *new_propid;
   double dlat, dlon;
   void  free_and_alloc();


/* read each station in file ... */

    count = 0;
    
    while ((error = get_station(file, &hdr, &station)) == 0) {


  /* check if basic properties pr, te, and sa are available ... */

      pts_avail = (available(PR, &hdr) && available(TE, &hdr) 
                 && available(SA, &hdr));

/* expand hdr.prop_id to hold new properties */
      new_propid = (int *) calloc((size_t)MAXPROP, sizeof(int));
      input_s0 = 0;
      input_th = 0;
      
      for (i=0; i < hdr.nprops; ++i) {
         new_propid[i] = hdr.prop_id[i];
	 if (hdr.prop_id[i] == (int) S0)
	     input_s0 = 1;
	 if (hdr.prop_id[i] == (int) TH)
	     input_th = 1;
      } 
      free(hdr.prop_id);
      hdr.prop_id = new_propid;
      
 /* compute appropriate properties at each level in station ... */

/* !**! Special cases for individual properties... */

   for (i = 0; i < MAXPROP; ++i) {
       if (!available((enum property)i, &hdr) && prop_needed[i] && pts_avail) {
          switch ((enum property) i) {
             case OX:
               if (available(O2, &hdr)) {
                    free_and_alloc(&station.observ[i], hdr.nobs);
                    for (j=0; j < hdr.nobs; ++j) {
                       station.observ[i][j] = ox_kg2l(station.observ[(int)O2][j], station.observ[(int)PR][j], station.observ[(int)TE][j],station.observ[(int)SA][j]);
                    }
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               }
               break;
             case O2:
               if (available(OX, &hdr)) {
                    free_and_alloc(&station.observ[i], hdr.nobs);
                    for (j=0; j < hdr.nobs; ++j) {
                       station.observ[i][j] = ox_l2kg(station.observ[(int)OX][j], station.observ[(int)PR][j], station.observ[(int)TE][j],station.observ[(int)SA][j]);
                    }
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               }
               break;
             case TH:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_theta(hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA] );
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case TP:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_dtdp(station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], hdr.nobs, window, w_incr);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case TZ:
               free_and_alloc(&station.observ[i], hdr.nobs);
               dlat = (double) hdr.lat;
               compute_dtdz(station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], hdr.nobs, window, w_incr, dlat);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case HC:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_heat_capacity(hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA] );
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case S0:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(0., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case S1:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(1000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case S2:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(2000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case S3:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(3000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case S4:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(4000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case S_:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(s_pref, hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case HT:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_height(hdr.nobs, station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], ht_pref, station.observ[i]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case PE:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_energy(hdr.nobs, station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], pe_pref, station.observ[i]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case SV:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sp_vol( hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;

             case VA:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_svan( hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;

             case DR:
                   free_and_alloc(&station.observ[i], hdr.nobs);
		   for (j = 0; j < hdr.nobs; ++j ) {
                      station.observ[i][j] = hb_ratio(station.observ[(int)PR][j], station.observ[(int)TE][j],station.observ[(int)SA][j]);
		   }
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case AL:
                   free_and_alloc(&station.observ[i], hdr.nobs);
		   for (j = 0; j < hdr.nobs; ++j ) {
                      station.observ[i][j] = hb_alpha(station.observ[(int)PR][j], station.observ[(int)TE][j],station.observ[(int)SA][j]);
		   }
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case BE:
                   free_and_alloc(&station.observ[i], hdr.nobs);
		   for (j = 0; j < hdr.nobs; ++j ) {
                      station.observ[i][j] = hb_beta(station.observ[(int)PR][j], station.observ[(int)TE][j],station.observ[(int)SA][j]);
		   }
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case PV: 
               free_and_alloc(&station.observ[i], hdr.nobs);
               dlat = (double) hdr.lat;
               buoy_freq(station.observ[i], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], hdr.nobs, window, w_incr);
               for (j = 0; j < hdr.nobs; ++j) {
                 if (station.observ[i][j] > -99.) 
                   station.observ[i][j] *= station.observ[i][j];
               }
               po_vort(station.observ[i], station.observ[i], hdr.nobs, dlat);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;

               break;

             case RR: 
               free_and_alloc(&station.observ[i], hdr.nobs);
               dlat = (double) hdr.lat;
	       compute_approx_rossby_radius(station.observ[i], hdr.nobs, hdr.pdr, station.observ[(int)DE], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], dlat, window, w_incr);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;

               break;

	       
             case BF:
               free_and_alloc(&station.observ[i], hdr.nobs);
               buoy_freq(station.observ[i], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], hdr.nobs, window, w_incr);
               for (j = 0; j < hdr.nobs; ++j) {
                    if (station.observ[i][j] > -9998.)   /* convert the units */
                        station.observ[i][j] *= 1.0e5;
               }
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;

             case GN:
	       if (!prop_req[(int)GE]) {
                  free_and_alloc(&station.observ[i], hdr.nobs);
                  dlat = (double) hdr.lat;
                  dlon = (double) hdr.lon;
	          compute_gamma_n(&ginfo, hdr.nobs, station.observ[i], station.observ[(int)PR],station.observ[(int)TE],station.observ[(int)SA], dlon, dlat);
	       }
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
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
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
		    hdr.prop_id[hdr.nprops] = (int)GN;
		    ++hdr.nprops;
	       break;
	       
             default:
               break;
          } /* end switch */
       } /* end if */
     } /* end for */
     
 /* determine availability of requested properties... */

   nreq = 0;
    for (i = 0; i < MAXPROP; ++i) {
       if (station.observ[i] != NULL )
          ++nreq;
    }
    
 /* check for values in min/max range */
    found = 0;
    if (available( (enum property) propindex, &hdr)) { 
       for (i = 0; i < hdr.nobs; ++i) {
          if ((station.observ[(int)PR][i] >= zmin) && (station.observ[(int)PR][i] <= zmax )) {
             if ((station.observ[(int)TH][i] >= tmin) && (station.observ[(int)TH][i] <= tmax )) {
               if ((station.observ[(int)S0][i] >= sigmin) && (station.observ[(int)S0][i] <= sigmax )) {
                  if ((station.observ[propindex][i] >= xmin) && (station.observ[propindex][i] <= xmax))
                    found = 1;
	        }
              }
           }
       }
    }

       
       if ((propindex != (int)TH) && ! input_th) {
          --nreq;
          free(station.observ[(int)TH]);
          station.observ[(int)TH] = NULL;
       }
 
        if ( (propindex != (int) S0) && ! input_s0) {
          --nreq;
          free(station.observ[(int)S0]);
          station.observ[(int)S0] = NULL;
       }
       hdr.nprops = station.nprops = nreq;
       free(hdr.prop_id);
       hdr.prop_id = (int *) malloc(nreq * sizeof(int));
       j = 0;
       for (i = 0; i < MAXPROP; ++i) {
          if (station.observ[i] != NULL)  {
            hdr.prop_id[j++] = i;
          }
       }
       
       
   if (found) { 
      ++count;   
       /* output the station */
       write_hydro_station(STDOUT, &hdr, &station);
    } /* end if found */
    else {
     if (o_flag) {
        write_hydro_station(outcastfile, &hdr, &station);
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

